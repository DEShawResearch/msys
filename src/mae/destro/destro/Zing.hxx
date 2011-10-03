/* @COPYRIGHT@ */

#ifndef DESTRO_ZING_HXX
#define DESTRO_ZING_HXX
#include <vector>
#include <string>
#include <stdint.h>

namespace desres {
  class ZingPool {

    bool m_inserts_allowed;
    std::vector< std::vector<char> > m_pool;
    std::vector< std::vector<unsigned char> > m_pool_active;
    std::vector< std::string> m_overflow;
    std::vector< unsigned char >  m_overflow_active;

  public:
    static void check_assumptions();
    void dump() const;

    /*! Number of bins in the pool */
    static const size_t s_size;

    /*! Creates a string pool for use with zing's.
     */
    ZingPool();

    /*! Destroy a zing state pool.
     */
    ~ZingPool();

    /*! Add a string to the appropriate pool.

      The routine may return the offset of a previous string.
     */
    size_t insert(const std::string& string);

    /*! Get a string from the pool.
     */
    std::string retrieve(size_t n,size_t length) const;


    /*! Mark all strings as "not-in-use"

      This routine is used to start a pool cleanup step.  At
      this point, all active strings must be touched before
      another insert takes place.  So, after calling clearbits,
      call touch() on each string, then call finish_compaction()
      to allow inserts.
    */
    void clearbits();

    /*! Touch a string to mark that it is still in use.

      Must be called on all active strings after calling clearbits().
    */
    void touch(size_t n,size_t length);

    /*! Call this to allow inserts again after calling clearbits()

      Inserts are not allowed when in compaction mode... Call this
      to signal that the compaction is complete.
    */
    void finish_compaction();

    /*! Compute size of storage for pool
     */
    size_t size() const;
  };

  class Zing {
    uint32_t m_bits;

    static const std::string s_alphabet_amino;
    static const std::string s_alphabet_5bit;
    static const std::string s_alphabet_4bit;
    static const std::string s_alphabet_bcd;

    static const uint32_t s_two_to_the_22nd;
    static const uint32_t s_two_to_the_26th;
    static const uint32_t s_two_to_the_27th;

    static const uint32_t mask_float;
    static const uint32_t control_float;
    static const uint32_t mask_5bit;
    static const uint32_t control_5bit;
    static const uint32_t mask_bcd;
    static const uint32_t control_bcd;
    static const uint32_t mask_7bit;
    static const uint32_t control_7bit;
    static const uint32_t mask_size;
    static const uint32_t control_size;
    static const uint32_t mask_any;
    static const uint32_t control_any;
    static const uint32_t mask_4bit;
    static const uint32_t control_4bit;
    static const uint32_t mask_amino;
    static const uint32_t control_amino;

    static uint32_t add_control(uint32_t x, uint32_t mask, uint32_t control);

    static const uint32_t encode_nbit(const std::string& string,unsigned nbits, unsigned nchars, uint32_t mask, uint32_t control, const std::string& alphabet);

    static const uint32_t encode_float(const std::string& string);
    static const std::string decode_float(uint32_t x);

    static const uint32_t encode_5bit(const std::string& string);
    static const std::string decode_5bit(uint32_t x);

    static const uint32_t encode_4bit(const std::string& string);
    static const std::string decode_4bit(uint32_t x);

    static const uint32_t encode_amino(const std::string& string);
    static const std::string decode_amino(uint32_t x);

    static const uint32_t encode_any(const std::string& string,ZingPool& pool);
    static const std::string decode_any(uint32_t x,const ZingPool& pool);

    static const uint32_t encode_size(const std::string& string,ZingPool& pool);
    static const std::string decode_size(uint32_t x,const ZingPool& pool);

    static const uint32_t encode_bcd(const std::string& string);
    static const std::string decode_bcd(uint32_t x);

    static const uint32_t encode_7bit(const std::string& string);
    static const std::string decode_7bit(uint32_t x);

  public:

    enum category_t {
      ZING_EMPTY=0,
      ZING_7BIT=1,
      ZING_5BIT=2,
      ZING_4BIT=3,
      ZING_AMINO=4,
      ZING_BCD=5,
      ZING_FLOAT=6,
      ZING_SIZE=7,
      ZING_ANY=8};

    static void check_assumptions();

    Zing();
    Zing(const Zing& other);
    Zing(const std::string& string, ZingPool& pool);
    ~Zing();

    Zing& operator=(const Zing& other);
    bool operator==(const Zing& other) const;
    bool operator!=(const Zing& other) const;

    std::string string(const ZingPool& pool) const;
    uint32_t integer() const;

    //! \brief Size is sometimes encoded.
    size_t size() const;

    category_t category() const;

    bool is_empty() const;

    void touch(ZingPool& pool) const;

  };
    
}
#endif
