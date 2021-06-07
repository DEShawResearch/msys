/* @COPYRIGHT@ */

#ifndef DESTRO_DESTRO_HXX
#define DESTRO_DESTRO_HXX

#include "Zing.hxx"
#include <msys/types.hxx>

#include <map>
#include <deque>
#include <vector>
#include <string>
#include <iostream>
#include <sys/types.h>

#ifdef WIN32
#ifdef _WIN64
 typedef __int64 ssize_t;
#else
 typedef int ssize_t;
#endif

#endif


namespace desres { namespace msys {
  class DestroArray;

  class Destro {

  public:
    //! \brief The maximum precision we allow (>75 can segv sprintf)
    static const short s_max_precision;
    class Tokenizer;

  /*!
   * \brief A temporary that holds an attribute pending cast or assigment
   *
   * Get_attr() calls from Destro return these intermediary
   * Attribute objects.  These are simply convenience classes
   * that leverage casting, assignment, and dereference to
   * "do the right thing" and maintain type safety.
   */
    class Attribute {
      static DestroArray s_empty_array;

      Destro* m_block;
      std::string m_attr;

      void compatible_schema(char type);

      template<class T> void assigni(T value);

      void assign(const char* value);
      void assign(const std::string& value);
      void assign(int value);
      void assign(short value);
      void assign(long value);
      void assign(long long value);
      void assign(unsigned int value);
      void assign(unsigned short value);
      void assign(unsigned long value);
      void assign(bool value);
      void assign(double value);
      void assign(float value);

      Destro* owner();
      const Destro* owner() const;

      Attribute();
    public:

      //! \brief Holds an (immutable) empty attribute value
      static Attribute s_empty_attribute;

      Attribute(Destro* block,const std::string& attr);
      ~Attribute();

      template<class T>
      inline Attribute& operator=(const T& value);

      Attribute& operator=(const char* value);

      //! \brief The type (b, i, r, s, x, e) of this attribute
      char type() const;

      //! \brief The raw string value of this attribute
      std::string value() const;

      //! \brief The doc string associated with this attribute
      std::string doc() const;

      //! \brief Set the doc string for this attribute
      void doc(const std::string&);

      //! \brief Number of digits of floating point precision
      short precision() const;

      //! \brief Set the precision for this attribute
      void precision(short);

      //! \brief True iff attribute is the <> token
      bool is_empty() const;

      //! \brief Directly update the value
      void update(const std::string& value);

      // Conversions
      operator std::string() const;
      operator double() const;
      operator float() const;
      operator long() const;
      operator int() const;
      operator short() const;
      operator unsigned long() const;
      operator unsigned int() const;
      operator unsigned short() const;
      operator bool() const;

      Destro* as_block();
      const Destro* as_block() const;

      //! \brief If the attribute is empty, use the alternate.
      template<class T>
      T or_else(const T& alternate) const {
        if (is_empty()) return alternate;
        return *this;
      }

      //! \brief Can't cast to const char* , so specialize or_else.
      std::string or_else(const char* alternate) const {
        if (is_empty()) return alternate;
        return *this;
      }


      //! \brief comparison
      template<class T>
      bool operator==(T operand) const {
        T left_operand = *this;
        return (left_operand == operand);
      }
    
      //! \brief Comparison specialized for const char*
      bool operator==(const char* operand) const {
        std::string left_operand = *this;
        return (left_operand == operand);
      }

      //! \brief Just returns ! of operator==()
      template<class T>
      bool operator!=(T operand) const {
        return ! operator==(operand);
      }


      //! \brief Get sub-attribute or sub-block
      Attribute operator[](const std::string&);

      //! \brief Get sub-attribute or sub-block
      Attribute operator[](const std::string&) const;

      //! \brief Get sub-attribute or sub-block
      Attribute operator[](const char*);

      //! \brief Get sub-attribute or sub-block
      Attribute operator[](const char*) const;

      //! \brief Get sub-block reference
      Destro& operator[](size_t);

      //! \brief Converts to block reference at end of a [][] chain
      const Destro& operator[](size_t) const;

      //! \brief Get sub-attribute or sub-block
      Attribute operator()(const std::string&) const;

      //! \brief Get sub-attribute or sub-block
      Attribute operator()(const char*);

      //! \brief Get sub-attribute or sub-block
      Attribute operator()(const char*) const;

      //! \brief Converts to block pointer at end of a [][] chain
      Destro* operator->();

      //! \brief Get sub-block pointer
      const Destro* operator->() const;

      //! \brief Convert to Destro (block only)
      operator Destro&() const;

      //! \brief Convert to Destro (block only)
      operator DestroArray&() const;
    };

  protected:
    Destro* m_parent;

    //! \brief Basic interface
    Destro(Destro* parent=NULL);

    static std::string indent(int level,int& newlevel);

    static std::string schema_out(char type, Zing key, Zing doc, const ZingPool& pool);

    //! \brief Make a string suitable for output (add quotes or escapes)
    //static std::string escape(const char* s);

    static void fill_named(Destro& block,Tokenizer& tokenizer);
    static void fill_nameless(Destro& block,Tokenizer& tokenizer);

    static std::string quotify(Zing z,const ZingPool& zpool);

    static Zing zingify(const std::string& value, ZingPool& zpool);

    //! \brief Set an attribute from a bool value
    virtual void set(const std::string& attr,const bool& value);

    //! \brief Set an attribute from a long value
    virtual void set(const std::string& attr,const long& value);

    //! \brief Set an attribute from a int value
    virtual void set(const std::string& attr,const int& value);

    //! \brief Set an attribute from a double value
    virtual void set(const std::string& attr,const double& value);

    //! \brief Set an attribute from a string value
    virtual void set(const std::string& attr,const std::string&value);

    //! \brief Set an attribute from a char* value
    virtual void set(const std::string& attr,const char* value);

    //! \brief Set an attribute from an existing Attribute
    virtual void set(const std::string& attr,const Attribute& value);


    // OK to allow compiler created copy c'tor
    //Destro(const Destro&);

  public:
    static void test_private_parts();

    /*!
     * \brief Takes a stream and returns maestro tokens.
     * This tokenizer is built on streams and uses a small, tight
     * finite state automata to construct a token
     */
    class Tokenizer {

      /*! \brief Actions for the DFA token builder
       *
       */
      typedef enum {
        DONE = 0,
        SKIPWHITE,
        INCOMMENT,
        CHOOSEKIND,
        SINGLECHAR,
        STARTSTRING,
        INSTRING,
        ESCAPE,
        STARTOTHER,
        CONTINUEOTHER
      } action_t;

      //! \brief The current character
      char m_c;

      //! \brief A stream buffer
      std::streambuf *m_buffer;

      //! \brief The actual input stream
      std::istream* m_input;

      //! \brief Delete stream only if owned
      bool m_input_owned;

      //! \brief The current token
      char * m_token;

      //! \brief number of malloc'ed bytes in m_token
      ssize_t max_token_size;

      //! \brief True iff the token is already read
      bool m_isfresh;

      //! \brief Current line in file
      unsigned m_line;

      //! \brief Current location in file
      size_t m_point;

      //! \brief Line where token starts
      unsigned m_tokenline;

      //! \brief Offset where token starts
      size_t m_tokenpoint;

      //! \brief Get current character
      char peek();

      //! \brief Read a new character
      char read();

      //! \brief Get the 1st character in place
      void prime_the_pump();

      //! \brief True if character is a token
      static bool issingle(char c);

      Tokenizer(const Tokenizer&); // No copy c'tor

    public:

      //! \brief Tokenizer based on a stream
      Tokenizer(std::istream& input);

      //! \brief Tokenizer based on a FILE*
      Tokenizer(FILE* input);

      //! \brief Tokenizer based on a string (stringstream)
      Tokenizer(const std::string& input="");

      //! \brief Clean up and release buffers.
      ~Tokenizer();

      //! \brief Current token under cursor
      const char * token(bool ignore_single_character_tokens=false);

      //! \brief Advance to next token
      void next();

      //! \brief Line number associated with current token
      unsigned line() const;

      //! \brief File seek point for current token
      size_t point() const;

      //! \brief Grab a # comment \n from input stream
      std::string optional_comment();

      //! \brief Predict a particular token
      const char * predict(const char * match="");

      //! \brief Predict a value token
      const char * predict_value();

      //! \brief For while(not_a(match)) loops
      bool not_a(const char * match=END_OF_FILE);

      //! \brief True if at end of file
      bool eof() const;

      //! \brief End of file token
      static const char * END_OF_FILE;
    };

    struct schema_t { /*GCOV-IGNORE*/
      char type;
      std::string attr;
      std::string doc;
    };

    //! \brief standard d'tor.
    virtual ~Destro();

    //! \brief Parent object for this destro.
    Destro* parent();

    //! \brief Parent object for this destro.
    const Destro* parent() const;

    //! \brief Progenitor object for this destro.
    Destro& progenitor();

    //! \brief Progenitor object for this destro.
    const Destro& progenitor() const;

    //! \brief Do we have this named attr or block?
    bool contains(const std::string& name) const;

    //! \brief Get zing pool associated with this object
    virtual ZingPool& pool();

    //! \brief Get zing pool associated with this object or return mutablity error
    virtual ZingPool& mutable_pool();

    //! \brief Get zing pool associated with this object
    virtual const ZingPool& pool() const;

    //! \brief Create a string representation
    operator std::string() const;

    //! \brief Add a schema and add a value
    virtual void add_schema_and_value(char type,const std::string& attr,const std::string& doc,const std::string& value);

    //! \brief Grab a attribute with a simpler interface
    virtual Attribute operator[](const std::string&);

    //! \brief Grab a attribute with a simpler interface
    virtual const Attribute operator[](const std::string&) const;

    template<class T>
    inline Attribute operator()(T fetch) const;

    //! \brief True if this is an array
    virtual bool is_array() const;

    //! \brief This makes writing [] accessors more uniform
    virtual Destro* operator->();

    //! \brief This makes writing [] accessors more uniform
    virtual const Destro* operator->() const;

    //! \brief Add a new block (or row) at end
    virtual Destro& append(const Destro& source);

    //! \brief Add a new block (or row) at end
    virtual Destro& append();

    //! \brief Fetch a array by index
    virtual DestroArray& array(size_t i);

    //! \brief Fetch a array by index
    virtual const DestroArray& array(size_t i) const;

    // ---- Pure virtuals...

    //! \brief Number of sub-blocks
    virtual size_t size() const = 0;

    //! \brief Name of this block (may be "")
    virtual std::string name() const = 0;

    //! \brief Rename this block
    virtual void name(const std::string& name) = 0;

    //! \brief Add a schema without adding a value
    virtual void add_schema(char type,const std::string& attr,const std::string& doc="") = 0;

    //! \brief Attributes known by this block
    virtual std::map<std::string,schema_t> schemas() const = 0;

    //! \brief Copy of schema as an inorder vector
    virtual std::vector<schema_t> ordered_schema() const = 0;

    //! \brief Look up the value of an attribute as a string
    virtual std::string get_value(const std::string& attr) const = 0;

    //! \brief Look up the value of an attribute as a string
    virtual char get_type(const std::string& attr) const = 0;

    //! \brief Look up the float precision of an attribute
    virtual int  get_precision(const std::string& attr) const = 0;
    virtual void set_precision(const std::string& attr,int precision) = 0;

    //! \brief Look up the docstring of an attribute
    virtual std::string get_doc(const std::string& attr) const = 0;
    virtual void set_doc(const std::string& attr,const std::string& doc) = 0;

    //! \brief Look up the value of an attribute
    virtual Attribute get_attr(const std::string& attr) = 0;

    //! \brief Look up the value of an attribute
    virtual const Attribute get_attr(const std::string& attr) const = 0;

    //! \brief This is a generic adaptor for set_attr
    template<class T> void setattr(const std::string& attr,T x) {
      set(attr,x);
    }

    //! \brief Unsafe set
    virtual void set_unsafe(const std::string& attr,const std::string& value) = 0;

    //! \brief Create from raw string value and schema
    virtual void set_attr(char type,const std::string& attr,const std::string& value) = 0;

    //! \brief Set an attribute to the <> value
    virtual void clear(const std::string& attr) = 0;

    //! \brief Delete an attribute or block
    virtual void del(const std::string& attr) = 0;

    //! \brief Delete a block by offset
    virtual void del(size_t blockno) = 0;

    //! \brief Create a new subblock
    virtual Destro& new_block(const std::string& name) = 0;

    //! \brief Create a new subarray
    virtual DestroArray& new_array(const std::string& name,size_t num_elements=0) = 0;

    //! \brief Grab a block by offset
    virtual Destro& operator[](ssize_t) = 0;

    //! \brief Grab a block by offset
    virtual const Destro& operator[](ssize_t) const = 0;

    //! \brief Do we have this named subblock?
    virtual bool has_block(const std::string& name) const = 0;

    //! \brief Do we have this named attr?
    virtual bool has_attr(const std::string& name) const = 0;

    //! \brief Does this attr exist *and* have a non-empty value
    virtual bool has_value(const std::string& name) const = 0;

    //! \brief Fetch a block by index
    virtual Destro& block(size_t i) = 0;

    //! \brief Fetch a block by index
    virtual const Destro& block(size_t i) const = 0;

    //! \brief Fetch a subblock by name
    virtual Destro& block(const std::string& name) = 0;

    //! \brief Fetch a subblock by name
    virtual const Destro& block(const std::string& name) const = 0;

    //! \brief Write a destro block to a stream
    virtual void write(std::ostream& os, int level=0) const = 0;
  };

  class DestroNamedBlock : public Destro {
    //! \brief Block name
    Zing m_name;

    // Do not allow compiler created copy c'tor
    DestroNamedBlock(const DestroNamedBlock&);

    
  public:
    DestroNamedBlock(Destro* parent=NULL);
    ~DestroNamedBlock();

    virtual std::string name() const;
    virtual void name(const std::string& name);
  };

  class DestroBlock : public DestroNamedBlock {
    struct key_value_t {
      Zing key;
      Zing value;
      Zing doc;
      char type;
      int8_t precision;

      key_value_t() : type('e'), precision(-1) {}
    };

    //! \brief List of subblocks
    std::vector<Destro*> m_subblocks;

    //! \brief key value pairs
    std::vector<key_value_t> m_data;

    //! \brief schema lookup
    key_value_t* find_schema(const std::string& attr,const ZingPool& pool);
    const key_value_t* find_schema(const std::string& attr,const ZingPool& pool) const;

    ssize_t offset_of_block(const std::string& name) const;

    // Do not allow compiler created copy c'tor
    DestroBlock(const DestroBlock&);


  protected:
    Destro& add_block(const std::string& name);

    void raw_update(key_value_t* location,ZingPool& zpool,std::string value);

  public:
    DestroBlock(Destro* parent=NULL);
    ~DestroBlock();

    virtual size_t size() const;

    virtual void add_schema(char type,const std::string& attr,const std::string& doc="");
    virtual void add_schema_and_value(char type,const std::string& attr,const std::string& doc,const std::string& value);
    virtual std::map<std::string,schema_t> schemas() const;
    virtual std::vector<schema_t> ordered_schema() const;
    virtual std::string get_value(const std::string& attr) const;
    virtual char get_type(const std::string& attr) const;
    virtual int  get_precision(const std::string& attr) const;
    virtual void set_precision(const std::string& attr,int precision);
    virtual std::string get_doc(const std::string& attr) const;
    virtual void set_doc(const std::string& attr,const std::string& doc);
    virtual Attribute get_attr(const std::string& attr);
    virtual const Attribute get_attr(const std::string& attr) const;
    virtual void set_unsafe(const std::string& attr,const std::string& value);
    virtual void set_attr(char type,const std::string& attr,const std::string& value);
    virtual void clear(const std::string& attr);
    virtual void del(const std::string& attr);
    virtual void del(size_t blockno);
    virtual Destro& new_block(const std::string& name);
    virtual DestroArray& new_array(const std::string& name,size_t num_elements=0);
    virtual Attribute operator[](const std::string&);
    virtual const Attribute operator[](const std::string&) const;
    virtual Destro& operator[](ssize_t);
    virtual const Destro& operator[](ssize_t) const;
    virtual bool has_block(const std::string& name) const;
    virtual bool has_attr(const std::string& name) const;
    virtual bool has_value(const std::string& name) const;
    virtual Destro& block(size_t i);
    virtual const Destro& block(size_t i) const;
    virtual Destro& block(const std::string& name);
    virtual const Destro& block(const std::string& name) const;
    virtual void write(std::ostream& os, int level=0) const;
  };

  class DestroTop : public DestroBlock {
    ZingPool m_pool;

    // Do not allow compiler created copy c'tor
    DestroTop(const DestroTop&);

  public:
    DestroTop();

    DestroTop(Tokenizer& tokenizer);

    template<class T>
    DestroTop(T& input) {
      Tokenizer tokenizer(input);
      if (std::string(tokenizer.token()) == "{") {
        fill_nameless(*this,tokenizer);
      } else {
        fill_named(*this,tokenizer);
      }
    }

    virtual ~DestroTop();

    virtual ZingPool& pool();
    virtual const ZingPool& pool() const;
  };

  class Maeff : public DestroTop {
    void init(Tokenizer& tokenizer);

    DestroBlock m_meta;

    static std::string adjustname(const std::string& name);

    // Do not allow compiler created copy c'tor
    Maeff(const Maeff&);

  public:
    Maeff();

    Maeff(Tokenizer& tokenizer);

    template<class T>
    Maeff(T& input)
      : DestroTop(), m_meta(this)
    {
      Tokenizer tokenizer(input);
      init(tokenizer);
    }

    virtual ~Maeff();

    Destro& meta();
    const Destro& meta() const;

    virtual Destro& new_block(const std::string& name);
    virtual void add_schema(char type,const std::string& attr,const std::string& doc="");

    // For block name tricks (f_m_ct == m_ct == ct)
    virtual bool has_block(const std::string& name) const;
    virtual Destro& block(size_t i);
    virtual const Destro& block(size_t i) const;
    virtual Destro& block(const std::string& name);
    virtual const Destro& block(const std::string& name) const;

    virtual void write(std::ostream& os, int level=0) const;
  };

  class DestroArray : public DestroNamedBlock {
    class DestroRow : public Destro {
      size_t m_row;

      DestroArray* owner();
      const DestroArray* owner() const;

      // OK to allow compiler created copy c'tor
      // DestroRow(const DestroRow&);
    public:
      DestroRow(Destro* parent=NULL,size_t row=0);
      ~DestroRow();

      void setrow(size_t i);

      virtual size_t size() const;
      virtual std::string name() const;
      virtual void name(const std::string& name);
      virtual void add_schema(char type,const std::string& attr,const std::string& doc="");
      virtual void add_schema_and_value(char type,const std::string& attr,const std::string& doc,const std::string& value);
      virtual std::map<std::string,schema_t> schemas() const;
      virtual std::vector<schema_t> ordered_schema() const;
      virtual std::string get_value(const std::string& attr) const;
      virtual char get_type(const std::string& attr) const;
      virtual int  get_precision(const std::string& attr) const;
      virtual void set_precision(const std::string& attr,int precision);
      virtual std::string get_doc(const std::string& attr) const;
      virtual void set_doc(const std::string& attr,const std::string& doc);
      virtual Attribute get_attr(const std::string& attr);
      virtual const Attribute get_attr(const std::string& attr) const;
      virtual void set_unsafe(const std::string& attr,const std::string& value);
      virtual void set_attr(char type,const std::string& attr,const std::string& value);
      virtual void clear(const std::string& attr);
      virtual void del(const std::string& attr);
      virtual void del(size_t blockno);
      virtual Destro& new_block(const std::string& name);
      virtual DestroArray& new_array(const std::string& name,size_t num_elements=0);
      virtual Attribute operator[](const std::string&);
      virtual const Attribute operator[](const std::string&) const;
      virtual Destro& operator[](ssize_t);
      virtual const Destro& operator[](ssize_t) const;
      virtual bool has_block(const std::string& name) const;
      virtual bool has_attr(const std::string& name) const;
      virtual bool has_value(const std::string& name) const;
      virtual Destro& block(size_t i);
      virtual const Destro& block(size_t i) const;
      virtual Destro& block(const std::string& name);
      virtual const Destro& block(const std::string& name) const;
      virtual void write(std::ostream& os, int level=0) const;
    };

    struct key_value_t {
      Zing key;
      Zing doc;
      char type;
      int8_t precision;
      std::vector<Zing> values;

      key_value_t() : type('e'), precision(-1) {}
    };

    //! \brief key value pairs
    std::vector<key_value_t> m_data;

    //! \brief schema lookup
    key_value_t* find_schema(const std::string& attr,const ZingPool& pool);
    const key_value_t* find_schema(const std::string& attr,const ZingPool& pool) const;

    //! \brief Convert string to integer
    static ssize_t integer(const std::string& string);

    std::deque<DestroRow> m_rows;

    // Do not allow compiler created copy c'tor
    DestroArray(const DestroArray&);

  public:
    DestroArray(Destro* parent=NULL, size_t num_elements=0);
    ~DestroArray();

    Destro& bulk_row(const std::string& offset, const std::vector<std::string>& row);

    void resize(size_t n);

    template<class T>
    void set_column(const std::string& attr, const T& container);

    template<class T>
    void set_column(const std::string& attr, const T* begin, size_t size, size_t stride=1);

    template<class T>
    void column(const std::string& attr, std::vector<T>& values) const;

    template<class T>
    void column(const std::string& attr, std::vector<T>& values, const T& defval) const;

    virtual size_t size() const;
    virtual bool is_array() const;

    virtual void add_schema(char type,const std::string& attr,const std::string& doc="");
    virtual std::map<std::string,schema_t> schemas() const;
    virtual std::vector<schema_t> ordered_schema() const;
    virtual std::string get_value(const std::string& attr) const;
    virtual char get_type(const std::string& attr) const;
    virtual int  get_precision(const std::string& attr) const;
    virtual void set_precision(const std::string& attr,int precision);
    virtual std::string get_doc(const std::string& attr) const;
    virtual void set_doc(const std::string& attr,const std::string& doc);
    virtual Attribute get_attr(const std::string& attr);
    virtual const Attribute get_attr(const std::string& attr) const;
    virtual void set_unsafe(const std::string& attr,const std::string& value);
    virtual void set_attr(char type,const std::string& attr,const std::string& value);
    virtual void clear(const std::string& attr);
    virtual void del(const std::string& attr);
    virtual void del(size_t blockno);
    virtual Destro& new_block(const std::string& name);
    virtual DestroArray& new_array(const std::string& name,size_t num_elements=0);
    virtual Attribute operator[](const std::string&);
    virtual const Attribute operator[](const std::string&) const;
    virtual Destro& operator[](ssize_t);
    virtual const Destro& operator[](ssize_t) const;
    virtual bool has_block(const std::string& name) const;
    virtual bool has_attr(const std::string& name) const;
    virtual bool has_value(const std::string& name) const;
    virtual Destro& block(size_t i);
    virtual const Destro& block(size_t i) const;
    virtual Destro& block(const std::string& name);
    virtual const Destro& block(const std::string& name) const;
    virtual void write(std::ostream& os, int level=0) const;
  };


  //! \brief Write a destro to a stream
  std::ostream& operator<<(std::ostream& os,const Destro& M);

  //! \brief Write a destro to a stream
  std::ostream& operator<<(std::ostream& os,const Destro::Attribute& A);

  template<class T>
  Destro::Attribute& Destro::Attribute::operator=(const T& value) {
    assign(value);
    return *this;
  }

  /*!
   * Return a special "empty" value.
   * @param attr The attribute or block name
   * @return The associated attribute if it exists (otherwise an empty)
   */
  template<class T>
  Destro::Attribute Destro::operator()(T attr) const {
    if (!contains(attr)) return Attribute::s_empty_attribute;
    return operator[](attr);
  }


}}
#endif
