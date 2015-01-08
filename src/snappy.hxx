#ifndef desres_msys_snappy_hxx
#define desres_msys_snappy_hxx

#include <algorithm>
#include <vector>

#include <iosfwd>                       // streamsize
#include <boost/iostreams/concepts.hpp> // multichar_input_filter
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/char_traits.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <snappy.h>

namespace bio = boost::iostreams;

namespace desres { namespace msys { 

enum snappy_state {
    s_start = 0,
    s_done,
    s_have_data,
};

struct snappy_header {
    uint64_t        size;
};

class snappy_decompressor : public bio::multichar_input_filter {
public:
    explicit snappy_decompressor(size_t chunk_size = 1024 * 1024) :
    s_state(s_start),
    m_chunk(chunk_size),
    m_dec_offset(0)
    {
    }

    template<typename Source>
    std::streamsize read(Source& src, char *s, std::streamsize n) {
        std::streamsize total = 0;
        std::streamsize tmp;

        while (total < n) {
            if (s_state == s_have_data) {
                tmp = copy(s + total, n);
                n -= tmp;
                total += tmp;

                if (n == 0)
                    break;
            }

            struct snappy_header header;
            std::streamsize h = bio::read(src, (char *)&header, sizeof(struct snappy_header));
            if (h < 0) {
                if (!total)
                    total = -1;
                break;
            }

            m_chunk.resize(header.size);
            std::streamsize have = bio::read(src, (char *)m_chunk.data(), header.size);
            if (have == -1) {
                if (!total)
                    total = -1;
                break;
            }

            bool good = ::snappy::Uncompress(m_chunk.data(), have, &m_dec);
            if (!good)
                return -1;

            m_dec_offset = 0;
            s_state = s_have_data;

            tmp = copy(s + total, n);
            n -= tmp;
            total += tmp;
        }

        return total;
    }

    template<typename Source>
    void close(Source &) {
        s_state = s_start;
    }

private:
    snappy_state s_state;
    std::vector<char> m_chunk;
    std::string m_dec;
    std::streamsize m_dec_offset;

    std::streamsize copy(char *s, std::streamsize have_space) {
        std::streamsize sz = std::min<std::streamsize>(have_space, m_dec.size() - m_dec_offset);

        memcpy(s, m_dec.data() + m_dec_offset, sz);
        m_dec_offset += sz;

        if (m_dec_offset == (std::streamsize)m_dec.size()) {
            s_state = s_start;
            m_dec_offset = 0;
        }

        return sz;
    }
};

class snappy_compressor : public bio::multichar_output_filter {
public:
    explicit snappy_compressor(size_t chunk_size = 1024 * 1024) :
    s_state(s_start),
    m_chunk(chunk_size),
    m_chunk_size(0),
    m_compr_offset(0)
    {
    }

    template<typename Sink>
    std::streamsize write(Sink& dst, const char* s, std::streamsize n) {
        std::streamsize consumed = 0;
        std::streamsize tmp;

        while (consumed < n) {
            if (s_state == s_start) {
                if (m_chunk_size + n < (std::streamsize)m_chunk.size()) {
                    memcpy((char *)m_chunk.data() + m_chunk_size, s, n);
                    m_chunk_size += n;
                    consumed += n;
                } else {
                    compress(dst);
                }
            }

            if (s_state == s_have_data) {
                tmp = copy<Sink>(dst);
                if (tmp < 0) {
                    if (consumed)
                        return consumed;

                    return -1;
                }
            }
        }

        return consumed;
    }

    template<typename Sink>
    void close(Sink &dst) {
        if (s_state == s_have_data)
            copy<Sink>(dst);

        if ((s_state == s_start) && (m_chunk_size > 0)) {
            compress(dst);
            copy<Sink>(dst);
        }

        s_state = s_start;
    }

private:
    snappy_state s_state;
    std::vector<char> m_chunk;
    std::streamsize m_chunk_size;
    std::string m_compr;
    std::streamsize m_compr_offset;

    template<typename Sink>
    void compress(Sink &dst) {
        ::snappy::Compress(m_chunk.data(), m_chunk_size, &m_compr);

        m_compr_offset = 0;
        s_state = s_have_data;
        m_chunk_size = 0;

        struct snappy_header header;
        memset(&header, 0, sizeof(struct snappy_header));

        header.size = m_compr.size();
        bio::write(dst, (char *)&header, sizeof(struct snappy_header));
    }

    template<typename Sink>
    std::streamsize copy(Sink &dst) {
        std::streamsize written = bio::write(dst, m_compr.data() + m_compr_offset, m_compr.size() - m_compr_offset);
        if (written < 0)
            return written;

        m_compr_offset += written;
        if (m_compr_offset == (std::streamsize)m_compr.size()) {
            s_state = s_start;
            m_compr_offset = 0;
        }

        return written;
    }
};


}}

#endif /* __SNAPPY_HPP */
