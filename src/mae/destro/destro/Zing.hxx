/* @COPYRIGHT@ */

#ifndef DESTRO_ZING_HXX
#define DESTRO_ZING_HXX
#include <vector>
#include <string>
#include <stdint.h>
#include <unordered_map>

namespace desres { namespace msys {
    class ZingPool {
        std::unordered_map<std::string, uint32_t> m_map;
        std::vector<const std::string*> m_names;

    public:
        ZingPool();
        uint32_t insert(const std::string& string);
        std::string const& retrieve(uint32_t n) const;
    };

    struct Zing {
        uint32_t index = 0;
        Zing() = default;
        Zing(Zing const&) = default;
        Zing(std::string const& s, ZingPool& pool) : index(pool.insert(s)) {}

        std::string const& string(ZingPool const& pool) const {
            return pool.retrieve(index);
        }

        bool is_empty() const {
            return index==0;
        }
    };
}}
#endif
