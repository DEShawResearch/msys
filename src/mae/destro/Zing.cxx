#include "destro/Zing.hxx"

namespace desres { namespace msys {

    ZingPool::ZingPool() {
        insert("<>");
    }

    uint32_t ZingPool::insert(const std::string& name) {
        auto pair = m_map.emplace(name, m_names.size());
        if (pair.second) {
            m_names.push_back(&pair.first->first);
        }
        return pair.first->second;
    }

    std::string const& ZingPool::retrieve(uint32_t n) const {
        return *m_names.at(n);
    }
}}
