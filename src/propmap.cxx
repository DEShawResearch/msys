#include "propmap.hxx"

using namespace desres::msys;

PropertyMap::~PropertyMap() {
    for (Map::iterator i=_map.begin(), e=_map.end(); i!=e; ++i) {
        if (i->second.first==StringType) {
            free(i->second.second.s);
        }
    }
}

std::vector<String> PropertyMap::keys() const {
    std::vector<String> result;
    result.reserve(_map.size());
    for (Map::const_iterator i=_map.begin(), e=_map.end(); i!=e; ++i) {
        result.push_back(i->first);
    }
    return result;
}

bool PropertyMap::has(String const& key) const {
    Map::const_iterator i=_map.find(key);
    return i!=_map.end();
}

void PropertyMap::del(String const& key) {
    Map::iterator i=_map.find(key);
    if (i==_map.end()) MSYS_FAIL("No such property: " << key);
    _map.erase(i);
}

ValueRef PropertyMap::get(String const& key) {
    Map::iterator i=_map.find(key);
    if (i==_map.end()) MSYS_FAIL("No such property: " << key);
    return ValueRef(i->second.first, i->second.second);
}

ValueRef PropertyMap::set(String const& key, ValueType type) {
    Map::iterator i=_map.find(key);
    if (i==_map.end()) {
        Property& p = (_map[key] = Property(type,Value()));
        return ValueRef(p.first, p.second);
    } else if (i->second.first != type) {
        if (i->second.first==StringType) {
            free(i->second.second.s);
            i->second.second.s = NULL;
        }
    }
    return ValueRef(i->second.first, i->second.second);
}

