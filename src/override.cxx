#include "override.hxx"

using namespace desres::msys;

OverrideTable::OverrideTable(ParamTablePtr target)
: _target(target), _params(ParamTable::create()) {
}

OverrideTablePtr OverrideTable::create(ParamTablePtr target) {
    return OverrideTablePtr(new OverrideTable(target));
}

void OverrideTable::clear() {
    for (OverrideMap::const_iterator it=_map.begin(); it!=_map.end(); ++it) {
        _target->decref(it->first.first);
        _target->decref(it->first.second);
        _params->decref(it->second);
    }
    _map.clear();
}


Id OverrideTable::get(IdPair params) const {
    if (params.first>params.second) std::swap(params.first, params.second);
    OverrideMap::const_iterator it=_map.find(params);
    if (it!=_map.end()) return it->second;
    return BadId;
}

void OverrideTable::del(IdPair params) {
    if (params.first>params.second) std::swap(params.first, params.second);
    OverrideMap::iterator it=_map.find(params);
    if (it!=_map.end()) {
        _target->decref(it->first.first);
        _target->decref(it->first.second);
        _params->decref(it->second);
        _map.erase(it);
    }
}

void OverrideTable::set(IdPair params, Id param) {
    if (!_target->hasParam(params.first) ||
        !_target->hasParam(params.second)) {
        MSYS_FAIL("Invalid params: " << params.first << ", " << params.second);
    }
    if (!_params->hasParam(param)) {
        MSYS_FAIL("Invalid override param " << param);
    }
    if (params.first>params.second) std::swap(params.first, params.second);
    std::pair<OverrideMap::iterator,bool> ret;
    ret = _map.insert(std::make_pair(params, param));
    OverrideMap::iterator it = ret.first;
    if (ret.second) {
        _target->incref(it->first.first);
        _target->incref(it->first.second);
    } else {
        _params->decref(it->second);
    }
    _params->incref(param);
    it->second = param;
}

Id OverrideTable::count() const {
    return _map.size();
}

std::vector<IdPair> OverrideTable::list() const {
    std::vector<IdPair> p;
    for (OverrideMap::const_iterator it=_map.begin(); it!=_map.end(); ++it) {
        p.push_back(it->first);
    }
    return p;
}

