/* @COPYRIGHT@ */

#include "json.hxx"
#include "JSON_parser.h"

#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <cassert>

using namespace desres::msys::fastjson;

static const Json invalid_json;

void Json::clear() {
    switch (_kind) {
        case Array: 
            for (int i=0; i<u.a.num; i++) {
                u.a.elems[i].clear();
            }
            free(u.a.elems);
            break;
        case Object: 
            for (int i=0; i<u.o.num; i++) {
                free(u.o.elems[i].key);
                u.o.elems[i].val.clear();
            }
            free(u.o.elems);
            break;
        case String: 
            free(u.s); 
            break;
        default: ;
    }
    _kind = Invalid;
}

Json& Json::to_array() {
    clear();
    _kind = Array;
    u.a.num=0;
    u.a.max=4;
    u.a.elems = (Json *)calloc(u.a.max, sizeof(Json));
    return *this;
}

Json& Json::to_object() {
    clear();
    _kind = Object;
    u.o.num=0;
    u.o.max=4;
    u.o.elems = (keyval_t *)calloc(u.o.max, sizeof(keyval_t));
    return *this;
}

Json& Json::append(Json& js) {
    if (_kind != Array) {
        throw std::invalid_argument("not an array");
    }
    if (u.a.num==u.a.max) {
        u.a.max *= 2;
        u.a.elems = (Json *)realloc( u.a.elems, u.a.max * sizeof(Json));
        memset(u.a.elems+u.a.num, 0, (u.a.max-u.a.num)*sizeof(Json));
    }
    u.a.elems[u.a.num++].swap(js);
    return *this;
}

Json& Json::append(const char * key, Json& val) {
    if (_kind != Object) {
        throw std::invalid_argument("not an array");
    }
    if (u.o.num==u.a.max) {
        u.o.max *= 2;
        u.o.elems = (keyval_t *)realloc(
                u.o.elems, u.o.max * sizeof(*u.o.elems)); 
        memset(u.o.elems+u.o.num, 0, (u.o.max-u.o.num)*sizeof(*u.o.elems));
    }
    u.o.elems[u.o.num].key = strdup(key);
    u.o.elems[u.o.num++].val.swap(val);
    return *this;
}

const char * Json::kindstr() const {
#define MYCASE(x) case x : return #x 
    switch (_kind) {
        MYCASE(Invalid);
        MYCASE(Array);
        MYCASE(Object);
        MYCASE(Int);
        MYCASE(Float);
        MYCASE(String);
        MYCASE(Null);
        MYCASE(Bool);
    }
#undef MYCASE
    assert(false);
    return NULL;
}

int Json::size() const {
    if (kind()==Array) return u.a.num;
    if (kind()==Object) return u.o.num;
    throw std::invalid_argument("not array or object");
}

const Json& Json::elem(int i) const {
    if (i<0) throw std::invalid_argument("negative index");
    if (kind()==Array) {
        if (i<u.a.num) return u.a.elems[i];
        throw std::invalid_argument("index out of bounds");
    }
    if (kind()==Object) {
        if (i<u.o.num) return u.o.elems[i].val;
        throw std::invalid_argument("index out of bounds");
    }
    if (kind()!=Invalid) {
        throw std::invalid_argument("not array or object");
    }
    return invalid_json;
}

Json& Json::elem(int i) {
    if (i<0) throw std::invalid_argument("negative index");
    if (kind()==Array) {
        if (i<u.a.num) return u.a.elems[i];
        throw std::invalid_argument("index out of bounds");
    }
    if (kind()==Object) {
        if (i<u.o.num) return u.o.elems[i].val;
        throw std::invalid_argument("index out of bounds");
    }
    throw std::invalid_argument("not array or object");
}

const char * Json::key(int i) const {
    if (i<0) throw std::invalid_argument("negative index");
    if (kind()==Object) {
        if (i<u.o.num) return u.o.elems[i].key;
        throw std::invalid_argument("index out of bounds");
    }
    throw std::invalid_argument("not an object");
}

const Json& Json::get(const char * key) const {
    if (kind()==Object) {
        for (int i=0; i<u.o.num; i++) {
            if (!strcmp(key, u.o.elems[i].key)) {
                return u.o.elems[i].val;
            }
        }
        return invalid_json;
    }
    if (kind()!=Invalid) {
        throw std::invalid_argument("not an object");
    }
    return invalid_json;
}

void Json::swap(Json& js) {
    Json tmp;
    memcpy(&tmp.u, &u, sizeof(u));
    memcpy(&u, &js.u, sizeof(u));
    memcpy(&js.u, &tmp.u, sizeof(u));
    kind_t t = _kind;
    _kind = js._kind;
    js._kind = t;
}

Json& Json::copy(const Json& js) {
    switch (js.kind()) {
        case Int:
            to_int(js.as_int());
            break;
        case Float:
            to_float(js.as_float());
            break;
        case String:
            to_string(js.as_string());
            break;
        case Bool:
            to_bool(js.as_bool());
            break;
        case Null:
            to_null();
            break;
        case Array:
            {
                to_array();
                int i,n=js.size();
                for (i=0; i<n; i++) {
                    Json tmp;
                    tmp.copy(js.elem(i));
                    append(tmp);
                }
            }
            break;
        case Object:
            {
                to_object();
                int i,n=js.size();
                for (i=0; i<n; i++) {
                    Json tmp;
                    tmp.copy(js.elem(i));
                    append(js.key(i), tmp);
                }
            }
            break;
        default:
            throw std::invalid_argument("invalid json");
    }
    return *this;
}

bool Json::equals(const Json& js) const {
    if (kind()!=js.kind()) return false;

    switch (js.kind()) {
        case Int:   
            return as_int()==js.as_int();
        case Float: 
            return as_float()==js.as_float();
        case String:
            return !strcmp(as_string(), js.as_string());
        case Bool:  
            return as_bool()==js.as_bool();
        case Null:  
            return true;
        case Array:
            {
                int i,n=js.size();
                if (size()!=n) return false;
                for (i=0; i<n; i++) {
                    if (!elem(i).equals(js.elem(i))) return false;
                }
            }
            return true;
        case Object:
            {
                int i,n=js.size();
                if (size()!=n) return false;
                for (i=0; i<n; i++) {
                    const char * k = key(i);
                    const Json& lhs = elem(i);
                    const Json& rhs = js.get(k);
                    /* if rhs is invalid, the kind() test will fail */
                    if (!lhs.equals(rhs)) return false;
                }
            }
            return true;
            break;
        default:
            throw std::invalid_argument("invalid json");
    }
    /* never gets here */
    return false;
}

const char * Json::as_string() const {
    if (kind()==String) return u.s;
    throw std::invalid_argument("not a string");
}

bool Json::as_bool() const {
    if (kind()==Int || kind()==Bool) return u.i;
    throw std::invalid_argument("not bool");
}

int Json::as_int() const {
    if (kind()==Int || kind()==Bool) return u.i;
    throw std::invalid_argument("not an int");
}

double Json::as_float() const {
    if (kind()==Int) return u.i;
    if (kind()==Float) return u.f;
    throw std::invalid_argument("not a double");
}
