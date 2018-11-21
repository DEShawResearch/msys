#include "value.hxx"
#include <stdexcept>
#include <stdlib.h>
#include <string.h>

using namespace desres::msys;

Int ValueRef::asInt() const {
    if (_type==IntType) return _val.i;
    if (_type==FloatType) return (Int)_val.f;
    MSYS_FAIL("Cannot convert ValueRef string '" << _val.s << "' to int");
}

Float ValueRef::asFloat() const {
    if (_type==IntType) return _val.i;
    if (_type==FloatType) return _val.f;
    throw std::runtime_error("Cannot convert ValueRef to float");
}
String ValueRef::asString() const {
    if (_type==StringType) return _val.s ? _val.s : "";
    throw std::runtime_error("Cannot convert ValueRef to string");
}
const char* ValueRef::c_str() const {
    if (_type==StringType) return _val.s ? _val.s : "";
    throw std::runtime_error("Cannot convert ValueRef to string");
}

void ValueRef::fromInt(const Int& i) {
    if (_type==IntType) _val.i = i;
    else if (_type==FloatType) _val.f = i;
    else throw std::runtime_error("cannot assign int to string prop");
    if (_cb) _cb->valueChanged();
}

void ValueRef::fromFloat(const Float& i) {
    if (_type==IntType) _val.i = (Int)i;
    else if (_type==FloatType) _val.f = i;
    else throw std::runtime_error("cannot assign float to string prop");
    if (_cb) _cb->valueChanged();
}

void ValueRef::fromString(const String& i) {
    if (_type==IntType) 
        throw std::runtime_error("cannot assign string to int prop");
    else if (_type==FloatType) 
        throw std::runtime_error("cannot assign string to float prop");
    else {
        if (_val.s) free(_val.s);
        _val.s = strdup(i.c_str());
    }
    if (_cb) _cb->valueChanged();
}

bool ValueRef::operator==(const ValueRef& rhs) const {
    switch (_type) {
        case IntType:
            if (rhs._type==IntType) return _val.i==rhs._val.i;
            if (rhs._type==FloatType) return _val.i==rhs._val.f;
            break;
        case FloatType: 
            if (rhs._type==IntType) return _val.f==rhs._val.i;
            if (rhs._type==FloatType) return _val.f==rhs._val.f;
            break;
        default:
        case StringType:
            if (rhs._type==StringType) 
                return strcmp(_val.s ?     _val.s : "", 
                          rhs._val.s ? rhs._val.s : "")==0;
            break;
    }
    return false;
}

int ValueRef::compare(const ValueRef& rhs) const {
    switch (_type) {
        case IntType:
            if (rhs._type==IntType) 
                return _val.i<rhs._val.i ? -1 : 
                       _val.i>rhs._val.i ?  1 : 
                                            0 ;

            if (rhs._type==FloatType) 
                return _val.i<rhs._val.f ? -1 :
                       _val.i>rhs._val.f ?  1 : 
                                            0 ;
            break;
        case FloatType:
            if (rhs._type==FloatType) 
                return _val.f<rhs._val.f ? -1 :
                       _val.f>rhs._val.f ?  1 :
                                            0;

            if (rhs._type==IntType) 
                return _val.f<rhs._val.i ? -1 :
                       _val.f<rhs._val.i ?  1 :
                                            0 ;
            break;
        default:
        case StringType:
            if (rhs._type==StringType) 
                return strcmp(_val.s ?     _val.s : "", 
                          rhs._val.s ? rhs._val.s : "");
            break;
    }
    return 0;
}
