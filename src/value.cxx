#include "value.hxx"
#include <stdexcept>
#include <stdlib.h>
#include <string.h>

using namespace desres::msys;

Int ValueRef::asInt() const {
    if (_type==IntType) return _val.i;
    if (_type==FloatType) return (Int)_val.f;
    throw std::runtime_error("Cannot convert Params value to int");
}

Float ValueRef::asFloat() const {
    if (_type==IntType) return _val.i;
    if (_type==FloatType) return _val.f;
    throw std::runtime_error("Cannot convert Params value to float");
}
String ValueRef::asString() const {
    if (_type==StringType) return _val.s ? _val.s : "";
    throw std::runtime_error("Cannot convert Params value to string");
}

void ValueRef::fromInt(const Int& i) {
    if (_type==IntType) _val.i = i;
    else if (_type==FloatType) _val.f = i;
    else throw std::runtime_error("cannot assign int to string prop");
}

void ValueRef::fromFloat(const Float& i) {
    if (_type==IntType) _val.i = (Int)i;
    else if (_type==FloatType) _val.f = i;
    else throw std::runtime_error("cannot assign float to string prop");
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
}
