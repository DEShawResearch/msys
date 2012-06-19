/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"
#include <sstream>
#include <cstdio>

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::s_empty_attribute;

desres::msys::DestroArray desres::msys::Destro::Attribute::s_empty_array;

desres::msys::Destro::Attribute::Attribute()
  : m_block(NULL),m_attr(" e m p t y ")
{
}

desres::msys::Destro::Attribute::Attribute(Destro* block, const std::string& attr)
  : m_block(block), m_attr(attr)
{
  if (m_block == NULL) throw desres::msys::dessert("invalid block for attribute");
  if (!m_block->contains(attr)) throw desres::msys::dessert("Attribute error: "+attr);
}

desres::msys::Destro::Attribute::~Attribute() {
}

desres::msys::Destro* desres::msys::Destro::Attribute::owner() {
  if (!m_block) throw desres::msys::dessert("Immutable attribute cannot be updated");
  return m_block;
}

const desres::msys::Destro* desres::msys::Destro::Attribute::owner() const {
  if (!m_block) throw desres::msys::dessert("Empty attribute has no value");
  return m_block;
}

void desres::msys::Destro::Attribute::assign(const std::string& value) {
  compatible_schema('s');
  update(value);
}

void desres::msys::Destro::Attribute::assign(const char* value) {
  assign(std::string(value));
}

template<class T>
void desres::msys::Destro::Attribute::assigni(T value) {
  switch(type()) {
  case 'e':
    if (m_block == NULL) throw desres::msys::dessert("immutable");
    m_block->add_schema('i',m_attr); /*GCOV-IGNORE*/
    /* fall through */
  case 'i': {
    std::ostringstream ss; ss << value; update(ss.str());
    break;
  }
  case 'r':
    assign((double)value);
    break;
  default:
    throw desres::msys::dessert(m_attr+": LHS is not an integer");
  }
}

//template<> void desres::msys::Destro::Attribute::assigni<int>(int value);

void desres::msys::Destro::Attribute::assign(int value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(short value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(long value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(long long value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(unsigned int value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(unsigned short value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(unsigned long value) {
  assigni(value);
}

void desres::msys::Destro::Attribute::assign(bool value) {
  compatible_schema('b');

  std::ostringstream ss;
  ss << (int)value;
  update(ss.str());
}

void desres::msys::Destro::Attribute::compatible_schema(char desired) {
  char atype = type();

  if (atype == desired) return;
  if (atype == 'e') {
    if (m_block == NULL) throw desres::msys::dessert("immutable");/*GCOV-IGNORE*/
    m_block->add_schema(desired,m_attr);/*GCOV-IGNORE*/
    return;/*GCOV-IGNORE*/
  }

  throw desres::msys::dessert(m_attr+" is type " + atype + " and does not conform to type "+desired);
}

void desres::msys::Destro::Attribute::assign(double value) {
  compatible_schema('r');

  int precision = m_block->get_precision(m_attr);

  //char string[32+s_max_precision];
  char string[1024];
  if (precision < 0) precision = desres::msys::Destro::s_default_double_precision;
  ::sprintf(string,"%.*g",precision,value);

  update(string);
}

void desres::msys::Destro::Attribute::assign(float value) {
  compatible_schema('r');

  int precision = m_block->get_precision(m_attr);

  //char string[32+s_max_precision];
  char string[1024];
  if (precision < 0) precision = desres::msys::Destro::s_default_float_precision;
  ::sprintf(string,"%.*g",precision,value);

  update(string);
}

desres::msys::Destro::Attribute&
desres::msys::Destro::Attribute::operator=(const char* rhs) {
  assign(rhs);
  return *this;
}

char desres::msys::Destro::Attribute::type() const {
  if (m_block == NULL) return 'e';
  if (m_block->has_attr(m_attr)) return m_block->get_type(m_attr);
  if (m_block->has_block(m_attr)) return 'x';
  return 'e';/*GCOV-IGNORE*/
}

std::string desres::msys::Destro::Attribute::doc() const {
  if (m_block == NULL) throw desres::msys::dessert("invalid block for attribute");
  return m_block->get_doc(m_attr);
}

void desres::msys::Destro::Attribute::doc(const std::string& doc) {
  if (m_block == NULL) throw desres::msys::dessert("invalid block for attribute");
  return m_block->set_doc(m_attr,doc);
}


std::string desres::msys::Destro::Attribute::value() const {
  if (m_block == NULL) throw desres::msys::dessert("invalid block for attribute");
  return m_block->get_value(m_attr);
}

short
desres::msys::Destro::Attribute::precision() const {
  return owner()->get_precision(m_attr);
}

void
desres::msys::Destro::Attribute::precision(short precision) {
  owner()->set_precision(m_attr,precision);
}

bool
desres::msys::Destro::Attribute::is_empty() const {
  if (m_block == NULL) return true;
  if (type() == 'x') return false;
  if (!m_block->has_value(m_attr)) return true;
  std::string val = m_block->get_value(m_attr);
 
  return val == "<>";
}

void desres::msys::Destro::Attribute::update(const std::string& value) {
  owner()->set_unsafe(m_attr,value);
}

desres::msys::Destro::Attribute::operator std::string() const {
  if (!m_block) return "<>";
  if (type() != 's') throw desres::msys::dessert(m_attr+" not a string");
  return m_block->get_value(m_attr);
}

desres::msys::Destro::Attribute::operator double() const {
  if (type() != 'r') throw desres::msys::dessert(m_attr+" not a float");
  std::string val = owner()->get_value(m_attr);
  const char* bytes = val.c_str();
  char* end_byte = NULL;
  double result = strtod(bytes,&end_byte);
  if (end_byte && *end_byte != '\0') throw desres::msys::dessert(m_attr+" not a float");
  return result;
}

desres::msys::Destro::Attribute::operator float() const {
  return operator double();
}

desres::msys::Destro::Attribute::operator long() const {
  if (!m_block) throw desres::msys::dessert("converting empty");
  if (type() != 'i') throw desres::msys::dessert(m_attr+" not an integer");
  std::string val = owner()->get_value(m_attr);
  const char* bytes = val.c_str();
  char* end_byte = NULL;
  long result = strtol(bytes,&end_byte,10);
  if (end_byte && *end_byte != '\0') throw desres::msys::dessert(m_attr+" not an integer");
  return result;
}

desres::msys::Destro::Attribute::operator int() const {
  return operator long();
}

desres::msys::Destro::Attribute::operator short() const {
  return operator long();
}

desres::msys::Destro::Attribute::operator unsigned long() const {
  return operator long();
}

desres::msys::Destro::Attribute::operator unsigned int() const {
  return operator long();
}

desres::msys::Destro::Attribute::operator unsigned short() const {
  return operator long();
}

desres::msys::Destro::Attribute::operator bool() const {
  if (!m_block) throw desres::msys::dessert("converting empty attribute");
  if (type() != 'b') throw desres::msys::dessert(m_attr+" not a bool");
  std::string val = owner()->get_value(m_attr);
  if (val == "0") return false;
  if (val == "1") return true;
  if (val == "<>") throw desres::msys::dessert("converting empty value");
  throw desres::msys::dessert(m_attr+" not a bool");
}

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator[](char const* attr) {
  return this->operator[](std::string(attr));
}

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator[](char const* attr) const {
  return this->operator[](std::string(attr));
}

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator[](const std::string& attr) {
  if (!m_block) throw desres::msys::dessert("Attribute error: " + attr);

  // If we refer to a block, forward the request...
  if (m_block->has_block(m_attr)) {
    return m_block->block(m_attr)->operator[](attr);
  }

  // If we don't refer to a block, we can't pull a subsequent attr
  throw desres::msys::dessert(m_attr+" is a leaf attribute and cannot contain " + attr);
}

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator[](const std::string& attr) const {
  if (!m_block) throw desres::msys::dessert("Attribute error: " + attr);

  // If we refer to a block, forward the request...
  if (m_block->has_block(m_attr)) {
    return m_block->block(m_attr)->operator[](attr);
  }

  // If we don't refer to a block, we can't pull a subsequent attr
  throw desres::msys::dessert(m_attr+" is a leaf attribute and cannot contain " + attr);
}

desres::msys::Destro&
desres::msys::Destro::Attribute::operator[](size_t block) {
  if (!m_block) throw desres::msys::dessert("subblock reference error");

  // If we refer to a block, forward the request...
  if (m_block->has_block(m_attr)) {
    return m_block->block(m_attr)->block(block);
  }

  // If we don't refer to a block, we can't pull a subsequent attr
  throw desres::msys::dessert(m_attr+" is a leaf attribute and does not contain subblocks");
}

const desres::msys::Destro&
desres::msys::Destro::Attribute::operator[](size_t block) const {
  if (!m_block) throw desres::msys::dessert("subblock reference error");

  // If we refer to a block, forward the request...
  if (m_block->has_block(m_attr)) {
    return m_block->block(m_attr)->block(block);/*GCOV-IGNORE*/
  }

  // If we don't refer to a block, we can't pull a subsequent attr
  throw desres::msys::dessert(m_attr+" is a leaf attribute and does not contain subblocks");
}

#if 0
desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator()(const std::string& attr) {

  // Empty attributes just return a new empty
  if (!m_block) return *this;

  Destro* block = as_block();

  // If this attribute is a block, then we build a new sub-attribute
  if (block) {
    return Attribute(block,attr);
  }
  
  // Otherwise, just return the empty
  return s_empty_attribute;
}
#endif

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator()(const std::string& attr) const {
  // Empty attributes just return a new empty
  if (!m_block) return *this;

  const Destro* block = as_block();

  // If this attribute is a block, then we build a new sub-attribute
  if (block) {
    return Attribute(const_cast<Destro*>(block),attr);
  }
  
  // Otherwise, just return the empty
  return s_empty_attribute;/*GCOV-IGNORE*/
}

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator()(char const* attr) {
  return this->operator()(std::string(attr));
}

desres::msys::Destro::Attribute
desres::msys::Destro::Attribute::operator()(const char* attr) const {
  return this->operator()(std::string(attr));
}

desres::msys::Destro*
desres::msys::Destro::Attribute::as_block() {
  if (!m_block) return &s_empty_array;

  if (m_block->has_block(m_attr)) {
    return &(m_block->block(m_attr));
  }

  return NULL;
}

const desres::msys::Destro*
desres::msys::Destro::Attribute::as_block() const {
  if (!m_block) return &s_empty_array;

  if (m_block->has_block(m_attr)) {
    return &(m_block->block(m_attr));
  }

  return NULL;
}

desres::msys::Destro*
desres::msys::Destro::Attribute::operator->() {
  Destro* block = as_block();
  if (!block) throw desres::msys::dessert("attribute " + m_attr + " is not a block");
  return block;
}


const desres::msys::Destro*
desres::msys::Destro::Attribute::operator->() const {
  const Destro* block = as_block();
  if (!block) throw desres::msys::dessert("attribute " + m_attr + " is not a block");
  return block;
}


desres::msys::Destro::Attribute::operator desres::msys::Destro&() const {
  const Destro* block = as_block();
  return *(const_cast<Destro*>(block));
}

desres::msys::Destro::Attribute::operator desres::msys::DestroArray&() const {
  const Destro* block = as_block();
  if (!m_block) return s_empty_array;

  const DestroArray* array = dynamic_cast<const DestroArray*>(block);
  if (array == NULL) throw desres::msys::dessert(m_attr+" not an array");
  return *const_cast<DestroArray*>(array);
}

std::ostream& desres::msys::operator<<(std::ostream& os,const Destro::Attribute& A) {
  const Destro* block = A.as_block();
  if (block) return (os << *block);

  if (&A == &Destro::Attribute::s_empty_attribute) return os << "<>";

  
  return os << A.value();
}
