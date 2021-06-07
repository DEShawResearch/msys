/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include <cstdio>
#include <msys/fastjson/print.hxx>

using desres::msys::fastjson::floatify;

// -----------------------------------------------
//  P R O T E C T E D
// -----------------------------------------------

desres::msys::DestroArray::key_value_t*
desres::msys::DestroArray::find_schema(const std::string& attr,const desres::msys::ZingPool& pool) {
  for(std::vector<key_value_t>::iterator p = m_data.begin(),
      en = m_data.end(); p != en; ++p) {
    if ((*p).key.string(pool) == attr) {
      return &(*p);
    }
  }
  return NULL;
}

const desres::msys::DestroArray::key_value_t*
desres::msys::DestroArray::find_schema(const std::string& attr,const desres::msys::ZingPool& pool) const {
  for(std::vector<key_value_t>::const_iterator p = m_data.begin(),
      en = m_data.end(); p != en; ++p) {
    if ((*p).key.string(pool) == attr) return &(*p);
  }
  return NULL;
}



ssize_t desres::msys::DestroArray::integer(const std::string& string) {
  size_t length = string.size();
  if (length == 0) {
    MSYS_FAIL("string not convertable to integer");
  }

  // Optional sign
  ssize_t sign = 1;
  std::string::const_iterator p = string.begin();
  std::string::const_iterator en = string.end();
  if (*p == '-') {
    sign = -1;
    ++p;
  }

  ssize_t total = 0;
  for(; p != en; ++p) {
    char c = *p;
    if (!isdigit(c)) {
      MSYS_FAIL("string not convertable to integer");
    }
    total = total*10 + (c-'0');
  }

  return sign*total;
}

// -----------------------------------------------
//  P U B L I C
// -----------------------------------------------
desres::msys::DestroArray::DestroArray(Destro* parent, size_t num_elements)
  : DestroNamedBlock(parent)
{
  for(size_t i=0; i<num_elements; ++i) {
    m_rows.push_back( DestroRow(this,i) );
  }
}

desres::msys::DestroArray::~DestroArray() {
}

size_t desres::msys::DestroArray::size() const {
  return m_rows.size();
}


desres::msys::Destro& desres::msys::DestroArray::bulk_row(const std::string& offset, const std::vector<std::string>& elements) {
  ZingPool& zpool = mutable_pool(); // check that array is mutable

  ssize_t rowid;
  try {
    rowid = integer(offset);
  } catch (...) {
    MSYS_FAIL("invalid row number " + offset);
  }

  size_t width = m_data.size();
  if ((size_t)rowid != size()+1) {
    MSYS_FAIL("out of order row " + offset);
  }

  if (elements.size() != width) {
    MSYS_FAIL("invalid number of elements in bulk row");
  }

  desres::msys::Destro& row = new_block(offset);

  for(unsigned i=0;i<width;++i) {
    key_value_t& location = m_data[i];

    const std::string& value = elements[i];

    // If we have <> or an empty string, then we don't set a value
    if (value == "<>" || value == "") continue;

    // Strip quotes if present
    Zing zval;
    if (value[0] == '"' && value[value.size()-1]) {
      zval = Zing(value.substr(1,value.size()-2),zpool);
    } else {
      zval = Zing(value,zpool);
    }

    location.values[rowid-1] = zval;
  }

  return row;
}

void desres::msys::DestroArray::resize(size_t n) {
  mutable_pool(); // check that array is mutable

  size_t sz = size();

  if (n < sz) {
    for(std::vector<key_value_t>::iterator p=m_data.begin(),
          en=m_data.end(); p != en; ++p) {
      std::vector<Zing>& values = (*p).values;
      values.erase(values.begin()+n,values.end());
    }
    m_rows.erase(m_rows.begin()+n,m_rows.end());
  } else {
    Zing empty;
    for(std::vector<key_value_t>::iterator p=m_data.begin(),
          en=m_data.end(); p != en; ++p) {
      std::vector<Zing>& values = (*p).values;
      values.insert(values.begin()+sz,n-sz,empty);
    }
    
    for(size_t i=sz; i<n; ++i) {
      m_rows.push_back( DestroRow(this,i) );
    }
  }
}



void desres::msys::DestroArray::add_schema(char type,const std::string& attr,const std::string& doc) {
  ZingPool& zpool = mutable_pool();

  // Perhaps we just need to update the one we have
  key_value_t* location = find_schema(attr,zpool);
  if (location) {
    if (location->type != type) MSYS_FAIL("Invalid override of schema for " + attr);
    if (doc.size() > 0) location->doc = Zing(doc,zpool);
    return;
  }
 
  // We don't have that key, so we add a new one
  m_data.push_back(key_value_t());
  key_value_t& keyvalue = m_data[m_data.size()-1];

  keyvalue.key = Zing(attr,zpool);
  Zing empty;
  size_t n = size();
  keyvalue.values.clear();
  keyvalue.values.reserve(n);
  for(size_t i=0;i<n;++i) {
    keyvalue.values.push_back(empty);
  }
  keyvalue.doc = Zing("",zpool);
  keyvalue.type = (unsigned)type;
  keyvalue.precision = -1;
}

std::map<std::string,desres::msys::Destro::schema_t> desres::msys::DestroArray::schemas() const {
  const ZingPool& zpool = pool();
  std::map<std::string,schema_t> result;
  for(std::vector<key_value_t>::const_iterator p = m_data.begin(),
        en = m_data.end(); p != en; ++p) {
    const key_value_t& description = *p;
    std::string attr = description.key.string(zpool);
    schema_t& schema = result[attr];
    schema.type = description.type;
    schema.attr = attr;
    schema.doc  = description.doc.string(zpool);
  }
  return result;/*GCOV-IGNORE*/
}

std::vector<desres::msys::Destro::schema_t> desres::msys::DestroArray::ordered_schema() const {
  const ZingPool& zpool = pool();
  std::vector<desres::msys::Destro::schema_t> result;
  for(std::vector<key_value_t>::const_iterator p = m_data.begin(),
        en = m_data.end(); p != en; ++p) {
    const key_value_t& description = *p;
    schema_t schema;
    schema.type = description.type;
    schema.attr = description.key.string(zpool);
    schema.doc  = description.doc.string(zpool);

    result.push_back(schema);
  }
  return result;/*GCOV-IGNORE*/
}

std::string desres::msys::DestroArray::get_value(const std::string&) const {
  MSYS_FAIL("get_value not meaningful for arrays");
}

char desres::msys::DestroArray::get_type(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    MSYS_FAIL("Attribute Error: "+attr); /*GCOV-IGNORE*/
  }
  return location->type;
}

int desres::msys::DestroArray::get_precision(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    MSYS_FAIL("Attribute Error: "+attr);
  }

  return location->precision;
}

void desres::msys::DestroArray::set_precision(const std::string& attr,int precision) {
  mutable_pool();

  if (precision > s_max_precision) precision = s_max_precision;
  const ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    MSYS_FAIL("Attribute Error: "+attr);
  }
  location->precision = precision;
}

std::string desres::msys::DestroArray::get_doc(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    MSYS_FAIL("Attribute Error: "+attr);
  }
  return location->doc.string(zpool);
}

void desres::msys::DestroArray::set_doc(const std::string& attr, const std::string& doc) {
  mutable_pool();

  ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    MSYS_FAIL("Attribute Error: "+attr);
  }
  location->doc = Zing(doc,zpool);
}

desres::msys::Destro::Attribute desres::msys::DestroArray::get_attr(const std::string& attr) {
  return Attribute(this,attr);
}

const desres::msys::Destro::Attribute desres::msys::DestroArray::get_attr(const std::string& attr) const {
  return Attribute(const_cast<DestroArray*>(this),attr);
}

void desres::msys::DestroArray::set_unsafe(const std::string& attr,const std::string& value) {
  ZingPool& zpool = mutable_pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    MSYS_FAIL("Attribute Error: "+attr);
  }

  Zing zvalue(value,zpool);
  size_t n = size();
  location->values.clear();
  location->values.reserve(n);
  for(size_t i=0;i<n;++i) {
    location->values.push_back(zvalue);
  }
}

void desres::msys::DestroArray::set_attr(char type,const std::string& attr,const std::string& value) {
  ZingPool& zpool = mutable_pool();
  key_value_t* location = find_schema(attr,zpool);

  // TODO: This can probably be improved to lookup and fetch
  // at the same time
  if (location == NULL) {
    add_schema(type,attr);
    location = find_schema(attr,zpool);
    if (location == NULL) {
      MSYS_FAIL("Attribute Error: " << attr);
    }
  }
  location->type = type;

  Zing zvalue = zingify(value,zpool);
  size_t n = size();
  for(size_t i=0;i<n;++i) {
    location->values.at(i) = zvalue;
  }

}

void desres::msys::DestroArray::clear(const std::string& attr) {
  ZingPool& zpool = mutable_pool();
  key_value_t* location = find_schema(attr,zpool);

  if (location == NULL) return;

  Zing zvalue;
  size_t n = size();
  for(size_t i=0;i<n;++i) {
    location->values.at(i) = zvalue;
  }
}

void desres::msys::DestroArray::del(const std::string& attr) {
  // Perhaps it is an attribute....
  ZingPool& zpool = mutable_pool();
  for(std::vector<key_value_t>::iterator p = m_data.begin(),
      en = m_data.end(); p != en; ++p) {
    if ((*p).key.string(zpool) == attr) {
      m_data.erase(p);
      return;
    }
  }

  // Otherwise, silently ignore
}

void desres::msys::DestroArray::del(size_t blockno) {
  if (blockno >= 1 && blockno <= size()) {

    // Erase row's entry in each column data
    for(std::vector<key_value_t>::iterator p = m_data.begin(),
          en = m_data.end(); p != en; ++p) {
      std::vector<Zing>& column = (*p).values;
      column.erase(column.begin()+(blockno-1));
    }

    // Delete the appropriate row object
    std::deque<DestroRow>::iterator p = m_rows.begin()+(blockno-1);
    p = m_rows.erase(p);

    // Renumber the remaining rows...
    for(;p != m_rows.end(); ++blockno,++p) {
      (*p).setrow(blockno);
    }
  }
}

desres::msys::Destro& desres::msys::DestroArray::new_block(const std::string& name) {
  mutable_pool();
  size_t n = size()+1;
  resize(n);
  return block(n);
}

desres::msys::DestroArray& desres::msys::DestroArray::new_array(const std::string& name,size_t num_elements) {
  MSYS_FAIL("new_array not meaningful for arrays");
}

desres::msys::Destro::Attribute desres::msys::DestroArray::operator[](const std::string& name) {
  return Attribute(this,name);
}

const desres::msys::Destro::Attribute desres::msys::DestroArray::operator[](const std::string& name) const {
  return Attribute(const_cast<DestroArray*>(this),name);
}

desres::msys::Destro& desres::msys::DestroArray::operator[](ssize_t i) {
  return block(i);
}

const desres::msys::Destro& desres::msys::DestroArray::operator[](ssize_t i) const {
  return block(i);
}

bool desres::msys::DestroArray::is_array() const {
  return true;
}

bool desres::msys::DestroArray::has_block(const std::string& blockname) const {
  ssize_t ival;
  try {
    ival = integer(blockname);
  } catch(...) {
    return false;
  }
  if (ival < 0) return false;
  size_t uval = static_cast<size_t>(ival);
  return (uval >= 1 && uval <= size());
}

bool desres::msys::DestroArray::has_attr(const std::string& attr) const {
  const key_value_t* location = find_schema(attr,pool());
  return (location != NULL);
}

bool desres::msys::DestroArray::has_value(const std::string&) const {
  MSYS_FAIL("has_value not meaningful for arrays");
}

desres::msys::Destro& desres::msys::DestroArray::block(size_t i) {
  return m_rows.at(i-1);
}

const desres::msys::Destro& desres::msys::DestroArray::block(size_t i) const {
  return m_rows.at(i-1);
}

desres::msys::Destro& desres::msys::DestroArray::block(const std::string& name) {
  ssize_t i = integer(name);
  return block(i);
}

const desres::msys::Destro& desres::msys::DestroArray::block(const std::string& name) const {
  ssize_t i = integer(name);
  return block(i);
}

void desres::msys::DestroArray::write(std::ostream& os, int level) const {
  std::string aname = name();
  //if (aname == "") MSYS_FAIL("Cannot write unnamed array");
  if (aname == "") aname = "anonymous";

  size_t rows = size();
  if (rows == 0) return; // do not write out empty DestroArray

  int newlevel = 0;
  std::string indentation = indent(level,newlevel);
  int dummy = 0;
  std::string indentation2 = indent(newlevel,dummy);
  const char* lineterm = (level >= 0)?"\n":"";

  os << indentation << aname << '[' << rows << "] {" << lineterm;

  const ZingPool& zpool = pool();

  // Write the attributes...
  for(std::vector<key_value_t>::const_iterator p=m_data.begin(),
        en = m_data.end(); p != en; ++p) {
    const key_value_t& schema = *p;
    os << indentation2 << schema_out(schema.type,schema.key,schema.doc,zpool) << lineterm;
    // width = max(zingvalue.size());
  }
  os << indentation2 << ":::" << lineterm;

  // Write the rows...
  size_t n = size();
  for(size_t i=0; i<n; ++i) {
    os << indentation2 << i+1;
    for(std::vector<key_value_t>::const_iterator p=m_data.begin(),
          en = m_data.end(); p != en; ++p) {
      const key_value_t& schema = *p;
      os << ' ' << quotify(schema.values[i],zpool);
    }
    os << lineterm;
  }

  os << indentation2 << ":::" << lineterm;

  os << indentation << '}' << lineterm;

}

// -----------------------------------------------
// Unpack string to value
// -----------------------------------------------
template<class T>
static T unpack(const std::string& x) {
  MSYS_FAIL("Not implemented");  /*GCOV-IGNORE*/
}

template<> long unpack<long>(const std::string& x) {
  if (x == "<>") MSYS_FAIL("unpacking empty value");
  return strtol(x.c_str(),NULL,10);
}
template<> bool unpack<bool>(const std::string& x) {return unpack<long>(x);}
template<> int unpack<int>(const std::string& x) {return unpack<long>(x);}
template<> short unpack<short>(const std::string& x) {return unpack<long>(x);}
template<> unsigned int unpack<unsigned int>(const std::string& x) {return unpack<long>(x);}
template<> unsigned short unpack<unsigned short>(const std::string& x) {return unpack<long>(x);}
template<> unsigned long unpack<unsigned long>(const std::string& x) {return unpack<long>(x);}
template<> double unpack<double>(const std::string& x) {
  if (x == "<>") MSYS_FAIL("unpack empty token");
  return strtod(x.c_str(),NULL);
}
template<> float unpack<float>(const std::string& x) {return unpack<double>(x);}
template<>
std::string unpack<std::string>(const std::string& x) {
  return x;
}

// -----------------------------------------------
// Unpack string to value (with default)
// -----------------------------------------------
template<class T>
T unpack(const std::string& x, T defval) {
  if (x == "<>") return defval;
  return unpack<T>(x);
}
template int unpack(const std::string&,int);
template short unpack(const std::string&,short);
template long unpack(const std::string&,long);
template unsigned int unpack(const std::string&,unsigned int);
template unsigned short unpack(const std::string&,unsigned short);
template unsigned long unpack(const std::string&,unsigned long);

template<class T>
void desres::msys::DestroArray::column(const std::string& attr, std::vector<T>& values) const {
  values.resize(0);

  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (!location) return;

  size_t n = size();
  values.resize(n);
  for(size_t i=0; i<n; ++i) {
    Zing z = location->values.at(i);
    if (z.is_empty()) MSYS_FAIL("Cannot unpack empty value");
    values[i] = unpack<T>(z.string(zpool));
  }

}

template void desres::msys::DestroArray::column<bool>(const std::string& attr, std::vector<bool>& values) const;
template void desres::msys::DestroArray::column<int>(const std::string& attr, std::vector<int>& values) const;
template void desres::msys::DestroArray::column<short>(const std::string& attr, std::vector<short>& values) const;
template void desres::msys::DestroArray::column<long>(const std::string& attr, std::vector<long>& values) const;
template void desres::msys::DestroArray::column<unsigned int>(const std::string& attr, std::vector<unsigned int>& values) const;
template void desres::msys::DestroArray::column<unsigned short>(const std::string& attr, std::vector<unsigned short>& values) const;
template void desres::msys::DestroArray::column<unsigned long>(const std::string& attr, std::vector<unsigned long>& values) const;
template void desres::msys::DestroArray::column<float>(const std::string& attr, std::vector<float>& values) const;
template void desres::msys::DestroArray::column<double>(const std::string& attr, std::vector<double>& values) const;
template void desres::msys::DestroArray::column<std::string>(const std::string& attr, std::vector<std::string>& values) const;

template<class T>
void desres::msys::DestroArray::column(const std::string& attr, std::vector<T>& values, const T& defval) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);

  size_t n = size();
  values.resize(n);
  for(size_t i=0; i<n; ++i) {
    if (location) {
      Zing z = location->values.at(i);
      values[i] = z.is_empty()?defval:unpack<T>(z.string(zpool),defval);
    } else {
      values[i] = defval;
    }
  }

}

template void desres::msys::DestroArray::column<bool>(const std::string& attr, std::vector<bool>& values,const bool&) const;
template void desres::msys::DestroArray::column<int>(const std::string& attr, std::vector<int>& values,const int&) const;
template void desres::msys::DestroArray::column<short>(const std::string& attr, std::vector<short>& values,const short&) const;
template void desres::msys::DestroArray::column<long>(const std::string& attr, std::vector<long>& values,const long&) const;
template void desres::msys::DestroArray::column<unsigned int>(const std::string& attr, std::vector<unsigned int>& values,const unsigned int&) const;
template void desres::msys::DestroArray::column<unsigned short>(const std::string& attr, std::vector<unsigned short>& values,const unsigned short&) const;
template void desres::msys::DestroArray::column<unsigned long>(const std::string& attr, std::vector<unsigned long>& values,const unsigned long&) const;
template void desres::msys::DestroArray::column<float>(const std::string& attr, std::vector<float>& values,const float&) const;
template void desres::msys::DestroArray::column<double>(const std::string& attr, std::vector<double>& values,const double&) const;
template void desres::msys::DestroArray::column<std::string>(const std::string& attr, std::vector<std::string>& values,const std::string&) const;

static std::string precision_string(const bool* val,char type, int precision) {
  if (type != 'b') MSYS_FAIL("cannot store a bool here");
  if (*val) return "1";
  return "0";
}
static std::string precision_string(const double* val,char type, int precision) {
  if (type != 'r') MSYS_FAIL("cannot store a double here");
  char buf[32];
  floatify(*val, buf);
  return buf;
}

template<class T>
static std::string iprecision_string(const T* val,char type, int precision) {
  switch(type) {
  case 'i': {
    std::ostringstream ss; ss << *val;
    return ss.str();
  }
  case 'r': {
    double xval = (double)*val;
    return precision_string(&xval,type,precision);
  }
  default:
    MSYS_FAIL("cannot store an int here");
  }
}
static std::string precision_string(const int* val,char type, int precision) {
  return iprecision_string(val,type,precision);
}
static std::string precision_string(const short* val,char type, int precision) {
  return iprecision_string(val,type,precision);
}
static std::string precision_string(const long* val,char type, int precision) {
  return iprecision_string(val,type,precision);
}
static std::string precision_string(const unsigned int* val,char type, int precision) {
  return iprecision_string(val,type,precision);
}
static std::string precision_string(const unsigned short* val,char type, int precision) {
  return iprecision_string(val,type,precision);
}
static std::string precision_string(const unsigned long* val,char type, int precision) {
  return iprecision_string(val,type,precision);
}
static std::string precision_string(const float* val,char type, int precision) {
  if (type != 'r') MSYS_FAIL("cannot store a float here");
  char buf[32];
  floatify(*val, buf);
  return buf;
}

template<class T>
void desres::msys::DestroArray::set_column(const std::string& attr, const T* begin, size_t n, size_t stride) {
  ZingPool& zpool = mutable_pool();
  key_value_t* location = find_schema(attr,zpool);
  if (!location) MSYS_FAIL("attibute error: "+attr);
  if (!begin) return; // Ignore NULL pointers

  // Must be the right size
  if (n != size()) MSYS_FAIL("size mismatch in column");

  std::vector<Zing>& column = location->values;
  for(size_t i=0; i<n; ++i) {
    std::string rval = precision_string(begin+i*stride,location->type,location->precision);
    Zing zz(rval,zpool);
    column[i] = zz;
  }
}

template void desres::msys::DestroArray::set_column<bool>(const std::string& attr, const bool* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<int>(const std::string& attr, const int* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<short>(const std::string& attr, const short* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<long>(const std::string& attr, const long* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<unsigned int>(const std::string& attr, const unsigned int* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<unsigned short>(const std::string& attr, const unsigned short* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<unsigned long>(const std::string& attr, const unsigned long* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<float>(const std::string& attr, const float* begin, size_t size, size_t stride);
template void desres::msys::DestroArray::set_column<double>(const std::string& attr, const double* begin, size_t size, size_t stride);

template<class T>
void desres::msys::DestroArray::set_column(const std::string& attr, const T& container) {
  ZingPool& zpool = mutable_pool();
  key_value_t* location = find_schema(attr,zpool);
  if (!location) MSYS_FAIL("attibute error: "+attr);

  // Must be the right size
  if (container.size() != size()) MSYS_FAIL("size mismatch in column");

  std::vector<Zing>& column = location->values;
  size_t i = 0;
  for(typename T::const_iterator p = container.begin(), en=container.end();
      p != en; ++p,++i) {
    typename T::value_type value = *p;
    std::string rval = precision_string(&value,location->type,location->precision);
    Zing zz(rval,zpool);
    column[i] = zz;
  }
}

template void desres::msys::DestroArray::set_column<std::vector<bool> >(const std::string& attr, const std::vector<bool>& container);
template void desres::msys::DestroArray::set_column<std::vector<int> >(const std::string& attr, const std::vector<int>& container);
template void desres::msys::DestroArray::set_column<std::vector<short> >(const std::string& attr, const std::vector<short>& container);
template void desres::msys::DestroArray::set_column<std::vector<long> >(const std::string& attr, const std::vector<long>& container);
template void desres::msys::DestroArray::set_column<std::vector<unsigned int> >(const std::string& attr, const std::vector<unsigned int>& container);
template void desres::msys::DestroArray::set_column<std::vector<unsigned short> >(const std::string& attr, const std::vector<unsigned short>& container);
template void desres::msys::DestroArray::set_column<std::vector<unsigned long> >(const std::string& attr, const std::vector<unsigned long>& container);
template void desres::msys::DestroArray::set_column<std::vector<float> >(const std::string& attr, const std::vector<float>& container);
template void desres::msys::DestroArray::set_column<std::vector<double> >(const std::string& attr, const std::vector<double>& container);
