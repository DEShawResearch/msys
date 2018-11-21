/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"

// -----------------------------------------------
// P R I V A T E
// -----------------------------------------------
desres::msys::DestroBlock::key_value_t*
desres::msys::DestroBlock::find_schema(const std::string& attr,const desres::msys::ZingPool& pool) {
  for(std::vector<key_value_t>::iterator p = m_data.begin(),
      en = m_data.end(); p != en; ++p) {
    if ((*p).key.string(pool) == attr) return &(*p);
  }
  return NULL;
}

const desres::msys::DestroBlock::key_value_t*
desres::msys::DestroBlock::find_schema(const std::string& attr,const desres::msys::ZingPool& pool) const {
  for(std::vector<key_value_t>::const_iterator p = m_data.begin(),
      en = m_data.end(); p != en; ++p) {
    if ((*p).key.string(pool) == attr) return &(*p);
  }
  return NULL;
}


ssize_t desres::msys::DestroBlock::offset_of_block(const std::string& name) const {
  for(std::vector<Destro*>::const_iterator p=m_subblocks.begin(),
      en = m_subblocks.end(); p != en; ++p) {
    if ((*p)->name() == name) {
      return p-m_subblocks.begin();
    }
  }
  return -1;
}

// -----------------------------------------------
// P R O T E C T E D
// -----------------------------------------------
desres::msys::Destro& desres::msys::DestroBlock::add_block(const std::string& name) {
  Destro* block = NULL;
  try {
    block = new DestroBlock(this);
    block->name(name);
    m_subblocks.push_back(block);
    return *block;
  } catch(...) { /*GCOV-IGNORE*/
    delete block; /*GCOV-IGNORE*/
    throw; /*GCOV-IGNORE*/
  }
}

void desres::msys::DestroBlock::touch(ZingPool& zpool) const {
  // Touch my name
  DestroNamedBlock::touch(zpool);

  // Touch my subblocks
  for(std::vector<Destro*>::const_iterator p = m_subblocks.begin(),
        en = m_subblocks.end(); p != en; ++p) {
    (*p)->touch(zpool);
  }

  // Touch my data
  for(std::vector<key_value_t>::const_iterator p = m_data.begin(),
        en = m_data.end(); p != en; ++p) {
    const key_value_t data = *p;
    data.key.touch(zpool);
    data.value.touch(zpool);
    data.doc.touch(zpool);
  }

}

// -----------------------------------------------
//  P U B L I C
// -----------------------------------------------
desres::msys::DestroBlock::DestroBlock(Destro* parent)
  : DestroNamedBlock(parent)
{
}

desres::msys::DestroBlock::~DestroBlock()
{
  // Clean up subblocks
  for(std::vector<Destro*>::const_iterator p=m_subblocks.begin(),
      en = m_subblocks.end(); p != en; ++p) {
    delete *p;
  }
}

// size_t size() const;
size_t desres::msys::DestroBlock::size() const {
  return m_subblocks.size();
}


void desres::msys::DestroBlock::add_schema(char type,const std::string& attr,const std::string& doc) {
  ZingPool& zpool = pool();

  key_value_t* location = find_schema(attr,zpool);
  if (location) {
    if (location->type != type) throw dessert("schema exists with different type");
    if (doc.size() > 0) location->doc = Zing(doc,zpool);
  }
  
  // We don't have that key, so we add a new one
  else {
    key_value_t keyvalue;
    keyvalue.key = Zing(attr,zpool);
    keyvalue.value = Zing();
    keyvalue.doc = Zing(doc,zpool);
    keyvalue.type = (unsigned)type;
    keyvalue.precision = -1;
    m_data.push_back(keyvalue);
  }
}

void desres::msys::DestroBlock::add_schema_and_value(char type,const std::string& attr,const std::string& doc,const std::string& value) {
  add_schema(type,attr,doc);
  set_attr(type,attr,value);
}

std::map<std::string,desres::msys::Destro::schema_t> desres::msys::DestroBlock::schemas() const {
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

std::vector<desres::msys::Destro::schema_t> desres::msys::DestroBlock::ordered_schema() const {
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

std::string desres::msys::DestroBlock::get_value(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute Error: "+attr);
  }
  if (location->value.is_empty()) {
    throw dessert(attr+" is empty");
  }
  return location->value.string(zpool);
}

char desres::msys::DestroBlock::get_type(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute Error: "+attr);
  }
  return location->type;
}

int  desres::msys::DestroBlock::get_precision(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute Error: "+attr);
  }
  return location->precision;
}

void desres::msys::DestroBlock::set_precision(const std::string& attr,int precision) {
  const ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute Error: "+attr);
  }
  if (precision > s_max_precision) precision = s_max_precision;
  location->precision = precision;
}

std::string desres::msys::DestroBlock::get_doc(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute Error: "+attr);
  }
  return location->doc.string(zpool);
}

void desres::msys::DestroBlock::set_doc(const std::string& attr, const std::string& doc) {
  ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute Error: "+attr);
  }
  location->doc = Zing(doc,zpool);
}

desres::msys::Destro::Attribute desres::msys::DestroBlock::get_attr(const std::string& attr) {
  return Attribute(this,attr);
}

const desres::msys::Destro::Attribute desres::msys::DestroBlock::get_attr(const std::string& attr) const {
  return Attribute(const_cast<DestroBlock*>(this),attr);
}

void desres::msys::DestroBlock::set_unsafe(const std::string& attr,const std::string& value) {
  ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  raw_update(location,zpool,value);
}

void desres::msys::DestroBlock::set_attr(char type,const std::string& attr,const std::string& value) {
  ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location == NULL) {
    add_schema(type,attr);
    location = find_schema(attr,zpool);
  } else {
    location->type = type;
  }
  raw_update(location,zpool,value);
}

void desres::msys::DestroBlock::clear(const std::string& attr) {
  ZingPool& zpool = pool();
  key_value_t* location = find_schema(attr,zpool);
  if (location != NULL) {
    location->value = Zing();
  }
}

void desres::msys::DestroBlock::del(const std::string& attr) {

  // Perhaps it is an attribute....
  ZingPool& zpool = pool();
  for(std::vector<key_value_t>::iterator p = m_data.begin(),
      en = m_data.end(); p != en; ++p) {
    if ((*p).key.string(zpool) == attr) {
      m_data.erase(p);
      return;
    }
  }

  // Perhaps it is a block
  for(std::vector<Destro*>::iterator p = m_subblocks.begin(),
      en = m_subblocks.end(); p != en; ++p) {
    if ((*p)->name() == attr) {
      delete *p;
      m_subblocks.erase(p);
      return;
    }
  }

  // Otherwise, silently ignore
}

void desres::msys::DestroBlock::del(size_t blockno) {
  if (blockno >= 1 && blockno <= m_subblocks.size()) {
    m_subblocks.erase(m_subblocks.begin()+(blockno-1));
  }

  // Silently ignore non-existant blocks
}

desres::msys::Destro& desres::msys::DestroBlock::new_block(const std::string& name) {
  // Might already exist
  ssize_t offset = offset_of_block(name);
  if (offset >= 0) {
    if (m_subblocks[offset]->is_array()) {
      throw dessert("attempting to convert an array with new_block()");
    }
    return *m_subblocks[offset];
  }

  return add_block(name);
}

desres::msys::DestroArray& desres::msys::DestroBlock::new_array(const std::string& name, size_t num_elements) {
  // Might already exist
  ssize_t offset = offset_of_block(name);
  if (offset >= 0) {
    DestroArray* block = dynamic_cast<DestroArray*>(m_subblocks[offset]);
    if (block) {
      block->resize(num_elements);
      return *block;
    }
    throw dessert("attempting to convert a block with new_array()");
  }

  DestroArray* arr = NULL;
  try {
    arr = new DestroArray(this,num_elements);
    arr->name(name);
    m_subblocks.push_back(arr);
    return *arr;
  } catch(...) { /*GCOV-IGNORE*/
    delete arr;/*GCOV-IGNORE*/
    throw;/*GCOV-IGNORE*/
  }
}

desres::msys::Destro::Attribute desres::msys::DestroBlock::operator[](const std::string& attr) {
  return Destro::operator[](attr);
}

const desres::msys::Destro::Attribute desres::msys::DestroBlock::operator[](const std::string& attr) const {
  return Destro::operator[](attr);
}

desres::msys::Destro& desres::msys::DestroBlock::operator[](ssize_t i) {
  return Destro::operator[](i);
}

const desres::msys::Destro& desres::msys::DestroBlock::operator[](ssize_t i) const { /*GCOV-IGNORE*/
  return Destro::operator[](i);/*GCOV-IGNORE*/
}

bool desres::msys::DestroBlock::has_block(const std::string& name) const {
  return offset_of_block(name) >= 0;
}

bool desres::msys::DestroBlock::has_attr(const std::string& attr) const {
  const key_value_t* location = find_schema(attr,pool());
  return (location != NULL);
}

bool desres::msys::DestroBlock::has_value(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = find_schema(attr,zpool);
  if (!location) return false;
  return !(location->value.is_empty() || location->value.string(zpool) == "<>");
}

desres::msys::Destro& desres::msys::DestroBlock::block(size_t i) {
  return *m_subblocks.at(i-1);
}

const desres::msys::Destro& desres::msys::DestroBlock::block(size_t i) const {
  return *m_subblocks.at(i-1);
}

desres::msys::Destro& desres::msys::DestroBlock::block(const std::string& name) {

  for(std::vector<Destro*>::iterator p = m_subblocks.begin(),
        en=m_subblocks.end(); p != en; ++p) {
    //if ((*p)->name() == name) {
    //std::cerr << "Found block " << name << " at " << *p << std::endl;
    //}
    if ((*p)->name() == name) return **p;
  }

  throw dessert("Attribute error: " + name);
}

const desres::msys::Destro& desres::msys::DestroBlock::block(const std::string& name) const {

  for(std::vector<Destro*>::const_iterator p = m_subblocks.begin(),
        en=m_subblocks.end(); p != en; ++p) {
    if ((*p)->name() == name) return **p;
  }

  throw dessert("Attribute error: " + name);
}

void desres::msys::DestroBlock::write(std::ostream& os, int level) const {
  int newlevel = 0;
  std::string indentation = indent(level,newlevel);
  int dummy = 0;
  std::string indentation2 = indent(newlevel,dummy);

  const char* lineterm = (level >= 0)?"\n":"";

  os << indentation;

  // Write the name (if any) followed by a space
  std::string blockname = name();
  if (blockname.size() > 0) os << blockname << ' ';


  os << '{' << lineterm;

  const ZingPool& zpool = pool();

  // Write the attributes...
  for(std::vector<key_value_t>::const_iterator p=m_data.begin(),
        en = m_data.end(); p != en; ++p) {
    const key_value_t& schema = *p;
    os << indentation2 << schema_out(schema.type,schema.key,schema.doc,zpool) << lineterm;
  }
  os << indentation2 << ":::" << lineterm;
  for(std::vector<key_value_t>::const_iterator p=m_data.begin(),
        en = m_data.end(); p != en; ++p) {
    const key_value_t& schema = *p;
    os << indentation2 << quotify(schema.value,zpool) << lineterm;
  }

  // Write the sub-blocks
  for(std::vector<Destro*>::const_iterator p=m_subblocks.begin(),
        en = m_subblocks.end(); p != en; ++p) {
    (*p)->write(os,newlevel);
  }
  os << indentation << '}' << lineterm;
}

size_t desres::msys::DestroBlock::footprint() const {
  size_t foot = DestroNamedBlock::footprint() +
    sizeof(m_subblocks) + sizeof(m_data);

  // Add in the space for the sub-blocks
  for(std::vector<Destro*>::const_iterator p = m_subblocks.begin(),
        en = m_subblocks.end(); p != en; ++p) {
    foot += sizeof(Destro*);
    foot += (*p)->footprint();
  }

  // Add in the space for the key-value pairs
  foot += m_data.size() * sizeof(m_data[0]);

  return foot;
}


// std::string name() const;
// void name(const std::string& name);
// desres::msys::ZingPool& pool();
// const desres::msys::ZingPool& pool() const;

// void add_schema(char type,const std::string& attr,const std::string& doc="");

// std::map<std::string,schema_t>& schemas() const;
// std::vector<schema_t> ordered_schema() const;

// std::string get_value(const std::string& attr) const;



// Attribute get_attr(const std::string& attr);
// const Attribute get_attr(const std::string& attr) const;

void desres::msys::DestroBlock::raw_update(key_value_t* location,ZingPool& zpool,std::string value) {
  if (location == NULL) throw dessert("no location provided");

  location->value = zingify(value,zpool);
}

/*GCOV-IGNORE*/
// void set_attr(char type,const std::string& attr,const std::string& value);


// void clear(const std::string& attr);
// void del(const std::string& attr);
// void del(size_t blockno);
// Destro& append(const Destro& source);
// Destro& append();

// Destro& new_block(const std::string& name);

// DestroArray& new_array(const std::string& name,size_t num_elements=0);
