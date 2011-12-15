/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"
#include <sstream>

/*
  A thought to cut space... you may be able to remove the
  m_row member because you can compute the row number given
  just a reference to the row object (look it up in the
  deque, do an offset analysis to figure out where it is
  in a vector!).

  This would cut the footprint by 8 bytes per row
 */

// -----------------------------------------------
// P R I V A T E
// -----------------------------------------------

desres::msys::DestroArray* desres::msys::DestroArray::DestroRow::owner() {
  DestroArray* array = dynamic_cast<DestroArray*>(m_parent);
  if (array == NULL) {
    throw dessert("invalid parent for row",DESSERT_LOC); /*GCOV-IGNORE*/
  }
  return array;
}

const desres::msys::DestroArray* desres::msys::DestroArray::DestroRow::owner() const {
  const DestroArray* array = dynamic_cast<const DestroArray*>(m_parent);
  if (array == NULL) {
    throw dessert("invalid parent for row",DESSERT_LOC); /*GCOV-IGNORE*/
  }
  return array;
}


// -----------------------------------------------
// P U B L I C
// -----------------------------------------------
desres::msys::DestroArray::DestroRow::DestroRow(Destro* parent,size_t row)
  : Destro(parent),m_row(row)
{
}

desres::msys::DestroArray::DestroRow::~DestroRow() {
}

void desres::msys::DestroArray::DestroRow::setrow(size_t i) {
  m_row = i-1;
}


/*!
 * Rows always have zero size (no subblocks)
 */
size_t desres::msys::DestroArray::DestroRow::size() const {
  return 0;
}

/*!
 * The row name is it's one-based offset
 */
std::string desres::msys::DestroArray::DestroRow::name() const {
  std::ostringstream ss;
  ss << (m_row+1);
  return ss.str();
}

/*!
 * The row name is unsettable
 */
void desres::msys::DestroArray::DestroRow::name(const std::string& name) {
  throw dessert("cannot set row name");
}

/*!
 * Adding a schema to a row adds it to the enclosing array
 */
void desres::msys::DestroArray::DestroRow::add_schema(char type,const std::string& attr,const std::string& doc) {
  owner()->add_schema(type,attr,doc);
}

/*!
 * Adding a schema to a row adds it to the enclosing array.
 * Then we set just this row's value.
 */
void desres::msys::DestroArray::DestroRow::add_schema_and_value(char type,const std::string& attr,const std::string& doc,const std::string& value) {
  owner()->add_schema(type,attr,doc);

  set_attr(type,attr,value);
}

/*!
 * The enclosing array knows all the schemas
 */
std::map<std::string,desres::msys::Destro::schema_t> desres::msys::DestroArray::DestroRow::schemas() const {
  return owner()->schemas();
}

/*!
 * Get schema list from enclosing array.
 */
std::vector<desres::msys::Destro::schema_t> desres::msys::DestroArray::DestroRow::ordered_schema() const {
  return owner()->ordered_schema();
}

/*!
 * The enclosing array has the schema/values.  Look for the
 * value set and then extract this row's value.
 */
std::string desres::msys::DestroArray::DestroRow::get_value(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = owner()->find_schema(attr,zpool);
  if (location == NULL) {
    std::string message("Attribute Error: ");
    throw dessert(message+attr);
  }

  if (location->values.at(m_row).is_empty()) {
    throw dessert(attr+" is empty");
  }

  return location->values.at(m_row).string(zpool);
}

/*!
 * We get attribute type info from enclosing array.
 */
char desres::msys::DestroArray::DestroRow::get_type(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = owner()->find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Enclosing array Attribute Error: "+attr);
  }

  return location->type;
}


/*!
 * Precision is pulled from the enclosing array
 */
int desres::msys::DestroArray::DestroRow::get_precision(const std::string& attr) const {
  return owner()->get_precision(attr);
}

/*!
 * Precision is set in the enclosing array (effects all rows)
 */
void desres::msys::DestroArray::DestroRow::set_precision(const std::string& attr,int precision) {
  owner()->set_precision(attr,precision);
}

std::string desres::msys::DestroArray::DestroRow::get_doc(const std::string& attr) const {
  return owner()->get_doc(attr);
}

void desres::msys::DestroArray::DestroRow::set_doc(const std::string& attr, const std::string& doc) {
  owner()->set_doc(attr,doc);
}

desres::msys::Destro::Attribute desres::msys::DestroArray::DestroRow::get_attr(const std::string& attr) {
  return Attribute(this,attr);
}

const desres::msys::Destro::Attribute desres::msys::DestroArray::DestroRow::get_attr(const std::string& attr) const {
  return Attribute(const_cast<DestroRow*>(this),attr);
}

void desres::msys::DestroArray::DestroRow::set_unsafe(const std::string& attr,const std::string& value) {
  ZingPool& zpool = pool();
  key_value_t* location = owner()->find_schema(attr,zpool);
  if (location == NULL) {
    throw dessert("Attribute error: " + attr);
  }

  Zing zval(value,zpool);
  location->values[m_row] = zval;
}

/*!
 * We look up the data descriptor in the enclosing array and then
 * set our row's data value.
 */
void desres::msys::DestroArray::DestroRow::set_attr(char type,const std::string& attr,const std::string& value) {
  ZingPool& zpool = pool();
  key_value_t* location = owner()->find_schema(attr,zpool);
  if (location == NULL) {
    owner()->set_attr(type,attr,value);
    return;
  }

  Zing zval(value,zpool);
  location->values[m_row] = zval;
}

void desres::msys::DestroArray::DestroRow::clear(const std::string& attr) {
  owner()->clear(attr);
}

/*!
 * Deleting an attribute deletes it through the whole array
 */
void desres::msys::DestroArray::DestroRow::del(const std::string& attr) {
  owner()->del(attr);
}

/*!
 * No blocks to delete.  We silently return (since that block
 * clearly doesn't exist anymore!).
 */
void desres::msys::DestroArray::DestroRow::del(size_t blockno) {
}

/*!
 * Cannot add blocks to a row.
 */
desres::msys::Destro& desres::msys::DestroArray::DestroRow::new_block(const std::string& name) {
  throw dessert("new_block() not applicable in rows");
}

/*!
 * Cannot add blocks to a row.
 */
desres::msys::DestroArray& desres::msys::DestroArray::DestroRow::new_array(const std::string& name,size_t num_elements) {
  throw dessert("new_array() not applicable in rows");
}

desres::msys::Destro::Attribute desres::msys::DestroArray::DestroRow::operator[](const std::string& attr) {
  return Attribute(this,attr);
}

const desres::msys::Destro::Attribute desres::msys::DestroArray::DestroRow::operator[](const std::string& attr) const {/*GCOV-IGNORE*/
  return Attribute(const_cast<DestroRow*>(this),attr);/*GCOV-IGNORE*/
}

/*!
 * No blocks, so no operator[int]
 */
desres::msys::Destro& desres::msys::DestroArray::DestroRow::operator[](ssize_t) {
  throw dessert("operator[] not applicable in rows");
}

/*!
 * No blocks, so no operator[int]
 */
const desres::msys::Destro& desres::msys::DestroArray::DestroRow::operator[](ssize_t) const {
  throw dessert("operator[] not applicable in rows");
}

/*!
 * Rows contain no blocks
 */
bool desres::msys::DestroArray::DestroRow::has_block(const std::string& blockname) const {
  return false;
}

/*!
 * The row has only the attributes of the enclosing array
 */
bool desres::msys::DestroArray::DestroRow::has_attr(const std::string& name) const {
  return owner()->has_attr(name);
}

/*!
 * Find value set in enclosing array and test this row's value
 */
bool desres::msys::DestroArray::DestroRow::has_value(const std::string& attr) const {
  const ZingPool& zpool = pool();
  const key_value_t* location = owner()->find_schema(attr,zpool);
  if (!location) return false;
  return !(location->values.at(m_row).is_empty() ||
           location->values.at(m_row).string(zpool) == "<>");
}

/*!
 * Rows have no subblocks
 */
desres::msys::Destro& desres::msys::DestroArray::DestroRow::block(size_t i) {
  throw dessert("block() not applicable to rows");
}

/*!
 * Rows have no subblocks
 */
const desres::msys::Destro& desres::msys::DestroArray::DestroRow::block(size_t i) const {
  throw dessert("block() not applicable to rows");
}

/*!
 * Rows have no subblocks
 */
desres::msys::Destro& desres::msys::DestroArray::DestroRow::block(const std::string& name) {
  throw dessert("block() not applicable to rows");
}

/*!
 * Rows have no subblocks
 */
const desres::msys::Destro& desres::msys::DestroArray::DestroRow::block(const std::string& name) const {
  throw dessert("block() not applicable to rows");
}

/*!
 * Rows always write horizontally
 */
void desres::msys::DestroArray::DestroRow::write(std::ostream& os, int level) const {
  //const char* lineterm = (level >= 0)?"\n":"";
  const ZingPool& zpool = pool();

  os << m_row+1;
  for(std::vector<key_value_t>::const_iterator p=owner()->m_data.begin(),
        en = owner()->m_data.end(); p != en; ++p) {
    const key_value_t& schema = *p;
    os << ' ' << quotify(schema.values[m_row],zpool);
  }
}

/*!
 * Rows have no deep storage since their values are held by
 * the enclosing array.
 */
size_t desres::msys::DestroArray::DestroRow::footprint() const {
  return Destro::footprint() + sizeof(m_row);
}

/*!
 * Row's don't directly contain any Zing data...  And in any case,
 * DestroArray doesn't call this, so nobody should every invoke
 * the function.
 */
void desres::msys::DestroArray::DestroRow::touch(ZingPool& zpool) const { /*GCOV-IGNORE*/
}

