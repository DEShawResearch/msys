/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"
#include <sstream>
#include <cassert>

// -----------------------------------------------
// P A R S E R
// -----------------------------------------------

/*FORWARD*/
static void
predict_block(desres::Destro& M, desres::Destro::Tokenizer& tokenizer);

/**************************************************************************/
/* LOCAL  **************       predict_schema      ************************/
/**************************************************************************/
/* Use the tokenizer to return blocks looking for the "schema" pattern    */
/* We predict that we get new schema names (b_ i_ r_ or s_) until we see  */
/* the ::: token that marks the end of this section.  Note the use of the */
/* optional_comment() method from tokenizer to look for the doc string    */
/*                                                                        */
/* This is factored out since we find these sections at the beginning of  */
/* normal blocks (followed by their values) and array sections (followed  */
/* by rows)                                                               */
/**************************************************************************/
static std::vector<desres::Destro::schema_t>
predict_schema(desres::Destro::Tokenizer& tokenizer) {
  std::vector<desres::Destro::schema_t> schemas;
  while(tokenizer.not_a(":::")) {
    desres::Destro::schema_t schema;
    std::string token = tokenizer.token();
    if (token[0] != 'b' && token[0] != 'i' && token[0] != 'r' && token[0] != 's') {
      std::stringstream str;
      str << "Line " << tokenizer.line() << " predicted a schema, but " 
          << token << " didn't start b_ i_ r_ or s_ ";
      throw desres::dessert(str.str(),DESSERT_LOC);
    }
    schema.type = token[0];
    schema.attr = token.substr(2);
    schema.doc  = tokenizer.optional_comment();
    schemas.push_back(schema);
    tokenizer.next();
  }
  return schemas;
}

/**************************************************************************/
/* LOCAL  ************** predict_schema_and_values ************************/
/**************************************************************************/
/* Normal blocks start with a list of schema followed ::: followed by a   */
/* matching number of data values.  We read them without interpretation   */
/**************************************************************************/
static void
predict_schema_and_values(desres::Destro& M, desres::Destro::Tokenizer& tokenizer) {
  std::vector<desres::Destro::schema_t> schema = predict_schema(tokenizer);
  tokenizer.predict(":::");
  for(unsigned i=0;i<schema.size();++i) {
    std::string value = tokenizer.predict_value();
    if (value == "<>" || value == "") {
      M.add_schema(schema[i].type,schema[i].attr,schema[i].doc);
    } else {
      // Strip quotes if present
      if (value[0] == '"' && value[value.size()-1]) {
        value = value.substr(1,value.size()-2);
      }
      M.add_schema_and_value(schema[i].type,schema[i].attr,schema[i].doc,value);
    }
  }
}


/**************************************************************************/
/* LOCAL  **************     predict_blockbody     ************************/
/**************************************************************************/
/* A blockbody looks like {, schema_and_values, zero or more blocks, and  */
/* a closing }.  This is factored out so that Maeff below can read just   */
/* the opening version block (which has no name).                         */
/**************************************************************************/
static void predict_blockbody(desres::Destro& subblock, desres::Destro::Tokenizer& tokenizer) {
  tokenizer.predict("{");
  predict_schema_and_values(subblock,tokenizer);
  while(tokenizer.not_a("}")) {
    predict_block(subblock,tokenizer);
  }
  tokenizer.predict("}");
}

/**************************************************************************/
/* LOCAL  **************     predict_arraybody     ************************/
/**************************************************************************/
/* */
/**************************************************************************/
static void predict_arraybody(desres::DestroArray& subarray, desres::Destro::Tokenizer& tokenizer) {

  // Read header
  tokenizer.predict("[");
  tokenizer.predict();
  tokenizer.predict("]");
  tokenizer.predict("{");

  // Read schema
  std::vector<desres::Destro::schema_t> schema = predict_schema(tokenizer);
  size_t width = schema.size();
  for(unsigned i=0;i<width;++i) {
    subarray.add_schema(schema[i].type,schema[i].attr,schema[i].doc);
  }
  tokenizer.predict(":::");

  // Read rows
  std::vector<std::string> elements(width);
  while(tokenizer.not_a(":::")) {
    std::string offset = tokenizer.predict();
    for(unsigned i=0;i<width;++i) {
      elements[i] = tokenizer.predict_value();
    }
    //desres::Destro& row = subarray.bulk_row(offset,elements);
    subarray.bulk_row(offset,elements);
  }

  tokenizer.predict(":::");

  tokenizer.predict("}");
}

/**************************************************************************/
/* LOCAL  **************   predict_nameless_block  ************************/
/**************************************************************************/
/* Assume the name from elsewhere, just read the block or array body      */
/**************************************************************************/
static void predict_nameless_block(std::string name,desres::Destro& M, desres::Destro::Tokenizer& tokenizer) {
  // -----------------------------------------------
  // May be an array
  // -----------------------------------------------
  if (std::string(tokenizer.token()) == "[") {
    desres::DestroArray& subarray = M.new_array(name);
    predict_arraybody(subarray,tokenizer);
  } 

  // -----------------------------------------------
  // Otherwise just a block
  // -----------------------------------------------
  else {
    desres::Destro& subblock = M.new_block(name);
    predict_blockbody(subblock,tokenizer);
  }
}

/**************************************************************************/
/* LOCAL  **************         check_name        ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
static void
check_name(const desres::Destro::Tokenizer& tokenizer,const std::string& name) {
  if (name.size() > 0 && !(isalpha(name[0]) || name[0] == '_')) {
    std::stringstream str;
    str << "Line " << tokenizer.line() << " predicted a block name have " << name << std::endl;
    throw desres::dessert(str.str());
  }
}

/**************************************************************************/
/* LOCAL  **************       predict_block       ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
static void
predict_block(desres::Destro& M, desres::Destro::Tokenizer& tokenizer) {
  std::string name = tokenizer.predict();
  check_name(tokenizer,name);
  predict_nameless_block(name,M,tokenizer);
}

// -----------------------------------------------
// P R O T E C T E D
// -----------------------------------------------
void desres::Destro::fill_named(Destro& block,Tokenizer& tokenizer) {
  std::string name = tokenizer.predict();
  check_name(tokenizer,name);
  predict_blockbody(block,tokenizer);
  block.name(name);
}

void desres::Destro::fill_nameless(Destro& block,Tokenizer& tokenizer) {
  predict_blockbody(block,tokenizer);
}

/*!
 * 
 * @param attr The attribute to set
 * @param value The boolean value to create
 */
void desres::Destro::set(const std::string& attr,const bool& value) {
  add_schema('b',attr);
  get_attr(attr) = value;
}

/*!
 * 
 * @param attr The attribute to set
 * @param value The integer value to create
 */
void desres::Destro::set(const std::string& attr,const long& value) {
  add_schema('i',attr);
  get_attr(attr) = value;
}

/*!
 * 
 * @param attr The attribute to set
 * @param value The integer value to create
 */
void desres::Destro::set(const std::string& attr,const int& value) {
  add_schema('i',attr);
  get_attr(attr) = value;
}

/*!
 * 
 * @param attr The attribute to set
 * @param value The floating point value to create
 */
void desres::Destro::set(const std::string& attr,const double& value) {
  add_schema('r',attr);
  get_attr(attr) = value;
}

/*!
 * 
 * @param attr The attribute to set
 * @param value The string value to create
 */
void desres::Destro::set(const std::string& attr,const std::string& value) {/*GCOV-IGNORE*/
  set(attr,value.c_str());/*GCOV-IGNORE*/
}/*GCOV-IGNORE*/

/*!
 * @param attr The attribute to set
 * @param value The string value to create
 */
void desres::Destro::set(const std::string& attr,const char* value) {
  //std::cerr << "char* set(" << attr << ")\n";
  add_schema('s',attr);
  get_attr(attr) = value;
}

/*!
 * @param attr The attribute to set
 * @param value The attribute value to copy
 */
void desres::Destro::set(const std::string& attr,const Attribute& value) { /*GCOV-IGNORE*/
  if (value.type() == 'e') throw dessert("invalid empty rhs"); /*GCOV-IGNORE*/
  add_schema(value.type(),attr); /*GCOV-IGNORE*/
  get_attr(attr).update( value.value().c_str() ); /*GCOV-IGNORE*/
} /*GCOV-IGNORE*/

desres::Zing desres::Destro::zingify(const std::string& value, ZingPool& zpool) {
  if (value == "<>") {
    return Zing();
  } else if (value == "\"<>\"") {
    return Zing("<>",zpool);
  }

  return Zing(value,zpool);
}

// -----------------------------------------------
//  P U B L I C
// -----------------------------------------------
const short desres::Destro::s_max_precision = 75;
int desres::Destro::s_default_float_precision = 7;
int desres::Destro::s_default_double_precision = 7;

void desres::Destro::test_private_parts() {
  ZingPool pool;
  Zing z1;
  Zing z2("",pool);
  Zing z3("\n",pool);
  assert(quotify(z1,pool) == "<>");
  assert(quotify(z2,pool) == "\"\"");
  try {
    quotify(z3,pool);
    assert(false); /*GCOV-IGNORE*/
  } catch(...) {}

  Zing z4 = zingify("<>",pool);
  Zing z5 = zingify("\"<>\"",pool);
  Zing z6 = zingify("abc",pool);

  assert(z4.is_empty());
  assert(z5.string(pool) == "<>");
  assert(z6.string(pool) == "abc");
}

desres::Destro::Destro(Destro* parent)
  : m_parent(parent)
{
}

desres::Destro::~Destro() {
}

std::string desres::Destro::quotify(desres::Zing z, const desres::ZingPool& zpool) {
  // uninitialized --> quoteless <>
  if (z.is_empty()) return "<>";

  std::string raw = z.string(zpool);

  // empty string --> quoted ""
  if (raw == "") return "\"\"";


  // Check for non-printable characters and "
  for(std::string::iterator p=raw.begin(), en=raw.end();
      p != en; ++p) {
    if (isspace(*p) || !isprint(*p) || *p == '"' || *p == '<' || *p == '\\') {
      std::string escaped(raw.begin(),p);
      //std::cerr << '?' << raw << "?" << escaped << '?' << std::endl;
      for(;p!=en;++p) {
        // We only support space and tab
        if (isspace(*p) && !(*p == ' ' || *p == '\t')) {
          throw dessert("unprintable whitespace in '" + raw + '\'');
        }

        if (*p == '"') {
          escaped += "\\\"";
        } else if (*p == '\\') {
          escaped += "\\\\";
        } else {
          escaped += *p;
        }

      }
      raw = '"' + escaped + '"';
      break;
    }
  }

  return raw;
}

std::string desres::Destro::indent(int level, int& newlevel) {
  if (level < 0) {
    newlevel = -1;
    return " ";
  }

  newlevel = level + 1;
  return std::string(level*2,' ');
}

std::string desres::Destro::schema_out(char type, Zing key, Zing doc, const ZingPool& pool) {
  std::string stype(&type,1);
  std::string skey = key.string(pool);
  std::string schema = stype + '_' + skey;

  std::string sdoc = doc.string(pool);
  if (sdoc.size() > 0) {
    std::string comment(" # ");

    for(std::string::iterator p = sdoc.begin(), en=sdoc.end();
        p != en; ++p) {
      if (*p != '#' && *p != '\n') comment += *p;
    }

    comment += " #";
    schema += comment;
  }

  return schema;
}

desres::Destro* desres::Destro::parent() {
  return m_parent;
}

const desres::Destro* desres::Destro::parent() const {
  return m_parent;
}

desres::Destro& desres::Destro::progenitor() {
  if (!m_parent) return *this;
  return m_parent->progenitor();
}

const desres::Destro& desres::Destro::progenitor() const {
  if (!m_parent) return *this;
  return m_parent->progenitor();
}

desres::ZingPool& desres::Destro::pool() {
  Destro& top = progenitor();

  if (&top == this) {
    throw dessert("pool not defined",DESSERT_LOC);
  }
  return top.pool();
}

desres::ZingPool& desres::Destro::mutable_pool() {
  try {
    return pool();
  } catch (const std::exception&) {
  }

  throw dessert("item is immutable",DESSERT_LOC);
}

const desres::ZingPool& desres::Destro::pool() const {
  const Destro& top = progenitor();

  if (&top == this) {
    throw dessert("pool not defined",DESSERT_LOC); /*GCOV-IGNORE*/
  }
  return top.pool();
}

desres::Destro::operator std::string() const {
  std::ostringstream ss;
  write(ss,-1);
  return ss.str();
}

void desres::Destro::add_schema_and_value(char type,const std::string& attr,const std::string& doc,const std::string& value) {
  add_schema(type,attr,doc);
  set_attr(type,attr,value);
}


desres::Destro::Attribute desres::Destro::operator[](const std::string& attr) {
  if (contains(attr)) return Attribute(this,attr);
  throw dessert("Attribute error: "+attr);
}

const desres::Destro::Attribute desres::Destro::operator[](const std::string& attr) const {
  return Attribute(const_cast<Destro*>(this),attr);
}

desres::Destro& desres::Destro::operator[](ssize_t blockno) {
  return block(blockno);
}

const desres::Destro& desres::Destro::operator[](ssize_t blockno) const { /*GCOV-IGNORE*/
  return block(blockno);/*GCOV-IGNORE*/
}

bool desres::Destro::is_array() const {
  return false;
}

desres::Destro* desres::Destro::operator->() {
  return this;
}

const desres::Destro* desres::Destro::operator->() const {
  return this;
}

bool desres::Destro::contains(const std::string& attr) const {
  return has_attr(attr) || has_block(attr);
}

size_t desres::Destro::footprint() const {
  return sizeof(Destro);
}

void desres::Destro::compact() {
  ZingPool& zpool = pool();

  zpool.clearbits();
  progenitor().touch(zpool);
  zpool.finish_compaction();
}

desres::Destro& desres::Destro::append(const Destro& source) {

  // -----------------------------------------------
  // Horrid ugliness.... Just leaving it be for now
  // but really needs to be factored back into the
  // DestroBlock and DestroArray classes
  // -----------------------------------------------

  // -----------------------------------------------
  // A R R A Y S
  // -----------------------------------------------
  if (dynamic_cast<const DestroArray*>(&source)) {
    DestroArray* newarray = &new_array(source.name());

    // Yes, we can try to copy onto ourself!
    if (newarray == this) return *this;

    // Copy over just the schemas
    std::vector<Destro::schema_t> schemas = source.ordered_schema();
    for(std::vector<Destro::schema_t>::iterator p = schemas.begin();
        p != schemas.end(); ++p) {
      Destro::schema_t& schema = *p;
      newarray->add_schema(schema.type,schema.attr,schema.doc);
      newarray->set_precision(schema.attr,source.get_precision(schema.attr));
    }

    // Add matching rows
    for(size_t i = 1; i <= source.size(); ++i) {
      Destro& newrow = newarray->new_block("");
      const Destro& oldrow = source[i];
      for(std::vector<Destro::schema_t>::iterator p = schemas.begin();
          p != schemas.end(); ++p) {
        Destro::schema_t& schema = *p;
        if (oldrow.has_value(schema.attr)) {
          newrow.set_attr(schema.type,schema.attr,
                          oldrow.get_value(schema.attr));
        }
      }
      //newblock->append(source.block(i));
    }

    return *newarray;
  }


  // -----------------------------------------------
  // B L O C K S
  // -----------------------------------------------
  Destro* newblock = &new_block(source.name());

  // Yes, we can try to copy onto ourself!
  if (newblock == this) return *this;

  // Copy over the schema and values
  std::vector<Destro::schema_t> schemas = source.ordered_schema();
  for(std::vector<Destro::schema_t>::iterator p = schemas.begin();
      p != schemas.end(); ++p) {
    Destro::schema_t& schema = *p;
    newblock->add_schema(schema.type,schema.attr,schema.doc);
    if (source.has_value(schema.attr)) {
      newblock->get_attr(schema.attr).update(source.get_attr(schema.attr).value().c_str());
    }
  }
    
  // Copy over the subblocks
  for(size_t i = 1; i <= source.size(); ++i) {
    newblock->append(source.block(i));
  }

  return *newblock;
}

desres::Destro& desres::Destro::append() {
  return new_block("");
}

desres::DestroArray& desres::Destro::array(size_t i) {
  Destro& ablock = block(i);
  DestroArray* arr = dynamic_cast<DestroArray*>(&ablock);
  if (!arr) {
    throw dessert(ablock.name()+" is not an array");
  }
  return *arr;
}

const desres::DestroArray& desres::Destro::array(size_t i) const {
  const Destro& ablock = block(i);
  const DestroArray* arr = dynamic_cast<const DestroArray*>(&ablock);
  if (!arr) {
    throw dessert(ablock.name()+" is not an array");
  }
  return *arr;
}


// -----------------------------------------------
// global functions
// -----------------------------------------------
std::ostream& desres::operator<<(std::ostream& os,const desres::Destro& M) {
  M.write(os);
  return os;
}
