/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"
#include <sstream>

// -----------------------------------------------
// P R I V A T E
// -----------------------------------------------
void desres::Maeff::init(Tokenizer& tokenizer) {
  fill_nameless(m_meta,tokenizer);

  while(tokenizer.not_a()) {
    std::string name = tokenizer.predict();
    fill_nameless(new_block(name),tokenizer);
  }
}

// -----------------------------------------------
//  P U B L I C
// -----------------------------------------------
desres::Maeff::Maeff()
  : DestroTop(), m_meta(this)
{
  m_meta.add_schema_and_value('s',"m_m2io_version","","2.0.0");
}

desres::Maeff::Maeff(Tokenizer& tokenizer)
  : DestroTop(), m_meta(this)
{
  init(tokenizer);
}

desres::Maeff::~Maeff()
{
}

desres::Destro& desres::Maeff::meta() {
  return m_meta;
}

const desres::Destro& desres::Maeff::meta() const {
  return m_meta;
}

std::string desres::Maeff::adjustname(const std::string& name) {
  // Name must be m_ct or ct
  if (name == "m_ct" || name == "ct" || name == "") {
    return "f_m_ct";
  }
  return name;
}

desres::Destro& desres::Maeff::new_block(const std::string& name) {
  std::string realname = adjustname(name);
  if (realname != "f_m_ct") {
    std::stringstream str;
    str << "Maeff block names must be f_m_ct, m_ct, or ct, not " << name << std::endl;
    throw dessert(str.str());
  }

  return add_block(realname);
}

void desres::Maeff::add_schema(char type,const std::string& attr,const std::string& doc) {
  throw dessert("use the meta() block to modify m_m2io_version");
}


bool desres::Maeff::has_block(const std::string& name) const {
  std::string realname = adjustname(name);
  return DestroTop::has_block(realname);
}

desres::Destro& desres::Maeff::block(size_t i) {
  return DestroTop::block(i);
}

const desres::Destro& desres::Maeff::block(size_t i) const {
  return DestroTop::block(i);
}

desres::Destro& desres::Maeff::block(const std::string& name) {
  std::string realname = adjustname(name);
  return DestroTop::block(realname);
}

const desres::Destro& desres::Maeff::block(const std::string& name) const {
  std::string realname = adjustname(name);
  return DestroTop::block(realname);
}

void desres::Maeff::write(std::ostream& os, int level) const {
  m_meta.write(os,level);
  size_t n = size();
  for(size_t i=1; i<=n; ++i) {
    this->block(i).write(os,level);
  }
}

size_t desres::Maeff::footprint() const {
  return DestroTop::footprint() + m_meta.footprint();
}

void desres::Maeff::touch(ZingPool& zpool) const {
  m_meta.touch(zpool);
  desres::DestroTop::touch(zpool);
}
