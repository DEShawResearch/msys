/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"
#include <sstream>

// -----------------------------------------------
// P R I V A T E
// -----------------------------------------------
void desres::msys::Maeff::init(Tokenizer& tokenizer) {
  fill_nameless(m_meta,tokenizer);

  while(tokenizer.not_a()) {
    std::string name = tokenizer.predict();
    fill_nameless(new_block(name),tokenizer);
  }
}

// -----------------------------------------------
//  P U B L I C
// -----------------------------------------------
desres::msys::Maeff::Maeff()
  : DestroTop(), m_meta(this)
{
  m_meta.add_schema_and_value('s',"m_m2io_version","","2.0.0");
}

desres::msys::Maeff::Maeff(Tokenizer& tokenizer)
  : DestroTop(), m_meta(this)
{
  init(tokenizer);
}

desres::msys::Maeff::~Maeff()
{
}

desres::msys::Destro& desres::msys::Maeff::meta() {
  return m_meta;
}

const desres::msys::Destro& desres::msys::Maeff::meta() const {
  return m_meta;
}

std::string desres::msys::Maeff::adjustname(const std::string& name) {
  // Name must be m_ct or ct
  if (name == "m_ct" || name == "ct" || name == "") {
    return "f_m_ct";
  }
  return name;
}

desres::msys::Destro& desres::msys::Maeff::new_block(const std::string& name) {
  std::string realname = adjustname(name);
  if (realname != "f_m_ct") {
    std::stringstream str;
    str << "Maeff block names must be f_m_ct, m_ct, or ct, not " << name << std::endl;
    throw dessert(str.str());
  }

  return add_block(realname);
}

void desres::msys::Maeff::add_schema(char type,const std::string& attr,const std::string& doc) {
  throw dessert("use the meta() block to modify m_m2io_version");
}


bool desres::msys::Maeff::has_block(const std::string& name) const {
  std::string realname = adjustname(name);
  return DestroTop::has_block(realname);
}

desres::msys::Destro& desres::msys::Maeff::block(size_t i) {
  return DestroTop::block(i);
}

const desres::msys::Destro& desres::msys::Maeff::block(size_t i) const {
  return DestroTop::block(i);
}

desres::msys::Destro& desres::msys::Maeff::block(const std::string& name) {
  std::string realname = adjustname(name);
  return DestroTop::block(realname);
}

const desres::msys::Destro& desres::msys::Maeff::block(const std::string& name) const {
  std::string realname = adjustname(name);
  return DestroTop::block(realname);
}

void desres::msys::Maeff::write(std::ostream& os, int level) const {
  m_meta.write(os,level);
  size_t n = size();
  for(size_t i=1; i<=n; ++i) {
    this->block(i).write(os,level);
  }
}

size_t desres::msys::Maeff::footprint() const {
  return DestroTop::footprint() + m_meta.footprint();
}

void desres::msys::Maeff::touch(ZingPool& zpool) const {
  m_meta.touch(zpool);
  desres::msys::DestroTop::touch(zpool);
}
