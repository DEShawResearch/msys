/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"

// -----------------------------------------------
// D E S T R O N A M E D B L O C K
// -----------------------------------------------
desres::msys::DestroNamedBlock::DestroNamedBlock(Destro* parent)
  : Destro(parent)
{
}

desres::msys::DestroNamedBlock::~DestroNamedBlock() {
}

std::string desres::msys::DestroNamedBlock::name() const {
  if (m_name.is_empty()) return "";
  return m_name.string(pool());
}

void desres::msys::DestroNamedBlock::name(const std::string& name) {
  Zing zname(name,pool());
  m_name = zname;
}

size_t desres::msys::DestroNamedBlock::footprint() const {
  return Destro::footprint() + sizeof(m_name);
}

void desres::msys::DestroNamedBlock::touch(ZingPool& zpool) const {
  m_name.touch(zpool);
}
