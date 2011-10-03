/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"

// -----------------------------------------------
// D E S T R O N A M E D B L O C K
// -----------------------------------------------
desres::DestroNamedBlock::DestroNamedBlock(Destro* parent)
  : Destro(parent)
{
}

desres::DestroNamedBlock::~DestroNamedBlock() {
}

std::string desres::DestroNamedBlock::name() const {
  if (m_name.is_empty()) return "";
  return m_name.string(pool());
}

void desres::DestroNamedBlock::name(const std::string& name) {
  Zing zname(name,pool());
  m_name = zname;
}

size_t desres::DestroNamedBlock::footprint() const {
  return Destro::footprint() + sizeof(m_name);
}

void desres::DestroNamedBlock::touch(ZingPool& zpool) const {
  m_name.touch(zpool);
}
