/* @COPYRIGHT@ */

#include "destro/Destro.hxx"

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

