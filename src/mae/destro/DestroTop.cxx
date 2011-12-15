/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"

// -----------------------------------------------
//  D E S T R O T O P
// -----------------------------------------------
desres::msys::DestroTop::DestroTop()
  : DestroBlock()
{
}

desres::msys::DestroTop::DestroTop(Tokenizer& tokenizer) /*GCOV-IGNORE*/
  : DestroBlock() /*GCOV-IGNORE*/
{
  fill_named(*this,tokenizer); /*GCOV-IGNORE*/
} /*GCOV-IGNORE*/

desres::msys::DestroTop::~DestroTop()
{
}



desres::msys::ZingPool& desres::msys::DestroTop::pool() {
  return m_pool;
}

const desres::msys::ZingPool& desres::msys::DestroTop::pool() const {
  return m_pool;
}

size_t desres::msys::DestroTop::footprint() const {
  return DestroBlock::footprint() + sizeof(m_pool) + m_pool.size();
}
