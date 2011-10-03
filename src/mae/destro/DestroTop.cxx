/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"

// -----------------------------------------------
//  D E S T R O T O P
// -----------------------------------------------
desres::DestroTop::DestroTop()
  : DestroBlock()
{
}

desres::DestroTop::DestroTop(Tokenizer& tokenizer) /*GCOV-IGNORE*/
  : DestroBlock() /*GCOV-IGNORE*/
{
  fill_named(*this,tokenizer); /*GCOV-IGNORE*/
} /*GCOV-IGNORE*/

desres::DestroTop::~DestroTop()
{
}



desres::ZingPool& desres::DestroTop::pool() {
  return m_pool;
}

const desres::ZingPool& desres::DestroTop::pool() const {
  return m_pool;
}

size_t desres::DestroTop::footprint() const {
  return DestroBlock::footprint() + sizeof(m_pool) + m_pool.size();
}
