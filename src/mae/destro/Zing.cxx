/* @COPYRIGHT@ */

#include "destro/Zing.hxx"
#include "dessert/dessert.hpp"
#include <sstream>
#include <cstring>  // memset

#include <iostream> // Debug

#ifdef WIN32
#ifdef _WIN64
 typedef __int64 ssize_t;
#else
 typedef int ssize_t;
#endif

#endif


/* Frequency of letters in English text...
   http://www.askoxford.com/asktheexperts/faq/aboutwords/frequency

E  	1.1607%  	56.88
A 	8.4966% 	43.31 	
R 	7.5809% 	38.64 	
I 	7.5448% 	38.45 	
O 	7.1635% 	36.51 	
T 	6.9509% 	35.43 	
N 	6.6544% 	33.92 	
S 	5.7351% 	29.23 	
L 	5.4893% 	27.98 	
C 	4.5388% 	23.13 	
U 	3.6308% 	18.51 	
D 	3.3844% 	17.25 	
P 	3.1671% 	16.14 	
M  	3.0129%  	15.36
H 	3.0034% 	15.31
G 	2.4705% 	12.59
B 	2.0720% 	10.56
F 	1.8121% 	9.24
Y 	1.7779% 	9.06
W 	1.2899% 	6.57
K 	1.1016% 	5.61
V 	1.0074% 	5.13
X 	0.2902% 	1.48
Z 	0.2722% 	1.39
J 	0.1965% 	1.00
Q 	0.1962% 	(1)

EARIOTNSLCUDPMHGBFYWKVXZJQ
eariotnslcudpmhgbfywkvxzjq
*/

const uint32_t desres::msys::Zing::s_two_to_the_27th = 1 << 27;
const uint32_t desres::msys::Zing::s_two_to_the_26th = 1 << 26;
const uint32_t desres::msys::Zing::s_two_to_the_22nd = 1 << 22;

void desres::msys::ZingPool::check_assumptions() {
}

void desres::msys::Zing::check_assumptions() {
  if (sizeof(desres::msys::Zing) != 4) {
    throw dessert("Zing changed size",DESSERT_LOC);
  }

  if (desres::msys::Zing::s_alphabet_5bit.size() != 32) {
    throw dessert("bad alphabet",DESSERT_LOC);/*GCOV-IGNORE*/
  }
  if (desres::msys::Zing::s_alphabet_5bit[31] != 0) {
    throw dessert("bad alphabet doesn't end in NUL",DESSERT_LOC);/*GCOV-IGNORE*/
  }

  if (desres::msys::Zing::s_alphabet_amino.size() != 32) {
    throw dessert("bad alphabet",DESSERT_LOC);/*GCOV-IGNORE*/
  }
  if (desres::msys::Zing::s_alphabet_amino[31] != 0) {
    throw dessert("bad alphabet doesn't end in NUL",DESSERT_LOC);/*GCOV-IGNORE*/
  }

  if (desres::msys::Zing::s_alphabet_4bit.size() != 16) {
    throw dessert("bad alphabet",DESSERT_LOC);/*GCOV-IGNORE*/
  }
  if (desres::msys::Zing::s_alphabet_4bit[15] != 0) {
    throw dessert("bad alphabet doesn't end in NUL",DESSERT_LOC);/*GCOV-IGNORE*/
  }

  if (desres::msys::Zing::s_alphabet_bcd.size() != 16) {
    throw dessert("bad alphabet",DESSERT_LOC);/*GCOV-IGNORE*/
  }
  if (desres::msys::Zing::s_alphabet_bcd[15] != 0) {
    throw dessert("bad alphabet doesn't end in NUL",DESSERT_LOC);/*GCOV-IGNORE*/
  }
}

const size_t desres::msys::ZingPool::s_size = 16;

// Fixed point
// 0ooo oooo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_float    = 0x1 << 31;
const uint32_t desres::msys::Zing::control_float = 0x0 << 31;

// 6 chars (5-bit alphabet)
// 10oo oooo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_5bit     = 0x3 << 30;
const uint32_t desres::msys::Zing::control_5bit  = 0x2 << 30; 

// 7 chars (BCD 4-bit alphabet)
// 1100 oooo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_bcd      = 0xf << 28;
const uint32_t desres::msys::Zing::control_bcd   = 0xc << 28;

// 4 chars (7-bit ASCII)
// 1101 oooo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_7bit     = 0xf << 28;
const uint32_t desres::msys::Zing::control_7bit  = 0xd << 28; 

// Pool (1..poolsize)
// 1110 oooo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_size     = 0xf << 28;
const uint32_t desres::msys::Zing::control_size  = 0xe << 28;

// Pool (overflow)
// 1111 0ooo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_any      = 0x1f << 27;
const uint32_t desres::msys::Zing::control_any   = 0x1e << 27;

// 6 chars (specialized 4-bit)
// 1111 10oo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_4bit     = 0x3f << 26;
const uint32_t desres::msys::Zing::control_4bit  = 0x3e << 26; 

// 5 chars (specialized 5 bit for amino acids)
// 1111 11oo oooo oooo oooo oooo oooo oooo
const uint32_t desres::msys::Zing::mask_amino    = 0x3f << 26;
const uint32_t desres::msys::Zing::control_amino = 0x3f << 26; 

desres::msys::ZingPool::ZingPool()
  : m_inserts_allowed(true),
    m_pool(s_size),
    m_pool_active(s_size)
{
}

desres::msys::ZingPool::~ZingPool()
{
}

void desres::msys::ZingPool::clearbits() {
  m_inserts_allowed = false;

  for(std::vector< std::vector<unsigned char> >::iterator p = m_pool_active.begin(),
        en = m_pool_active.end(); p != en; ++p) {
    std::vector<unsigned char>& v = *p;
    size_t n = v.size();
    if (n > 0) memset(&v[0],0,n);
  }
  
  size_t n = m_overflow_active.size();
  if (n > 0) memset(&m_overflow_active[0],0,n);
}

void desres::msys::ZingPool::touch(size_t n, size_t length) {
  if (m_inserts_allowed) throw dessert("pool is in insertion mode.  Must call clearbits() first",DESSERT_LOC);

  if (length == 0) {
    // Nothing to do
  } else if (length <= s_size) {
    m_pool_active.at(length-1).at(n) = 1;
  } else {
    m_overflow_active[n] = 1;
  }
}

void desres::msys::ZingPool::finish_compaction() {
  m_inserts_allowed = true;
}


size_t desres::msys::ZingPool::size() const {
  size_t size = 0;
  for(std::vector< std::vector<char> >::const_iterator p = m_pool.begin(),
        en=m_pool.end(); p != en; ++p) {
    size += (*p).size();
  }
  for(std::vector< std::vector<unsigned char> >::const_iterator p = m_pool_active.begin(),
        en=m_pool_active.end(); p != en; ++p) {
    size += (*p).size();
  }

  for(std::vector< std::string >::const_iterator p = m_overflow.begin(),
        en=m_overflow.end(); p != en; ++p) {
    size += (*p).size();
  }
  size += m_overflow_active.size();
  return size;
}

size_t desres::msys::ZingPool::insert(const std::string& string) {
  if (!m_inserts_allowed) {
    throw dessert("The pool is locked.  Inserts not allowed",
                           DESSERT_LOC);
  }

  size_t length = string.size();

  // Zero length strings are always encoded as 0
  if (length == 0) {
    return 0;
  }

  // If we have a bin defined, we look there for a duplicate
  else if (length <= s_size) {
    // Note length-1 to get the i'th pool (no need to store
    // zero length strings)
    std::vector<char>& storage = m_pool[length-1];
    std::vector<unsigned char>& active = m_pool_active[length-1];

    // Look to see if we have this string already
    unsigned n = storage.size()/length;
    ssize_t empty_slot = -1;
    //for(size_t i=0; i<storage.size(); i += length) {
    for(size_t i=0; i<n; ++i) {
      // Don't check empty slots
      if (!active[i]) {
        empty_slot = i;
        continue;
      }

      bool match = true;
      for(size_t j=0; j<length; ++j) {
        if (string[j] != storage[i*length+j]) {
          match = false;
          break;
        }
      }
      if (match) return i;
    }

    // If we didn't match it, we can copy it into an empty
    // slot (if we have one).
    if ( empty_slot >= 0 ) {
      for(size_t j=0; j<length; ++j) {
        storage[empty_slot*length+j] = string[j];
      }
      active[empty_slot] = 1;
      return empty_slot;
    }

    // Have to append this string (and mark active)
    for(size_t j=0; j<length; ++j) {
      storage.push_back(string[j]);
    }
    active.push_back(1);
    return (storage.size()-length)/length;
  }
   
  // It is possible that we already have the string in
  // our storage.  A hash here would be faster, but that
  // would take more space and that is what we wish to
  // optimize.
  ssize_t empty_slot = -1;
  size_t n = m_overflow.size();
  for(unsigned i=0;i<n;++i) {
    if (!m_overflow_active[i]) {
      empty_slot = i;
      continue;
    }
    if (m_overflow[i] == string) return i;
  }

  // If we have an empty slot, use it...
  if ( empty_slot >= 0) {
    m_overflow[empty_slot] = string;
    m_overflow_active[empty_slot] = 1;
    return empty_slot;
  }

  // OK, this is a new string.  Push it onto the end.
  m_overflow.push_back(string);
  m_overflow_active.push_back(1);
  return m_overflow.size()-1;
}

std::string desres::msys::ZingPool::retrieve(size_t n,size_t length) const {
  if (length == 0) {
    return "";
  } else if (length <= s_size) {
    size_t start = n*length;
    size_t end = start+length;
    
    // Deal with "lost" strings 
    const std::vector<char>& storage = m_pool[length-1];
    if (end > storage.size()) throw dessert("string lost");

    if (!m_pool_active[length-1].at(n)) {
      throw dessert("The zing was marked dead",DESSERT_LOC);
    }

    // Assumes vector<> contiguous storage
    return std::string(&storage[start],length);
  }
   
  if (n >= m_overflow.size()) throw dessert("string lost");
  if (!m_overflow_active.at(n)) {
    throw dessert("The zing was marked dead",DESSERT_LOC); /*GCOV-IGNORE*/
  }
  return m_overflow[n];
}

void desres::msys::ZingPool::dump() const {
  std::cerr << "DUMP POOL =============================" << std::endl;
  for(unsigned ii=0; ii<s_size; ++ii) {
    size_t n = m_pool[ii].size()/(ii+1);
    std::cerr << ii+1 << "->" << n << std::endl;
    for(unsigned jj=0;jj<n;++jj) {
      if (m_pool_active[ii][jj]) {
        for(unsigned kk=0; kk < (ii+1); ++kk) std::cerr << ' ';
      } else {
        for(unsigned kk=0; kk < (ii+1); ++kk) std::cerr << 'X';
      }
    }
    std::cerr << std::endl;
    if (m_pool[ii].size() > 0) {
      std::string x(&m_pool[ii][0],m_pool[ii].size());
      std::cerr << x << std::endl;
    }
  }
  for(unsigned ii=0; ii<m_overflow.size(); ++ii) {
    if (m_overflow_active[ii]) {
      std::cerr << "o";
    } else {
      std::cerr << "x";
    }
    std::cerr << ii << " = " << m_overflow[ii] << std::endl;
  }
  std::cerr << "=======================================" << std::endl;
}

uint32_t desres::msys::Zing::add_control(uint32_t x,uint32_t mask,uint32_t control) {
  if (x & mask) throw dessert("overflows into mask",DESSERT_LOC);
#if PEDANTIC
  assert( (~mask) & control == 0 );
#endif
  return x | control;
}

const uint32_t desres::msys::Zing::encode_float(const std::string& string) {
  size_t length = string.size();
  
  // Let BCD handle short strings...
  if (length < 8) return 0;

  // We have problems with numbers like 00.1234567 and
  // 01.456789 because of the integer encoding scheme
  // we use.  Just punt and let "any" catch them.
  if (string[0] == '0' && string[1] != '.') return 0;

  // Optional lead sign
  unsigned i = 0;
  int sign = 0;
  if (string[0] == '-') {
    ++i;
    sign = 1;
  }

  uint32_t ival = 0;
  int eee = -1;
  for(;i<length;++i) {
    char c = string[i];

    if (c == '.') {
      if (eee >= 0) return 0; // Multiple . in string
      eee = 0;
      continue;
    }
    if (eee >= 0) {
      if (++eee >= 8) return 0; // Exponent out of range
    }

    if (c < '0' || c > '9') return 0;
    uint32_t digit = c-'0';
    ival = ival*10 + digit;
    if (ival >= s_two_to_the_27th) return 0;
  }
  if (eee <= 0) return 0;  // No digits after the decimal point

  uint32_t encode = (((ival << 1) | sign) << 3) | (eee-1);
  uint32_t result = add_control(encode,mask_float,control_float);
  return result;
}
const std::string desres::msys::Zing::decode_float(uint32_t x) {
#ifdef PEDANTIC
  uint32_t control = x & mask_float;
  if (control != control_float) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif

  x = x & ~mask_float;

  uint32_t eee = (x & 0x7) + 1;
  x = x >> 3;

  uint32_t sign = x & 0x1;
  x = x >> 1;

  uint32_t ival = x;

  std::ostringstream ss;
  ss << ival;
  std::string sval = ss.str();

  // Pat with zeros to get to exponent length
  while(sval.size() < eee) sval = '0'+sval;

  // Stick the . back between prefix and suffix
  sval.insert(sval.size()-eee,".");

  // Canonicalize the number for fractionals
  if ( sval[0] == '.' ) sval.insert(0,"0");

  // Add the sign back in
  if (sign) sval.insert(0,"-");

  return sval;
}

const std::string
desres::msys::Zing::s_alphabet_5bit("eariotnslcudpmhgbfywkvxzjq_012\"\0",32);

const std::string
desres::msys::Zing::s_alphabet_amino("MVCINUOGHTYEPARSLDW 01234567890\0",32);

const std::string
desres::msys::Zing::s_alphabet_4bit("0123456789 CONH\0",16);

const uint32_t desres::msys::Zing::encode_nbit(const std::string& string, unsigned nbits, unsigned nchars, uint32_t mask, uint32_t control, const std::string& alphabet) {
  size_t length = string.size();
  if (length == 0 || length > nchars) return 0;

  // Fill with NUL characters as needed
  uint32_t encode = 0;
  unsigned nul = (1<<nbits)-1;
  for(unsigned i=length; i<nchars; ++i) {
    encode = (encode << nbits) | nul;
  }

  // We merrily encode characters, but if we find an
  // embedded null or non-alphabet character, we
  // simply return '0' to show we can't encode it
  // this way
  for(std::string::const_reverse_iterator p = string.rbegin(),
        en=string.rend(); p != en; ++p) {
    char c = *p;

    std::string::size_type pos = alphabet.find(c);
    if (pos == std::string::npos) return 0; // Not in alphabet

    encode = (encode << nbits) | pos;
  }

  uint32_t result = add_control(encode,mask,control);
  return result;
}

// 00 ..... ..... ..... ..... .....
const uint32_t desres::msys::Zing::encode_5bit(const std::string& string) {
  return encode_nbit(string,5,6,mask_5bit,control_5bit,s_alphabet_5bit);
#if 0
  size_t length = string.size();
  if (length == 0 || length > 6) return 0;

  // Fill with NUL characters as needed
  uint32_t encode = 0;
  for(unsigned i=length; i<6; ++i) {
    encode = (encode << 5) | 31;
  }

  // We merrily encode characters, but if we find an
  // embedded null or non-alphabet character, we
  // simply return '0' to show we can't encode it
  // this way
  for(std::string::const_reverse_iterator p = string.rbegin(),
        en=string.rend(); p != en; ++p) {
    char c = *p;

    std::string::size_type pos = s_alphabet_5bit.find(c);
    if (pos == std::string::npos) return 0; // Not in alphabet

    encode = (encode << 5) | pos;
  }

  uint32_t result = add_control(encode,mask_5bit,control_5bit);
  return result;
#endif
}

const std::string desres::msys::Zing::decode_5bit(uint32_t x) {
#ifdef PEDANTIC
  uint32_t control = x & mask_5bit;
  if (control != control_5bit) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif

  x = x & ~mask_5bit;

  char result[7];
  result[0] = s_alphabet_5bit[x & 0x1f]; x = x >> 5;
  result[1] = s_alphabet_5bit[x & 0x1f]; x = x >> 5;
  result[2] = s_alphabet_5bit[x & 0x1f]; x = x >> 5;
  result[3] = s_alphabet_5bit[x & 0x1f]; x = x >> 5;
  result[4] = s_alphabet_5bit[x & 0x1f]; x = x >> 5;
  result[5] = s_alphabet_5bit[x & 0x1f]; x = x >> 5;
  result[6] = 0;

  return result;
}

const uint32_t desres::msys::Zing::encode_4bit(const std::string& string) {
  return encode_nbit(string,4,6,mask_4bit,control_4bit,s_alphabet_4bit);
}

const std::string desres::msys::Zing::decode_4bit(uint32_t x) {
#ifdef PEDANTIC
  uint32_t control = x & mask_4bit;
  if (control != control_4bit) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif

  x = x & ~mask_4bit;

  char result[7];
  result[0] = s_alphabet_4bit[x & 0xf]; x = x >> 4;
  result[1] = s_alphabet_4bit[x & 0xf]; x = x >> 4;
  result[2] = s_alphabet_4bit[x & 0xf]; x = x >> 4;
  result[3] = s_alphabet_4bit[x & 0xf]; x = x >> 4;
  result[4] = s_alphabet_4bit[x & 0xf]; x = x >> 4;
  result[5] = s_alphabet_4bit[x & 0xf]; x = x >> 4;
  result[6] = 0;

  return result;
}

const uint32_t desres::msys::Zing::encode_amino(const std::string& string) {
  return encode_nbit(string,5,5,mask_amino,control_amino,s_alphabet_amino);
}

const std::string desres::msys::Zing::decode_amino(uint32_t x) {
#ifdef PEDANTIC
  uint32_t control = x & mask_amino;
  if (control != control_amino) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif

  x = x & ~mask_amino;

  char result[6];
  result[0] = s_alphabet_amino[x & 0x1f]; x = x >> 5;
  result[1] = s_alphabet_amino[x & 0x1f]; x = x >> 5;
  result[2] = s_alphabet_amino[x & 0x1f]; x = x >> 5;
  result[3] = s_alphabet_amino[x & 0x1f]; x = x >> 5;
  result[4] = s_alphabet_amino[x & 0x1f]; x = x >> 5;
  result[5] = 0;

  return result;
}

const uint32_t desres::msys::Zing::encode_any(const std::string& string,ZingPool& pool) {
  size_t length = string.size();
  if (length <= pool.s_size) return 0; // Use length string for small ones

  uint32_t offset = pool.insert(string);
  if (offset > s_two_to_the_26th) return 0;
  uint32_t result = add_control(offset,mask_any,control_any);

  return result;
}

const std::string desres::msys::Zing::decode_any(uint32_t x,const ZingPool& pool) {
#ifdef PEDANTIC
  uint32_t control = x & mask_any;
  if (control != control_any) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif
  x = x & ~mask_any;

  return pool.retrieve(x,pool.s_size+1);
}

// 0111 oooooooooooooooooooooooo ssss
const uint32_t desres::msys::Zing::encode_size(const std::string& string,ZingPool& pool) {
  size_t length = string.size();
  if (length == 0) return 0;  // Can't encode length 0 string this way
  if (length > pool.s_size) return 0; // Use "any" for long strings

  uint32_t offset = pool.insert(string);
  if (offset > s_two_to_the_22nd) return 0;
  uint32_t result = add_control((offset <<4) | (length-1) ,mask_size,control_size);

  return result;
}
const std::string desres::msys::Zing::decode_size(uint32_t x,const ZingPool& pool) {
#ifdef PEDANTIC
  uint32_t control = x & mask_size;
  if (control != control_size) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif
  x = x & ~mask_size;

  size_t length = x & 0xf;
  x = x >> 4;

  return pool.retrieve(x,length+1);
}

const std::string desres::msys::Zing::s_alphabet_bcd("0123456789.+-eE\0",16);

// 0100 .... .... .... .... .... .... ....
const uint32_t desres::msys::Zing::encode_bcd(const std::string& string) {
  size_t length = string.size();
  if (length == 0 || length > 7) return 0;

  // Fill with NUL characters as needed
  uint32_t encode = 0;
  for(unsigned i=length; i<7; ++i) {
    encode = (encode << 4) | 15;
  }

  // We merrily encode characters, but if we find an
  // embedded null or non-alphabet character, we
  // simply return '0' to show we can't encode it
  // this way
  for(std::string::const_reverse_iterator p = string.rbegin(),
        en=string.rend(); p != en; ++p) {
    char c = *p;

    std::string::size_type pos = s_alphabet_bcd.find(c);
    if (pos == std::string::npos) return 0; // Not in alphabet

    encode = (encode << 4) | pos;
  }

  uint32_t result = add_control(encode,mask_bcd,control_bcd);
  return result;
}

const std::string desres::msys::Zing::decode_bcd(uint32_t x) {
#ifdef PEDANTIC
  uint32_t control = x & mask_bcd;
  if (control != control_bcd) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif

  x = x & ~mask_bcd;

  char result[8];
  result[0] = s_alphabet_bcd[x & 0xf]; x = x >> 4;
  result[1] = s_alphabet_bcd[x & 0xf]; x = x >> 4;
  result[2] = s_alphabet_bcd[x & 0xf]; x = x >> 4;
  result[3] = s_alphabet_bcd[x & 0xf]; x = x >> 4;
  result[4] = s_alphabet_bcd[x & 0xf]; x = x >> 4;
  result[5] = s_alphabet_bcd[x & 0xf]; x = x >> 4;
  result[6] = s_alphabet_bcd[x & 0xf];
  result[7] = 0;

  return result;
}

// 0101 ------- ------- ------- -------
const uint32_t desres::msys::Zing::encode_7bit(const std::string& string) {
  size_t length = string.size();

  // I only have room for four 7-bit letters
  if (length > 4) return 0;

  // Fill with NUL characters as needed
  uint32_t encode = 0;
  for(unsigned i=length; i<4; ++i) {
    encode = (encode << 7);
  }

  // We fill the string backwards to speed unpacking
  for(std::string::const_reverse_iterator p = string.rbegin(),
        en=string.rend(); p != en; ++p) {
    char c = *p;
    if (c & 0x80) return 0; // Isn't an 7-bit ascii character
    encode = (encode << 7) | c;
  }

  uint32_t result = add_control(encode,mask_7bit,control_7bit);

  return result;
}

const std::string desres::msys::Zing::decode_7bit(uint32_t x) {
#ifdef PEDANTIC
  uint32_t control = x & mask_7bit;
  if (control != control_7bit) {
    throw dessert("illegal decode",DESSERT_LOC);
  }
#endif
  x = x & ~mask_7bit;

  char result[5];
  result[0] = x & 0x7f; x = x >> 7;
  result[1] = x & 0x7f; x = x >> 7;
  result[2] = x & 0x7f; x = x >> 7;
  result[3] = x & 0x7f; x = x >> 7;
  result[4] = 0;

  return result;
}


desres::msys::Zing::Zing()
  : m_bits(0) // The default gets EMPTY encoding
{
}

desres::msys::Zing::Zing(const Zing& other)
  : m_bits(other.m_bits)
{
}

desres::msys::Zing& desres::msys::Zing::operator=(const desres::msys::Zing& other)
{
  m_bits = other.m_bits;
  return *this;
}

bool desres::msys::Zing::operator==(const desres::msys::Zing& other) const {
  return m_bits == other.m_bits;
}

bool desres::msys::Zing::operator!=(const desres::msys::Zing& other) const {
  return m_bits != other.m_bits;
}

desres::msys::Zing::Zing(const std::string& string, desres::msys::ZingPool& pool)
  : m_bits(0)
{
  // The empty string "" gets standard 7BIT encoding

  if (!((m_bits = encode_7bit(string)) ||
        (m_bits = encode_bcd(string)) ||
        (m_bits = encode_5bit(string)) ||
        (m_bits = encode_4bit(string)) ||
        (m_bits = encode_amino(string)) ||
        (m_bits = encode_float(string)) ||
        (m_bits = encode_size(string,pool)) ||
        (m_bits = encode_any(string,pool)))) {
    // This could only happen if we exhausted storage for
    // 'any' strings or if we had an encode error
    throw dessert("out of storage",DESSERT_LOC);/*GCOV-IGNORE*/
  }
}


desres::msys::Zing::~Zing() {
}

std::string desres::msys::Zing::string(const desres::msys::ZingPool& pool) const {
  uint32_t bits = m_bits;

  if (bits == 0) throw dessert("Uninitialized zing",DESSERT_LOC);

  if ((bits & mask_7bit) == control_7bit) return decode_7bit(bits);
  if ((bits & mask_bcd) == control_bcd) return decode_bcd(bits);
  if ((bits & mask_5bit) == control_5bit) return decode_5bit(bits);
  if ((bits & mask_4bit) == control_4bit) return decode_4bit(bits);
  if ((bits & mask_amino) == control_amino) return decode_amino(bits);
  if ((bits & mask_float) == control_float) return decode_float(bits);
  if ((bits & mask_size) == control_size) return decode_size(bits,pool);
  if ((bits & mask_any) == control_any) return decode_any(bits,pool);

  throw dessert("Cannot decipher control bits",DESSERT_LOC);/*GCOV-IGNORE*/
}

uint32_t desres::msys::Zing::integer() const {
  return m_bits;
}

desres::msys::Zing::category_t desres::msys::Zing::category() const {
  uint32_t bits = m_bits;

  if (bits == 0) return ZING_EMPTY;
  if ((bits & mask_7bit) == control_7bit) return ZING_7BIT;
  if ((bits & mask_bcd) == control_bcd) return ZING_BCD;
  if ((bits & mask_5bit) == control_5bit) return ZING_5BIT;
  if ((bits & mask_4bit) == control_4bit) return ZING_4BIT;
  if ((bits & mask_amino) == control_amino) return ZING_AMINO;
  if ((bits & mask_float) == control_float) return ZING_FLOAT;
  if ((bits & mask_size) == control_size) return ZING_SIZE;
  if ((bits & mask_any) == control_any) return ZING_ANY;
  
  throw dessert("Invalid category",DESSERT_LOC);/*GCOV-IGNORE*/
}


bool desres::msys::Zing::is_empty() const {
  return m_bits == 0;
}

void desres::msys::Zing::touch(desres::msys::ZingPool& pool) const {
  uint32_t bits = m_bits;

  if (bits == 0) return;
  if ((bits & mask_7bit) == control_7bit) return;
  if ((bits & mask_bcd) == control_bcd) return;
  if ((bits & mask_5bit) == control_5bit) return;
  if ((bits & mask_4bit) == control_4bit) return;
  if ((bits & mask_amino) == control_amino) return;
  if ((bits & mask_float) == control_float) return;

  if ((bits & mask_size) == control_size) {
    bits = bits & ~mask_size;

    size_t length = bits & 0xf;
    bits = bits >> 4;

    pool.touch(bits,length+1);
    return;
  }

  if ((bits & mask_any) == control_any) {
     bits = bits & ~mask_any;
     pool.touch(bits,pool.s_size+1);
     return;
  }
  
  throw dessert("Invalid category",DESSERT_LOC);/*GCOV-IGNORE*/
}
