/* @COPYRIGHT@ */

#include "destro/Destro.hxx"
#include "dessert/dessert.hpp"
#include "destro/FILEbuf.hxx"
#include <fstream>
#include <sstream>

/*!
 * Returns the current character in the file
 * @return current char
 */
char desres::Destro::Tokenizer::peek() {
  return m_c;
}

/*!
 * Read into the current character and update line
 * and point information.
 * @return (new) current char
 */
char desres::Destro::Tokenizer::read() {
  m_c = m_input->get(); 
  if (m_c == '\n') m_line++;
  m_point++;
  return m_c;
}

/*!
 * Test for end-of-file (stream,string, ...)
 * @return true iff at end of file
 */
bool desres::Destro::Tokenizer::eof() const {
  return m_input->eof();
}

/*!
 * A few Maestro tokens are just 1 char long.  We use
 * this test to short-circuit the tokenizer when at
 * [, ], {, or }.
 */
bool desres::Destro::Tokenizer::issingle(char c) {
  return (c == '[' || c == ']' || c == '{' || c == '}');
}

/*!
 * Build from an istream
 * @param input The stream to parse
 */
desres::Destro::Tokenizer::Tokenizer(std::istream& input)
  : m_c(0),
    m_buffer(NULL),
    m_input(&input),
    m_input_owned(false),
    m_token(NULL),
    m_isfresh(false),
    m_line(1),
    m_point(0),
    m_tokenline(1),
    m_tokenpoint(0)
{
  max_token_size = 16;
  m_token = (char *)malloc(max_token_size);
  prime_the_pump();
}

/*!
 * Make sure stream is ready to read and grab 1st token.
 * We require that the stream is not at EOF to start.
 */
void desres::Destro::Tokenizer::prime_the_pump() {
  if (m_input->eof()) throw dessert("cannot parse empty destro");
  read();
}

/*!
 * When we have a string, we build a temporary stringstream
 * and use that as our input.  The stringstream is deleted
 * in the destructor.
 * @param input The string to parse
 */
desres::Destro::Tokenizer::Tokenizer(const std::string& input)
  : m_c(0),
    m_buffer(NULL),
    m_input(new std::stringstream(input)),
    m_input_owned(true),
    m_token(NULL),
    m_isfresh(false),
    m_line(1),
    m_point(0),
    m_tokenline(1),
    m_tokenpoint(0)
{
  max_token_size = 16;
  m_token = (char *)malloc(max_token_size);
  prime_the_pump();
}

/*!
 * When we have a FILE*, we use the adaptor defined above to
 * create a C++ style input buffer object and then a new
 * istream (both deallocated on destruction).  File is not
 * closed.
 * @param input The FILE* to parse
 */
desres::Destro::Tokenizer::Tokenizer(FILE* input)
  : m_c(0),
    m_buffer(new FILEbuf(input)),
    m_input(new std::istream(m_buffer)),
    m_input_owned(true),
    m_token(NULL),
    m_isfresh(false),
    m_line(1),
    m_point(0),
    m_tokenline(1),
    m_tokenpoint(0)
{
  if (!input) throw dessert("file was not open");
  max_token_size = 16;
  m_token = (char *)malloc(max_token_size);
  prime_the_pump();
}

/*!
 * The destructor cleans up any heap allocated temporaries created
 * during construction.
 */
desres::Destro::Tokenizer::~Tokenizer() {
  if (m_token) free(m_token);
  delete m_buffer;
  if (m_input_owned) delete m_input;
}

/*!
 * This routine assembles a token character-by-character.
 * At its heart is a little DFA that looks at a character
 * to determine the next state.  For instance, the DFA
 * starts in a SKIPWHITE state and stays there unless
 * it finds a comment (#) or other character.  In the
 * INCOMMENT state, it looks for end-of-line before
 * returning to SKIPWHITE.  Similarly, all tokens are
 * defined using this one-character lookahead.
 * @return The current (possibly new) token 
 */
const char * desres::Destro::Tokenizer::token(bool ignore_single_character_tokens) {
  // -----------------------------------------------
  // Keep returning the same token until next()
  // is called.
  // -----------------------------------------------
  if (m_isfresh) return m_token;
  
  char * ptr = m_token;
  m_isfresh = true;
  
  unsigned state = SKIPWHITE;
  char c = peek();
  bool good = false;
  ssize_t diff;
  while(state != DONE && c >= 0) {
    // make sure we have space in m_token for 2 more characters
    if ((diff = ptr-m_token) >= max_token_size-1) {
      m_token = (char *)realloc( m_token, 2*max_token_size );
      ptr = m_token + diff;
      max_token_size *= 2;
    }
    switch(state) {
    case SKIPWHITE:
      // -----------------------------------------------
      // Skip whitespace and see if its a token or a
      // comment
      // -----------------------------------------------
      if (::isspace(c)) {
        c = read();
      } else if (c == '#') {
        state = INCOMMENT;
        c = read();
      } else {
        state = CHOOSEKIND;
      }
      break;
    case INCOMMENT:
      // -----------------------------------------------
      // On a comment, read until end of line
      // -----------------------------------------------
      if (c == '\n' || c == '#') state = SKIPWHITE;
      c = read();
      break;
    case CHOOSEKIND:
      // -----------------------------------------------
      // We []{} are single character tokens,
      // Strings start with "
      // Everything else starts with some other character
      // -----------------------------------------------
      if (issingle(c)) {
        state = SINGLECHAR;
      } else if (c == '"') {
        state = STARTSTRING;
      } else {
        state = STARTOTHER;
      }
      break;
    case SINGLECHAR:
      good = true;
      m_tokenline = m_line;
      m_tokenpoint = m_point;
      *ptr++ = c;
      *ptr++ = '\0';
      read();
      state = DONE;
      break;
    case STARTSTRING:
      good = true;
      m_tokenline = m_line;
      m_tokenpoint = m_point;
      *ptr++ = c;
      read(); // Skip opening quote
      c = peek();
      state = INSTRING;
      break;
    case INSTRING:
      if ( c == '"' ) {
        *ptr++ = c;
        *ptr++ = '\0';
        state = DONE;
      } else if ( c == '\\' ) {
        state = ESCAPE;
      } else {
        *ptr++ = c;
      }
      c = read();
      break;
    case ESCAPE:
      *ptr++ = c;
      state = INSTRING;
      c = read();
      break;
    case STARTOTHER:
      good = true;
      m_tokenline = m_line;
      m_tokenpoint = m_point;
      state = CONTINUEOTHER;
      break;
    case CONTINUEOTHER:
      if ( (!ignore_single_character_tokens && issingle(c)) || ::isspace(c) || c == '#' || c == '"' ) {
        *ptr++ = '\0';
        state = DONE;
      } else {
        *ptr++ = c;
        c = read();
      }
      break;
    }
  }
  
  // -----------------------------------------------
  // Maybe we just read trailing whitespace...
  // -----------------------------------------------
  if (!good) *m_token = '\0';
  
  return m_token;
}

/*!
 * Set state to read a new token on the next request.
 */
void desres::Destro::Tokenizer::next() {
  m_isfresh = false;
}

/*!
 * Line associated with current token
 * @return The line number
 */
unsigned desres::Destro::Tokenizer::line() const {
  return m_tokenline;
}

/*!
 * File offset associated with current token.
 * @return The offset
 */
size_t desres::Destro::Tokenizer::point() const {
  return m_tokenpoint;
}

/*!
 * Destro has occasion to look for the doc string
 * that follows an attribute discription.  This isn't
 * done with a strict parse (i.e. there is no
 * predict_optional_comment() routine) because that
 * would promote a comment to a full fledged token
 * class that could occur anywhere. So, if we think
 * a comment is there, we use this "mini-parser"
 * to grab it.
 * @return The comment (or empty string if no comment)
 */
std::string desres::Destro::Tokenizer::optional_comment() {
  std::string comment;
  char c;
  // -----------------------------------------------
  // Skip leading whitespace
  // -----------------------------------------------
  for(c=peek(); c >= 0 && ::isspace(c); c = peek()) {
    read();
  }
  
  // -----------------------------------------------
  // If we have a comment, read until end of line
  // or a #
  // -----------------------------------------------
  if (c == '#') {
    // Skip whitespace after the #
    for(c = read(); c >= 0;c = read()) {
      if (c == '\n') return comment; // Blank, naughty single #
      if (!::isspace(c)) break;
    }

    
    // Read body of comment until # or \n
    for(c = peek();
        c != '\n' && c != '#' && c >= 0;
        c = read()) {
      comment += c;
    }
    if (c >= 0) read(); // Consume final # or \n

    // Remove trailing whitespace
    while(comment.size() > 0 && ::isspace(comment[comment.size()-1])) {
      comment.erase(comment.size()-1);
    }
  }
  
  return comment; /*GCOV-IGNORE*/
}

/*!
 * The predictive, recursive-descent parsers I use do a lot
 * of "I expect the next token to look like this" calls.
 * (e.g. I expect a "{" here).  This simplifies the logic
 * for that.
 * @param match
 * @return The matching token body
 */
const char * desres::Destro::Tokenizer::predict(const char * match) {
  const char * tok = token();
  if (strcmp(match, "") && strcmp(tok, match)) {
    std::stringstream str;
    str << "Line " << line() << " predicted '" << match << "' have '"
        << (isprint(tok[0])?tok:"<unprintable>")
        << "'" << std::endl;
    throw dessert(str.str());
  }
  next();
  return tok;
}

/*!
 * Once I've bored down to a block or array block, I need
 * to grab a whitespace delimited token and have to ignore
 * []{} and the like.  This is particularly true for SMILES
 * names used in block (which embed [ and {) without quotes.
 * @return The matching token body
 */
const char * desres::Destro::Tokenizer::predict_value() {
  const char * tok = token(true);  // Ignore single char tokens here
  if ( (tok[0] == '\0') || (strcmp(tok,":::") == 0) || (strcmp(tok,"}") == 0)) {
    std::stringstream str;
    str << "Line " << line() << " predicted a value token, but I have a '"
        << (isprint(tok[0])?tok:"<unprintable>")
        << "'" << std::endl;
    throw dessert(str.str());
  }
  next();
  return tok;
}

/*!
 * Another common pattern is replication.  So, while (not_a("}")) { ... }
 * This function makes that easy.
 * @param match Token body to try to match
 * @return True on a match, False on EOF or a non-match
 */
bool desres::Destro::Tokenizer::not_a(const char * match) {
  const char * tok = token();
  if (!strcmp(tok, END_OF_FILE)) return false; // EOF always quits
  return strcmp(tok, match);
}

/*!
 * The special end of file token (string of length 1 containing a null byte)
 * Normal parsing will not create this character sequence, so it makes a
 * good special token.
 */
const char * desres::Destro::Tokenizer::END_OF_FILE = "";
