/* @COPYRIGHT@ */

#include "mae.hxx"
#include "../types.hxx"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifdef WIN32
#ifdef _WIN64
 typedef __int64 ssize_t;
#else
 typedef int ssize_t;
#endif

#endif



using desres::msys::fastjson::Json;

/*!
 * \brief Takes a stream and returns maestro tokens.
 * This tokenizer is built on streams and uses a small, tight
 * finite state automata to construct a token
 */

namespace desres { namespace msys { namespace mae {

    struct tokenizer {

        char buf[256];
  
        /*! \brief The current character */
        char m_c;
  
        std::streamsize bufpos;
        std::streamsize bufsize;
  
        /*! \brief the stream for the file we're parsing */
        std::istream * m_input;
  
        /*! \brief The current token */
        char * m_token;
  
        /*! \brief number of malloc'ed bytes in m_token */
        ssize_t max_token_size;
  
        /*! \brief True iff the token is already read */
        int m_isfresh;
  
        /*! \brief Current line in file */
        unsigned m_line;

        /* current position in file */
        std::streamsize m_offset;
  
        /*! \brief Line where token starts */
        unsigned m_tokenline;
    };
}}}

using namespace desres::msys;
using desres::msys::mae::tokenizer;

/*! \brief Get current character */
/*!
* Returns the current character in the file
* @return current char
*/
static inline char tokenizer_peek(const tokenizer * tk) { return tk->m_c; }

/*! \brief Read a new character */
/*!
* Read into the current character and update line
* and point information.
* @return (new) current char
*/
static inline char tokenizer_read(tokenizer * tk) {
  if (tk->bufpos==tk->bufsize) {
      tk->m_offset += tk->bufsize;
      //printf("tokenizer_read %lu\n", sizeof(tk->buf));
      tk->m_input->read(tk->buf, sizeof(tk->buf));
      tk->bufsize = tk->m_input->gcount();
      tk->bufpos = 0;
      //printf("  gcount %ld\n", tk->bufsize);
  } 
  tk->m_c = tk->bufsize ? tk->buf[tk->bufpos++] : -1;
  if (tk->m_c == '\n') tk->m_line++;
  return tk->m_c;
}

/*!
 * Build from an istream
 * @param input The stream to parse
 */

static void tokenizer_init( tokenizer * tk, std::istream& input );

/*!
 * The destructor cleans up any heap allocated temporaries created
 * during construction.
 */
static void tokenizer_release( tokenizer * tk);

/*!
 * Set state to read a new token on the next request.
 */
static inline void tokenizer_next(tokenizer * tk) {
  tk->m_isfresh = 0;
}

/*!
 * Line associated with current token
 * @return The line number
 */
static inline unsigned tokenizer_line(const tokenizer * tk) {
  return tk->m_tokenline;
}


/*!
 * The special end of file token (string of length 1 containing a null byte)
 * Normal parsing will not create this character sequence, so it makes a
 * good special token.
 */
static const char * END_OF_FILE = "";

static inline int issingle(char c) {
    return (c == '[' || c == ']' || c == '{' || c == '}');
}


/*! \brief Actions for the DFA token builder
 *
 */
enum ACTION {
  DONE = 0,
  SKIPWHITE,  /* 1  */
  INCOMMENT,  /*  2 */
  CHOOSEKIND,  /* 3 */
  SINGLECHAR,  /* 4 */
  STARTSTRING,  /* 5  */
  INSTRING,  /* 6 */
  ESCAPE,  /* 7  */
  STARTOTHER,  /* 8  */
  CONTINUEOTHER  /* 9 */
};

void tokenizer_init( tokenizer * tk, std::istream& input ) {
    memset(tk,0,sizeof(*tk));
    tk->m_input = &input;
    tk->m_line = 1;
    tk->m_tokenline = 1;
    tk->max_token_size = 16;
    tk->m_token = (char *)malloc(tk->max_token_size);

    /* grab 1st token */
    tokenizer_read(tk);
}

/*!
 * The destructor cleans up any heap allocated temporaries created
 * during construction.
 */
void tokenizer_release( tokenizer * tk) {
    if (tk->m_token) free(tk->m_token);
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
static inline const char * tokenizer_token(tokenizer * tk, int ignore_single) {
  /* -----------------------------------------------
  // Keep returning the same token until next()
  // is called.
  // -----------------------------------------------
  */

  char * ptr;
  char c;
  int good;
  ssize_t diff;
  unsigned state = SKIPWHITE;
  if (tk->m_isfresh) return tk->m_token;
  
  /* -----------------------------------------------
  // End of file simply returns an empty string
  // -----------------------------------------------
  */

  /* begin at start of token space */
  ptr = tk->m_token;
  tk->m_isfresh = 1;
  
  c = tokenizer_peek(tk);
  good = 0;
  while(state != DONE && c != 0 && c != -1) {
    /* make sure we have space in m_token for 2 more characters */
    if ((diff = ptr-tk->m_token) >= tk->max_token_size-1) {
      tk->m_token = (char *)realloc( tk->m_token, 2*tk->max_token_size );
      ptr = tk->m_token + diff;
      tk->max_token_size *= 2;
    }
    switch(state) {
    case SKIPWHITE:
      /* -----------------------------------------------
      // Skip whitespace and see if its a token or a
      // comment
      // ----------------------------------------------- */
      if (isspace(c)) {
        c = tokenizer_read(tk);
      } else if (c == '#') {
        state = INCOMMENT;
        c = tokenizer_read(tk);
      } else {
        state = CHOOSEKIND;
      }
      break;
    case INCOMMENT:
      /* -----------------------------------------------
      // On a comment, read until end of line
      // ----------------------------------------------- */
      if (c == '\n' || c == '#') state = SKIPWHITE;
      c = tokenizer_read(tk);
      break;
    case CHOOSEKIND:
      /* -----------------------------------------------
      // We []{} are single character tokens,
      // Strings start with "
      // Everything else starts with some other character
      // ----------------------------------------------- */
      if (issingle(c)) {
        if (ignore_single)
          state = STARTOTHER;
        else  
          state = SINGLECHAR;
      } else if (c == '"') {
        state = STARTSTRING;
      } else {
        state = STARTOTHER;
      }
      break;
    case SINGLECHAR:
      good = 1;
      tk->m_tokenline = tk->m_line;
      *ptr++ = c;
      *ptr++ = '\0';
      tokenizer_read(tk);
      state = DONE;
      break;
    case STARTSTRING:
      good = 1;
      tk->m_tokenline = tk->m_line;
      *ptr++ = c;
      tokenizer_read(tk); /* Skip opening quote */
      c = tokenizer_peek(tk);
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
      c = tokenizer_read(tk);
      break;
    case ESCAPE:
      *ptr++ = c;
      state = INSTRING;
      c = tokenizer_read(tk);
      break;
    case STARTOTHER:
      good = 1;
      tk->m_tokenline = tk->m_line;
      state = CONTINUEOTHER;
      break;
    case CONTINUEOTHER:
      //if ( (!ignore_single && issingle(c)) || isspace(c) || c == '#' || c == '"' ) {
      if (ignore_single) {
        if (isspace(c) || c == '\n') {
          *ptr++ = '\0';
          state = DONE;
        } else {
          *ptr++ = c;
          c = tokenizer_read(tk);
        }
      } else {
        if (issingle(c) || isspace(c) || c == '"') {
          *ptr++ = '\0';
          state = DONE;
          /* allow hashes in tokens; comments require preceding space */
        } else if (c=='#' && isspace(ptr[-1])) {
          *ptr++ = '\0';
          state = DONE;
        } else {
          *ptr++ = c;
          c = tokenizer_read(tk);
        }
      }
      break;
    }
  }
  
  /* -----------------------------------------------
  // Maybe we just read trailing whitespace...
  // ----------------------------------------------- */
  if (!good) *tk->m_token = '\0';
  
  return tk->m_token;
}

#define MAE_ERROR3(fmt, arg1, arg2, arg3) do { \
    char buf[4096]; \
    sprintf(buf, fmt, arg1, arg2, arg3); \
    throw std::runtime_error(buf); \
} while (0)

#define MAE_ERROR2(fmt, arg1, arg2) do { \
    char buf[4096]; \
    sprintf(buf, fmt, arg1, arg2); \
    throw std::runtime_error(buf); \
} while (0)

#define MAE_ERROR1(fmt, arg1) do { \
    char buf[4096]; \
    sprintf(buf, fmt, arg1); \
    throw std::runtime_error(buf); \
} while (0)


/*!
 * The predictive, recursive-descent parsers I use do a lot
 * of "I expect the next token to look like this" calls.
 * (e.g. I expect a "{" here).  This simplifies the logic
 * for that.
 * @param match
 * @return The matching token body
 */
static inline const char *
tokenizer_predict(tokenizer * tk, const char * match) {
  const char * tok = tokenizer_token(tk,0);
  if (strcmp(match, "") && strcmp(tok, match)) {
      MAE_ERROR3("Line %d predicted '%s', have '%s'",
                  tokenizer_line(tk) ,
                  match ,
                  (isprint(tok[0])?tok:"<unprintable>"));
  }
  tokenizer_next(tk);
  return tok;
}

static inline const char *
tokenizer_predict_value(tokenizer * tk) {
  const char * tok = tokenizer_token(tk,1);
  if(tok[0]=='\0') {
      MSYS_FAIL("Premature end of file at line " << tokenizer_line(tk));
  } else if (!strcmp(tok,":::") || !strcmp(tok,"}")) {
      MSYS_FAIL("Line " << tokenizer_line(tk) << " unexpected characters '" << tok << "'");
  } else {
    tokenizer_next(tk);
  }
  return tok;
}

/*!
 * Another common pattern is replication.  So, while (not_a("}")) { ... }
 * This function makes that easy.
 * @param match Token body to try to match
 * @return True on a match, False on EOF or a non-match
 */
static inline int tokenizer_not_a(tokenizer * tk, const char * match) {
  const char * tok = tokenizer_token(tk,0);
  if (!strcmp(tok, END_OF_FILE)) return 0; /* EOF always quits */
  return strcmp(tok, match);
}


/************ 
 * parser 
 ************
 */

/* forward declarations */
static void predict_blockbody( Json& js, tokenizer * tk );
static int predict_schema( Json& js, tokenizer * tk );

static void check_name( const tokenizer * tk, const char * name ) {
    if (strlen(name) && !(isalpha(*name) || *name == '_')) {
        MSYS_FAIL("Line " << tokenizer_line(tk) << " predicted a block name, have '" << name << "'");
    }
}

static inline void parse_value( Json& js, Json::kind_t kind, const char * val ) {
    char *s;
    unsigned len;

    if (!strcmp(val, "<>")) {
        js.to_null();
        return;
    }

    /* strip quotes if present */
    len = strlen(val);
    if (len>1 && val[0] == '"' && val[len-1]=='"') {
        s=strdup(val+1);
        s[len-2]='\0';
    } else {
        s=strdup(val);
    }
    switch (kind) {
        case Json::Bool:
            js.to_bool(atoi(s));
            free(s);
            break;
        case Json::Int:
            js.to_int(atoi(s));
            free(s);
            break;
        case Json::Float:
            js.to_float(atof(s));
            free(s);
            break;
        case Json::String:
            js.to_string(s);
            free(s);
            break;
        default:
            MAE_ERROR1("Unexpected kind %d", kind);
    }
}


static void predict_arraybody( Json& js, tokenizer * tk ) {
    /* read header */
    tokenizer_predict(tk, "[");
    tokenizer_predict(tk, END_OF_FILE); /* FIXME: check item count */
    tokenizer_predict(tk, "]");
    tokenizer_predict(tk, "{");

    /* read schema */
    int ncols = predict_schema( js, tk );
    tokenizer_predict(tk, ":::");

    /* at this point the values in the js are atomic types, rather than
     * lists.  Convert them, saving the type as we go. */
    std::vector<Json::kind_t> kinds(ncols);
    for (int i=0; i<ncols; i++) {
        Json& val = js.elem(i);
        kinds[i] = val.kind();
        val.to_array();
    }

    /* read rows */
    int nrows=0;
    while (tokenizer_not_a(tk, ":::")) {
        tokenizer_predict(tk, END_OF_FILE); /* throw away row index */
        for (int i=0; i<ncols; i++) {
            Json& col = js.elem(i);
            const char * tok = tokenizer_predict_value(tk);
            Json val;
            parse_value( val, kinds[i], tok );
            col.append(val);
        }
        ++nrows;
    }
    tokenizer_predict(tk, ":::");
    // DESRESCode#4634 - allow one level of nested array
    const char* tok = tokenizer_token(tk, false);
    tokenizer_next(tk);
    if (*tok != '}') {
        // got start of a new array body
        MSYS_WARN("WARNING: Skipping nested array '" << tok << "'");
        Json subblock;
        subblock.to_object();
        predict_arraybody( subblock, tk );
        tokenizer_predict(tk, "}");
    }
    /* add row count as __size__ member of js */
    Json size;
    size.to_int(nrows);
    js.append("__size__", size);
    // does this work?  js.append("__size__", Json().to_int(nrows))
}

static void predict_nameless_block( Json& js, const char * name,
                                    tokenizer * tk ) {
    char * key = strdup(name);
    const char * tok;
    Json subblock;
    subblock.to_object();
    /* may be an array */
    tok = tokenizer_token(tk,0);
    if (!strcmp(tok, "[")) {
        predict_arraybody( subblock, tk );
        Json jname;
        jname.to_string(key);
        subblock.append( "__name__", jname );
    }
    /* otherwise just a block */
    else {
        predict_blockbody( subblock, tk );
    }
    js.append(key, subblock);
    free(key);
}

static void predict_block(Json& js, tokenizer * tk ) {
    const char * name = tokenizer_predict(tk, END_OF_FILE);
    check_name(tk,name);
    predict_nameless_block(js, name, tk);
}

static void strip_comments(char* buf, size_t len) {
    for (size_t i=0; i<len; i++) {
        if (buf[i]==' ' && buf[i+1]=='#') {
            buf[i] = '\0';
            break;
        }
    }
}

/* append keyvals to the object.  Return how many were added. */
static int predict_schema( Json& js, tokenizer * tk ) {
    int n = js.size();
    char buf[256];
    while (tokenizer_not_a(tk, ":::")) {
        Json attr;
        const char * token = tokenizer_token(tk,1);
        size_t len = strlen(token);
        if (len+1 >= sizeof(buf)) {
            MAE_ERROR1("schema token '%s' is too long", token);
        }
        memcpy(buf, token, len);
        while ((buf[len++] = tokenizer_peek(tk)) != '\n') {
            if (len == sizeof(buf)) {
                MAE_ERROR1("schema too long at line '%d'", tokenizer_line(tk));
            }
            tokenizer_read(tk);
        }
        buf[len-1] = '\0';
        strip_comments(buf, len-1);
        switch (*buf) {
            case 'b':
                attr.to_bool(false); break;
            case 'i':
                attr.to_int(0); break;
            case 'r':
                attr.to_float(0); break;
            case 's':
                attr.to_string(""); break;
            default:
                MAE_ERROR2("Line %d predicted schema, but '%s' is invalid",
                tokenizer_line(tk), buf);
        }
        js.append(buf+2, attr);
        tokenizer_next(tk);
    }
    return js.size()-n;
}

static void predict_schema_and_values( Json& js, tokenizer * tk ) {
    int istart = js.size();
    int nvalues = predict_schema( js, tk );
    tokenizer_predict(tk, ":::");
    for (int i=0; i<nvalues; i++) {
        Json& elem = js.elem(istart+i);
        parse_value( elem, elem.kind(), tokenizer_predict_value(tk) );
    }
}

static void predict_blockbody( Json& js, tokenizer * tk ) {
    tokenizer_predict(tk, "{");
    predict_schema_and_values(js, tk);
    while (tokenizer_not_a(tk, "}")) {
        predict_block(js, tk);
    }
    tokenizer_predict(tk, "}");
}

static void fill_nameless( Json& js, const char * name, tokenizer * tk ) {
    Json blockname;
    blockname.to_string(name);
    js.append("__name__", blockname);
    predict_blockbody(js,tk);
}

#define INDENT do { int j; for (j=0; j<depth; j++) putc(' ', fd); } while (0)
#define START_BLOCK do { fprintf(fd, "{\n"); depth += 2; } while (0)
#define END_BLOCK   do { depth -=2; INDENT; fprintf(fd, "}\n"); } while (0)

/* scan the key names in the given object.  Infer the type and thus the
 * prefix for the key by examining the type of value, or the type of the
 * first element if the value is a list.  If the value is an object,
 * we ignore it.
 */
static void write_meta( const Json& ct, int depth, FILE * fd ) {
    for (int i=0; i<ct.size(); i++) {
        const char * key = ct.key(i);
        const Json& val = ct.elem(i);
        const char * prefix = NULL;
        if (!strcmp(key, "__name__")) continue;
        if (!strcmp(key, "__size__")) continue;
        Json::kind_t kind = val.kind();
        if (kind == Json::Array) {
            if (!val.size()) continue; /* skipping empty array blocks */
            kind = val.elem(0).kind();
        }
        switch (kind) {
            case Json::Int:    prefix="i_"; break;
            case Json::Bool:   prefix="i_"; break;
            case Json::Float:  prefix="r_"; break;
            case Json::String: prefix="s_"; break;
            default: ;
        }
        if (prefix) {
            INDENT; fprintf(fd, "%s%s\n", prefix, key);
        }
    }
    /* close the labels section */
    INDENT; fprintf(fd, ":::\n");
}

/* Return NULL if the string can be printed as-is; otherwise return malloc'ed
 * version. */
static char * quotify( const char * s ) {
    const char * p;
    std::string str;
    int needs_escape = 0;
    if (!s) return strdup("<>");
    if (!strlen(s)) return strdup("\"\"");
    for (p=s; *p; ++p) {
        if (isspace(*p) || !isprint(*p) || *p=='"' || *p=='<' || *p=='\\') {
            needs_escape = 1;
            //if (isspace(*p) && !(*p==' ' || *p=='\t'))
                //MAE_ERROR1("unprintable whitespace in <%s>", s);
        }
    }
    if (!needs_escape) return NULL;
    str += '"';
    for (p=s; *p; ++p) {
        if (*p=='"') {
            str += "\\\"";
        } else if (*p=='\\') {
            str += "\\\\";
        } else {
            str += *p;
        }
    }
    str += '"';
    return strdup(str.c_str());
}

static void write_array( const Json& arr, 
                         int size, 
                         int depth, 
                         FILE * fd ) {

    const char * s;
    char * q;
    for (int i=0; i<size; i++) {
        INDENT;
        /* first column is row number */
        fprintf(fd, "%d ", i+1);
        for (int j=0; j<arr.size(); j++) {
            const Json& col = arr.elem(j);
            if (col.kind()!=Json::Array) continue;
            const Json& elem = col.elem(i);
            switch(elem.kind()) {
                case Json::Int:
                    fprintf(fd, "%d ", elem.as_int()); break;
                case Json::Bool:
                    fprintf(fd, "%d ", elem.as_bool()); break;
                case Json::Float:
                    fprintf(fd, "%g ", elem.as_float()); break;
                case Json::String:
                    s = elem.as_string();
                    q = quotify(s);
                    fprintf(fd, "%s ", q ? q : s);
                    if (q) free(q);
                    break;
                default:
                    ;
            }
        }
        fprintf(fd, "\n");
    }
    INDENT; fprintf(fd, ":::\n");
}

static void write_values( const Json& ct, int depth, FILE * fd ) {
    const char * s;
    char * q;
    for (int i=0; i<ct.size(); i++) {
        const char * key = ct.key(i);
        const Json& val = ct.elem(i);
        if (!strcmp(key, "__name__")) continue;
        INDENT;
        switch (val.kind()) {
            case Json::Int:
                fprintf(fd, "%d\n", val.as_int()); break;
            case Json::Bool:
                fprintf(fd, "%d\n", val.as_bool()); break;
            case Json::Float:
                fprintf(fd, "%g\n", val.as_float()); break;
            case Json::String:
                s = val.as_string();
                q = quotify(s);
                fprintf(fd, "%s\n", q ? q : s);
                if (q) free(q);
                break;
            case Json::Array:
                throw std::runtime_error("Unexected array type in json");
                break;
            case Json::Object:
                {
                fprintf(fd, "%s", key);
                const Json& size = val.get("__size__");
                if (size.valid()) {
                    fprintf(fd, "[%d]", size.as_int());
                }
                fprintf(fd, " ");
                START_BLOCK;
                write_meta( val, depth, fd );
                if (size.kind()) write_array( val, size.as_int(), depth, fd );
                else             write_values( val, depth, fd );
                END_BLOCK;
                }
                break;
            default: ;
        }
    }
}

namespace desres { namespace msys { namespace mae {

    import_iterator::import_iterator(std::istream& file) 
    : in(maybe_compressed_istream(file)), tk(), _offset() {
        tk = new tokenizer;
        tokenizer_init(tk, *in);
    }

    import_iterator::~import_iterator() {
        tokenizer_release(tk);
        delete tk;
    }

    bool import_iterator::next(Json& block) {
        _offset = tk->m_offset + tk->bufpos;
        /* eat the meta block, if any */
        while (!strcmp("{", tokenizer_token(tk,0))) {
            tokenizer_predict(tk, "{");
            block.to_object();
            predict_schema_and_values(block, tk);
            tokenizer_predict(tk, "}");
            _offset = tk->m_offset + tk->bufpos;
        }
        if (tokenizer_not_a(tk, END_OF_FILE)) {
            block.to_object();
            const char * name = tokenizer_predict(tk, END_OF_FILE);
            fill_nameless( block, name, tk );
            return true;
        }
        return false;
    }

    void import_mae( std::istream& input, Json& js ) {
        Json block;
        js.to_array();
        import_iterator it(input);
        while (it.next(block)) js.append(block);
    }

}}}
