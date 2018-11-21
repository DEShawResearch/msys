/* @COPYRIGHT@ */

#include "parse.hxx"
#include "JSON_parser.h"

#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <fstream>

using desres::msys::fastjson::Json;

/*

   Start with a stack of size 1.

   When an array starts, stuff that data in, then push the stack.
   When an object starts, start that data in, then push the stack.

   When any value comes in, either atom or composite, stuff that
   data into the top of the stack.  Check for previous stack value.
   If it's non-NULL, then it must be composite, and we need to
   append top value in top onto that composite.

   When a key comes in, copy the context key.  Appending to an
   object steals the key and sets it back to NULL.

   When an array or object ends, just pop.

   When we're all done, we should have a stack of size 1 again, and
   that's our return value.

   Oops!  If it's an array, we can't append it to prev until we're
   all completed because the elems pointer might move around.  So,
   we don't have a newval for the array until we end it.

*/

#define PUSH() do { \
    list * node = new list; \
    node->next = head[0]; \
    head[0] = node; \
} while (0)

#define POP() do { \
    list * node = head[0]->next; \
    delete head[0]; \
    head[0] = node; \
} while (0)

namespace {

    struct list {
        list *  next;
        char *  key;
        Json    json;

        list() : next(NULL), key(NULL) {}
        ~list() { if (key) free(key); }
    };

    int build(void *ctx_, int type, const JSON_value *value) {
        list ** head = (list **)ctx_;
        //json_parser_context * ctx = (json_parser_context *)ctx_;

        Json * newval = NULL;
        Json & json = head[0]->json;
        list * prev = head[0]->next;

        switch (type) {
            case JSON_T_ARRAY_BEGIN:
                json.to_array();
                PUSH();
                break;
            case JSON_T_ARRAY_END:
                POP();
                newval = &head[0]->json;
                prev = head[0]->next;
                break;
            case JSON_T_OBJECT_BEGIN:
                json.to_object();
                PUSH();
                break;
            case JSON_T_OBJECT_END:
                POP();
                newval = &head[0]->json;
                prev = head[0]->next;
                break;
            case JSON_T_KEY:
                head[0]->key = strdup(value->vu.str.value);
                break;
            case JSON_T_INTEGER:
                newval = &json.to_int(value->vu.integer_value);
                break;
            case JSON_T_STRING:
                newval = &json.to_string(value->vu.str.value);
                break;
            case JSON_T_FLOAT:
                newval = &json.to_float(value->vu.float_value);
                break;
            case JSON_T_NULL:
                newval = &json.to_null();
                break;
            case JSON_T_TRUE:
                newval = &json.to_bool(true);
                break;
            case JSON_T_FALSE:
                newval = &json.to_bool(false);
                break;
            default:
                return 0;
        }

        if (newval && prev) {

            if (prev->json.kind() == Json::Array) {
                prev->json.append(*newval);

            } else if (prev->json.kind() == Json::Object) {
                prev->json.append(head[0]->key, *newval);
                free(head[0]->key);
                head[0]->key=NULL;

            } else {
                return 0;
            }
        }

        return 1;
    }

    struct vjc {
        list head;
        list * ctx;
        JSON_config config;
        JSON_parser_struct * jc;

        vjc() {
            ctx = &head;
            init_JSON_config(&config);
            config.depth                  = 20;
            config.callback               = &build;
            config.callback_ctx           = &ctx;
            config.allow_comments         = 1;
            config.handle_floats_manually = 0;

            jc = new_JSON_parser(&config);
        }
        ~vjc() {
            if (jc) delete_JSON_parser(jc);
            list* node = ctx;
            while (node != &head) {
                list* tmp = node->next;
                delete node;
                node = tmp;
            }
        }
    };

}

namespace desres { namespace msys { namespace fastjson {

    void parse_json( std::istream& in, Json& js ) {
        int col=1;
        int line=1;
        char buf[4096];
        vjc v[1];

        if (in.peek()==EOF || !in.good()) {
            throw std::runtime_error("bad input stream");
        }

        while (!in.eof()) {
            in.read(buf,sizeof(buf));
            int len=in.gcount();
            for (int i=0; i<len; i++) {
                if (!JSON_parser_char(v->jc, buf[i])) {
                    std::stringstream ss;
                    ss << "Parser failed near line " 
                        << line << " column " << col;
                    throw std::runtime_error(ss.str());
                }
                ++col;
                if (buf[i]=='\n') {
                    ++line;
                    col=1;
                }
            }
        }
        if (!JSON_parser_done(v->jc)) {
            throw std::runtime_error("Parser failed at end of input");
        }
        if (v->head.next) {
            throw std::runtime_error("Unexpected stack size > 1");
        }
        js.swap( v->head.json );
    }

    void parse_json( const char * path, Json& js ) {
       std::ifstream in(path);
       parse_json(in,js);
       in.close();
    }
}}}

std::istream& operator>>(std::istream& in, desres::msys::fastjson::Json& js) {
    parse_json( in, js );
    return in;
}
