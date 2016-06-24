#include "token.hxx"

using namespace desres::msys::atomsel;

extern void* atomselParseAlloc(void* (mallocProc)(size_t));
extern void atomselParseFree(void*, void (freeProc)(void*));
extern void atomselParse(void*, int, desres::msys::atomsel::Token, desres::msys::atomsel::Query*);
extern void atomselParseTrace(FILE*, char*);

static void fail_with_error(std::string const& sel, int start, int stop) {
    std::stringstream ss;
    ss << "Parse failed:\n" << sel << "\n";
    for (int i=0; i<start; i++) ss << " ";
    ss << "^";
    for (int i=start; i<stop-2; i++) ss << "-";
    ss << "^\n";
    MSYS_FAIL(ss.str());
}

void Query::parse(std::string const& sel) {
    if (sel.empty()) MSYS_FAIL("empty selection");
    void* p = atomselParseAlloc(malloc);
    std::shared_ptr<void> dtpr(p, [](void *v) { atomselParseFree(v,free); });
    atomsel::Tokenizer tk(sel.data());
    if (getenv("ATOMSEL_DEBUG")) {
        atomselParseTrace(stderr, (char *)"atomsel: ");
    }

    int tokenId;
    do {
        atomsel::Token token;
        int start = tk.location();
        tokenId = tk.next(&token, mol);
        int stop = tk.location();
        if (tokenId<0) {
            fail_with_error(sel,start,stop);
        }
        if (tokenId==MACRO) {
            parse(token.str());
        }
        atomselParse(p, tokenId, token, this);
        if (!mol) {
            fail_with_error(sel,start,stop);
        }
    } while (tokenId>0);
}

