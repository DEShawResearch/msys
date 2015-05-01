#ifndef desres_msys_inchi_hxx
#define desres_msys_inchi_hxx

#include "system.hxx"

namespace desres { namespace msys {


    /* data from an inchi invocation */
    class InChI {
        int _rc;
        String _string;
        String _auxinfo;
        String _message;

        InChI(int rc, const char* s, const char* a, const char* m)
        : _rc(rc),
          _string (s ? s : ""),
          _auxinfo(a ? a : ""),
          _message(m ? m : "")
          {}

    public:
        enum Flags {
            /* structure perception (compatible with standard inchi) */
            DoNotAddH = 1 << 0,
            SNon      = 1 << 1,

            /* InChI creation options (generates non-standard InChI) */
            FixedH    = 1 << 2 
        };

        /* create inchi */
        static InChI create(SystemPtr mol, unsigned options=0);

        String const& string() const { return _string; }
        String const& auxinfo() const { return _auxinfo; }
        String const& message() const { return _message; }
        bool ok() const { return _rc==0; }

        /* compute inchi key */
        String key() const;
    };

}}

#endif
