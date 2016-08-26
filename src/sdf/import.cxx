#include <cmath>

#include "../sdf.hxx"
#include "elements.hxx"
#include "../append.hxx"

#include <stdio.h>
#include <errno.h>
#include <string>

#include <fstream>
#include <sstream>

using namespace desres::msys;

namespace {
    void add_typed_keyval(String const& key, String const& val, 
                          component_t& ct) {
        try {
            int v = stringToInt(val);
            ct.add(key,IntType);
            ct.value(key)=v;
            return;
        } catch (Failure& e) {
        }
        try {
            double v = stringToDouble(val);
            if (std::isfinite(v)) {
                ct.add(key,FloatType);
                ct.value(key)=v;
                return;
            }
        } catch (Failure& e) {
        }
        ct.add(key,StringType);
        ct.value(key)=val;
    }

    class iterator : public LoadIterator {
    protected:
        char buf[1024];
        int line = 0;
        virtual bool getline() = 0;
        virtual bool eof() = 0;

    private:
        std::string skip_to_end() {
            std::string current = std::to_string(line)+": " + buf;
            while (getline()) {
                if (!strncmp(buf, "$$$$", 4)) break;
            }
            return current;
        }

        static const short BAD_COUNT = 9999;
        // parse text of the form __z, _yz, or xyz, where underscore
        // denotes leading space and x,y,z are digits.
        static short parse_count(const char* s) {
            int n=0;
            unsigned short result=0;
            unsigned short sign = 1;
            if (isspace(*s)) { ++n; ++s; }
            if (isspace(*s)) { ++n; ++s; }
            if (*s == '-') {
                sign = -sign;
                ++s;
                ++n;
            }
            for (; n<3; n++, s++) {
                if (!isdigit(*s)) return BAD_COUNT;
                unsigned short d = *s - '0';
                switch (n) {
                    case 0: result += 100*d; break;
                    case 1: result +=  10*d; break;
                    case 2: result +=     d; break;
                };
            }
            return sign*result;
        }

        // parse xxxxx.yyyy as float
        static float parse_coord(const char* s) {
            float scale = 10000.0f;
            unsigned result = 0;
            static const float fail = std::numeric_limits<float>::max();
            int n=0;
            for (; n<4; n++, s++) {
                if (!isspace(*s)) break;
            }
            if (*s == '-') {
                ++s;
                ++n;
                scale = -scale;
            }
            for (; n<5; n++, s++) {
                if (!isdigit(*s)) return fail;
                unsigned d = *s - '0';
                switch (n) {
                    case 0: result += 100000000*d; break;
                    case 1: result +=  10000000*d; break;
                    case 2: result +=   1000000*d; break;
                    case 3: result +=    100000*d; break;
                    case 4: result +=     10000*d; break;
                }
            }
            if (*s != '.') return fail;
            result += 1000*(s[1]-'0');
            result +=  100*(s[2]-'0');
            result +=   10*(s[3]-'0');
            result +=      (s[4]-'0');
            return float(result)/scale;
        }

        static char parse_element(const char* s) {
            if (isspace(*s)) ++s;
            if (isspace(*s)) {
                // one-character element
                ++s;
                switch (*s) {
                    case 'C': return 6;
                    case 'H': return 1;
                    case 'N': return 7;
                    case 'O': return 8;
                    case 'F': return 9;
                    case 'P': return 15;
                    case 'S': return 16;
                    case 'K': return 19;
                    case 'B': return 5;
                    case 'V': return 23;
                    case 'Y': return 39;
                    case 'I': return 53;
                    case 'W': return 74;
                    case 'U': return 92;
                    default:  return 0;
                };
            }
            char buf[3] = {s[0],s[1],'\0'};
            return ElementForAbbreviationSlow(buf);
        }

    public:
        SystemPtr next() {
            SystemPtr ptr;

            // three lines for header block, then one for counts
            if (!getline()) return ptr;
            std::string name = buf;
            name.pop_back();    // remove trailing newline
            getline();
            getline();
            getline();
            if (eof()) return ptr;
            auto natoms = parse_count(buf);
            auto nbonds = parse_count(buf+3);
            if (natoms==BAD_COUNT || nbonds==BAD_COUNT) {
                MSYS_FAIL("Bad counts line: " << skip_to_end());
            }

            ptr = System::create();
            System& mol = *ptr;
            mol.addChain();
            mol.addResidue(0);
            mol.name = name;
            mol.ct(0).setName(name);

            // atoms
            for (unsigned short i=0; i<natoms; i++) {
                if (!getline()) {
                    skip_to_end();
                    MSYS_FAIL("Missing expected atom record");
                }
                auto sz = strlen(buf);
                if (sz<34) {
                    MSYS_FAIL("Malformed atom line: " << skip_to_end());
                }
                auto& atm = mol.atomFAST(mol.addAtom(0));
                atm.x = parse_coord(buf   );
                atm.y = parse_coord(buf+10);
                atm.z = parse_coord(buf+20);
                atm.atomic_number = parse_element(buf+31);
                atm.name = AbbreviationForElement(atm.atomic_number);

                if (sz>=39) {
                    auto q = parse_count(buf+36);
                    if (q==BAD_COUNT) {
                        MSYS_FAIL("Bad charge: " << skip_to_end());
                    }
                    atm.formal_charge = q==0 ? 0 : 4-q;
                }
                if (sz>=42) {
                    auto s = parse_count(buf+39);
                    if (s==BAD_COUNT) {
                        MSYS_FAIL("Bad stereo: " << skip_to_end());
                    }
                    atm.stereo_parity = s;
                }
            }

            // bonds
            for (unsigned short i=0; i<nbonds; i++) {
                if (!getline()) {
                    skip_to_end();
                    MSYS_FAIL("Missing expected bond record");
                }
                auto sz = strlen(buf);
                if (sz<9) {
                    MSYS_FAIL("Malformed bond line: " << skip_to_end());
                }
                auto ai = parse_count(buf  )-1;
                auto aj = parse_count(buf+3)-1;
                auto& bnd = mol.bondFAST(mol.addBond(ai,aj));
                bnd.order = parse_count(buf+6);
                bnd.stereo = parse_count(buf+9);
                if ((bnd.stereo==1 || bnd.stereo==6) && ai>aj) {
                    // msys stores the bonds with ai<aj, but flipping the
                    // order of the atoms necessitates flipping the
                    // absolute stereo specification.
                    bnd.stereo = -bnd.stereo;
                }
                if (bnd.order == 4) {
                    bnd.order = 1;
                    bnd.resonant_order = 1.5;       
                    bnd.aromatic = 1;
                } else if (bnd.order<1 || bnd.order>4) {
                    MSYS_FAIL("Unsupported bond type in bond record: " << skip_to_end());
                }
            }

            // M entries
            bool cleared_charges = false;
            while (getline()) {
                if (!strncmp(buf, "M  ", 3)) {
                    if (!strncmp(buf+3, "END", 3)) {
                        break;
                    } else if (!strncmp(buf+3, "CHG", 3)) {
                        if (!cleared_charges) {
                            cleared_charges = true;
                            for (int i=0; i<natoms; i++) {
                                mol.atomFAST(i).formal_charge = 0;
                            }
                        }
                        int n = parse_count(buf+6);
                        if (n==BAD_COUNT) {
                            MSYS_FAIL("Malformed CHG line: " << skip_to_end());
                        }
                        for (int i=0; i<n; i++) {
                            short aid = parse_count(buf+10+8*i);
                            short chg = parse_count(buf+14+8*i);
                            if (aid==BAD_COUNT || chg==BAD_COUNT) {
                                MSYS_FAIL("Malformed CHG line: " << skip_to_end());
                            }
                            mol.atom(aid-1).formal_charge = chg;
                        }
                    }
                } else if (!strncmp(buf, "A  ", 3)) {
                    getline();
                } else if (!strncmp(buf, "G  ", 3)) {
                    getline();
                } else if (!strncmp(buf, "V  ", 3)) {
                    // ignore
                } else {
                    MSYS_FAIL("Malformed properties line: " << skip_to_end());
                }
            }

            // data fields.  
            bool needline = true;
            for (;;) {
                if (needline && !getline()) {
                    MSYS_FAIL("Unexpected end of file");
                }
                needline = true;
                if (!strncmp(buf, "$$$$", 4)) {
                    break;
                } else if (!strncmp(buf, "> ", 2)) {
                    const char* langle = strchr(buf+2,'<');
                    const char* rangle = strchr(buf+3,'>');
                    std::string key;
                    if (langle && rangle) {
                        key = std::string(langle+1,rangle);
                    }
                    std::string val;
                    // support multiline values containing newlines
                    for (;;) {
                        if (!getline()) {
                            MSYS_FAIL("Unexpected end of file");
                        }
                        if (!strncmp(buf, "$$$$", 4) ||
                            !strncmp(buf, "> ", 2)) {
                            if (*key.c_str()) {
                                auto sz = val.size();
                                if (sz==1) val.clear(); /* just a newline */
                                else if (sz>1) val.resize(sz-2);
                                add_typed_keyval(key,val,mol.ct(0));
                            }
                            needline = false;
                            break;
                        }
                        val += buf;
                    }
                } else if (buf[0]=='\n') {
                    // allow blank line
                } else {
                    MSYS_FAIL("Malformed data field line: " << skip_to_end());
                }
            }
            return ptr;
        }
    };

    class text_iterator : public iterator {
        std::string data;
        size_t pos = 0;

    protected:
        bool getline() {
            if (pos>=data.size()) return false;
            ++line;
            size_t len = strcspn(&data[pos], "\n");
            if (len>sizeof(buf)) MSYS_FAIL("Line " << line << " too long");
            memcpy(buf, &data[pos], len+1);
            buf[len+1]='\0';
            pos += len+1;
            return true;
        }
        bool eof() {
            return pos>=data.size();
        }

    public:
        text_iterator(std::string const& d) : data(d) {}
    };

    class file_iterator : public iterator {
        FILE* fp = nullptr;
        int (*closer)(FILE *);

    protected:
        bool getline() {
            ++line;
            return fgets(buf, sizeof(buf), fp)!=NULL;
        }
        bool eof() {
            return feof(fp);
        }

    public:
        file_iterator(std::string const& path) {
            fp = fopen(path.data(), "rb");
            if (!fp) MSYS_FAIL(strerror(errno));
            closer = fclose;

            // handle gzipped files
            unsigned char buf[2];
            size_t rc = fread(buf, 1, 2, fp);
            if (rc<2) {
                if (feof(fp)) return;
                MSYS_FAIL(strerror(errno));
            }
            fseek(fp, 0, SEEK_SET);
            if (buf[0]==0x1f && buf[1]==0x8b) {
                fclose(fp);
                fp = NULL;
                std::string cmd("gzip -dc \"");
                cmd += path;
                cmd += "\"";
                fp = popen(cmd.data(), "r");
                if (!fp) {
                    MSYS_FAIL("popen failed: " << strerror(errno));
                }
                closer = pclose;
            }
        }
        ~file_iterator() {
            if (fp) closer(fp);
        }
    };
}

SystemPtr desres::msys::ImportSdf(std::string const& path) {
    file_iterator iter(path);
    SystemPtr ct, mol = System::create();
    mol->name = path;
    while ((ct = iter.next())) {
        AppendSystem(mol, ct);
    }
    return mol;
}

LoadIteratorPtr desres::msys::SdfFileIterator(std::string const& path) {
    return LoadIteratorPtr(new file_iterator(path));
}
LoadIteratorPtr desres::msys::SdfTextIterator(std::string const& data) {
    return LoadIteratorPtr(new text_iterator(data));
}

