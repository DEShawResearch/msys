#include <cmath>

#include "../sdf.hxx"
#include "../elements.hxx"
#include "../append.hxx"

#include <stdio.h>
#include <errno.h>
#include <string>

#include <fstream>
#include <sstream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

using namespace desres::msys;

namespace {
    void add_typed_keyval(String const& key, String& val, 
                          component_t& ct) {
        if (key.empty()) return;
        auto sz = val.size();
        if (sz==1) val.clear(); /* just a newline */
        else if (sz>1) val.resize(sz-2);
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
        virtual size_t offset() const = 0;
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
                    bnd.aromatic = 1;
                } else if (bnd.order<0 || bnd.order>4) {
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
                    } else if (!strncmp(buf+3, "ISO", 3)) {
                        int n = parse_count(buf+6);
                        if (n==BAD_COUNT) {
                            MSYS_FAIL("Malformed ISO line: " << skip_to_end());
                        }
                        auto id = mol.addAtomProp("isotope", IntType);
                        for (int i=0; i<n; i++) {
                            short aid = parse_count(buf+10+8*i);
                            short iso = parse_count(buf+14+8*i);
                            if (aid==BAD_COUNT || iso==BAD_COUNT) {
                                MSYS_FAIL("Malformed ISO line: " << skip_to_end());
                            }
                            mol.atomPropValue(aid-1, id) = iso;
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
            if (!getline()) return ptr;
            bool needline = false;
            for (;;) {
                if (needline && !getline()) {
                    MSYS_FAIL("Unexpected end of file");
                }
                if (eof()) break;
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
                            add_typed_keyval(key,val,mol.ct(0));
                            needline = false;
                            break;
                        }
                        if (!strncmp(buf, "$$$$", 4) ||
                            !strncmp(buf, "> ", 2)) {
                            add_typed_keyval(key,val,mol.ct(0));
                            needline = false;
                            break;
                        }
                        val += buf;
                    }
                } else if (buf[0]=='\n' || (buf[0]=='\r' && buf[1]=='\n')) {
                    // allow blank line
                } else {
                    MSYS_FAIL("Malformed data field line: " << skip_to_end());
                }
            }
            return ptr;
        }
    };

    class buffer_iterator : public iterator {
        const char* data = nullptr;
        size_t size = 0;
        size_t pos = 0;

    protected:
        bool getline() {
            if (pos>=size) return false;
            ++line;
            size_t len = strcspn(&data[pos], "\n");
            if (len>sizeof(buf)) MSYS_FAIL("Line " << line << " too long");
            memcpy(buf, &data[pos], len+1);
            buf[len+1]='\0';
            pos += len+1;
            return true;
        }
        bool eof() {
            return pos>=size;
        }

    public:
        size_t offset() const {
            return pos;
        }
        buffer_iterator(const char* d, size_t sz) : data(d), size(sz) {}
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
        size_t offset() const {
            return pos;
        }
        text_iterator(std::string const& d) : data(d) {}
    };

    class file_iterator : public iterator {
        FILE* fp = nullptr;
        int (*closer)(FILE *);

    protected:
        bool getline() {
            ++line;
            bool rc = fgets(buf, sizeof(buf), fp)!=NULL;
            return rc;
        }
        bool eof() {
            return feof(fp);
        }
    public:
        size_t offset() const {
            return ftello(fp);
        }

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
    mol->updateFragids();
    return mol;
}

// return offset in bytes of each sdf entry
static std::vector<size_t> sdf_offsets(std::string const& path) {
    file_iterator iter(path);
    SystemPtr ct;
    std::vector<size_t> offsets;
    while ((ct = iter.next())) {
        offsets.push_back(iter.offset());
        //if (!(offsets.size() % 100000)) printf("%lu\n", offsets.size());
    }
    return offsets;
}

namespace {
    // index file format: 
    // byte 0: version = 0x01
    // byte 1-7: unused
    // byte 8-15: num_entries
    // remainder: 8-byte offsets

    class IndexedSdfLoader : public IndexedFileLoader {
        int sdf_fd = 0;
        int idx_fd = 0;
        size_t _size = 0;
        std::string _path;

        void parse_header()
        {
            unsigned char buf[16];
            if (::read(idx_fd, buf, sizeof(buf)) != sizeof(buf)) {
                MSYS_FAIL("Parsing idx file header: " << strerror(errno));
            }
            if (buf[0] != 0x01) {
                MSYS_FAIL("Bad version in header: got " << buf[0] << " want " << 0x01);
            }
            _size = ((size_t *)buf)[1];
        }

    public:
        IndexedSdfLoader(std::string const& sdf_path, 
                         std::string const& idx_path)
        : _path(sdf_path)
        {
            idx_fd = ::open( idx_path.data(), O_RDONLY);
            if (idx_fd<=0) {
                MSYS_FAIL(idx_path << ": " << strerror(errno));
            }
            sdf_fd = ::open(sdf_path.data(), O_RDONLY);
            if (sdf_fd <= 0) {
                MSYS_FAIL("Opening sdf file: " << strerror(errno));
            }
            parse_header();
        }

        ~IndexedSdfLoader() {
            if (sdf_fd > 0) close(sdf_fd);
            if (idx_fd > 0) close(idx_fd);
        }


        std::string const& path() const { return _path; }
        size_t size() const { return _size; }
        SystemPtr at(size_t i) const {
            if (i>=_size) MSYS_FAIL("Invalid index " << i << " >= " << _size);

            // get range of needed data
            size_t range[2];
            if (i==0) {
                range[0] = 0;
                if (pread(idx_fd, &range[1], 8, 16)!=8) {
                    MSYS_FAIL("Reading index entry 0: " << strerror(errno));
                }
            } else {
                if (pread(idx_fd, &range[0], 16, 8+8*i)!=16) {
                    MSYS_FAIL("Reading index entry " << i << ": " << strerror(errno));
                }
            }

            size_t size = range[1] - range[0];
            std::vector<char> buf(size);
            ssize_t rc = pread(sdf_fd, buf.data(), size, range[0]);
            if (rc!=(ssize_t)(size)) {
                MSYS_FAIL("Reading entry " << i << " from sdf " << _path << ": " << strerror(errno));
            }
            return buffer_iterator(buf.data(), buf.size()).next();
        }
    };
}

void desres::msys::CreateIndexedSdf(std::string const& sdf_path, 
                                    std::string const& idx_path) {
    auto offsets = sdf_offsets(sdf_path);
    FILE* fp = fopen(idx_path.data(), "wb");
    if (!fp) MSYS_FAIL(idx_path << ": " << strerror(errno));
    if (offsets.empty()) {
        fclose(fp);
        return;
    }
    auto n = offsets.size();
    char header[16];
    header[0] = 0x01;
    ((size_t*)header)[1] = n;
    fwrite(header, 1, sizeof(header), fp);
    if (fwrite(&offsets[0], sizeof(offsets[0]), n, fp) != n) {
        std::string err = strerror(errno);
        fclose(fp);
        unlink(idx_path.data());
        MSYS_FAIL(idx_path << ": failed to write complete index: " << err);
    }
    fclose(fp);
}

std::shared_ptr<IndexedFileLoader>
desres::msys::OpenIndexedSdf(std::string const& path, std::string const& idx) {
    return std::make_shared<IndexedSdfLoader>(path, idx);
}

LoadIteratorPtr desres::msys::SdfFileIterator(std::string const& path) {
    return LoadIteratorPtr(new file_iterator(path));
}
LoadIteratorPtr desres::msys::SdfTextIterator(std::string const& data) {
    return LoadIteratorPtr(new text_iterator(data));
}

