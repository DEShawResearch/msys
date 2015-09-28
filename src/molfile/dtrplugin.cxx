//
// Version info for VMD plugin tree:
//   $Id$
//
// Version info for last sync with D. E. Shaw Research:
//  //depot/desrad/main/sw/libs/molfile/plugins/dtrplugin.cxx#30
//

/*
Copyright 2009, D. E. Shaw Research, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research, LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "dtrplugin.hxx"
#include "dtrframe.hxx"
#include "dtrutil.hxx"

using namespace desres::molfile;
using namespace desres::molfile::dtr;

#include <sstream>
#include <ios>
#include <iomanip>
#include <math.h>
#include <errno.h>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>
#include <ios>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "vmddir.h"

#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace bfs = boost::filesystem;
namespace bio = boost::iostreams;

static const char SERIALIZED_VERSION[] = "0008";

const char * desres::molfile::dtr_serialized_version() {
    return SERIALIZED_VERSION;
}

#ifndef DESRES_WIN32
static const char s_sep = '/';

#include <netinet/in.h> /* for htonl */
#if defined(_AIX)
#include <fcntl.h>
#else
#include <sys/fcntl.h>
#endif

#if defined(__sun)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#else
/// windows version

#define M_PI (3.1415926535897932385)
#define M_PI_2 (1.5707963267948966192)

#ifndef S_ISREG
#define S_ISREG(x) (((x) & S_IFMT) == S_IFREG)
#endif

#ifndef S_ISDIR
#define S_ISDIR(x) (((x) & S_IFMT) == S_IFDIR)
#endif

static const char s_sep = '\\';

#endif

static const uint32_t magic_timekey = 0x4445534b;

namespace {

  const double PEAKmassInAmu = 418.4;

  double sfxp_ulp32flt(int32_t x) {
    return ldexp(((double) x),-31);
  }

  typedef struct key_prologue {
    uint32_t magic;           /* Magic number for frames */
    uint32_t frames_per_file; /* Number of frames in each file */
    uint32_t key_record_size; /* The size of each key record */
  } key_prologue_t;


  /*!
   * The byte order associated with this machine.  We use
   * 1234 for little endian, 4321 for big endian, and
   * 3412 for the unlikely PDB endianism.
   */

  // unused(!)
#if 0
  uint32_t machineEndianism() {
#if __BYTE_ORDER == __LITTLE_ENDIAN
    uint32_t byteorder = 1234;
#else
#if __BYTE_ORDER == __BIG_ENDIAN
    uint32_t byteorder = 4321;
#else
#ifdef PDB_ENDIAN
#if __BYTE_ORDER == __PDB_ENDIAN
    uint32_t byteorder = 3412;
#endif
#endif
#endif
#endif
    // If we get a compile error here, then __BYTE_ORDER
    // has an unexpected value.
    return byteorder;
  }
#endif

  bool isfile(const std::string &name) {
    struct stat statbuf;
    return (stat(name.c_str(),&statbuf) == 0 && S_ISREG(statbuf.st_mode));
  }

  bool exists(const std::string& path) {
    struct stat statbuf;
#ifdef DESRES_WIN32
    // Use ::stat instead of ::lstat on windows since there are no symlinks
    return (stat(path.c_str(),&statbuf) == 0);
#else
    return (::lstat(path.c_str(),&statbuf) == 0);
#endif
  }

  /*!
   * Remove a file or directory.  For directories,
   * we recurse through subfiles and remove those
   * before attempting the ::rmdir();
   */
  void recursivelyRemove(std::string path) {
    struct stat statbuf;

    // -----------------------------------------------
    // Only try to unlink if the file exists
    // We recurse through directories and unlink
    // other files.
    // -----------------------------------------------

#ifdef DESRES_WIN32
    // Use ::stat instead of ::lstat on windows since there are no symlinks
    if (stat(path.c_str(),&statbuf) == 0) {
#else
    if (::lstat(path.c_str(),&statbuf) == 0) {
#endif
      if (!S_ISDIR(statbuf.st_mode)) {
        if (::unlink(path.c_str()) != 0) {
            throw std::runtime_error(strerror(errno));
        }
      } else {
        VMDDIR* directory = NULL;
        try {
          directory = vmd_opendir(path.c_str());
          if (directory) {
            // Remove subfiles
            char * entry;
            while( (entry=vmd_readdir(directory)) != NULL ) {
              // Don't unlink . or ..
              if (entry[0] == '.') {
                if (entry[1] == 0) continue;
                if (entry[1] == '.' && entry[2] == 0) continue;
              }
              recursivelyRemove(path + s_sep + entry);
            }
            vmd_closedir(directory);
            directory = NULL;

            // Remove the actual directory
            if (::rmdir(path.c_str()) != 0) {
              throw std::runtime_error(strerror(errno));
            }
          }
        } catch(...) {
          if (directory) vmd_closedir(directory);
          throw;
        }
      }
    }
  }


  ////////
  // CRC
  ////////
  
  typedef uint32_t crc;
  
  #define POLYNOMIAL 0x04C11DB7
  #define WIDTH  (8 * sizeof(crc))
  #define TOPBIT (1 << (WIDTH - 1))
  #define FINAL_XOR_VALUE 0xFFFFFFFF

  crc processByte( crc remainder, char msg ) {
        remainder ^= (msg << (WIDTH - 8));
        for (uint8_t bit = 8; bit > 0; --bit)
        {
            if (remainder & TOPBIT) {
                remainder = (remainder << 1) ^ POLYNOMIAL;
            } else {
                remainder = (remainder << 1);
            }
        }
        return remainder;
  }

  crc processBytes(const char *message, int nBytes) {
    crc  remainder = 0;        
    for (int byte = 0; byte < nBytes; ++byte) {
        remainder = processByte( remainder, message[byte] );
    }
    return remainder;
  } 

  int32_t cksum(const std::string &s) {
    ssize_t len = s.size();
    int32_t result = processBytes( s.c_str(), len );
  
    for ( ; len; len >>= 8) {
      result = processByte( result, len & 0xff );
    }
    return result ^ FINAL_XOR_VALUE;
  }

}

namespace {
    struct defer_close {
        int _fd;
        explicit defer_close(int fd) : _fd(fd) {}
        ~defer_close() { close(_fd); }
    };
}

/* read bytes from the given file at the given offset.  If buffer.size()==0,
 * read all the remaining bytes after the offset and resize buffer.  Otherwise,
 * read only buffer.size() bytes into buffer. */
static void read_file( std::string const& path, 
                       std::vector<char>& buffer,
                       off_t offset = 0) {

    int fd = open( path.c_str(), O_RDONLY|O_BINARY);
    if (fd<=0) {
        DTR_FAILURE("Unable to read file at " << path << "\n  " << strerror(errno));
    }
    defer_close _(fd);
    if (buffer.size()==0) {
        struct stat statbuf;
        if (fstat(fd,&statbuf)!=0) {
            DTR_FAILURE("Could not stat file: " << strerror(errno));
        }
        if (statbuf.st_size < offset) {
            DTR_FAILURE("offset " << offset << " is greater than size " << statbuf.st_size << " of file " << path);
        }
        buffer.resize(statbuf.st_size-offset);
    }

    ssize_t rc = pread(fd, &buffer[0], buffer.size(), offset);
    if (rc<0) {
        DTR_FAILURE("reading " << path << ": " << strerror(errno));
    }
    if (size_t(rc) != buffer.size()) {
        DTR_FAILURE("unexpected short read of file " << path);
    }
}

void Timekeys::init(const std::string& path ) {
    std::string timekeys_path = path;
    timekeys_path += s_sep;
    timekeys_path += "timekeys";
    std::vector<char> buffer;
    try {
        read_file(timekeys_path, buffer);
        initWithBytes(buffer.size(), &buffer[0]);
    } catch (std::exception& e) {
        DTR_FAILURE("failed reading timekeys at " << timekeys_path << ": " << e.what());
    }
}

void Timekeys::initWithBytes(size_t tksize, void* bytes) {
  
    const bool verbose = getenv("DTRPLUGIN_VERBOSE");
    key_prologue_t* prologue = static_cast<key_prologue_t*>(bytes);
    /* check the magic number */
    if (tksize<(sizeof(key_prologue_t))) {
        DTR_FAILURE("timekeys file is too small");
    }
    prologue->magic = htonl(prologue->magic);
    if (prologue->magic != magic_timekey) {
        char buf[100];
        sprintf(buf, "timekeys magic number %x doesn't match %x", 
                          prologue->magic, magic_timekey);
        DTR_FAILURE(buf);
    }
  
    /* get frames per file and key record size */
    prologue->frames_per_file = ntohl( prologue->frames_per_file );
    prologue->key_record_size = ntohl( prologue->key_record_size );
    m_fpf = prologue->frames_per_file;

    if (prologue->key_record_size != sizeof(key_record_t)) {
        /* technically we should be able to handle this; in practice we don't */
        DTR_FAILURE("timekeys record size " << prologue->key_record_size
                << " doesn't match expected size " << sizeof(key_record_t));
    }
  
    /* read all key records */
    size_t nframes = (tksize-sizeof(key_prologue_t))/sizeof(key_record_t);
    size_t nbytes = nframes * sizeof(key_record_t);
    if (nbytes + sizeof(key_prologue_t) != tksize) {
        DTR_FAILURE("Invalid size of timekeys records");
    }
  
    keys.resize(nframes);
    if (nframes>0) {
        memcpy(&keys[0], prologue+1, nbytes);
    }

    /* Check that we didn't get zero-length frames; this would be a strong
     * indicator of file corruption! */
    size_t i;
    for (i=0; i<nframes; i++) {
        if (keys[i].size()==0) {
            DTR_FAILURE("timekeys frame " << i << " had 0 size");
        }

        if ((i > 0) && (keys[i].time() <= keys[(i-1)].time())) {
            DTR_FAILURE("timekeys frame " << i << " had time " << keys[i].time() << " but frame " << (i-1) << " had time " << keys[(i-1)].time());
        }
    }

    m_size = m_fullsize = keys.size();
    if (!keys.size()) return;

    m_first = keys[0].time();
    m_framesize = keys[0].size();
    if (keys.size()==1) {
        m_interval=0;
        keys.clear();
    } else if (keys.size()>1) {
        m_interval=keys[1].time()-keys[0].time();
    }

    int time_interval_warnings = 0;
    for (i=1; i<keys.size(); i++) {

        /* constant frame size */
        if (keys[i].size() != m_framesize) {
            if (verbose) {
                std::cerr << "non-constant framesize at frame " << i << "\n";
                std::cerr << "  size " << keys[i].size() 
                          << " framesize " << m_framesize << "\n\n";
            }
            m_framesize = 0;
            return;
        }
        /* constant offset */
        if (keys[i].offset() != m_framesize*( i % m_fpf)) {
            if (verbose) {
                std::cerr << "unexpected offset for frame " << i << "\n";
            }
            m_framesize = 0;
            return;
        }
        /* constant time interval */
        if (m_interval>0) {
          if (((keys[i].time()-keys[i-1].time())-m_interval) > 1e-3) {
            if (++time_interval_warnings<=10 && verbose)
                printf("non-constant time interval at frame %ld\n", i);
          }
          /* make sure that m_interval precisely equals the given time */
          double computed_time = m_first + i * m_interval;
          if (keys[i].time() != computed_time) {
            /* keep checking for constant framesize, but record that
            * the interval is not constant */
            if (verbose)
              printf("frame %lu time %g != computed %g\n",
                i, keys[i].time(), computed_time);
            m_interval=0;
          }
        }
    }
    /* If we still have good m_interval and m_framesize, we don't need
     * the explicit key records anymore. */
    if (m_interval>0 && m_framesize>0) {
      if (verbose) {
        printf("all times computable from interval %g\n", m_interval);
      }
      keys.clear();
    }
}

key_record_t Timekeys::operator[](uint64_t i) const {
    if (i>=m_fullsize) {
        DTR_FAILURE("Frame index " << i << " must be less than " << m_fullsize);
    }
    if (keys.size()) return keys.at(i);

    key_record_t timekey;
#if defined(_MSC_VER)
    double time = m_first + ((__int64) i)*m_interval;
#else
    double time = m_first + i*m_interval;
#endif
    uint64_t offset = (i % m_fpf) * m_framesize;

    timekey.time_lo = htonl(lobytes(time));
    timekey.time_hi = htonl(hibytes(time));
    timekey.offset_lo = htonl(lobytes(offset));
    timekey.offset_hi = htonl(hibytes(offset));
    timekey.framesize_lo = htonl(lobytes(m_framesize));
    timekey.framesize_hi = htonl(hibytes(m_framesize));
    return timekey;
}

namespace {
    template <typename T> 
    void rawdump(std::ostream& out, const T& v) {
        out.write((char *)&v, sizeof(v));
    }

    template <typename T> 
    void rawload(std::istream& in, T& v) {
        in.read((char *)&v, sizeof(v));
    }
}

void Timekeys::dump(std::ostream& out) const {
    rawdump(out, m_first);
    rawdump(out, m_interval);
    rawdump(out, m_framesize);
    rawdump(out, m_size);
    rawdump(out, m_fullsize);
    rawdump(out, m_fpf);
    rawdump(out, keys.size());
    if (keys.size()) {
        out.write((const char *)&keys[0], keys.size()*sizeof(keys[0]));
    }
}

void Timekeys::load(std::istream& in) {
    size_t sz;
    rawload(in, m_first);
    rawload(in, m_interval);
    rawload(in, m_framesize);
    rawload(in, m_size);
    rawload(in, m_fullsize);
    rawload(in, m_fpf);
    rawload(in, sz);
    if (sz) {
        keys.resize(sz);
        in.read((char *)&keys[0], keys.size()*sizeof(keys[0]));
    }
}

static inline std::string addslash(const std::string& s){
    return (s.rbegin()[0] == '/') ? s : s + "/";
}

#define DD_RELPATH_MAXLEN (9) 
static std::string 
DDreldir(const std::string& fname, int ndir1, int ndir2){
    if( fname.find('/', 0) != std::string::npos ) {
        DTR_FAILURE("filename '" << fname << "' must not contain '/'");
    }

    uint32_t hash = cksum(fname);
    uint32_t d1, d2;
    char answer[DD_RELPATH_MAXLEN];
    if(ndir1 > 0){
        d1 = hash%ndir1;
        if(ndir2 > 0){
            d2 = (hash/ndir1)%ndir2;
            sprintf(answer, "%03x/%03x/", d1, d2);
        }else{
            sprintf(answer, "%03x/", d1);
        }
    }else{
        sprintf(answer, "./");
    }
    return std::string(answer);
}

namespace {
  class DDException : public std::runtime_error{
  public:
      int eno;
      DDException(const std::string &text, int _eno=0) 
      : std::runtime_error(text + strerror(_eno)), eno(_eno){}
  };
}

void DDmkdir(const std::string &dirpath, mode_t mode, int ndir1, int ndir2){
    std::string dpslash(addslash(dirpath));

    mode_t openmode = mode | 0300; // make sure we can write into the directory
    if( mkdir(dpslash.c_str(), openmode) < 0 )
        throw DDException("mkdir", errno);
        
    if( mkdir((dpslash + "not_hashed").c_str(), openmode) < 0 )
        throw DDException("mkdir not_hashed subdirectory", errno);

    FILE *fp = fopen((dpslash + "not_hashed/.ddparams").c_str(), "w");
    if(fp == NULL)
        throw DDException("fopen( .ddparams, \"w\" )", errno);
    if( fprintf(fp, "%d %d\n", ndir1, ndir2) < 0 ){
        fclose(fp);
        throw DDException("fprintf(.ddparams ...)", errno);
    }
    if( fclose(fp) )
        throw DDException("fclose(.ddparams)", errno);

    for(int i=0; i<ndir1; ++i){
        char sub[6];
        sprintf(sub, "%03x/", i);
        std::string dirsub = dpslash + sub;
        {
            if( mkdir(dirsub.c_str(), openmode) < 0 )
                throw DDException("mkdir " + dirsub, errno);
        }
        for(int j=0; j<ndir2; ++j){
            char subsub[6];
            sprintf(subsub, "%03x", j);
            std::string dirsubsub = dirsub + subsub;
            if( mkdir(dirsubsub.c_str(), mode) < 0 ) // NOT openmode!
                throw DDException("mkdir " + dirsubsub, errno);
        }
        if( mode != openmode ){
            // change the mode back to what the user requested now
            // that we're done creating stuff...
            if( chmod(dirsub.c_str(), mode) < 0 )
                throw DDException("chmod " + dirsub, errno);
        }
    }
    if( mode != openmode ){
        // change the mode back to what the user requested now
        // that we're done creating stuff...
        if( chmod(dpslash.c_str(), mode) < 0 )
            throw DDException("chmod " + dpslash, errno);
        if( chmod((dpslash + "not_hashed").c_str(), mode) < 0 )
          throw DDException("chmod " + dpslash + "not_hashed", errno);
    }
}


static void 
DDgetparams(const std::string& dirpath, int *ndir1, int *ndir2) {
  // get ddparams, or assume (0,0) and let the frame file opens fail.
  *ndir1 = *ndir2 = 0;
  std::string dirslash(addslash(dirpath));
  // New convention - .ddparams is in not_hashed/.
  FILE *fp = fopen((dirslash + "not_hashed/.ddparams").c_str(), "r");
  // Allow the old convention of placing .ddparams in the top-level.
  if( fp == NULL && errno == ENOENT ) {
      fp = fopen((dirslash + ".ddparams").c_str(), "r");
  }
  if(fp != NULL) {
    if( fscanf(fp, "%d%d", ndir1, ndir2) != 2 ) 
      DTR_FAILURE("Failed to parse .ddparams");
    if( fclose(fp) ) {
      DTR_FAILURE("Warning: Failed to close .ddparams file: " << strerror(errno));
    }
  }
}

static std::string framefile( const std::string &dtr,
                              size_t frameno, 
                              size_t frames_per_file,
                              int ndir1,
                              int ndir2) {
  unsigned frame_file = frameno / frames_per_file;
  std::ostringstream filename;
  filename << "frame" << std::setfill('0') << std::setw(9)
           << frame_file;
  std::string fname = filename.str();

  std::string fullpath(dtr);
  fullpath += "/";
  fullpath += DDreldir(fname, ndir1, ndir2);
  fullpath += fname;
  return fullpath;
}

uint64_t key_record_t::size() const {
  return assemble64(ntohl(framesize_lo), ntohl(framesize_hi));
}
uint64_t key_record_t::offset() const {
  return assemble64(ntohl(offset_lo), ntohl(offset_hi));
}
double key_record_t::time() const {
  return assembleDouble(ntohl(time_lo), ntohl(time_hi));
}

void DtrReader::read_meta() {
    
    if (metap) {
        return;
    }

    std::string metafile = get_path() + s_sep + "metadata";

    try {
        std::vector<char> file_buffer;
        read_file(metafile, file_buffer);
        metap.reset(new metadata(file_buffer.data(), file_buffer.size(), &jobstep_id));
    } catch (std::exception &e) {
        DTR_FAILURE("unable to read metadata frame " << metafile);
    }
}

bool StkReader::recognizes(const std::string &path) {
      return path.size()>4 && 
             path.substr(path.size()-4)==".stk" &&
             isfile(path);
}

static std::string filename_to_cache_location_v8(std::string const& stk) {
  std::stringstream ss;
  bfs::path fullpath = bfs::system_complete(stk);

  //
  // Lots of chemists create symlinks to stk files in the
  // workdir, but we want to reuse stk cache info for all
  // users regardless of where their respective symlinks
  // live, so traverse all symlinks to get to the real
  // file.
  //
  bfs::path canonicalpath = bfs::canonical(fullpath);

  ss << canonicalpath.parent_path().string() << "/." << canonicalpath.filename().string() << ".cache.0008";
  return ss.str();
}

static std::string filename_to_cache_location_v7(std::string const& stk) {
  std::stringstream ss;
  bfs::path fullpath = bfs::system_complete(stk);
  ss << fullpath.parent_path().string() << "/." << fullpath.filename().string() << ".cache.0007";
  return ss.str();
}
 
static std::string filename_to_old_cache_location(std::string const& stk) {
  const char*        cache_prefix = getenv("MOLFILE_STKCACHE_LOCATION");
  if (!cache_prefix) cache_prefix = "/d/en/cache-0/molfile_stk/";

  std::stringstream ss;
  ss << cache_prefix;
  ss << "0007"; // V7 cache file

  bfs::path fullpath = bfs::system_complete(stk);
  for (bfs::path::iterator it=fullpath.begin(); it!=fullpath.end(); ++it) {
    std::string s = it->string();
    if (s=="/") continue; /* for first component */
    ss << "#" << s;
  }
  ss << ".dump";
  return ss.str();
}

bool StkReader::has_velocities() const {
    if (!framesets.size()) return false;
    return framesets[0]->has_velocities();
}

uint32_t StkReader::natoms() const {
    if (!framesets.size()) return false;
    return framesets[0]->natoms();
}

bool StkReader::read_stk_cache_file(const std::string &cachepath, bool verbose, bool v8) {
    if (verbose) {
        printf("StkReader: cachepath %s\n", cachepath.c_str());
    }

    bool read_cache = false;

    /* attempt to read the cache file */
    std::ifstream file(cachepath.c_str());
    if (file) {
        if (verbose) printf("StkReader: cache file %s found\n", cachepath.c_str());
        bio::filtering_istream in;
        in.push(bio::gzip_decompressor());
        in.push(file);

        try {
            if (v8) {
                load_v8(in);
            } else {
                load_v7(in);
            }
        }
        catch (std::exception& e) {
            if (verbose) {
                printf("StkReader: Reading cache file failed, %s\n", e.what());
            }
            in.setstate(std::ios::failbit);
        }

        if (!in) {
            if (verbose) printf("StkReader: reading cache file failed.\n");
            /* reading failed for some reason.  Clear out anything we
             * might have loaded and start over */
            framesets.clear();
        } else {
            if (verbose) {
                printf("StkReader: reading cache file suceeded.\n");
            }
            read_cache = true;
        }
    }

    return read_cache;
}

void StkReader::init(int* changed) {

    curframeset=0;
    if (changed) {
        *changed = 0;
    }
    bool stale_ddparams = false;
    const bool verbose = getenv("DTRPLUGIN_VERBOSE");
    const bool use_cache = getenv("DESRES_ROOT") && 
                          !getenv("MOLFILE_STKCACHE_DISABLE");

    //
    // In the beginning, desres programmers created frameset.
    // Frameset was without much form and was void of an index.
    //
    // And desres programmers said "let there be an index",
    // and it was so.
    //
    // And desres programmers saw the indices, that it was
    // good, and divided the data into index files and timeseries
    // files, and called the indices "timekeys" and the timeseries
    // files framesXXXXXXXXX.  And it was good.
    //
    // And desres programmers said "let there be a cache of
    // the timekeys" to divide the fast processing of trajectories
    // from the slow - and it was so.  And on the 7th iteration,
    // SERIALIZED_VERSION was set to 0007, and it contained the
    // time indices and the INVMASS from the meta frame.
    //
    // And desres programmers called the cached timekeys data
    // "molfile stk cache files" and placed them in a single
    // large directory "/d/en/cache-0/molfile_stk/".
    //
    // And desres programmers saw that the "molfile stk cache
    // files" were good, and they multiplied, and in the
    // 2015th year, desres programmers said "behold, it is a
    // miracle that one directory can hold over 800K files",
    // and they were moved to a new location next to the ".stk"
    // files (see filename_to_*_cache_location).
    //
    // And desres programmers formed Time Squared Sampling from
    // the bits and the bytes, and breathed into it the ability
    // to create many large text files.  And it was good (mostly).
    //
    // And desres programmers said "it is not good that TSS should
    // alone write large text files over NFS" and I will modify
    // the frame format for DTRs so that it may contain these data
    // sets as well, and I will call them ETRs.  And so it was.
    // 
    // And desres programmers realized that the smaller framed
    // ETRs were not able to use the existing 7th version of the
    // "molfile stk cache" format, and an 8th format would need
    // to be created, that contained the entirety of the meta
    // frame.  And so it was.
    //
    // To not pay a huge penalty on startup with the conversion
    // to v8, find any possible v7 cache file and get the time
    // indices from it, so that only processing of a single meta
    // frame is needed.
    //
    // And as a crude hack, look at the first and last meta frames
    // in the stk, and only read all meta frames if those two are
    // not identical.
    //
    
    bool found_v7_cache = false;
    bool found_v8_cache = false;

    /* use an stk cache if running in DESRES */
    if (use_cache) {

        std::string cachepath = filename_to_cache_location_v8(dtr);
        found_v8_cache = read_stk_cache_file(cachepath, verbose, true);

        if (verbose) {
            if (found_v8_cache) {
                printf("StkReader: found v8 stk cache.\n");
            } else {
                printf("StkReader: no v8 stk cache found.\n");
            }
        }

        if (found_v8_cache == false) {
            //
            // If no v8 file is found, look for a v7 file in 2 locations.
            //
            cachepath = filename_to_cache_location_v7(dtr);
            found_v7_cache = read_stk_cache_file(cachepath, verbose, false);

            if (verbose) {
                if (found_v7_cache) {
                    printf("StkReader: found v7 stk cache at %s\n", cachepath.c_str());
                } else {
                    printf("StkReader: no v7 stk cache found at %s\n", cachepath.c_str());
                }
            }
        }

        if ((found_v8_cache == false) && (found_v7_cache == false)) {
            //
            // One more place to look for an old v7 cache entry.
            //
            cachepath = filename_to_old_cache_location(dtr);
            found_v7_cache = read_stk_cache_file(cachepath, verbose, false);

            if (verbose) {
                if (found_v7_cache) {
                    printf("StkReader: found v7 stk cache at %s\n", cachepath.c_str());
                } else {
                    printf("StkReader: no v7 stk cache found at %s\n", cachepath.c_str());
                }
            }
        }
    }

    /* read stk file */
    std::ifstream input(path().c_str());
    if (!input) {
        DTR_FAILURE("Could not open stk file");
    }
    std::string fname;
    std::vector<std::string> fnames;
    // Prep stk absolute path, in case it might be needed.
    bfs::path dirname = bfs::path(path());
    if (bfs::is_symlink(dirname)) {
        dirname = bfs::read_symlink(dirname);
    }
    dirname = dirname.parent_path();
    while (std::getline(input, fname)) {
        /* strip whitespace?  ignore comment lines? */
        /* Fix for relative paths here */
        if (fname.size()) {
            if (bfs::path(fname).is_absolute()) {
                fnames.push_back(fname);
            } else {               
                fnames.push_back((dirname/bfs::path(fname)).c_str());
            }
        }
    }
    input.close();

    if (!fnames.size()) {
        DTR_FAILURE("Empty stk file");
    }

    /* find which framesets have already been loaded */
    unsigned i=0;
    for (; i<fnames.size(); i++) {

        if (i==framesets.size() || fnames[i]!=framesets[i]->path()) {
            break;
        }

        if (verbose) {
            printf("StkReader: Reusing dtr at %s\n", fnames[i].c_str());
        }

        /* previous versions did lazy loading of the ddparams.  We now
         * wish we hadn't done that because previously cached dtrs will
         * never get updated.  So we force the update here, and check
         * whether anything changed, necessitating a cache update. */
        stale_ddparams |= framesets[i]->update_ddparams();
    }

    /* delete any remaining framesets */
    for (unsigned j=i; j<framesets.size(); j++) {
        delete framesets[j];
    }
    framesets.erase(framesets.begin() + i, framesets.end());

    /* delete the filenames we've already loaded */
    fnames.erase(fnames.begin(), fnames.begin()+i);

    /* The set of overlapping frames may have changed!  Restore the keys
     * to their full, non-overlapping glory.  */
    for (i=0; i<framesets.size(); i++) {
        DtrReader * r = framesets[i];
        r->keys.restore_full_size();
    }

    /* read the unread timekeys files */
    std::vector<Timekeys> timekeys(fnames.size());
    for (unsigned i=0; i<timekeys.size(); i++) {
        if (verbose) {
            printf("StkReader: Loading timekeys from dtr at %s\n", fnames[i].c_str());
        }
        timekeys[i].init(fnames[i]);
    }

    if (changed) {
        *changed = fnames.size();
    }

    append(fnames, timekeys);

    if ((fnames.size() || stale_ddparams || !found_v8_cache) && use_cache) {

        /* update the cache */
        if (verbose) printf("StkReader: updating cache\n");

        /* create temporary file.  We don't care about errors writing
         * the file because we'll detect any error on rename */
        std::string cachepath = filename_to_cache_location_v8(dtr);
        std::string tmpfile(bfs::unique_path(cachepath+"-%%%%-%%%%").string());
        int fd = open(tmpfile.c_str(), O_WRONLY|O_CREAT|O_BINARY, 0666);
        if (fd<0) {
          if (verbose)
            printf("StkReader: warning, creation of cache tmpfile at %s failed: %s\n",
            tmpfile.c_str(), strerror(errno));
        } else {
          bio::filtering_ostream out;
          bio::file_descriptor_sink file(fd, bio::close_handle);
          out.push(bio::gzip_compressor());
          out.push(file);
          dump(out);
        }

        if (fd >= 0) {
            /* do the rename, check for failure */
            boost::system::error_code ec;
            bfs::rename(bfs::path(tmpfile), bfs::path(cachepath), ec);
            if (ec) {
                if (verbose)
                    /* FIXME: this should probably be more noticeable... */
                    printf("StkReader: rename of tmpfile to %s failed: %s\n",
                           cachepath.c_str(), strerror(errno));
            } else {
                if (verbose)
                    printf("StkReader: cache update succeeded.\n");

                //
                // This chmod is needed to deal with umasks that
                // might create the file with more restrictive
                // permissions than 0666 in the open O_CREAT.
                //
                int rc = chmod(cachepath.c_str(), 0666);
                if (verbose) {
                    if (rc == 0) {
                        printf("StkReader: cache file %s successfully chmod to 0666.", cachepath.c_str());
                    } else {
                        printf("StkReader: unable to chmod cache file %s to 0666.", cachepath.c_str());
                    }
                }
            }
        }
    }
}

void StkReader::append(std::vector<std::string> fnames,
                       std::vector<Timekeys> const& timekeys ) {

    uint32_t starting_framesets = framesets.size();

    /* instantiate dtr readers */
    for (unsigned i=0; i<fnames.size(); i++) {
        DtrReader *reader = new DtrReader(fnames[i], _access);
        framesets.push_back(reader);
    }

    process_meta_frames();

    /* intialize dtr readers */
    DtrReader* first = NULL;
    for (unsigned i=0; i<framesets.size(); i++) {
        if (framesets[i]->natoms()>0) {
            first=framesets[i];
            break;
        }
    }
    for (unsigned i=0; i<fnames.size(); i++) {
        auto reader = framesets[i + starting_framesets];
        if (first) {
            reader->set_natoms(first->natoms());
            reader->set_has_velocities(first->has_velocities());
        }

        try {
            reader->initWithTimekeys(timekeys[i]);
        } catch (std::exception &e) {
            delete reader;
            DTR_FAILURE("Failed opening frameset at " << fnames[i + starting_framesets] << ": " << e.what());
        }
        if (first==NULL && reader->natoms()>0) {
            first = reader;
            framesets[0]->set_natoms(first->natoms());
        }
    }

    // now remove overlaps
    while (!framesets.empty() && !framesets.back()->size()) {
        delete framesets.back();
        framesets.pop_back();
    }

    if (framesets.size()) {
        double first=framesets.back()->keys[0].time();
        size_t i=framesets.size()-1;
        while (i--) {
            /* find out how many frames to keep in frameset[i] */
            Timekeys& cur = framesets[i]->keys;
            size_t n = cur.size();
            if (n==0) {
                continue;
            }

            if (first < cur[0].time()) {
                DTR_FAILURE(
                "Frameset " << framesets[i+1]->path() <<
                " in stk " << path() << 
                " has initial time " << first << " which is earlier than any time in the preceding frameset, thus superseding it entirely.  This is probably an erroneously generated stk file.");
            }

            while (n && cur[n-1].time() >= first) {
                --n;
            }

            cur.truncate( n );
            if (cur.size()) {
                double c0t = cur[0].time();
                first = (first < c0t) ? first : c0t;
            }
        }
    }
}

ssize_t StkReader::size() const {
  ssize_t result=0;
  for (size_t i=0; i<framesets.size(); i++) 
    result += framesets[i]->keys.size();
  return result;
}

bool StkReader::next(molfile_timestep_t *ts) {
  int rc=false;
  while (curframeset < framesets.size() && 
         (rc=framesets[curframeset]->next(ts))==false) {
    ++curframeset;
  }
  return rc;
}

const DtrReader * StkReader::component(ssize_t &n) const {
  for (size_t i=0; i<framesets.size(); i++) {
    ssize_t size = framesets[i]->size();
    if (n < size) return framesets[i];
    n -= size;
  }
  return NULL;
}

dtr::KeyMap StkReader::frame(ssize_t n, molfile_timestep_t *ts,
                            void ** bufptr) const {
  const DtrReader *comp = component(n);
  if (!comp) DTR_FAILURE("Bad frame index " << n);
  return comp->frame(n, ts, bufptr);
}

StkReader::~StkReader() {
  for (size_t i=0; i<framesets.size(); i++) 
    delete framesets[i];
}

std::string DtrReader::framefile(ssize_t n) const {
  return ::framefile( dtr, n, framesperfile(), ndir1(), ndir2() );
}

void DtrReader::initWithTimekeys(Timekeys const& tk) {
  keys = tk;
  init_common();
}

void DtrReader::init(int* changed) {
  keys.init(path());
  if (changed) *changed=1;
  init_common();
}

void DtrReader::init_common() {

  bool with_momentum = false;

  /* update the ddparams and meta if we haven't read them already */
  update_ddparams();
  read_meta();

  std::string format;

  if (metap->get_frame_map()->find("FORMAT") != metap->get_frame_map()->end()) {
      Key kv = metap->get_frame_map()->at("FORMAT");
      format = (char *) kv.data;
  }

  if (format != "ETR_V1") {
      // read the first frame to see how many atoms there are, and whether 
      // there are any velocities.
      // Do this only if n_atoms isn't already set
      if (keys.size()>0 && _natoms==0) {
          if (getenv("DTRPLUGIN_VERBOSE")) {
              printf("DtrReader: reading first frame to get atom count\n");
          }
          std::string fname=::framefile(dtr, 0, keys.framesperfile(), ndir1(), ndir2());
          std::vector<char> buffer;


          // If we know the framesize, there is no need to read the entire file, just the
          // first frame's worth.  So, resize buffer such that read_file will only read
          // the first frame.
          if (keys.framesize() > 0) {
              buffer.resize(keys.framesize());
          }

          read_file(fname, buffer);
          bool swap;
          KeyMap blobs = ParseFrame(buffer.size(), &buffer[0], &swap);
          with_momentum = blobs.find("MOMENTUM")!=blobs.end();

          // I'm aware of these sources of atom count: 
          //  "POSN" (the original frameset format)
          //  "POSITION" (the wrapped frameset formats)
          //  "POS" (anton trajectories)
          //  "FORCES" (forces)
          const char *posnames[] = { "POSN", "POSITION", "POS", "FORCES" };
          for (unsigned i=0; i<sizeof(posnames)/sizeof(posnames)[0]; i++) {
              if (blobs.find(posnames[i])!=blobs.end()) {
                  _natoms = blobs[posnames[i]].count / 3;
                  break;
              }
          }
          // similar for velocities
          const char *velnames[] = { "MOMENTUM", "VELOCITY" };
          for (int i=0; i<2; i++) {
              if (blobs.find(velnames[i])!=blobs.end()) {
                  with_velocity=true;
                  break;
              }
          }
      }
  }

  if (with_momentum && (metap->get_frame_map()->find("INVMASS") == metap->get_frame_map()->end())) {
    if (getenv("DTRPLUGIN_VERBOSE")) {
      printf("DtrReader: dtr %s contains MOMENTUM but no RMASS, so no velocities will be read\n", dtr.c_str());
    }
    with_velocity = false;
  }
}

bool DtrReader::update_ddparams() {
  if (m_ndir1<0 || m_ndir2<0) {
      if (getenv("DTRPLUGIN_VERBOSE")) {
          printf("DtrReader: reading ddparams for %s\n", dtr.c_str());
      }
      DDgetparams(dtr, &m_ndir1, &m_ndir2);
      return true;
  }
  return false;
}

ssize_t DtrReader::times(ssize_t start, ssize_t count, double *t) const {
    ssize_t remaining = keys.size()-start;
    count = (count < remaining) ? count : remaining;
    for (ssize_t j=0; j<count; j++) {
        t[j]=keys[start++].time();
    }
    return count;
}

ssize_t StkReader::times(ssize_t start, ssize_t count, double *t) const {
    ssize_t nread=0;
    size_t i=0,n=framesets.size();
    if (start<0) return 0;
    if (count<=0) return 0;
    /* Find the first frameset containing frames in the desired range */
    /* FIXME: could do this using a binary search... */
    for (; i<n; i++) {
        ssize_t sz = framesets[i]->size();
        if (start<sz) break;
        start -= sz;
    }
    /* Read times from framesets until count times are read. */
    for (; i<n; i++) {
        ssize_t sz = framesets[i]->times(start, count, t+nread);
        nread += sz;
        count -= sz;
        start=0;
        if (!count) break;
    }
    return nread;
}

static void read_homebox( const double *box,
                          molfile_timestep_t *ts ) {

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      ts->unit_cell[3*i+j] = box[3*j+i];
    }
  }
}

void write_homebox( const molfile_timestep_t * ts,
                    double * box ) {
  /* our on-disk format stores the box vectors as columns... */
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      box[3*i+j] = ts->unit_cell[3*j+i];
    }
  }
}

static void read_scalars(
    KeyMap &blobs,
    molfile_timestep_t *ts ) {

#if defined(DESRES_READ_TIMESTEP2)
  if (blobs.find("ENERGY")!=blobs.end()) {
      blobs["ENERGY"].get(&ts->total_energy);
  }
  if (blobs.find("POT_ENERGY")!=blobs.end()) {
      blobs["POT_ENERGY"].get(&ts->potential_energy);
  }
  if (blobs.find("POTENTIALENERGY")!=blobs.end()) {
      blobs["POTENTIALENERGY"].get(&ts->potential_energy);
  }
  if (blobs.find("KIN_ENERGY")!=blobs.end()) {
      blobs["KIN_ENERGY"].get(&ts->kinetic_energy);
  }
  if (blobs.find("KINETICENERGY")!=blobs.end()) {
      blobs["KINETICENERGY"].get(&ts->kinetic_energy);
  }
  if (blobs.find("EX_ENERGY")!=blobs.end()) {
      blobs["EX_ENERGY"].get(&ts->extended_energy);
  }
  if (blobs.find("PRESSURE")!=blobs.end()) {
      blobs["PRESSURE"].get(&ts->pressure);
  }
  if (blobs.find("TEMPERATURE")!=blobs.end()) {
      blobs["TEMPERATURE"].get(&ts->temperature);
  }
  if (blobs.find("PRESSURETENSOR")!=blobs.end()) {
      blobs["PRESSURETENSOR"].get(ts->pressure_tensor);
  }
  if (blobs.find("VIRIALTENSOR")!=blobs.end()) {
      blobs["VIRIALTENSOR"].get(ts->virial_tensor);
  }
#endif
}

static void handle_etr_v1(uint32_t len, const void *buf, molfile_timestep_t *ts, dtr::KeyMap meta_blobs, dtr::KeyMap *frame_blobsp, bool swap) {

    //
    // Iteratre through the meta blobs and add entries to the
    // frames map pointing into the correct offset within the
    // single true frame blob named "_D"
    //
    auto blob = frame_blobsp->find("_D");
    if (blob == frame_blobsp->end()) {
        DTR_FAILURE("etr_v1 frame has no _D blob");
    }
    
    char *base_datap = (char *) (blob->second.data);

    for (auto it = meta_blobs.begin(); it != meta_blobs.end(); ++it) {
        uint32_t *blobp = (uint32_t *) it->second.data;
        uint32_t type = blobp[0];
        uint32_t offset = blobp[1];
        uint32_t count = blobp[2];
        void *addr = base_datap + offset;
        (*frame_blobsp)[it->first] = Key(addr, count, type, swap);
    }

    //
    // Now delete the "_D" blob since that is an internal encoding
    // hack that should not be exposed via a keymap to users.
    //
    frame_blobsp->erase("_D");
}

void handle_force_v1(
    KeyMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    molfile_timestep_t *ts ) {

    if (blobs.find("FORCES")==blobs.end()) {
        DTR_FAILURE("Missing FORCES field in frame");
    }
    Key frc=blobs["FORCES"];
    if (frc.count != 3*natoms) {
        DTR_FAILURE("Expected " << 3*natoms  << " elements in FORCES; got " << frc.count);
    }
  if (ts->dcoords) frc.get(ts->dcoords);
  if (ts->coords)  frc.get(ts->coords);

  read_scalars(blobs, ts);
}

static void handle_wrapped_v2(
    KeyMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    molfile_timestep_t *ts ) {

  KeyMap::const_iterator iter;

  // just read POSITION in either single or double precision
  if (blobs.find("POSITION")==blobs.end()) {
      DTR_FAILURE("Missing POSITION field in frame");
  }
  Key pos=blobs["POSITION"];
  if (pos.count != 3*natoms) {
      DTR_FAILURE("Expected " << 3*natoms  << " elements in POSITION; got " << pos.count);
  }
  if (ts->dcoords) pos.get(ts->dcoords);
  if (ts->coords)     pos.get(ts->coords);

  if (with_velocity && (ts->velocities || ts->dvelocities) 
                    && blobs.find("VELOCITY")!=blobs.end()) {
    Key vel=blobs["VELOCITY"];
    if (vel.count != 3*natoms) {
      DTR_FAILURE("Expected " << 3*natoms  << " elements in VELOCITY; got " << vel.count);
    }
    if (ts->dvelocities) vel.get(ts->dvelocities);
    if (ts->velocities)     vel.get(ts->velocities);
  }

  if ((iter=blobs.find("UNITCELL"))!=blobs.end()) {
    double box[9];
    iter->second.get(box);
    read_homebox( box, ts );
  }

  if (ts->gids && (iter=blobs.find("GID"))!=blobs.end()) {
    iter->second.get(ts->gids);
  }

  read_scalars(blobs, ts);
}

namespace {

  inline void
  compute_center(int partition,
                 int nx, int ny, int nz,
                 float b0, float b1, float b2,
                 float b3, float b4, float b5,
                 float b6, float b7, float b8,
                 float* cx, float* cy, float* cz) {
    double nu, nv, nw, mu, mv, mw;
    double xc, yc, zc;

    // -----------------------------------------------
    // Map the partition number to its "mesh" position
    // (see define_mesh_collective in topology.c)
    // -----------------------------------------------
    int hmx = partition;
    int hmy  = hmx / nx;     /* y = y + ny*( z + nz*r ) */
    int hmz  = hmy / ny;     /* z = z + nz*r */
    hmx -= hmy * nx;         /* x = x */
    hmy -= hmz * ny;         /* y = y */

    nu = (double)nx;
    nv = (double)ny;
    nw = (double)nz;

    // -----------------------------------------------
    // Code adapted from configure_global_cell in
    // topology.c
    // -----------------------------------------------
    mu = -0.5*(nu-1) + (double)hmx;
    mv = -0.5*(nv-1) + (double)hmy;
    mw = -0.5*(nw-1) + (double)hmz;

    // We used to do FORCE_PRECISION(xc,float) here, but that
    // seems unnecessary in the context of trajectory writing.
    xc = b0*mu + b1*mv + b2*mw; 
    yc = b3*mu + b4*mv + b5*mw; 
    zc = b6*mu + b7*mv + b8*mw; 

    *cx = xc;
    *cy = yc;
    *cz = zc;
  }

  inline void
  posn_momentum_v_1(int32_t nx, int32_t ny, int32_t nz,
                    uint64_t nparticles,
                    const double  * home_box,
                    const uint32_t* gid,
                    const uint32_t* npp,
                    const float   * rm, // reciprocal mass
                    const float* posn, const float* momentum,
                    /* returns */
                    float *position, float *velocity, double *box) {

    // bounding box is a straight multiple of the home box
    if (box) {
      box[0] = home_box[0]*nx;
      box[1] = home_box[1]*ny;
      box[2] = home_box[2]*nz;
        
      box[3] = home_box[3]*nx;
      box[4] = home_box[4]*ny;
      box[5] = home_box[5]*nz;

      box[6] = home_box[6]*nx;
      box[7] = home_box[7]*ny;
      box[8] = home_box[8]*nz;
    }


    int partition = 0;
    int remaining = 0;
    float cx = 0;
    float cy = 0;
    float cz = 0;
    float ux = home_box[0];
    float vx = home_box[1];
    float wx = home_box[2];
    float uy = home_box[3];
    float vy = home_box[4];
    float wy = home_box[5];
    float uz = home_box[6];
    float vz = home_box[7];
    float wz = home_box[8];

    for(uint64_t i=0; i<nparticles; ++i) {
      if (remaining == 0) {
        do {
          remaining = npp[partition];
          ++partition;
        } while (!remaining); // skip empty partitions
          compute_center(partition-1, nx,ny,nz, ux,vx,wx, uy,vy,wy, uz,vz,wz,
                        &cx,&cy,&cz);
      }
      uint32_t id = gid[i];
      if (id >= nparticles) {
          DTR_FAILURE("non-contiguous particles");
      }

      if (posn) {
        float x = posn[3*i+0];
        float y = posn[3*i+1];
        float z = posn[3*i+2];

        position[3*id+0] = ux*x + vx*y + wx*z + cx;
        position[3*id+1] = uy*x + vy*y + wy*z + cy;
        position[3*id+2] = uz*x + vz*y + wz*z + cz;
      }

      if (velocity && momentum && rm) {
        velocity[3*id+0] = momentum[3*i+0]*rm[id];
        velocity[3*id+1] = momentum[3*i+1]*rm[id];
        velocity[3*id+2] = momentum[3*i+2]*rm[id];
      } else if (velocity) {
        velocity[3*id+0] = 0.0;
        velocity[3*id+1] = 0.0;
        velocity[3*id+2] = 0.0;
      }
      --remaining;
    }
  }
}

static void handle_posn_momentum_v1(
    KeyMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    const float * rmass,
    molfile_timestep_t *ts ) {

  if (ts->dcoords) {
      DTR_FAILURE("Sorry, can't load double precision coordinates from this frameset format");
  }
  int32_t nx, ny, nz;
  double home_box[9], box[9];
  blobs["HOME_BOX"].get(home_box);
  blobs["NX"].get(&nx);
  blobs["NY"].get(&ny);
  blobs["NZ"].get(&nz);
  
  std::vector<uint32_t> gid, npp;
  std::vector<float> pos, mtm;
  Key gidblob=blobs["GID"];
  Key nppblob=blobs["NPP"];
  Key posblob=blobs["POSN"];
  Key mtmblob=blobs["MOMENTUM"];

  if (gidblob.count != natoms) {
      DTR_FAILURE("Missing GID field");
  }
  if (posblob.count != 3*natoms) {
      DTR_FAILURE("Missing POSN field");
  }
  gid.resize(gidblob.count);
  npp.resize(nppblob.count);
  pos.resize(posblob.count);
  mtm.resize(mtmblob.count);

  gidblob.get(&gid[0]);
  nppblob.get(&npp[0]);
  posblob.get(&pos[0]);

  if (rmass && with_velocity) mtmblob.get(&mtm[0]);

  posn_momentum_v_1( nx, ny, nz, natoms, home_box, 
                     &gid[0], &npp[0], rmass,
                     &pos[0],
                     &mtm[0],
                     ts->coords,
                     ts->velocities,
                     box );

  read_homebox( box, ts );
}

static void handle_wrapped_v1(
    KeyMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    molfile_timestep_t *ts ) {

  if (ts->dcoords) {
      DTR_FAILURE("Sorry, can't load double precision coordinates from this frameset format");
  }
  {
    // homebox
    double home_box[9], box[9];
    int32_t nx, ny, nz;
    blobs["HOME_BOX"].get(home_box);
    blobs["NX"].get(&nx);
    blobs["NY"].get(&ny);
    blobs["NZ"].get(&nz);
    box[0] = home_box[0]*nx;
    box[1] = home_box[1]*ny;
    box[2] = home_box[2]*nz;
      
    box[3] = home_box[3]*nx;
    box[4] = home_box[4]*ny;
    box[5] = home_box[5]*nz;
 
    box[6] = home_box[6]*nx;
    box[7] = home_box[7]*ny;
    box[8] = home_box[8]*nz;
    read_homebox( box, ts );
  }

  Key posblob=blobs["POSN"];
  Key velblob=blobs["VELOCITY"];

  // get positions
  if (posblob.count != 3*natoms) {
      DTR_FAILURE("Missing POSN field");
  }
  posblob.get(ts->coords);
  
  // if required, get velocities
  if (ts->velocities && velblob.count > 0) {
    if (velblob.count != 3*natoms) {
      DTR_FAILURE("Expected " << 3*natoms  << " elements in VELOCITY; got " << velblob.count);
    }
    velblob.get(ts->velocities);
  }
}

static void handle_anton_sfxp_v3(
    KeyMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    const float * rmass,
    molfile_timestep_t *ts ) {

  double positionScale=0, momentumScale=0;
  // position scale...
  {
    Key blob = blobs["POSITIONSCALE"];
    if (blob.count != 1) {
        DTR_FAILURE("Missing POSITIONSCALE field");
    }
    blob.get(&positionScale);
  }
  // momentum scale
  if (ts->velocities || ts->dvelocities) {
    if (!rmass) {
        DTR_FAILURE("Cannot read anton_sfxp_v3 velocities without rmass");
    }

    Key blob = blobs["MOMENTUMSCALE"];
    if (blob.count != 1) {
        DTR_FAILURE("Missing MOMENTUMSCALE field");
    }
    blob.get(&momentumScale);
    momentumScale *= PEAKmassInAmu;
  }

  // box
  {
    double box[9] = { 0,0,0, 0,0,0, 0,0,0 };
    uint32_t anton_box[3];
    Key boxblob = blobs["BOX"];
    if (boxblob.count != 3) {
        DTR_FAILURE("Missing BOX field");
    }
    if (boxblob.type==Key::TYPE_UINT32) {
        boxblob.get(anton_box);
    } else if (boxblob.type==Key::TYPE_INT32) {
        int32_t tmp[3];
        boxblob.get(tmp);
        memcpy(anton_box, tmp, sizeof(anton_box));
    } else {
        DTR_FAILURE("BOX field in sfxp_v3 trajectory has unsupported type");
    }

    box[0] = sfxp_ulp32flt(anton_box[0])*positionScale;
    box[4] = sfxp_ulp32flt(anton_box[1])*positionScale;
    box[8] = sfxp_ulp32flt(anton_box[2])*positionScale;
    read_homebox( box, ts );
  }

  // velocities
  std::vector<int32_t> vel;
  if (ts->velocities || ts->dvelocities) {
    Key velblob = blobs["MOMENTUM"];
    if (velblob.count != 3*natoms) {
        DTR_FAILURE("Missing or wrong size MOMENTUM field");
    }
    vel.resize(3*natoms);
    velblob.get(&vel[0]);
  }

  // positions
  std::vector<int32_t> pos(3*natoms);
  {
    Key posblob = blobs.count("POS") ? blobs["POS"]
                                      : blobs["POSITION"]
                                      ;
    if (posblob.count != 3*natoms) {
        DTR_FAILURE("Missing or wrong size POS/POSITION field");
    }
    posblob.get(&pos[0]);
  }
  // convert and read into supplied storage
  if (ts->coords) for (unsigned i=0; i<natoms; i++) {
    ts->coords[3*i  ] = sfxp_ulp32flt(pos[3*i+0])*positionScale;
    ts->coords[3*i+1] = sfxp_ulp32flt(pos[3*i+1])*positionScale;
    ts->coords[3*i+2] = sfxp_ulp32flt(pos[3*i+2])*positionScale;
    if (ts->velocities) {
      const double rm = rmass[i] * momentumScale; // includes PEAKmassInAmu
      ts->velocities[3*i  ] = (float)(rm * sfxp_ulp32flt(vel[3*i  ]));
      ts->velocities[3*i+1] = (float)(rm * sfxp_ulp32flt(vel[3*i+1]));
      ts->velocities[3*i+2] = (float)(rm * sfxp_ulp32flt(vel[3*i+2]));
    }
  }
  if (ts->dcoords) for (unsigned i=0; i<natoms; i++) {
    ts->dcoords[3*i  ] = sfxp_ulp32flt(pos[3*i+0])*positionScale;
    ts->dcoords[3*i+1] = sfxp_ulp32flt(pos[3*i+1])*positionScale;
    ts->dcoords[3*i+2] = sfxp_ulp32flt(pos[3*i+2])*positionScale;
  }
  if (ts->dvelocities) for (unsigned i=0; i<natoms; i++) {
    const double rm = rmass[i] * momentumScale; // includes PEAKmassInAmu
    ts->dvelocities[3*i  ] = (rm * sfxp_ulp32flt(vel[3*i  ]));
    ts->dvelocities[3*i+1] = (rm * sfxp_ulp32flt(vel[3*i+1]));
    ts->dvelocities[3*i+2] = (rm * sfxp_ulp32flt(vel[3*i+2]));
  }
}

bool DtrReader::next(molfile_timestep_t *ts) {

  if (eof()) return false;
  if (!ts) {
    ++m_curframe;
    return true;
  }
  ssize_t iframe = m_curframe;
  ++m_curframe;
  frame(iframe, ts);
  return true;
}

namespace {
    struct FdCloser {
        int _fd;
        FdCloser(int fd) : _fd(fd) {}
        ~FdCloser() { if (_fd>=0) close(_fd); }
    };
}

dtr::KeyMap DtrReader::frame(ssize_t iframe, molfile_timestep_t *ts, void ** bufptr) const {

    if (iframe<0 || ((size_t)iframe)>=keys.full_size()) {
        DTR_FAILURE("dtr " << dtr << " has no frame " << iframe << ": nframes=" << keys.full_size());
    }
    key_record_t key = keys[iframe];
    off_t offset = assemble64( ntohl(key.offset_lo), 
                               ntohl(key.offset_hi) );
    ssize_t framesize = assemble64( ntohl(key.framesize_lo), 
                                    ntohl(key.framesize_hi) );
    if (ts) ts->physical_time = key.time();

    /* use realloc'ed buffer if bufptr is given, otherwise use temporary 
     * space. */
    std::vector<char> tmp;
    void * buffer = NULL;
    if (bufptr) {
        buffer = *bufptr = realloc(*bufptr, framesize);
    } else {
        tmp.resize(framesize);
        buffer = &tmp[0];
    }

    /* get file descriptor for framefile */
    std::string fname=::framefile(dtr, iframe, framesperfile(), ndir1(), ndir2());
    int fd = -1;
    if (_access==RandomAccess) {
        fd = open(fname.c_str(), O_RDONLY|O_BINARY);
    } else { // SequentialAccess
        if (_last_path != fname) {
            if (_last_fd) close(_last_fd);
            _last_fd = 0;
        }
        if (_last_fd<=0) {
            _last_fd = open(fname.c_str(), O_RDONLY|O_BINARY);
            _last_path = fname;
        }
        fd = _last_fd;
    }
    if (fd<0) {
        DTR_FAILURE("Error opening " << fname << ": " << strerror(errno));
    }
    /* ensure that the file descriptor gets closed when we go out of scope,
     * unless it's being cached as _last_fd. */
    FdCloser _(fd==_last_fd ? -1 : fd);

    if (pread(fd, buffer, framesize, offset) != framesize) {
        std::string err = strerror(errno);
        close(fd);
        _last_fd=0;
        DTR_FAILURE("Error reading " << fname << " with offset " << offset << " size " << framesize << ": " << strerror(errno));
    }

    KeyMap map = frame_from_bytes(buffer, framesize, ts);

    if (!bufptr) {
        map.clear();
    }

    return map;
}

KeyMap DtrReader::frame_from_bytes(const void *buf, uint64_t len, 
                                molfile_timestep_t *ts) const {

    bool swap;
    KeyMap blobs = ParseFrame(len, buf, &swap);

    if (ts) {

        const float * rmass = NULL;

        std::string key = "INVMASS";
        if (metap->get_frame_map()->find("INVMASS") != metap->get_frame_map()->end()) {
            Key blob = metap->get_frame_map()->at("INVMASS");
            rmass = (float *) blob.data;
        }

        // Now, dispatch to routines based on format, which can be
        // defined in either the meta frame or the frame.
        std::string format;
        auto p = metap->get_frame_map()->find("FORMAT");
        if (p != metap->get_frame_map()->end()) {
            format += (char *) p->second.data;
        }
        
        if (format == "") {
            format = blobs["FORMAT"].toString();
        }

        if (format=="WRAPPED_V_2" || format == "DBL_WRAPPED_V_2") {
            handle_wrapped_v2(blobs, _natoms, with_velocity, ts);

        } else if (format=="POSN_MOMENTUM_V_1" || format=="DBL_POSN_MOMENTUM_V_1") {
            handle_posn_momentum_v1(blobs, _natoms, with_velocity, rmass, ts);

        } else if (format=="WRAPPED_V_1" || format == "DBL_WRAPPED_V_1") {
            handle_wrapped_v1(blobs, _natoms, with_velocity, ts);

        } else if (format=="ANTON_SFXP_V3") {
            handle_anton_sfxp_v3(blobs, _natoms, with_velocity, rmass, ts);

        } else if (format=="FORCE_V_1" || format == "DBL_FORCE_V_1") {
            handle_force_v1(blobs, _natoms, with_velocity, ts);

        } else if (format=="ETR_V1") {
            handle_etr_v1(len, buf, ts, *metap->get_frame_map(), &blobs, swap);

        } else {
            DTR_FAILURE("can't handle format " << format);
        }
    }

    return blobs;
}

void write_all( int fd, const char * buf, ssize_t count ) {
    while (count) {
        ssize_t n = ::write(fd, buf, count);
        if (n<0) {
            if (errno==EINTR) continue;
            throw std::runtime_error(strerror(errno));
        }
        buf += n;
        count -= n;
    }
}

void DtrWriter::init(const std::string &path, DtrWriter::Mode mode) {

    dtr=path;
    m_directory=path;
    this->mode = mode;
    char cwd[4096];

    while(m_directory.size() > 0 && m_directory[m_directory.size()-1] == s_sep) {
      m_directory.erase(m_directory.size()-1);
    }

    if ( m_directory[0] != s_sep) {
      if (! ::getcwd(cwd,sizeof(cwd))) {
        throw std::runtime_error(strerror(errno));
      }
      m_directory = std::string(cwd) + s_sep + m_directory;
    }

    std::string timekeys_path = dtr + s_sep + "timekeys";

    if (mode==NOCLOBBER && exists(m_directory)) {
        DTR_FAILURE("path '" << m_directory << "' already exists and mode is NOCLOBBER");
    } else if (mode==APPEND && exists(m_directory)) {
        /* read existing timekeys */
        Timekeys tk;
        tk.init(dtr);
        frames_per_file = tk.framesperfile();
        nwritten = tk.size();
        if (nwritten==0) {
            last_time=0;
            framefile_offset = 0;
            frame_fd = 0;
        } else {
            key_record_t last = tk[nwritten-1];
            last_time = last.time();
            framefile_offset = last.offset() + last.size();
            std::string filepath=framefile(dtr, nwritten, frames_per_file,0,0);
            frame_fd = open(filepath.c_str(),O_WRONLY|O_APPEND|O_BINARY,0666);
        }
        timekeys_file = fopen(timekeys_path.c_str(), "a+b");
        if (!timekeys_file) {
            DTR_FAILURE("Opening timekeys failed: " << strerror(errno));
        }

    } else {
        /* start afresh */
        recursivelyRemove(m_directory);
        ::DDmkdir(m_directory,0777,0, 0);

        // craft an empty metadata frame
        {
          std::string metadata_file = m_directory + s_sep + "metadata";
          FILE *fd = fopen(metadata_file.c_str(), "wb");
          if (!fd) DTR_FAILURE("Failed opening metadata frame file: " << 
                  strerror(errno));
          KeyMap map;
          void * buf = NULL;
          size_t framesize = ConstructFrame(map, &buf);
          fwrite( buf, framesize, 1, fd );
          fclose(fd);
          free(buf);
        }

        // start writing timekeys file */
        timekeys_file = fopen( timekeys_path.c_str(), "a+b" );
        if (!timekeys_file) {
            DTR_FAILURE("Opening timekeys failed: " << strerror(errno));
        } else {
          key_prologue_t prologue[1];
          prologue->magic = htonl(magic_timekey);
          prologue->frames_per_file = htonl(frames_per_file);
          prologue->key_record_size = htonl(sizeof(key_record_t));
          fwrite( prologue, sizeof(key_prologue_t), 1, timekeys_file );
          fflush(timekeys_file);
        }
    }
}

void DtrWriter::truncate(double t) {
    rewind(timekeys_file);
    key_prologue_t prologue[1];
    key_record_t record[1];
    if (fread(prologue, sizeof(prologue), 1, timekeys_file)!=1) { 
        DTR_FAILURE("Reading prologue from timekeys file: " << strerror(errno));
    }
    
    for (;;) {
        if (fread(record, sizeof(record), 1, timekeys_file)!=1) {
            if (feof(timekeys_file)) break;
            DTR_FAILURE("Reading a timekeys record: " << strerror(errno));
        }
        if (record->time() <= t) continue;
        /* found a record past the specified time.  rewind to just before
         * this record, and truncate the rest */
        if (fseek(timekeys_file, -sizeof(record), SEEK_CUR)) {
            DTR_FAILURE("Seeking to position for truncation: " 
                    << strerror(errno));
        }
        if (ftruncate(fileno(timekeys_file), ftell(timekeys_file))) {
            DTR_FAILURE("Truncating timekeys file: "
                    << strerror(errno));
        }
        sync();
        break;
    }
}


int DtrWriter::sync() {
    int frc, trc;
    if (timekeys_file) fflush(timekeys_file);
#if defined(_MSC_VER)
    frc = frame_fd>0 ?    _commit(frame_fd) : 0;
    trc = timekeys_file ? _commit(fileno(timekeys_file)) : 0;
#else
    frc = frame_fd>0 ?    fsync(frame_fd) : 0;
    trc = timekeys_file ? fsync(fileno(timekeys_file)) : 0;
#endif
    return !(frc==0 && trc==0);
}

void DtrWriter::next(const molfile_timestep_t *ts) {

    static const char *title = "written by molfile/VMD";
    const double time = ts->physical_time;
    /* require increasing times (DESRESCode#1053) */
    if (last_time != HUGE_VAL && time <= last_time) {
        DTR_FAILURE("framesets require increasing times. previous " << last_time << " current " << time);
    }

    KeyMap map;
    map["TITLE"].set(title, strlen(title));
    map["CHEMICAL_TIME"].set(&time, 1);
    double box[9];
    const char* format=NULL;

    if (mode==FORCE) {
        format = ts->coords ? "FORCE_V_1" : 
                 ts->dcoords ? "DBL_FORCE_V_1" :
                 NULL;
        if (!format) DTR_FAILURE("Either coords or dcoords must be set");

        if (ts->coords) {
            map["FORCES"].set(ts->coords, 3*natoms);
        } else if (ts->dcoords) {
            map["FORCES"].set(ts->dcoords, 3*natoms);
        }

        if (ts->velocities) {
            map["ENERGIES"].set(ts->velocities, 1*natoms);
        } else if (ts->dvelocities) {
            map["ENERGIES"].set(ts->dvelocities, 1*natoms);
        }

    } else {
        format = ts->coords ? "WRAPPED_V_2" 
                            : "DBL_WRAPPED_V_2";
        write_homebox( ts, box );
        map["UNITCELL"].set(box, 9);
        if      (ts->coords)  map["POSITION"].set(ts->coords,  3*natoms);
        else if (ts->dcoords) map["POSITION"].set(ts->dcoords, 3*natoms);

        if      (ts->velocities)  map["VELOCITY"].set(ts->velocities, 3*natoms);
        else if (ts->dvelocities) map["VELOCITY"].set(ts->dvelocities,3*natoms);
    }
    map["FORMAT"].set(format, strlen(format));
    if (ts->gids) map["GID"].set(ts->gids, natoms);

#if defined(DESRES_READ_TIMESTEP2)
    map["ENERGY"].set(&ts->total_energy, 1);
    map["POT_ENERGY"].set(&ts->potential_energy, 1);
    map["KIN_ENERGY"].set(&ts->kinetic_energy, 1);
    map["EX_ENERGY"].set(&ts->extended_energy, 1);
    map["TEMPERATURE"].set(&ts->temperature, 1);
    map["PRESSURE"].set(&ts->pressure, 1);
    map["PRESSURETENSOR"].set(ts->pressure_tensor, 9);
    map["VIRIALTENSOR"].set(ts->virial_tensor, 9);
#endif

    append(time, map);
}

void DtrWriter::append(double time, KeyMap const& map) {

    uint64_t framesize = ConstructFrame(map, &framebuffer);
    uint64_t keys_in_file = nwritten % frames_per_file;

    if (!keys_in_file) {
      if (frame_fd>0) {
          sync();
          ::close(frame_fd);
      }
      framefile_offset = 0;
      std::string filepath=framefile(dtr, nwritten, frames_per_file, 0, 0);
      frame_fd = open(filepath.c_str(),O_WRONLY|O_CREAT|O_BINARY,0666);
      if (frame_fd<0) throw std::runtime_error(strerror(errno));
    }

    // write the data to disk
    write_all( frame_fd, (const char *)framebuffer, framesize );

    // add an entry to the keyfile list
    key_record_t timekey;
    timekey.time_lo = htonl(lobytes(time));
    timekey.time_hi = htonl(hibytes(time));
    timekey.offset_lo = htonl(lobytes(framefile_offset));
    timekey.offset_hi = htonl(hibytes(framefile_offset));
    timekey.framesize_lo = htonl(lobytes(framesize));
    timekey.framesize_hi = htonl(hibytes(framesize));

    if (fwrite(&timekey, sizeof(timekey), 1, timekeys_file)!=1) {
        DTR_FAILURE("Writing timekey failed: " << strerror(errno));
    }

    ++nwritten;
    framefile_offset += framesize;
}

DtrWriter::~DtrWriter() {
  sync();
  if (frame_fd>0) ::close(frame_fd);
  if (timekeys_file) fclose(timekeys_file);
  if (framebuffer) free(framebuffer);
}

/* Write out the size and then the bytes */
std::ostream& operator<<(std::ostream& out, const boost::shared_ptr < metadata > &mp) {
    out << mp->get_frame_size() << ' ';
    if (mp->get_frame_size()) {
        out.write((const char *)mp->get_frame_data(), mp->get_frame_size());
        out << ' ';
    }
    return out;
}

std::istream& operator>>(std::istream& in, boost::shared_ptr < metadata > &mp) {
    char c;
    uint32_t frame_size;
    in >> frame_size;
    in.get(c);
    if (frame_size) {
        std::vector<char> buffer(frame_size);
        in.read((char *)buffer.data(), frame_size);
        mp.reset(new metadata(buffer.data(), frame_size, NULL));
    }

    return in;
}

std::istream& operator>>(std::istream& in, std::vector<float> &inv_mass) {
    uint32_t sz;
    char c;
    in >> sz;
    in.get(c);
    inv_mass.resize(sz);
    if (sz) {
        in.read((char *)inv_mass.data(), sz*sizeof(float));
    }
    return in;
}

std::ostream& DtrReader::dump(std::ostream &out) const {
    out << dtr << ' '
        << _natoms << ' '
        << with_velocity << ' '
        << m_ndir1 << ' '
        << m_ndir2 << ' '
        << metap->get_hash() << ' ';
    
    if (jobstep_id == "") {
        out << "(NULL)";
    } else {
        out << jobstep_id;
    }

    out << ' ';

    /* write raw m_ndir values so that we don't read them from .ddparams
     * if they haven't been read yet */
    keys.dump(out);
    return out;
}

std::istream& DtrReader::load_v7(std::istream &in) {
  char c;
  std::string version;
  bool has_meta;
  bool owns_meta;

  in >> version;

  if (version != "0007") {
      DTR_FAILURE("Bad version string: " << version);
  }

  in >> dtr
     >> _natoms
     >> with_velocity
     >> owns_meta
     >> has_meta;

  if (owns_meta && has_meta) {
      //
      // V7 format had a serialized INVMASS vector of floats.
      // This is now superceded by V8 holding the entire meta
      // frame, but we still need to parse out correctly the
      // INVMASS right now.
      //
      std::vector<float> inv_mass;
      in.get(c);
      in >> inv_mass;
  }

  in >> m_ndir1
     >> m_ndir2;

  in.get(c);
  keys.load(in);

  return in;
}

std::istream& DtrReader::load_v8(std::istream &in) {
  char c;
  uint64_t meta_hash;

  in >> dtr
     >> _natoms
     >> with_velocity
     >> m_ndir1
     >> m_ndir2
     >> meta_hash
     >> jobstep_id;

  in.get(c);
  keys.load(in);
  return in;
}

std::ostream& StkReader::dump(std::ostream &out) const {

    out << dtr << ' ';

    //
    // Build a set of unique meta frame hashes in use.  Since we
    // support adding and deleting entries in the stk, we don't
    // want an orphaned set of meta frames to gunk up the cache
    // file, so make sure to only write out those that are still
    // in use.
    //
    std::set<uint64_t> meta_hashes;
    for (size_t i=0; i<framesets.size(); i++) {
        meta_hashes.insert(framesets[i]->get_meta()->get_hash());
    }

    out << meta_hashes.size() << ' ';
    for (auto hash : meta_hashes) {
        auto metap = meta_data_map.at(hash);
        out << metap;
    }

    out << framesets.size() << ' ';
    for (size_t i=0; i<framesets.size(); i++) {
        out << framesets[i]->get_meta()->get_hash();
        out << ' ';
        framesets[i]->dump(out);
    }

    return out;
}

std::istream& StkReader::load_v8(std::istream &in) {
  uint32_t meta_hashes_count;
  // The old path isn't useful - we already have the path specified by
  // the user, and it doesn't matter if all that happened is that the
  // directory containing the stk was renamed.  One could argue that
  // we shouldn't use the cache file if the stk was moved, on the
  // suspicion that something has happened to the dtr files as well.
  // But if the dtr contents have changed we're screwed anyway because
  // the positions may already be cached by zendo; all we'd do is
  // possibly wind up with time values which are different from the
  // zendo-cached time values, which would be even more confusing.
  std::string old_path;
  in >> old_path;
  char c;
  in.get(c);
  in >> meta_hashes_count;
  in.get(c);

  for (unsigned int i = 0; i < meta_hashes_count; i++) {
      boost::shared_ptr<metadata> mp;
      in >> mp;
      meta_data_map[mp->get_hash()] = mp;
  }

  size_t size; in >> size; framesets.resize(size);
  in.get(c);
  for (size_t i=0; i<framesets.size(); i++) {
    delete framesets[i];

    uint64_t meta_hash;
    in >> meta_hash;
    in.get(c);
    framesets[i] = new DtrReader("<no file>", _access);
    framesets[i]->load_v8(in);
    framesets[i]->set_meta(meta_data_map[meta_hash]);
    if (framesets[0]->natoms()==0 && framesets[i]->natoms()>0) {
        framesets[0]->set_natoms(framesets[i]->natoms());
    }
  }
  return in;
}

void StkReader::process_meta_frames() {

    //
    // Read the first and last meta frames.  If they have
    // the same hash value, then skip reading the meta frame
    // for each dtr in the stk, assuming that all meta frames
    // in between are also identical.
    //
    //
    // First meta frame
    //
    framesets[0]->read_meta();
    meta_data_map[framesets[0]->get_meta()->get_hash()] = framesets[0]->get_meta();

    if (framesets.size() == 1) {
        //
        // This is the first and last meta frame.  Nothing
        // else to do.
        //
        return;
    }

    //
    // Last meta frame
    //
    uint32_t last_index = (framesets.size() - 1);
    framesets[last_index]->read_meta();

    //
    // See if this meta frame is already in the Stk's map.  If not, add it.
    //
    if (meta_data_map.find(framesets[last_index]->get_meta()->get_hash()) == meta_data_map.end()) {
        meta_data_map[framesets[last_index]->get_meta()->get_hash()] = framesets[last_index]->get_meta();
                                                                                     
    } else {
        //
        // If it's already in the Stk's map, then take that pointer
        // instead of the new one, and the new one will automatically
        // get freed.
        //
        framesets[last_index]->set_meta(meta_data_map[framesets[last_index]->get_meta()->get_hash()]);
    }

    if (meta_data_map.size() == 1) {
        //
        // The common case, assume all intermediate meta frame
        // hashes would be identical.
        //
        for (size_t i = 1; i < (framesets.size() - 1); i++) {
            framesets[i]->set_meta(framesets[0]->get_meta());
            framesets[i]->set_jobstep_id("UNREAD");
        }

    } else {

        for (size_t i = 1; i < (framesets.size() - 1); i++) {

            framesets[i]->read_meta();

            //
            // See if this meta frame is already in the Stk's map.  If not, add it.
            //
            if (meta_data_map.find(framesets[i]->get_meta()->get_hash()) == meta_data_map.end()) {
                meta_data_map[framesets[i]->get_meta()->get_hash()] = framesets[i]->get_meta();
            } else {
                //
                // If it's already in the Stk's map, then take that pointer
                // instead of the new one, and the new one will automatically
                // get freed.
                //
                framesets[i]->set_meta(meta_data_map[framesets[i]->get_meta()->get_hash()]);
            }
        }
    }
}

std::istream& StkReader::load_v7(std::istream &in) {

  std::string old_path;
  in >> old_path;
  size_t size; in >> size; framesets.resize(size);
  char c; in.get(c);
  for (size_t i=0; i<framesets.size(); i++) {
    delete framesets[i];
    framesets[i] = new DtrReader("<no file>", _access);
    framesets[i]->load_v7(in);
    if (framesets[0]->natoms()==0 && framesets[i]->natoms()>0) {
        framesets[0]->set_natoms(framesets[i]->natoms());
    }
  }

  process_meta_frames();

  return in;
}

///////////////////////////////////////////////////////////////////
//
// Plugin Interface
//
// ////////////////////////////////////////////////////////////////

static void *
open_file_read( const char *filename, const char *filetype, int *natoms ) {

  FrameSetReader *h = NULL;

  // check for .stk file
  if (StkReader::recognizes(filename)) {
    h = new StkReader(filename);

  } else {
    // check for "clickme.dtr"
      std::string fname=filename;
    std::string::size_type pos = fname.rfind( "clickme.dtr" );
    if (pos != std::string::npos) {
      fname.resize(pos);
      filename = fname.c_str();
    }
    h = new DtrReader(fname);
  }

  try {
      h->init();
  }
  catch (std::exception &e) {
      delete h;
      fprintf(stderr, "Failed to open %s file at %s: %s\n", 
              filetype, filename, e.what());
      return NULL;
  }
  *natoms = h->natoms();
  return h;
}

static int 
read_timestep_metadata(void *v, molfile_timestep_metadata *m) {
  FrameSetReader* handle = reinterpret_cast<FrameSetReader*>(v);
  m->has_velocities = handle->has_velocities();
  m->count = handle->size();
  m->supports_double_precision = true;
  return MOLFILE_SUCCESS;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  try {
    if (h->next(ts)) return MOLFILE_SUCCESS;
  }
  catch (std::exception &e) {
      fprintf(stderr, "Failed reading timestep: %s\n", e.what());
  }
  return MOLFILE_ERROR;
}

#if defined(DESRES_READ_TIMESTEP2)
static int read_timestep2(void *v, ssize_t n, molfile_timestep_t *ts) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  try {
    h->frame(n, ts);
    return MOLFILE_SUCCESS;
  }
  catch (std::exception &e) {
      fprintf(stderr, "Failed reading timestep: %s\n", e.what());
  }
  return MOLFILE_ERROR;
}

static molfile_ssize_t read_times(void *v, 
                                  molfile_ssize_t start, 
                                  molfile_ssize_t count,
                                  double * times) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  return h->times(start, count, times);
}
#endif

static void close_file_read( void *v ) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  delete h;
}

static void *open_file_write(const char *path, const char *type, int natoms) {
  DtrWriter::Mode mode;
  if (!strcmp(type, "dtr") || !strcmp(type, "dtr_clobber")) {
      mode = DtrWriter::CLOBBER;
  } else if (!strcmp(type, "dtr_append")) {
      mode = DtrWriter::APPEND;
  } else if (!strcmp(type, "dtr_noclobber")) {
      mode = DtrWriter::NOCLOBBER;
  } else if (!strcmp(type, "force_dtr")) {
      mode = DtrWriter::FORCE;
  } else {
      fprintf(stderr, "Unrecognized filetype '%s'\n", type);
      return NULL;
  }

  DtrWriter *h = new DtrWriter(natoms);

  try {
      h->init(path, mode);
  }
  catch (std::exception &e) {
      delete h;
      h = NULL;
      fprintf(stderr, "Failed to open file of type %s at %s for writing: %s\n",
              type, path, e.what());
  }
  return h;
}

static int write_timestep(void *v, const molfile_timestep_t *ts) {
  DtrWriter *h = reinterpret_cast<DtrWriter *>(v);
  try {
      h->next(ts);
      return MOLFILE_SUCCESS;
  }
  catch (std::exception &e) {
      fprintf(stderr, "Failed writing timestep: %s\n", e.what());
  }
  return MOLFILE_ERROR;
}

static void close_file_write( void * v ) {
  DtrWriter *h = reinterpret_cast<DtrWriter *>(v);
  delete h;
}

static int sync_file_write( void * v ) {
  DtrWriter *h = reinterpret_cast<DtrWriter *>(v);
  if (h) return h->sync();
  return MOLFILE_ERROR;
}

static int truncate_file_write(void *v, double t) {
  DtrWriter *h = reinterpret_cast<DtrWriter *>(v);
  try {
      h->truncate(t);
      return MOLFILE_SUCCESS;
  }
  catch (std::exception &e) {
      fprintf(stderr, "Failed truncating file: %s\n", e.what());
  }
  return MOLFILE_ERROR;
}

static molfile_plugin_t desmond;
static molfile_plugin_t dtr_append;
static molfile_plugin_t dtr_clobber;
static molfile_plugin_t dtr_noclobber;
static molfile_plugin_t force_dtr;

extern "C"
int msys_dtrplugin_init() {
  /* Plugin for desmond trajectory files */
  ::memset(&desmond,0,sizeof(desmond));
  desmond.abiversion = vmdplugin_ABIVERSION;
  desmond.type = MOLFILE_PLUGIN_TYPE;
  desmond.name = "dtr";
  desmond.prettyname = "DESRES Trajectory, clobber";
  desmond.author = "D.E. Shaw Research";
  desmond.majorv = 4;
  desmond.minorv = 0;
  desmond.is_reentrant = VMDPLUGIN_THREADUNSAFE;

  desmond.open_file_read = open_file_read;
  desmond.read_timestep_metadata = read_timestep_metadata;
  desmond.read_next_timestep = read_next_timestep;
  desmond.read_timestep2 = read_timestep2;
  desmond.read_times = read_times;
  desmond.close_file_read = close_file_read;

  desmond.open_file_write = open_file_write;
  desmond.write_timestep = write_timestep;
  desmond.close_file_write = close_file_write;
  desmond.sync_file_write = sync_file_write;

  memcpy(&dtr_append, &desmond, sizeof(desmond));
  memcpy(&dtr_clobber, &desmond, sizeof(desmond));
  memcpy(&dtr_noclobber, &desmond, sizeof(desmond));
  memcpy(&force_dtr, &desmond, sizeof(desmond));
  dtr_append.name="dtr_append";
  dtr_append.prettyname = "DESRES Trajectory, append";
  dtr_append.truncate_file_write = truncate_file_write;
  dtr_clobber.name="dtr_clobber";
  dtr_noclobber.name="dtr_noclobber";
  dtr_noclobber.prettyname = "DESRES Trajectory, no-clobber";

  // Make sure that the regular, clobbering dtr plugin is the one that gets 
  // returned by molfile.guess_filetype by initializing this string only
  // after the other plugins have been copied from the desmond template.
  desmond.filename_extension = "dtr,dtr/,stk,atr,atr/,etr,etr/";

  force_dtr.name="force_dtr";
  force_dtr.prettyname = "DESRES Force Trajectory";

  return VMDPLUGIN_SUCCESS;
}

extern "C"
int msys_dtrplugin_register(void *v, vmdplugin_register_cb cb) {
  cb(v,reinterpret_cast<vmdplugin_t*>(&desmond));
  cb(v,reinterpret_cast<vmdplugin_t*>(&dtr_append));
  cb(v,reinterpret_cast<vmdplugin_t*>(&dtr_clobber));
  cb(v,reinterpret_cast<vmdplugin_t*>(&dtr_noclobber));
  cb(v,reinterpret_cast<vmdplugin_t*>(&force_dtr));
  return VMDPLUGIN_SUCCESS;
}

