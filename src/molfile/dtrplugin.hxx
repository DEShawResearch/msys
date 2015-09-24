//
// Version info for VMD plugin tree:
//   $Id$
//
// Version info for last sync with D. E. Shaw Research:
//  //depot/desrad/main/sw/libs/vmd_plugins,DIST/dtrplugin.cxx#3
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

#ifndef MOLFILE_DTRPLUGIN_HXX
#define MOLFILE_DTRPLUGIN_HXX

#if defined(_MSC_VER)
#ifndef DESRES_WIN32
#define DESRES_WIN32
#endif
#endif

#include <math.h>
#include <stdio.h>
#ifdef DESRES_WIN32
#include <io.h>
#include <direct.h>
#include <fcntl.h>
#include <windows.h>

typedef int int32_t;
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
#if 1
typedef unsigned __int64 uint64_t;    // This also works with MVSC6
#else
typedef unsigned long long uint64_t;
#endif
typedef unsigned short uint16_t;
typedef unsigned int ssize_t;
typedef int mode_t;
#define mkdir(a,b) _mkdir(a)
#define rmdir(a)   _rmdir(a)
#define ftello(a)  ftell(a)
#else
#define O_BINARY 0
#include <inttypes.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#endif

#include "molfile_plugin.h"

#include <vector>
#include <string>
#include <stdexcept>

#include "dtrframe.hxx"
#include <boost/shared_ptr.hpp>
#include <ThreeRoe/ThreeRoe.hpp>

namespace desres { namespace molfile {
  
  const char * dtr_serialized_version();

  struct key_record_t {
    uint32_t time_lo;       /* Time associated with this file (low bytes). */
    uint32_t time_hi;       /* Time associated with this file (high bytes). */
    uint32_t offset_lo;     /* Zero in the 1 frame/file case (low bytes) */
    uint32_t offset_hi;     /* Zero in the 1 frame/file case (high bytes) */
    uint32_t framesize_lo;  /* Number of bytes in frame (low bytes) */
    uint32_t framesize_hi;  /* Number of bytes in frame (high bytes) */

    double time() const;
    uint64_t offset() const;
    uint64_t size() const;

  };

  class Timekeys {

      double    m_first;        /* last time in the key list */
      double    m_interval;     /* representative time interval between keys */
      uint64_t  m_framesize;    /* size of a frames */
      size_t    m_size;         /* number of non-overlapping frames */
      size_t    m_fullsize;     /* total number of frames */
      uint32_t  m_fpf;          /* frames per file */

      /* storage for keys until compressed, or if not compressible */
      std::vector<key_record_t> keys;

    public:
      Timekeys() 
      : m_first(0), m_interval(0), m_framesize(0), 
        m_size(0), m_fullsize(0), m_fpf(0) {}

      void init( const std::string& path );

      /* initialize from timekeys bytes.  Writable but owned by caller */
      void initWithBytes( size_t nbytes, void* tkbytes );

      uint32_t framesperfile() const { return m_fpf; }

      /* size of a frame, 0 if the size varies between frames, or
       * if there are no frames in the dtr */
      uint64_t framesize() const { return m_framesize; }

      /* interval between frames, 0 if the interval varies between frames
       * or if there are fewer than 2 frames in the dtr */
      double interval() const { return m_interval; }

      bool is_compact() const { return keys.size()==0; }

      size_t size() const { return m_size; }
      size_t full_size() const { return m_fullsize; }

      void truncate( size_t nframes ) { m_size = nframes; }

      void restore_full_size() { m_size = m_fullsize; }

      key_record_t operator[](uint64_t i) const;

      void dump( std::ostream& out ) const;
      void load( std::istream& in  );
  };

  class DtrReader;

  class FrameSetReader {
  protected:
    std::string dtr;

  public:
    virtual ~FrameSetReader() {}

    virtual bool has_velocities() const = 0;

    virtual uint32_t natoms() const = 0;

    const std::string &path() const { return dtr; }

    // initialize all members from frameset at given path.  
    // If changed is provided, set to true/false if timekeys were/were not
    // reloaded.  Throws on failure.
    virtual void init(int * changed = NULL) = 0;

    // number of frames
    virtual ssize_t size() const = 0;
    virtual ssize_t total_bytes() const = 0;

    // read the next frame.  If ts is NULL, skip it.  If there are no more
    // frames to read, return MOLFILE_EOF.  Otherwise, fill in coordinates
    // in ts->coords, velocities in ts->velocities if ts->velocities is 
    // non-NULL, and return MOLFILE_SUCCESS if all goes well.
    virtual bool next(molfile_timestep_t *ts) = 0;

    // Get the DtrReader component corresponding to frame n.  Change
    // n to the local index within the returned dtr.
    virtual const DtrReader * component(ssize_t &n) const = 0;

    // number of framesets
    virtual ssize_t nframesets() const = 0;

    // nth frameset 
    virtual const DtrReader * frameset(ssize_t n) const = 0;

    // read a specified frame.  If a pointer to a buffer is supplied, it 
    // will be realloc'ed to a size sufficient to hold the frame contents, 
    // and a // KeyMap will be returned pointing into the supplied buffer.  
    // If no buffer pointer is supplied, only molfile_timestep_t information 
    // will be filled in, and the returned KeyMap will be empty.
    virtual dtr::KeyMap frame(ssize_t n, molfile_timestep_t *ts,
                              void ** bufptr = NULL) const = 0;

    // read up to count times beginning at index start into the provided space;
    // return the number of times actually read.
    virtual ssize_t times(ssize_t start, ssize_t count, double * times) const = 0;
  };

  class metadata {
  public:
      metadata(const void *bufptr, ssize_t n, std::string *jobstep_id = NULL) {
          
          //
          // We want to put metadata objects into the stack
          // cache file, but that means we only want to copy
          // them in once.  Unfortunately, each meta frame
          // has a JOBSTEP_ID entry, which changes with each
          // new dtr.  So, we now create a copy of everything
          // in the meta frame except for key JOBSTEP_ID tags.
          // 
          // If jobstep_id is not NULL, then the value of the
          // key JOBSTEP_ID is stored there.
          //
          bool swap;
          auto full_frame_map = dtr::ParseFrame(n, bufptr, &swap);

          frame_data = calloc(n, 1);
          memcpy(frame_data, bufptr, n);
          frame_size = n;

          ThreeRoe hasher;

          for (auto kv : full_frame_map) {

              //
              // Iterate through the keys, copying all non JOB*
              // keys, updating the data pointers to the newly
              // copied data region, and computing a ThreeRoe
              // hash of the data as we go.
              //
              uint32_t v_size = kv.second.get_element_size();
              if (kv.first != "JOBSTEP_ID") {
                  void *copied_data_address = (void *) (((char *)kv.second.data - (char *)bufptr) +
                                                        (char *) frame_data);
                  frame_map[kv.first] = dtr::Key(copied_data_address, kv.second.count, kv.second.type, swap);
                  hasher.Update(kv.second.data, v_size);
              } else {
                  if (jobstep_id) {
                      (*jobstep_id) += (char *) kv.second.data;
                  }
              }
          }

          hash = hasher.Final().first;
      }

      ~metadata() {
          if (frame_size) {
              free(frame_data);
          }
      }

      uint64_t get_hash() const {
          return(hash);
      }

      dtr::KeyMap *get_frame_map() {
          return &frame_map;
      }
 
      uint32_t get_frame_size() const {
          return frame_size;
      }

      void *get_frame_data() {
          return frame_data;
      }

  private:
      uint32_t frame_size;
      void *frame_data;
      uint64_t hash;
      dtr::KeyMap frame_map;
  };

  class DtrReader : public FrameSetReader {
    uint32_t _natoms;
    bool with_velocity;

    int m_ndir1;
    int m_ndir2;
    ssize_t m_curframe;

    boost::shared_ptr < metadata > metap;

    bool eof() const { return m_curframe >= (ssize_t)keys.size(); }

    void init_common();
    const unsigned _access;

    // if doing SequentialAccess, cache the last used file descriptor
    // and file path.
    mutable int _last_fd;
    mutable std::string _last_path;

    std::string jobstep_id;

  public:
    enum {
        RandomAccess
      , SequentialAccess    /* WARNING: MAKES frame() NOT REENTRANT */
    };

    // initializing 
    DtrReader(std::string const& path, unsigned access = RandomAccess) 
    : _natoms(0), with_velocity(false), m_ndir1(-1), m_ndir2(-1), m_curframe(0),
      _access(access), _last_fd(0), _last_path("")
    { 
        dtr = path;
    }

    void set_meta(boost::shared_ptr < metadata > p) {
        metap = p;
    }

    void set_jobstep_id(std::string id) {
        jobstep_id = id;
    }

    boost::shared_ptr <metadata> get_meta() {
        return metap;
    }

    virtual ~DtrReader() {
      if (_last_fd>0) close(_last_fd);
    }

    Timekeys keys;

    // DDparams
    int ndir1() const { return m_ndir1; }
    int ndir2() const { return m_ndir2; }

    bool has_velocities() const { return with_velocity; }
    uint32_t natoms() const { return _natoms; }

    /* set by StkReader */
    void set_has_velocities(bool b) { with_velocity = b; }
    void set_natoms(uint32_t n) { _natoms = n; }

    /* update ddparams.  Return true if they were stale and updated.  */
    bool update_ddparams();

    void read_meta();

    ssize_t curframe() const { return m_curframe; }

    uint32_t framesperfile() const { return keys.framesperfile(); }

    void initWithTimekeys(Timekeys const& tk);

    virtual void init(int * changed=NULL);
    virtual ssize_t size() const { return keys.size(); }
    virtual ssize_t total_bytes() const {
        ssize_t total = 0;

        for (unsigned int i = 0; i < keys.size(); i++) {
            total += keys[i].size();
        }

        return(total);
    }

    virtual ssize_t times(ssize_t start, ssize_t count, double * times) const;

    virtual bool next(molfile_timestep_t *ts);

    virtual const DtrReader * component(ssize_t &n) const {
      return this;
    }

    virtual ssize_t nframesets() const { return 1; }
    virtual const DtrReader * frameset(ssize_t n) const {
        if (n!=0) throw std::runtime_error("bad index");
        return this;
    }

    /* WARNING: this method is reentrant only when using RandomAccess */
    virtual dtr::KeyMap frame(ssize_t n, molfile_timestep_t *ts,
                              void ** bufptr = NULL) const;

    // path for frame at index.  Empty string on not found.
    std::string framefile(ssize_t n) const;

    // parse a frame from supplied bytes
    dtr::KeyMap frame_from_bytes( const void *buf, uint64_t len,
                             molfile_timestep_t *ts ) const;

    std::ostream& dump(std::ostream &out) const;
    std::istream& load_v7(std::istream &in);
    std::istream& load_v8(std::istream &in);
    const std::string get_path() {
        return dtr;
    }
  };

  struct DtrWriter {
    std::string dtr;
    std::string m_directory;
    const uint32_t natoms;
    int frame_fd;        // for writing
    uint32_t frames_per_file;
    uint64_t framefile_offset;
    uint64_t nwritten;
    double last_time;
    FILE * timekeys_file;
    int mode;
    void* framebuffer;

    explicit DtrWriter(uint32_t natoms_, uint32_t fpf = 0) 
    : natoms(natoms_), frame_fd(0), framefile_offset(0),
      nwritten(0), last_time(HUGE_VAL), timekeys_file(NULL),
      framebuffer()
    {
        if (fpf > 0) {
            frames_per_file = fpf;
        } else if (natoms==0) {
            frames_per_file = 1024; /* arbitrary */
        } else {
            /* try to achieve 100MB frame files.  Unfortunately for those 
             * trying to write just a few atoms per frame, the frameset 
             * format requires a frame size of no less than 4k. */
            static const uint32_t maxsize = 100 * 1024 * 1024;
            static const uint32_t maxcount = maxsize / 4096;
            /* we assume single precision and no velocities.  If we end
             * up with 200MB frames, that's ok. */
            frames_per_file = maxsize / (12 * natoms);
            if (frames_per_file < 1) frames_per_file = 1;
            if (frames_per_file > maxcount) frames_per_file = maxcount;
        }
    }

    ~DtrWriter();

    enum Mode {
        CLOBBER, 
        NOCLOBBER, 
        APPEND,
        FORCE
    };

    // initialize for writing at path
    void init(const std::string &path, Mode mode=CLOBBER);

    // write another frame.  
    void next(const molfile_timestep_t *ts);

    // write an arbitrary set of keyvals
    void append(double time, dtr::KeyMap const& keyvals);

    // commit timekeys current frame file to disk.  0 on success.
    int sync();

    // remove timekeys with times strictly greater than the given time
    void truncate(double after_time);
  };

  class StkReader : public FrameSetReader {
    std::vector<DtrReader*> framesets;
    size_t curframeset;
    const unsigned _access;
    std::map<uint64_t, boost::shared_ptr < metadata > > meta_data_map;

  public:
    explicit StkReader(std::string const& path, 
                       unsigned access = DtrReader::RandomAccess) 
    : curframeset(0), _access(access) 
    { dtr = path; 
    }
    ~StkReader();

    void append(std::vector<std::string> fnames, 
                std::vector<Timekeys> const& timekeys);

    virtual bool has_velocities() const;
    virtual uint32_t natoms() const;

    virtual void init(int * changed=NULL);
    virtual ssize_t size() const;
    virtual ssize_t total_bytes() const { return 0; };
    virtual ssize_t times(ssize_t start, ssize_t count, double * times) const;
    virtual bool next(molfile_timestep_t *ts);
    virtual dtr::KeyMap frame(ssize_t n, molfile_timestep_t *ts,
                              void ** bufptr = NULL) const;

    virtual const DtrReader * component(ssize_t &n) const;

    virtual ssize_t nframesets() const { return framesets.size(); }
    virtual const DtrReader * frameset(ssize_t n) const {
        return framesets.at(n);
    }

    static bool recognizes(const std::string &path);

    bool read_stk_cache_file(const std::string &cachepath, bool verbose, bool v8);

    std::ostream& dump(std::ostream &out) const;
    std::istream& load_v7(std::istream &in);
    std::istream& load_v8(std::istream &in);
    void process_meta_frames();
  };
} }

#endif
