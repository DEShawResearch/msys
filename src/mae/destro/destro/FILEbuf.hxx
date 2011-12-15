/* @COPYRIGHT@ */

#ifndef _desSTL__FILEbuf_dot_hxx
#define _desSTL__FILEbuf_dot_hxx

#include <iostream>
#include <streambuf>
#include <cstdio>

namespace desres { namespace msys {


  /*!
    \brief Adapt FILE* and file descriptors to a stream.

    Sometimes you want a streambuf that is synced with a FILE*.  With
    such a thing, you can easily make an iostream that is synced with a
    FILE*.  desSTL::FILEbuf is such a streambuf.  Use it like this:

    <pre>
    FILE *fp;
    FILEbuf ibuf(FILEbuf(fp));
    istream in(&ibuf);

    int fd;
    FILEbuf obuf(FILEbuf(fd, "w"));
    ostream out(&obuf);
    </pre>

    This implementation has no state of its own.  All buffering is
    handled by the underlying FILE*.

    Constructors are provided for FILE* and integer file descriptors.
    For a FILE, all I/O is synced with the original FILE.  The FILE is
    not closed when the FILEbuf is destroyed.
  
    In the case of the file descriptor, the argument is dup'ed, and a new
    FILE* is created with fdopen.  The FILEbuf is associated with this
    FILE/fd, which are both closed when the FILEbuf is destroyed.  Since
    we dup'ed the original fd, it remains open.  I/O is not "sync"ed with
    the original fd, which may be considered either a bug or a feature.

    See Matt Austern's article "The Standard Librarian: IOStreams and Stdio"
    http://www.ddj.com/cpp/184401305?pgno=1

    TODO
    - add seekpos and seekoff to facilitate repositioning.
    - add iFILEstream, oFILEstream and FILEstream wrappers
    by analogy with ifstream, ofstream and fstream.

    @param fp_ File pointer to adapt to a buffer

  */
  class FILEbuf : public std::streambuf{
  private:
    FILE *fp;
    bool owned;

  public:
    /*!
      \brief Adapt a FILE* or filedescriptor to an std::stream

    */
    FILEbuf(FILE* fp_) : std::streambuf(), fp(fp_), owned(false) {}

    /*!
     * \brief Destructor
     */
    virtual ~FILEbuf(){
      if(owned)
        fclose(fp); /*GCOV-IGNORE*/
    }

  protected:
    /*!
     * \brief Put back a character on overflow
     */
    virtual int overflow(int c=EOF){
      return c!=EOF ? fputc(c, fp) : EOF;
    }

    /*!
     * \brief get the next char *without* advancing the get pointer.
     *
     * @return the next character in the stream
     */
    virtual int underflow(){ /*GCOV-IGNORE*/
      int c = getc(fp); /*GCOV-IGNORE*/
      if(c!=EOF) /*GCOV-IGNORE*/
        ungetc(c, fp); /*GCOV-IGNORE*/
      return c; /*GCOV-IGNORE*/
    }
    
    /*!
     * \brief get the next char *and* advance the get pointer.
     *
     * @return the next character in the stream
     */
    virtual int uflow(){
      return getc(fp);
    }


    /*!
     * \brief pbackfail 
     */
    virtual int pbackfail(int c=EOF){ /*GCOV-IGNORE*/
      return c!=EOF ? ungetc(c, fp) : EOF; /*GCOV-IGNORE*/
    }

    /*!
     * \brief flush to disk
     *
     * @return 0 on success otherwise EOF
     */
    virtual int sync(){ return fflush(fp); } /*GCOV-IGNORE*/
  };
}}

#endif
