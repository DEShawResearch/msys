#include "dms.hxx"
#include "../istream.hxx"
#include <sqlite3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <cerrno>
#include <iomanip>
#include <unistd.h>

#include <sys/stat.h>
#include <fcntl.h>

using namespace desres::msys;

namespace {
    struct dms_file : sqlite3_file {
        char * contents;
        sqlite3_int64 size;
        sqlite3_int64 capacity;
        char * path;

        void write() const {
            if (!path) return;
            /* if path is non-NULL, we need to write the contents */
            int fd=open(path, O_WRONLY | O_CREAT, 0644);
            if (fd<0) {
                MSYS_FAIL("Error opening " << path 
                    << " for writing: " << strerror(errno));
            }
            sqlite3_int64 sz = size;
            char* ptr = contents;
            while (sz) {
                errno = 0;
                ssize_t rc = ::write(fd, ptr, sz);
                if (rc<0 || (rc==0 && errno!=0)) {
                    std::string errmsg = strerror(errno);
                    close(fd);
                    MSYS_FAIL("Error writing DMS contents of size " << size
                        << " bytes to " << path << " : " << errmsg);
                }
                sz -= rc;
                ptr += rc;
            }
            if (close(fd)!=0) {
                std::string errmsg = strerror(errno);
                MSYS_FAIL("Error closing DMS file at " << path
                        << " : " << errmsg);
            }
        }
    };
}

namespace {
    int dms_xClose(sqlite3_file *file) {
        dms_file* dms = static_cast<dms_file*>(file);
        if (dms->contents) {
            free(dms->contents);
            dms->contents = NULL;
        }
        if (dms->path) {
            free(dms->path);
            dms->path = NULL;
        }
        return SQLITE_OK;
    }

    int dms_xRead(sqlite3_file *file, void *pBuf, int iAmt, sqlite3_int64 offset) {
        dms_file *dms = (dms_file *)file;
        const char *ptr = dms->contents;
        sqlite3_int64 size=dms->size;

        sqlite3_int64 max=size-offset;  // max allowable read
        if (max<0) max=0;
        sqlite3_int64 amt=iAmt < max ? iAmt : max;
        memcpy(pBuf, ptr+offset, amt);
        if (amt < max) return SQLITE_OK;
        /* Unread parts of the buffer must be zero-filled */
        memset(&((char *)pBuf)[amt], 0, amt-max);
        return SQLITE_IOERR_SHORT_READ;
    }

    int dms_xWrite(sqlite3_file*file, const void*pBuf, int iAmt, sqlite3_int64 offset) {
        dms_file *dms = (dms_file *)file;
        sqlite3_int64 last=offset+iAmt;
        if (dms->capacity < last) {
            dms->capacity = last * 1.5;
            dms->contents = (char *)realloc(dms->contents, dms->capacity);
        }
        dms->size = std::max(dms->size, offset + iAmt);
        memcpy(dms->contents+offset, pBuf, iAmt);
        return SQLITE_OK;
    }

    int dms_xLock(sqlite3_file*file, int lockType) {
        return SQLITE_OK;
    }
    int dms_xUnlock(sqlite3_file*file, int) {
        return SQLITE_OK;
    }

    int dms_xSync(sqlite3_file* file, int flags) {
        return SQLITE_OK;
    }

    int dms_xFileSize(sqlite3_file *file, sqlite3_int64 *pSize) {
        dms_file *dms = (dms_file *)file;
        *pSize = dms->size;
        return SQLITE_OK;
    }

    int dms_xDeviceCharacteristics(sqlite3_file *file) {
        return 0;
    }
    int dms_xFileControl(sqlite3_file* file, int op, void *pArg) {
        return SQLITE_NOTFOUND;
    }
  
    sqlite3_io_methods iomethods = {
        1, //int iVersion;
        dms_xClose,
        dms_xRead,
        dms_xWrite,
        0, // int (*xTruncate)(sqlite3_file*, sqlite3_int64 size);
        dms_xSync,
        dms_xFileSize,
        dms_xLock,
        dms_xUnlock, // int (*xUnlock)(sqlite3_file*, int);
        0, // int (*xCheckReservedLock)(sqlite3_file*, int *pResOut);
        dms_xFileControl, // int (*xFileControl)(sqlite3_file*, int op, void *pArg);
        0, // int (*xSectorSize)(sqlite3_file*);
        dms_xDeviceCharacteristics //int (*xDeviceCharacteristics)(sqlite3_file*);
    };

    /* if buf looks like gzipped data, decompress it, and update *sz.  
     * Return buf, which may now point to new space. */
    char* maybe_decompress(char* buf, sqlite3_int64 *sz) {
        if (*sz<2) return buf;
        /* check for gzip magic number */
        unsigned char s0 = buf[0];
        unsigned char s1 = buf[1];
        if (s0==0x1f && s1==0x8b) {
            std::istringstream file;
            //this ought to work but doesn't, at least on MacOS.
            //file.rdbuf()->pubsetbuf(buf, *sz);
            file.str(std::string(buf, buf+*sz));    // extra copy :(
            std::unique_ptr<istream> in(istream::wrap(file));
            char tmp[16384];
            std::vector<char> contents;
            for (;;) {
                auto rc = in->read(tmp, sizeof(tmp));
                if (rc<=0) break;
                contents.insert(contents.end(), tmp, tmp+rc);
            }
            *sz = contents.size();
            buf = (char *)realloc(buf, contents.size());
            memcpy(buf, contents.data(), contents.size());
        }
        return buf;
    }

    struct dms_vfs : sqlite3_vfs {

        dms_vfs() {
            zName = "dms";
            xOpen = dms_xOpen;
            xDelete = dms_xDelete;
            xAccess = dms_xAccess;
            xFullPathname = dms_xFullPathname;
            xAccess = dms_xAccess;
            mxPathname=1024;
            szOsFile=sizeof(dms_file);
        }

        static int dms_xOpen(sqlite3_vfs* self, const char *zName, 
                sqlite3_file*file, int flags, int *pOutFlags) {
            dms_file *dms = (dms_file *)file;
            dms->pMethods = &iomethods;
            dms->path = NULL;
            if (flags & SQLITE_OPEN_CREATE) {
                dms->contents = NULL;
                dms->size = 0;
                if (flags & SQLITE_OPEN_MAIN_DB) {
                    dms->path = strdup(zName);
                }
            }
            return SQLITE_OK;
        }

        static int dms_xDelete(sqlite3_vfs*, const char *zName, int syncDir) {
            return SQLITE_OK;
        }

        static int dms_xAccess(sqlite3_vfs*, const char *zName, int flags, 
                int *pResOut) {
            switch (flags) {
                case SQLITE_ACCESS_EXISTS:    *pResOut=0; break;
                case SQLITE_ACCESS_READWRITE: *pResOut=0; break;
                case SQLITE_ACCESS_READ:      *pResOut=0; break;
                default:
                                              return SQLITE_ERROR;
            }
            return SQLITE_OK;
        }

        static int dms_xFullPathname(sqlite3_vfs*, const char *zName, int nOut, 
                char *zOut) {
            strncpy(zOut, zName, nOut);
            return 0;
        }
    } vfs[1];
}

const char* Sqlite::errmsg() const {
    return sqlite3_errmsg(_db.get());
}

std::string Sqlite::contents() const {
    dms_file* dms;
    sqlite3_file_control(_db.get(), "main", SQLITE_FCNTL_FILE_POINTER, &dms);
    return std::string(dms->contents, dms->contents+dms->size);
}

void Sqlite::finish() {
    if (_unbuffered) return;
    dms_file* dms;
    sqlite3_file_control(_db.get(), "main", SQLITE_FCNTL_FILE_POINTER, &dms);
    dms->write();
}

Sqlite Sqlite::read(std::string const& path, bool unbuffered) {
    sqlite3* db;
    if (unbuffered) {
        int rc = sqlite3_open_v2( path.data(), &db, SQLITE_OPEN_READONLY, NULL);
        if (rc!=SQLITE_OK) MSYS_FAIL(sqlite3_errmsg(db));
        return Sqlite(std::shared_ptr<sqlite3>(db, sqlite3_close), unbuffered);
    }

    sqlite3_vfs_register(vfs, 0);

    int fd=open(path.c_str(), O_RDONLY);
    if (fd<0) {
        MSYS_FAIL("Reading DMS file at '" << path << "': " << strerror(errno));
    }
    struct stat statbuf[1];
    if (fstat(fd, statbuf)!=0) {
        int _errno = errno;
        ::close(fd);
        MSYS_FAIL("Getting size of DMS file at '" << path << "': " << strerror(_errno));
    }

    ssize_t tmpsize = statbuf->st_size;
    if (tmpsize==0) {
        close(fd);
        MSYS_FAIL("DMS file at '" << path << "' has zero size");
    }
    char* tmpbuf = (char *)malloc(tmpsize);
    if (!tmpbuf) {
        close(fd);
        MSYS_FAIL("Failed to allocate read buffer for DMS file at '" << path
           << "' of size " << tmpsize);
    }
    ssize_t sz = tmpsize;
    char* ptr = tmpbuf;
    while (sz) {
        errno = 0;
        ssize_t rc = ::read(fd, ptr, sz);
        if (rc<0 || (rc==0 && errno!=0)) {
            std::string errmsg = strerror(errno);
            close(fd);
            free(tmpbuf);
            MSYS_FAIL("Error reading DMS contents at " << path 
                    << ": " << errmsg);
        }
        sz -= rc;
        ptr += rc;
    }
    close(fd);
    int rc = sqlite3_open_v2( "::dms::", &db, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, 
            vfs->zName);
    if (rc!=SQLITE_OK) {
        free(tmpbuf);
        MSYS_FAIL(sqlite3_errmsg(db));
    }
    dms_file* dms;
    sqlite3_file_control(db, "main", SQLITE_FCNTL_FILE_POINTER, &dms);
    dms->size = tmpsize;
    dms->contents = maybe_decompress(tmpbuf, &dms->size);
    return std::shared_ptr<sqlite3>(db, sqlite3_close);
}

Sqlite Sqlite::read_bytes(const void* bytes, int64_t len ) {
    sqlite3* db;
    sqlite3_vfs_register(vfs, 0);

    ssize_t tmpsize = len;
    char* tmpbuf = (char *)malloc(len);
    if (!tmpbuf) MSYS_FAIL("Failed to allocate read buffer for DMS file of size " << len);
    memcpy(tmpbuf, bytes, len);
    int rc = sqlite3_open_v2( "::dms::", &db, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, 
            vfs->zName);
    if (rc!=SQLITE_OK) {
        free(tmpbuf);
        MSYS_FAIL(sqlite3_errmsg(db));
    }
    dms_file* dms;
    sqlite3_file_control(db, "main", SQLITE_FCNTL_FILE_POINTER, &dms);
    dms->size = tmpsize;
    dms->contents = maybe_decompress(tmpbuf, &dms->size);
    return std::shared_ptr<sqlite3>(db, sqlite3_close);
}

Sqlite Sqlite::write(std::string const& path, bool unbuffered) {
    sqlite3* db;
    sqlite3_vfs_register(vfs, 0);

    int rc = sqlite3_open_v2(path.c_str(), &db, 
            SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, 
            unbuffered ? NULL : vfs->zName);
    if (rc!=SQLITE_OK) MSYS_FAIL(sqlite3_errmsg(db));

    return Sqlite(std::shared_ptr<sqlite3>(db, sqlite3_close), unbuffered);
}

void Sqlite::exec(std::string const& sql) {
    if (sqlite3_exec(_db.get(), sql.c_str(), NULL, NULL, NULL))
        MSYS_FAIL("Error executing SQL '" << sql << "': " << errmsg());
}

bool Sqlite::has(std::string const& table) const {
    int rc;
    sqlite3_stmt * stmt;
    char * sql = sqlite3_mprintf("select rowid from sqlite_master where name=%Q", table.c_str());
    if (sqlite3_prepare_v2(_db.get(), sql, -1, &stmt, NULL))
        MSYS_FAIL(errmsg());
    sqlite3_free(sql);
    rc=sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    return rc==SQLITE_ROW;
}

int Sqlite::size(std::string const& table) const {
    char * sql;
    sqlite3_stmt * stmt;
    int result;
    sql=sqlite3_mprintf("select count() from %q", table.c_str());
    if (sqlite3_prepare_v2(_db.get(), sql, -1, &stmt, NULL)) 
        MSYS_FAIL(errmsg());
    sqlite3_free(sql);
    sqlite3_step(stmt);
    result=sqlite3_column_int64(stmt,0);
    sqlite3_finalize(stmt);
    return result;
}

Reader Sqlite::fetch(std::string const& table, bool strict) const {
    return Reader(_db, table, strict);
}

const char* Reader::errmsg() const {
    return sqlite3_errmsg(_db.get());
}

static void close_stmt(sqlite3_stmt* stmt) {
    if (sqlite3_finalize(stmt)) 
        MSYS_FAIL("Error finalizing statement: " << sqlite3_errmsg(sqlite3_db_handle(stmt)));
}

Reader::Reader(std::shared_ptr<sqlite3> db, std::string const& table,
               bool strict) 
: _db(db), _table(table), _strict_types(strict) {

    sqlite3_stmt * stmt;

    char* sql = sqlite3_mprintf("select * from %q", table.c_str());
    int rc = sqlite3_prepare_v2(_db.get(), sql, -1, &stmt, NULL);
    sqlite3_free(sql);
    if (rc) return; /* no such table */
    _stmt.reset(stmt, close_stmt);

    /* We used to get the table schema from the first row, but that makes
     * it impossible to distinguish between an empty table and no table
     * at all.  Also, Sqlite lets you store types of different values
     * in the same column, so just reading the types of the first row
     * doesn't necessarily tell you how you should read the others.  The
     * only sensible way to proceed seems to be to parse the table_info
     * pragma and use the "declared" type of the table. */

    /* msys/3.0.1: fall back to the old way when strict_types is false,
     * because otherwise we have no way to handle things like
     * "create table particle (id integer primary key, ... grp_energy)"
     * with no declared type.  */
    if (!_strict_types) {

        _cols.resize(sqlite3_column_count(stmt));
        for (unsigned i=0; i<_cols.size(); i++) {
            _cols[i].first = (const char *)(sqlite3_column_name(stmt,i));
        }

        next();
        if (done()) return;

        for (unsigned i=0; i<_cols.size(); i++) {
            int type = sqlite3_column_type(stmt,i);
            switch(type) {
                default:
                case SQLITE_TEXT:    _cols[i].second = StringType; break;
                case SQLITE_INTEGER: _cols[i].second = IntType; break;
                case SQLITE_FLOAT:   _cols[i].second = FloatType; break;
            }
        }
        return;
    }

    /* get table schema from table_info pragma */
    sql = sqlite3_mprintf("pragma table_info(%q)", table.c_str());
    rc = sqlite3_prepare_v2(_db.get(), sql, -1, &stmt, NULL);
    sqlite3_free(sql);
    if (rc) MSYS_FAIL("Error getting schema for " << table << errmsg());
    /* we expect name and type to be columns 1 and 2, respectively */
    std::shared_ptr<sqlite3_stmt> tmp(stmt, close_stmt);
    while (sqlite3_step(stmt)==SQLITE_ROW) {
        std::string name = (const char *)sqlite3_column_text(stmt,1);
        std::string type = (const char *)sqlite3_column_text(stmt,2);
        ValueType t;
        to_upper(type);
        /* see http://www.sqlite.org/datatype3.html, section 2.1 */
        if (strstr(type.c_str(), "INT")) {
            t = IntType;
        } else if (strstr(type.c_str(), "CHAR") ||
                   strstr(type.c_str(), "CLOB") ||
                   strstr(type.c_str(), "STRING") || /* msys mistake! */
                   strstr(type.c_str(), "TEXT")) {
            t = StringType;
        } else if (strstr(type.c_str(), "BLOB") ||
                   type.size()==0) {
            MSYS_FAIL("Blob type for column " << name << " in table " << table << " not supported.");
        } else if (strstr(type.c_str(), "REAL") ||
                   strstr(type.c_str(), "FLOA") ||
                   strstr(type.c_str(), "DOUB")) {
            t = FloatType;
        } else {
            MSYS_FAIL("Numeric type for column " << name << " in table " << table << " not supported.");
        }
        _cols.push_back(std::make_pair(name,t));
    }
    /* advance to the first row */
    next();
}

void Reader::next() {
    int rc=sqlite3_step(_stmt.get());
    if (rc==SQLITE_ROW) {
        /* nothing to do */
    } else if (rc==SQLITE_DONE) {
        _stmt.reset();
    } else {
        MSYS_FAIL(errmsg());
    }
}

ValueType Reader::current_type(int col) const {
    int type = sqlite3_column_type(_stmt.get(),col);
    switch(type) {
        default:
        case SQLITE_TEXT:    return StringType;
        case SQLITE_INTEGER: return IntType;
        case SQLITE_FLOAT:   return FloatType;
    }
    return StringType;
}

std::string Reader::name(int col) const {
    if (col<0 || col>=size()) MSYS_FAIL("no such column " << col);
    return _cols[col].first;
}

ValueType Reader::type(int col) const {
    if (col<0 || col>=size()) MSYS_FAIL("no such column " << col);
    return _cols[col].second;
}

ValueType Reader::numeric_type(int col) const {
    if (col<0 || col>=size()) MSYS_FAIL("no such column " << col);
    auto type = sqlite3_value_numeric_type(sqlite3_column_value(_stmt.get(), col));
    switch(type) {
        default:
        case SQLITE_TEXT:    return StringType;
        case SQLITE_INTEGER: return IntType;
        case SQLITE_FLOAT:   return FloatType;
    }
    return StringType;
}

int Reader::column(std::string const& name) const {
    for (unsigned i=0; i<_cols.size(); i++) {
        if (name==_cols[i].first) return i;
    }
    return -1;
}

int Reader::get_int(int col) const {
    if (_strict_types && type(col)!=IntType) 
        MSYS_FAIL("Type error reading int from column " << _cols[col].first << " in table " << _table);
    return sqlite3_column_int(_stmt.get(),col);
}

double Reader::get_flt(int col) const {
    if (_strict_types && type(col)!=FloatType) 
        MSYS_FAIL("Type error reading float from column " << _cols[col].first << " in table " << _table);
    return sqlite3_column_double(_stmt.get(),col);
}

const char * Reader::get_str(int col) const {
    if (_strict_types && type(col)!=StringType) 
        MSYS_FAIL("Type error reading string from column " << _cols[col].first << " in table " << _table);
    const char* ret = (const char *)sqlite3_column_text(_stmt.get(),col);
    return ret ? ret : "";
}

Writer Sqlite::insert(std::string const& table) const {
    return Writer(_db, table);
}

Writer::Writer(std::shared_ptr<sqlite3> db, std::string const& table) 
: _db(db) {

    sqlite3_stmt * stmt;
    int i, n=0;
    char * insert_prefix;
    char * sql;

    insert_prefix = sqlite3_mprintf("insert into %q values (", table.c_str());
    sql = sqlite3_mprintf("pragma table_info(%q)", table.c_str());

    std::vector<std::string> cols;
    if (sqlite3_prepare_v2(_db.get(), sql, -1, &stmt, NULL)) {
        sqlite3_free(sql);
        sqlite3_free(insert_prefix);
        MSYS_FAIL(sqlite3_errmsg(_db.get()));
    }
    sqlite3_free(sql);
    while (sqlite3_step(stmt)==SQLITE_ROW) {
        cols.push_back((const char * )sqlite3_column_text(stmt, 1));
        ++n;
    }
    sqlite3_finalize(stmt);

    std::string insert_sql(insert_prefix);
    sqlite3_free(insert_prefix);
    for (i=0; i<n; i++) {
        insert_sql += "?";
        insert_sql += i==n-1 ? ")" : ",";
    }
    if (sqlite3_prepare_v2(_db.get(), insert_sql.c_str(), -1, &stmt, NULL))
        MSYS_FAIL(sqlite3_errmsg(_db.get()));

    _stmt.reset(stmt, close_stmt);
}

void Writer::bind_int(int col, int v) {
    sqlite3_bind_int(_stmt.get(), col+1, v); 
}

void Writer::bind_flt(int col, double v) {
    sqlite3_bind_double(_stmt.get(), col+1, v); 
}

void Writer::bind_str(int col, std::string const& v) {
    sqlite3_bind_text(_stmt.get(), col+1, v.data(), v.size(), SQLITE_TRANSIENT);
}

void Writer::next() {
    if (sqlite3_step(_stmt.get()) != SQLITE_DONE) {
        MSYS_FAIL(sqlite3_errmsg(sqlite3_db_handle(_stmt.get())));
    }
    sqlite3_reset(_stmt.get());
}

