#include "dms.hxx"
#include <sqlite3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

namespace DM = desres::msys;

#define THROW_FAILURE( args ) do { \
    std::stringstream ss; \
    ss << args; \
    throw std::runtime_error(ss.str()); \
} while (0)

struct dms_reader {
    sqlite3_stmt * stmt;
    int ncols;
    char ** cols;
    DM::ValueType * types;
};

struct dms_writer {
    sqlite3_stmt * stmt;
    std::vector<char *> cols;
};

dms_t * dms_read( const char * path ) {
    sqlite3* db;
    if (sqlite3_open_v2(path, &db, SQLITE_OPEN_READONLY, NULL)!=SQLITE_OK)
        THROW_FAILURE(sqlite3_errmsg(db));
    return new dms(db);
}

dms_t * dms_write( const char * path ) {
    sqlite3* db;
    if (sqlite3_open(path, &db)!=SQLITE_OK)
        THROW_FAILURE(sqlite3_errmsg(db));
    return new dms(db);
}

void dms_close( dms_t * dms ) {
    sqlite3_close(dms->db);
    delete dms;
}

void dms_exec( dms_t * dms, const char * sql ) {
    if (sqlite3_exec(dms->db, sql, NULL, NULL, NULL))
        THROW_FAILURE("Error executing SQL '" << sql << "': " 
                << sqlite3_errmsg(dms->db));
}

int dms_has_table( dms_t * dms, const char * name ) {
    int rc;
    sqlite3_stmt * stmt;
    char * sql = sqlite3_mprintf("select rowid from sqlite_master where name=%Q", name);
    if (sqlite3_prepare_v2(dms->db, sql, -1, &stmt, NULL))
        THROW_FAILURE(sqlite3_errmsg(dms->db));
    sqlite3_free(sql);
    rc=sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    return rc==SQLITE_ROW;
}

int dms_table_size( dms_t * dms, const char * name ) {
    char * sql;
    sqlite3_stmt * stmt;
    int result;
    sql=sqlite3_mprintf("select count() from %q", name);
    if (sqlite3_prepare_v2(dms->db, sql, -1, &stmt, NULL))
        THROW_FAILURE(sqlite3_errmsg(dms->db));
    sqlite3_free(sql);
    sqlite3_step(stmt);
    result=sqlite3_column_int64(stmt,0);
    sqlite3_finalize(stmt);
    return result;
}

int dms_fetch( dms_t * dms, const char * name, dms_reader_t ** pr ) {
    sqlite3_stmt * stmt;
    char * sql;
    dms_reader_t * r;
    int i, rc;

    *pr = NULL;
    sql = sqlite3_mprintf("select * from %q", name);
    rc = sqlite3_prepare_v2(dms->db, sql, -1, &stmt, NULL);
    sqlite3_free(sql);
    if (rc) return 0;

    r=new dms_reader_t;
    r->stmt=stmt;
    r->ncols=sqlite3_column_count(stmt);
    r->cols = new char *[r->ncols];
    r->types = new DM::ValueType[r->ncols];
    for (i=0; i<r->ncols; i++) {
        r->cols[i] = strdup(sqlite3_column_name(stmt,i));
    }
    dms_reader_next(&r);
    if (!r) return 0;

    *pr = r;
    for (i=0; i<r->ncols; i++) {
        int type = sqlite3_column_type(stmt,i);
        switch(type) {
            default:
            case SQLITE_TEXT:    r->types[i]=DM::StringType; break;
            case SQLITE_INTEGER: r->types[i]=DM::IntType; break;
            case SQLITE_FLOAT:   r->types[i]=DM::FloatType; break;
        }
    }
    return 1;
}

void dms_reader_next( dms_reader_t ** pr ) {
    dms_reader_t * r = *pr;
    int rc;
    if (!r) return;
    rc=sqlite3_step(r->stmt);
    if (rc==SQLITE_ROW) {
        /* nothing to do */
    } else if (rc==SQLITE_DONE) {
        /* end of table.  free and set *pr to NULL */
        dms_reader_free(r);
        *pr = NULL;
    } else {
        THROW_FAILURE(sqlite3_errmsg(sqlite3_db_handle(r->stmt)));
    }
}
void dms_reader_free( dms_reader_t * r ) {
    if (!r) return;
    sqlite3_finalize(r->stmt);
    while (r->ncols) free(r->cols[--r->ncols]);
    delete [] r->cols;
    delete [] r->types;
    delete r;
}

int dms_reader_column_count( dms_reader_t * r ) {
    return r->ncols;
}

/* get the name of the given column */
const char * dms_reader_column_name( dms_reader_t * r, int col ) {
    if (col<0 || col>=r->ncols)
        THROW_FAILURE("no such column " << col);
    return r->cols[col];
}

DM::ValueType dms_reader_column_type( dms_reader_t * r, int col ) {
    if (col<0 || col>=r->ncols)
        THROW_FAILURE("no such column " << col);
    return r->types[col];
}

int dms_reader_column( dms_reader_t * r, const char * col ) {
    int i,n=r->ncols;
    for (i=0; i<n; i++) {
        if (!strcmp(r->cols[i],col)) return i;
    }
    return -1;
}

int dms_reader_get_int( dms_reader_t * r, int col ) {
    if (col<0 || col>=r->ncols)
        THROW_FAILURE("no such column " << col);
    return sqlite3_column_int(r->stmt,col);
}

double dms_reader_get_double( dms_reader_t * r, int col ) {
    if (col<0 || col>=r->ncols)
        THROW_FAILURE("no such column " << col);
    return sqlite3_column_double(r->stmt,col);
}

const char * dms_reader_get_string( dms_reader_t * r, int col ) {
    if (col<0 || col>=r->ncols)
        THROW_FAILURE("no such column " << col);
    return (const char *)sqlite3_column_text(r->stmt,col);
}

void dms_insert( dms_t * dms, const char * table, dms_writer_t ** pw ) {
    dms_writer_t * w;
    sqlite3_stmt * stmt;
    int i, n=0;
    char * insert_prefix;
    char * sql;

    insert_prefix = sqlite3_mprintf("insert into %q values (", table);
    sql = sqlite3_mprintf("pragma table_info(%q)", table);

    std::vector<char*> cols;
    if (sqlite3_prepare_v2(dms->db, sql, -1, &stmt, NULL))
        THROW_FAILURE(sqlite3_errmsg(dms->db));
    sqlite3_free(sql);
    while (sqlite3_step(stmt)==SQLITE_ROW) {
        cols.push_back(strdup((const char * )sqlite3_column_text(stmt, 1)));
        ++n;
    }
    sqlite3_finalize(stmt);

    std::string insert_sql(insert_prefix);
    for (i=0; i<n; i++) {
        insert_sql += "?";
        insert_sql += i==n-1 ? ")" : ",";
    }
    if (sqlite3_prepare_v2(dms->db, insert_sql.c_str(), -1, &stmt, NULL))
        THROW_FAILURE(sqlite3_errmsg(dms->db));

    w = new dms_writer_t;
    w->stmt=stmt;
    w->cols=cols;
    *pw = w;
}

void dms_writer_bind_null( dms_writer_t * w, int col) {
    sqlite3_bind_null(w->stmt, col+1);
}
void dms_writer_bind_int( dms_writer_t * w, int col, int v ) {
    sqlite3_bind_int(w->stmt, col+1, v); 
}

void dms_writer_bind_double( dms_writer_t * w, int col, double v ) {
    sqlite3_bind_double(w->stmt, col+1, v); 
}

void dms_writer_bind_string( dms_writer_t * w, int col, const char* v ) {
    sqlite3_bind_text(w->stmt, col+1, v, -1, SQLITE_TRANSIENT);
}

void dms_writer_next( dms_writer_t * w ) {
    if (sqlite3_step(w->stmt) != SQLITE_DONE) {
        THROW_FAILURE(sqlite3_errmsg(sqlite3_db_handle(w->stmt)));
    }
    sqlite3_reset(w->stmt);
}

void dms_writer_free( dms_writer_t * w ) {
    if (!w) return;
    sqlite3_finalize(w->stmt);
    delete w;
}

