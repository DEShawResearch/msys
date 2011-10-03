#ifndef desres_msys_dms_dms_hxx 
#define desres_msys_dms_dms_hxx 

#include "../value.hxx"

struct sqlite3;

typedef struct dms {
    struct sqlite3 * db;
    dms( sqlite3* _db) : db(_db) {}
} dms_t;

/* opaque handles */
typedef struct dms_reader dms_reader_t;
typedef struct dms_writer dms_writer_t;

/* Initialize the dms handle */
dms_t * dms_read( const char * path );

/* Initialize dms handle for writing */
dms_t * dms_write( const char * path );

/* Release dms resources */
void dms_close( dms_t * dms );

/* execute some raw SQL */
void dms_exec( dms_t * dms, const char * sql );

/* does the table exist? */
int dms_has_table( dms_t * dms, const char * name );

/* size of table; 0 if it doesn't exist */
int dms_table_size( dms_t * dms, const char * name );

/* fetch all entries in the table.  NULL if doesn't exist.  Return true
 * if table has at least one entry, false if not. */
int dms_fetch( dms_t * dms, const char * name, dms_reader_t ** r );

/* advance to next row.  If we're done, free r and set to NULL */
void dms_reader_next( dms_reader_t ** r );

/* free the reader */
void dms_reader_free( dms_reader_t * r );

/* get the number of columns */
int dms_reader_column_count( dms_reader_t * r );

/* get the name of the given column */
const char * dms_reader_column_name( dms_reader_t * r, int col );

/* what kind of data was declared for the given column */
desres::msys::ValueType dms_reader_column_type( dms_reader_t * r, int col );

/* get the column for the given name, or -1 if not present */
int dms_reader_column( dms_reader_t * r, const char * name );

/* get an int from the given column */
int dms_reader_get_int( dms_reader_t * r, int col );

/* get a double from the given column */
double dms_reader_get_double( dms_reader_t * r, int col );

/* get a volatile string from the given column */
const char * dms_reader_get_string( dms_reader_t * r, int col );

/* prepare to insert records */
void dms_insert( dms_t * dms, const char * table, dms_writer_t ** w );

/* bind values to columns */
void dms_writer_bind_int( dms_writer_t * w, int col, int v );
void dms_writer_bind_double( dms_writer_t * w, int col, double v );
void dms_writer_bind_string( dms_writer_t * w, int col, const char * v );
void dms_writer_bind_null( dms_writer_t* w, int col);

/* execute the write step */
void dms_writer_next( dms_writer_t * w );

/* free the write */
void dms_writer_free( dms_writer_t * w );

#endif
