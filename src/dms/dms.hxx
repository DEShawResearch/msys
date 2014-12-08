#ifndef desres_msys_dms_dms_hxx 
#define desres_msys_dms_dms_hxx 

#include "../value.hxx"
#include <boost/shared_ptr.hpp>

struct sqlite3;
struct sqlite3_stmt;

namespace desres { namespace msys {

    class Reader;
    class Writer;

    class Sqlite {
        boost::shared_ptr<sqlite3> _db;

        const char* errmsg() const;

    public:
        Sqlite() {}

        Sqlite(boost::shared_ptr<sqlite3> db)
        : _db(db) {}

        static Sqlite read(std::string const& path);
        static Sqlite read_bytes(const char* bytes, int64_t len);
        static Sqlite write(std::string const& path, bool unbuffered = false);

        void exec(std::string const& sql);

        bool has(std::string const& table) const;
        int size(std::string const& table) const;
        Reader fetch(std::string const& table, bool strict_types=true) const;
        Writer insert(std::string const& table) const;
    };

    class Reader {
        boost::shared_ptr<sqlite3> _db;
        boost::shared_ptr<sqlite3_stmt> _stmt;
        std::string _table;

        typedef std::pair<std::string, ValueType> column_t;
        std::vector<column_t> _cols;
        bool _strict_types;

        const char* errmsg() const;

    public:
        Reader(boost::shared_ptr<sqlite3> db, std::string const& table,
               bool strict_types);
        bool strict_types() const { return _strict_types; }
        void next();
        bool done() const { return !_stmt; }
        operator bool() const { return !done(); }
        int size() const { return _cols.size(); }
        std::string name(int col) const;
        ValueType type(int col) const;

        /* Get the type actually stored for the current row */
        ValueType current_type(int col) const;
        int column(std::string const& name) const;
        int get_int(int col) const;
        double get_flt(int col) const;
        const char* get_str(int col) const;
    };

    class Writer {
        boost::shared_ptr<sqlite3> _db;
        boost::shared_ptr<sqlite3_stmt> _stmt;

    public:
        Writer(boost::shared_ptr<sqlite3> db, std::string const& table);
        void next();
        void bind_int(int col, int v);
        void bind_flt(int col, double v);
        void bind_str(int col, std::string const& v);
    };
}}

#endif
