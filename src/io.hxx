#ifndef msys_io_hxx
#define msys_io_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Supported file formats */
    enum FileFormat {
        UnrecognizedFileFormat = 0,
        DmsFileFormat          = 1,
        MaeFileFormat          = 2,
        PdbFileFormat          = 3,
        ParmTopFileFormat      = 4,
        Mol2FileFormat         = 5,
        XyzFileFormat          = 6,
        SdfFileFormat          = 7,
        WebPdbFileFormat       = 8,
        PsfFileFormat          = 9

    };

    /* Guess file format for the given path.  Returns UnrecognizedFileFormat
     * if unable to make a guess based on the filename.  */
    FileFormat GuessFileFormat(std::string const& path);

    /* Returns a name for the given file format, for logging purposes.
     * For unrecognized the name will be "UNRECOGNIZED". */
    std::string FileFormatAsString(FileFormat format);

    /* Inverse of previous: file format given a string. */
    FileFormat FileFormatFromString(std::string const& name);

    SystemPtr LoadWithFormat(std::string const& path, FileFormat format,
                             bool structure_only,
                             bool without_tables);

    /* Load using the specified format with default options.  Returns
     * NULL SystemPtr if format is UnrecognizedFileFormat */
    inline
    SystemPtr LoadWithFormat(std::string const& path, FileFormat format,
                             bool structure_only = false) {
        return LoadWithFormat(path, format, structure_only, structure_only);
    }

    SystemPtr Load(std::string const& path,
                   FileFormat* opt_format,
                   bool structure_only,
                   bool without_tables);

    /* Try to guess the file type from the path, and load the system
     * using default options.  Returns NULL SystemPtr if file type could
     * not be guessed by GuessFileFormat (in this way one can distinguish
     * between errors due to filename and errors due to failure to import
     * the file).  If optional format pointer is supplied, stores 
     * format which was guessed based on the path.  */
    inline
    SystemPtr Load(std::string const& path, 
                   bool structure_only = false,
                   FileFormat* opt_format = NULL) {

        return Load(path, opt_format, structure_only, structure_only);
    }


    /* Interface class to iterate over structures */
    class LoadIterator;
    typedef std::shared_ptr<LoadIterator> LoadIteratorPtr;

    class LoadIterator {
    public:
        virtual ~LoadIterator() {}

        static LoadIteratorPtr create(
                std::string const& path,
                FileFormat* opt_format,
                bool structure_only,
                bool without_tables);

        static inline LoadIteratorPtr create(
                std::string const& path,
                bool structure_only = false,
                FileFormat* opt_format = NULL) {
            return create(path, opt_format, structure_only, structure_only);
        }

        virtual SystemPtr next() = 0;
    };

    // IndexedFileLoader provides random access to multi-structure files.
    class IndexedFileLoader {
    public:
        virtual ~IndexedFileLoader() {}
        virtual std::string const& path() const = 0;
        virtual size_t size() const = 0;
        virtual SystemPtr at(size_t zero_based_entry) const = 0;

        // open an indexed file loader, inferring the type based on
        // file name.  The index file must already exist and is expected
        // to be placed at $path.idx; you may optionally specify your own 
        // index file path.
        static std::shared_ptr<IndexedFileLoader> open(
                std::string const& path,
                std::string const& idx_path = "");

        // create an index for the given file, inferring type based on
        // file name.  Place index file at idx_path, defaulting to 
        // $path.idx.
        static void index(std::string const& path,
                          std::string const& idx_path = "");

        // convenience: open an indexed file loader, creating the 
        // index file if necessary at $path.idx
        static std::shared_ptr<IndexedFileLoader> create(
                std::string const& path);

    };

    /* Flags for Save. */
    struct SaveOptions {
        enum Flags {
            Default         = 0,
            Append          = 1 << 0,
            StructureOnly   = 1 << 1
        };
    };

    /* Save the system with the specified format. */
    void SaveWithFormat(SystemPtr mol, 
                        std::string const& path, 
                        Provenance const& prov,
                        FileFormat format,
                        unsigned flags);

    /* Save using format guessed from GuessFileFormat.  */
    void Save(SystemPtr mol, 
              std::string const& path, 
              Provenance const& prov,
              unsigned flags);


}}

#endif
