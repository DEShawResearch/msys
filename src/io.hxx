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
        SdfFileFormat          = 7

    };

    /* Guess file format for the given path.  Returns UnrecognizedFileFormat
     * if unable to make a guess based on the filename.  */
    FileFormat GuessFileFormat(std::string const& path);

    /* Returns a name for the given file format, for logging purposes.
     * For unrecognized the name will be "UNRECOGNIZED". */
    std::string FileFormatAsString(FileFormat format);

    /* Inverse of previous: file format given a string. */
    FileFormat FileFormatFromString(std::string const& name);

    /* Load using the specified format with default options.  Returns
     * NULL SystemPtr if format is UnrecognizedFileFormat */
    SystemPtr LoadWithFormat(std::string const& path, FileFormat format,
                             bool structure_only = false);

    /* Try to guess the file type from the path, and load the system
     * using default options.  Returns NULL SystemPtr if file type could
     * not be guessed by GuessFileFormat (in this way one can distinguish
     * between errors due to filename and errors due to failure to import
     * the file).  If optional format pointer is supplied, stores 
     * format which was guessed based on the path.  */
    SystemPtr Load(std::string const& path, 
                   bool structure_only = false,
                   FileFormat* opt_format = NULL);


    /* Interface class to iterate over structures */
    class LoadIterator;
    typedef boost::shared_ptr<LoadIterator> LoadIteratorPtr;

    class LoadIterator {
    public:
        virtual ~LoadIterator() {}
        static LoadIteratorPtr create(
                std::string const& path,
                bool structure_only = false,
                FileFormat* opt_format = NULL);

        virtual SystemPtr next() = 0;
    };
}}

#endif
