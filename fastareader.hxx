#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <fstream>

using std::iostream;
using std::endl;
using std::cout;
using std::ifstream;
using std::map;
using std::string;
using std::vector;
using std::logic_error;
using std::out_of_range;
using std::runtime_error;
using std::pair;
using std::vector;
using std::exception;
using std::ifstream;

namespace tkbio {
    typedef unsigned int uint;
    class fastareader {
        static const uint MAGIC_NUMBER;
        static const uint FORMAT_VERSION;
        static unsigned int CACHE_SIZE;
        static char LINE_SEPARATOR;

        string _filename;
        ifstream* _filehandler;
        map<string,vector<pair<size_t,uint> > > _pointers;

        char* _cache;
        string _cache_chromosome;
        int _cache_size;
        int _cache_start;
        int _cache_stop;

        static string get_cache_filename(const char* filename);
        void serialize(const char* filename) const throw (exception);
        static fastareader* deserialize(const char* filename) throw (exception);

        void load_section(const string& chromosome, int start, int stop);
        void initialize();
        void load_cache(const char* filename) throw (exception);

        const fastareader& operator = (const fastareader& rhs);
        fastareader(const fastareader& rhs);
    public:
        fastareader();
        ~fastareader();
        static fastareader* load_genome(const char* filename) throw (exception);
        string get_sequence(const string& chromosome, int start, int end) const throw (out_of_range);
        int get_length(const string& chromosome) const;
        vector<string> sequence_names() const;
    };
}
