#include <string>
#include <stdexcept>
#include <vector>

using std::string;
using std::exception;
using std::vector;

typedef unsigned int uint;

namespace tkbio {
    class motif_counter;
    class moccs_result {
    private:
        moccs_result();
        const moccs_result& operator =(const moccs_result& rhs);
        moccs_result(const moccs_result& rhs);
    private:        
        int _motifsize;
        uint _motif;
        int _observed;
        double _score;
        int _datasize;
        double _expected;
        bool _includes_complementary;
        double _pvalue;
        int* _counts;

        void initialize();
    public:
        moccs_result(int motifsize, unsigned int motif, int distance, int const* counts);
        //moccs_result(int index, unsigned int motif, int observed, double score, int distance, int const* counts);
        ~moccs_result();
        void set_parameters(int observed, double expected, double score, double pvalue, bool complemtnray);
        bool includes_complementary() const {
            return _includes_complementary;
        }
        //int index() const { return _index; }
        string motif(bool reverse=false) const;
        double score() const { return _score; }
        double expected() const { return _expected; }
        int size() const { return _datasize; }
        const int* counts() const { return _counts; }
        

        string to_string() const;
        
        // moccs_result() {
        //     throw std::runtime_error("moccs_result can be instanciated with an array");
        //     // this->index = 0;
        //     // this->motif = 0;
        //     // this->score = 0;
        //     // this->observed = 0;
        //     // this->datasize = 0;
        //     // this->counts = NULL;
        // }

        // const moccs_result& operator =(const moccs_result& rhs) {
        //     if (this == &rhs) return *this;
        //     delete[] counts;
        //     this->index = rhs.index;
        //     this->motif = rhs.motif;
        //     this->score = rhs.score;
        //     this->observed = rhs.observed;
        //     this->datasize = rhs.datasize;
        //     this->counts = new int[rhs.datasize];
        //     memcpy(this->counts, rhs.counts, sizeof(int) * rhs.datasize);
        //     return *this;
        // }
        
        // moccs_result(const moccs_result& rhs) {
        //     this->index = rhs.index;
        //     this->motif = rhs.motif;
        //     this->score = rhs.score;
        //     this->observed = rhs.observed;
        //     this->datasize = rhs.datasize;
        //     this->counts = new int[rhs.datasize];
        //     memcpy(this->counts, rhs.counts, sizeof(int) * rhs.datasize);
        // }
        
        // // moccs_result(int index, unsigned int motif, int observed, double score) {
        // //     this->index = index;
        // //     this->motif = motif;
        // //     this->score = score;
        // //     this->observed = observed;
        // //     this->datasize = 0;
        // //     this->counts = NULL;
        // // }
        // ~moccs_result() {
        //     delete[] counts;
        // }
        // moccs_result(int index, unsigned int motif, int observed, double score, int distance, int const* counts) {
        //     this->index = index;
        //     this->motif = motif;
        //     this->score = score;
        //     this->observed = observed;
        //     this->datasize = distance;
        //     this->counts = new int[distance];
        //     memcpy(this->counts, counts, sizeof(int) * distance);
        // }
        // static bool compare_score(const moccs_result& lhs, const moccs_result& rhs) {
        //     return lhs._score > rhs._score;
        // }
        // static bool compare_score(const moccs_result* lhs, const moccs_result* rhs) {
        //     return lhs._score > rhs._score;
        // }
        static bool compare_score(const moccs_result* lhs, const moccs_result* rhs);
        friend class motif_counter;
    };
    
    class motif_counter {
    public:
        enum MODE {DEFAULT=0, STRAND_CONCIOUS=1, PARTIAL_MOTIFS=2, ENRICHED=4};
    private:
        //string _motif;
        int _motifsize;

        int _reserved;
        int _size;
        int _distance;
        unsigned int* _patterns;
        int** _count;
        bool _bidirectional;
        bool _whole;
        bool _sorted;
        int* _nucleotides; // counts of ACGT
        int* _samples;     // number of tested
        static int default_span;
    private:
        void initialize_buffer(int size);
        void expand_buffer();
        void sort_patterns();

        int find_index(unsigned int code) const;

        motif_counter& operator = (const motif_counter& rhs);
        motif_counter(const motif_counter& rhs);
        motif_counter();
        
    public:
        motif_counter(int motif_size, int distance=0, MODE mode=DEFAULT, int reserved=0);
        //motif_counter(const string& motif, int dist=0);
        ~motif_counter();

        void add_motif(unsigned int code) throw (exception);
        void add_motif(const string& motif) throw (exception);

        void set_sequence(int span, const char* sequence, int center=-1) throw (std::exception);
        static string decode_sequence(int length, int code, bool complementary=false);
        static unsigned int encode_sequence(const string& sequence, int size=0);
        static unsigned int encode_sequence(int length, const char* sequence);
        static unsigned int generate_complementary(int length, unsigned int code);

        vector<moccs_result*> get_results(int max_num, double score_threshold=0.0, double excess=0.0) const;
//        std::string to_string(int max_num=0) const;
    };
}
