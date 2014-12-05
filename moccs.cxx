#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iomanip>

using std::vector;
using std::pair;
using std::make_pair;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;

#include <moccs.hxx>

using namespace tkbio;

motif_counter::motif_counter(int motif_size, int distance, MODE mode, int buffer_size) {
    _motifsize = motif_size;
    _distance = distance;
    _bidirectional = (mode & motif_counter::STRAND_CONCIOUS) == 0;
    _whole = (mode & motif_counter::PARTIAL_MOTIFS) == 0;
    //_reserved = 0;
    //_size = 0;
    _patterns = NULL;
    _count = NULL;
    _nucleotides = new int[distance];
    for (int i = 0; i < _distance; i++) {
        _nucleotides[i] = 0;
    }
    if (buffer_size == 0) {
        initialize_buffer(1 << (motif_size * 2));
        _size = _reserved;
        for (int i = 0; i < _size; i++) {
            _patterns[i] = i;
        }
        _whole = true;
    } else {
        initialize_buffer(128);
        _whole = false;
        _sorted = false;
    }
}

void motif_counter::initialize_buffer(int size) {
    _reserved = size;
    _patterns = new unsigned int[_reserved];
    _count = new int*[_reserved];
    _samples = new int[_reserved];
    for (int i = 0 ; i < _reserved; i++) {
        _count[i] = new int[_distance];
        _samples[i] = 0;
        for (int j = 0; j < _distance; j++) _count[i][j] = 0;
        //cout << "buffer " << i << " : " << _distance << endl;
    }
    _nucleotides = new int[4];
    for (int i = 0; i < 4; i++) _nucleotides[i] = 0;
    _size = 0;
}

motif_counter::~motif_counter() {
    for (int i = 0; i < _reserved; i++) {
        delete[] _count[i];
    }
    delete[] _count;
    delete[] _patterns;
    delete[] _nucleotides;
    delete[] _samples;
}

// motif_counter& motif_counter::operator = (const motif_counter& rhs) {
//     if (this == &rhs) {
//         return *this;
//     }
//     return *this;
// }

// motif_counter(const motif_counter& rhs);
// motif_counter();

void motif_counter::expand_buffer() {
    if (_reserved <= 0) {
        initialize_buffer(128);
    } else {
        int r2 = _reserved * 2;
        unsigned int* pat = new unsigned int[r2];
        memcpy(pat, _patterns, sizeof(unsigned int) * _size);
        delete[] _patterns;
        _patterns = pat;
        
        int** cnt = new int*[r2];
        memcpy(cnt, _count, sizeof(int*) * _size);
        for (int i = _size; i < r2; i++) {
            cnt[i] = new int[_distance];
            for (int j = 0; j < _distance; j++) {
                cnt[i][j] = 0;
            }
        }
        delete[] _count;
        _count = cnt;
        _reserved = _reserved * 2;
    }
}

void motif_counter::add_motif(unsigned int code) throw (std::exception) {
    if (_whole) {
        throw std::runtime_error("cannot add motif for completed set");
    }
    while (_reserved == _size) {
        expand_buffer();
    }
    _patterns[_size] = code;
    for (int i = 0; i < _distance; i++) {
        _count[_size][i] = 0;
    }
    _size++;
    _sorted = false;
}

std::string motif_counter::decode_sequence(int size, int code, bool complementary) {
    string pat;
    if (complementary) {
        for (int i = 0; i < size; i++) {
            int num = (code >> (i * 2)) & 0x03;
            switch (num) {
            case 0:
                pat += 'T'; break;
            case 1:
                pat += 'G'; break;
            case 2:
                pat += 'C'; break;
            case 3:
                pat += 'A'; break;
            default:
                pat += '.'; break;
                break;
                // useless
            }
        }
    } else {
        for (int i = 0; i < size; i++) {
            int num = (code >> ((size - 1 - i) * 2)) & 0x03;
            switch (num) {
            case 0:
                pat += 'A'; break;
            case 1:
                pat += 'C'; break;
            case 2:
                pat += 'G'; break;
            case 3:
                pat += 'T'; break;
            default:
                pat += '.'; break;
                break;
                // useless
            }
        }
    }
    return pat;
}

namespace {
    unsigned int INVALID_MOTIF = 0x7fffffff;
}

unsigned int motif_counter::encode_sequence(int length, const char* sequence) {
    unsigned int code = 0;
    for (int i = 0; i < length; i++) {
        char c = sequence[i];
        code <<= 2;
        switch (c) {
        case 'A': break;
        case 'C': code |= 1; break;
        case 'G': code |= 2; break;
        case 'T': code |= 3; break;
        default:
            return INVALID_MOTIF;//(unsigned int)0xffffffff;
        }
    }
    return code;
}

unsigned int motif_counter::encode_sequence(const string& sequence, int size) {
    if (size <= 0) size = sequence.size();
    unsigned int code = 0;
    for (int i = 0; i < size; i++) {
        char c = sequence.c_str()[i];
        code <<= 2;
        switch (c) {
        case 'A': break;
        case 'C': code |= 1; break;
        case 'G': code |= 2; break;
        case 'T': code |= 3; break;
        default:
            return (unsigned int)0xffffff;
        }
    }
    return code;
}

void motif_counter::sort_patterns() {
    if (_sorted) return;
    std::sort(_patterns, _patterns + _size);
    for (int i = 0; i < _size; i++) {
        for (int j = 0; j < _distance; j++) {
            _count[i][j] = 0;
        }
    }
    _sorted = true;
}

unsigned int motif_counter::generate_complementary(int size, unsigned int code) {
    unsigned int cmp = 0;
    for (int i = 0; i < size; i++) {
        cmp <<= 2;
        cmp |= (((code >> (i * 2)) & 3) ^ 3);
    }
    return cmp;
}

void motif_counter::set_sequence(int size, const char* sequence, int center) throw (exception) {
    sort_patterns();
    if (center < 0) {
        center = size / 2;
    }
    //cout << size << ", " << center << " "<<  strlen(sequence) << std::endl;
    int half = _motifsize / 2;
    int loops = _bidirectional ? 2 : 1;
    for (int i = 0; i < size; i++) {
        char c = sequence[i];
        int index;
        switch (c) {
        case 'A': index = 0; break;
        case 'C': index = 1; break;
        case 'G': index = 2; break;
        case 'T': index = 3; break;
        default:
            continue;
        }
        _nucleotides[index]++;
    }
    for (int i = 0; i < size - _motifsize; i++) {
        //cout << string(sequence + i, _motifsize) << endl;
        unsigned int code = encode_sequence(_motifsize, sequence + i);
        if (code == INVALID_MOTIF) {
            i += _motifsize;
            continue;
        }
        for (int loop = 0; loop < loops; loop++) {
            int pos = i < center - half ? center - i - half: i - center + half;
            //cout << "\n";
            //cout << i << " : " << center << " " << pos << " " << std::hex << code << std::dec << " / " << size << ", " << _motifsize << " / " << size << endl;//", " << sequence << endl;
            
            if (pos < 0 || pos >= _distance) continue;
            _samples[pos] ++;
            //cout << _whole << " " << std::hex << code << std::dec << " " << pos << " / " << _distance << endl;
            if (_whole) {
                _count[code][pos] ++;
            } else {
                int left = 0;
                int right = _size;
                int pivot;
                for (;;) {
                    pivot = (left + right) / 2;
                    unsigned int pat = _patterns[pivot];
                    if (pat < code) {
                        left = pivot + 1;
                    } else if (pat > code) {
                        right = pivot;
                    } else {
                        _count[pivot][pos]++;
                        break;
                    }
                    if (left == right) break;
                }
            }
            if (_bidirectional) {
                code = generate_complementary(_motifsize, code);
            }
        }
    }
}

namespace {
    template<typename T1,typename T2> bool sort_by_second(const std::pair<T1,T2>& lhs, const std::pair<T1,T2>& rhs) {
        return lhs.second > rhs.second;
    }
}
    

int motif_counter::find_index(unsigned int code) const {
    const_cast<motif_counter*>(this)->sort_patterns();
    int left = 0;
    int right = _size;
    for (;;) {
        if (left == right) {
            return -1;
        }
        int center = (left + right) / 2;
        unsigned int pat = _patterns[center];
        if (pat < code) {
            left = center + 1;
        } else if (pat > code) {
            right = center;
        } else {
            return center;
        }
    }
    return -1;
}

namespace {
    double calculate_expected(int size, uint motif, int N, double bases[4]) {
        double ratio = 1.0;
        for (int i = 0; i < size; i++) {
            uint b = motif & 3;
            ratio *= bases[b];
            motif >>= 2;
        }
        return N * ratio;
    }
}

//std::string motif_counter::to_string(int max_num) const {
vector<moccs_result*> motif_counter::get_results(int max_num, double score_threshold, double oe_threshold) const {
    if (max_num <= 0 || max_num > _size) max_num = _size;
    std::vector<moccs_result*> scores;
    //std::vector<std::pair<int,double> > scores;
    int* bufint = NULL;
    if (_bidirectional) {
        bufint = new int[_distance];
    }
    double composition[4];
    int total = _nucleotides[0] + _nucleotides[1] + _nucleotides[2] + _nucleotides[3];
    if (total <= 0) {
        return scores;
    }
    for (int i = 0; i < 4; i++) {
        composition[i] = (double)_nucleotides[i] / total;
    }

    double* normalized = new double[_distance];
    for (int i = 0; i < _size; i++) {
        unsigned int code = _patterns[i];
        const int* row = _count[i];
        if (_bidirectional) {
            unsigned int rev = generate_complementary(_motifsize, code);
            if (code < rev) {
                int ic = find_index(rev);
                if (ic >= 0) {
                    for (int j = 0; j < _distance; j++) {
                        bufint[j] = _count[i][j] + _count[ic][j];
                    }
                    row = bufint;
                }
            } else if (code > rev) {
                continue;
            }
        }
        //cerr << i << " " << code << endl;

        // evaluation
        int sum = 0;
        double score = 0.0;
        double accum = 0.0;
        double area = 0.0;
        double expected = calculate_expected(_motifsize, code, total, composition);
        double pvalue = 1.0;
        
        for (int j = 0; j < _distance; j++) {
            sum += row[j];
        }
        if (sum <= 0) {
            continue;
        }
        if (sum < oe_threshold * expected) {
            continue;
        }
        for (int j = 0; j < _distance; j++) {
            normalized[j] = (double)row[j] / sum;
        }
        for (int j = 0; j < _distance; j++) {
            accum += normalized[j];
            area += accum;
        }
        score = area - _distance * 0.5;

        if (score >= score_threshold) {
            moccs_result* item = new moccs_result(_motifsize, _patterns[i], _distance, row);
            item->set_parameters(sum, expected, score, pvalue, _bidirectional);
            scores.push_back(item);
            //sum, score, _distance, row));
        }
    }
    //cerr << "exit loop\n";
    delete[] bufint;
    delete[] normalized;
    sort(scores.begin(), scores.end(), moccs_result::compare_score);
    if (max_num > 0 && (int)scores.size() > max_num) {
        for (int i = max_num; i < (int)scores.size(); i++) {
            delete scores[i];
        }
        scores.erase(scores.begin() + max_num, scores.end());
    }
    return scores;
}


moccs_result::moccs_result() {
    throw std::runtime_error("moccs_result can be instanciated with an array");
}

// const moccs_result& moccs_result::operator =(const moccs_result& rhs) {
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

moccs_result::moccs_result(const moccs_result& rhs) {
    _motifsize = rhs._motifsize;
    _motif = rhs._motif;
    _score = rhs._score;
    _observed = rhs._observed;
    _expected = rhs._expected;
    _datasize = rhs._datasize;
    _pvalue = rhs._pvalue;
    _includes_complementary = rhs._includes_complementary;
    _counts = new int[rhs._datasize];
    memcpy(_counts, rhs._counts, sizeof(int) * rhs._datasize);
}

        // moccs_result(int index, unsigned int motif, int observed, double score) {
        //     this->index = index;
        //     this->motif = motif;
        //     this->score = score;
        //     this->observed = observed;
        //     this->datasize = 0;
        //     this->counts = NULL;
        // }
moccs_result::~moccs_result() {
    delete[] _counts;
}

//moccs_result::moccs_result(int index, unsigned int motif, int observed, double score, int distance, int const* counts) {
moccs_result::moccs_result(int motifsize, unsigned int motif, int distance, int const* counts) {
    _motifsize = motifsize;
    _motif = motif;
    _datasize = distance;
    _counts = new int[distance];
    _score = 0.0;
    _observed = 0;
    _includes_complementary = true;
    _expected = 0.0;
    _pvalue = 1.0;
    memcpy(_counts, counts, sizeof(int) * distance);
}

void moccs_result::set_parameters(int observed, double expected, double score, double pvalue, bool complementary) {
    _observed = observed;
    _expected = expected;
    _score = score;
    _pvalue = pvalue;
    _includes_complementary = true;
}

bool moccs_result::compare_score(const moccs_result* lhs, const moccs_result* rhs) {
    return lhs->_score > rhs->_score;
}

std::string moccs_result::to_string() const {
    string contents;
    std::stringstream ss;
    ss << motif_counter::decode_sequence(_motifsize, _motif);
    if (_includes_complementary) ss << "/" << motif_counter::decode_sequence(_motifsize, _motif, true);
    ss << "\t" << _observed << "/" << std::setprecision(3) << _expected;
    ss << "\t" << _score;
    for (int i = 0; i < _datasize; i++) {
        ss << (i == 0 ? "\t" : ",") << _counts[i];
    }
    return ss.str();
}
