#include <fastareader.hxx>

using namespace tkbio;
using std::ofstream;
using std::ifstream;
using std::make_pair;

const unsigned int fastareader::MAGIC_NUMBER = 4096;
const unsigned int fastareader::FORMAT_VERSION = 100;
unsigned int fastareader::CACHE_SIZE = 8192;
char fastareader::LINE_SEPARATOR = '/';

fastareader::~fastareader() {
    delete[] _cache;
    if (_filehandler != NULL) {
        _filehandler->close();
        delete _filehandler;
    }
}

string fastareader::get_cache_filename(const char* filename) {
    string fn = filename;
    size_t pos = fn.rfind(LINE_SEPARATOR);
    string cache_file;
    if (pos != string::npos) {
        cache_file = fn.substr(0, pos + 1);
        pos++;
    } else {
        pos = 0;
        cache_file = "";
    }
    cache_file += "." + fn.substr(pos, fn.size() - pos) + ".cache";
    return cache_file;
}

void fastareader::initialize() {
    //_filename = "";
    //_chromosomes;
    //_pointers;
    _cache = NULL;//new char[CACHE_SIZE];
    _cache_size = 0;
    _cache_start = 0;
    _cache_stop = 0;
    _filehandler = NULL;
}

fastareader::fastareader() {
    initialize();
}

/**
   0-4 : magic number
   4-8 : version
   8-12: filename size
   12- : filename
   
   number of chromosomes
   pointers

   8-12: number of sequences
   12- : pointers to chromosomes
   12-16: filename size
   16- : filename
   
   sequence
   0-4:number of elements
   
 */

void fastareader::serialize(const char* filename) const throw (exception) {
    string fn = get_cache_filename(filename);
    ofstream fo(fn.c_str());
    if (fo.is_open() == false) {
        throw runtime_error(string("cannot open file ") + fn);
    }
    fo.write((char*)&MAGIC_NUMBER, sizeof(uint));
    fo.write((char*)&FORMAT_VERSION, sizeof(uint));
    vector<pair<string, size_t> > cachesize;
    for (map<string,vector<pair<size_t, uint> > >::const_iterator it = _pointers.begin(); it != _pointers.end(); it++) {
        size_t size_seq = sizeof(int) + it->first.size();
        size_seq += it->second.size() * (sizeof(size_t) + sizeof(uint)) + sizeof(int) + it->first.size() + sizeof(uint);
        cachesize.push_back(make_pair(it->first, size_seq));
    }
    int fnlen = _filename.size();
    fo.write(reinterpret_cast<const char*>(&fnlen), sizeof(int));
    fo.write(_filename.c_str(), sizeof(char) * fnlen);

    // number of sequencess
    int seqnum = cachesize.size();
    fo.write(reinterpret_cast<const char*>(&seqnum), sizeof(int));
    uint seq_top = 1024;
    uint pointer = seq_top;
    for (int i = 0; i < (int)cachesize.size(); i++) {
        fo.write(reinterpret_cast<const char*>(&pointer), sizeof(uint));
        pointer += cachesize[i].second;
    }
    fo.seekp(seq_top);
    for (int i = 0; i < (int)cachesize.size(); i++) {
        uint nlen = cachesize[i].first.size();
        fo.write(reinterpret_cast<const char*>(&nlen), sizeof(uint));
        fo.write(cachesize[i].first.c_str(), cachesize[i].first.size());
        const string& name = cachesize[i].first;
        map<string,vector<pair<size_t, uint> > >::const_iterator it = _pointers.find(name);
        const vector<pair<size_t, uint> >& ptr = it->second;
        uint clen = ptr.size();//_pointers[].size();
        fo.write(reinterpret_cast<const char*>(&clen), sizeof(uint));
        for (int i = 0; i < (int)ptr.size(); i++) {
            fo.write(reinterpret_cast<const char*>(&ptr[i].first), sizeof(size_t));
            fo.write(reinterpret_cast<const char*>(&ptr[i].second), sizeof(uint));
        }
    }

    fo.close();
}

fastareader* fastareader::deserialize(const char* filename) throw (exception) {
    string fn = get_cache_filename(filename);
    ifstream fi(fn.c_str());
    if (fi.is_open() == false) {
        throw runtime_error(string("cannot open ") + filename);
    }
    try {
        uint magic;
        fi.read(reinterpret_cast<char*>(&magic), sizeof(uint));
        uint version;
        fi.read(reinterpret_cast<char*>(&version), sizeof(uint));
        {
            uint len;
            fi.read(reinterpret_cast<char*>(&len), sizeof(uint));
            char* buf = new char[len + 1];
            fi.read(buf, len);
            buf[len] = '\0';
            string name = buf;
            delete[] buf;
            if (strncmp(filename, name.c_str(), name.size()) != 0) {
                throw runtime_error(string("cache is incompatible with ") + filename);
            }
        }
        int seqnum;
        fi.read(reinterpret_cast<char*>(&seqnum), sizeof(int));
        vector<size_t> pointers;
        for (int i = 0; i < seqnum; i++) {
            size_t pos;
            fi.read(reinterpret_cast<char*>(&pos), sizeof(size_t));
            pointers.push_back(pos);
        }
        fastareader* reader = new fastareader();
        reader->_filename = filename;
        for (int i = 0; i < (int)pointers.size(); i++) {
            int len;
            char* buf;
            fi.seekg(pointers[i]);
            fi.read(reinterpret_cast<char*>(&len), sizeof(uint));
            buf = new char[len + 1];
            fi.read(buf, len);
            buf[len] = '\0';
            string name = buf;
            fi.read(reinterpret_cast<char*>(&len), sizeof(uint));
            vector<pair<size_t, uint> > pointers;
            for (int j = 0; j < len; j++) {
                size_t fpos;
                uint spos;
                fi.read(reinterpret_cast<char*>(&fpos), sizeof(size_t));
                fi.read(reinterpret_cast<char*>(&spos), sizeof(uint));
                pointers.push_back(make_pair(fpos, spos));
            }
            reader->_pointers[name] = pointers;
        }
        
        fi.close();
        //reader->_pointers[name] = points;
        //reader->_pointers = pointers;
        return reader;
    } catch (exception& e) {
        fi.close();
        throw;
    }
}

fastareader* fastareader::load_genome(const char* filename) throw (exception) {
    // try {
    //     //fastareader* reader = NULL;
    //     fastareader* reader = deserialize(filename);
    //     return reader;
    // } catch (exception& e) {
    // }
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw runtime_error(string("cannot open ") + filename);
    }
    //string chromosome;
    fastareader* fr = new fastareader();
    vector<pair<size_t,uint> > positions;
    uint position = 0;
    string name;
    uint next_position = 0;
    uint skip = CACHE_SIZE / 2;
    size_t fpos = 0;
    while (fi.is_open() == false) {
        string line;
        //size_t 
        fpos = fi.tellg();
        getline(fi, line);
        if (line.c_str()[0] == '>') {
            if (positions.size() > 0 && name != "") {
                positions.push_back(std::make_pair(fpos, position));
                fr->_pointers[name] = positions;
            }
            positions.erase(positions.begin(), positions.end());
            position = 0;
            next_position = skip;
            name = line.substr(1, line.size() - 1);
            std::cerr << "NAME : " << name << " at " << fpos << endl;
            positions.push_back(make_pair(fi.tellg(), 0));
        } else {
            if (position >= next_position) {
                positions.push_back(make_pair(fpos, position));
                next_position += skip;
            }
            for (int i = 0; i < (int)line.size(); i++) {
                char c = line.c_str()[i];
                if ((c >= 'A' && c <= 'Z') && (c >= 'a' && c <= 'z')) {
                    position ++;
                }
            }
        }
    }
    if (positions.size() > 0 && name != "") {
        positions.push_back(std::make_pair(fpos, position));
        fr->_pointers[name] = positions;
    }
    fi.close();
    try {
        fr->serialize(filename);
    } catch (exception& e) {
        std::cerr << e.what() << " : cannot save cache file" << endl;
    }
    return fr;
}

string fastareader::get_sequence(const string& chromosome, int start, int stop) const throw (out_of_range) {
    const_cast<fastareader*>(this)->load_section(chromosome, start, stop);
    int offset = start - _cache_start;
    int length = stop - start;
    if (offset < 0) { // out of range
        length += offset;
        offset = 0;
    } 
    if (length < 0) length = 0;
    if (length > _cache_stop - _cache_start) {
        length = _cache_stop - _cache_start;
    }
    return string(_cache + offset, length);
}

void fastareader::load_section(const string& chromosome, int start, int stop) {
// already cached
    if (_cache_chromosome == chromosome && _cache_start <= start && stop <= _cache_stop) {
        return;
    }
    if (_filehandler != NULL) {
        _filehandler = new ifstream(_filename.c_str());
        if (_filehandler->is_open() == false) {
            throw runtime_error(string("cannot open fasta file ") + _filename);
        }
    }
    _cache_chromosome = chromosome;

    if (start < 0) start = 0;
    const map<string,vector<pair<size_t, uint> > >::const_iterator it = _pointers.find(chromosome);
    if (it == _pointers.end()) { // no chromosome
        _cache_start = _cache_stop = 0;
        return;
        //return string("");
    }
    const vector<pair<size_t, uint> >& pointers = it->second;//_pointers[chromosome];
    if (pointers.size() == 0) { // no pointers
        _cache_start = _cache_stop = 0;
        return;
        //return string("");
    }
    // set in range
    if (start < 0) start = 0;
    if (stop >= pointers[_pointers.size() - 1].second) {
        stop = pointers[_pointers.size() - 1].second;
    }
    
    int left = 0;
    int right = pointers.size();
    int index = -1;
    for (;;) {
        int center = (left + right) / 2;
        const pair<size_t, uint>& pos = pointers[center];
        if (pos.first > stop) {
            right = center;
        } else if (pos.first < start - CACHE_SIZE) {
            left = center + 1;
        } else {
            index = center;
            break;
        }
    }
    if (index < 0) {
        throw out_of_range("");
    }

    int index_begin = index;
    int index_end = index;
    while (index >= 0) {
        const pair<size_t, uint>& pos = pointers[index];
        if (pos.first > start) {
            index_begin = index;
            break;
        }
        index--;
    }
    index = index + 1;
    while (index < (int)pointers.size()) {
        const pair<size_t, uint>& pos = pointers[index];
        if (pos.first > stop) {
            index_end = index;
            break;
        }
        index++;
    }
    int required_buffer_size = pointers[index_end].second - pointers[index_begin].second;
    if (required_buffer_size > _cache_size) {
        delete[] _cache;
        _cache = new char[required_buffer_size];
        _cache_size = required_buffer_size;
    }
    _filehandler->seekg(pointers[index_begin].second);
    _filehandler->read(_cache, required_buffer_size);
    int ptr = 0;
    for (int i = 0; i < required_buffer_size; i++) {
        char c = _cache[i];
        if (c >= 'a' && c <= 'z') {
            c -= (char)('a' - 'A');
        } else if (c < 'A' or c > 'Z') {
            continue;
        }
        _cache[ptr++] = c;
    }
    _cache_chromosome = chromosome;
    _cache_start = pointers[index_begin].first;
    _cache_stop = pointers[index_end].second;
}

