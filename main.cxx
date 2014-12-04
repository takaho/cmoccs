#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdexcept>

using std::istream;
using std::ostream;
using std::iostream;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::cout;
using std::cin;
using std::ifstream;
using std::map;
using std::string;
using std::vector;
using std::logic_error;
using std::out_of_range;
using std::runtime_error;
using std::cerr;

#include <moccs.hxx>
#include <tktools.hxx>
#include <fastareader.hxx>

using namespace tkbio;
using namespace tktools;
using namespace tktools::util;
using namespace tktools::io;

typedef unsigned long long ullong;

namespace {
    void show_help() {
        cerr << "OPTIONS\n";
        cerr << " -i <filename>   : filename\n";
        cerr << " -g <filename>   : genome\n";
        cerr << " -o <filename>   : output\n";
        cerr << " -w <number>     : window size\n";
        cerr << " -m <number>     : motif size\n";
        cerr << " -n <number>     : number of results (default:0=all)\n";
        cerr << " -e <float>      : excess per expected (default:1)\n";
        cerr << " -s <float>      : score threshold (default:0)\n";
        cerr << "--strand        : strand specific search\n";
    }
}

int main(int argc, char** argv) {
    // if (1) {
    //     fastareader* r = fastareader::load_genome("/Data/Genomes/M_musculus/mm9/mm9.mfa");
    //     cout << r->get_sequence("chr8", 110011000, 110012000) << endl;
    //     delete r;
    //     exit(0);
    // }
    
    try {
        enum file_mode {BED=1, FASTA=2};
        if (has_option(argc, argv, "h")) {
            show_help();
            return 0;
        }
        const char* filename_input = get_argument_string(argc, argv, "i", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* filename_genome = get_argument_string(argc, argv, "g", NULL);
        double excess = get_argument_float(argc, argv, "e", 1.0);
        double score_threshold = get_argument_float(argc, argv, "s", 0.0);
        bool verbose = has_option(argc, argv, "verbose");
        int window_size = get_argument_integer(argc, argv, "w", 250);
        int motif_size = get_argument_integer(argc, argv, "m", 6);
        bool ignore_complementary = has_option(argc, argv, "-strand");
        int num_display = get_argument_integer(argc, argv, "n", 0);
        if (verbose) {
            cerr << "Motif size      : " << motif_size << endl;
            cerr << "Sequence window : " << window_size << endl;
            cerr << "Output          : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
            cerr << "Input           : " << (filename_input == NULL ? "stdin" : filename_input) << endl;
            cerr << "Strand specific : " << (ignore_complementary ? "yes" : "no") << endl;
            cerr << "O/E             : " << excess << endl;
            cerr << "Display         : " << num_display << endl;
            if (filename_genome != NULL) {
                cerr << "Genome          : " << filename_genome << endl;
            }
        }

        motif_counter::MODE mmode = (ignore_complementary ? motif_counter::STRAND_CONCIOUS : motif_counter::DEFAULT);
        motif_counter* counter = new motif_counter(motif_size, window_size, mmode);
        fastareader* reader = NULL;

        file_mode mode = FASTA;
        if (filename_input != NULL) {
            if (file_exists(filename_input) == false) {
                throw logic_error("cannot open input file");
            }
            string ext = get_file_extension(filename_input);
            if (ext == "bed") {
                mode = BED;
                reader = fastareader::load_genome(filename_genome);
            } else if (ext == "fa" || ext == "mfa" || ext == "fasta") {
                mode = FASTA;
            }
        }

        istream *ist = &cin;
        if (filename_input != NULL) {
            ifstream* file_in = new ifstream(filename_input);
            if (file_in->is_open() == false) {
                throw logic_error(string("cannot open ") + filename_input);
            }
            ist = file_in;
        }
        ostream *ost = &cout;
        if (filename_output != NULL) {
            ofstream* file_out = new ofstream(filename_output);
            if (file_out->is_open() == false) {
                throw logic_error(string("cannot open ") + filename_input);
            }
            ost = file_out;
        }
        
        string line;
        if (mode == BED) {
            vector<ullong> peaks;
            for (;;) {
                std::ios_base::iostate state = getline(*ist, line);
                //cout << line << endl;
                //cout << state << " " << std::ios_base::eofbit << " " << std::ios_base::failbit << " " << std::ios_base::badbit << endl;
                if ((state & (std::ios_base::eofbit | std::ios_base::failbit)) != 0) {// | std::ios_base::badbit)) != 0) {
                    break;
                }
                vector<string> items = split_items(line, '\t');
                string chromosome;
                int position = -1;
                if (items.size() == 2) {
                    chromosome = items[0];
                    position = atoi(items[1].c_str());
                } else if (items.size() >= 3) {
                    chromosome = items[0];
                    int start = atoi(items[1].c_str());
                    int end = atoi(items[2].c_str());
                    position = (start + end) / 2;
                }
                if (position > 0) {
                    int chrm = tktools::bio::convert_chromosome_to_code(items[0].c_str());
                    if (chrm > 0) {
                        peaks.push_back(((ullong)chrm << 32) | (ullong)position);
                    }
                }
            }
            sort(peaks.begin(), peaks.end());
            vector<pair<int,int> > regions;
            for (int i = 0; i < (int)peaks.size(); i++) {
                ullong prev, post, p;
                if (i > 0) {
                    prev = peaks[i-1];
                } else {
                    prev = 0;
                }
                if (i < (int)peaks.size() - 1) {
                    post = peaks[i + 1];
                } else {
                    post = 0;
                }
                p = peaks[i];
                int start, end;
                if (((prev ^ p) & 0xff00000000) != 0 || p - prev > window_size) {
                    start = (int)(p & 0xffffffffLL) - window_size;
                } else {
                    start = (p - prev) / 2;
                }
                if (((post & p) & 0xff00000000) != 0 || post - p > window_size) {
                    end = (int)(p & 0xffffffffLL) + window_size;
                } else {
                    end = (post - p) / 2;
                }
                regions.push_back(std::make_pair(start, end));
            }
            for (int i = 0; i < (int)peaks.size(); i++) {
                pair<int,int> pos = regions[i];
                string seq = reader->get_sequence(tktools::bio::convert_code_to_chromosome((int)(peaks[i] >> 32)), pos.first, pos.second);
                counter->set_sequence(seq.size(), seq.c_str(), (int)(peaks[i] & 0x8ffffff) - pos.first);
            }
        } else if (mode == FASTA) {
            string sequence;
            string name;
            for (;;) {
                std::ios_base::iostate state = getline(*ist, line);
                //cout << line << endl;
                if ((state & (std::ios_base::eofbit | std::ios_base::failbit)) != 0 || line.size() == 0) {//| std::ios_base::badbit)) != 0) {
                    break;
                }
                //cout << line.c_str()[0] << endl;
                if (line.c_str()[0] == '>') {
                    //cout << "--------------\n";
                    if (sequence.size() > 0) {
                        //cout << sequence.size() << " : " << sequence << endl;
                        //cerr << state << " : " << name << ":" << sequence.size() << "             \r";
                        counter->set_sequence(sequence.size(), sequence.c_str(), sequence.size() / 2);
                        sequence = "";
                    }
                    name = line.substr(1, line.size() - 1);
                } else {
                    sequence += line;
                    //cerr << sequence.size() << endl;
                }
            }
            if (sequence.size() > 0) {
                counter->set_sequence(sequence.size(), sequence.c_str(), sequence.size() / 2);
            }
        }

        vector<moccs_result*> results = counter->get_results(num_display, score_threshold, excess);
        for (int i = 0; i < (int)results.size(); i++) {
            *ost << (i + 1) << "\t" << results[i]->to_string() << endl;
            delete results[i];
        }

        //*ost << counter->to_string(num_display) << endl;

        delete counter;
        
        if (filename_input != NULL) {
            dynamic_cast<ifstream*>(ist)->close();
            delete ist;
        }
        if (filename_output != NULL) {
            dynamic_cast<ofstream*>(ost)->close();
            delete ost;
        }
        return 0;
    } catch (exception& e) {
        show_help();
        cerr << e.what() << endl;
        return -1;
    }
}
