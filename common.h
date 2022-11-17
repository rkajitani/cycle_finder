#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <bitset>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <utility>
#include <list>
#include <iomanip>
#include <stdio.h>
#include <array>
#include <omp.h>
#include <cstdint>
#include <chrono>
using namespace std;

#define PROGRAM_NAME "cycle_finder"
#define VERSION "1.0.0"
#define K_MER 17
#define ROOT_PATH "/data1/kajitani/misc/cycle_finder_v1.0/cycle_finder"
#define CDHIT_IDENTITY_THR 0.8  //cd-hitで用いられるidentity
#define CDHIT_COVERAGE_THR_TANDEM 0.7  //全長に対してアライメント長が超えるべき割合(tandemでつないでいるとき)
#define CDHIT_COVERAGE_THR_NOT_TANDEM 0.8  //(tandemでつないでないとき)



vector<string> split(const string &str, const char delim);
bitset<K_MER*2> str_to_bit(const string &str);
string bit_to_str(bitset<K_MER*2> BIT);
vector<string> split_str(const string& s, const string& delim);
pair<bitset<K_MER*2>, bitset<K_MER*2>> str_to_bit_com(string str);
bitset<K_MER*2> bit_to_bit_com(bitset<K_MER*2> bit_k);
string complement(string str);
uint64_t file_exist(string file);
void file_write(string file, string log);
vector<string> option_multi_file(int argc,char *argv[],int *num);
void makeblastdb(string reference, string log_file);
void stopwatch(std::chrono::system_clock::time_point start, bool out);

