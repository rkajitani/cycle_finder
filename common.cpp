#include "common.h"

#ifndef COMMON_H
#define COMMON_H


vector<string> split(const string &str, const char delim){
  istringstream iss(str); string tmp; vector<string> res;
  while(getline(iss, tmp, delim)) res.push_back(tmp);
  return res;
}

bitset<K_MER*2> str_to_bit(const string &str){
  bitset<K_MER*2> BIT;
  for(int i = 0; i < K_MER; i++){
    switch(str[i]){
    case 'A':
      BIT <<= 2;
      break;
    case 'G':
      (BIT <<= 2) |=1;
      break;
    case 'C':
      (BIT <<= 2) |=2;
      break;
    case 'T':
      (BIT <<= 2) |=3;
      break;
    default:
      break;
    }
  }
  return BIT;
}

string bit_to_str(bitset<K_MER*2> BIT){
  string str;
  bitset<K_MER*2> BIT_tmp;
  for(int i = 0; i < K_MER; i++){
    BIT_tmp = BIT;
    BIT_tmp <<= i*2;
    BIT_tmp >>= (K_MER-1)*2;
    switch(BIT_tmp.to_ulong()){
    case 0:
      str+="A";
      break;
    case 1:
      str+="G";
      break;
    case 2:
      str+="C";
      break;
    case 3:
      str+="T";
      break;
    default:
      cout << "ERROR" << endl;
      break;
    }
  }
  return str;
}

vector<string> split_str(const string& s, const string& delim){
  vector<string> result;
  result.clear();
  size_t pos = 0;
  while(pos != string::npos){
    size_t p = s.find(delim, pos);
    if(p == string::npos){
      result.push_back(s.substr(pos));
      break;
    } else {
      result.push_back(s.substr(pos, p - pos));
    }
    pos = p + delim.size();
  }
  return result;
}

pair<bitset<K_MER*2>, bitset<K_MER*2>> str_to_bit_com(string str){
  pair<bitset<K_MER*2>,bitset<K_MER*2>> BIT_P;
  bitset<K_MER*2> BIT,BIT_COM;
  for(int i = 0; i < K_MER; i++){
    switch(str[i]){
    case 'A':
      BIT <<= 2;
      break;
    case 'G':
      (BIT <<= 2) |=1;
      break;
    case 'C':
      (BIT <<= 2) |=2;
      break;
    case 'T':
      (BIT <<= 2) |=3;
      break;
    default:
      break;
    }
  }
  BIT_P.first = BIT;
  BIT = ~BIT;
  for(int i = 1; i <= K_MER; i++){
    BIT_COM = ((BIT << ((K_MER - i) * 2)) >> (K_MER * 2 - 2)) | (BIT_COM << 2);
  }
  BIT_P.second = BIT_COM;
  return BIT_P;
}

bitset<K_MER*2> bit_to_bit_com(bitset<K_MER*2> bit_k){
  bitset<K_MER*2> bit_k_com;
  bit_k = ~bit_k;
  for(int i = 1; i <= K_MER; i++){
    bit_k_com = ((bit_k << ((K_MER-i)*2)) >> (K_MER*2-2)) | (bit_k_com << 2);
  }
  return bit_k_com;
}

string complement(string str){
  string com;
  for(int j = str.size()-1; j >= 0; j--){
    switch(str[j])
      {
      case  'A':
        com = com + 'T';
        break;
      case 'G':
        com = com + 'C';
        break;
      case 'C':
        com = com + 'G';
        break;
      case 'T':
        com = com + 'A';
        break;
      default:
        break;
      }
  }
  return com;
}

uint64_t file_exist(string file){
  char buff[1024];
  uint64_t file_size;
  FILE *fp;
  fp = fopen(file.c_str(),"r");
  if(fp == NULL){
    fclose(fp);
    return 0;
  }
  while(fgets(buff, sizeof(buff), fp)){
    file_size = strlen(buff);
    if(file_size > 0){
      break;
    }
  }
  if(file_size == 0){
    fclose(fp);
    return 0;
  } else {
    fclose(fp);
    return 1;
  }
}

void file_write(string file, string log){
  ofstream ofs(file, ios::app);
  ofs << log << endl;
}

vector<string> option_multi_file(int argc,char *argv[],int *num){
  vector<string> files;
  for(int i = *num + 1; i <= argc; ++i){
    ++*num;
    if(i == argc){ //オプションの最後が-fだったときのために
      break;
    }
    if(*argv[i] == '-'){
      break;
    }
    files.push_back(argv[i]);
  }
  return files;
}

void makeblastdb(string reference, string log_file){
  string root_path = ROOT_PATH;
  string cmd = root_path + "/makeblastdb -in " + reference + " -dbtype nucl -hash_index >> " + log_file;
  system(cmd.c_str());
}

void stopwatch(std::chrono::system_clock::time_point start, bool out){
  auto end = std::chrono::system_clock::now();
  auto dur = end - start;
  auto sec = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
  int min, sec_out, min_out, hour_out;
  int SIXTY = 60;
  min = sec / SIXTY;
  sec_out = sec % SIXTY;
  hour_out = min / SIXTY;
  min_out = min % SIXTY;
  if(out) std::cout << hour_out << " h " << min_out << " min " << sec_out << " sec\n";
}


#endif
