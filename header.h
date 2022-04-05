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
#include <sys/time.h>
using namespace std;

#define K_MER 17

inline vector<string> split(const string &str, const char delim){
  istringstream iss(str); string tmp; vector<string> res;
  while(getline(iss, tmp, delim)) res.push_back(tmp);
  return res;
}

inline vector<string> split_str(const string& s, const string& delim){
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

double summation(vector<double> v, int start, int end){
  double sum = 0;
  for(int i = start; i < end + 1; i++){
    sum += v[i];
  }
  return sum; 
}

//20文字、20文字のpairを入れると、配列を返す
string unsl_to_str(pair<unsigned long long,unsigned long long> sequence){ //42ビット使っている
  string mer;
  unsigned long long m = 0;
  for(int i = 0; i < 20; i++)
    {
      m = (sequence.first >> (40 - i * 2)) & 3;
      switch(m)
	{
	case 0:
	  mer += 'A';
	  break;
	case 1:
	  mer += 'G';
	  break;
	case 2:
	  mer += 'C';
	  break;
	case 3:
	  mer += 'T';
	  break;
	default:
	  break;
	}
    }
  for(int i  = 0; i < 20; i++)
    {
      m = (sequence.second >> (40 - i * 2)) & 3;
      switch(m)
	{
	case 0:
	  mer += 'A';
	  break;
	case 1:
	  mer += 'G';
	  break;
	case 2:
	  mer += 'C';
	  break;
	case 3:
	  mer += 'T';
	  break;
	default:
	  break;
	}
    }
  return mer;
}

pair<unsigned long long,unsigned long long> str_to_pair(string line40){
  unsigned long long sequence1,sequence2;
  pair<unsigned long long,unsigned long long> sequence;
  for(int j = 0; j < 20; ++j)
    {
      switch(line40[j])
	{
	case 'A':
	  sequence1 = sequence1 | 0;
	  sequence1 = sequence1 << 2;
	  break;
	case 'G':
	  sequence1 = sequence1 | 1;
	  sequence1 = sequence1 << 2;
	  break;
	case 'C':
	  sequence1 = sequence1 | 2;
	  sequence1 = sequence1 << 2;
	  break;
	case 'T':
	  sequence1 = sequence1 | 3;
	  sequence1 = sequence1 << 2;
	  break;
	default:
	  break;
	}
    }
  sequence.first = sequence1;
  for(int j = 20; j < 40; ++j)
    {
      switch(line40[j])
	{
	case 'A':
	  sequence2 = sequence2 | 0;
	  sequence2 = sequence2 << 2;
	  break;
	case 'G':
	  sequence2 = sequence2 | 1;
	  sequence2 = sequence2 << 2;
	  break;
	case 'C':
	  sequence2 = sequence2 | 2;
	  sequence2 = sequence2 << 2;
	  break;
	case 'T':
	  sequence2 = sequence2 | 3;
	  sequence2 = sequence2 << 2;
	  break;
	default:
	  break;
	}
    }
  sequence.second = sequence2;
  return sequence;
}

//相補鎖のstring.ver
string complement(string str){
  string com;
    for(int j = str.size()-1; j >= 0; j--)
  //  for(int j = 0; j < str.size(); j++)
    {
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
	case 'N':
	  com = com + 'N';
	  break;
	default:
	  break;
	}
    }
  return com;
}

struct Comparer {
  bool operator() (const bitset<K_MER*2> &b1, const bitset<K_MER*2> &b2) const {
    return b1.count() < b2.count();
  }
  };


struct Comparer1 {
  bool operator() (const std::bitset<K_MER*2>& x, const std::bitset<K_MER*2>& y){
    for (int i = K_MER*2-1; i >= 0; i--) {
      if (x[i] ^ y[i]) return y[i];
    }
    return false;
  }
  };

//bit演算 
inline string bit_to_str(bitset<K_MER*2> BIT){
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

inline bitset<K_MER*2> str_to_bit(string str){
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

//input:string, output:pair<bitset,bitset> stringから相補鎖含め両方のbitが欲しい時に用いる
inline pair<bitset<K_MER*2>, bitset<K_MER*2>> str_to_bit_com(string str){
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
    BIT_COM = ((BIT << (K_MER-i)*2) >> K_MER*2-2) | (BIT_COM << 2);
  }
  BIT_P.second = BIT_COM;
  return BIT_P;
}

//bitを相補鎖にする input:bit output:bit
bitset<K_MER*2> bit_to_bit_com(bitset<K_MER*2> bit_k){
  bitset<K_MER*2> bit_k_com;
  bit_k = ~bit_k;
  for(int i = 1; i <= K_MER; i++){
    bit_k_com = ((bit_k << (K_MER-i)*2) >> K_MER*2-2) | (bit_k_com << 2);
  }  
  return bit_k_com;
}


//intput:vector<int> output:average
double vector_to_average(vector<double> vec){
  double sum;
  for(int i = 0; i < vec.size(); i++){
    sum += vec[i];
  }
  return sum/vec.size();
}

//intput:vector<int> output:max
double vector_to_max(vector<double> vec){
  double num = vec[0];
  for(int i = 1; i < vec.size(); i++){
    if(num < vec[i]){
      num = vec[i];
    }
  }
  return num;
}

//intput:vector<int> output:min
double vector_to_min(vector<double> vec){
  double num = vec[0];
  for(int i = 1; i < vec.size(); i++){
    if(num > vec[i]){
      num = vec[i];
    }
  }
  return num;
}

inline string insert_newline(string &str){ //80文字で改行入れる
  int c = 70;
  string str_out;
  for(int i = 0; i < str.size()/c; ++i){
    str_out += str.substr(i*c,c) + "\n";
  }
  str_out += str.substr((str.size()/c)*c,str.size()-(str.size()/c)*c);
  if(str.size()-(str.size()/c)*c == 0){ //c以下の最後のつけたしがなければ（cの倍数なら）最後の改行文字を削る
    str_out = str_out.substr(0,str_out.size()-1);
  }
  return str_out;
}

inline bool alignment(string &str1, string &str2, double perc){ //db:str1, query:str2 str1の方が長さは長い想定 percに総長の割合のスコアの閾値を設ける
  int score = 0;
  int tuple_size = 32; //word_size
  int match = 1; //match_score
  int mismatch = -3; //mis match score
  int start = 0;
  string k_tuple = str2.substr(0,tuple_size);
  bool conclusion = false;
  if(str1.find(k_tuple) != string::npos){ //k_tupleが見つかれば
    start = str1.find(k_tuple);
    for(int i = 0; i < str2.length(); ++i){ //1文字ずつ見ていく
      if(str1[start+i] == str2[i]){ //matchしてたらスコアをあげる
	score += match;
      } else {
	score += mismatch; //mismatchならスコアを下げる
      }
    }
    if(score > str2.size()*perc){
      conclusion = true;
    }
  }
  return conclusion;
}

inline bool alignment_smith(string &str1, string &str2, double perc){
  int tuple_size = 32; //word_size
  bool result = false;
  int start = 0;
  int step_size = 10;
  while(result == false){
    string k_tuple = str2.substr(start,tuple_size);
    if(str1.find(k_tuple) != string::npos){ //k_tupleが見つかれば
      int sz1 = str1.size()+1;
      int sz2 = str2.size()+1;
      int a[sz1][sz2];
      for(int i = 0; i < sz1; ++i){
	for(int j = 0; j < sz2; ++j){
	  a[i][j] = 0;
	}
      }
      int score_match = 3,score_mis = 2,score_gap = 1;
      for(int i = 1; i < sz1; ++i){
	for(int j = 1; j < sz2; ++j){
	  if(str1[i] == str2[j]){
	    a[i][j] = a[i-1][j-1] + score_match;
	  } else if(a[i-1][j] != 0 || a[i][j-1] != 0 || a[i-1][j-1] != 0){
	    if(a[i-1][j] >= a[i-1][j-1] && a[i-1][j] >= a[i][j-1] && a[i-1][j] > score_gap){ //gap
	      a[i][j] = a[i-1][j-1] - score_mis;
	    } else if(a[i][j-1] >= a[i-1][j-1] && a[i][j-1] >= a[i-1][j] && a[i][j-1] > score_gap){ //gap
	      a[i][j] = a[i][j-1] - score_gap;
	    } else if(a[i-1][j-1] >= a[i][j-1] && a[i-1][j-1] >= a[i-1][j] && a[i-1][j-1] > score_mis){ //mismatch
	      a[i][j] = a[i-1][j] - score_gap;
	    }
	  }
	}
      }
      int max = 0,max_pos[2];
      for(int i = 0; i < sz1; ++i){
	for(int j = 0; j < sz2; ++j){
	  if(max < a[i][j]){
	    max = a[i][j];
	    max_pos[0] = i,max_pos[1] = j;
	  }
	}
      }
      int i = max_pos[0],j = max_pos[1];
      while(a[i][j] != 0){
	if(a[i-1][j-1] >= a[i][j-1] && a[i-1][j-1] >= a[i-1][j]){
	  --i,--j;
	}else if(a[i-1][j] >= a[i][j-1] && a[i-1][j] >= a[i-1][j-1]){
	  --i;
	} else if(a[i][j-1] >= a[i-1][j] && a[i][j-1] >= a[i-1][j-1]){
	  --j;
	}
      }
      if(max_pos[1]-j >= str2.size()*perc){
	result = true;
      }
    }
    start += step_size;
    if(start > str2.size() * (1-perc)){
      break;
    }
  }
  return result;
}

//k-merの配列から1塩基違いの配列をvectorにpush_back(k*3個)
inline void snp_kmer(string &str, vector<string> &v){
  bitset<K_MER*2> bit,test_bit,test_bit_c,bit_out1,bit_out2,bit_out3;
  test_bit = str_to_bit(str);
  test_bit_c = test_bit; //固定
  for(int i = 1; i <= K_MER; ++i){
    test_bit = test_bit_c;
    bit = ((test_bit << (i - 1)*2 ) >> (K_MER-1)*2);
    if(bit == 0){ //00
      //      bit_out1 = test_bit.set((K_MER-i)*2 , true); //01
      //      bit_out2 = test_bit.set((K_MER-i)*2 +1, true); //11
      //      bit_out3 = test_bit.set((K_MER-i)*2 , false); //10
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 , true)) );
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 +1, true) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 , false) ));
    } else if(bit == 1){ //01
      //      bit_out1 = test_bit.set((K_MER-i)*2 , false); //00
      //      bit_out2 = test_bit.set((K_MER-i)*2 + 1 , true); //10
      //      bit_out3 = test_bit.set((K_MER-i)*2 , true); //11
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 , false) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 + 1 , true) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 , true) ));
    } else if(bit == 2){ //10
      //      bit_out1 = test_bit.set((K_MER-i)*2 + 1 , false); //00
      //      bit_out2 = test_bit.set((K_MER-i)*2 , true); //01
      //      bit_out3 = test_bit.set((K_MER-i)*2 +1, true); //11
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 + 1 , false) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 , true) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 +1, true) ));
    } else if(bit == 3){ //11
      //      bit_out3 = test_bit.set((K_MER-i)*2 + 1 , false); //01
      //      bit_out1 = test_bit.set((K_MER-i)*2 , false); //00
      //      bit_out2 = test_bit.set((K_MER-i)*2 + 1 , true); //10
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 + 1 , false) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 , false) ));
      v.push_back(bit_to_str( test_bit.set((K_MER-i)*2 + 1 , true) ));
    }
  }
}

//k-merの配列から2塩基違いの配列をvectorにpush_back(kC2 * 3 * 3個)
inline void snp2_kmer(string &str, vector<string> &v){
  bitset<K_MER*2> bit,bit2,test_bit,test_bit2,test_bit_c,bit_out1,bit_out2,bit_out3;
  test_bit = str_to_bit(str);
  test_bit_c = test_bit;
  for(int i = 1; i <= K_MER; ++i){
    test_bit = test_bit_c;
    bit = ((test_bit << (i - 1)*2 ) >> (K_MER-1)*2);
    if(bit == 0){ //00
      bit_out1 = test_bit.set((K_MER-i)*2 , true); //01
      bit_out2 = test_bit.set((K_MER-i)*2 +1, true); //11
      bit_out3 = test_bit.set((K_MER-i)*2 , false); //10
    } else if(bit == 1){ //01
      bit_out1 = test_bit.set((K_MER-i)*2 , false); //00
      bit_out2 = test_bit.set((K_MER-i)*2 + 1 , true); //10
      bit_out3 = test_bit.set((K_MER-i)*2 , true); //11
    } else if(bit == 2){ //10
      bit_out1 = test_bit.set((K_MER-i)*2 + 1 , false); //00
      bit_out2 = test_bit.set((K_MER-i)*2 , true); //01
      bit_out3 = test_bit.set((K_MER-i)*2 +1, true); //11
    } else if(bit == 3){ //11
      bit_out1 = test_bit.set((K_MER-i)*2 + 1 , false); //01
      bit_out2 = test_bit.set((K_MER-i)*2 , false); //00
      bit_out3 = test_bit.set((K_MER-i)*2 + 1 , true); //10
    }
    bitset<K_MER*2> bit_list[3] = {bit_out1,bit_out2,bit_out3};
    for(int l = 0; l < 3; ++l){
      for(int j = i + 1; j <= K_MER; ++j){
	test_bit2 = bit_list[l];
	bit2 = ((test_bit << (j - 1)*2 ) >> (K_MER-1)*2);
	if(bit2 == 0){ //00
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 , true)) ); //01
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 +1, true)) ); //11
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 , false)) ); //10
	} else if(bit2 == 1){ //01
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 , false)) ); //00
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 + 1 , true)) ); //10
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 , true)) ); //11
	} else if(bit2 == 2){ //10
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 + 1 , false)) ); //00
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 , true)) ); //01
	  v.push_back( bit_to_str(test_bit.set((K_MER-j)*2 +1, true)) ); //11
	} else if(bit2 == 3){ //11
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 + 1 , false)) ); //01
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 , false)) ); //00
	  v.push_back( bit_to_str(test_bit2.set((K_MER-j)*2 + 1 , true)) ); //10
	}
      }

    }
  }
}
