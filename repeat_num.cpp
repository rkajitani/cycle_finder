#include "cycle.h"
#include "common.h"

namespace std{
  template <>
  class hash<std::bitset<K_MER*2> >
  {
  public:
    size_t operator()(const std::bitset<K_MER*2>& x) const
    {
      size_t value=0;
      for(int i = 0; i < K_MER*2; ++i){
        if(x[i] == 1){
          value = ((value << 1) | 1);
        } else {
          value = value << 1;
        }
      }
      return value;
    }
  };
}

void Cycle::repeat_num_exe(){
  int repeat_type = this->repeat_type;
  string file_repeat = this->o;
  int method = this->method;
  string file_kmer = this->fa;
  bool data_type = this->d;
  int PEAK = 0;
  int PEAK1 = 0,PEAK2 = 0;
  string output;
  if(method == 1){
    output = this->repeat_num_out;
  } else if(method == 0){
    output = this->repeat_num_out_min;
  } else {
    cout << "error in repeat_num.cpp" << endl;
  }
  if(this->d == 1){
    PEAK = this->p;
  } else if(this->d == 0){
    PEAK1 = this->p1;
    PEAK2 = this->p2;
  } else {
    cerr << "error : data type is incorrect" << endl;
    return;
  }
  FILE *fp1,*fp2;
  char str[16384];
  ifstream del, norm;
  vector<int> n,list;
  unordered_multimap<string, pair<unsigned long long, unsigned long long> > mp;
  string line, open_norm, open_del, name, line40,FILE,column1,column2,fasta;
  pair<unsigned long long, unsigned long long> sequence, key;
  vector<bitset<K_MER*2> >SEQUENCE;
  fp1 = fopen(file_kmer.c_str(),"r");
  fp2 = fopen(file_repeat.c_str(),"r");
  if(fp1 == NULL || fp2 == NULL || (data_type != 1 && data_type != 0)){
    if(fp1 == NULL){
      cerr << "can't open  " << file_kmer << endl;
    } else if(fp2 == NULL){
      cerr << "can't open " << file_repeat << endl;
    } else {
      cerr << "data type is wrong" << endl;
    }
    printf("CANTOPEN\n");
    return;
  }

  //k-mer fastaファイル
  ofstream ofs(output);
  unsigned long long number_s = 0,number_t = 0,number = 0;
  unordered_map<bitset<K_MER*2>, int> NORM, DEL;
  vector<string> a;
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>' && data_type == 0){
	  number_t = stol( line.substr(7,line.find("Somatic")-7));
	  number_s = stol( line.substr(line.find("Somatic")+7,line.size()-6-line.find("Somatic") ) );
	  //	  cout << line << "\t" << number_t << "\t" << number_s << endl;
    }
    if(line[0] == '>' && data_type == 1){
      number = stol(line.substr(1, line.size()-1));
    }
    else if(data_type == 0){
      NORM[str_to_bit(line)] =  number_t;
      NORM[bit_to_bit_com(str_to_bit(line))] =  number_t;      
      DEL[str_to_bit(line)] = number_s;
      DEL[bit_to_bit_com(str_to_bit(line))] = number_s;
    }
    else if(data_type == 1){
      NORM[str_to_bit(line)] = number;
      NORM[bit_to_bit_com(str_to_bit(line))] =  number;
    }
  }
  fclose(fp1);

  //fastaファイル
  vector<string> est;
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] != '>'){
      est.push_back(line);
    }
  }
  fclose(fp2);

  string mer;
  long long numsum = 0,num = 0,numsums = 0,nummin=0,nummins=0,est_size;
  bitset<K_MER*2> bit_kmer,bit_kmer_com;
  vector<string> replace(3);
  pair<int,int> pair_result;
  unordered_map<string, pair<int, int> > mp_result;
  for(unsigned z = 0; z < est.size(); ++z){ //各リピートについて
    line.clear();
    if(repeat_type == 1){ //tandem repeatなら直列につなげる
      if(est[z].size() < K_MER){ //k文字以下のとき
	for(unsigned i_copy = 0; i_copy <  (K_MER / est[z].size())+2; i_copy++){ //k文字以上になるように足す
	  line += est[z]; 
	}
      } else { 
	line = est[z] + est[z];
      }
      est_size = est[z].size(); //部分文字列は文字数分
    } else { //intersperse repeatなら
      line = est[z];
      est_size = est[z].size() - K_MER + 1; //部分文字列はk-merがとれるだけ
    }
    for(int i = 0; i < est_size; ++i){ //lineからk文字部分配列
      SEQUENCE.push_back(str_to_bit(line.substr(i, K_MER)));
    }
    for(unsigned j = 0; j < SEQUENCE.size(); j++){ //それぞれのk-merの出現回数を調べる
      bit_kmer = SEQUENCE[j];
      numsum += NORM[bit_kmer];
      numsums += DEL[bit_kmer];
      if(nummin > NORM[bit_kmer] || nummin == 0){ //最小値になりうる値なら
	nummin = NORM[bit_kmer];
      }
      num++;
    }
    if(method == 0){
      pair_result.first = nummin;
      pair_result.second = nummins;
    } else {
      pair_result.first = numsum/SEQUENCE.size();
      pair_result.second = numsums/SEQUENCE.size();
    }
    mp_result[est[z]] = pair_result;
    SEQUENCE.clear();
    numsum = 0,num = 0,numsums = 0,nummin=0,nummins=0;
  }
  unordered_multimap<int,string> mp_sort; //長さをkey、配列をvalueに
  unordered_set<int> mp_set; //長さ格納
  string repeat;
  for(auto it = mp_result.begin(); it != mp_result.end(); ++it){
    mp_sort.insert(make_pair((it->first).size(),it->first));
    mp_set.insert(it->first.size());
  }
  vector<int> v_sort;
  for(auto it = mp_set.begin(); it != mp_set.end(); ++it){
    v_sort.push_back(*it);
  }
  sort(v_sort.begin(),v_sort.end());
  for(auto it = v_sort.begin(); it != v_sort.end(); ++it){
    auto range = mp_sort.equal_range(*it);
    for(auto it_multi = range.first; it_multi != range.second; ++it_multi){
      repeat = it_multi->second;
      nummin = mp_result[repeat].first;
      nummins = mp_result[repeat].second;
      if(data_type == 0){
	ofs << repeat << "\t" << repeat.size() << "\t" << "copy1_copy2_dif" << "\t" << (nummin / PEAK1) - (nummins / PEAK2) << "\t" << nummin/PEAK1 << "\t" << nummins/PEAK2  << "\t" << "base1_base2_dif" << "\t" << repeat.size() * nummin/PEAK1 - repeat.size()*nummins/PEAK2 << "\t" << repeat.size()*nummin/PEAK1 << "\t" <<  repeat.size()*nummins/PEAK2 << endl;
      } else if(data_type == 1){
	ofs << repeat << "\t" << repeat.size() << "\t" << "copy1_copy2_dif" << "\t" << (nummin / PEAK) << "\t" << nummin / PEAK << "\t" << 0 << "\t" << "base1_base2_dif" << "\t" << repeat.size() * nummin/PEAK  << "\t" << repeat.size()*nummin/PEAK << "\t" <<  0 << "\t" << fixed << setprecision(1)<< 0  <<  endl;	
      }
    }
  }
}

void Cycle::repeat_num(){
  this->method = 1;
  repeat_num_exe();
  this->method = 0;
  repeat_num_exe();
}
