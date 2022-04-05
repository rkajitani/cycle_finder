#include "cluster.h"

namespace std {
  template <>
  class hash<std::pair<int, int>> {
  public:
    size_t operator()(const std::pair<int, int>& x) const{
      return hash<int>()(x.first) ^ hash<int>()(x.second);
    }
  };
}

void Cluster::blast_to_copynum(){
  int PEAK1,PEAK2,PEAK;
  ifstream file1,file2,file3;
  int data_type = 0;
  unordered_map<int,double> mp_c,mp_m,mp_min;
  unordered_map<int,double> mp_t,mp_s;
  string line,fasta,result,clstr,column1,column2,column3,column4,column5,column6,consensus,column2_pre,clstr2,column7,column8,column9,column10,FILE1,FILE2,FILE3;
  string output;
  output = this->o + "_blst_family";
  FILE1 = this->o + "_blst.blastn";
  FILE2 = this->f3; //result_min
  FILE3 = this->o + "_clst.fa";
  data_type = this->d;
  PEAK1 = this->p1;
  PEAK2 = this->p2;
  PEAK = this->p;
  if(data_type == 0){
    cerr << "blast_to_copy_number2.cpp -f1 " << FILE1 << " -f2 " << FILE2 << " -f3 " << FILE3 << " -t " << data_type << " -p1 " << PEAK1 << " -p2 " << PEAK2 << " -o " << output << endl;
  } else if(data_type == 1){
    //    cerr << "blast_to_copy_number2.cpp -f1 " << FILE1 << " -f2 " << FILE2 << " -f3 " << FILE3 << " -t " << data_type << " -p " << PEAK << " -o " << output << endl;
  }

  //fastaファイルopen
  vector<string> a,b,c;
  double appear_t = 0, appear_s = 0,appear_x=0;
  long long family_id = 0,family_id_pre = -1,repeat_id,counter = 0,sum_min_t,sum_min_s;
  double score_tmp_pre = 0, score_tmp;
  vector<int> Family_id,Med_tmp;
  unordered_multiset<string> Column2;
  unordered_map<int,pair<long long,long long> > mp_rep2freq;
  unordered_map<int,double> mp_rep2sc_max;
  unordered_map<int,long long> mp_rep2sum_t,mp_rep2sum_s;
  map<pair<int,int>,double > mp_famrep2sc; //unordered_mapだとkeyがpairのhash関数がないのかコンパイルエラーがでる
  unordered_multimap<int,int> mp_fam2rep,mp_rep2fam;
  unordered_map<string,int> mp_str2cp_t,mp_str2cp_s;
  unordered_map<int,string> mp_fam2ref;
  pair<int,int> pr_freq,pr_famrep;
  unordered_set<int> set_tes;
  double score;
  FILE *fp1 = fopen(FILE1.c_str(),"r");
  FILE *fp2 = fopen(FILE2.c_str(),"r");
  FILE *fp3 = fopen(FILE3.c_str(),"r");
  char str[16384];
  unordered_set<int> set_dupli,set_dupli_,set_rep;

  //resultファイルから key:配列, value:コピー数
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    mp_str2cp_t.insert(make_pair(a[0],stoi(a[4])));
    mp_str2cp_s.insert(make_pair(a[0],stoi(a[5])));
  }
  //  cout << "result読み込み完了" << endl;
  //clstファイルから key:family_id, value:代表配列
  while((fgets(str,sizeof(str),fp3)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line.find("*") != string::npos){ //代表配列だけ格納
      b = split(line,'_');
      fgets(str,sizeof(str),fp3);	
      str[strlen(str)-1] = '\0';
      line = str;
      mp_fam2ref.insert(make_pair(stoi(b[1]),line)); //key:family_id, value:配列
    }
  }
  //  cout << "clst読み込み完了" << endl;
  unordered_set<pair<int,int> > mp_fam2rep_set;
  pair<int,int> mp_fam2rep_pr;
  //blastnファイル
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    b = split(a[0],'_'); //b[1]:family_id
    c = split(a[1],'_'); //c[1]:repeat_id
    family_id = stoi(b[1]);
    repeat_id = stoi(c[1]);
    //changed
    mp_fam2rep_pr.first = family_id; 
    mp_fam2rep_pr.second = repeat_id;
    mp_fam2rep_set.insert(mp_fam2rep_pr); //pairをsetに格納して、重複を除去
    if(data_type == 0){
      pr_freq.first = stoll(a[1].substr(6,a[1].find("Somatic")-6));
      pr_freq.second = stoll(a[1].substr(a[1].find("Somatic")+7,a[1].find("_")-a[1].find("Somatic")-7));
    } else if(data_type == 1){
      pr_freq.first = stoll(a[1].substr(0,a[1].find("_")));
      pr_freq.second = 0;
    }
    mp_rep2freq.insert(make_pair(repeat_id,pr_freq)); //key:repeat_id, value:pair(frequency)
    pr_famrep.first = family_id;
    pr_famrep.second = repeat_id;
    //    score = stod(a[11]);
    score = stod(a[2]);
    mp_famrep2sc.insert(make_pair(pr_famrep,score)); //key:(family_id,repeat_id), value:bit_score
    //重複k-mer格納
    if(set_rep.find(repeat_id) != set_rep.end()){ //repeat_idがこれまでに出てきていれば(重複k-mer)
      set_dupli.insert(repeat_id);
    }
    set_rep.insert(repeat_id);
    //family_idが変われば 長さも格納
    if(family_id != family_id_pre){ 
      family_id_pre = family_id;
      Family_id.push_back(family_id);
      mp_c[family_id] = stoi(b[2]); //とりあえず1つめのものを代表配列の*がなくなってた時用に格納する
    }
    //代表配列の長さを格納
    if(a[0].find("*") != string::npos){
      mp_c[family_id] = stoi(b[2]); //key:family_id, value:repeat_length
    }
  }
  //  cout << set_rep.size() << "\t" << set_dupli.size() << endl;
  //changed
  for(auto it = mp_fam2rep_set.begin(); it != mp_fam2rep_set.end(); ++it){
    mp_fam2rep_pr = *it;
    mp_fam2rep.insert(make_pair(mp_fam2rep_pr.first,mp_fam2rep_pr.second));
  }

  for(auto it = set_dupli.begin(); it != set_dupli.end(); ++it){ //ユニークk-merセットを作成
    set_rep.erase(*it);
  }
  //  cout << set_rep.size() << "\t" << set_dupli.size() + set_rep.size() << endl;
  
  //重複k-merからfamily_idを呼び出せるように
  unordered_set<pair<int,int> > mp_rep2fam_set;
  pair<int,int> mp_rep2fam_pr;
  for(unsigned i = 0; i < Family_id.size(); ++i){ //family_id全てまわす
    family_id = Family_id[i];
    mp_rep2fam_pr.second = family_id;
    auto range = mp_fam2rep.equal_range(family_id); //repeat_id全てまわす(multimap)
    for(auto it = range.first; it != range.second; ++it){
      repeat_id = it->second;
      if(set_dupli.find(repeat_id) != set_dupli.end()){ //重複k-mer 
	mp_rep2fam_pr.first = repeat_id;
	mp_rep2fam_set.insert(mp_rep2fam_pr); //pairをsetに格納して、重複を除去
      } else { //ユニークk-merのコピー数を足す
	appear_t += mp_rep2freq[repeat_id].first;
	appear_s += mp_rep2freq[repeat_id].second;
	//	appear_x += mp_rep2freq[repeat_id].first/PEAK1 - mp_rep2freq[repeat_id].second/PEAK2;
	++counter;
	set_tes.insert(repeat_id); //hoge
      }
    }
    //family_id毎にユニークk-mer由来のfrequencyを格納
    mp_t[family_id] = appear_t; 
    mp_s[family_id] = appear_s;
    appear_t = 0;
    appear_s = 0;
  }
  //  cout << mp_fam2rep.size() << endl;
  //changed
  for(auto it = mp_rep2fam_set.begin(); it != mp_rep2fam_set.end(); ++it){
    mp_rep2fam_pr = *it;
    mp_rep2fam.insert(make_pair(mp_rep2fam_pr.first,mp_rep2fam_pr.second));
  }
  
  //重複k-merをもつrepeat_idからfamily_idと最大scoreを出せるように    
  long long repeat_id_max = 0;
  double score_tmp_max = 0;
  for(auto it = set_dupli.begin(); it != set_dupli.end(); ++it){ //重複k-mer全体をまわす
    repeat_id = *it;
    auto range = mp_rep2fam.equal_range(repeat_id); //重複k-merをkeyにfamily_id全体をまわす
    for(auto it1 = range.first; it1 != range.second; ++it1){
      family_id = it1->second;
      pr_famrep.first = family_id;
      pr_famrep.second = repeat_id;
      score_tmp = mp_famrep2sc[pr_famrep];
      if(score_tmp > score_tmp_pre){ //重複k-merをbit-scoreが最大のものを特定
	repeat_id_max = repeat_id;
	score_tmp_max = score_tmp;
	score_tmp_pre = score_tmp;
      }
    }
    score_tmp_pre = 0;
    //最後にmaxのscoreのものを格納
    mp_rep2sc_max.insert(make_pair(repeat_id_max,score_tmp_max)); //repeat_idの中で最大のscoreを格納
  }
  //  cout << "重複k-merの最大score" << "\t" << mp_rep2sc_max.size() << endl;
  
  //重複k-merを回して、分配する
  for(auto it = set_dupli.begin(); it != set_dupli.end(); ++it){ //重複k-mer全体をまわす
    repeat_id = *it;
    auto range = mp_rep2fam.equal_range(repeat_id);
    for(auto it1 = range.first; it1 != range.second; ++it1){ //family全体をまわす
      family_id = it1->second;
      pr_famrep.first = family_id;
      pr_famrep.second = repeat_id;
      if(mp_famrep2sc[pr_famrep] == mp_rep2sc_max[repeat_id]){ //scoreがmaxのfamilyに絞ってminの総和を求める
	sum_min_t += mp_str2cp_t[mp_fam2ref[family_id]]; //family_id -> 代表配列 -> copy数(result_min)
	sum_min_s += mp_str2cp_s[mp_fam2ref[family_id]];
      }
    }
    if(sum_min_t == 0){
      cerr << "*********error:testisでコピー数0を出力しています*****************"  << endl;
      exit(1);
    }
    mp_rep2sum_t.insert(make_pair(repeat_id,sum_min_t)); //repeat_idをkeyにmin総和を格納
    if(sum_min_s == 0){ //somaは0のときがある。このときは全て0なので、分母をとりあえず0以外にしておく
      sum_min_s = 1;
    }
    mp_rep2sum_s.insert(make_pair(repeat_id,sum_min_s)); //repeat_idをkeyにmin総和を格納
    sum_min_t = 0;
    sum_min_s = 0;
  }

  unordered_set<int> set_dupli_done;
  for(auto it = set_dupli.begin(); it != set_dupli.end(); ++it){ //もう一度重複k-mer全体をまわす
    repeat_id = *it;
    auto range = mp_rep2fam.equal_range(repeat_id);
    for(auto it1 = range.first; it1 != range.second; ++it1){ //もう一度family全体をまわす
      family_id = it1->second;
      pr_famrep.first = family_id;
      pr_famrep.second = repeat_id;
      if(mp_famrep2sc[pr_famrep] == mp_rep2sc_max[repeat_id]){ //scoreがmaxのfamilyに絞る
	mp_t[family_id] += (double(mp_rep2freq[repeat_id].first * mp_str2cp_t[mp_fam2ref[family_id]]) / mp_rep2sum_t[repeat_id]);
	mp_s[family_id] += (double(mp_rep2freq[repeat_id].second * mp_str2cp_t[mp_fam2ref[family_id]]) / mp_rep2sum_t[repeat_id]);
	//appear_x += (double(mp_rep2freq[repeat_id].first * mp_str2cp_t[mp_fam2ref[family_id]]) / mp_rep2sum_t[repeat_id]);
	appear_x += (double)mp_rep2freq[repeat_id].first * mp_str2cp_t[mp_fam2ref[family_id]] / (mp_rep2sum_t[repeat_id]*PEAK1) - double(mp_rep2freq[repeat_id].second * mp_str2cp_t[mp_fam2ref[family_id]]) / (mp_rep2sum_t[repeat_id] * PEAK2);
	//	cout << double(mp_rep2freq[repeat_id].first * mp_str2cp_t[mp_fam2ref[family_id]]) / mp_rep2sum_t[repeat_id] << "\t" << mp_rep2freq[repeat_id].first << "\t" << mp_str2cp_t[mp_fam2ref[family_id]] << "\t" << mp_rep2sum_t[repeat_id] << endl;
	set_dupli_done.insert(repeat_id);
      }
    }
  }
  for(auto it = set_dupli_done.begin(); it != set_dupli_done.end(); ++it){
    set_dupli.erase(*it);
  }
  //  cerr << "出力されてないk-merは" << set_dupli.size() << endl;
  //  cerr << uint64_t(appear_x) << endl;
  
  sort(Family_id.begin(),Family_id.end());
  ofstream ofs(output);
  if(data_type == 0){
    for(unsigned i = 0; i < Family_id.size(); i++){
      ofs << "family_" << Family_id[i] << "_" << mp_c[Family_id[i]] << "\t" << uint64_t(round(mp_t[Family_id[i]]/PEAK1 - mp_s[Family_id[i]]/PEAK2)) << "\t" << uint64_t(round(mp_t[Family_id[i]]/PEAK1)) << "\t" << uint64_t(round(mp_s[Family_id[i]]/PEAK2))  << "\t" << uint64_t(round((mp_t[Family_id[i]]/PEAK1 - mp_s[Family_id[i]]/PEAK2)/mp_c[Family_id[i]])) << endl;
    }
  } else if(data_type == 1){
    for(unsigned i = 0; i < Family_id.size(); i++){
      ofs << "family_" << Family_id[i] << "_" << mp_c[Family_id[i]] << "\t" <<  uint64_t(round(mp_t[Family_id[i]]/PEAK)) << "\t" << uint64_t(round(mp_t[Family_id[i]]/PEAK)) << "\t" << 0 << "\t" << uint64_t(round((mp_t[Family_id[i]]/PEAK)/mp_c[Family_id[i]])) << endl;
    }    
  }
}
