#include "cycle.h"

#define LENGTH_THR_STR 2000
#define COVERAGE 0.5
#define READ_LEN 100

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

//pairのsecondでsortできるように
template <class T1, class T2, class Pred = std::less<T2> >
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
    Pred p;
    return p(left.second, right.second);
  }
};

//intput:bitset<K_MER*2>,map, output:array<bitset<K_MER*2>> ;kmer入れるとmapに含まれる次のkmer全てを格納したvectorを出力
inline vector<bitset<K_MER*2> > next_kmer(const bitset<K_MER*2> bit_k,const unordered_map<bitset<K_MER*2>,int>  mp_k){
  vector<bitset<K_MER*2> > Next;
  bitset<K_MER*2> bit_k_next;
  for(int i = 0; i < 4; i++){
    bit_k_next = ((bit_k << 2) |= i);
    if(mp_k.find(bit_k_next) != mp_k.end()){
      Next.push_back(bit_k_next);
    }
  }
  return Next;
}

//intput:bitset<K_MER*2>,map, output:int ;kmer入れるとそのkmerに入るedgeの数を出力
inline int before_kmer(const bitset<K_MER*2> bit_k,const unordered_map<bitset<K_MER*2>,int> mp_k){
  vector<bitset<K_MER*2> > Before;
  bitset<K_MER*2> bit_k_before;
  int edge_num = 0;
  bit_k_before = (bit_k >> 2);
  if(mp_k.find(bit_k_before) != mp_k.end()){
    ++edge_num;
  }
  bit_k_before.set(K_MER*2-2);
  if(mp_k.find(bit_k_before) != mp_k.end()){
    ++edge_num;
  }
  bit_k_before.set(K_MER*2-1);
  if(mp_k.find(bit_k_before) != mp_k.end()){
    ++edge_num;
  }
  bit_k_before.set(K_MER*2-2,false);
  if(mp_k.find(bit_k_before) != mp_k.end()){
    ++edge_num;
  }
  return edge_num;
}
//intput:bitset<K_MER*2>,map, output:int ;kmer入れるとそのkmerから出るedgeの数を出力
inline int after_kmer(const bitset<K_MER*2> bit_k,const unordered_map<bitset<K_MER*2>,int> mp_k){
  vector<bitset<K_MER*2> > After;
  bitset<K_MER*2> bit_k_after;
  int edge_num = 0;
  for(int i = 0; i < 4; i++){
    bit_k_after = ((bit_k << 2) |= i);
    if(mp_k.find(bit_k_after) != mp_k.end()){
      ++edge_num;
    }
  }
  return edge_num;
}

//input:vector<bit>, output:string  初めのbitにK_MER文字目を付け足していく
inline string vec_bit_to_str(const vector<bitset<K_MER*2> > v_bit){
  string str_out;
  bitset<K_MER*2> bit_tmp;
  str_out = bit_to_str(v_bit[0]);
  for(uint64_t i = 1; i < v_bit.size(); i++){
    switch(((v_bit[i] << (K_MER*2-2)) >> (K_MER*2-2)).to_ulong()){
    case 0:
      str_out += 'A';
      break;
    case 1:
      str_out += 'G';
      break;
    case 2:
      str_out += 'C';
      break;
    case 3:
      str_out += 'T';
      break;
    default:
      break;
    }
  }
  return str_out;
}

//input:vector<bit>,vector<bit>開始のbit, output:string  初めのbitにK_MER文字目を付け足していく
inline string vec_bit_to_str_mid(const vector<bitset<K_MER*2> > v_bit, const bitset<K_MER*2> bit_mid){
  string str_out;
  bool start_num=false;
  str_out = bit_to_str(bit_mid);
  for(uint64_t i = 1; i < v_bit.size(); i++){
    if(v_bit[i] == bit_mid){ //閉路開始点になったら
      start_num = true;
      continue;
    }
    if(start_num == true){ 
      switch(((v_bit[i] << (K_MER*2-2)) >> (K_MER*2-2)).to_ulong()){
      case 0:
	str_out += 'A';
	break;
      case 1:
	str_out += 'G';
	break;
      case 2:
	str_out += 'C';
	break;
      case 3:
	str_out += 'T';
	break;
      default:
	break;
      }
    }
  }
  return str_out;
}

//再帰関数による閉路探索
//input:kmer(bit),map(bit,int),vector<bit>
//output:int(通過node数)
inline void cycle_find(bitset<K_MER*2> &start_node, unordered_map<bitset<K_MER*2>,int > &mp_graph, unordered_map<bitset<K_MER*2>,bool > mp_visit, vector<bitset<K_MER*2> > Seq_bit, vector<string> &Repeat, int &depth, unsigned long long &node_num, bool &limit, const int &LENGTH_THR, const int &NODE_NUM_THR, const int &DEPTH_THR){
  ++depth; ++node_num;
  string str_answer; //出力用配列
  str_answer = vec_bit_to_str(Seq_bit);
  vector<bitset<K_MER*2> > Next;
  bitset<K_MER*2> bit_k_next;
  Next.clear();
  for(int i = 0; i < 4; ++i){
    bit_k_next = ((Seq_bit[Seq_bit.size()-1] << 2) |= i);
    if(mp_graph.find(bit_k_next) != mp_graph.end()){
      Next.push_back(bit_k_next);
    }
  }
  if(int(Next.size()) == 0 || int(Seq_bit.size()) > (LENGTH_THR - K_MER + 1)){ //長くなりすぎるorつなげなくなる:伸長終わり
    --depth;
    return;
  } else if(int(Next.size()) == 1){ //枝分かれなし1本
    while(int(Next.size()) == 1 && int(Seq_bit.size()) <= LENGTH_THR - K_MER + 1){ //長さが200以下で枝分かれなしの間
      if(Next[0] == start_node || mp_visit.find(Next[0]) != mp_visit.end()){ //loopを作っている
	if(Next[0] == start_node){ //startに戻っていれば
	  str_answer = vec_bit_to_str(Seq_bit);
	  //	  cout << str_answer.size() << "\t" << str_answer << "\t" << str_answer.substr(0,str_answer.size()-K_MER+1) << "\t" << bit_to_str(start_node) << "\t" << "直鎖"  << endl;
	  Repeat.push_back(str_answer);
	}
	--depth;
	return;
      } else { //startに戻らなければvisitつけて次のnodeへ(再帰)
	mp_visit.insert(make_pair(Next[0],true));
	Seq_bit.push_back(Next[0]); //bitを格納
	Next.clear();
	for(int i = 0; i < 4; ++i){
	  bit_k_next = ((Seq_bit[Seq_bit.size()-1] << 2) |= i);
	  if(mp_graph.find(bit_k_next) != mp_graph.end()){
	    Next.push_back(bit_k_next);
	  }
	}
	++node_num;
      }
    }
  }

  //Next[]をsort Next.size()の大きさで場合分け
  switch (int(Next.size())){
  case 2:
    if(mp_graph[Next[0]] < mp_graph[Next[1]]){
      bitset<K_MER*2> tmp1 = Next[0];
      Next[0] = Next[1];
      Next[1] = tmp1;
    }
    break;
  case 3:
    if(mp_graph[Next[0]] < mp_graph[Next[1]]){
      if(mp_graph[Next[1]] < mp_graph[Next[2]]){ //2,1,0の順 ->0と2を入れ替え
	bitset<K_MER*2> tmp1 = Next[0];
	Next[0] = Next[2];
	Next[2] = tmp1;
      } else if(mp_graph[Next[0]] < mp_graph[Next[2]]){ //1,2,0の順 -> 0,1,2を入れ替え
	bitset<K_MER*2> tmp1 = Next[0];
	Next[0] = Next[1];
	Next[1] = Next[2];
	Next[2] = tmp1;
      } else { //1,0,2の順 -> 0と1を入れ替え
	bitset<K_MER*2> tmp1 = Next[0];
	Next[0] = Next[1];
	Next[1] = tmp1;
      }
    } else if(mp_graph[Next[0]] < mp_graph[Next[2]]){ //2,0,1の順 ->0,1,2を入れ替え
	bitset<K_MER*2> tmp1 = Next[0];
	Next[0] = Next[2];
	Next[2] = Next[1];
	Next[1] = tmp1;
    } else if(mp_graph[Next[2]] > mp_graph[Next[1]]){ //0,2,1の順 ->1,2の入れ替え
      bitset<K_MER*2> tmp1 = Next[1];
      Next[1] = Next[2];
      Next[2] = tmp1;
    } //0,1,2の順だと何もしない
    break;
  case 4:
    int tmp_int = 0;
    bitset<K_MER*2> tmp1;
    for(int j = 1; j < 4; j++){ //Next[0]に出現頻度の高いものを格納
      if(mp_graph[Next[0]] < mp_graph[Next[j]]){ //Next[0]よりも大きいのがあればNext[0]に代入
	tmp_int = j;
      }
    }
    if(tmp_int != 0){ //0より大きいのがあれば交換
      tmp1 = Next[0];
      Next[0] = Next[tmp_int];
      Next[tmp_int] = Next[0];
      tmp_int = 0;
    }
    for(int j = 2; j < 4; j++){ //Next[0]に出現頻度の高いものを格納
      if(mp_graph[Next[1]] < mp_graph[Next[j]]){ //Next[1]よりも大きいのがあればNext[0]に代入
	tmp_int = j;
      }
      if(tmp_int != 0){ //1より大きいのがあれば交換
	tmp1 = Next[1];
	Next[1] = Next[tmp_int];
	Next[tmp_int] = Next[1];
	tmp_int = 0;
      }
    }
    if(mp_graph[Next[3]] > mp_graph[Next[2]]){
      tmp1 = Next[2];
      Next[2] = Next[3];
      Next[3] = tmp1;
    }
    break;
  }
  

  if(Next.size() > 1){ //枝分かれ2本以上
    int next_size;
    //    if(depth >= DEPTH_THR || node_num  > NODE_NUM_THR){ //depthが閾値を超えると
    if(depth >= DEPTH_THR){ //depthが閾値を超えると(nodeの閾値はなくす)
      for(uint64_t j = 1; j < Next.size(); j++){ //Next[0]に出現頻度の高いものを格納
	if(mp_graph[Next[0]] < mp_graph[Next[j]]){ //Next[0]よりも大きいのがあればNext[0]に代入
	    Next[0] = Next[j];
	}
      }
       next_size = 1; //Next[0]しか考えない
    } else { //閾値を超えてないなら全ての枝分かれを考える
      next_size = Next.size();
    }
    //////////i=0のときにloopすると、Seq_bitにpush_backされず、i=1以降でpop_backすると削ってしまう...
    for(int i = 0; i < next_size; ++i){
      if(Next[i] == start_node){ //loopをつくっている
	if(i != 0){ //初めの枝分かれを格納してるので、次からはそこを変更する必要あり、出力する場合は消す
	  Seq_bit.pop_back();
	}
	str_answer = vec_bit_to_str(Seq_bit);
	//	ofs << str_answer.size() << "\t" << str_answer  << "\t" << str_answer.substr(0,str_answer.size()-K_MER+1) << "\t" << bit_to_str(start_node) << "\t" << "枝分かれ" << i << endl;
	Repeat.push_back(str_answer);
	--depth;
	return;
      } else if(mp_visit.find(Next[i]) != mp_visit.end()){ //内部ループ
	--depth;
	return;
      } else { //startに戻らなければvisitつけて次のnodeへ(再帰)
	mp_visit.insert(make_pair(Next[i],true));
	if(i != 0){
	  Seq_bit[Seq_bit.size()-1] = Next[i];
	} else {
	  Seq_bit.push_back(Next[i]); //bitを格納
	}
	cycle_find(start_node,mp_graph,mp_visit,Seq_bit,Repeat,depth,node_num,limit,LENGTH_THR,NODE_NUM_THR,DEPTH_THR); //再帰呼び出し
      }
    }
  }
  --depth;
  return;
}


int Cycle::cycle_find_parent(){
  FILE *fp;
  char str[1024];
  string line,file_kmer;
  file_kmer = this->fa;
  string output = this->o;
  int data_type = this->d;
  int length_thr = this->max_l;
  int node_num_thr = this->max_n;
  int depth_thr = this->max_d;
  int num_threads = this->t;
  num_threads = 1; //hoge
  fp = fopen(file_kmer.c_str(),"r");
  if(fp == NULL || (data_type != 1 && data_type != 0)){
    if(fp == NULL){
      cout << "fp:" << file_kmer << endl;
    }
    if(data_type !=1 && data_type !=0){
      cout << data_type << endl;
    }
    cerr << "cant open in cycle_find.cpp" << endl;
    return -2;
  }
  long long number_t,number_s;
  int number_dif = 0;
  vector<string> a1,a2,a3;
  unordered_map<bitset<K_MER*2>,int> mp;
  pair<bitset<K_MER*2>,bitset<K_MER*2>> bit_pr;
  vector<bitset<K_MER*2> > start_bit;
  while((fgets(str,sizeof(str),fp)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      if(line[0] == '>' && data_type == 0){ //配列の出現回数を取得 TestisとSomaの差を考える
	a1 = split_str(line,"Testis");
	a2 = split_str(a1[1],"Somatic");
	a3 = split(a2[1],'_');
	number_t = stoll(a2[0]);
	number_s = stoll(a3[0]);
	number_dif = number_t - number_s;
      } else if(line[0] == '>' && data_type == 1) { //出現回数のみを考える
	a1 = split(line,'_');
	number_dif = stoi(a1[0].substr(1,line.size()-1));
      } else { //配列を取得、配列をkeyに出現回数を格納
	bit_pr = str_to_bit_com(line);
	mp[bit_pr.first] = number_dif;
	mp[bit_pr.second] = number_dif;
	start_bit.push_back(bit_pr.first); //閉路探索は相補鎖の片方のみで良いので格納
      }
  }

  //k-merの分類(edgeの数によるnodeの分類) 
  vector<bitset<K_MER*2> > Before,start_bit_junction;
  int in_num,out_num,one_one=0,counter_one=0;
  bitset<K_MER*2> bit_tmp1,bit_k_before,bit_k_after;
  vector<bitset<K_MER*2> > Next,Next_str,Seq_bit;
  int counter = 0,edge_num=0;
  unsigned long long node_num_=0;
  vector<string> Repeat_tmp,Repeat;

  for(auto it = start_bit.begin(); it != start_bit.end(); ++it){
     bit_k_before = (*it >> 2);
     edge_num = 0;
     if(mp.find(bit_k_before) != mp.end()){
       ++edge_num;
     }
     bit_k_before.set(K_MER*2-2);
     if(mp.find(bit_k_before) != mp.end()){
       ++edge_num;
     }
     bit_k_before.set(K_MER*2-1);
     if(mp.find(bit_k_before) != mp.end()){
       ++edge_num;
     }
     bit_k_before.set(K_MER*2-2,false);
     if(mp.find(bit_k_before) != mp.end()){
       ++edge_num;
     }
     in_num = edge_num;
     edge_num = 0;
     for(int i = 0; i < 4; i++){
       bit_k_after = ((*it << 2) |= i);
       if(mp.find(bit_k_after) != mp.end()){
   	++edge_num;
       }
     }
     out_num = edge_num;
     counter_one++;
     if(in_num * out_num > 1){ //candidate for start node
       Next.clear();
       start_bit_junction.push_back(*it);
       one_one++;
     }
   }

   fclose(fp);
   Next.reserve(4);
   unordered_map<bitset<K_MER*2>, bool> mp_visit,mp_visit_entire;
   cout << "\t" << "cycle finding parameter" << "\n" 
	<< "\t" << "length threshold:" << length_thr << "\n" 
	<< "\t" << "nodes threshold :" << node_num_thr << "\n" 
	<< "\t" << "depth threshold :" << depth_thr  << endl;
   int depth = 0;
   bool limit;
  omp_set_num_threads(num_threads);
#pragma omp parallel
  {
#pragma omp for firstprivate(Repeat_tmp,node_num_,depth,limit,counter,Seq_bit,mp_visit)
    for(uint64_t i = 0; i < start_bit_junction.size(); ++i){
#pragma omp critical (one)
          {
            Seq_bit.clear();
            mp_visit.clear();
            Seq_bit.push_back(start_bit_junction[i]);
            depth = 0;
            node_num_ = 0;
            limit = false;
            counter++;
          }
	  cycle_find(start_bit_junction[i],mp,mp_visit,Seq_bit,Repeat_tmp,depth,node_num_,limit,length_thr,node_num_thr,depth_thr);
	  Repeat.insert(Repeat.end(), Repeat_tmp.begin(), Repeat_tmp.end());
	  Repeat_tmp.clear();
#pragma omp critical (two)
        mp.erase(start_bit_junction[i]); //johnson's algorithmより探索したnodeは消去
    }
  }
  //  cout << Repeat.size() << endl; //hoge

  //枝分かれのない閉路探索を行う ここでの目的は枝分かれのない簡素な閉路を回収とパラメータによる閉路探索では得られなかった長いリピートの検出
  bool at_first,next_exist;
  bitset<K_MER*2> bit_k_next_max,start_node,bit_k_next;
  string str_answer,str_answer_;
  //  cout << "max length(no branch path):" << LENGTH_THR_STR  << endl;
  counter = 0;
  //  cout << start_bit.size() << endl;//hoge
  for(auto it = start_bit.begin(); it != start_bit.end(); ++it){ 
    start_node = *it;
    if(mp_visit_entire.find(start_node) == mp_visit_entire.end()){ //これまでにvisitしてないnodeのみを対象
      Seq_bit.clear();
      Seq_bit.push_back(start_node);
      mp_visit.clear();
      while(true){
	at_first = true;
	next_exist = false;
	for(int i = 0; i < 4; ++i){
	  bit_k_next = ((Seq_bit[Seq_bit.size()-1] << 2) |= i);
	  if(mp.find(bit_k_next) != mp.end()){
	    if(at_first == true){ //一つ目ならとりあえず格納
	      bit_k_next_max = bit_k_next; //最大のものに
	      at_first = false; //一つ目をfalseに
	      next_exist = true; //次があったのでtrueに
	    } else if(mp[bit_k_next_max] < mp[bit_k_next]){ //二つ目以降なら比較
	      bit_k_next_max = bit_k_next;
	    }
	  }
	}
	if(next_exist == false || Seq_bit.size() > LENGTH_THR_STR - K_MER + 1){ //case1:伸長終わり(次のkmerがない or 長くなりすぎる)
	  break;
	} else if(bit_k_next_max == start_node && next_exist == true){                            //case2:startに戻る
	  str_answer = vec_bit_to_str(Seq_bit);
	  Repeat.push_back(str_answer);
	  break;
	} else if(mp_visit.find(bit_k_next_max) != mp_visit.end() && next_exist == true){         //case3:start_nodeを含まずに閉路をつくる
	  str_answer = vec_bit_to_str_mid(Seq_bit,bit_k_next_max); //第二引数から最後までのstring
	  str_answer_ = vec_bit_to_str(Seq_bit);
	  //Repeat.push_back(str_answer);
	  break;
	} else if(next_exist == true){                                                            //case4:次のkmerがある
	  mp_visit.insert(make_pair(bit_k_next_max,true));
	  mp_visit_entire.insert(make_pair(bit_k_next_max,true));
	  mp_visit_entire.insert(make_pair(bit_to_bit_com(bit_k_next_max),true));
	  Seq_bit.push_back(bit_k_next_max); //bitを格納
	}
      } //while終わり
    }
  }

  //以下長さによるsort
  ofstream ofs(output);
  unordered_multimap<int,string> mp_sort; //長さをkey、配列をvalueに
  unordered_set<int> mp_set; //長さ格納
  for(auto it = Repeat.begin(); it != Repeat.end(); ++it){
    mp_sort.insert(make_pair((*it).size()-K_MER+1,(*it).substr(0,(*it).size() - K_MER + 1)));
    mp_set.insert((*it).size()-K_MER+1);
  }
  vector<int> v_sort;
  for(auto it = mp_set.begin(); it != mp_set.end(); ++it){
    v_sort.push_back(*it);
  }
  sort(v_sort.begin(),v_sort.end());
  counter = 0;
  for(auto it = v_sort.begin(); it != v_sort.end(); ++it){
    auto range = mp_sort.equal_range(*it);
    for(auto it_multi = range.first; it_multi != range.second; ++it_multi){
      ++counter;
      ofs << ">seq" << counter << "_len" << it_multi->first << "\n" << it_multi->second << endl;
    }
  }
  cout << "detected tandem repeats by cycle finding:" << Repeat.size() << endl;
  if(Repeat.size() == 0){
    return -1;
  } else {
    return 1;
  }
}

/*
void trf(string *file,string file_out, string file_out_tandem){
  string output = *file;
  char command_trf[16384];
  int trf_pr[7] = {2,7,7,80,10,50,500};
  snprintf(command_trf,sizeof(command_trf),"%s %d %d %d %d %d %d %d %s",("trf "+output).c_str(),trf_pr[0],trf_pr[1],trf_pr[2],trf_pr[3],trf_pr[4],trf_pr[5],trf_pr[6],(" -h -ngs |grep '@'|awk '{print substr($1,2);}' > " + output+"_tmp_trf").c_str());
  system(command_trf);
  FILE *fp1,*fp2;
  ifstream file_tmp,file0;
  fp1 = fopen((output+"_tmp_trf").c_str(),"r");
  fp2 = fopen(output.c_str(),"r");
  if(fp1 == NULL || fp2 == NULL){
    cerr << "can't open" << endl;
    return;
  }
  string line,line_out,seq_name,line_tandem;
  char str[16384];
  int length,length_pre = 0;
  uint64_t read_length = 100;
  unordered_set<string> set,set_exist;
  //  string output1 = output + "_remove_trf";
  //  string output2 = output + "_remove_trf_tandem";
  ofstream ofs(file_out);
  ofstream ofs_t(file_out_tandem);
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    set.insert(line);
  }
  fclose(fp1);

  remove((output+"_tmp_trf").c_str()); //trfの結果ファイルを削除

  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      seq_name = line.substr(1,line.size()-1);
    } else {
      length = line.size();
      if(length_pre == length){
	if(set_exist.find(line) == set_exist.end() && set_exist.find(complement(line)) == set_exist.end() && set.find(seq_name) == set.end()){
	  ofs << ">" << seq_name << "\n" << line << endl;
	  line_out = line + line;
	  while(line_out.size() < read_length){
	    line_out += line;
	  }
	  ofs_t << ">" << seq_name << "\n" << line_out << endl;
	  set_exist.clear();
	  for(int i = 0; i < length; ++i){
	    set_exist.insert(line.substr(length-1-i,1+i)+line.substr(0,length-1-i));
	  }
	}
      } else if(set.find(seq_name) == set.end()){
	ofs << ">" << seq_name << "\n" << line << endl;
	line_out = line + line;
	while(line_out.size() < read_length){
	  line_out += line;
	}
	ofs_t << ">" << seq_name << "\n" << line_out << endl;	
	set_exist.clear();
        for(int i = 0; i < length; ++i){
          set_exist.insert(line.substr(length-1-i,1+i)+line.substr(0,length-1-i));
        }
	length_pre = line.size();
      } else if(set.find(seq_name) != set.end() && length_pre != length){
	set_exist.clear();
	length_pre = line.size();
      }
    }
  }
  fclose(fp2);
}
*/
