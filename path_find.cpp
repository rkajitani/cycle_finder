#include "intersperse.h"

#define MIN_LEN 0

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


//pairのsecondでsortできるように
template <class T1, class T2, class Pred = std::less<T2> >
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
    Pred p;
    return p(left.second, right.second);
  }
};

//再帰関数によるパス探索
//input:kmer(bit),map(bit,int),vector<bit>
//output:int(通過node数)
inline void path_find(bitset<K_MER*2> &start_node, unordered_map<bitset<K_MER*2>,int > &mp_graph, unordered_map<bitset<K_MER*2>,bool > mp_visit, vector<bitset<K_MER*2> > Seq_bit, vector<string> &Repeat, int &depth, long long &node_num, bool &limit, const int &NODE_NUM_THR, const int &DEPTH_THR){
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
  ++node_num;
  if(Next.size() == 0){ //つなげなくなる:伸長終わり
    if(Seq_bit.size() > MIN_LEN){
      str_answer = vec_bit_to_str(Seq_bit);
      Repeat.push_back(str_answer);
    }
    --depth;
    return;
  } else if(Next.size() == 1){ //枝分かれなし1本
    while(Next.size() == 1){ //長さが200以下で枝分かれなしの間
      if(Next[0] == start_node || mp_visit.find(Next[0]) != mp_visit.end()){ //loopを作っている
        --depth;
        return;
      } else { //startに戻らなければvisitつけて次のnodeへ(再帰)
        mp_visit[Next[0]] = true; //visit
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
  if(Next.size() > 1){ //枝分かれ2本以上
    int next_size;
    if(depth >= DEPTH_THR || node_num  > NODE_NUM_THR){ //depthが閾値を超えると
      if(limit == false){ //まだ文を出力してなければ
        limit = true;
      }
      for(unsigned j = 1; j < Next.size(); j++){ //Next[0]に出現頻度の高いものを格納
        if(mp_graph[Next[0]] < mp_graph[Next[j]]){ //Next[0]よりも大きいのがあればNext[0]に代入
	  Next[0] = Next[j];
        }
      }
      next_size = 1; //Next[0]しか考えない
    } else { //閾値を超えてないなら全ての枝分かれを考える
      next_size = Next.size();
    }
    for(int i = 0; i < next_size; ++i){
      if(Next[i] == start_node || mp_visit.find(Next[i]) != mp_visit.end()){ //loopをつくってる
	--depth;
	return;
      } else { //startに戻らなければvisitつけて次のnodeへ(再帰)
        mp_visit[Next[i]] = true; //visit
        if(i == 0){
          Seq_bit.push_back(Next[i]); //bitを格納
        } else {
          Seq_bit[Seq_bit.size()-1] = Next[i];
        }
        path_find(start_node,mp_graph,mp_visit,Seq_bit,Repeat,depth,node_num,limit,NODE_NUM_THR,DEPTH_THR); //再帰呼び出し
      }
    }
  }
  --depth;
  return;
}

int Intersperse::path_find_parent(){
  FILE *fp;
  char str[1024];
  string line,file_kmer;

  int number = 0;
  string file_kmer_i = this->kmer_for_path_find;
  string output = this->o;
  int data_type = this->d;
  int length_thr = this->max_l;
  int node_num_thr = this->max_n;
  int depth_thr = this->max_d;
  int num_threads = this->t;
  int length_min_thr = this->length_min_thr;
  cout << "\t" << "path finding parameter" << endl;
  cout << "\t" << "length threshold:" << length_thr << endl;
  cout << "\t" << "node threshold  :" << node_num_thr << endl;
  cout << "\t" << "depth threshold :" << depth_thr << endl;
  fp = fopen(file_kmer_i.c_str(),"r");
  if(fp == NULL || (data_type != 1 && data_type != 0)){
    printf("CANTOPEN\n");
    return -1;
  }
  vector<string> a;
  unordered_map<bitset<K_MER*2>,int> mp;
  pair<bitset<K_MER*2>,bitset<K_MER*2>> bit_pr;
  vector<bitset<K_MER*2> > start_bit;
  int counter = 0;
  char str1[1024],str2[1024];
  if(data_type == 1){
    while((fgets(str,256,fp)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      if(str[0] == '>'){
        a = split(line,'\t');
        number = stoi(a[0].substr(1,line.size()-1));
      } else {
        bit_pr = str_to_bit_com(line);
        mp.insert(make_pair(bit_pr.first,number));
        mp.insert(make_pair(bit_pr.second,number));
        start_bit.push_back(bit_pr.first); //閉路探索は相補鎖の片方のみで良いので格納
      }
    }
  } else if(data_type == 0){
    while((fgets(str,256,fp)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      if(str[0] == '>'){
        while(str[counter] != 'S'){
          ++counter;
        }
        for(int i = 0; i < counter; i++){
          str1[i] = str[i+7];
        }
        counter += 7;
        for(unsigned i = 0; i < strlen(str)-counter; i++){
          str2[i] = str[i+counter];
        }
        number = atoi(str1) - atoi(str2);
        counter = 0;
        memset(str1,0,sizeof(str1));
        memset(str2,0,sizeof(str2));
      } else {
        bit_pr = str_to_bit_com(line);
        mp.insert(make_pair(bit_pr.first,number));
        mp.insert(make_pair(bit_pr.second,number));
        start_bit.push_back(bit_pr.first); //閉路探索は相補鎖の片方のみで良いので格納
      }
    }
  }
  //  return 0; //hoge

  //k-merの分類(edgeの数によるnodeの分類) hoge
  vector<bitset<K_MER*2> > Before,start_bit_junction;
  int in_num,one_one=0,counter_one=0;
  vector<bitset<K_MER*2> > Next,Next_str,Seq_bit;
  bitset<K_MER*2> bit_tmp1,bit_k_before,bit_k_after;
  vector<string> Repeat_tmp;
  int edge_num=0;
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
    counter_one++;
    if(in_num ==  0){ //start nodeとなりうるnode
      Next.clear();
      start_bit_junction.push_back(*it);
      one_one++;
    }
  }
  fclose(fp);
  ofstream ofs(output);
  Next.reserve(4);
  unordered_map<bitset<K_MER*2>, bool> mp_visit,mp_visit_entire;
  int depth = 0;
  long long node_num = 0;
  vector<string> Repeat;
  bool limit;
  counter = 0;
  omp_set_num_threads(num_threads);
#pragma omp parallel
  {
#pragma omp for firstprivate(Repeat_tmp,node_num,depth,limit,counter,Seq_bit,mp_visit)
    for(unsigned i = 0; i < start_bit_junction.size(); ++i){
#pragma omp critical (one)
      {
        Seq_bit.clear();
        mp_visit.clear();
        //    Seq_bit.push_back(*it);
        Seq_bit.push_back(start_bit_junction[i]);
        depth = 0;
        node_num = 0;
        limit = false;
        ++counter;
      }
      path_find(start_bit_junction[i],mp,mp_visit,Seq_bit,Repeat_tmp,depth,node_num,limit,node_num_thr,depth_thr);
#pragma omp critical (two)
      Repeat.insert(Repeat.end(), Repeat_tmp.begin(), Repeat_tmp.end());
      Repeat_tmp.clear();
      mp.erase(start_bit_junction[i]);
    }
  }

  //以下長さによるsort
  unordered_multimap<int,string> mp_sort; //長さをkey、配列をvalueに
  unordered_set<int> mp_set; //長さ格納
  for(auto it = Repeat.begin(); it != Repeat.end(); ++it){
    mp_sort.insert(make_pair((*it).size(),*it));
    mp_set.insert((*it).size());
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
      if(it_multi->first >= length_min_thr){
        ++counter;
        ofs << ">seq" << counter << "_len" << it_multi->first << "\n" << it_multi->second << endl;
      }
    }
  }
  cout << "\t" << "detected interspersed repeats :" << counter << endl;
  if(counter == 0){
    return -1;
  } else {
    return 0;
  }
}
