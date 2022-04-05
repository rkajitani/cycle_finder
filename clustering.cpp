#include "cluster.h"
#define TANDEM 2

struct Cluster_str {
  unordered_set<int> clst_mem; //構成するクラスターのfamily_id
  int max_mem; //最大コピー数をもつfamily_id
  int max_copy; //最大のコピー数
};

void Cluster::tandem(string repeat_num_file, string file_out){
  FILE *fp1;
  fp1 = fopen(repeat_num_file.c_str(),"r");
  if(fp1 == NULL){
    cerr << "error at tandem() in clustering.cpp; cant open " << repeat_num_file << endl;
  }
  char str[16384];
  string line,open_file1,tandem_str;
  int counter = 0;
  vector<string> a;
  ofstream ofs(file_out);
  while((fgets(str,sizeof(str),fp1)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      a = split(line,'\t');
      ++counter;
      ofs << ">" << counter << "_" << a[0].size() << endl;
      tandem_str = a[0];
      for(unsigned i = 0; i < K_MER/a[0].size(); ++i){
	tandem_str += a[0];
      }
      for(int i = 0; i < TANDEM; ++i){
	ofs << tandem_str;
      }
      ofs << "\n";
    }  
}

void Cluster::tandem_undo(string repeat_num_file, string file_out){
  FILE *fp1;
  fp1 = fopen(repeat_num_file.c_str(),"r");
  if(fp1 == NULL){
    printf("CANTOPEN\n");
  }
  char str[16384];
  string line,open_file1;
  vector<string> a;
  ofstream ofs(file_out);
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      ofs << line << endl;
    } else {
      ofs << line.substr(0,line.size()/TANDEM) << endl;
    }
  }
}

void Cluster::cdhit_integrage(string repeat_num_file, string cdhit1_file, string cdhit2_file, string file_out){
  FILE *fp1,*fp2,*fp3;
  char str[16384];
  fp1 = fopen(repeat_num_file.c_str(),"r");
  fp2 = fopen(cdhit1_file.c_str(),"r");
  fp3 = fopen(cdhit2_file.c_str(),"r");
  if(fp1 == NULL || fp2 == NULL || fp3 == NULL){
    printf("CANTOPEN\n");
  }
  string line,open_file1,open_file2,open_file3,open_file4;
  int number,consensus_id,cluster_num;
  vector<string> a;
  string tandem_str;

  vector<int> v_id; //クラスター内リピートの各idを一時的に格納し、mp_clstrにまとめて格納
  unordered_multimap<int,int> mp_clstr; //代表配列のid(行番号)をkeyにそのクラスターのidをvalueとして格納
  map<int,int> mp_copy;
  map<int,string> mp_str;
  int line_number = 0;
  ofstream ofs(file_out);
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    ++line_number;
    mp_copy[line_number] = stoi(a[3]);
    mp_str[line_number] = a[0];
  }
  fclose(fp1);
  //clusteringファイル1を開く
  fgets(str,sizeof(str),fp2);
  str[strlen(str)-1] = '\0';
  line = str;
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      for(auto it = v_id.begin(); it != v_id.end(); it++){
	if(*it != consensus_id){
	  mp_clstr.insert(make_pair(consensus_id,*it));
	}
      }
      v_id.clear();
    } else {
      number = stoi(line.substr(line.find(">")+1,line.find("_")-line.find(">")-1)); //何行目にあるやつかを格納
      v_id.push_back(number);
      if(line.find("*") != string::npos){ //*があればコンセンサス配列に
	consensus_id = number;
      }
      cluster_num++;
    }
  }
  for(auto it = v_id.begin(); it != v_id.end(); it++){
    if(*it != consensus_id){
      mp_clstr.insert(make_pair(consensus_id,*it));
    }
  }
  fclose(fp2);

  //clusteringファイル2を開く    
  v_id.clear();
  int copy_dif,copy_dif_pre = 0,column2_pre;
  string str0,str_tandem;
  unordered_set<string> Str_40mer;
  string file_out1 = file_out + "_1";
  ofstream ofs1(file_out1); //tandemにしないものを出力
  fgets(str,sizeof(str),fp3); //先に1行読む
  str[strlen(str)-1] = '\0';
  line = str;
  a = split(line,' ');
  column2_pre = stoi(a[1]); //cluster番号
  while((fgets(str,sizeof(str),fp3)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,' ');
    if(line[0] == '>'){
      str0 = mp_str[consensus_id];
      str_tandem.clear();
      for(unsigned i = 0; i < K_MER/str0.size()+2; ++i){ //tandemにつなげる
	str_tandem += str0;
      }
      ofs << ">family_" << column2_pre << "_" << mp_str[consensus_id].size() << "_"<< consensus_id << "*" << "\n" << str_tandem << endl; //代表配列のみfamily_の形で出力
      ofs1 << ">family_" << column2_pre << "_" << mp_str[consensus_id].size() << "_"<< consensus_id << "*" << "\n" << str0 << endl; //代表配列のみfamily_の形で出力

      for(auto it = v_id.begin(); it != v_id.end(); ++it){
	if(*it != consensus_id){
	  str0 = mp_str[*it];
	  str_tandem.clear();
	  for(unsigned i = 0; i < K_MER/str0.size()+2; ++i){ //tandemにつなげる
	    str_tandem += str0;
	  }
	  ofs << ">family_" << column2_pre << "_" << mp_str[*it].size() << "_"<< *it << "\n" << str_tandem << endl; 
	  ofs1 << ">family_" << column2_pre << "_" << mp_str[*it].size() << "_"<< *it << "\n" << str0 << endl;
	}
      }
      column2_pre = stoi(a[1]);
      v_id.clear();
      copy_dif_pre = 0;
    } else {
      number = stoi(line.substr(line.find(">")+1,line.find("_")-line.find(">")-1)); //何行目にあるやつかを格納
      copy_dif = mp_copy[number];
      if(copy_dif_pre < copy_dif){ //より大きいコピー差のものがきたら
	consensus_id = number;
	copy_dif_pre = copy_dif;
      }
      v_id.push_back(number); //各代表配列
      auto range = mp_clstr.equal_range(number);
      for(auto it = range.first; it != range.second; ++it){
	v_id.push_back(it->second); //クラスター全体
	copy_dif = mp_copy[it->second];
	if(copy_dif_pre < copy_dif){ //より大きいコピー差のものがきたら
	  consensus_id = it->second;
	  copy_dif_pre = copy_dif;
	}
      }
    }
  }
  str0 = mp_str[consensus_id];
  str_tandem.clear();
  for(unsigned i = 0; i < K_MER/str0.size()+2; ++i){ //tandemにつなげる
    str_tandem += str0;
  }
  ofs << ">family_" << column2_pre << "_" << mp_str[consensus_id].size() << "_"<< consensus_id << "*" << "\n" << str_tandem << endl; //代表配列のみfamily_の形で出力
  ofs1 << ">family_" << column2_pre << "_" << mp_str[consensus_id].size() << "_"<< consensus_id << "*" << "\n" << str0 << endl; //代表配列のみfamily_の形で出力
  for(auto it = v_id.begin(); it != v_id.end(); it++){
    if(*it != consensus_id){
      str0 = mp_str[*it];
      str_tandem.clear();
      for(unsigned i = 0; i < K_MER/str0.size()+2; ++i){ //tandemにつなげる
	str_tandem += str0;
      }
      ofs << ">family_" << column2_pre << "_" << mp_str[*it].size() << "_"<< *it  << "\n" << str_tandem << endl; //代表配列のみfamily_の形で出力
      ofs1 << ">family_" << column2_pre << "_" << mp_str[*it].size() << "_"<< *it  << "\n" << str0 << endl; //代表配列のみfamily_の形で出力
    }
  }
  fclose(fp3);
}

void Cluster::max_copy_fasta(string repeat_num_file, string cdhit_file, string file_out){
  FILE *fp1,*fp2;
  char str[16384];
  fp1 = fopen(repeat_num_file.c_str(),"r");
  fp2 = fopen(cdhit_file.c_str(),"r");
  if(fp1 == NULL || fp2 == NULL){
    printf("CANTOPEN\n");
  }
  string line,open_file1,open_file2,open_file3,open_file4;
  int number,consensus_id,cluster_num;
  vector<string> a;
  string tandem_str;

  vector<int> v_id; //クラスター内リピートの各idを一時的に格納し、mp_clstrにまとめて格納
  unordered_multimap<int,int> mp_clstr; //代表配列のid(行番号)をkeyにそのクラスターのidをvalueとして格納
  int line_number = 0;
  ofstream ofs(file_out);
  map<int,int> mp_copy;
  map<int,string> mp_str;
  unordered_multimap<int,int> mp_out;
  unordered_set<int> set_out;
  vector<int> v_out;
  int copy_dif_pre = 0,copy_dif = 0,cluster_num1=0;
  cluster_num = 0;
  string file_out1 = file_out + ".clstr";
  ofstream ofs1(file_out1);
  while((fgets(str,sizeof(str),fp1)) != NULL){ //resultファイルからcopy数,配列を格納
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    ++line_number;
    mp_copy[line_number] = stoi(a[3]); //copy数を格納
    mp_str[line_number] = a[0]; //配列を格納
  }
  fclose(fp1);

  fgets(str,sizeof(str),fp2); //1行読んでおく
  str[strlen(str)-1] = '\0';
  line = str;
  while((fgets(str,sizeof(str),fp2)) != NULL){ //.clstrファイルのクラスターの中で最もcopy数多いものを選ぶ
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>' && consensus_id != 0){
      mp_out.insert(make_pair(mp_str[consensus_id].size(),consensus_id));
      set_out.insert(mp_str[consensus_id].size());
      for(auto it = v_id.begin(); it != v_id.end(); ++it){ //consensus_idをkeyにクラスター内のidを格納
	mp_clstr.insert(make_pair(consensus_id,*it));
      }
      v_id.clear();
      copy_dif_pre = 0;
    } else if(line[0] != '>'){
      number = stoi(line.substr(line.find(">")+1,line.find("_")-line.find(">")-1)); //何行目にあるやつかを格納      
      v_id.push_back(number);
      copy_dif = mp_copy[number];
      if(copy_dif_pre < copy_dif){ //より大きいコピー差のものがきたら
	consensus_id = number;
	copy_dif_pre = copy_dif;
      }
    }
  }
  mp_out.insert(make_pair(mp_str[consensus_id].size(),consensus_id));
  set_out.insert(mp_str[consensus_id].size());
  for(auto it = v_id.begin(); it != v_id.end(); ++it){ //consensus_idをkeyにクラスター内のidを格納
    mp_clstr.insert(make_pair(consensus_id,*it));
  }
  v_id.clear();
    
  for(auto it = set_out.begin(); it != set_out.end(); ++it){ //長さ(重複無し)でsortするためにvectorに格納
    v_out.push_back(*it);
  }

  sort(v_out.begin(), v_out.end()); //sort
  string str_tandem;
  for(auto it = v_out.begin(); it != v_out.end(); ++it){ //長さ順に出力
    auto range = mp_out.equal_range(*it);
    for(auto it1 = range.first; it1 != range.second; ++it1){ //長さからconsensus_idを取り出し、consensus配列を出力
      str_tandem = mp_str[it1->second];
      for(unsigned i = 0; i < K_MER/mp_str[it1->second].size(); ++i){
	str_tandem += mp_str[it1->second];
      }
      ofs << ">" << it1->second << "_" << *it << "\n" << str_tandem << endl;
      auto range1 = mp_clstr.equal_range(it1->second); 
      for(auto it2 = range1.first; it2 != range1.second; ++it2){ //consensus_idからクラスター内のidを取り出す
	v_id.push_back(it2->second);
      }
      ofs1 << ">Cluster " << cluster_num << endl;
      for(auto it3 = v_id.begin(); it3 != v_id.end(); ++it3){ //全てのid
	if(*it3 != it1->second){
	  ofs1 << cluster_num1 << "\t" << mp_str[*it3].size() << "nt, >" << *it3 << "_" << mp_str[*it3].size() << "... at" << endl;
	} else {
	  ofs1 << cluster_num1 << "\t" << mp_str[*it3].size() << "nt, >" << *it3 << "_" << mp_str[*it3].size() << "... *" << endl;	    
	}
	++cluster_num1; //Clusterの中での番号
      }
      v_id.clear();
      ++cluster_num; //Clusterの番号
      cluster_num1 = 0;
    }
  }
  fclose(fp2);
}


void Cluster::cluster_by_blast(string repeat_num_file, string family_file, string blast_self_file, string max_copy_fasta, string file_out,int repeat_type){
  FILE *fp1,*fp2,*fp3,*fp4;
  fp1 = fopen(repeat_num_file.c_str(),"r");
  fp2 = fopen(family_file.c_str(),"r");
  fp3 = fopen(blast_self_file.c_str(),"r");
  fp4 = fopen(max_copy_fasta.c_str(),"r");    
  char str[16384];
  string line,open_file1,open_file2,open_file3,open_file4;

  if(fp1 == NULL || fp2 == NULL || fp3 == NULL || fp4 == NULL){
    printf("CANTOPEN\n");
  }
  vector<string> a;
  string tandem_str;

  vector<int> v_id; //クラスター内リピートの各idを一時的に格納し、mp_clstrにまとめて格納
  unordered_multimap<int,int> mp_clstr; //代表配列のid(行番号)をkeyにそのクラスターのidをvalueとして格納
  map<int,int> mp_copy;
  map<int,string> mp_str;
  int line_number = 0;
  ofstream ofs(file_out);
  string file_out1 = file_out + ".fa";
  ofstream ofs1(file_out1);
  string file_out2 = file_out + ".element";
  ofstream ofs2(file_out2);
  vector<string> b;
  if(repeat_type == 0){ //散在性ではコピー数ではなく塩基数を参考に代表配列を決める
    while((fgets(str,sizeof(str),fp1)) != NULL){ //resultファイルからcopy数,配列を格納
      str[strlen(str)-1] = '\0';
      line = str;
      a = split(line,'\t');
      ++line_number;
      //mp_copy[line_number] = stoi(a[3]); //copy数を格納
      //2017/12/10変更してみる(コピー数じゃなく塩基数を参考に代表配列を決める)(散在性でシミュレーションしているときに都合が悪かったので)
      mp_copy[line_number] = stoi(a[7]);
      mp_str[line_number] = a[0]; //配列を格納
    }
  } else if(repeat_type == 1){ //タンデムでは塩基数ではなくコピー数を参考に代表配列を決める
    while((fgets(str,sizeof(str),fp1)) != NULL){ //resultファイルからcopy数,配列を格納
      str[strlen(str)-1] = '\0';
      line = str;
      a = split(line,'\t');
      ++line_number;
      mp_copy[line_number] = stoi(a[3]); //copy数を格納
      //2017/12/10変更してみる(コピー数じゃなく塩基数を参考に代表配列を決める)(散在性でシミュレーションしているときに都合が悪かったので)
      mp_str[line_number] = a[0]; //配列を格納
    }
  }
  unordered_map<int,uint64_t> mp_copy_tes,mp_copy_soma;
  pair<int,int> tes_soma;
  unordered_map<int,string> mp_family;
  while((fgets(str,sizeof(str),fp2)) != NULL){ //familyファイルからfamily_idをkeyに総copy数(tes,soma)を格納
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    b = split(a[0],'_');
    mp_copy_tes.insert(make_pair(stoi(b[1]),stoll(a[2]))); //family_idをkey, (tes)をvalue
    mp_copy_soma.insert(make_pair(stoi(b[1]),stoll(a[3]))); //family_idをkey, (soma)をvalue
    mp_family[stoi(b[1])] = line; //family_idをkeyにline全体をvalue(blastnでhitしてないもの対象)
  }
  unordered_map<int,int> mp_fa_id;
  unordered_map<int,string> mp_fa_str;
  while((fgets(str,sizeof(str),fp4)) != NULL){ //fastaファイルからfamily_idをkeyにconsensus_id,配列を格納
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>' && line.find("*") != string::npos){
      a = split(line,'_');
      mp_fa_id[stoi(a[1])] = stoi(a[3].substr(0,a[3].size()-1));
      fgets(str,sizeof(str),fp4);
      str[strlen(str)-1] = '\0';
      line = str;
      mp_fa_str[stoi(a[1])] = line;
    }
  }
  vector<Cluster_str> Cluster_v; //形成するクラスター
  vector<string> c; 
  fgets(str,sizeof(str),fp3); //1行読む
  str[strlen(str)-1] = '\0';
  line = str;
  a = split(line,'\t');
  b = split(a[1],'_');
  c = split(a[2],'_');
  int clm2 = stoi(b[1]);
  int clm3 = stoi(c[1]);
  uint64_t copy_add_tes = mp_copy_tes[clm2] + mp_copy_tes[clm3];
  uint64_t copy_add_soma = mp_copy_soma[clm2] + mp_copy_soma[clm3];

  int Family_num = 0,copy_dif = 0;
  int copy_dif1 = mp_copy[stoi(b[3].substr(0,b[3].size()-1))]; //column2のコピー数差(resultファイルより)
  int copy_dif2 = mp_copy[stoi(c[3].substr(0,c[3].size()-1))]; //column3のコピー数差
  int max_id;
  v_id.push_back(clm2);
  v_id.push_back(clm3);
  unordered_set<int> used; //これまでに出てきたfamily_idを格納
  used.insert(clm2);
  used.insert(clm3);
  if(copy_dif1 > copy_dif2){
    max_id = clm2;
    copy_dif = copy_dif1;
  } else {
    max_id = clm3;
    copy_dif = copy_dif2;
  }
  Cluster_str tmp;
  unordered_map<int,int> mp_id2cl; //family_idをkeyに所属しているクラスターのCluster_vでの添字をvalue
  tmp.clst_mem.insert(clm2);
  tmp.clst_mem.insert(clm3); 
  tmp.max_mem = max_id;
  tmp.max_copy = copy_dif;
  Cluster_v.push_back(tmp);
  mp_id2cl.insert(make_pair(clm2,Cluster_v.size()-1));
  mp_id2cl.insert(make_pair(clm3,Cluster_v.size()-1));
  while((fgets(str,sizeof(str),fp3)) != NULL){ //blastファイルからクラスタリング(column1にfamily_idでsortしたもの、column2,3にfamily名)
    str[strlen(str)-1] = '\0';
    line = str;
    a = split(line,'\t');
    b = split(a[1],'_');
    c = split(a[2],'_');
    clm2 = stoi(b[1]);
    clm3 = stoi(c[1]);
    copy_dif1 = mp_copy[stoi(b[3].substr(0,b[3].size()-1))]; //column2のコピー数差(resultファイルより)
    copy_dif2 = mp_copy[stoi(c[3].substr(0,c[3].size()-1))]; //column3のコピー数差
    //2,3カラムが既出かどうかの2*2通り
    if(used.find(clm2) != used.end() && used.find(clm3) == used.end()){ //2カラムが既出、3カラムが新 
      //3カラムをクラスターに格納
      Cluster_v[mp_id2cl[clm2]].clst_mem.insert(clm3);
      used.insert(clm3);
      if(Cluster_v[mp_id2cl[clm2]].max_copy < copy_dif2){ //3カラムの方が大きければ
	Cluster_v[mp_id2cl[clm2]].max_copy = copy_dif2;
	Cluster_v[mp_id2cl[clm2]].max_mem = clm3;
      }
      mp_id2cl.erase(clm3);
      mp_id2cl.insert(make_pair(clm3,mp_id2cl[clm2]));
    } else if(used.find(clm2) != used.end() && used.find(clm3) != used.end() && mp_id2cl[clm2] != mp_id2cl[clm3]){ //両方既出で異なるクラスター
      //2カラムと3カラムのクラスターを統合(3カラムのクラスターを2カラムのクラスターに全て格納)
      if(Cluster_v[mp_id2cl[clm2]].max_copy < Cluster_v[mp_id2cl[clm3]].max_copy){ //3カラムの方が大きければ
	Cluster_v[mp_id2cl[clm2]].max_mem = Cluster_v[mp_id2cl[clm3]].max_mem;
	Cluster_v[mp_id2cl[clm2]].max_copy = Cluster_v[mp_id2cl[clm3]].max_copy;
      }
      unordered_set<int> tmp_set = Cluster_v[mp_id2cl[clm3]].clst_mem; 
      for(auto it = tmp_set.begin(); it != tmp_set.end(); ++it){ //3カラムのクラスターを全走
	Cluster_v[mp_id2cl[clm2]].clst_mem.insert(*it); //2カラムのクラスターに格納
	mp_id2cl.erase(*it);
	mp_id2cl.insert(make_pair(*it,mp_id2cl[clm2])); //3カラムのクラスター要素から2カラムのクラスターの添字にリンク
      }
    } else if(used.find(clm2) == used.end() && used.find(clm3) == used.end()){ //両方、新
      //2カラムから新しいクラスターをつくる
      if(copy_dif1 > copy_dif2){
	max_id = clm2;
	copy_dif = copy_dif1;
      } else {
	max_id = clm3;
	copy_dif = copy_dif2;
      }
      Cluster_str tmp1;
      tmp1.clst_mem.insert(clm2);
      tmp1.clst_mem.insert(clm3); 
      tmp1.max_mem = max_id;
      tmp1.max_copy = copy_dif;
      Cluster_v.push_back(tmp1);
      mp_id2cl.insert(make_pair(clm2,Cluster_v.size()-1));
      mp_id2cl.insert(make_pair(clm3,Cluster_v.size()-1));
      used.insert(clm2);
      used.insert(clm3);
    } else if(used.find(clm2) == used.end() && used.find(clm3) != used.end()){ //2カラムが新、3カラムが既出
      //2カラムをクラスターに格納
      Cluster_v[mp_id2cl[clm3]].clst_mem.insert(clm2);
      used.insert(clm2);
      if(Cluster_v[mp_id2cl[clm3]].max_copy < copy_dif1){ //2カラムの方がコピー数大きければ
	Cluster_v[mp_id2cl[clm3]].max_copy = copy_dif1;
	Cluster_v[mp_id2cl[clm3]].max_mem = clm2;
      }
      mp_id2cl.erase(clm2);
      mp_id2cl.insert(make_pair(clm2,mp_id2cl[clm3]));
    }
  }

  //出力
  unordered_set<int> tmp_set1;
  for(auto it = mp_id2cl.begin(); it != mp_id2cl.end(); ++it){
    if(tmp_set1.find(it->second) == tmp_set1.end()){
      tmp_set1.insert(it->second);
      tmp = Cluster_v[it->second];
      ++Family_num;
      copy_add_tes = 0,copy_add_soma = 0;
      ofs2 << ">Family_" << Family_num << "_" << mp_fa_str[tmp.max_mem].size() << "_" << mp_fa_id[tmp.max_mem] << endl;
      for(auto it1 = tmp.clst_mem.begin(); it1 != tmp.clst_mem.end(); ++it1){ //クラスターの要素を全走
	copy_add_tes += mp_copy_tes[*it1];
	copy_add_soma += mp_copy_soma[*it1];
	ofs2 << *it1 << endl;
      }
      ofs << "Family_" << Family_num  << "_" << mp_fa_str[tmp.max_mem].size() << "_" << mp_fa_id[tmp.max_mem] << "\t" << copy_add_tes - copy_add_soma << "\t" << copy_add_tes << "\t" << copy_add_soma << "\t" << double(copy_add_tes - copy_add_soma)/mp_fa_str[tmp.max_mem].size() << endl;
      ofs1 << ">Family_" << Family_num << "_" << mp_fa_str[tmp.max_mem].size() << "_" << mp_fa_id[tmp.max_mem] << "\n" << mp_fa_str[tmp.max_mem] << endl;
    }
  }
  //usedに格納されたもの全てをmp_copy_tesから削除
  for(auto it = used.begin(); it != used.end(); ++it){
    mp_copy_tes.erase(*it);
  }
  for(auto it = mp_copy_tes.begin(); it != mp_copy_tes.end(); ++it){ //これまでの手順で出力されなかったものを出力
    ++Family_num;
    line = mp_family[it->first];
    a = split(line,'\t');
    b = split(a[0],'_');
    ofs << "Family_" << Family_num << "_" << mp_fa_str[stoi(b[1])].size() << "_" << mp_fa_id[stoi(b[1])] << "\t" << a[1] << "\t" << a[2] << "\t" << a[3] << "\t" << stod(a[1])/mp_fa_str[stoi(b[1])].size() << endl;
    ofs1 << ">Family_" << Family_num << "_" << mp_fa_str[stoi(b[1])].size() << "_" << mp_fa_id[stoi(b[1])] << "\n" << mp_fa_str[stoi(b[1])] << endl;
  }
  fclose(fp1);
}
