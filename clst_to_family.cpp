//#include "clst_to_family.h"
#include "cluster.h"
//cd-hit-estの結果をfamilyの名前にして、すべての配列をfastaで出力

void Cluster::clst_to_family(string open_file1, string open_file2, string output){
  FILE *fp1,*fp2;
  char str[16384]; //1行16384文字までしか読み込めない
  string line;
  fp1 = fopen(open_file1.c_str(),"r");
  fp2 = fopen(open_file2.c_str(),"r");
  if(fp1 == NULL || fp2 == NULL){
    if(fp1 == NULL) cerr << "error at clst_to_family() in clst_to_family.cpp; cant open " << open_file1 << endl;
    if(fp2 == NULL) cerr << "error at clst_to_family() in clst_to_family.cpp; cant open " << open_file2 << endl;
    return;
  }

  //fastaファイルを格納
  unordered_map<string,string> mp_fa;
  string seq_name;
  while((fgets(str,sizeof(str),fp1)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      if(line[0] == '>'){
	seq_name = line;
      } else {
	mp_fa.insert(make_pair(seq_name,line));
      }
  }
  fclose(fp1);

  //cd-hit outputファイル
  int family_num;
  vector<string> a,Seq_name;
  string seq_represent;
  ofstream ofs(output);
  //1行読む
  fgets(str,sizeof(str),fp2);
  str[strlen(str)-1] = '\0';
  line = str;
  a = split(line,' ');
  family_num = stoi(a[1]);
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      ofs << ">family_" << family_num << "_" << seq_represent.substr(seq_represent.find("_")+1) << "_" << seq_represent.substr(1,seq_represent.find("_")-1) << "*" << "\n" << mp_fa[seq_represent] << endl; //代表配列
      for(auto it = Seq_name.begin(); it != Seq_name.end(); ++it){
	ofs << ">family_" << family_num << "_" << (*it).substr((*it).find("_")+1) << "_" << (*it).substr(1,(*it).find("_")-1) << "\n" << mp_fa[*it] << endl;
      }
      a = split(line,' ');
      family_num = stoi(a[1]);
      Seq_name.clear();
    } else {
      a = split_str(line.substr(line.find(">")),"...");
      if(line.find("*") != string::npos){ //代表配列だけ別
	seq_represent = a[0]; 
      } else { //代表配列以外はvectorに格納
	Seq_name.push_back(a[0]);
      }
    }
  }
  //最終行
      ofs << ">family_" << family_num << "_" << seq_represent.substr(seq_represent.find("_")+1) << "_" << seq_represent.substr(1,seq_represent.find("_")-1) << "*" << "\n" << mp_fa[seq_represent] << endl; //代表配列
  for(auto it = Seq_name.begin(); it != Seq_name.end(); ++it){
    ofs << ">family_" << family_num << "_" << (*it).substr((*it).find("_")+1) << "_" << (*it).substr(1,(*it).find("_")-1) << "\n" << mp_fa[*it] << endl;
  }

  fclose(fp2);
  
}
