#include "cycle.h"

struct REGION{
  int start,end,len;
  string family;
};

bool comp(const REGION& left, const REGION& right){ //AREAの比較関数:startでsort
  return left.start < right.start;
}

inline REGION line_to_struct(string &line){ //lineをタブ区切りしてそれぞれALIの構造体にいれていく
  vector<string> a,b;
  REGION reg;
  a = split(line,'\t');
  b = split_str(a[0],"len");
  reg.family = a[0];
  reg.len = stoi(b[1]);
  reg.start = stoi(a[6]);
  reg.end = stoi(a[7]);
  return reg;
}


int Cycle::blastn2fa(){
  string file_blastn = this->blast_read_out;
  string file_fa = this->o;
  int repeat_type = this->repeat_type;
  string file_out = this->fa_filter;
  FILE *fp1,*fp2,*fp3;
  char str[16384]; //1行16384文字までしか読み込めない
  string line,open_file1,open_file2;
  fp1 = fopen(file_blastn.c_str(),"r");
  fp2 = fopen(file_fa.c_str(),"r");  
  fp3 = fopen(file_out.c_str(),"w");  
  if(fp1 == NULL ||fp2 == NULL){
    printf("CANTOPEN\n");
    return -1;
  }
  vector<string> Column,Column1;
  string family_name,family_name_pre,seq_name;
  REGION region;
  vector<REGION> Region;
  unordered_map<string,vector<REGION> > mp; //family_nameをkeyにvector<REGION>を呼び出す
  unordered_set<string> Family_set;

  //fastaファイル読み込む
  unordered_map<string,string> mp_fa;
  while((fgets(str,sizeof(str),fp2)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    if(line[0] == '>'){
      seq_name = line.substr(1,line.size()-1);
    } else {
      mp_fa.insert(make_pair(seq_name,line));
    }
  }

  //1行読む
  fgets(str,sizeof(str),fp1);
  str[strlen(str)-1] = '\0';
  line = str;
  region = line_to_struct(line);
  family_name = region.family;
  Family_set.insert(family_name);
  family_name_pre = family_name;
  Region.push_back(region);
  //~1行読む
  while((fgets(str,sizeof(str),fp1)) != NULL){
    str[strlen(str)-1] = '\0';
    line = str;
    region = line_to_struct(line);
    family_name = region.family;
    if(family_name == family_name_pre){ //上の行と同じfamily(repeat)なら
      Region.push_back(region);
    } else {
      Family_set.insert(family_name);
      sort(Region.begin(),Region.end(),comp); //startでsort
      mp.insert(make_pair(family_name_pre,Region)); //familyをkeyにfamilyにアライメントされた位置の構造体のvectorをvalueに
      Region.clear();
    }
    family_name_pre = family_name;
  }
  sort(Region.begin(),Region.end(),comp); //startでsort
  mp.insert(make_pair(family_name_pre,Region)); //familyをkeyにfamilyにアライメントされた位置の構造体のvectorをvalueに
  Region.clear();
  fclose(fp1);

  int range_min,range_max;
  vector<string> seq_name_len;
  string seq_str; //出力する配列
  //tandemの場合
  if(repeat_type == 1){
    for(auto it = Family_set.begin(); it != Family_set.end(); ++it){
      Region = mp[*it];
      range_min = Region[0].start;
      range_max = Region[0].end;
      if(range_min == 1){
    for(unsigned i = 0; i < Region.size(); ++i){
      if(Region[i].start > range_max){ //一部の領域を飛ばすようにアライメントされていればその配列はミスアセンブルとみなす
        break;
      } else if(Region[i].end > range_max){ //領域のmaxの更新
        range_max = Region[i].end;
      }
    }
    if(range_max >= Region[0].len*2){
      fprintf(fp3,">%s\n%s\n",(*it).c_str(),mp_fa[*it].c_str()); //%d:number, %s:文字列
    }
      }
    }
  }
  //tandemでない場合
  else if(repeat_type == 0){
    for(auto it = Family_set.begin(); it != Family_set.end(); ++it){
      Region = mp[*it];
      range_min = Region[0].start;
      range_max = Region[0].end;
      if(range_min <= ceil(Region[0].len * 0.05)){ //長さの0.05倍からのstartなら認める
		for(unsigned i = 0; i < Region.size(); ++i){
		  if(Region[i].start > range_max){ //一部の領域を飛ばすようにアライメントされていればその配列はミスアセンブルとみなす
			break;
		  } else if(Region[i].end > range_max){ //領域のmaxの更新
			range_max = Region[i].end;
		  }
		}
		if(range_max >= floor(Region[0].len * 0.95)){ //長さの0.95倍以上なら出力
		  if(range_min != 0 || range_max < Region[0].len){ //完全長のアライメントがとれてなければ配列名の長さを変更
			seq_name_len = split_str(*it,"len");
			seq_name = seq_name_len[0] + to_string(static_cast<long long>(range_max - range_min + 1));
			seq_str = mp_fa[*it].substr(range_min, range_max - range_min + 1);
		  } else {
			seq_name = *it;
			seq_str = mp_fa[*it];
		  }
		  fprintf(fp3,">%s\n%s\n",seq_name.c_str(),seq_str.c_str()); //%d:number, %s:文字列
		}
      }
    }
  }
  fclose(fp2);
  fclose(fp3);

  if(file_exist(file_out) == 0){ //出力したファイルがあるか
    cerr << repeat_type << " ERROR : COULDN'T MAP READS; no " << file_out << endl;
    return -1;
  } else {
    return 0;
  }
}
