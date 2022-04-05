#include "extract.h"

inline uint64_t peak_detect_exe(string file){
  FILE *fp1;
  char str[16384]; //1行16384文字までしか読み込めない
  string line,output;
  fp1 = fopen(file.c_str(),"r");
  vector<string> a;
  vector<long long> v;
  long long num1,num2,num3,def1,def2,max=0,max_x=0;
  if(fp1 == NULL){
    printf("CANTOPEN\n");
    return -1;
  }

  //まず2行よむ
  fgets(str,sizeof(str),fp1);
  str[strlen(str)-1] = '\0';
  line = str;
  a = split(line,' ');
  v.push_back(stol(a[1]));
  num1 = stol(a[1]);
  long long min = stol(a[1]);
  //2行目
  fgets(str,sizeof(str),fp1);
  str[strlen(str)-1] = '\0';
  line = str;
  a = split(line,' ');
  v.push_back(stol(a[1]));
  num2 = stol(a[1]);
  def1 = num2 - num1;
  while((fgets(str,sizeof(str),fp1)) != NULL){
      str[strlen(str)-1] = '\0';
      line = str;
      a = split(line,' ');
      num3 = stol(a[1]);
      v.push_back(num3);
      def2 = num3 - num2;
      if(def1*def2 > 0){ //単調
	
      } else if(def1 >= 0 && def2 <= 0 && max < num2){ //極大値候補
	max = num2;
	max_x = stol(a[0])-1;
      } else if(def1 <= 0 && def2 >= 0 && min > num2){ //極小値候補
	min = num2;
      }
      def1 = def2;
      num1 = num2;
      num2 = num3;
  }
  return max_x;
  fclose(fp1);
}

void Extract::write_peak_to_file(){
  ofstream ofs(this->kmer_peak_file);
  if(this->d == 1){
    ofs << "peak:" << this->peak << endl;
    cerr << "peak:" << this->peak << endl;    
  } else if(this->d == 0){
    ofs << "peak1:" << this->peak1 << endl;
    ofs << "peak2:" << this->peak2 << endl;
    cerr << "peak1:" << this->peak1 << endl;
    cerr << "peak2:" << this->peak2 << endl;
  }
}

void Extract::peak_detect(){
  string file, file1, file2;
  if(this->d == 1){
    file = this->histo_file;
    this->peak = peak_detect_exe(file);
  } else if(this->d == 0){
    file1 = this->histo_file1;
    file2 = this->histo_file2;
    this->peak1 = peak_detect_exe(file1);
    this->peak2 = peak_detect_exe(file2);
  } else {
    cout << "error" << endl;
  }
  write_peak_to_file();
}
