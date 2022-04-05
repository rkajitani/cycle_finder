#include "cycle.h"

void Cycle::trf_filter(){
  string output = this->o;
  string file_out = this->trf_out;
  string file_out_tandem = this->trf_out_tandem;
  char command_trf[16384];
  int trf_pr[7] = {2,7,7,80,10,50,500};  
  string root_path = ROOT_PATH;
  if(file_exist(output) != 0){
    snprintf(command_trf,sizeof(command_trf),"%s %d %d %d %d %d %d %d %s",(root_path + "/trf "+output).c_str() ,trf_pr[0],trf_pr[1],trf_pr[2],trf_pr[3],trf_pr[4],trf_pr[5],trf_pr[6],("-h -ngs |grep '@'|awk \'{print substr($1,2);}' > " + output+"_tmp_trf").c_str());
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
}
