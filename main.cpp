#include "common.h"
#include "bll.h"

void print_usage(void){
  cerr << "Usage: cycle_finder [Command] [option]" << endl;
  cerr << "------------------------------------------------" << endl;
  cerr << "Command                                         " << endl;
  cerr << "all         : Run the whole pipeline" << endl;
  cerr << "extract     : Extract k-mer" << endl;
  cerr << "cycle       : Find cycles from de Bruijn graph" << endl;
  cerr << "cluster     : Cluster repeats" << endl;
  cerr << "intersperse : Detect intersperse repeats" << endl;
  cerr << "------------------------------------------------" << endl;
}

void print_version(void){
  cout << "ROOT PATH: " << ROOT_PATH << endl;
  cout << "Version  : " << VERSION << endl;
  cout << "k-mer    : " << K_MER << endl;
}

int main(int argc, char **argv){
  print_version();  
  for(int i = 0; i < argc; ++i){
    string ss = argv[i];
    //------------------------- "all" -------------------------//
    if(ss == "all"){
      All all;
      all.all_exe(argc, argv);
    }
    //----------------------- "extract" -----------------------//
    if(ss == "extract"){
      Extract extract;
      if(extract.option_extract_parse(argc, argv) == -1){
	extract.print_extract_usage();
	return -1;
      }
      extract.extract_exe();
    }

    //----------------------- "cycle" -----------------------//
    if(ss == "cycle"){
      Cycle cycle;
      if(cycle.option_cycle_parse(argc, argv) == -1){
	cycle.print_cycle_usage();
	return -1;
      }
      cycle.cycle_exe();
    }
    //----------------------- "cluster" -----------------------//
    if(ss == "cluster"){
      Cluster cluster;
      if(cluster.option_cluster_parse(argc,argv) == -1){
	cluster.print_cluster_usage();
	return -1;
      }
      cluster.cluster_exe(); //"cluster"
    }

    //----------------------- "intersperse" -----------------------//
    if(ss == "intersperse"){
      Intersperse intersperse;
      if(intersperse.option_intersperse_parse(argc, argv) == -1){
	intersperse.print_intersperse_usage();
	return -1;
      }
      intersperse.intersperse_exe();
    }

    //just print help message
    if(argc < 2 || ss == "-h" ||  ss == "-help" || ss == "--help"){
      print_usage();
      return 0;
    }

    //when input other commands
    if((i == 1 && ss != "extract" && ss != "cycle" && ss != "cluster" && ss != "all" && ss != "intersperse")){
      print_usage();
      return -1;
    }
  }
  return 0;
}
