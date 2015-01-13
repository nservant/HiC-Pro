// g++ -std=c++0x -o cutsite_trimming cutsite_trimming.cpp
//./cutsite_trimming -fastq fastq -cutsite AGCTT


#include <iostream>     // std::cout
#include <stdlib.h>
#include <string.h>
#include <fstream>

static const char* prog;

static int usage(int ret=1)
{
  std::cerr << "usage: " << prog << " --fastq FASTQFILE --cutsite CUTSITE --out OUTFILE [--rmuntrim] \n";
  std::cerr << "usage: " << prog << " --help\n";
  return ret;
}

static int get_options(int argc, char* argv[], std::string& fastqFile, std::string& cutSite, std::string& output, bool& rmuntrim)
{
  prog = argv[0];
  if (argc == 1){
    exit(usage());
  }
  for (int ac = 1; ac < argc; ++ac) {
    const char* opt = argv[ac];
    if (*opt == '-') {
      if (!strcmp(opt, "--fastq")) {
	fastqFile = std::string(argv[++ac]);
      } else if (!strcmp(opt, "--cutsite")) {
	cutSite = std::string(argv[++ac]);
      } 
      else if (!strcmp(opt, "--out")) {
        output = std::string(argv[++ac]);
      }
      else if (!strcmp(opt, "--rmuntrim")) {
        rmuntrim = true;
      }
    }else {
      std::cerr << prog << ": unknown option " << opt << std::endl;
      return usage();
    } 
  }
  return 0;
}

static int trim_fastq(std::string& fastqFile, std::string& cutSite, std::string& outFile, bool& rmuntrim)
{

  int trim_count=0;
  std::string ID;
  std::ifstream ifs (fastqFile);
  std::ofstream ofs (outFile);

  if (ifs.is_open()){
    while (getline(ifs, ID)) {
      std::string seq;
      std::string dummy;
      std::string qual;
      
      getline(ifs, seq);
      getline(ifs, dummy);
      getline(ifs, qual);
      size_t pos = seq.find(cutSite);
      
      if (pos != std::string::npos) {
	trim_count++;
	ofs << ID << '\n';
	ofs << seq.substr(0, pos) << '\n';
	ofs << "+\n";
	ofs << qual.substr(0, pos) << '\n';
      } else {
	if (!rmuntrim){
	  ofs << ID << '\n';
       	  ofs << seq << '\n';
	  ofs << "+\n";
	  ofs << qual << '\n';
	}
      }
    }
  }else{
    std::cerr << "Error : Cannot open file : " << fastqFile;
  }
  return trim_count;
}

int main(int argc, char* argv[])
{
  
  std::string fastqFile;
  std::string cutSite;
  std::string outFile;
  bool rmuntrim = false;

  int ret = get_options(argc, argv, fastqFile, cutSite, outFile, rmuntrim);
  printf("##Fastq file: %s\n", fastqFile.c_str());
  printf("##Restriction site: %s\n", cutSite.c_str());
  printf("##Output File: %s\n", outFile.c_str());

  if (fastqFile.empty() || cutSite.empty() || outFile.empty()){
    usage();
    exit(ret);
  }

  int trim_count=trim_fastq(fastqFile, cutSite, outFile, rmuntrim);
  printf("\n##Trimmed reads: %d\n", trim_count);
  return(0);
 }



