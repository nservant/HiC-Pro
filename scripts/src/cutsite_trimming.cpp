// HiC-Pro
// Copyright 2019 Institut Curie                               
// Author(s): Nicolas Servant
// Contact: nicolas.servant@curie.fr
// This software is distributed without any guarantee under the terms of the BSD-3 licence

// g++ -std=c++0x -o cutsite_trimming cutsite_trimming.cpp
//./cutsite_trimming -fastq fastq -cutsite AGCTT


#include <iostream>     // std::cout
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>

static const char* prog;

static int usage(int ret=1)
{
  std::cerr << "usage: " << prog << " --fastq FASTQFILE --cutsite CUTSITE --out OUTFILE [--rmuntrim] \n";
  std::cerr << "usage: " << prog << " --help\n";
  return ret;
}

/*
Get options from command line
*/
static int get_options(int argc, char* argv[], std::string& fastqFile,
                       std::vector<std::string>& cutSites, std::string& output, bool& rmuntrim)
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

        std::string cutSitesSequence;
        cutSitesSequence = std::string(argv[++ac]);
        size_t pos = cutSitesSequence.find(",");
        size_t begin = 0;
        while(pos != std::string::npos){
          cutSites.push_back(cutSitesSequence.substr(begin, pos - begin));
          begin = pos + 1;
          pos = cutSitesSequence.find(",", begin + 1);
        }
        cutSites.push_back(cutSitesSequence.substr(begin, pos));

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

/*
Replace 'N' in motifs
*/
static std::vector<std::string> check_cutSites(std::vector<std::string>& cutSites)
{
  std::vector<std::string> outCutSites;
  for(std::vector<std::string>::iterator it = cutSites.begin(); it != cutSites.end(); ++it){
    std::string motif = *it;
    size_t pos = motif.find('N');

    if( pos != std::string::npos){
      // replace N by A,C,G,T
      outCutSites.push_back(motif.replace(pos, 1, "A"));
      outCutSites.push_back(motif.replace(pos, 1, "C"));
      outCutSites.push_back(motif.replace(pos, 1, "G"));
      outCutSites.push_back(motif.replace(pos, 1, "T"));

      // Get the next occurrence from the current position
      pos = motif.find('N', pos);
      if( pos != std::string::npos){
	outCutSites = check_cutSites(outCutSites);
      }
    }else{
      outCutSites.push_back(motif);
    }
  }
  return (outCutSites);
}


static int trim_fastq(std::string& fastqFile,
                      std::vector<std::string>& cutSites,
                      std::string& outFile, bool& rmuntrim)
{

  int trim_count=0;
  std::string ID;
  std::ifstream ifs (fastqFile);
  std::ofstream ofs (outFile);

  if (ifs.is_open()){
    // foreach reads
    while (getline(ifs, ID)) {
      std::string seq;
      std::string dummy;
      std::string qual;
      
      getline(ifs, seq);
      getline(ifs, dummy);
      getline(ifs, qual);

      bool find_pos = false;
      size_t pos = std::string::npos;
      for (std::vector<std::string>::iterator it = cutSites.begin(); it != cutSites.end(); ++it){
	// return the first position detected
	size_t tmp_pos = seq.find(*it);
        if (tmp_pos != std::string::npos) {
          // If find_pos is alread True, there are two cut sites in the same read).
	  // Cut at the left-most position
          if (find_pos == true){
            if(tmp_pos < pos) {
              pos = tmp_pos;
            }
          } else {
            find_pos = true;
            pos = tmp_pos;
          }
        }
      }
      
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
      find_pos = false;
    }
  }else{
    std::cerr << "Error : Cannot open file : " << fastqFile;
  }
  return trim_count;
}

int main(int argc, char* argv[])
{
  
  std::string fastqFile;
  std::vector<std::string> cutSites;
  std::string outFile;
  bool rmuntrim = false;

  int ret = get_options(argc, argv, fastqFile, cutSites, outFile, rmuntrim);
  if (fastqFile.empty() || cutSites.size() == 0 || outFile.empty()){
    usage();
    exit(ret);
  }

  printf("##Fastq file: %s\n", fastqFile.c_str());

  cutSites = check_cutSites(cutSites);
  printf("##Restriction sites (n=%d):\n", (int)cutSites.size());
  for(std::vector<std::string>::iterator it = cutSites.begin(); it != cutSites.end(); ++it){
    std::cout << *it << std::endl;
  }
  printf("##Output File: %s\n", outFile.c_str());

  int trim_count=trim_fastq(fastqFile, cutSites, outFile, rmuntrim);
  printf("\n##Trimmed reads: %d\n", trim_count);
  return(0);
 }



