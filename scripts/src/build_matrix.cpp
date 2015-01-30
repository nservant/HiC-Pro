// g++ -std=c++0x -o build_matrix build_matrix.cpp
//./build_matrix --binsize 10000 --chrsizes ../annotation/chrom.sizes --ifile ../hic_results/data/mAST-Rad21-WT-rep1/SRR941305.mm9.interaction --oprefix /tmp/zoo

// ./build_matrix --binsize 10000 --chrsizes ../annotation/chrom.sizes --ifile ../../Hi-C_pipelintest_ev_all/hic_results/data/mAST-Rad21-WT-rep1/SRR941305.mm9.interaction --oprefix /tmp/zoo --chrA chr1:chr2 --chrB chr1:chr2:chr8:chr10


#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

static const int SPARSE_FMT = 0x1;
static const int BED_FMT = 0x2;
static const char* prog;

typedef unsigned int chrsize_t;
class Chromosome {

private:
  static std::unordered_map<std::string, Chromosome*> chr_map;
  void computeSizes(chrsize_t ori_binsize, chrsize_t step, bool binadjust);

  std::string name;
  chrsize_t chrsize;
  chrsize_t binsize;
  chrsize_t stepsize;
  chrsize_t bincount;

public:
  Chromosome(const std::string& name, chrsize_t chrsize, chrsize_t ori_binsize, chrsize_t step, bool binadjust) : name(name), chrsize(chrsize) {
    computeSizes(ori_binsize, step, binadjust);
    assert(chr_map.find(name) == chr_map.end());
    chr_map[name] = this;
  }

  void adjustBinsize(chrsize_t ori_binsize, const chrsize_t step);

  const std::string& getName() const {return name;}
  chrsize_t getChrsize() const {return chrsize;}
  chrsize_t getBinsize() const {return binsize;}
  chrsize_t getStepsize() const {return stepsize;}
  chrsize_t getBincount() const {return bincount;}

  static chrsize_t getCount() {
    return chr_map.size();
  }

  static Chromosome* getByName(const std::string& name) {
    return chr_map[name];
  }
};

class AxisChromosome {
  int idx; // really needed ?
  const Chromosome* chr;
  chrsize_t binstart;
  chrsize_t binend;

public:
  AxisChromosome(int binoffset, const Chromosome* chr, const AxisChromosome* lastAxisChr) : chr(chr) {
    if (lastAxisChr != NULL) {
      binstart = lastAxisChr->getBinend();
    } else {
      binstart = binoffset;
    }
    binend = binstart + chr->getBincount();
    std::cerr << "AxisChromosome: " << chr->getName() << " " << binstart << " " << binend << " " << chr->getBincount() << std::endl;
  }

  chrsize_t getBinstart() const {return binstart;}
  chrsize_t getBinend() const {return binend;}

  chrsize_t getChrsize() const {return chr->getChrsize();}
  chrsize_t getBinsize() const {return chr->getBinsize();}
  chrsize_t getStepsize() const {return chr->getStepsize();}
  chrsize_t getBincount() const {return chr->getBincount();}

  const Chromosome* getChromosome() const {return chr;}
};

class Matrix {

  std::vector<AxisChromosome*> axis_chr_abs;
  std::vector<AxisChromosome*> axis_chr_ord;
  std::unordered_map<std::string, AxisChromosome*> axis_chr_abs_map;
  std::unordered_map<std::string, AxisChromosome*> axis_chr_ord_map;

  std::map<chrsize_t, std::map<chrsize_t, chrsize_t> > mat;
  void addAxisChromosome(const std::vector<const Chromosome*>& chr_v, std::vector<AxisChromosome*>& axis_chr, std::unordered_map<std::string, AxisChromosome*>& axis_chr_map);

  const AxisChromosome* getAxisChromosome(const std::string& chrname, const std::unordered_map<std::string, AxisChromosome*>& axis_chr_map) const {
    std::unordered_map<std::string, AxisChromosome*>::const_iterator iter = axis_chr_map.find(chrname);
    if (iter == axis_chr_map.end()) {
      return NULL;
    }
    return (*iter).second;
  }

  void displayBed(std::ostream& ofs, const std::vector<AxisChromosome*>& axis_chr) const {
    std::vector<AxisChromosome*>::const_iterator begin = axis_chr.begin();
    std::vector<AxisChromosome*>::const_iterator end = axis_chr.end();
    while (begin != end) {
      const AxisChromosome* axis_chr = *begin;
      const std::string& name = axis_chr->getChromosome()->getName();
      chrsize_t binstart = axis_chr->getBinstart();
      chrsize_t binend = axis_chr->getBinend();
      chrsize_t binsize = axis_chr->getBinsize();
      chrsize_t chrsize = axis_chr->getChrsize();
      binend -= binstart;
      for (chrsize_t bin = 0; bin < binend; ++bin) {
	// bed are 0-based begin, 1-based end
	chrsize_t beg = bin * binsize;
	chrsize_t end = beg + binsize - 1;
	if (end > chrsize) {
	  end = chrsize-1;
	}
	ofs << name << '\t' << beg << '\t' << (end+1) << '\t' << (bin+binstart) << '\n';
      }
      ++begin;
    }
  }

  int binoffset;

public:
  Matrix(int binoffset) : binoffset(binoffset) {}

  void addXAxisChromosome(const std::vector<const Chromosome*>& chr_v);
  void addYAxisChromosome(const std::vector<const Chromosome*>& chr_v);

  const AxisChromosome* getXAxisChromosome(const std::string& chrname) const {
    return getAxisChromosome(chrname, axis_chr_abs_map);
  }

  const AxisChromosome* getYAxisChromosome(const std::string& chrname) const {
    return getAxisChromosome(chrname, axis_chr_ord_map);
  }

  void add(chrsize_t abs_bin, chrsize_t ord_bin) {
    std::map<chrsize_t, std::map<chrsize_t, chrsize_t> >::iterator iter = mat.find(abs_bin);
    if (iter == mat.end()) {
      mat[abs_bin] = std::map<chrsize_t, chrsize_t>();
      mat[abs_bin][ord_bin] = 1;
    } else {
      (*iter).second[ord_bin]++;
    }
  }

  void displayMatrix(std::ostream& ofs) const {
    std::map<chrsize_t, std::map<chrsize_t, chrsize_t> >::const_iterator begin = mat.begin();
    std::map<chrsize_t, std::map<chrsize_t, chrsize_t> >::const_iterator end = mat.end();
    while (begin != end) {
      chrsize_t abs = (*begin).first;
      const std::map<chrsize_t, chrsize_t>& line = (*begin).second;
      std::map<chrsize_t, chrsize_t>::const_iterator bb = line.begin();
      std::map<chrsize_t, chrsize_t>::const_iterator ee = line.end();
      while (bb != ee) {
	ofs << abs << '\t' << (*bb).first << '\t' << (*bb).second << '\n';
	++bb;
      }
      ++begin;
    }
  }

  void displayXBed(std::ostream& ofs) const {
    displayBed(ofs, axis_chr_abs);
  }

  void displayYBed(std::ostream& ofs) const {
    displayBed(ofs, axis_chr_ord);
  }
};

void Matrix::addAxisChromosome(const std::vector<const Chromosome*>& chr_v, std::vector<AxisChromosome*>& axis_chr, std::unordered_map<std::string, AxisChromosome*>& axis_chr_map)
{
  std::vector<const Chromosome*>::const_iterator begin = chr_v.begin();
  std::vector<const Chromosome*>::const_iterator end = chr_v.end();

  const AxisChromosome* lastAxisChr = NULL;
  while (begin != end) {
    const Chromosome* chr = *begin;
    AxisChromosome* axisChr = new AxisChromosome(binoffset, chr, lastAxisChr);
    axis_chr.push_back(axisChr);
    axis_chr_map[chr->getName()] = axisChr;
    lastAxisChr = axisChr;
    ++begin;
  }
}

void Matrix::addXAxisChromosome(const std::vector<const Chromosome*>& chr_v)
{
  addAxisChromosome(chr_v, axis_chr_abs, axis_chr_abs_map);
}

void Matrix::addYAxisChromosome(const std::vector<const Chromosome*>& chr_v)
{
  addAxisChromosome(chr_v, axis_chr_ord, axis_chr_ord_map);
}

std::unordered_map<std::string, Chromosome*> Chromosome::chr_map;

enum Format {
  SPARSE_IND_FMT = SPARSE_FMT,
  SPARSE_BED_FMT = SPARSE_FMT|BED_FMT,
  EXPANDED_FMT = 0x4
};

void Chromosome::adjustBinsize(chrsize_t ori_binsize, const chrsize_t step)
{
  bincount = 1 + (chrsize_t)floor( (double)(chrsize-ori_binsize) / (ori_binsize/step));
  binsize = chrsize / bincount;
  stepsize = binsize / step;
}

void Chromosome::computeSizes(chrsize_t ori_binsize, chrsize_t step, bool binadjust)
{
  if (chrsize < ori_binsize) {
    binsize = chrsize;
    stepsize = chrsize;
    bincount = 1;
  } else if (binadjust) {
    adjustBinsize(ori_binsize, step);
  } else {
    binsize = ori_binsize;
    stepsize = (chrsize_t)floor(ori_binsize/step);
    chrsize_t remainder = (chrsize - ori_binsize) % stepsize;
    chrsize_t tmp_bincount = 1 + (chrsize_t)floor(chrsize-ori_binsize)/stepsize;
    bincount = remainder > 0 ? tmp_bincount+1 : tmp_bincount;
  }
  std::cerr << name << " sizes: " << chrsize << " " << binsize << " " << stepsize << " " << bincount << std::endl;
}

static int usage(int ret = 1)
{
  std::cerr << "usage: " << prog << " --binsize BINSIZE --chrsizes FILE --ifile FILE\n";
  //  std::cerr << "       --oprefix PREFIX [--fmt sparse_bed|sparse_ind|expanded] [--bed-prefix PREFIX]\n";
  std::cerr << "       --oprefix PREFIX [--binadjust] [--step STEP] [--binoffset OFFSET]\n";
  std::cerr << "       [--matrix-format asis|upper|lower|complete][--chrA CHR... --chrB CHR...]\n";
  std::cerr << "       [--legacy-data]\n";
  std::cerr << "\nusage: " << prog << " --help\n";
  return ret;
}

static int help()
{
  // TBD: complete help
  std::cerr << "\n";
  (void)usage();
  std::cerr << "\nOPTIONS\n\n";
  std::cerr << "  --binsize BINSIZE      : bin size\n";
  std::cerr << "  --chrsizes FILE        : file containing chromosome sizes\n";
  std::cerr << "  --ifile FILE           : input interaction file\n";
  std::cerr << "  --oprefix PREFIX       : output prefix of generated files (matrix and bed)\n";
  //  std::cerr << "  --fmt FORMAT             : \n";
  //  std::cerr << "  --bed-prefix PREFIX      :\n";
  std::cerr << "  --binadjust            : [optional] adjust bin sizes, default is false\n";
  std::cerr << "  --step STEP            : [optional] step size, default is 1\n";
  std::cerr << "  --binoffset OFFSET     : [optional] starting bin offset, default is 1\n";
  std::cerr << "  --matrix-format FORMAT : [optional] FORMAT may be:\n";
  std::cerr << "                           - asis: matrix is generated according to input data (default)\n";
  std::cerr << "                           - upper: only the upper matrix is generated\n";
  std::cerr << "                           - lower: only the lower matrix is generated\n";
  std::cerr << "                           - complete: generate both parts of the matrix (upper and lower);\n";
  std::cerr << "                             input data must contain only one part (upper or lower) \n";
  std::cerr << "  --chrA CHR             : [optional] colon separated list of abscissa chromosomes; default is all chromosomes\n";
  std::cerr << "  --chrB CHR             : [optional] colon separated list of ordinate chromosomes; default is all chromosomes\n";
  //  std::cerr << "  --organism ORGANISM      :\n";
  std::cerr << "  --legacy-data          : [optional] use for compatibility mode for old format \n";
  return -1;
}

enum MatrixFormat {
  ASIS_MATRIX = 1,
  UPPER_MATRIX,
  LOWER_MATRIX,
  COMPLETE_MATRIX
};
  
static int get_options(int argc, char* argv[], chrsize_t& binsize, const char*& chrsize_file, const char*& ifile, const char*& oprefix, Format& format, std::string& bed_prefix, bool& binadjust, MatrixFormat& matrix_format, chrsize_t& step, bool& whole_genome, bool& public_data, int& binoffset, const char*& chrA, const char*& chrB)
{
  prog = argv[0];
  for (int ac = 1; ac < argc; ++ac) {
    const char* opt = argv[ac];
    if (*opt == '-') {
      if (!strcmp(opt, "--binadjust")) {
	binadjust = true;
      } else if (!strcmp(opt, "--matrix-format")) {
	if (ac == argc-1) {
	  return usage();
	}
	std::string matrix_format_str = argv[++ac];
	if (matrix_format_str == "asis") {
	  matrix_format = ASIS_MATRIX;
	} else if (matrix_format_str == "upper") {
	  matrix_format = UPPER_MATRIX;
	} else if (matrix_format_str == "lower") {
	  matrix_format = LOWER_MATRIX;
	} else if (matrix_format_str == "complete") {
	  matrix_format = COMPLETE_MATRIX;
	} else {
	  return usage();
	}
      } else if (!strcmp(opt, "--step")) {
	if (ac == argc-1) {
	  return usage();
	}
	step = atoi(argv[++ac]);
      } else if (!strcmp(opt, "--binsize")) {
	if (ac == argc-1) {
	  return usage();
	}
	binsize = atoi(argv[++ac]);
      } else if (!strcmp(opt, "--public-data")) {
	public_data = true;
      } else if (!strcmp(opt, "--legacy-data")) {
	public_data = false;
      } else if (!strcmp(opt, "--binoffset")) {
	if (ac == argc-1) {
	  return usage();
	}
	binoffset = atoi(argv[++ac]);
      } else if (!strcmp(opt, "--ifile")) {
	if (ac == argc-1) {
	  return usage();
	}
	ifile = argv[++ac];
      } else if (!strcmp(opt, "--oprefix")) {
	if (ac == argc-1) {
	  return usage();
	}
	oprefix = argv[++ac];
      } else if (!strcmp(opt, "--chrsizes")) {
	if (ac == argc-1) {
	  return usage();
	}
	chrsize_file = argv[++ac];
	/*
      } else if (!strcmp(opt, "--fmt")) {
	if (ac == argc-1) {
	  return usage();
	}
	const char* fmt = argv[++ac];
	if (!strcasecmp(fmt, "sparse_bed")) {
	  format = SPARSE_BED_FMT;
	} else if (!strcasecmp(fmt, "sparse_ind")) {
	  format = SPARSE_IND_FMT;
	} else if (!strcasecmp(fmt, "expanded")) {
	  format = EXPANDED_FMT;
	} else {
	  return usage();
	}
      } else if (!strcmp(opt, "--bed-prefix")) {
	if (ac == argc-1) {
	  return usage();
	}
	bed_prefix = argv[++ac];
	*/
      } else if (!strcmp(opt, "--chrA")) {
	if (ac == argc-1) {
	  return usage();
	}
	chrA = argv[++ac];
	whole_genome = false;
      } else if (!strcmp(opt, "--chrB")) {
	if (ac == argc-1) {
	  return usage();
	}
	chrB = argv[++ac];
	whole_genome = false;
      } else if (!strcmp(opt, "--help")) {
	return help();
      } else {
	std::cerr << prog << ": unknown option " << opt << std::endl;
	return usage();
      }
    }
  }

  return 0;
}

static void split_in_vect(const std::string& str, std::vector<const Chromosome*>& vect)
{
  size_t last_pos = 0;
  while (size_t pos = str.find(':', last_pos)) {
    std::string chrname;
    bool last = pos == std::string::npos;
    if (last) {
      chrname = str.substr(last_pos);
    } else {
      chrname = str.substr(last_pos, pos-last_pos);
    }
    const Chromosome* chr = Chromosome::getByName(chrname);
    if (!chr) {
      std::cerr << prog << ": unknown chromosome " << chrname << std::endl;
      exit(1);
    }
    vect.push_back(chr);
    if (last) {
      break;
    }
    last_pos = pos+1;
  }
}

static int build_matrix_init(Matrix& matrix, const char* ifile, std::ifstream& ifs, const std::string& oprefix, std::ofstream& matfs, std::ofstream& xbedfs, std::ofstream& ybedfs, const char* chrsize_file, bool whole_genome, const char* chrA, const char* chrB, chrsize_t ori_binsize, chrsize_t step, bool binadjust)
{
  ifs.open(ifile);
  if (ifs.bad() || ifs.fail()) {
    std::cerr << prog << " cannot open interaction file: " << ifile << " for reading\n";
    return 1;
  }

  std::ifstream chrsizefs;
  chrsizefs.open(chrsize_file);
  if (chrsizefs.bad()) {
    std::cerr << prog << " cannot open chrsizes file: " << chrsize_file << " for reading\n";
    return 1;
  }

  std::string matfile = oprefix + ".matrix";
  matfs.open(matfile);
  if (matfs.bad() || matfs.fail()) {
    std::cerr << prog << " cannot open file: " << matfile << " for writing\n";
    return 1;
  }

  std::string xbedfile = oprefix + "_abs.bed";
  xbedfs.open(xbedfile);
  if (xbedfs.bad() || xbedfs.fail()) {
    std::cerr << prog << " cannot open file: " << xbedfile << " for writing\n";
    return 1;
  }

  std::string ybedfile = oprefix + "_ord.bed";
  if (whole_genome) {
    std::string xbedlink;
    size_t pos = xbedfile.rfind('/');
    if (pos != std::string::npos) {
      xbedlink = xbedfile.substr(pos+1);
    } else {
      xbedlink = xbedfile;
    }
    unlink(ybedfile.c_str());
    if (symlink(xbedlink.c_str(), ybedfile.c_str())) {
      std::cerr << prog << " cannot created link: " << ybedfile << "\n";
      return 1;
    }
  } else {
    ybedfs.open(ybedfile);
    if (ybedfs.bad() || ybedfs.fail()) {
      std::cerr << prog << " cannot open file: " << ybedfile << " for writing\n";
      return 1;
    }
  }


  std::vector<const Chromosome*> all_chr_v;
  while (!chrsizefs.eof()) {
    std::string buffer;
    getline(chrsizefs, buffer);

    chrsize_t chrsize;
    /*
    char name[256];
    if (sscanf(buffer.c_str(), "%s %u", name, &chrsize) == 2) {
      Chromosome* chromosome = new Chromosome(name, chrsize, ori_binsize, step, binadjust);
      all_chr_v.push_back(chromosome);
    }
    */
    std::istringstream istr(buffer);
    std::string name;
    istr >> name >> chrsize;
    if (!istr.fail()) {
      Chromosome* chromosome = new Chromosome(name, chrsize, ori_binsize, step, binadjust);
      all_chr_v.push_back(chromosome);
    }
  }

  chrsizefs.close();

  if (chrA) {
    assert(chrB != NULL);
    std::vector<const Chromosome*> chrA_v;
    std::vector<const Chromosome*> chrB_v;
    split_in_vect(chrA, chrA_v);
    split_in_vect(chrB, chrB_v);
    matrix.addXAxisChromosome(chrA_v);
    matrix.addYAxisChromosome(chrB_v);
  } else {
    matrix.addXAxisChromosome(all_chr_v);
    matrix.addYAxisChromosome(all_chr_v);
  }

  return 0;
}

static int interaction_parse2(char* buffer, char*& lchr, chrsize_t& lstart, char*& rchr, chrsize_t& rstart)
{
  char c;
  char* str;
  while ((c = *buffer++) != 0) {
    if (c == '\t') {
      lchr = buffer;
      break;
    }
  }
  while ((c = *buffer) != 0) {
    if (c == '\t') {
      *buffer++ = 0;
      str = buffer;
      break;
    }
    buffer++;
  }

  while ((c = *buffer) != 0) {
    if (c == '\t') {
      *buffer++ = 0;
      lstart = atoi(str);
      break;
    }
    buffer++;
  }

  while ((c = *buffer++) != 0) {
    if (c == '\t') {
      rchr = buffer;
      break;
    }
  }

  while ((c = *buffer) != 0) {
    if (c == '\t') {
      *buffer++ = 0;
      str = buffer;
      break;
    }
    buffer++;
  }

  while ((c = *buffer) != 0) {
    if (c == '\t') {
      *buffer++ = 0;
      rstart = atoi(str);
      break;
    }
    buffer++;
  }

  return 0;
}

static int interaction_parse(const std::string& str, std::string& mark, std::string& org, std::string& chr, chrsize_t& start, int &dist)
{
  size_t lastpos = 0;
  size_t pos = str.find('|', lastpos);
  if (pos == std::string::npos) {
    return 1;
  }
  //mark = str.substr(lastpos, pos);

  lastpos = pos+1;
  pos = str.find('|', lastpos);
  if (pos == std::string::npos) {
    return 1;
  }
  //org = str.substr(lastpos, pos-lastpos);

  lastpos = pos+1;
  pos = str.find(':', lastpos);
  if (pos == std::string::npos) {
    return 1;
  }
  chr = str.substr(lastpos, pos-lastpos);

  lastpos = pos+1;
  pos = str.find('-', lastpos);
  if (pos == std::string::npos) {
    return 1;
  }
  start = atoi(str.substr(lastpos, pos-lastpos).c_str());

  lastpos = pos+1;
  pos = str.find('@', lastpos);
  if (pos == std::string::npos) {
    return 1;
  }
  dist = atoi(str.substr(pos+1).c_str());
  return 0;
}

static void get_left_right(const char buffer[], char left[], char right[])
{
  const char* start = buffer;
  for (;;) {
    char c = *start;
    if (c != ' ' && c != '\t') {
      break;
    }
    start++;
  }
  const char* end = start+1;
  for (;;) {
    char c = *end;
    if (c == ' ' || c == '\t' || !c) {
      break;
    }
    end++;
  }
  size_t size = end-start;
  strncpy(left, start, size);
  left[size] = 0;
  start = end + 1;
  for (;;) {
    char c = *start;
    if (c != ' ' && c != '\t') {
      break;
    }
    start++;
  }
  end = start+1;
  for (;;) {
    char c = *end;
    if (c == ' ' || c == '\t' || !c) {
      break;
    }
    end++;
  }
  size = end-start;
  strncpy(right, start, size);
  right[size] = 0;
}

#if 0
static void adjust_binsize(const chrsize_t chrsize, const chrsize_t binsize, const chrsize_t step, chrsize_t& binsize_ret, chrsize_t& stepsize_ret, chrsize_t& bincount_ret)
{
  chrsize_t bincount = 1 + (chrsize_t)floor( (double)(chrsize-binsize) / (binsize/step));
  //chrsize_t bincount = (chrsize_t)floor( (double)(chrsize) / (binsize/step));
#if 1
  // EV code
  binsize_ret = chrsize / bincount;
  stepsize_ret = binsize_ret / step;
  bincount_ret = bincount;
#else
  // original code !
  chrsize_t stepsize = (chrsize_t)floor(binsize/step);
  chrsize_t cur_remainder = (chrsize-binsize) % stepsize;
  chrsize_t binsize_mod = binsize+floor((double)cur_remainder/bincount);
  chrsize_t stepsize_mod = floor((double)binsize_mod/step);
  if (binsize == binsize_mod){
    binsize_ret = binsize_mod;
    stepsize_ret = stepsize_mod;
    bincount_ret = bincount;
    return;
  }
  adjust_binsize(chrsize, binsize_mod, step, binsize_ret, stepsize_ret, bincount_ret);
#endif
}

static void get_sizes(const std::unordered_map<std::string, chrsize_t>& chrsize_map, std::unordered_map<std::string, chrsize_t>& binsize_map, std::unordered_map<std::string, chrsize_t>& stepsize_map, std::unordered_map<std::string, chrsize_t>& bincount_map, chrsize_t ori_binsize, chrsize_t step, bool binadjust)
{
  std::unordered_map<std::string, chrsize_t>::const_iterator begin = chrsize_map.begin();
  std::unordered_map<std::string, chrsize_t>::const_iterator end = chrsize_map.end();
  while (begin != end) {
    const std::string& chr = (*begin).first;
    chrsize_t chrsize = (*begin).second;
    if (chrsize < ori_binsize) {
      binsize_map[chr] = chrsize;
      stepsize_map[chr] = chrsize;
      bincount_map[chr] = 1;
    } else if (binadjust) {
      chrsize_t binsize, stepsize, bincount;
      adjust_binsize(chrsize, ori_binsize, step, binsize, stepsize, bincount);
      binsize_map[chr] = binsize;
      stepsize_map[chr] = stepsize;
      bincount_map[chr] = bincount;
    } else {
      binsize_map[chr] = ori_binsize;
      chrsize_t stepsize = (chrsize_t)floor(ori_binsize/step);
      stepsize_map[chr] = stepsize;
      chrsize_t remainder = (chrsize - ori_binsize) % stepsize;
      chrsize_t bincount = 1 + (chrsize_t)floor(chrsize-ori_binsize)/stepsize;
      bincount_map[chr] = remainder > 0 ? bincount+1 : bincount;
    }
    std::cout << chr << " sizes: " << chrsize << " " << binsize_map[chr] << " " << stepsize_map[chr] << " " << bincount_map[chr] << '\n';
    ++begin;
  }
  std::cout << std::endl;
}
#endif

/*
sub assign_bin{
    my ($chr,$st,$dist,$genome)=@_;
    my $loc=$st+$dist-1;
    my $binsize=$chrm_bin_def{$chr}{"binsize"};
    my $stepsize=$chrm_bin_def{$chr}{"stepsize"};
    my $cur_binidx=1+ceil(($loc-$binsize)/$stepsize);
    my $cur_binst=$stepsize*($cur_binidx-1)+1;
    my $cur_binend=$cur_binst+$binsize-1;
    $cur_binend=$chrm_size{$chr} if ($cur_binend>$chrm_size{$chr});
    my $cur_hit="HIC_" . $chr . "_" . $cur_binidx . "|" . $genome . "|" . $chr . ":" . $cur_binst . "-" . $cur_binend;
    return ($cur_hit,$cur_binidx);
}
*/

// starting port

static chrsize_t assign_bin(const std::string& org, const AxisChromosome* chr, chrsize_t start, int dist)
{
  int loc = start + dist-1;
  int binsize = chr->getBinsize();
  int stepsize = chr->getStepsize();
  int cur_binidx = 1 + ceil((double)(loc-binsize)/stepsize);
  int cur_binbeg = stepsize * (cur_binidx-1)+1;
  int cur_binend = cur_binbeg + binsize-1;
  int chrsize = chr->getChrsize();
  if (cur_binend > chrsize) {
    cur_binend = chrsize;
  } 
  return cur_binidx + chr->getBinstart() - 1; // warning: should depends on step... no actually, cur_binidx already depends on step
}

static bool is_empty_line(const char* buffer)
{
  while (char c = *buffer++) {
    if (c != ' ' || c != '\n' || c != '\t') {
      return false;
    }
  }
  return true;
}

static int build_matrix(bool public_data, int binoffset, chrsize_t ori_binsize, const char* chrsize_file, const char* ifile, const char* oprefix, Format _dummy_format, const std::string& _dummy_bed_prefix, bool binadjust, MatrixFormat matrix_format, chrsize_t step, bool whole_genome, const char* chrA, const char* chrB)
{
  std::ifstream ifs;
  std::ofstream matfs, xbedfs, ybedfs;

  Matrix matrix(binoffset);
  if (int ret = build_matrix_init(matrix, ifile, ifs, oprefix, matfs, xbedfs, ybedfs, chrsize_file, whole_genome, chrA, chrB, ori_binsize, step, binadjust)) {
    return ret;
  }

  char buffer[4096];
  size_t line_cnt = 1;
  size_t line_num = 0;
  std::string lmark, rmark, lorg, rorg, lchr, rchr;
  chrsize_t lstart, rstart;
  if (public_data) {
    char* lchr = NULL;
    char* rchr = NULL;
    while (!ifs.eof()) {
      ifs.getline(buffer, sizeof(buffer)-1);
      line_num++;
      if (is_empty_line(buffer)) {
	continue;
      }
      interaction_parse2(buffer, lchr, lstart, rchr, rstart);
      const AxisChromosome* abs_chr = matrix.getXAxisChromosome(lchr);
      if (!abs_chr) {
	continue;
      }
      const AxisChromosome* ord_chr = matrix.getYAxisChromosome(rchr);
      if (!ord_chr) {
	continue;
      }
      chrsize_t abs_bin = assign_bin(lorg, abs_chr, lstart, 1);
      chrsize_t ord_bin = assign_bin(rorg, ord_chr, rstart, 1);
      switch(matrix_format) {

      case ASIS_MATRIX:
	matrix.add(abs_bin, ord_bin);
	break;

      case UPPER_MATRIX:
	if (abs_bin < ord_bin) {
	  matrix.add(abs_bin, ord_bin);
	} else {
	  matrix.add(ord_bin, abs_bin);
	}
	break;

      case LOWER_MATRIX:
	if (abs_bin > ord_bin) {
	  matrix.add(abs_bin, ord_bin);
	} else {
	  matrix.add(ord_bin, abs_bin);
	}
	break;

      case COMPLETE_MATRIX:
	matrix.add(abs_bin, ord_bin);
	if (abs_bin != ord_bin) {
	  matrix.add(ord_bin, abs_bin);
	}
	break;
      }
      line_cnt++;
      if ((line_cnt % 100000) == 0) {
	std::cerr << line_cnt << std::endl;
      }
    }
  } else {
    // legacy format
    int ldist, rdist;
    char left[512];
    char right[512];
    while (!ifs.eof()) {
      ifs.getline(buffer, sizeof(buffer)-1);
      if (is_empty_line(buffer)) {
	continue;
      }
      get_left_right(buffer, left, right);
      line_num++;
      if (interaction_parse(left, lmark, lorg, lchr, lstart, ldist)) {
	std::cerr << "invalid line #" << line_num << " [" << buffer << "]\n";
	continue;
      }
      if (ldist < 0) {
	continue;
      }
      if (interaction_parse(right, rmark, rorg, rchr, rstart, rdist)) {
	std::cerr << "invalid line #" << line_cnt << " [" << buffer << "]\n";
	continue;
      }
      if (rdist < 0) {
	continue;
      }

      const AxisChromosome* abs_chr = matrix.getXAxisChromosome(lchr);
      if (!abs_chr) {
	continue;
      }
      const AxisChromosome* ord_chr = matrix.getYAxisChromosome(rchr);
      if (!ord_chr) {
	continue;
      }
      chrsize_t abs_bin = assign_bin(lorg, abs_chr, lstart, ldist);
      chrsize_t ord_bin = assign_bin(rorg, ord_chr, rstart, rdist);
      matrix.add(abs_bin, ord_bin);
      /*
      if (symetric && abs_bin != ord_bin) {
	matrix.add(ord_bin, abs_bin);
      }
      */
      line_cnt++;
      if ((line_cnt % 100000) == 0) {
	std::cerr << line_cnt << std::endl;
      }
    }
  }
  matrix.displayXBed(xbedfs);
  if (!whole_genome) {
    matrix.displayYBed(ybedfs);
  }
  matrix.displayMatrix(matfs);
  xbedfs.close();
  ybedfs.close();
  matfs.close();
  return 0;
}

int main(int argc, char* argv[])
{
  chrsize_t step = 1;
  bool binadjust = false;
  MatrixFormat matrix_format = ASIS_MATRIX;
  chrsize_t binsize = 0;
  const char* ifile = NULL;
  const char* oprefix = NULL;
  const char* chrA = NULL;
  const char* chrB = NULL;
  const char* chrsize_file = NULL;
  bool whole_genome = true;
  bool public_data = true; // EV: 2015-01-21: default becomes true
  int binoffset = 1;
  std::string bed_prefix;
  Format format = SPARSE_BED_FMT;

  if (int ret = get_options(argc, argv, binsize, chrsize_file, ifile, oprefix, format, bed_prefix, binadjust, matrix_format, step, whole_genome, public_data, binoffset, chrA, chrB)) {
    if (ret < 0) {
      return 0;
    }
    return ret;
  }

  if (!binsize || !chrsize_file || !ifile || !oprefix) {
    return usage();
  }

  if (whole_genome && (chrA || chrB)) { // should never happen
    return usage();
  }

  if ((chrA && !chrB) || (!chrA && chrB)) {
    return usage();
  }

  return build_matrix(public_data, binoffset, binsize, chrsize_file, ifile, oprefix, format, bed_prefix, binadjust, matrix_format, step, whole_genome, chrA, chrB);
}
