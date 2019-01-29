// HiC-Pro
// Copyright 2015 Institut Curie                               
// Author(s): Eric Viara
// Contact: nicolas.servant@curie.fr
// This software is distributed without any guarantee under the terms of the BSD-3 License

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>


static const int SPARSE_FMT = 0x1;
static const int BED_FMT = 0x2;
static const char* prog;
static bool progress = false;
static bool detail_progress = false;
static bool quiet = false;

static bool NO_DICHO = getenv("NO_DICHO") != NULL;

typedef unsigned int chrsize_t;

const std::string VERSION = "1.2 [2015-10-20]";

const static chrsize_t BIN_NOT_FOUND = (chrsize_t)-1;

class AxisChromosome;

static bool is_empty_line(const char* buffer)
{
  while (char c = *buffer++) {
    if (c != ' ' || c != '\n' || c != '\t') {
      return false;
    }
  }
  return true;
}

static int bed_line_parse(char* buffer, char chr[], chrsize_t& start, chrsize_t& end, const std::string& bedfile, size_t line_num)
{
  if (sscanf(buffer, "%s %u %u", chr, &start, &end) != 3) {
    std::cerr << "bed file \"" << bedfile << "\" at line #" << line_num << " format error\n";
    return 1;
  }
  return 0;
}

struct Interval {
  chrsize_t start;
  chrsize_t end;

  Interval(chrsize_t start = 0, chrsize_t end = 0) : start(start), end(end) { }
};
 
class ChrRegions {

  std::vector<std::string> chr_v;
  std::map<std::string, std::vector<Interval>* > intervals;

public:
  ChrRegions() { }

  int readBedfile(const std::string& bedfile) {
    std::ifstream ifs(bedfile.c_str());
    if (ifs.bad() || ifs.fail()) {
      std::cerr << prog << " cannot open bed file: " << bedfile << " for reading\n";
      return 1;
    }
    char buffer[4096];
    size_t line_num = 0;
    chrsize_t lastend = 0;
    char lastchr[2048] = {0};
    while (!ifs.eof()) {
      ifs.getline(buffer, sizeof(buffer)-1);
      line_num++;
      if (is_empty_line(buffer)) {
	continue;
      }
      chrsize_t start = 0;
      chrsize_t end = 0;
      char chr[2048];
      if (bed_line_parse(buffer, chr, start, end, bedfile, line_num)) {
	return 1;
      }
      if (intervals.find(chr) == intervals.end()) {
	intervals[chr] = new std::vector<Interval>();
	chr_v.push_back(chr);
      }
      /*
      if (lastend != 0 && !strcmp(lastchr, chr) && start != lastend) {
	std::cerr << "warning: discontinuous segment for chromosome " << chr << " at position " << start << " " << end << std::endl;
      }
      */
      if (*lastchr && strcmp(lastchr, chr)) {
	lastend = 0;
      }

      if (lastend != 0 && start < lastend) {
	std::cerr << "error: bedfile not sorted at line #" << line_num << std::endl;
	exit(1);
      }
      strcpy(lastchr, chr);
      lastend = end;
      intervals[chr]->push_back(Interval(start, end));
      if (progress && (line_num % 100000) == 0) {
	std::cerr << '.' << std::flush;
      }
    }
    if (progress) {
      std::cerr << std::endl;
    }
    return 0;
  }

  void displayBed(std::ostream& ofs, const std::vector<AxisChromosome*>& axis_chr) const {
    std::vector<std::string>::const_iterator begin = chr_v.begin();
    std::vector<std::string>::const_iterator end = chr_v.end();
    unsigned int num = 1;
    while (begin != end) {
      const std::string& chrname = *begin;
      std::map<std::string, std::vector<Interval>* >::const_iterator iter = intervals.find(chrname);
      assert(iter != intervals.end());
      const std::vector<Interval>* itv_vect = (*iter).second;
      std::vector<Interval>::const_iterator itv_begin = itv_vect->begin();
      std::vector<Interval>::const_iterator itv_end = itv_vect->end();
      while (itv_begin != itv_end) {
	const Interval& itv = (*itv_begin);
	ofs << chrname << '\t' << itv.start << '\t' << itv.end << '\t' << num << '\n';
	if (progress && (num % 100000) == 0) {
	  std::cerr << '.' << std::flush;
	}
	num++;
	++itv_begin;
      }
      ++begin;
    }
    if (progress) {
      std::cerr << std::endl;
    }
  }

  const std::vector<Interval>* getIntervalsFromChr(const std::string& chr) const {
    std::map<std::string, std::vector<Interval>* >::const_iterator iter = intervals.find(chr);
    if (iter != intervals.end()) {
      return (*iter).second;
    }
    return NULL;
  }
};

class Dichotomic {

  int min, max;
  const std::vector<Interval>& intervals;

public:
  Dichotomic(const std::vector<Interval>& intervals) : intervals(intervals) {
    //min = middle(intervals[0]);
    //max = middle(intervals[intervals.size()-1]);
    min = 0;
    max = intervals.size()-1;
  }

  static chrsize_t middle(const Interval& itv) {
    return (itv.start+1 + itv.end) / 2;
  }

  int find(chrsize_t value) {
    int l = min;
    int r = max;
    int n = 0;
    while (l <= r) {
      n = (l + r) >> 1;
      const Interval& itv = intervals[n];
      if (value >= itv.start+1 && value <= itv.end) {
	return n;
      }

      int x = middle(itv) - value;
      
      if (x < 0) {
	l = n + 1;
      } else {
	r = n - 1;
      }
      //std::cout << "l: " << l << '\n';
      //std::cout << "r: " << r << '\n';
    }

    return -1;
  }
};

class Chromosome {

private:
  static std::unordered_map<std::string, Chromosome*> chr_map;

  void computeSizes(chrsize_t ori_binsize, chrsize_t step, bool binadjust, const ChrRegions* chr_regions);

  std::string name;

  chrsize_t chrsize;

  chrsize_t binsize;
  chrsize_t stepsize;
  chrsize_t bincount;

  const ChrRegions* chr_regions;

public:
  Chromosome(const std::string& name, chrsize_t chrsize, chrsize_t ori_binsize, chrsize_t step, bool binadjust, const ChrRegions* chr_regions) : name(name), chrsize(chrsize), chr_regions(chr_regions) {
    computeSizes(ori_binsize, step, binadjust, chr_regions);
    assert(chr_map.find(name) == chr_map.end());
    chr_map[name] = this;
  }

  void adjustBinsize(chrsize_t ori_binsize, const chrsize_t step);

  const std::string& getName() const {return name;}
  chrsize_t getChrsize() const {return chrsize;}
  chrsize_t getBinsize() const {return binsize;}
  chrsize_t getStepsize() const {return stepsize;}
  chrsize_t getBincount() const {return bincount;}

  const ChrRegions* getChrRegions() const {return chr_regions;}

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
    /*
    if (verbose) {
      std::cerr << "AxisChromosome: " << chr->getName() << " " << binstart << " " << binend << " " << chr->getBincount() << std::endl;
    }
    */
  }

  chrsize_t getBinstart() const {return binstart;}
  chrsize_t getBinend() const {return binend;}
  chrsize_t getChrsize() const {return chr->getChrsize();}
  chrsize_t getBinsize() const {return chr->getBinsize();}
  chrsize_t getStepsize() const {return chr->getStepsize();}
  chrsize_t getBincount() const {return chr->getBincount();}

  const Chromosome* getChromosome() const {return chr;}

  chrsize_t assign_bin(const std::string& org, chrsize_t start) const {
    const ChrRegions* chr_regions = chr->getChrRegions();
    if (chr_regions != NULL) {
      const std::vector<Interval>* intervals = chr_regions->getIntervalsFromChr(chr->getName());
      assert(intervals != NULL);

      if (!NO_DICHO) {
	Dichotomic dicho(*intervals);
	int where = dicho.find(start);
	if (where < 0) {
	  if (!quiet) {
	    std::cerr << "warning: no bin at position " << chr->getName() << ":" << start << std::endl;
	  }
	  return BIN_NOT_FOUND;
	}
	return where + getBinstart();
      }

      std::vector<Interval>::const_iterator begin = intervals->begin();
      std::vector<Interval>::const_iterator end = intervals->end();

      chrsize_t binidx = 1;
      while (begin != end) {
	const Interval& itv = *begin;
	if (start >= itv.start+1 && start <= itv.end) {
	  break;
	}
	++binidx;
	++begin;
      }
      
      return binidx + getBinstart() - 1;
    }

    int loc = (int)start;
    int binsize = getBinsize();
    int stepsize = getStepsize();
    int cur_binidx = 1 + ceil((double)(loc-binsize)/stepsize);
    int cur_binbeg = stepsize * (cur_binidx-1)+1;
    int cur_binend = cur_binbeg + binsize-1;
    int chrsize = getChrsize();
    if (cur_binend > chrsize) {
      cur_binend = chrsize;
    } 
    return cur_binidx + getBinstart() - 1;
  }
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
    size_t line_total = 0;
    if (progress) {
      while (begin != end) {
	const std::map<chrsize_t, chrsize_t>& line = (*begin).second;
	line_total += line.size();
	++begin;
      }
      begin = mat.begin();
    }

    size_t line_cnt = 1;
    if (progress) {
      std::cerr << "\n=================\n";
      std::cerr << " Dumping matrix\n";
      std::cerr << "=================\n\n";
    }
    size_t modulo = line_total / 1000;
    while (begin != end) {
      chrsize_t abs = (*begin).first;
      const std::map<chrsize_t, chrsize_t>& line = (*begin).second;
      std::map<chrsize_t, chrsize_t>::const_iterator bb = line.begin();
      std::map<chrsize_t, chrsize_t>::const_iterator ee = line.end();
      while (bb != ee) {
	if (progress && (line_cnt % modulo) == 0) {
	  double percent = (double(line_cnt)/line_total)*100;
	  std::cerr << "" << percent << "% " << line_cnt << " / " << line_total << std::endl;
	}
	ofs << abs << '\t' << (*bb).first << '\t' << (*bb).second << '\n';
	line_cnt++;
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

  const std::vector<AxisChromosome*>& getXAxisChromosomes() {return axis_chr_abs;}
  const std::vector<AxisChromosome*>& getYAxisChromosomes() {return axis_chr_ord;}
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

void Chromosome::computeSizes(chrsize_t ori_binsize, chrsize_t step, bool binadjust, const ChrRegions* chr_regions)
{
  if (NULL != chr_regions) {
    const std::vector<Interval>* intervals = chr_regions->getIntervalsFromChr(name);
    assert(intervals != NULL);
    bincount = intervals->size();
    /*
    if (verbose) {
      std::cerr << name << " bincount: " << bincount << std::endl;
    }
    */
  } else {
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
    /*
    if (verbose) {
      std::cerr << name << " sizes: " << chrsize << " " << binsize << " " << stepsize << " " << bincount << std::endl;
    }
    */
  }
}

static int usage(int ret = 1)
{
  std::cerr << "\nusage: " << prog << " --binsize BINSIZE|--binfile --chrsizes FILE --ifile FILE\n";
  std::cerr << "       --oprefix PREFIX [--binadjust] [--step STEP] [--binoffset OFFSET]\n";
  std::cerr << "       [--matrix-format asis|upper|lower|complete][--chrA CHR... --chrB CHR...] [--quiet] [--progress] [--detail-progress]\n";
  std::cerr << "\nusage: " << prog << " --version\n";
  std::cerr << "\nusage: " << prog << " --help\n";
  return ret;
}

static int help()
{
  (void)usage();
  std::cerr << "\nOPTIONS\n\n";
  std::cerr << "  --version              : display version\n";
  std::cerr << "  --binsize BINSIZE      : bin size\n";
  std::cerr << "  --binfile BEDFILE      : bed file containing bins (chr start end)\n";
  std::cerr << "  --chrsizes FILE        : file containing chromosome sizes\n";
  std::cerr << "  --ifile FILE           : input interaction file\n";
  std::cerr << "  --oprefix PREFIX       : output prefix of generated files (matrix and bed)\n";
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
  std::cerr << "  --quiet                : do not display any warning\n";
  std::cerr << "  --progress             : display progress\n";
  std::cerr << "  --detail-progress      : display detail progress (needs preliminary steps consuming time)\n";
  return -1;
}

enum MatrixFormat {
  ASIS_MATRIX = 1,
  UPPER_MATRIX,
  LOWER_MATRIX,
  COMPLETE_MATRIX
};
  
static int get_options(int argc, char* argv[], chrsize_t& binsize, const char*& binfile, const char*& chrsize_file, const char*& ifile, const char*& oprefix, Format& format, std::string& bed_prefix, bool& binadjust, MatrixFormat& matrix_format, chrsize_t& step, bool& whole_genome, int& binoffset, const char*& chrA, const char*& chrB)
{
  prog = argv[0];
  for (int ac = 1; ac < argc; ++ac) {
    const char* opt = argv[ac];
    if (*opt == '-') {
      if (!strcmp(opt, "--binadjust")) {
	binadjust = true;
      } else if (!strcmp(opt, "--version")) {
	std::cout << "build_matrix version " << VERSION << "\n";
	exit(0);
      } else if (!strcmp(opt, "--progress")) {
	progress = true;
      } else if (!strcmp(opt, "--quiet")) {
	quiet = true;
      } else if (!strcmp(opt, "--detail-progress")) {
	progress = true;
	detail_progress = true;
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
      } else if (!strcmp(opt, "--binfile")) {
	if (ac == argc-1) {
	  return usage();
	}
	binfile = argv[++ac];
      } else if (!strcmp(opt, "--binsize")) {
	if (ac == argc-1) {
	  return usage();
	}
	binsize = atoi(argv[++ac]);
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
	std::cerr << '\n' << prog << ": unknown option " << opt << std::endl;
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

static int interaction_parse(char* buffer, char*& lchr, chrsize_t& lstart, char*& rchr, chrsize_t& rstart)
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

static char p_buffer[512000];

static int build_matrix_init(Matrix& matrix, const char* ifile, std::ifstream& ifs, const std::string& oprefix, std::ofstream& matfs, std::ofstream& xbedfs, std::ofstream& ybedfs, const char* chrsize_file, bool whole_genome, const char* chrA, const char* chrB, chrsize_t ori_binsize, const char* binfile, chrsize_t step, bool binadjust, ChrRegions*& chr_regions, size_t& line_total)
{
  ifs.open(ifile);
  if (ifs.bad() || ifs.fail()) {
    std::cerr << prog << " cannot open interaction file: " << ifile << " for reading\n";
    return 1;
  }

  if (detail_progress) {
    if (progress) {
      std::cerr << "\n======================================\n";
      std::cerr << " Getting information for progress bar\n";
      std::cerr << "======================================\n\n";
    }
    std::cerr << std::setprecision(2) << std::fixed;
    int fd = open(ifile, O_RDONLY);
    struct stat st;
    assert(fstat(fd, &st) == 0);
    assert(fd >= 0);
    int nn;
    int cnt = 1;
    while ((nn = read(fd, p_buffer, sizeof(p_buffer))) > 0) {
      const char *p = p_buffer;
      while (nn-- > 0) {
	if (*p++ == '\n') {
	  line_total++;
	}
      }
      if ((cnt % 200) == 0) {
	std::cerr << '.' << std::flush;
      }
      cnt++;
    }
    std::cerr << std::endl;
    close(fd);
  }
  
  std::ifstream chrsizefs;
  chrsizefs.open(chrsize_file);
  if (chrsizefs.bad() || chrsizefs.fail()) {
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
  if (!whole_genome) {
    //std::string xbedlink;
    //size_t pos = xbedfile.rfind('/');
    //if (pos != std::string::npos) {
    //  xbedlink = xbedfile.substr(pos+1);
    //} else {
    //  xbedlink = xbedfile;
    //}
    //unlink(ybedfile.c_str());
    //if (symlink(xbedlink.c_str(), ybedfile.c_str())) {
    //  std::cerr << prog << " cannot created link: " << ybedfile << "\n";
    //  return 1;
    //}
    //} else {
    ybedfs.open(ybedfile);
    if (ybedfs.bad() || ybedfs.fail()) {
      std::cerr << prog << " cannot open file: " << ybedfile << " for writing\n";
      return 1;
    }
  }

  chr_regions = NULL;
  if (NULL != binfile) {
    chr_regions = new ChrRegions();
    if (progress) {
      std::cerr << "\n=================\n";
      std::cerr << " Reading binfile\n";
      std::cerr << "=================\n\n";
    }
    if (chr_regions->readBedfile(binfile)) {
      return 1;
    }
  }

  std::vector<const Chromosome*> all_chr_v;
  while (!chrsizefs.eof()) {
    std::string buffer;
    getline(chrsizefs, buffer);

    chrsize_t chrsize;
    std::istringstream istr(buffer);
    std::string name;
    istr >> name >> chrsize;
    if (!istr.fail()) {
      Chromosome* chromosome = new Chromosome(name, chrsize, ori_binsize, step, binadjust, chr_regions);
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

static int build_matrix(int binoffset, chrsize_t ori_binsize, const char* binfile, const char* chrsize_file, const char* ifile, const char* oprefix, Format _dummy_format, const std::string& _dummy_bed_prefix, bool binadjust, MatrixFormat matrix_format, chrsize_t step, bool whole_genome, const char* chrA, const char* chrB)
{
  std::ifstream ifs;
  std::ofstream matfs, xbedfs, ybedfs;

  Matrix matrix(binoffset);
  ChrRegions *chr_regions = NULL;
  size_t line_total = 0;
  if (int ret = build_matrix_init(matrix, ifile, ifs, oprefix, matfs, xbedfs, ybedfs, chrsize_file, whole_genome, chrA, chrB, ori_binsize, binfile, step, binadjust, chr_regions, line_total)) {
    return ret;
  }

  if (progress) {
    std::cerr << "\n=================\n";
    std::cerr << " Building matrix\n";
    std::cerr << "=================\n\n";
  }
  size_t line_cnt = 1;
  size_t line_num = 0;
  char buffer[4096];
  std::string lmark, rmark, lorg, rorg;
  while (!ifs.eof()) {
    ifs.getline(buffer, sizeof(buffer)-1);
    line_num++;
    if (is_empty_line(buffer)) {
      continue;
    }
    chrsize_t lstart = 0;
    chrsize_t rstart = 0;
    char* lchr = NULL;
    char* rchr = NULL;
    interaction_parse(buffer, lchr, lstart, rchr, rstart);
    const AxisChromosome* abs_chr = matrix.getXAxisChromosome(lchr);
    if (!abs_chr) {
      continue;
    }
    const AxisChromosome* ord_chr = matrix.getYAxisChromosome(rchr);
    if (!ord_chr) {
      continue;
    }
    chrsize_t abs_bin = abs_chr->assign_bin(lorg, lstart);
    if (abs_bin == BIN_NOT_FOUND) {
      continue;
    }
    chrsize_t ord_bin = ord_chr->assign_bin(rorg, rstart);
    if (ord_bin == BIN_NOT_FOUND) {
      continue;
    }
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
    if (progress && (line_cnt % 100000) == 0) {
      if (detail_progress) {
	double percent = (double(line_cnt)/line_total)*100;
	std::cerr << "" << percent << "% " << line_cnt << " / " << line_total << std::endl;
      } else {
	std::cerr << line_cnt << std::endl;
      }
    }
  }

  if (progress) {
    std::cerr << "\n==================\n";
    std::cerr << " Dumping bedfiles\n";
    std::cerr << "==================\n\n";
  }

  if (NULL != chr_regions) {
    chr_regions->displayBed(xbedfs, matrix.getXAxisChromosomes());
    if (!whole_genome) {
      chr_regions->displayBed(ybedfs, matrix.getYAxisChromosomes());
    }
  } else {
    matrix.displayXBed(xbedfs);
    if (!whole_genome) {
      matrix.displayYBed(ybedfs);
    }
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
  const char* binfile = NULL;
  bool whole_genome = true;
  int binoffset = 1;
  std::string bed_prefix;
  Format format = SPARSE_BED_FMT;

  if (int ret = get_options(argc, argv, binsize, binfile, chrsize_file, ifile, oprefix, format, bed_prefix, binadjust, matrix_format, step, whole_genome, binoffset, chrA, chrB)) {
    if (ret < 0) {
      return 0;
    }
    return ret;
  }

  if (!binsize && !binfile) {
    std::cerr << '\n';
    std::cerr << prog << ": missing --binsize or --binfile option\n";
    return usage();
  }

  if (!chrsize_file) {
    std::cerr << '\n';
    std::cerr << prog << ": missing --chrsizes option\n";
    return usage();
  }

  if (!ifile) {
    std::cerr << '\n';
    std::cerr << prog << ": missing --ifile option\n";
    return usage();
  }

  if (!oprefix) {
    std::cerr << '\n';
    std::cerr << prog << ": missing --oprefix option\n";
    return usage();
  }

  if ((chrA && !chrB) || (!chrA && chrB)) {
    std::cerr << '\n';
    std::cerr << prog << ": options --chrA and --chrB must be set simultanously\n";
    return usage();
  }

  if (binfile && binsize) {
    std::cerr << '\n';
    std::cerr << prog << ": options --binfile and --binsize cannot be set simultanously\n";
    return usage();
  }

  return build_matrix(binoffset, binsize, binfile, chrsize_file, ifile, oprefix, format, bed_prefix, binadjust, matrix_format, step, whole_genome, chrA, chrB);
}
