#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stack>
#include <map>
#include <malloc.h>
#include <getopt.h>


#include <errno.h>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "mmpriv.h"
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"

#include "khash.h"
#include "alignManage.h"
#include <stdexcept>


typedef struct {
  std::map<string,vector<uint32_t>> i;
  std::map<string,vector<uint32_t>> d;
  std::map<string,vector<uint32_t>> x;
}ssnp;
class correctRead{
public:
  aln_paf_v alnPaf;
  std::unique_ptr<char[]> lread;
  float threshold;
  int correctMode;
  int long_idx;
  int read_len;
  int n_qReads;
  std::unique_ptr<char []> newLread;
  std::map<uint32_t,ssnp> pos_snp;
  std::map<uint32_t,uint32_t> mutReads;
  std::map<uint32_t,vector<uint32_t>> map;

  correctRead(aln_paf_v ap,const Sequence& seq,int maxLength,float threshold);
  void prepare_snp();
  static void reverse(char *target, char *source, int len);
  static int msa(std::string pacbioFile,std::string illuminaFile,std::string outputFile,int msa_thread,float threshold,int splitNumber);

private:
  bool find_mut(int sidx,int matchcount);
  bool haveKey(std::map<string,vector<uint32_t>> &b,char* key,uint32_t sid,int mode);
  bool add_mut(int sid);
  void cigar_ite(uint64_t baseCode,uint16_t snp_len,char type,uint32_t &pos,uint32_t &srPos,uint32_t sid,char qt,int ite);
  void find_max(std::map<string,vector<uint32_t>> &b,int &maxdepth,int &type,string &base,int &rsnp_len,int t,int pos,int &num);
  
  uint16_t countContainsPos(uint32_t pos);
  void generate_mut();
  void update_lread();
  void update_snp(uint64_t* base,int lcigar,uint16_t* cigar,uint32_t &pos,uint32_t &srPos,int sid,char qt,int ite);
   void copy_upper_case(char *target, char *source, int len);
  
};