#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stack>
#include <map>
#include <sys/mman.h>
#include <fcntl.h>  
#include <unistd.h> 

#include <malloc.h>
#include <getopt.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <gatb/gatb_core.hpp>
#include "mmpriv.h"
#include "minimap.h"
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "khash.h"

#define MIN_KMER_LEN 4
#define MIN_SOLID_THR 1

#define MIN_ERROR_RATE 0.0
#define MAX_ERROR_RATE 0.5

#define MIN_TRIALS 1
#define MAX_TRIALS 100

#define MIN_NB_BRANCH 1
#define MAX_NB_BRANCH 10000

#define MIN_THREADS 0
#define MAX_THREADS 64


// other constants
//#define MAX_READ_LEN 500000
#define MAX_PATH_LEN 1000

// default values for correction
#define DEF_ERROR_RATE 0.40
#define DEF_MAX_BRANCH 200
#define DEF_TRIALS 5
#define DEF_THREADS 0
#define DEF_MAX_ABUNDANCE 2147483647

// alignment scores
#define ALIGN_MATCH 1
#define ALIGN_MISMATCH -3
#define ALIGN_INDEL -2

// Constants for trimming and splitting
#define DEF_MIN_REGION_LG 100

// Constants for statistics
#define STAT_FOUND 0
#define STAT_FOUND_LEN1 1
#define STAT_TOOLONG 2
#define STAT_EXPLOSION 3
#define STAT_NOPATH 4

#define STAT_TAIL "TAIL"
#define STAT_END2END "END2END"
#define STAT_GAPEXTEND "GAPEXTEND"

// File extensions
#define HDF5_EXT ".h5"
#define FASTA_EXT ".fa"
#define FASTQ_EXT ".fq"

// Utilities 

typedef struct aln_ls_paf{ 
	uint32_t index;
	uint32_t mut;
	char qt;
	uint32_t rs;
	uint32_t re;
  uint32_t qs;
	uint32_t qe;
	uint32_t matchcount;
	uint32_t n_cigar;
	uint32_t *cigar;

}aln_ls_paf;

typedef struct {
	size_t n, m; 
	aln_ls_paf *a;
} aln_paf_v;

typedef struct {
	size_t n, m; 
	kstring_t *a;
} Short_reads;

void init_aln_ls_paf(aln_ls_paf* aln, uint32_t _index, char _qt, mm_reg1_t reg) {
	aln->index = _index;
	aln->rs=reg.rs;
	aln->re=reg.re;
  aln->qs=reg.qs;
	aln->qe=reg.qe;
  aln->qt = _qt;
  aln->n_cigar = reg.p->n_cigar;
	aln->matchcount=reg.mlen;
	aln->mut=0;
  aln->cigar = (uint32_t *)malloc(sizeof(uint32_t) * aln->n_cigar);
  for (uint32_t i = 0; i < aln->n_cigar; ++i) {
      aln->cigar[i] = reg.p->cigar[i];
  }/**/
}
// min macro
#define MIN(a,b) ((a) < (b) ? (a) : (b))

// max macro
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// complement reverse a dna seq
void reverse(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    switch(source[len-i-1]) {
      case 'a':
        target[i] = 't';
        break;
      case 'c':
        target[i] = 'g';
        break;
      case 'g':
        target[i] = 'c';
        break;
      case 't':
        target[i] = 'a';
        break;
      case 'n':
        target[i] = 'n';
        break;
      case 'A':
        target[i] = 'T';
        break;
      case 'C':
        target[i] = 'G';
        break;
      case 'G':
        target[i] = 'C';
        break;
      case 'T':
        target[i] = 'A';
        break;
      case 'N':
        target[i] = 'N';
        break;
    }
  }
  target[len] = '\0';
}
//////////////////////////////////////////////////
// check wheter a file is present and readable 
bool is_readable( const std::string & file ) 
{ 
    std::ifstream fichier( file.c_str() ); 
    return !fichier.fail(); 
} 


typedef struct {
  int n,m;
  int* a;
}int_v;
struct LineOffset {
  size_t start;
  size_t end;
};


typedef struct {
  std::map<string,int_v> i;
  std::map<string,int_v> d;
  std::map<string,int_v> x;
}ssnp;


bool haveKey(std::map<string,int_v> &b,char* key,int sid,int mode){
  string base=key;
  auto search=b.find(base);
  if(mode==1){
    if(search==b.end()){
      int_v v;
      kv_init(v);
      kv_push(int,v,sid);
      b.insert({base,v}); 
      return false;  
    }
    int_v v=search->second;
    kv_push(int,v,sid);
    search->second=v;
  }
  else{
    kv_remove(int, search->second, sid);
    int_v ids=search->second;
    if (ids.n==0){
      b.erase(search);
      return true;
    }
  }
  return false; 
}

bool find_mut(std::map<int,int> id_num,int sidx,int matchcount,float threshold){
  auto search=id_num.find(sidx);
  if(search!=id_num.end()){
    float num=search->second;
    float fre=num/matchcount;
    if(fre>threshold)
      return true;
  }
  return false;
}

bool add_map(int pos,int pos2,std::map<int,int_v> &map){
  auto search=map.find(pos);
  if(search==map.end()){
    int_v array;
    kv_init(array);
    kv_push(int,array,pos2);
    map.insert({pos,array});
    return false;
  } 
  int_v array=map[pos];
  kv_push(int,array,pos2);
  search->second=array;
  return true;
}

bool add_mut(std::map<int,int> &mut,int sid){
  auto search=mut.find(sid);
  if(search==mut.end()){
    mut.insert({sid,1});
    return false;
  }
  search->second=search->second+1;
  return true;  
}

void getLineOffsets(const std::string& filePath,std::vector<LineOffset> &lineOffsets,int num) {
  int step;
  if (filePath.substr(filePath.size() - 3) == ".fa" || filePath.substr(filePath.size() - 6) == ".fasta") {
      step = 2; 
  } else if (filePath.substr(filePath.size() - 3) == ".fq" || filePath.substr(filePath.size() - 6) == ".fastq") {
      step = 4; 
  } else {
      std::cerr << "Invalid filePath: " << filePath << std::endl;

  }


  int fd = open(filePath.c_str(), O_RDONLY);
  if (fd == -1) {
      std::cerr << "Failed to open file: " << filePath << std::endl;

  }

  size_t fileSize = lseek(fd, 0, SEEK_END);
  char* fileData = (char*)mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
  if (fileData == MAP_FAILED) {
      std::cerr << "Failed to map file to memory: " << filePath << std::endl;
      close(fd);

  }


  
  int lineStart = 0;
  int lineCount = 0;

  for (int i = 0; i < fileSize; i++) {
      if (fileData[i] == '\n') {
          lineCount++;

          if (lineCount % step == 2) { 
              LineOffset offset;
              offset.start = lineStart;
              offset.end = i; 
              lineOffsets.push_back(offset);
              num++;
          }
          lineStart = i + 1; 
      }
  }

 
  if (lineStart < fileSize && lineCount % step == 1) {
      LineOffset offset;
      offset.start = lineStart;
      offset.end = fileSize;
      lineOffsets.push_back(offset);
  }

  munmap(fileData, fileSize);
  close(fd);

}
std::string readindex(std::string filePath, std::vector<LineOffset> lineOffsets, int index) {
  
  //cout<<"readindex:"<<index<<endl;
  if (index <0 || index >= lineOffsets.size()) {
      std::cerr << "Invalid index: " << index << std::endl;
      return "";
  }

  int fd = open(filePath.c_str(), O_RDONLY);
  size_t fileSize = lseek(fd, 0, SEEK_END);
  char* fileData = (char*)mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);

  const LineOffset& offset = lineOffsets[index];
  std::string line(fileData + offset.start, fileData + offset.end);

  munmap(fileData, fileSize);
  close(fd);
  return line;
}
void readShorts(char* sread,const char* tmp,int len,char qt){
  if(qt=='-'){
    char* buffer = new char[len+1];
    std::strcpy(buffer, tmp);
    reverse(sread, buffer, len);
    delete[] buffer;
  }
  else{
    std::strcpy(sread, tmp);
  }
  sread[len]='\0';
  
}

void cigar_ite(int snp_len,char type,int &pos,int &srPos,std::map<int,ssnp> &pos_snp,char *lread,char* sread,int sid,int ite){

  if(type=='=')
  {
    pos=pos+snp_len;
    srPos=srPos+snp_len;
    return ;
  }

  char* snp_base=new char[snp_len+1];

  for(int i=0;i<snp_len;i++){
    if(type=='I'||type=='X'){
      snp_base[i]=sread[srPos+i];
    }
    else{
      snp_base[i]=lread[srPos+i];
    }
  }  
  snp_base[snp_len]='\0'; 
  bool eras=false; 

  if(type=='I'){
    eras=haveKey(pos_snp[pos].i,snp_base,sid,ite);
    srPos=srPos+snp_len;
  }
  else if(type=='D'){
    eras=haveKey(pos_snp[pos].d,snp_base,sid,ite);
    pos=pos+snp_len;
  }
  else{
    eras=haveKey(pos_snp[pos].x,snp_base,sid,ite);
    srPos=srPos+snp_len;
    pos=pos+snp_len;
  }
  if(eras){
    int tmp=pos;
    if(type!='I') tmp=pos-snp_len;
    if(pos_snp[tmp].i.empty()&&pos_snp[tmp].d.empty()&&pos_snp[tmp].x.empty()){
      pos_snp.erase(pos_snp.find(tmp));

    }
  }
}   


void find_max(std::map<string,int_v> &b,int &maxdepth,int &type,string &base,int &rsnp_len,int t,char* lread,int pos,int &num){
  if(b.empty()) return;
  
  for (auto it = b.begin(); it != b.end(); ++it) {
   
    string snp_base=it->first;
    int snp_len=snp_base.length();
    int_v snp_reads=it->second;
    int depth=snp_reads.n;
    num+=depth;
    char* lbase=new char[snp_len+1];
    for(int i=0;i<snp_len;i++)
      lbase[i]=lread[pos+i];
    lbase[snp_len]='\0';  
    if(depth<maxdepth) continue;
    else if((depth==maxdepth&&(type==2||type==0)&&snp_base.c_str()==lbase)||depth>maxdepth){
      maxdepth=depth;
      type=t;
      base=lbase;
      rsnp_len=snp_len;
    }
  }

}
void generate_mut(char* lread,int read_len,std::map<int,ssnp> pos_snp,std::map<int,int> &mutReads,std::map<int,int_v> map){

  for(int pos=0;pos<read_len;pos++){
    int type=0;
    string base;
    int maxdepth=0,rsnp_len;
    int snp_num=0,map_num=0;
    float c=0.0;
 
    if(pos_snp.find(pos)!=pos_snp.end()){
      std::map<string,int_v> mis_snp=pos_snp[pos].x;
      find_max(mis_snp,maxdepth,type,base,rsnp_len,2,lread,pos,snp_num);  

      for(auto map_ite=map.begin();map_ite!=map.end();++map_ite){
        int start=map_ite->first;
        if(pos<start){break;}
        else{
          int_v end=map_ite->second;
          int i=0;
          for(int i=0;i<end.n;i++){
            if(pos<=end.a[i])
              map_num++;
          }
        }
      }
      c=float(snp_num)/float(map_num);  
      if(c<=0.5 or (map_num-snp_num)>5){
        std::map<string,int_v> mis_snp=pos_snp[pos].x;
        for(auto mis=mis_snp.begin();mis!=mis_snp.end();++mis){
          int_v snp_reads=mis->second;
          float snp_ratio=float(snp_reads.n)/float(map_num);
          if(snp_reads.n>5||snp_ratio>=0.5){
            for(int p=0;p<snp_reads.n;p++){
              add_mut(mutReads,snp_reads.a[p]);
            }
          }
        }        
      } 
    }
  }
}
void update_lread(char* lread,int read_len,char* newLread,std::map<int,ssnp> pos_snp,std::map<int,int_v> map){
  int newpos=0;

  for(int pos=0;pos<read_len;pos++){
    int type=0;
    string base;
    int maxdepth=0,rsnp_len;
    int snp_num=0,map_num=0;
    int c=0;

    if(pos_snp.find(pos)!=pos_snp.end()){
      find_max(pos_snp[pos].x,maxdepth,type,base,rsnp_len,2,lread,pos,snp_num);  
      find_max(pos_snp[pos].i,maxdepth,type,base,rsnp_len,0,lread,pos,snp_num);
      find_max(pos_snp[pos].d,maxdepth,type,base,rsnp_len,1,lread,pos,snp_num);
      
      for(auto map_ite=map.begin();map_ite!=map.end();++map_ite){
        int start=map_ite->first;
        if(pos<start){break;}
        else{
          int_v end=map_ite->second;
          int i=0;
          for(int i=0;i<end.n;i++){
            if(pos<=end.a[i])
              map_num++;
          }
        }
      }
      c=map_num-snp_num;  
    }
   
    if(c>=5||maxdepth<3){
      newLread[newpos]=lread[pos];
      newpos++;
      continue;
    }

    if(type==0){
      for(int i=0;i<rsnp_len;i++){
        newLread[newpos]=base[i];
        newpos++;
      }
      newLread[newpos]=lread[pos];
      newpos++;
    }
    else if(type==1){
      pos=pos+rsnp_len-1;
    }
    else{
      for(int i=0;i<rsnp_len;i++){
        newLread[newpos]=base[i];
        newpos++;
      }
      pos=pos+rsnp_len-1;
    }
  }
  newLread[newpos]='\0';  
}
void update_snp(int lcigar,uint32_t* cigar,std::map<int,ssnp> &pos_snp, aln_ls_paf *ls_detail,int &pos,int &srPos,char *lread,char *sread,int ite){

  for(int j=0;j<lcigar;++j){
    int snp_len=cigar[j]>>4;
    char type=MM_CIGAR_STR[cigar[j]&0xf];
    if(type!='='){
      if(pos_snp.find(pos)==pos_snp.end()){
        ssnp all_snp={{},{},{}};
        pos_snp[pos]=all_snp;
      }
    }
    cigar_ite(snp_len,type,pos,srPos,pos_snp,lread,sread,ls_detail->index,ite);
  }
  
}
void long_correct(Sequence seq,aln_paf_v alnPaf,char* newLread,std::map<int,int> lowPos,int threshold,std::string illuminaFile,std::vector<LineOffset> lineOffsets){
  char *lread = seq.getDataBuffer();
  int read_len = seq.getDataSize();
  int n_qReads=alnPaf.n;
  std::map<int,ssnp> pos_snp;
  std::map<int,int> have_short;
  std::map<int,int> mutReads;
  std::map<int,int_v> map;
  for(int i=0;i<n_qReads;i++){ 
    aln_ls_paf *ls_detail=&(alnPaf.a[i]);
    auto have=have_short.find(ls_detail->index);
    
    if(have!=have_short.end()){
      alnPaf.a[i].mut=1;
      continue;
    }
    have_short.insert({ls_detail->index,1});
    int lcigar=ls_detail->n_cigar;
    uint32_t* cigar=ls_detail->cigar;
    int srPos=0;
    int pos=ls_detail->rs;
    int pos2=ls_detail->re;
    std::string tmpSreads=readindex(illuminaFile,lineOffsets,ls_detail->index);
    int slen=tmpSreads.size();

    char* sread=new char[slen+1];

    auto search=map.find(pos);
    if(search==map.end()){
      int_v v;
      kv_init(v);
      kv_push(int,v,pos2);
      map.insert({pos,v});
    }
    else{
      int_v v=map[pos];
      kv_push(int,v,pos2);
      search->second=v;
    }
    readShorts(sread,tmpSreads.c_str(),slen,ls_detail->qt);
    update_snp(lcigar,cigar,pos_snp,ls_detail,pos,srPos,lread,sread,1);
  }

  generate_mut(lread,read_len,pos_snp,mutReads,map); 
  int nomutAlin=0;
  for(int i=0;i<n_qReads;i++){ 
    if(alnPaf.a[i].mut==1) {
      //cout<<alnPaf.a[i].index<<endl;
      continue;
    }
    aln_ls_paf* ls_detail=&(alnPaf.a[i]);
    int lcigar=ls_detail->n_cigar;
    uint32_t* cigar=ls_detail->cigar;
    int srPos=0;
    int pos=ls_detail->rs;
    int pos2=ls_detail->re;
    std::string tmpSreads=readindex(illuminaFile,lineOffsets,ls_detail->index);
    int slen=tmpSreads.size();
    char* sread=new char[slen+1];
    readShorts(sread,tmpSreads.c_str(),slen,ls_detail->qt);


    bool is_mut = find_mut(mutReads, ls_detail->index, ls_detail->matchcount,threshold);
    if (is_mut) {
      nomutAlin++;
      kv_remove(int,map[pos],pos2);
      if(map[pos].n==0){
        map.erase(map.find(pos));
      }

      for(int j=0;j<lcigar;++j){
          int snp_len=cigar[j]>>4;
          char type=MM_CIGAR_STR[cigar[j]&0xf];
          cigar_ite(snp_len,type,pos,srPos,pos_snp,lread,sread,ls_detail->index,2);
      }     
      //continue;
    }
  }
  //cout<<seq.getComment()<<" all:"<<n_qReads<<" filt:"<<nomutAlin<<endl;
  update_lread(lread,read_len,newLread,pos_snp,map);  

  pos_snp.clear();
  map.clear();
  have_short.clear();
}

