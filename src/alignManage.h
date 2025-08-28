#ifndef _Align_Manage_
#define _Align_Manage_

#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <filesystem>
#include <memory>
#include <cstdint>
#include "bseq.h"
#include "kseq.h"
#include <gatb/gatb_core.hpp>

#define BASE "ATGC"

using std::string;
using std::vector;


struct aln_ls_paf {
    uint32_t index;
    uint16_t mut;
    char qt;
    uint16_t rs;
    uint16_t re;
    uint16_t qs;
    uint16_t qe;
    int as;
    int matchcount;
    int n_cigar;  // 数组长度
    std::unique_ptr<uint16_t[]> cigar; 
    std::unique_ptr<uint64_t[]> base; 
};
  
struct aln_paf_v {
std::vector<std::shared_ptr<aln_ls_paf>> a;  // 自动管理生命周期
size_t n() const { return a.size(); }
size_t m() const { return a.capacity(); }

};


class rManage{
public:
    string rawLongReadsFilePath;
    string rawShortReadsFilePath;
    int rawSequnenceSplitNumber=1;
    int alnSplitNumber=1;
    int seqThreshold=30000;
    int anyFilesSeqNumber;
    int FilesSeqNumber;
    int longReadsNumber=0;
    int shortReadsNumber=0;
    int maxReadsLengh=0;
    string directoryPath;
    string* id_to_name;

    

    
    rManage(string s1,string s2,int ss1,int st);
    static void init_aln_ls_paf(std::shared_ptr<aln_ls_paf>& aln, std::vector<std::string>& reg);

    void ssplitFasta(const string inputFilename,int anyFilesSeqNumber, int numFiles) ;
    void splitFasta(const string& inputFilename, int numFiles ,string directoryPath);
    void countLongReadsNumber();
    void generateTmpFaSplitFile(int i ,string &tmplongFile,string &alnFilename);

    

private:
    void findDirectory();
    void countAlnSplitNumber();
    void runMinimap(const char *target_file,const char *query_file, string& AlnFilename);
};

#endif 