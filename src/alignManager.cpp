#include "alignManage.h"
#include "minimap.h"

rManage::rManage(string s1,string s2,int ss1,int st)
{
    cout<<"initialization..."<<endl;
    rawLongReadsFilePath=s1;
    rawShortReadsFilePath=s2;
    rawShortReadsFilePath=rawShortReadsFilePath.insert(s2.size()-3,"rep");
    rawSequnenceSplitNumber=ss1;
    alnSplitNumber=1;
    seqThreshold=st;
    longReadsNumber=0;
    shortReadsNumber=0;
    maxReadsLengh=0;
    if (!std::filesystem::exists(rawShortReadsFilePath)){
        string replace_cmd1="seqkit replace -p '.*' -r 'SR{nr}' -w 0 "+s2+" > "+rawShortReadsFilePath;
        system(replace_cmd1.c_str());
    }
    findDirectory();
    countLongReadsNumber();
    
    cout<<"split long reads..."<<endl;
    //cout<<rawLongReadsFilePath<<" "<<rawShortReadsFilePath<<" "<<rawSequnenceSplitNumber<<" "<<directoryPath<<endl;
    splitFasta(rawLongReadsFilePath,rawSequnenceSplitNumber,directoryPath);
    
}
void rManage::init_aln_ls_paf(std::shared_ptr<aln_ls_paf>& aln, std::vector<std::string>& reg) {
    aln->index = std::stoi(reg[1]);
    aln->qt = reg[2][0];
    aln->rs = std::stoi(reg[3]);
    aln->re = std::stoi(reg[4]);
    aln->qs = std::stoi(reg[5]);
    aln->qe = std::stoi(reg[6]);
    aln->as = std::stoi(reg[7]);
    aln->matchcount = std::stoi(reg[8]);

    aln->n_cigar = std::stoi(reg[9]);
    //aln->cigar = std::make_unique<uint16_t[]>(aln->n_cigar);
    aln->cigar = std::unique_ptr<uint16_t[]>(
        new uint16_t[aln->n_cigar],
        std::default_delete<uint16_t[]>()
    );
    aln->base = std::unique_ptr<uint64_t[]>(
        new uint64_t[aln->n_cigar],
        std::default_delete<uint64_t[]>()
    );
    //aln = std::make_shared<aln_ls_paf>(aln.n_cigar);
    int j=0;
    for (uint16_t i = 0; i < aln->n_cigar; ++i) {
        aln->cigar[i] = std::stoi(reg[i+j + 10]);
        if(MM_CIGAR_STR[aln->cigar[i]&0xf]=='I'||MM_CIGAR_STR[aln->cigar[i]&0xf]=='X'){
            j++;
            aln->base[i] = std::stoull(reg[i+j + 10]);
        }
        else{
            aln->base[i] =0;
        }
    }
}
void rManage::ssplitFasta(const std::string inputFilename,int anyFilesSeqNumber, int numFiles) {
    std::ifstream inputFile(inputFilename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }
    
    std::vector<std::ofstream> tmpFiles(numFiles);
    for (int i = 0; i < numFiles; ++i) { 
        std::string alnFilename=inputFilename;
        tmpFiles[i].open(alnFilename.insert(alnFilename.size()-3,"aln"+std::to_string(i)));
        if (!tmpFiles[i].is_open()) {
            std::cerr << "Error opening output file " << i << std::endl;
            return;
        }
        std::string line;
        int counters=anyFilesSeqNumber;
        int index=-1;
        while(--counters>=0&&std::getline(inputFile, line)){
            tmpFiles[i] << line << "\n";
            std::getline(inputFile, line);
            tmpFiles[i] << line << "\n";
        }
        tmpFiles[i].close();
    }
    inputFile.close();
    //std::filesystem::remove(inputFilename);
}
void rManage::splitFasta(const std::string& inputFilename, int numFiles ,string directoryPath) {
    countAlnSplitNumber();
    id_to_name=new string[rawSequnenceSplitNumber*FilesSeqNumber];
    std::ifstream inputFile(inputFilename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }
    std::vector<std::ofstream> tmpFiles(numFiles);
    for (int i = 0; i < numFiles; ++i) {
        tmpFiles[i].open(directoryPath+"tmp" + std::to_string(i) + ".fa");
        if (!tmpFiles[i].is_open()) {
            std::cerr << "Error opening output file " << i << std::endl;
            return;
        }
    }  
    std::string line;
    int fileIndex = 0;
    int index=0;
    string prefixL=">LR";
    string seqId=">LR"+std::to_string(0);
    while (std::getline(inputFile, line)) {
        
        if (line[0] == '>') {
            if (fileIndex == numFiles) {
                fileIndex = 0;
                index++;
                seqId=prefixL+std::to_string(index);
            }
            //cout<<line<<" "<<fileIndex<<" "<<FilesSeqNumber<<" "<<index<<" "<<endl;
            id_to_name[fileIndex*FilesSeqNumber+index]=line.substr(1);
            tmpFiles[fileIndex] <<seqId<<"\n";
            std::getline(inputFile, line);
            tmpFiles[fileIndex] << line << "\n";
            fileIndex++;
        } 
    
    }
    
    for (int i = 0; i < numFiles; ++i) {
        tmpFiles[i].close();
    }
    
    inputFile.close();
    
    
    for (int i = 0; i < numFiles; ++i) {
        ssplitFasta(directoryPath+"tmp"+std::to_string(i) + ".fa",anyFilesSeqNumber,alnSplitNumber);
    }
}
void splitAlnSortFile(const std::string inputFilename,int anyFilesSeqNumber) {
    std::ifstream inputFile(inputFilename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }
    std::ofstream tmpAlnFiles;
    int i=0;
    std::string alnFilename=inputFilename;
    int prefixLen=alnFilename.size()-4;
    tmpAlnFiles.open(alnFilename.insert(prefixLen,std::to_string(i)));
    if (!tmpAlnFiles.is_open()) {
        std::cerr << "Error opening output file " << i << std::endl;
        return;
    }
    std::string line;
    int head=0;
    int index=head;
    while(std::getline(inputFile, line)){
        index=std::stoi(line.substr(0, line.find(' ')));
        if(index-head>=anyFilesSeqNumber) {
            head=head+anyFilesSeqNumber;
            i++;
            tmpAlnFiles.close();
            tmpAlnFiles.open(alnFilename.replace(prefixLen,alnFilename.size()-prefixLen-4,std::to_string(i)));
        }
        tmpAlnFiles << line << "\n";
    }
    tmpAlnFiles.close();
}
void rManage::countLongReadsNumber(){
    BankFasta bsize(rawLongReadsFilePath);
    BankFasta::Iterator itSeqSize(bsize);
    size_t seqSize;
    for (itSeqSize.first(); !itSeqSize.isDone(); itSeqSize.next()) {
        seqSize = itSeqSize->getDataSize();
        if (seqSize > maxReadsLengh) {
            maxReadsLengh = seqSize;
        }
        longReadsNumber++;
    }
    if(longReadsNumber<alnSplitNumber){
        alnSplitNumber=longReadsNumber;
    }
    maxReadsLengh = maxReadsLengh * 1.25;
}
void rManage::generateTmpFaSplitFile(int i ,string &tmplongFile,string &alnFilename){
    tmplongFile=directoryPath+"tmp"+std::to_string(i) + ".fa";
    alnFilename=tmplongFile;
    alnFilename.replace(alnFilename.size()-2,2,"aln");
    cout<<"Target sequence path:"<<tmplongFile<<endl;
    cout<<"Query sequence path:"<<rawShortReadsFilePath<<endl;
    cout<<"Output file path:"<<alnFilename<<endl;
    runMinimap(tmplongFile.c_str(),rawShortReadsFilePath.c_str(),alnFilename);
    std::filesystem::remove(tmplongFile);
    splitAlnSortFile(alnFilename,anyFilesSeqNumber);
    alnFilename.insert(alnFilename.size()-4,"0");
    cout<<"sorted file path:"<<alnFilename<<endl;
    tmplongFile.insert(tmplongFile.size()-3,"aln0");
    
}


void rManage::findDirectory(){
    int pathpos = rawLongReadsFilePath.find_last_of("/");
    if (pathpos != std::string::npos) {
        directoryPath = rawLongReadsFilePath.substr(0, pathpos + 1);
    }
}
void rManage::countAlnSplitNumber(){
    FilesSeqNumber=longReadsNumber/rawSequnenceSplitNumber+1;
    if(FilesSeqNumber>seqThreshold){
        alnSplitNumber=FilesSeqNumber/seqThreshold+1;
        anyFilesSeqNumber=seqThreshold;   
    }
    else{
        anyFilesSeqNumber=FilesSeqNumber;
    }
}

