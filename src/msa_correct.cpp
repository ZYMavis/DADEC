#include "msa_correct.h"



correctRead::correctRead(aln_paf_v ap,const Sequence& seq,int maxLength,int thr){
    if (seq.getDataSize() > maxLength) {
        std::cout << "Too long read" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::map<uint32_t,uint32_t> lp;
    alnPaf=ap;
    long_idx=seq.getIndex();
    read_len = seq.getDataSize();
    n_qReads=alnPaf.n();
    threshold=thr;
    lread.reset(new char[maxLength]);
    newLread.reset(new char[maxLength]);
    copy_upper_case(lread.get(), seq.getDataBuffer(), read_len);
    lread[read_len] = '\0';
}
int correctRead::msa(std::string LongReadsFile,std::string ShortReadsFile,std::string outputFile,int msa_thread,float threshold,int splitNumber){
    try{
        rManage alrm(LongReadsFile,ShortReadsFile,splitNumber,200000);
        cout << "Found " << alrm.longReadsNumber << " reads.\n";
        cout << "Correcting reads...\n"; 
        BankFasta output(outputFile);
         
        for (int i = 0; i < alrm.rawSequnenceSplitNumber; ++i) {
            string tmplongFile;
            string alnFilename;
            alrm.generateTmpFaSplitFile(i,tmplongFile,alnFilename);
            
            int prefixLen=tmplongFile.size()-4;
            int prefixLen2=alnFilename.size()-5;

            for(int j=0;j<alrm.alnSplitNumber;j++){
                tmplongFile.replace(prefixLen,tmplongFile.size()-prefixLen-3,std::to_string(j));
                alnFilename.replace(prefixLen2,alnFilename.size()-prefixLen2-4,std::to_string(j));
                cout<<"sequence path:"<<tmplongFile<<endl;
                cout<<"alignment path:"<<alnFilename<<endl;
                std::ifstream alnFile(alnFilename);
                std::vector<aln_paf_v> alnPafs(alrm.longReadsNumber + 2);
                std::string line;
                int index=-1;
                while(std::getline(alnFile, line)){
                    bool isvalid = true;
                    std::vector<string> reg;
                    std::istringstream iss(line); 
                    string num;
                    while (iss >> num) {
                        if (num.length()>8) {
                            isvalid = false;
                            break;
                        }
                        reg.push_back(num);
                    }
                    if(!isvalid) {
                        cout << "Invalid number in line: " << line << endl;
                        continue; // Skip this line
                    }
                    index=std::stoi(reg[0]);           
                    if(index/alrm.anyFilesSeqNumber>j) {
                        cout<<"index/anyFilesSeqNumber>j "<<index<<" "<<alrm.anyFilesSeqNumber<<" "<<j<<endl;
                        break;
                    }
                    std::shared_ptr<aln_ls_paf> map_detail = std::make_shared<aln_ls_paf>();
                    rManage::init_aln_ls_paf(map_detail, reg);
                    alnPafs[index % alrm.anyFilesSeqNumber].a.push_back(map_detail);
                }
                alnFile.close();
                IBank *tmpBankPB = NULL;
                try{
                    tmpBankPB =  Bank::open(tmplongFile);
                }
                catch (gatb::core::system::Exception & bankPBExc){
                    std::cout << "Error message PacBio bank " << bankPBExc.getMessage () << std::endl;
                    return EXIT_FAILURE;
                }
                Iterator<Sequence> *itSeq= tmpBankPB->iterator();

                ISynchronizer *sync = System::thread().newSynchronizer();
                
                IDispatcher::Status status = Dispatcher(64).iterate(itSeq, [&] (const Sequence& seq) {
                    int long_index=seq.getIndex();
                    correctRead cread(alnPafs[long_index],seq,alrm.maxReadsLengh,threshold);
                    cread.prepare_snp();
                    int id=std::stoi(seq.getComment().substr(2));

                    Sequence s(cread.newLread.get());
                    s._comment = alrm.id_to_name[i*alrm.FilesSeqNumber+id];
                    {
                        LocalSynchronizer local(sync);
                        output.insert(s);
                    }
                });
                delete tmpBankPB;
                std::filesystem::remove(alnFilename);
            }
        
        }
        output.flush();
        return EXIT_SUCCESS;
    } catch (gatb::core::system::Exception & e){
        std::cout << "Error message " << e.getMessage () << std::endl;
        return EXIT_FAILURE;
    }
}
  
void correctRead::prepare_snp(){

    for(int i=0;i<n_qReads;i++){ 
        aln_ls_paf& ls_detail=*alnPaf.a[i];
        uint16_t* cigar=ls_detail.cigar.get();
        uint64_t* base=ls_detail.base.get();
        uint32_t srPos=0;
        uint32_t pos=ls_detail.rs;
        uint32_t pos2=ls_detail.re;

        auto search=map.find(pos);
        if(search==map.end()){
            map.emplace(pos, std::initializer_list<uint32_t>{pos2});
        }
        else{
            map[pos].push_back(pos2);
        }
        update_snp(base,ls_detail.n_cigar,cigar,pos,srPos,ls_detail.index,ls_detail.qt,1);
    }


    generate_mut();
    int nomutAlin=0;
    for(int i=0;i<n_qReads;i++){ 
        if(alnPaf.a[i]->mut==1) continue;
        aln_ls_paf& ls_detail=*alnPaf.a[i];
        uint16_t alnScore=ls_detail.as;
        uint16_t* cigar=ls_detail.cigar.get();
        uint64_t* base=ls_detail.base.get();
        uint32_t srPos=0;
        uint32_t pos=ls_detail.rs;
        uint32_t pos2=ls_detail.re;
        bool is_mut = find_mut( ls_detail.index, ls_detail.matchcount);
        if (is_mut) {

            nomutAlin++;
            auto it = std::find(map[pos].begin(), map[pos].end(), pos2);
            map[pos].erase(it);

            if(map[pos].empty()){
                map.erase(map.find(pos));
            }

            for(uint16_t j=0;j<ls_detail.n_cigar;++j){
                uint16_t snp_len=cigar[j]>>4;
                char type=MM_CIGAR_STR[cigar[j]&0xf];
                cigar_ite(base[j],snp_len,type,pos,srPos,ls_detail.index,ls_detail.qt,2);
            }     

        }
    }
    //cout<<long_idx<<" all:"<<n_qReads<<" filt:"<<nomutAlin<<endl;
    //cout<<"test1"<<endl;

    update_lread();
   
    pos_snp.clear();
    map.clear();
    mutReads.clear();

    }

bool correctRead::find_mut(int sidx,int matchcount){
    auto search=mutReads.find(sidx);
    if(search!=mutReads.end()){
        float num=search->second;
        float fre=num/matchcount;
        if(fre>threshold)
            return true;
    }
    return false;
}


bool correctRead::haveKey(std::map<string,vector<uint32_t>> &b,char* key,uint32_t sid,int mode){
    string base=key;
    auto search=b.find(base);
    if(mode==1){
        if(search==b.end()){
            b.emplace(base, std::initializer_list<uint32_t>{sid});
            return false;  
        }
        search->second.push_back(sid);
    }
    else{
        auto& vec = search->second;  
        auto key = std::find(vec.begin(), vec.end(), sid);
        if (key != vec.end()) {
            vec.erase(key);  

            if (vec.empty()) {
                b.erase(search); 
                return true;
            }
            return true;  
        }
    }
    return false; 
}


bool correctRead::add_mut(int sid){
    auto search=mutReads.find(sid);
    if(search==mutReads.end()){
        mutReads.insert({sid,1});
        return false;
    }
    search->second=search->second+1;
    return true;  
}


void correctRead::cigar_ite(uint64_t baseCode,uint16_t snp_len,char type,uint32_t &pos,uint32_t &srPos,uint32_t sid,char qt,int ite){
         
    if(type=='=')
    {
        pos=pos+snp_len;
        srPos=srPos+snp_len;
        return ;
    }

    std::unique_ptr<char[]> snp_base(new char[snp_len + 1]);
    
    for(int i=snp_len-1;i>=0;i--){
        if(type=='I'||type=='X'){
            snp_base[i]=BASE[baseCode&0x3];
            baseCode=baseCode>>2;
        }
        else{
            snp_base[i]=lread[srPos+i];
        }
    }  
    
    snp_base[snp_len]='\0'; 
    bool eras=false; 

    if(type=='I'){
        eras=haveKey(pos_snp[pos].i,snp_base.get(),sid,ite);
        srPos=srPos+snp_len;
    }
    else if(type=='D'){
        eras=haveKey(pos_snp[pos].d,snp_base.get(),sid,ite);
        pos=pos+snp_len;
    }
    else{
        eras=haveKey(pos_snp[pos].x,snp_base.get(),sid,ite);
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

void correctRead::find_max(std::map<string,vector<uint32_t>> &b,int &maxdepth,int &type,string &base,int &rsnp_len,int t,int pos,int &num){
    if(b.empty()) return;
    
    for (auto it = b.begin(); it != b.end(); ++it) {
    
        string snp_base=it->first;
        uint16_t snp_len=snp_base.length();
        vector<uint32_t> snp_reads=it->second;
        int depth=snp_reads.size();
        num+=depth;
        std::unique_ptr<char[]> lbase(new char[snp_len + 1]);

        for(uint16_t i=0;i<snp_len;i++)
            lbase[i]=lread[pos+i];
        lbase[snp_len]='\0';  
        if(depth<maxdepth) continue;
        else if((depth==maxdepth&&(type==2||type==0)&&snp_base.c_str()==lbase.get())||depth>maxdepth){
            maxdepth=depth;
            type=t;
            base=lbase.get();
            rsnp_len=snp_len;
        }
    }

}

uint16_t correctRead::countContainsPos(uint32_t pos){
    uint16_t map_num=0;
    for(auto map_ite=map.begin();map_ite!=map.end();++map_ite){
        uint32_t start=map_ite->first;
        if(pos<start){break;}
        else{
            vector<uint32_t> end=map_ite->second;

            for(int i=end.size()-1;i>=0;i--){
                if(pos<=end[i]) {
                    map_num++;
                }
            }
        }
    }
    return map_num;
}
void correctRead::generate_mut(){

    for(int pos=0;pos<read_len;pos++){
        int type=0;
        string base;
        int maxdepth=0,rsnp_len;
        int snp_num=0;
        float c=0.0;
   
        if(pos_snp.find(pos)!=pos_snp.end()){
            std::map<string,vector<uint32_t>> mis_snp=pos_snp[pos].x;
            find_max(mis_snp,maxdepth,type,base,rsnp_len,2,pos,snp_num);  

            int map_num=countContainsPos(pos);
            c=float(snp_num)/float(map_num);  
            if(c<=0.5 or (map_num-snp_num)>5){
                std::map<string,vector<uint32_t>> mis_snp=pos_snp[pos].x;
                for(auto mis=mis_snp.begin();mis!=mis_snp.end();++mis){
                    vector<uint32_t> snp_reads=mis->second;
                    float snp_ratio=float(snp_reads.size())/float(map_num);
                    if(snp_reads.size()>5||snp_ratio>=0.5){
                        for(int p=0;p<snp_reads.size();p++){
                            add_mut(snp_reads[p]);
                        }
                    }
                }        
            } 
        }
    }

}
void correctRead::update_lread(){
    int newpos=0;

    for(uint32_t pos=0;pos<read_len;pos++){
        int type=0;
        string base;
        int maxdepth=0,rsnp_len;
        int snp_num=0;
        int c=0;
 
        if(pos_snp.find(pos)!=pos_snp.end()){
            find_max(pos_snp[pos].x,maxdepth,type,base,rsnp_len,2,pos,snp_num);  
            find_max(pos_snp[pos].i,maxdepth,type,base,rsnp_len,0,pos,snp_num);
            find_max(pos_snp[pos].d,maxdepth,type,base,rsnp_len,1,pos,snp_num);
            
            int map_num=countContainsPos(pos);
            c=map_num-snp_num;  
        }
    
        if(c>=5||maxdepth<5){
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
void correctRead::update_snp(uint64_t* base,int lcigar,uint16_t* cigar,uint32_t &pos,uint32_t &srPos,int sid,char qt,int ite){
    
    for(int j=0;j<lcigar;++j){
        int snp_len=cigar[j]>>4;
        char type=MM_CIGAR_STR[cigar[j]&0xf];
        if(type!='='){
            if(pos_snp.find(pos)==pos_snp.end()){
                ssnp all_snp={{},{},{}};
                pos_snp.insert({pos,all_snp});
            }
        }
        cigar_ite(base[j],snp_len,type,pos,srPos,sid,qt,ite);
    }

}
void correctRead::copy_upper_case(char *target, char *source, int len) {
    uint32_t tag=0,start=0,end=0;
    for(uint32_t i = 0; i < len; i++) {
      if(tag==0&&islower(source[i])){
        start=i;
        tag=1;
      }
      else if(tag==1&&isupper(source[i])){
        end=i-1;
        tag=0;
      }
      target[i] = toupper(source[i]);
    }
}
void correctRead::reverse(char *target, char *source, int len) {
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
