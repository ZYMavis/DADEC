#ifndef dbgGraphbuild_h
#define dbgGraphbuild_h
#include <filesystem>
#include <gatb/gatb_core.hpp>
typedef GraphUnitigsTemplate<128> GraphType;

void bcalmDbgGraphbuild(string illuminaFile,int kmerSize,int AbundanceMin,int cores){
    cout<<"--------------build graph-------------------------"<<endl;
	IBank *shortReads = Bank::open(illuminaFile);
    string solid_kmer_thr_str="1";
    string illuminaGraph = illuminaFile + "_k" + std::to_string(kmerSize) + ".h5";
    string illuminaFa = illuminaFile + "_k" + std::to_string(kmerSize) + "_a" + std::to_string(AbundanceMin) + ".fa";
	std::string arg="-kmer-size "+std::to_string(kmerSize)+" -out "+illuminaGraph+" -nb-cores "+std::to_string(cores)+" -abundance-min "+std::to_string(AbundanceMin);
    cout<<"arg:"<<arg<<endl;
    cout<<"--------------end-------------------------"<<endl;
   try
    {
       GraphType graph = GraphType::create (shortReads,arg.c_str());
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }
    filesystem::copy("dummy.unitigs.fa",illuminaFa);

}

#endif