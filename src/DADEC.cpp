#include <iostream>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include "cxxopts.hpp"
#include "stream.hpp"
#include "ThreadReadAssertion.h"
#include "EValue.h"
#include "msa_correct.h"
#include <gatb/gatb_core.hpp>
#include "Aligner.h"


struct CorrectParams
{
	std::string shortFile;
	std::string longFile;
	std::string out;
	size_t numThreads;
	int kmerSize1;
	int kmerSize2;
	float msaThreshold;
	int splitNumber;
	int am1;
	int am2;
};
string VERSION="1.0";

void dbgCorrect(CorrectParams correctParams,int step){
    int kmerSize;
    int am;
    if(step==1){
        kmerSize=correctParams.kmerSize1;
        am=correctParams.am1;
    }
    else{
        kmerSize=correctParams.kmerSize2;
        am=correctParams.am2;
    }
	AlignerParams params;
	params.graphFile = correctParams.shortFile;
	params.fastqFiles={correctParams.longFile};
	params.outputGAMFile = "";
	params.outputJSONFile = "";
	params.outputGAFFile = "";
	params.outputCorrectedFile = correctParams.out;
	params.outputCorrectedClippedFile = "";
	params.numThreads = correctParams.numThreads;
	params.alignmentBandwidth = 5;
	params.dynamicRowStart = false;
	params.maxCellsPerSlice = 10000;
	params.verboseMode = false;
	params.mxmLength = 20;
	params.mumCount = 0;
	params.memCount = 0;
	params.seederCachePrefix = "";
	params.selectionECutoff = -1;
	params.compressCorrected = false;
	params.compressClipped = false;
	params.minimizerSeedDensity = 5;
	params.minimizerLength = 19;
	params.minimizerWindowSize = 30;
	params.seedClusterMinSize = 1;
	params.minimizerDiscardMostNumerousFraction = 0.0002;
	params.maxClusterExtend = 5;
	params.preciseClippingIdentityCutoff = 0.66;
	params.Xdropcutoff = 50;
	params.DPRestartStride = 0;
	params.multimapScoreFraction = 0.9;
	params.cigarMatchMismatchMerge = false;
	params.minAlignmentScore = 0;
	params.hpcCollapse = false;
	params.includeCigar = true;
	params.clipAmbiguousEnds = -1;
	params.maxTraceCount = 10;
	params.overlapIncompatibleCutoff = 0.3;
	params.realignFile = "";
	params.uniqueMemBonusFactor = 1.0;
	params.lowMemoryMEMIndexConstruction = false;
	params.MEMindexUsesWaveletTree = true;
	params.MEMwindowsize = 0;
	params.useDiploidHeuristic = false;
	params.diploidHeuristicCacheFile = "";
	params.keepSequenceNameTags = false;
	params.kmerSize = kmerSize;
	params.mode = am;
	params.dbgBuildMethod="bcalm2";
    alignReads(params);
	std::filesystem::remove("dummy.unitigs.fa");
}

void msaCorrect(CorrectParams correctParams) {
    try{
		correctRead::msa(correctParams.longFile,correctParams.shortFile,correctParams.out,correctParams.numThreads,correctParams.msaThreshold,correctParams.splitNumber);
	} catch (Exception & e){
        std::cout << "Error message " << e.getMessage () << std::endl;
  	}
}


int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

#ifndef NOBUILTINPOPCOUNT
	if (__builtin_cpu_supports("popcnt") == 0)
	{
		std::cerr << "CPU does not support builtin popcount operation" << std::endl;
		std::cerr << "recompile with -DNOBUILTINPOPCOUNT" << std::endl;
		std::abort();
	}
#endif

	struct sigaction act;
	act.sa_handler = ThreadReadAssertion::signal;
	sigemptyset(&act.sa_mask);
	act.sa_flags = 0;
	sigaction(SIGSEGV, &act, 0);

	cxxopts::Options options { "DADEC" };
	options.add_options("Mandatory parameters")
		("short-reads,s", "input short reads (fasta, uncompressed or gzipped)", cxxopts::value<std::string>())
		("long-reads,l", "input reads (fasta, uncompressed or gzipped)", cxxopts::value<std::string>())
		("out,o", "output corrected reads file (.fa)", cxxopts::value<std::string>())
	;
	options.add_options("General parameters")
		("help,h", "help message")
		("version,v", "print version")
		("threads,t", "number of threads (int) (default 1)", cxxopts::value<size_t>())
		("splitNumber,S",  "number of splits for long-read files", cxxopts::value<int>())
		;
	options.add_options("Build De Bruijn Graph")
	    ("kmer-size1,k", "k-mer size for the first DBG-based correction", cxxopts::value<int>())
		("kmer-size2,K", "k-mer size for the second DBG-based correction", cxxopts::value<int>())
		("msa-threshold,r", "Haplotype filtering threshold", cxxopts::value<float>())
		("abundance-min1,a", "Abundance threshold for the first DBG construction", cxxopts::value<int>())
		("abundance-min2,A",  "Abundance threshold for the second DBG construction", cxxopts::value<int>())
		;

	auto vm = options.parse(argc, argv);

	if (vm.count("help"))
	{
		std::cerr << options.help({"Mandatory parameters",  "General parameters", "Build De Bruijn Graph"}) << std::endl;
		std::cerr << "Preset parameters" << std::endl;
		std::exit(0);
	}
	if (vm.count("version"))
	{
		std::cout << "Version " << VERSION << std::endl;
		std::exit(0);
	}

	CorrectParams params;
	params.shortFile = "";
	params.longFile = "";
	params.out = "DADEC.fa";
	params.numThreads = 1;
	params.kmerSize1 = 39;
    params.kmerSize2 = 39;
	params.msaThreshold=0.008;
	params.splitNumber=10;
	params.am1 = 2;
	params.am2 = 1;

	try {
		if (vm.count("short-reads")) params.shortFile = vm["short-reads"].as<std::string>();
		if (vm.count("long-reads")) params.longFile = vm["long-reads"].as<std::string>();
		if (vm.count("out")) params.out = vm["out"].as<std::string>();
		if (vm.count("threads")) params.numThreads = vm["threads"].as<size_t>();
		if (vm.count("kmer-size1")) params.kmerSize1 = vm["kmer-size1"].as<int>();
		if (vm.count("kmer-size2")) params.kmerSize2 = vm["kmer-size2"].as<int>();
		if (vm.count("msa-threshold")) params.msaThreshold = vm["msa-threshold"].as<float>();
		if (vm.count("abundance-min1")) params.am1 = vm["abundance-min1"].as<int>();
		if (vm.count("abundance-min2")) params.am2 = vm["abundance-min2"].as<int>();
		if (vm.count("splitNumber")) params.splitNumber = vm["splitNumber"].as<int>();
	}
	catch (const cxxopts::exceptions::exception& e) {
		std::cerr << "parameters wrong: " << e.what() << std::endl;
		return 1;
	}
	
	if (params.shortFile == "")
	{
		std::cerr << "short reads file must be given" << std::endl;
	}
	if (params.longFile == "")
	{
		std::cerr << "long read file must be given" << std::endl;
	}
	if (params.numThreads < 1)
	{
		std::cerr << "number of threads must be >= 1" << std::endl;
	}


    string output=params.out;

    params.out="tmp_DADEC1.fa";
	cout<<"the first step with High-Confidence Error Elimination"<<endl;
	dbgCorrect(params,1);
	std::filesystem::remove("dummy.unitigs.fa");
	cout<<"-----------------end---------------------"<<endl;

	
	cout<<"the second step with Haplotype-Specific Refinement"<<endl;
    params.longFile=params.out;
    params.out="tmp_DADEC2.fa";
    msaCorrect(params);
	std::filesystem::remove(params.longFile);
	
	cout<<"-----------------end---------------------"<<endl;

	cout<<"the third step with Low-Abundance Information Recovery"<<endl;
	params.longFile=params.out;
    params.out=output;
    dbgCorrect(params,2);
	std::filesystem::remove(params.longFile);
	std::filesystem::remove("tmp*");
	cout<<"-----------------end---------------------"<<endl;


	return 0;
}

