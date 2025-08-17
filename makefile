# =====================
# 目录和路径配置
# =====================
GATB = gatb-core/gatb-core
MINIMAP = minimap2
GRAPHALIGNER = GraphAligner
ODIR =obj
BINDIR = bin
SRCDIR =src
GATB_LIB = $(GATB)/build/lib

PLATFORM = $(shell uname -s)
GPP=$(CXX)
MKDIR_P = mkdir -p

CXXFLAGS = -Wall -no-pie -Wextra -std=c++17 -O2 -g \
           -Wunused-function -Wignored-qualifiers -Wno-unused-parameter \
           -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE  

CPPFLAGS = `pkg-config --cflags protobuf libsparsehash zlib` \
           -Iconcurrentqueue -IBBHash -Izstr/src \
           -Iparallel-hashmap/parallel_hashmap/ \
           -I$(MINIMAP) -I$(GATB)/src \
           -I$(GATB)/build/include -I$(GATB)/thirdparty \
           -IMEMfinder/src -I`jemalloc-config --includedir` \
           -Icxxopts/include -I$(GRAPHALIGNER) -I$(SRCDIR) \
           -I`pkg-config --variable=includedir zlib`/bamtools

LDFLAGS = -L$(GATB_LIB) -L$(MINIMAP)
LDLIBS = -lgatbcore -lhdf5 -lminimap2 -ldl -lz -lm -lrt \
         `pkg-config --libs protobuf` -lsdsl -lbamtools \
         `pkg-config --libs libdivsufsort` `pkg-config --libs libdivsufsort64`

JEMALLOCFLAGS = -L`jemalloc-config --libdir` \
                -Wl,-rpath,`jemalloc-config --libdir` \
                -ljemalloc `jemalloc-config --libs`

LINKFLAGS = $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $(JEMALLOCFLAGS) -pthread

_DADECDEPS =  dbgGraphbuild.h msa_correct.h
DADECDEPS = $(patsubst %, $(SRCDIR)/%, $(_DADECDEPS))

_DADECOBJ =  msa_correct.o map.o alignManager.o
DADECOBJ = $(patsubst %, $(ODIR)/%, $(_DADECOBJ))


_DEPS = vg.pb.h GraphAlignerWrapper.h fastqloader.h \
        BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h \
        AlignmentGraph.h CommonUtils.h GfaGraph.h ReadCorrection.h \
        MinimizerSeeder.h AlignmentSelection.h EValue.h MEMSeeder.h \
        DNAString.h DiploidHeuristic.h
DEPS = $(patsubst %, $(GRAPHALIGNER)/%, $(_DEPS))


_OBJ = Aligner.o vg.pb.o  fastqloader.o\
       BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o \
       GraphAlignerWrapper.o GfaGraph.o ReadCorrection.o MinimizerSeeder.o \
       AlignmentSelection.o EValue.o MEMSeeder.o DNAString.o DiploidHeuristic.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) \
           commit $(shell git rev-parse HEAD) \
           $(shell git show -s --format=%ci)

$(shell $(MKDIR_P) $(BINDIR) $(ODIR))

DADEC: $(ODIR)/DADECMain.o $(OBJ) $(DADECOBJ) MEMfinder/lib/memfinder.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/GraphAlignerWrapper.o: $(GRAPHALIGNER)/GraphAlignerWrapper.cpp \
    $(GRAPHALIGNER)/GraphAligner.h $(GRAPHALIGNER)/NodeSlice.h $(GRAPHALIGNER)/WordSlice.h \
    $(GRAPHALIGNER)/ArrayPriorityQueue.h $(GRAPHALIGNER)/ComponentPriorityQueue.h \
    $(GRAPHALIGNER)/GraphAlignerVGAlignment.h $(GRAPHALIGNER)/GraphAlignerGAFAlignment.h \
    $(GRAPHALIGNER)/GraphAlignerBitvectorBanded.h \
    $(GRAPHALIGNER)/GraphAlignerBitvectorCommon.h $(GRAPHALIGNER)/GraphAlignerCommon.h \
    $(DEPS)

$(ODIR)/DADECMain.o: $(SRCDIR)/DADEC.cpp $(DADECDEPS) $(DEPS) $(OBJ) $(DADECOBJ) 
	$(GPP) -c -o $@ $< $(CPPFLAGS) $(CXXFLAGS)

$(ODIR)/%.o:$(SRCDIR)/%.cpp $(DADECDEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) $(CXXFLAGS)

$(ODIR)/%.o: $(GRAPHALIGNER)/%.cpp  $(DEPS) 
	$(GPP) -c -o $@ $< $(CPPFLAGS) $(CXXFLAGS)

$(ODIR)/vg.pb.o: $(GRAPHALIGNER)/vg.pb.cc
	$(GPP) -c -o $@ $< $(CPPFLAGS) $(CXXFLAGS)

$(GRAPHALIGNER)/%.pb.cc $(GRAPHALIGNER)/%.pb.h: $(GRAPHALIGNER)/%.proto
	protoc -I=$(GRAPHALIGNER) --cpp_out=$(GRAPHALIGNER) $<


MEMfinder/lib/memfinder.a:
	$(MAKE) -C MEMfinder lib DEBUGFLAG="-DNDEBUG"

.PHONY: all clean

build_gatb: 
	@echo "Building GATB-core..."
	@cd $(GATB) && mkdir -p build && cd build && cmake .. && make -j8 && cd .. && cd .. && cd ..

build_minimap2: 
	@echo "Building minimap2..."
	@cd $(MINIMAP) && make && cd ..

deps:build_gatb build_minimap2

all:deps DADEC 

clean:
	rm -rf $(ODIR) DADEC $(GRAPHALIGNER)/vg.pb.*
	$(MAKE) -C MEMfinder clean