CXX=g++
SAMTOOLS_ROOT=../vendor/samtools-0.1.19
CPP_ROOT=cpp

CXXFLAGS+=    -O2 -fopenmp
LDFLAGS+=    -L${SAMTOOLS_ROOT}
LIBS+=    -lm -lbam -lz -lpthread
INCLUDES+=    -I${SAMTOOLS_ROOT}
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic  
OBJS= $(patsubst %,%.o,$(SOURCE))

%.o:%.cpp
	        $(CXX) -std=c++11 $(CXXFLAGS) $(INCLUDES) -c $< -o $@ 

all: samtools msisensor-pro

samtools:
	        $(MAKE) -C ${SAMTOOLS_ROOT}

msisensor-pro: $(OBJS)
	        $(CXX) -std=c++11 $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ 

clean:
	        rm -f *.o msisensor-pro
			        $(MAKE) -C ${SAMTOOLS_ROOT}  clean

