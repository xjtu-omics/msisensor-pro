CXX=g++

# HTSLIB_ROOT=$(realpath ../vendor/htslib-1.11)
# HTSLIB_CPPFLAGS=-I$(realpath $(HTSLIB_ROOT))
# echo $(HTSLIB_CPPFLAGS)
# HTSLIB_LDFLAGS=-L/mnt/project/Project_Pindel/pindel_update/pindel/htslib/ -Wl,-rpath
CPP_ROOT=cpp
CXXFLAGS+=    -O2 -fopenmp
LDFLAGS+=    -L${HTSLIB}
INCLUDES+=    -I${HTSLIB}
LIBS+=    -lm -lz -lpthread -lhts
LD_LIBRARY_PATH+= -L${HTSLIB}
export LD_LIBRARY_PATH=${HTSLIB_ROOT}

SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic  
OBJS= $(patsubst %,%.o,$(SOURCE))

%.o:%.cpp
	        $(CXX) -L ${HTSLIB} -std=c++11 $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all:  msisensor-pro
#
htsfile:
	        $(MAKE) -C ${HTSLIB}

msisensor-pro: $(OBJS)
	        $(CXX)  -std=c++11 $^ $(CXXFLAGS) $(LDFLAGS)  $(LIBS) -Wl,-rpath=${HTSLIB} -o $@
# 	        $(CXX)  -std=c++11 $^ $(CXXFLAGS) $(INCLUDES)  $(LDFLAGS) $(LIBS) -o $@

clean:
	        rm -f *.o msisensor-pro
			        $(MAKE) -C ${HTSLIB}  clean

