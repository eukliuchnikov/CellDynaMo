# Makefile
# Generic Makefile for making cuda programs
#

# Makefile
# Generic Makefile for making cuda programs
#
BIN					:= cell
# flags
CUDA_INSTALL_PATH	:= /usr/local/cuda
OBJDIR				:= obj
INCLUDES			+= -I$(CUDA_INSTALL_PATH)/include -I.
LIBS				:= -L$(CUDA_INSTALL_PATH)/lib64
CFLAGS				:= -O0 -g
LDFLAGS				:= -lrt -lm -lcudart
# compilers
#NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc --compiler-options -fpermissive -arch sm_20 --ptxas-options=-v
NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc --compiler-options -fpermissive --ptxas-options=-v -use_fast_math 
CC					:= g++
LINKER				:= g++ -fPIC -g
# files
CPP_SOURCES			:= \
  src/main.cpp \
  src/TCO_Writers/dcdio.cpp \
  src/TCO_Writers/psfio.cpp \
  src/TCO_Writers/xyzio.cpp \
  src/Initializiation/init.cpp \
  src/Initializiation/kinesin.cpp \
  src/Initializiation/initSV.cpp \
  src/Initializiation/configreader.cpp \
  src/Initializiation/visualization.cpp \
  src/Initializiation/mt_init.cpp \
  src/Initializiation/kt_init.cpp \
  src/Initializiation/implicit.cpp \
  src/Initializiation/diff_update.cpp \
  src/Initializiation/parameters.cpp \
  src/Math/mat.cpp \
  src/Chemistry/react.cpp \
  src/Chemistry/diffusion.cpp \
  src/Chemistry/phos-dephos.cpp \
  src/Chemistry/attach-detach.cpp \
  src/Chemistry/mt_growth.cpp \
  src/Chemistry/mt_short.cpp \
  src/Chemistry/mt_cat.cpp \
  src/Chemistry/mt_resc.cpp \
  src/Chemistry/kinesin_react.cpp \
  src/Potentials/membrane.cpp

CU_SOURCES			:= \
  src/Potentials/langevin.cu \
  src/Potentials/pushing.cu \
  src/Potentials/pulling.cu \
  src/Potentials/harmonic.cu \
  src/Potentials/kt_ex_vol.cu \
  src/Potentials/cylinder.cu \
  src/Potentials/bending.cu
 
$(BIN):
	$(NVCC) -o $(BIN) $(CPP_SOURCES) $(CU_SOURCES) $(LDFLAGS) $(INCLUDES) $(LIBS)

run: $(BIN)
	LD_LIBRARY_PATH=$(CUDA_INSTALL_PATH)/lib ./$(BIN)
	
install:
	cp $(BIN) /usr/bin/$(BIN)
	
