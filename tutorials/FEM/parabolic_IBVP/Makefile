MODEL_INCLUDE =../../../
EIGEN_INCLUDE =/usr/local/include
usePARDISO = 0
MKL_PATH=/opt/intel/composer_xe_2013.3.171/mkl
##############################
IDIR += -I ./
IDIR += -I $(MODEL_INCLUDE)
IDIR += -I $(EIGEN_INCLUDE)

##############################
CFLAGS += -std=c++14  
CFLAGS += -O3 
CFLAGS += -fopenmp
CFLAGS += -Wall
CFLAGS += -Wextra

##############################
LIBS += 

##############################
# Settings for PARDISO solver
ifeq ($(usePARDISO), 1)
CFLAGS += -D _MODEL_PARDISO_SOLVER_
IDIR += -I $(MKL_PATH)/include
LIBS += -L $(MKL_PATH)/lib
LIBS += -lmkl_intel_lp64 
LIBS += -lmkl_intel_thread 
LIBS += -lmkl_core
endif


program:
	g++ main.cpp -o ibvp $(CFLAGS) $(IDIR)  $(LIBS)
	
clean:
		rm -f ibvp
		
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.132/mkl/lib:$DYLD_LIBRARY_PATH
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.132/compiler/lib:$DYLD_LIBRARY_PATH

# Hoffman2
#export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
