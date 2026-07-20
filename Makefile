
UNAME_S       := $(shell uname -s)
GSL_PREFIX    ?= $(shell if command -v brew >/dev/null 2>&1; then brew --prefix gsl 2>/dev/null; fi)
OPENMP        ?= 0
OMP_PREFIX    ?= $(shell if command -v brew >/dev/null 2>&1; then brew --prefix libomp 2>/dev/null; fi)

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -std=c++17
LD            = $(CXX)
LDFLAGS       = -O3

LIBS          = $(SYSLIBS) -lgsl -lgslcblas

ifneq ($(GSL_PREFIX),)
CXXFLAGS      += -I$(GSL_PREFIX)/include
LDFLAGS       += -L$(GSL_PREFIX)/lib
endif

ifneq ($(UNAME_S),Darwin)
THREAD_FLAGS  ?= -pthread
CXXFLAGS      += $(THREAD_FLAGS)
LDFLAGS       += $(THREAD_FLAGS)
endif

ifeq ($(OPENMP),1)
ifeq ($(UNAME_S),Darwin)
OMP_CXXFLAGS  ?= -Xpreprocessor -fopenmp
ifneq ($(OMP_PREFIX),)
OMP_CXXFLAGS  += -I$(OMP_PREFIX)/include
OMP_LDFLAGS   ?= -L$(OMP_PREFIX)/lib
endif
OMP_LIBS      ?= -lomp
else
OMP_CXXFLAGS  ?= -fopenmp
OMP_LDFLAGS   ?=
OMP_LIBS      ?= -fopenmp
endif
CXXFLAGS      += $(OMP_CXXFLAGS)
LDFLAGS       += $(OMP_LDFLAGS)
LIBS          += $(OMP_LIBS)
endif

vpath %.cpp src
ifeq ($(OPENMP),1)
objdir     = obj_omp
else
objdir     = obj
endif

SRC        = cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoCMF.cpp eoCMFe.cpp eoHadron.cpp eoAZH.cpp eoSmash.cpp \
			 trancoeff.cpp fld.cpp hdo.cpp s95p.cpp icurqmd.cpp ic.cpp ickw.cpp icPartUrqmd.cpp icPartSMASH.cpp \
			 icDynFlu.cpp main.cpp rmn.cpp cornelius.cpp icGlauber.cpp icGubser.cpp icGlissando.cpp icTrento.cpp \
			 icTrento3d.cpp icSuperMC.cpp vtk.cpp icTest.cpp particle.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = hlle_visc
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -rf obj obj_omp $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)

$(objdir)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
