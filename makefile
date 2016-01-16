# %I% %G%
# name of the project (executable)
#PROJ = current-0001
PROJ = Curxy
#PROJ = densityDens1
#CONFIG = ifort-mic
CONFIG = ifort-opt

ifeq ($(CONFIG),ifort-mic)
SUFOBJ = .o
SUFMOD = .mod
SUFEXE = .exe
DASHC = -c #
DASHO = -o #
DASHE = -o #
F90 = ifort
F90FLAGS += -free -openmp -mmic
LD = ifort
LDFLAGS += -static-intel -mmic -mkl
LDLIBS =
RM = rm -f
endif

ifeq ($(CONFIG),ifort-opt)
SUFOBJ = .o
SUFMOD = .mod
SUFEXE = .exe
DASHC = -c #
DASHO = -o #
DASHE = -o #
F90 = ifort
F90FLAGS += -free -xhost -openmp -heap-arrays
LD = ifort
LDFLAGS += -static-intel -mkl
LDLIBS =
RM = rm -f
endif

ifeq ($(CONFIG),ifl-debug)
SUFOBJ = .obj
SUFMOD = .mod
SUFEXE = .exe
DASHC = /c #
DASHO = /Fo#
DASHE = /Fe#
F90 = ifl
F90FLAGS += /Zi /FR /Od /Zd
LD = ifl
LDFLAGS +=
#LDLIBS += -llapack
#LDLIBS += -lblas
LDLIBS += mkl_c_dll.lib
RM = rm -f
endif

ifeq ($(CONFIG),ifl-opt)
SUFOBJ = .obj
SUFMOD = .mod
SUFEXE = .exe
DASHC = /c #
DASHO = /Fo#
DASHE = /Fe#
F90 = ifort
F90FLAGS += /Ox /FR
LD = ifort
LDFLAGS +=
#LDLIBS += -llapack
#LDLIBS += -lblas
LDLIBS += mkl_intel_c_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib
RM = rm -f
endif

ifeq ($(CONFIG),zahir-debug)
SUFOBJ = .o
SUFMOD = .mod
SUFEXE =
DASHC = -c #
DASHO = -o #
DASHE = -o #
#90 = mpxlf90_r
F90FLAGS += -g
LD = mpxlf90_r
LDFLAGS +=
#LDLIBS += -llapack
#LDLIBS += -lblas
LDLIBS += -lessl
RM = rm -f
endif

ifeq ($(CONFIG),zahir-opt)
SUFOBJ = .o
SUFMOD = .mod
SUFEXE =
DASHC = -c #
DASHO = -o #
DASHE = -o #
F90 = mpxlf90_r
F90FLAGS += -O4
LD = mpxlf90_r
LDFLAGS += -qipa
#LDLIBS += -llapack
#LDLIBS += -lblas
LDLIBS += -lessl
RM = rm -f
endif

.SUFFIXES: .f .F $(SUFOBJ) $(SUFMOD) $(SUFEXE)

%$(SUFOBJ) : %.f
	$(F90) $< $(DASHC) $(F90FLAGS)

#%$(SUFOBJ) : %.f90
#	$(F90) $< $(DASHC) $(F90FLAGS)

objs += params$(SUFOBJ)
mods += params$(SUFMOD)
objs += potential2d$(SUFOBJ)
mods += potential2d$(SUFMOD)
objs += currLattice$(SUFOBJ)
#objs += densityDens$(SUFOBJ)
#objs += currEQUIL$(SUFOBJ)
#objs += mainCONDE$(SUFOBJ)
#objs += comp_evLAPACK$(SUFOBJ)
objs += greenL$(SUFOBJ)
objs += greenR$(SUFOBJ)
objs += g11lead$(SUFOBJ)
#objs += mysort$(SUFOBJ)
#objs += limits$(SUFOBJ)
#objs += green$(SUFOBJ)
objs += greenTotal$(SUFOBJ)
#objs += zgeev$(SUFOBJ)
#zgeev$(SUFOBJ): F90FLAGS = /Zi /Od /Zd

exes += $(PROJ)$(SUFEXE)
all : $(exes)

$(PROJ)$(SUFEXE) : $(objs)
	$(LD) $(objs) $(F90FLAGS) $(LDFLAGS) $(DASHE)$@ $(LDLIBS)

clean:
	-$(RM) $(objs)
	-$(RM) $(mods)
	-$(RM) $(exes)



# dependences
potential2d$(SUFOBJ) : params$(SUFOBJ)
#mainCONDE$(SUFOBJ) : params$(SUFOBJ) greenL$(SUFOBJ) greenR$(SUFOBJ) greenTotal$(SUFOBJ) g11lead$(SUFOBJ) potential2d$(SUFOBJ) mysort$(SUFOBJ) limits$(SUFOBJ) green$(SUFOBJ)
#currEQUIL$(SUFOBJ) : params$(SUFOBJ) greenL$(SUFOBJ) greenR$(SUFOBJ) greenTotal$(SUFOBJ) g11lead$(SUFOBJ) potential2d$(SUFOBJ)
#densityDens$(SUFOBJ) : params$(SUFOBJ) greenL$(SUFOBJ) greenR$(SUFOBJ) greenTotal$(SUFOBJ) g11lead$(SUFOBJ) potential2d$(SUFOBJ)
currLattice$(SUFOBJ) : params$(SUFOBJ) greenL$(SUFOBJ) greenR$(SUFOBJ) greenTotal$(SUFOBJ) g11lead$(SUFOBJ) potential2d$(SUFOBJ)
#comp_evLAPACK$(SUFOBJ) : params$(SUFOBJ)
#green$(SUFOBJ) : params$(SUFOBJ)
g11lead$(SUFOBJ) : params$(SUFOBJ)
greenR$(SUFOBJ) : params$(SUFOBJ)
greenL$(SUFOBJ) : params$(SUFOBJ)
greenTotal$(SUFOBJ) : params$(SUFOBJ)
