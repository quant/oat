win32 {
 INCLUDEPATH += "$$(MKLROOT)/include"
 LIBS += -L"$$(MKLROOT)\ia32\lib"
 LIBS += -L"$$(MKLROOT)\..\lib\ia32"
 LIBS += mkl_intel_c_dll.lib
 LIBS += mkl_intel_thread_dll.lib
 LIBS += mkl_core_dll.lib
 LIBS += libiomp5md.lib
 INCLUDEPATH += "C:\Program Files\MPICH2\include"
 LIBS += "C:\Program Files\MPICH2\lib\mpi.lib"
}

unix {
  QMAKE_CC = mpicc -cc=icc
  QMAKE_CXX = mpicxx -cxx=icpc
  QMAKE_LINK = mpicxx -cxx=icpc
  LIBS += -lmkl_intel_lp64
  #LIBS += -lmkl_intel_thread
  LIBS += -lmkl_sequential
  LIBS += -lmkl_core
  LIBS += -liomp5
  LIBS += -lpthread
  QMAKE_CXXFLAGS += -xSSE4.2
  QMAKE_CXXFLAGS += -Wall -w2 -Werror
  QMAKE_CXXFLAGS += -wd2557 # comparison between signed and unsigned operands
  QMAKE_CXXFLAGS += -wd981 # operands are evaluated in unspecified order
  QMAKE_CXXFLAGS += -wd383 # value copied to temporary, reference to temporary used
  QMAKE_CXXFLAGS += -wd869 # parameter "ex" was never referenced
  QMAKE_CXXFLAGS += -wd2407 # the initialization of member "MainWindow::sigmaU" will be done before that of member "MainWindow::typeCond"
  QMAKE_CXXFLAGS += -wd177 # variable "AA" was declared but never referenced
  QMAKE_CXXFLAGS += -wd2415 # variable "NCUT" of static storage duration was declared but never referenced
  #CONFIG += debug
  #QMAKE_CXXFLAGS += -g -O0
  QMAKE_LFLAGS = -Wl,-Map,pdv.map
}

# Local variables:
# mode: makefile
# End:
