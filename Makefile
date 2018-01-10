PROJ = pdv2016c
MPICXX = mpiicpc
CXXFLAGS = -mkl

all: $(PROJ).exe
	qsub run.pbs

OBJ = mpidistributions1.o mpimainwindow.o percol2d.o PercolRectSides.o
HDR = $(wildcard *.h)

$(PROJ).exe: $(OBJ)
	$(MPICXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ): $(HDR)

clean:
	-rm *.o *.exe
