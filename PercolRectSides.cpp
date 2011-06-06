#include "PercolRectSides.h"

/*
o----o----o----o----o
     |    |    |
o----o----o----o----o
     |    |    |
o----o----o----o----o

rows: 3
cols: 5
defined nodes: V(6)
undefined nodes: W(9)
currents: 3*(5-1)+(3-1)*(5-2) = 18

indexing vertices (v-nodes):
0----6----9---12----3
     |    |    |
1----7---10---13----4
     |    |    |
2----8---11---14----5

indexing edges (currents):
o--0>o--3>o--6>o--9>o
     ^    ^    ^
    12   14   16
o--1>o--4>o--7>o-10>o
     ^    ^    ^
    13   15   17
o--2>o--5>o--8>o-11>o
*/


PercolRectSides::PercolRectSides(int _rows, int _cols)
: Percol2D(), rows(_rows), cols(_cols)
{
    V.resize(2*rows); // left and right sides of the grid
    W.resize(rows*cols - 2*rows); // inner nodes of the grid
    I.resize((rows-1)*(cols-2) + rows*(cols-1));
    difV.resize( I.size() );
    IdifV.resize(I.size());
    Sigma.resize(I.size());

    for (int r=0; r<rows; ++r)
    {
        V[r] = 1.0; V[r+rows] = -1.0;
    }
}

PercolRectSides::~PercolRectSides(void)
{
}

