#ifndef PERCOLRECTSIDES_H_INCLUDED
#define PERCOLRECTSIDES_H_INCLUDED
#include "percol2d.h"

// Rectangular percolation grid with source and drain at opposite
// sides (left and right).

class PercolRectSides : public Percol2D
{
    int rows, cols;
public:
    PercolRectSides(int _rows,int _cols);
    virtual ~PercolRectSides(void);

private:
    // Return id of vnode located at (row,col)
    int vnode(int row, int col) const
    {
        if (col == 0)
            return row;
        else if (col == cols-1)
            return row + rows;
        else
            return row + (col+1)*rows;
    }

    // Return (row,col) location of vnode v
    MYPAIR<int,int> rowcolv(int v) const
    {
        if (v < rows)
            return MYPAIR<int,int>(v % rows, 0);
        else if (v < 2*rows)
            return MYPAIR<int,int>(v % rows, cols-1);
        else
            return MYPAIR<int,int>(v % rows, (v/rows)-1);
    }

public:
    // Return a pair of v-nodes that current i flows <from,to>
    MYPAIR<int,int> ends(int i) const
    {
        int from,to;
        int total_horz_currents = rows*(cols-1);
        if (i < total_horz_currents)
        {
            // this is horisontal current
            int row = i % rows; // row of the current
            int col = i / rows; // col of the from-end of the current
            from = vnode(row,col);
            to   = vnode(row,col+1);
        }
        else
        {
            // this is vertical current
            int row = (i - total_horz_currents) % (rows-1); // row of to-end
            int col = (i - total_horz_currents) / (rows-1); // col of current
            to   = vnode(row,  col);
            from = vnode(row+1,col);
        }
        return MYPAIR<int,int>(from,to);
    }

    // Return -/+1 if current i comes from/to node v, 0 otherwise.
    int S(int i,int v) const
    {
        MYPAIR<int,int> rc = ends(i);
        if (v == rc.first) return -1;
        if (v == rc.second) return +1;
        return 0;
    }

    // Return list of currents coming out from vnode v
    MYARRAY<int> from(int v) const
    {
        MYPAIR<int,int> rc = rowcolv(v);
        int row = rc.first;
        int col = rc.second;
        if (col == cols-1)
        {
            // Right edge of the grid, no currents
            MYARRAY<int> res;
            return res;
        }
        else if (col == 0)
        {
            // Left edge of the grid, one h-current
            MYARRAY<int> res(1);
            res[0] = row;
            return res;
        }
        else if (row == 0)
        {
            // Top edge of the grid, one h-current, no ^-currents
            MYARRAY<int> res(1);
            res[0] = col*rows;
            return res;
        }
        else
        {
            // Middle of the grid, one h-current and one ^-current
            int total_horz_currents = rows*(cols-1);
            MYARRAY<int> res(2);
            res[0] = row + col*rows; // h-current
            res[1] = total_horz_currents + row-1 + (col-1)*(rows-1); // ^-current
            return res;
        }
    }
    MYARRAY<int> to(int v) const
    {
        MYPAIR<int,int> rc = rowcolv(v);
        int row = rc.first;
        int col = rc.second;
        if (col == 0)
        {
            // Left edge of the grid, no currents coming into it
            MYARRAY<int> res;
            return res;
        }
        else if (col == cols-1)
        {
            // Right edge of the grid, one h-current
            MYARRAY<int> res(1);
            res[0] = row + (col-1)*rows;
            return res;
        }
        else if (row == rows-1)
        {
            // Bottom side of the grid, one h-current and no ^-currents
            MYARRAY<int> res(1);
            res[0] = row + (col-1)*rows; // h-current
            return res;
        }
        else
        {
            // Middle of the grid, one h-current and one ^-current
            MYARRAY<int> res(2);
            int total_horz_currents = rows*(cols-1);
            res[0] = row + (col-1)*rows; // h-current
            res[1] = total_horz_currents + row + (col-1)*(rows-1); // ^-current
            return res;
        }
    }

    MYPAIR<double,double> xy(int v) const
    {
        MYPAIR<int,int> rc = rowcolv(v);
        return MYPAIR<double,double>(rc.second,rc.first);
    }

    int vnode(double x,double y) const
    {
        int r = int(y > 0 ? y+0.5 : y-0.5); //nearest int
        int c = int(x > 0 ? x+0.5 : x-0.5); //nearest int
        return vnode(r,c);
    }

    double xmin() const { return 0; }
    double xmax() const { return cols - 1; }
    double ymin() const { return 0; }
    double ymax() const { return rows - 1; }
};
#endif /*PERCOLRECTSIDES_H_INCLUDED*/
