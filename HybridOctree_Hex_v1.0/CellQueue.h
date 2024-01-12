#ifndef CELL_QUEUE_H
#define CELL_QUEUE_H

#include <stdio.h>
#include <math.h>
#include <memory>
#include <malloc.h>
#include <string>
#include <sys/types.h>
#include "Initialization.h"

class CellQueue {
public:
    inline CellQueue(int size = 10000);// constructor
    ~CellQueue();// destructor
    
    // add item to the queue
    inline void Add(unsigned int cell);

    // remove and return the first item in queue
    inline bool Get(int &cell);

    // return the first item in queue
    inline int Peek(int &cell);

    // remove the first item in queue
    inline void Pop();

    // reset to empty
    inline void Reset(void) { nEl = 0; }

    // check if queue is empty
    inline bool Empty(void) { return(nEl == 0); }

private:
    int nEl;
    int cellSize;// # of elements in cell array
    int start;
    unsigned int *cells;
};

// create a new cell queue with elements of specified size
inline CellQueue::CellQueue(int size)
{
    nEl = 0;
    start = 0;
    cellSize = size;
    cells    = (unsigned int *)malloc(sizeof(unsigned int) * cellSize);
}

// free storage
inline CellQueue::~CellQueue()
{
    if (cells != NULL)
        free(cells);
}

// add an item to the queue
inline void CellQueue::Add(unsigned int c)
{
    int n;
    int oldSize;
    int atEnd;

    n = nEl++;

    // resize the queue if needed
    if (nEl > cellSize) {
        oldSize = cellSize;
        cellSize *= 2;
        // ignore this warning, otherwise the computation cost increases
        cells = (unsigned int *)realloc(cells, sizeof(int) * cellSize);

        // move everything from 'start' to the 'end'
        if (start != 0) {
            atEnd = oldSize - start;
            memmove(&cells[cellSize - atEnd], &cells[start], sizeof(unsigned int)*atEnd);
            start = cellSize - atEnd;
        }
    }

    n += start;
    if (n >= cellSize)
        n -= cellSize;

    cells[n] = c;
}

// return the top item from the queue
inline bool CellQueue::Get(int &c)
{
    if (Peek(c) == -1)
        return false;
    Pop();
    return true;
}

// return the top item, but don't remove it
inline int CellQueue::Peek(int &c)
{
    if (nEl == 0)
        return(-1);

    c = cells[start];

    return(1);
}

// delete the top item in the queue
inline void CellQueue::Pop(void)
{
    start++;
    if (start == cellSize)
        start=0;
    nEl--;
}
#endif CELLQUEUE_H