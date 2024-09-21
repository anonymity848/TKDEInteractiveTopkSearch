#ifndef PRUNING_H
#define PRUNING_H

#include "../structure/data_struct.h"
#include "../structure/data_utility.h"
#include "../structure/point_set.h"
#include "../structure/Partition.h"

#include "operation.h"
#include "lp.h"
#include "../structure/rtree.h"
#include "../UH/frame.h"
#include "../qhull/io.h"
#include <queue>

using namespace std;

int halfspace(FILE* rPtr, FILE* wPtr);

// use the branch-and-bound skyline (BBS) algorithm for maintaining the candidate set
void rtree_pruning(point_set* P, vector<int>& C_idx, Partition* R);

#endif