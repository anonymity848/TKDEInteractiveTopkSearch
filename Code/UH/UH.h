#ifndef MAXUTILITY_H
#define MAXUTILITY_H

#include "structure/define.h"
#include "structure/data_struct.h"
#include "structure/data_utility.h"

#include <vector>
#include <algorithm>
#include "structure/rtree.h"
#include "Others/lp.h"
#include "Others/pruning.h"
#include "Others/operation.h"
#include <queue>
#define RANDOM 1
#define SIMPLEX 2

using namespace std;

int UHRandom(point_set *originalSet, point_t *u, int k, int s, int error);

int UHSimplex(point_set *originalSet, point_t *u, int k, int s, int error);

#endif