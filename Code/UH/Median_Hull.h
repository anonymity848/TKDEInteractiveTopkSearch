#ifndef RUN_MEDIAN_HULL_H
#define RUN_MEDIAN_HULL_H
#include "point_set.h"
#include "hyperplane_set.h"
#include "Partition.h"
#include "data_utility.h"

int Median(point_set* originalSet, point_t *u, int k, int s);

int Hull(point_set* originalSet, point_t *u, int k, int s);



#endif //RUN_MEDIAN_HULL_H
