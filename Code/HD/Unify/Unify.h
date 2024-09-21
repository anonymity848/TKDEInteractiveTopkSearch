#ifndef RUN_UNIFY_H
#define RUN_UNIFY_H
#include "point_set.h"
#include "hyperplane_set.h"
#include "Partition.h"
#include "choose_item.h"

int unify(point_set* pset, point_t *u, int k, int s, int CH, int HS);

int unify(point_set* pset, point_t *u, int k, int s, int CH, int HS, int error);

#endif //RUN_UNIFY_H
