#ifndef RUN_HDPI_H
#define RUN_HDPI_H
#include "point_set.h"
#include "hyperplane_set.h"
#include "Partition.h"
#include "choose_item.h"

void construct_partitions(point_set *pset, Partition *R, std::vector<point_set*> &topSet, std::vector<hyperplane*> &hyperSet,
                          std::vector<Partition*> &partitionSet, std::vector<int> &considered_halfset);

int build_choose_item_table(std::vector<Partition*> PartitionSet, std::vector<hyperplane*> hset,
                            std::vector<choose_item*> &choose_item_set);

int modify_choose_item_table(std::vector<choose_item*> &choose_item_set, std::vector<Partition*> PartitionSet,
                             std::vector<int> &considered_halfset, int H_num, bool direction);

int HDPI_sampling(point_set* pset, point_t *u, int k, int s);


int HDPI_grid(point_set* pset, point_t *u, int k, int s);

int HDPI_extreme(point_set* pset, point_t *u, int k, int s);










#endif //RUN_HDPI_H
