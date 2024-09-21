#ifndef RUN_CLUSTER_T_H
#define RUN_CLUSTER_T_H
#include "point_t.h"
#include "hyperplane.h"

class cluster_t
{
public:
    point_t* center;
    std::vector<hyperplane*> h_set;

    cluster_t();
    cluster_t(int dim);
};


#endif //RUN_CLUSTER_T_H
