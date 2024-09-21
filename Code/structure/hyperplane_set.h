#ifndef RUN_HYPERPLANE_SET_H
#define RUN_HYPERPLANE_SET_H
#include "hyperplane.h"

class hyperplane_set
{
public:
    std::vector<hyperplane*> hyperplanes;
    hyperplane_set();
    hyperplane_set(point_set *pset, point_t *exp, double RandRate);
};


#endif //RUN_HYPERPLANE_SET_H
