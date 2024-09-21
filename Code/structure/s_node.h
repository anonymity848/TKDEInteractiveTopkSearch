#ifndef RUN_S_NODE_H
#define RUN_S_NODE_H
#include "point_t.h"
#include "hyperplane.h"

class s_node
{
public:
    bool is_leaf;
    point_t* center;
    double angle;
    std::vector<s_node*> child;
    std::vector<hyperplane*> hyper;

    s_node(int dim);
};


#endif //RUN_S_NODE_H
