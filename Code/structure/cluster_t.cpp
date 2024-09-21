#include "cluster_t.h"


/**
 * @brief Constructor
 */
cluster_t::cluster_t()
{
    center = NULL;
}


/**
 * @brief Constructor
 * @param dim
 */
cluster_t::cluster_t(int dim)
{
    center = new point_t(dim);
}
