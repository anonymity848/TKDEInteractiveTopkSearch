#include "Groundtruth.h"

/**
 * @brief Find the nearest point of e
 * @param p_skyline The point set
 * @param e         The expected point
 */
void ground_truth(point_set *pset, point_t *u, int k)
{
    point_set *resultSet = new point_set();
    pset->findTopk(u, k, resultSet);
    //std::cout << "-----------------------------------------------------------------------------------\n";
    printf("|%15s |%15s |%15s |%15s |%10d |\n", "Ground Truth", "-", "-", "-", resultSet->points[0]->id);
    for(int i = 1; i < k; ++i)
        printf("|%15s |%15s |%15s |%15s |%10d |\n", "-", "-", "-", "-", resultSet->points[i]->id);
    std::cout << "-----------------------------------------------------------------------------------\n";
}
