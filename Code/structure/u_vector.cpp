//
// Created by 王伟程 on 2022/2/8.
//

#include "u_vector.h"

/**
 * @brief Constructor
 */
u_vector::u_vector(){}

/**
 * @brief Constructor
 * @param ut The x-axis of the utility vector
 */
u_vector::u_vector(double ut)
{
    x = ut;
}

/**
 * @brief Constructor
 * @param ut    The x-axis of the utility vector
 * @param up    The point from up to down
 * @param down  The point from down to up
 */
u_vector::u_vector(double ut, point_t *up, point_t *down)
{
    x = ut;
    point_up = up;
    point_down = down;
    point_up->surpass.push_back(point_down);
    point_up->surpassVector.push_back(ut);
}

/**
 * @brief insert the utility vector in the list based on descending order
 * @param lists     The list
 */
void u_vector::inserted(std::vector<u_vector *> &lists)
{
    int left = 0, right = lists.size() - 1, middle;
    while (left <= right)
    {
        middle = (left + right) / 2;
        if (lists[middle]->x > x)
            left = middle + 1;
        else
            right = middle - 1;
    }
    lists.insert(lists.begin() + left, this);
}




