#ifndef RUN_U_VECTOR_H
#define RUN_U_VECTOR_H
#include "point_t.h"

class u_vector
{
public:
    double x;
    int place;
    double time;
    point_t* point_up;
    point_t* point_down;

    u_vector();
    u_vector(double ut);
    u_vector(double ut, point_t *up, point_t *down);

    void inserted(std::vector<u_vector*> &lists);

};


#endif //RUN_U_VECTOR_H
