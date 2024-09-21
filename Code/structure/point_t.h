#ifndef POINT_T_H
#define POINT_T_H
#include "define.h"

class point_t
{
public:
    int	id, index, dim, place;
    double *attr, value;
    std::vector<point_t*> surpass;
    std::vector<double> surpassVector;
    int topk, result;
    std::vector<int> notRdominateID;
    int count, currentScanID;


    explicit point_t(int dim);
    explicit point_t(int d_order, int d_unorder);
    explicit point_t(point_t *p);
    point_t(point_t *p1, point_t *p2);
    ~point_t();
    void print(); //print the point

    bool is_same(point_t *p);//check whether the point is the same as p
    bool is_zero();
    bool is_positive();
    double dot_product(point_t* p);
    double dot_product(double *v);
    int compare(point_t* p1, point_t* p2, int error);
    point_t* sub(point_t* p);
    point_t* add(point_t* p);
    point_t* scale(double c);
    double cal_len();
    double distance(point_t* p);
    double distance(double *p);
    bool dominate(point_t *p);
    double bound(point_t* p, double y);
    bool is_changed(point_t *p);
    int countPassed(double ut);
    double tranBound(point_t *p);
    void printResult(char *name, int Qcount, timeval t1);
    void printResult(char *name, int Qcount, timeval t1, double preTime, long mem_baseline);
};










#endif //U_2_POINT_T_H
