#ifndef POINT_SET_H
#define POINT_SET_H

#include "point_t.h"
#include <string>


class point_set
{
public:
    std::vector<point_t*> points;

    point_set();
    explicit point_set(point_set *p_set);
    explicit point_set(const char* input);
    ~point_set();

    void print();
    void random(double RandRate);
    point_set* sort(point_t *u);
    void findTopk(point_t *u, int k, point_set *topSet);
    void findRanking(point_t *u, point_set *rankingSet);
    void find_convexHull(std::vector<point_set *> &topSet);
    void write(std::string fileName);
    void prunePt(point_t* p);
    bool isSame(point_set *pset);
    bool isSame_exact(point_set *pset);
    point_set* findsame(point_set *pset);
    bool checkExist(point_t *p);
    void sort_point(std::vector<point_t*> &return_point, point_t *u);
    void sort_point(std::vector<point_t*> &return_point);
    void printResult(char *name, int Qcount, int s,timeval t1, double preTime, double acc);
    void skyband(point_set *returnSet, int k);
    int findBest(point_t *u);
    void findTopk_sampling(std::vector<point_set*> &topSet, point_t *u, int k, int level, int used_seg);
    void findRanking_sampling(std::vector<point_set*> &rankSet, point_t *u, int level, int used_seg);
    void findTopk_sampling(std::vector<point_set*> &topSet, double *max, double *min, point_t *u, int k, int level, int used_seg);

};


#endif //U_2_POINT_SET_H
