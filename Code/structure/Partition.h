#ifndef HYPERPLANE_SET_H
#define HYPERPLANE_SET_H

#include "point_t.h"
#include "point_set.h"
#include "hyperplane.h"


class Partition
{
public:
    int ID;
    point_set *topSet; //the ordered dimension of the expected point
    std::vector<hyperplane*> hyperplanes;
    std::vector<point_t*> ext_pts;
    point_t *center;

    Partition();
    explicit Partition(int dim);
    Partition(Partition *Pr);
    ~Partition();

    //Prepare the file for computing the convex hull (the utility range R) via halfspace interaction
    void write(point_t* feasible_pt, char* filename);
    void print();//print the information of the hyperplane set
    bool buildPartition(point_set *pset, point_set *topSet, std::vector<point_set*> &cPointSet,
                                   std::vector<hyperplane*> &candHyper, std::vector<Partition*> &partitionSet);
    bool isExist(point_set *pset, point_set *topSet);
    bool set_ext_pts();
    bool set_ext_pts_withoutHull();
    point_t* find_feasible();
    double max_uq(hyperplane *h, point_t *q);
    bool R_dominate(point_t *p1, point_t *p2);
    void learning(point_t *est_u, double learning_rate, point_t* p1, point_t* p2, double ground);
    void nearestSkyline(point_set *pset);
    int check_relation(hyperplane *h);//check the relation between the hyperplane and the hyperplane set
    int check_relationlose(hyperplane *h);
    int check_positive(hyperplane *h);
    point_t* average_point();
    point_set* findTopk(point_set *pset, int k);
    point_set* findTopk(point_set *pset, int k, int s);
    bool is_prune(hyperplane *h);
    void findMinMax(int indexDimenion, double &min, double &max);
    double findL1Dis(int dim);
    void findendPts(point_t* lowerPts, point_t *upperPts, int dimInex);
    void find_possible_topK(point_set *pset, int k);
    point_set* find_Rdominate_topK(point_set *pset, int k, int s);
    void insertHyperplane(std::vector<hyperplane*> &hyperSet);
    void findMinMax(double *min, double *max);
    bool is_inside(point_t *u);
    bool isIn(point_t* p);
    void findTopk_sampling(point_set *pset, std::vector<point_set*> &topSet, double *max, double *min, point_t *u,
                           int k, int level, int used_seg);
    void findTopk_extreme(point_set *pset, std::vector<point_set*> &topSet, int k);

};


#endif //U_2_HYPERPLANE_SET_H
