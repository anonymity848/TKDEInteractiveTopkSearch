#include "Partition.h"
#include "../Others/lp.h"
#include "../Others/pruning.h"
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <algorithm>

/**
 * @brief Constructor
 */
Partition::Partition(){}

/**
 * @brief       Constructor: initial some hyperplanes to bound the utility space
 * @param       The point set
 */
Partition::Partition(int dim)
{
    if(dim > 1)
    {
        //initialize hyperplanes
        hyperplane *h;
        for(int i = 0; i < dim; ++i)
        {
            h = new hyperplane(dim);
            for (int j = 0; j < dim; ++j)
            {
                if (i == j)
                    h->norm[j] = -1;
                else
                    h->norm[j] = 0;
            }
            h->offset = 0;
            hyperplanes.push_back(h);
        }
        //set extreme points
        this->set_ext_pts();
        center = average_point();
    }
    else
    {
        cout<<"Error: The number of dimension is invalid.";
    }
}


/**
 * @brief Constructor
 * @param Pr The original Partition
 */
Partition::Partition(Partition *Pr)
{
    ID = Pr->ID;
    for(int i = 0; i < Pr->hyperplanes.size(); ++i)
    {
        hyperplanes.push_back(Pr->hyperplanes[i]);
    }
    for(int i = 0; i < Pr->ext_pts.size(); ++i)
    {
        ext_pts.push_back(new point_t(Pr->ext_pts[i]));
    }
}


/**
 * @brief Destructor
 *        Delete the hyperplanes and extreme points
 */
Partition::~Partition()
{
    //delete the hyperplanes
    int i = hyperplanes.size();
    while(i > 0)
    {
        delete hyperplanes[i-1];
        i--;
    }
    hyperplanes.clear();

    //delete the extreme points
    i = ext_pts.size();
    while(i > 0)
    {
        delete ext_pts[i-1];
        i--;
    }
    ext_pts.clear();

}


/**
    * @brief Prepare the file for computing the convex hull (the utility range R)
    *        via halfspace interaction
    * @param feasible_pt   A points inside R
    * @param filename      The name of the file written
    */
void Partition::write(point_t *feasible_pt, char *filename)
{
    //char filename[MAX_FILENAME_LENG];
    int dim = feasible_pt->dim, size = hyperplanes.size();
    ofstream wPtr;
    wPtr.open(filename, std::ios::out);
    wPtr.setf(ios::fixed, ios::floatfield);  // set as fixed model
    wPtr.precision(6);  // set precision to 6

    // write the feasible point
    wPtr << dim <<"\n" << 1 <<"\n";
    for(int i = 0; i < dim; i++)
        wPtr << feasible_pt->attr[i] << " ";
    wPtr << "\n";

    // write the hyperplane
    wPtr << dim + 1 <<"\n" << size << "\n";//record the offset as one dimension
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            wPtr << hyperplanes[i]->norm[j] <<" ";
        }
        wPtr << hyperplanes[i]->offset <<"\n";
    }
    wPtr.close();
}


/**
 * @brief   Print the information of the hyperplane set
 *          Including hyperplanes and extreme points of the intersection(convex hull)
 */
void Partition::print()
{
    //for(int i = 0; i < topSet->points.size(); ++i)
    //    topSet->points[i]->print();
    //print hyperplanes
    for (int i = 0; i < hyperplanes.size(); i++)
    {
        printf("hyperplane: ");
        hyperplanes[i]->print();
    }
    //print extreme points
    for (int i = 0; i < ext_pts.size(); i++)
    {
        printf("extreme point: ");
        ext_pts[i]->print();
    }
}

/**
 * @brief   Check whether the intersection of the hyperplane exists
 *          Set the extreme points of the hyperplane set.
 *          Refine the bounded hyperplanes
 * @return  Whether R is updated
 *          true  R is the different
 *          false R is the same as before
 */
bool Partition::set_ext_pts() {
    int M = hyperplanes.size(), size = 0;
    if (M <= 0) {
        printf("%s\n", "Error: No hyperplane in the set.");
        return false;
    }
    int dim = hyperplanes[0]->dim;
    char file1[MAX_FILENAME_LENG];
    sprintf(file1, "../output/hyperplane_data.txt");
    char file2[MAX_FILENAME_LENG];
    sprintf(file2, "../output/ext_pt.txt");

    hyperplane *h = new hyperplane(dim);
    for (int i = 0; i < dim; i++)
        h->norm[i] = 1;
    h->offset = -10000;
    hyperplanes.push_back(h);


    point_t *feasible_pt = find_feasible();//the feasible point

    if (feasible_pt == NULL)
    {
        hyperplanes.pop_back();
        cout << "The intersection is infeasible.\n";
        return false;
    }
    write(feasible_pt, file1);//write hyperplanes and the feasible point to file1,
    hyperplanes.pop_back();


    //conduct half space intersection and write results to file2
    FILE *rPtr, *wPtr;
    if ((rPtr = fopen(file1, "r")) == NULL) {
        fprintf(stderr, "Cannot open the data file.\n");
        exit(0);
    }
    wPtr = (FILE *) fopen(file2, "w");
    halfspace(rPtr, wPtr);
    fclose(rPtr);
    fclose(wPtr);

    //read extreme points in file2
    if ((rPtr = fopen(file2, "r")) == NULL) {
        fprintf(stderr, "Cannot open the data file %s.\n", file2);
        exit(0);
    }

    //update the set of extreme points
    fscanf(rPtr, "%i%i", &dim, &size);
    if (size > 0)
    {
        //delete all the original extreme points
        while (ext_pts.size() > 0) {
            point_t *pt = ext_pts[ext_pts.size() - 1];
            ext_pts.pop_back();
            delete pt;
        }
        ext_pts.clear();

        //input new extreme points
        for (int i = 0; i < size; ++i)
        {
            point_t *p = new point_t(dim);
            double dim_sum = 0;
            for (int j = 0; j < dim; ++j)
            {
                fscanf(rPtr, "%lf", &p->attr[j]);
                p->attr[j] /= 10000;
                dim_sum += p->attr[j];

            }
            if (dim_sum > 1e-4)
            {
                ext_pts.push_back(p);
                //p->print();
            }
        }

        // update the set of hyperplanes
        fscanf(rPtr, "%i", &size);
        int *hid = new int[size + 1];
        for (int i = 1; i <= size; ++i)
            fscanf(rPtr, "%i", &hid[i]);
        sort(hid, hid + size + 1);
        for (int i = M - 1, count = size; i > -1; --i) {
            if (hid[count] < i)
                hyperplanes.erase(hyperplanes.begin() + i);
            else
                --count;
        }
        delete[]hid;

        //for(int j = 0; j < hyperplanes.size(); ++j)
        //    hyperplanes[j]->print();
    }
    fclose(rPtr);
    return true;
}






/**
 * @brief   Check whether the intersection of the hyperplane exists
 *          Set the extreme points of the hyperplane set.
 *          Refine the bounded hyperplanes
 * @return  Whether R is updated
 *          true  R is the different
 *          false R is the same as before
 */
bool Partition::set_ext_pts_withoutHull() {
    int M = hyperplanes.size(), size = 0;
    if (M <= 0) {
        printf("%s\n", "Error: No hyperplane in the set.");
        return false;
    }
    int dim = hyperplanes[0]->dim;
    char file1[MAX_FILENAME_LENG];
    sprintf(file1, "../output/hyperplane_data.txt");
    char file2[MAX_FILENAME_LENG];
    sprintf(file2, "../output/ext_pt.txt");

    hyperplane *h = new hyperplane(dim);
    for (int i = 0; i < dim; i++)
        h->norm[i] = 1;
    h->offset = -1;
    hyperplanes.push_back(h);

    hyperplane *h2 = new hyperplane(dim);
    for (int i = 0; i < dim; i++)
        h2->norm[i] = -1;
    h2->offset = 0.99;
    hyperplanes.push_back(h2);

    point_t *feasible_pt = find_feasible();//the feasible point

    if (feasible_pt == NULL)
    {
        hyperplanes.pop_back();
        hyperplanes.pop_back();
        cout << "The intersection is infeasible.\n";
        return false;
    }
    //feasible_pt->print();
    write(feasible_pt, file1);//write hyperplanes and the feasible point to file1,
    hyperplanes.pop_back();
    hyperplanes.pop_back();


    //conduct half space intersection and write results to file2
    FILE *rPtr, *wPtr;
    if ((rPtr = fopen(file1, "r")) == NULL) {
        fprintf(stderr, "Cannot open the data file.\n");
        exit(0);
    }
    wPtr = (FILE *) fopen(file2, "w");
    halfspace(rPtr, wPtr);
    fclose(rPtr);
    fclose(wPtr);

    //read extreme points in file2
    if ((rPtr = fopen(file2, "r")) == NULL) {
        fprintf(stderr, "Cannot open the data file %s.\n", file2);
        exit(0);
    }

    //update the set of extreme points
    fscanf(rPtr, "%i%i", &dim, &size);
    if (size > 0)
    {
        //delete all the original extreme points
        while (ext_pts.size() > 0) {
            point_t *pt = ext_pts[ext_pts.size() - 1];
            ext_pts.pop_back();
            delete pt;
        }
        ext_pts.clear();

        //input new extreme points
        for (int i = 0; i < size; ++i)
        {
            point_t *p = new point_t(dim);
            double dim_sum = 0;
            for (int j = 0; j < dim; ++j)
            {
                fscanf(rPtr, "%lf", &p->attr[j]);
                dim_sum += p->attr[j];
            }
            if (dim_sum > 0.99)
            {
                ext_pts.push_back(p);
                //p->print();
            }
        }

        // update the set of hyperplanes
        fscanf(rPtr, "%i", &size);
        int *hid = new int[size + 1];
        for (int i = 1; i <= size; ++i)
            fscanf(rPtr, "%i", &hid[i]);
        sort(hid, hid + size + 1);
        for (int i = M - 1, count = size; i > -1; --i) {
            if (hid[count] < i)
                hyperplanes.erase(hyperplanes.begin() + i);
            else
                --count;
        }
        delete[]hid;

    }
    fclose(rPtr);
    return true;
}





/**
 * @brief build the partition based on the top-k points
 * @param pset              The dataset
 * @param topSet            The top-k points
 * @param cPointSet         The candidate top-k set
 * @param candHyper         The candidate hyperplanes (the boundries of partitions)
 * @param partitionSet      The partition set
 * @return  true    There exists a partition
 *          false   There does not exist a partition
 */
bool Partition::buildPartition(point_set *pset, point_set *topSet, std::vector<point_set*> &cPointSet,
                               std::vector<hyperplane*> &candHyper, std::vector<Partition*> &partitionSet)
{
    int dim = pset->points[0]->dim;
    //Initial a test
    int M = hyperplanes.size(), size = 0;
    char file1[MAX_FILENAME_LENG];
    sprintf(file1, "../output/hyperplane_data.txt");
    char file2[MAX_FILENAME_LENG];
    sprintf(file2, "../output/ext_pt.txt");

    Partition *cHull = new Partition();
    for (int i = 0; i < M; ++i)
        cHull->hyperplanes.push_back(new hyperplane(hyperplanes[i]));
    for(int i = 0; i < pset->points.size() * 0.05; ++i)
    {
        bool isTop = false;
        for (int j = 0; j < topSet->points.size(); ++j)
        {
            if (pset->points[i]->id == topSet->points[j]->id)
            {
                isTop = true;
                break;
            }
        }
        if(!isTop)
        {
            for (int j = 0; j < topSet->points.size(); ++j)
                cHull->hyperplanes.push_back(new hyperplane(pset->points[i], topSet->points[j]));
        }
    }
    if (cHull->set_ext_pts())
    {
        //cHull->print();
        for (int i = pset->points.size() * 0.05; i < pset->points.size(); ++i)
        {
            bool isTop = false;
            for (int j = 0; j < topSet->points.size(); ++j)
            {
                if (pset->points[i]->id == topSet->points[j]->id)
                {
                    isTop = true;
                    break;
                }
            }
            if(!isTop)
            {
                for (int j = 0; j < topSet->points.size(); ++j)
                {
                    hyperplane *h = new hyperplane(pset->points[i], topSet->points[j]);
                    if(cHull->check_relation(h) == 0)
                        cHull->hyperplanes.push_back(h);
                }
            }
        }
        if(!cHull->set_ext_pts())
            return false;
        partitionSet.push_back(cHull);

        for (int i = 0; i < cHull->hyperplanes.size(); ++i)
        {
            //check whether it is a useful boundary
            if (cHull->hyperplanes[i]->p_1 != NULL && cHull->hyperplanes[i]->p_2 != NULL)
            {
                //update the hyperplane set
                bool inH = false; //check whether it is in H
                for(int j = 0; j < candHyper.size(); ++j)
                {
                    if (candHyper[j]->p_1->id == cHull->hyperplanes[i]->p_1->id &&
                    candHyper[j]->p_2->id == cHull->hyperplanes[i]->p_2->id)
                    {
                        inH = true;
                        break;
                    }
                }
                if (!inH) //put the hyperplane in H
                    candHyper.push_back(
                            new hyperplane(cHull->hyperplanes[i]->p_1, cHull->hyperplanes[i]->p_2));
                //update the other partitions
                point_set *newps = new point_set();
                for(int j = 0; j < topSet->points.size(); ++j)
                {
                    if (cHull->hyperplanes[i]->p_2->id != topSet->points[j]->id)
                    {
                        newps->points.push_back(topSet->points[j]);
                    }
                    else
                    {
                        newps->points.push_back(cHull->hyperplanes[i]->p_1);
                    }
                }
                cPointSet.push_back(newps);
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}


/**
 * @brief Check whether there is a partition correponds to the top-k points in the set
 * @param pset      The datasets
 * @param topSet    The top-k points
 * @return
 */
bool Partition::isExist(point_set *pset, point_set *topSet)
{
    int M = hyperplanes.size(), size = 0, dim = hyperplanes[0]->dim;
    char file1[MAX_FILENAME_LENG];
    sprintf(file1, "../output/hyperplane_data.txt");
    char file2[MAX_FILENAME_LENG];
    sprintf(file2, "../output/ext_pt.txt");

    //list all the necessary hyperplanes
    Partition *cHull = new Partition();
    for(int i = 0; i < M; ++i)
        cHull->hyperplanes.push_back(new hyperplane(hyperplanes[i]));
    for (int i = pset->points.size(); i < pset->points.size(); ++i)
    {
        bool isTop = false;
        for (int j = 0; j < topSet->points.size(); ++j)
        {
            if (pset->points[i]->id == topSet->points[j]->id)
            {
                isTop = true;
                break;
            }
        }
        if(!isTop)
        {
            for (int j = 0; j < topSet->points.size(); ++j)
            {
                hyperplane *h = new hyperplane(pset->points[i], topSet->points[j]);
                if(cHull->check_relation(h) == 0)
                    cHull->hyperplanes.push_back(h);
            }
        }
    }
    point_t *feasible_pt = cHull->find_feasible();//the feasible point
    if (feasible_pt == NULL)
    {
        std::cout << "Warning: The partition cannot be found for the top-k points. \n";
        return false;
    }
    return true;
}

/**
 * @brief   Find a point inside the convex hull
 * @return  The point inside the convex hull
 */
point_t *Partition::find_feasible() {
    int M = hyperplanes.size();
    int D = hyperplanes[0]->dim;

    // D + 2variables: D for dim, 2 for additional var for feasible
    int* ia = new int[1 + (D + 2) * M];  //TODO: delete
    int* ja = new int[1 + (D + 2) * M];  //TODO: delete
    double* ar = new double[1 + (D + 2) * M];   //TODO: delete
    int i, j;

    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "find_feasible");
    glp_set_obj_dir(lp, GLP_MAX);


    glp_add_rows(lp, M);  // add D rows: q_1...q_D
    // Add rows q_1 ... q_D
    for (i = 1; i <= M; i++) {
        char buf[10];
        sprintf(buf, "q%d", i);
        glp_set_row_name(lp, i, buf);
        glp_set_row_bnds(lp, i, GLP_UP, 0, 0);
    }


    glp_add_cols(lp, D + 2);    // add D columns: v[1] ... v[D]
    // Add col v[1] ... v[D]
    for (i = 1; i <= D + 2; i++) {
        char buf[10];
        sprintf(buf, "v%d", i);

        glp_set_col_name(lp, i, buf);

        if(i <= D)
            glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0); // -infty <= v[i] < infty
        else if (i == D + 1)
            glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); // 0 <= v[i] < infty
        else
            glp_set_col_bnds(lp, i, GLP_UP, 0.0, D+1);

        if(i == D + 2)
            glp_set_obj_coef(lp, i, 1);  // objective: 0
        else
            glp_set_obj_coef(lp, i, 0.0);  // objective: 0
    }


    int counter = 1;
    // set value on row q1 ... qD
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= D + 2; j++) {

            ia[counter] = i; ja[counter] = j;

            if(j <= D)
            {
                ar[counter++] = hyperplanes[i-1]->norm[j-1];
                //printf("%lf ", hyperplane[i-1]->normal->coord[j-1]);
            }
            else if (j == D+1)
            {
                ar[counter++] = hyperplanes[i-1]->offset;
                //printf("%lf ", hyperplane[i-1]->offset);
            }
            else if (j == D+2)
            {
                ar[counter++] = 1;
                //printf("1.00000\n");
            }
        }
    }

    // loading data
    glp_load_matrix(lp, counter - 1, ia, ja, ar);

    // running simplex
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF; // turn off all message by glp_simplex

    glp_simplex(lp, &parm);


    point_t* feasible_pt = new point_t(D);
    double w1, w2;
    w1 = glp_get_col_prim(lp, D+1);
    w2 = glp_get_col_prim(lp, D+2);
    double status = glp_get_prim_stat(lp);

    if(w1 < EQN3 || w2 < EQN3 || isZero(w1) || isZero(w2))
    {
        //printf("LP feasible error.\n");
        return NULL;
    }
    for (i = 0; i < D; i++)
    {
        double v = glp_get_col_prim(lp, i + 1);
        //printf("w%d = %lf\n", i + 1, v);
        feasible_pt->attr[i] = v / w1;
    }
    feasible_pt->dim = D;

    glp_delete_prob(lp); // clean up
    delete[]ia;
    delete[]ja;
    delete[]ar;

    return feasible_pt;
}


double Partition::max_uq(hyperplane *h, point_t *q) {

    int M = hyperplanes.size();
    int D = hyperplanes[0]->dim;

    // D + 2variables: D for dim, 2 for additional var for feasible
    int* ia = new int[1 + D * (M + 2)];  //TODO: delete
    int* ja = new int[1 + D * (M + 2)];  //TODO: delete
    double* ar = new double[1 + D * (M + 2)];   //TODO: delete
    int i, j;

    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "max_uq");
    glp_set_obj_dir(lp, GLP_MAX);

    glp_add_rows(lp, M);  // add D rows: q_1...q_D
    // Add rows q_1 ... q_D
    for (i = 1; i <= M + 2; i++)
    {
        char buf[10];
        sprintf(buf, "q%d", i);
        glp_set_row_name(lp, i, buf);
        if(i <= M)
            glp_set_row_bnds(lp, i, GLP_UP, 0.0, -hyperplanes[i - 1]->offset); // qi <= offset
        else if (i == M + 1)
            glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0); // qi = offset
        else
            glp_set_row_bnds(lp, i, GLP_UP, 0.0, 10000.0); // qi = offset
    }

    glp_add_cols(lp, D);    // add D columns: v[1] ... v[D]
    // Add col v[1] ... v[D]
    for (i = 1; i <= D; i++)
    {
        char buf[10];
        sprintf(buf, "v%d", i);
        glp_set_col_name(lp, i, buf);
        glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0); // -infty <= v[i] < infty

        glp_set_obj_coef(lp, i, q->attr[i - 1]);  // objective
    }


    int counter = 1;
    // set value on row q1 ... qD
    for (i = 1; i <= M + 1; i++)
    {
        for (j = 1; j <= D; j++)
        {
            ia[counter] = i; ja[counter] = j;
            if(i <= M)
                ar[counter++] = hyperplanes[i-1]->norm[j-1];
            else if (i == M + 1)
                ar[counter++] = h->norm[j-1];
            else
                ar[counter++] = 1;
            //printf("%lf ", hyperplane[i-1]->normal->coord[j-1]);
        }
    }

    // loading data
    glp_load_matrix(lp, counter - 1, ia, ja, ar);

    // running simplex
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF; // turn off all message by glp_simplex

    glp_simplex(lp, &parm);

    double z = glp_get_obj_val(lp);
    point_t* feasible_u = new point_t(D);
    for (i = 0; i < D; i++)
    {
        double v = glp_get_col_prim(lp, i + 1);
        printf("v%d = %lf\n", i + 1, v);
        feasible_u->attr[i] = v;
    }

    glp_delete_prob(lp); // clean up
    delete[]ia;
    delete[]ja;
    delete[]ar;

    return z;
}







/**
 * @brief Check the relation between the hyperplane and the intersection of the hyperplane set
 *        Since the extreme points of the half_set can not be accurate enough, we set "Precision" to solve the error
 * @param h     The hyperplane
 * @return      1: partition on the positive side of the hyperplane
 *              -1: partition on the negative side of the hyperplane
 *              0: partition intersects with the hyperplane
 */
int Partition::check_relation(hyperplane *h)
{
    int M = ext_pts.size();
    if (M < 1)
    {
        printf("%s\n", "Warning: The extreme point set of the partition is empty.");
        return -2;
    }
    int positive = 0, negative = 0;
    for (int i = 0; i < M; i++)
    {
        int relation = h->check_position(ext_pts[i]);
        if (relation == -1)
            ++negative;
        else if (relation == 1)
            ++positive;
        if (positive > 0 && negative > 0)
        {
            return 0;
        }
    }
    if (negative > 0)
        return -1;
    else
        return 1;
}


/**
 * @brief Check whether the partition is on the positive side of h
 * @param h     The hyperplane
 * @return      1 it is on the positive side
 *             -1 it is not on the positive side
 */
int Partition::check_positive(hyperplane *h)
{
    int M = ext_pts.size();
    for (int i = 0; i < M; i++)
    {
        if (h->check_positive(ext_pts[i]) != 1)
            return -1;
    }
    return 1;
}



/**
 * @brief Check the relation between the hyperplane and the intersection of the hyperplane set
 *        Since the extreme points of the half_set can not be accurate enough, we set "Precision" to solve the error
 * @param h     The hyperplane
 * @return      1: partition on the positive side of the hyperplane
 *              -1: partition on the negative side of the hyperplane
 *              0: partition intersects with the hyperplane
 */
int Partition::check_relationlose(hyperplane *h)
{
    int M = ext_pts.size();
    if (M < 1)
    {
        printf("%s\n", "Warning: The extreme point set of the partition is empty.");
        return -2;
    }
    int positive = 0, negative = 0;
    for (int i = 0; i < M; i++)
    {
        int relation = h->check_positionlose(ext_pts[i]);
        if (relation == -1)
            ++negative;
        else if (relation == 1)
            ++positive;
        if (positive > 0 && negative > 0)
        {
            return 0;
        }
    }
    if (negative > 0)
        return -1;
    else
        return 1;
}




/**
 * @brief Check whether p1 R-dominates p2
 * @param p1 The first point
 * @param p2 The second point
 * @return   1 R-dominates
 *          -1 Does not R-dominate
 */
bool Partition::R_dominate(point_t *p1, point_t *p2)
{
    hyperplane *h = new hyperplane(p1, p2);
    int size = ext_pts.size();
    for (int i = 0; i < size; ++i)
    {
        if (h->check_positive(ext_pts[i]) != 1)
        {
            delete h;
            return false;
        }
    }
    delete h;
    return true;
}

/**
 * @brief Prune all the points which are R-dominated by at least one point
 *        There is a point which the closer to any expected point in R
 * @param The point set
 */
void Partition::nearestSkyline(point_set *pset)
{
    int i = 0;
    while (i < pset->points.size())
    {
        //cout << "ID: " <<pset->points[i]->id << "\n";
        //pset->points[i]->print();
        for(int j = 0; j < pset->points.size(); ++j)
        {
            //pset->points[j]->print();
            if(i != j && R_dominate(pset->points[j], pset->points[i]))
            {
                pset->points.erase(pset->points.begin() + i);
                --i;
                break;
            }
        }
        ++i;
    }
    //record the index of the point in the list
    for(int i = 0; i < pset->points.size(); ++i)
        pset->points[i]->index = i;
}


/**
 * @brief Calculate the average point of the extreme points
 * @param ap The average point
 */
point_t* Partition::average_point()
{
    int size = ext_pts.size();
    if(size <=0)
        return NULL;
    int dim = ext_pts[0]->dim;
    point_t* avgPoint = new point_t(dim);
    for(int j = 0; j < dim; ++j)
        avgPoint->attr[j] = 0;
    //calculate the sum
    for(int i = 0; i < size; ++i)
    {
        for (int j = 0; j < dim; ++j)
            avgPoint->attr[j] += ext_pts[i]->attr[j];
    }
    for(int j = 0; j < dim; ++j)
        avgPoint->attr[j] /= size;
    return avgPoint;
}


/**
 * @brief Find the k points which are the top-k points any u in R
 * @param pset     The point set
 * @param k        The parameter
 * @return         point set    The top-k points set
 *                 NULL         There do not exist such k points
 */
point_set* Partition::findTopk(point_set *pset, int k)
{
    int M = ext_pts.size();

    point_set *topSet1 = new point_set(), *topSet2 = new point_set();
    pset->findTopk(ext_pts[0], k, topSet1);
    //p_set->points[nearestPoint]->print();

    for(int i = 1; i < M; ++i)
    {
        pset->findTopk(ext_pts[i], k, topSet2);
        if(!topSet1->isSame(topSet2))
        {
            return NULL;
        }
    }
    return topSet1;
}


/**
 * @brief Find the s points which are the top-k points any u in R
 * @param pset     The point set
 * @param k        The parameter
 * @param s        The parameter
 * @return         point set    The top-k points set
 *                 NULL         There do not exist such k points
 */
point_set* Partition::findTopk(point_set *pset, int k, int s)
{
    if(ext_pts.size() <= 0)
        return NULL;

    int M = pset->points.size();
    point_set *topSet1 = new point_set(), *topSet2 = new point_set();
    pset->findTopk(ext_pts[0], k, topSet1);
    //p_set->points[nearestPoint]->print();
    int relation;

    for(int i = 0; i < topSet1->points.size(); ++i)
    {
        if(topSet1->points[i]->result == 1)
            topSet2->points.push_back(topSet1->points[i]);
        else
        {
            //previously scanned part
            for(int j = 0; j < topSet1->points[i]->notRdominateID.size(); ++j)
            {
                hyperplane *h = new hyperplane(topSet1->points[i],pset->points[topSet1->points[i]->notRdominateID[j]]);
                relation = check_positive(h);
                delete h;
                if (relation == 1)
                {
                    ++topSet1->points[i]->count;
                    if (topSet1->points[i]->count >= M - k)
                    {
                        topSet1->points[i]->result = 1;
                        topSet2->points.push_back(topSet1->points[i]);
                        break;
                    }
                    topSet1->points[i]->notRdominateID.erase(topSet1->points[i]->notRdominateID.begin() + j);
                    j--;
                }
            }

            if(topSet1->points[i]->notRdominateID.size() >= k)
                break;

            //unscanned part
            for (int j = topSet1->points[i]->currentScanID; j < M; ++j)
            {
                topSet1->points[i]->currentScanID = j + 1;
                if (topSet1->points[i]->id != pset->points[j]->id)
                {
                    hyperplane *h = new hyperplane(topSet1->points[i], pset->points[j]);
                    relation = check_positive(h);
                    delete h;
                    if (relation == 1)
                    {
                        ++topSet1->points[i]->count;
                        if (topSet1->points[i]->count >= M - k)
                        {
                            topSet1->points[i]->result = 1;
                            topSet2->points.push_back(topSet1->points[i]);
                            break;
                        }
                    }
                    else
                    {
                        topSet1->points[i]->notRdominateID.push_back(j);
                        if (topSet1->points[i]->notRdominateID.size() >= k)
                            break;
                    }
                }
            }
        }
        if(topSet2->points.size() >= s)
            break;
    }
    if(topSet2->points.size() >= s)
        return topSet2;
    else
        return NULL;
}


bool Partition::is_prune(hyperplane *h)
{
    int size = ext_pts.size();
    bool ispruned = true;
    for (int i = 0; i < size; ++i)
    {
        if (h->check_positive(ext_pts[i]) == -1)
        {
            ispruned = false;
            break;
        }
    }
    if(ispruned)
        return true;

    hyperplanes.push_back(h);
    if(!set_ext_pts())
        return true;
    else
        return false;
}


/**
 * @brief Find the upper and lower bound of the indexDimension-th dimension among the extreme points
 * @param indexDimenion     The index of the dimension
 * @param min               The maximum value
 * @param max               The minimum value
 */
void Partition::findMinMax(int indexDimenion, double &min, double &max)
{
    int size = ext_pts.size();
    min = INF; max = -INF;
    for(int i = 0; i < ext_pts.size(); ++i)
    {
        if(min > ext_pts[i]->attr[indexDimenion])
            min = ext_pts[i]->attr[indexDimenion];
        if(max < ext_pts[i]->attr[indexDimenion])
            max = ext_pts[i]->attr[indexDimenion];
    }
}

double Partition::findL1Dis(int dim)
{
    double L1 = 0, min, max;
    for(int i = 0; i < dim; ++i)
    {
        findMinMax(i, min, max);
        L1 = L1 + max - min;
    }
    return L1;
}



/**
 * @brief draw a line crossing the average points and calculate the intersecting points of the line and R
 * @param leftPts   One intersecting point
 * @param rightPts  One intersecting point
 * @param dimInex
 */
void Partition::findendPts(point_t *lowerPts, point_t *upperPts, int dimInex)
{
    int dim = hyperplanes[0]->dim;
    point_t* avgPts = average_point();
    for(int i = 0; i < dim; ++i)
    {
        lowerPts->attr[i] = avgPts->attr[i];
        upperPts->attr[i] = avgPts->attr[i];
    }


    lowerPts->attr[dimInex] = -INF;
    upperPts->attr[dimInex] = INF;

    for(int i = 0; i < hyperplanes.size(); ++i)
    {
        double sum = hyperplanes[i]->offset;
        for(int j = 0; j < dim; ++j)
        {
            if(j != dimInex)
                sum += hyperplanes[i]->norm[j] * avgPts->attr[j];
        }
        if(hyperplanes[i]->norm[dimInex] > 0)
        {
            sum = sum / (-hyperplanes[i]->norm[dimInex]);
            if(upperPts->attr[dimInex] > sum)
                upperPts->attr[dimInex] = sum;
        }
        else if(hyperplanes[i]->norm[dimInex] < 0)
        {
            sum = sum / (-hyperplanes[i]->norm[dimInex]);
            if(lowerPts->attr[dimInex] < sum)
                lowerPts->attr[dimInex] = sum;
        }
    }
}




/**
 * @brief Used to prune points which are not able to be the top-k based on utility space R
 * @param pset 	The dataset
 * @param k     The threshold top-k
 */
void Partition::find_possible_topK(point_set *pset, int k)
{
    for (int i = 0; i < pset->points.size(); ++i)
    {

        int num = 0;
        for (int j = 0; j < pset->points.size(); ++j)
        {
            if (!pset->points[i]->is_same(pset->points[j]))
            {
                hyperplane *h = new hyperplane(pset->points[j], pset->points[i]);
                if (check_relation(h) == 1)
                    num++;
                delete h;
            }
            if (num >= k)
            {
                pset->points.erase(pset->points.begin() + i);
                i--;
                break;
            }
        }


    }
}





/**
 * @brief Used to prune points which are not able to be the top-k based on utility space R
 * @param pset 	The dataset
 * @param k     The threshold top-k
 */
point_set* Partition::find_Rdominate_topK(point_set *pset, int k, int s)
{
    point_set *resultSet = new point_set();
    int M = pset->points.size();
    for (int i = 0; i < pset->points.size(); ++i)
    {
        int num = 0;
        for (int j = 0; j < pset->points.size(); ++j)
        {
            if (!pset->points[i]->is_same(pset->points[j]))
            {
                hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                if (check_relation(h) == 1)
                    num++;
                delete h;
            }
            if (num >= M - k)
            {
                resultSet->points.push_back(pset->points[i]);
                break;
            }
        }
        if(resultSet->points.size() >= s)
            return resultSet;
    }
}

/**
 * @brief Inset the hyperplanes which bounds the partition into hyperSet
 * @param hyperSet  The hyperplane set
 */
void Partition::insertHyperplane(std::vector<hyperplane*> &hyperSet)
{
    for(int i = 0; i < hyperplanes.size(); ++i)
    {
        if(hyperplanes[i]->p_1 != NULL && hyperplanes[i]->p_2 != NULL)
        {
            bool is_exist = false;
            for (int j = 0; j < hyperSet.size(); ++j)
            {
                if (hyperplanes[i]->is_same(hyperSet[j]))
                {
                    is_exist = true;
                    break;
                }
            }
            if (!is_exist)
                hyperSet.push_back(hyperplanes[i]);
        }
    }
}

/**
 * @brief Find the min and max value of each dimension in the partition
 * @param min  min values
 * @param max  max values
 */
void Partition::findMinMax(double *min, double *max)
{
    int dim = hyperplanes[0]->dim;
    for(int i = 0; i < dim; ++i)
    {
        min[i] = 99999;
        max[i] = -99999;
    }

    for(int i = 0; i < dim; ++i)
    {
        for(int j = 0; j < ext_pts.size(); ++j)
        {
            if(min[i] > ext_pts[j]->attr[i])
                min[i] = ext_pts[j]->attr[i];
            if(max[i] < ext_pts[j]->attr[i])
                max[i] = ext_pts[j]->attr[i];
        }
    }
}


/**
 * @brief Check whether the point is in the partition
 * @param u  The point
 * @return   1 It is inside
 *          -1 It is not inside
 */
bool Partition::is_inside(point_t *u)
{
    int dim = hyperplanes[0]->dim;
    for(int i = 0; i < hyperplanes.size(); ++i)
    {
        if(u->dot_product(hyperplanes[i]->norm) > 0)
            return false;
    }
    return true;
}


/**
 * @brief Check whether the point is in the partition
 * @param u  The point
 * @return   1 It is inside
 *          -1 It is not inside
 */
bool Partition::isIn(point_t* p)
{
    int dim = hyperplanes[0]->dim;
    for(int i = 0; i < hyperplanes.size(); ++i)
    {
        double sum = p->dot_product(hyperplanes[i]->norm);
        sum += hyperplanes[i]->offset;
        if(sum > 0)
            return false;
    }
    return true;
}




//@brief Use sampling to find all the points which is able to be top-1 at some utility vector
//@param p_set 		The point set containing all the points
//@param top_set	The returned point set containing all the possible top-1 point
//@param u 			The utility vector. For user, point_t* u = alloc_point(dim)
//@param level		The number of dimensions we have set. For user, only need to set level=0
//@param used_seg	The range which has been assigned to the u[i]. For user, set rest_seg=0
void Partition::findTopk_sampling(point_set *pset, std::vector<point_set*> &topSet, double *max, double *min, point_t *u, int k, int level, int used_seg)
{
    int dim = pset->points[0]->dim, M = pset->points.size();
    double segment = 3;
    if (level >= dim - 2)
    {
        for (int j = 0; j <= segment; j++)
        {
            u->attr[level] = min[level] + j * (max[level] - min[level]) / segment;
            u->attr[dim - 1] = 1;
            for (int i = 0; i < dim - 1; i++)
                u->attr[dim - 1] -= u->attr[i];

            //u->print();
            if(is_inside(u)) //Find the top-k set w.r.t u
            {
                point_set *topk = new point_set();
                pset->findTopk(u, k, topk);

                //Check if it is already in top_set
                bool is_inside = false;
                for (int i = 0; i < topSet.size(); i++)
                {
                    if (topSet[i]->isSame(topk))
                    {
                        is_inside = true;
                        break;
                    }
                }
                if (!is_inside)
                {
                    topSet.push_back(topk);
                    for (int i = 0; i < topk->points.size(); ++i)
                        topk->points[i]->topk = 1;

                }
            }
        }
    }
    else
    {
        for (int i = 0; i <= segment; i++)
        {
            u->attr[level] = min[level] + i * (max[level] - min[level]) / segment;
            findTopk_sampling(pset, topSet, max, min, u, k, level + 1, used_seg + i);
        }
    }
}






//@brief Use sampling to find all the points which is able to be top-1 at some utility vector
//@param p_set 		The point set containing all the points
//@param top_set	The returned point set containing all the possible top-1 point
//@param u 			The utility vector. For user, point_t* u = alloc_point(dim)
//@param level		The number of dimensions we have set. For user, only need to set level=0
//@param used_seg	The range which has been assigned to the u[i]. For user, set rest_seg=0
void Partition::findTopk_extreme(point_set *pset, std::vector<point_set*> &topSet, int k)
{
    int dim = pset->points[0]->dim, M = pset->points.size();
    point_t *u;
    for (int j = 0; j < ext_pts.size(); j++)
    {
        u = ext_pts[j];
        point_set *topk = new point_set();
        pset->findTopk(u, k, topk);

        //Check if it is already in top_set
        bool is_inside = false;
        for (int i = 0; i < topSet.size(); i++)
        {
            if (topSet[i]->isSame(topk))
            {
                is_inside = true;
                break;
            }
        }
        if (!is_inside)
        {
            topSet.push_back(topk);
            for (int i = 0; i < topk->points.size(); ++i)
                topk->points[i]->topk = 1;
        }
    }
}

void Partition::learning(point_t *est_u, double learning_rate, point_t* p1, point_t* p2, double ground)
{
    double estResult = 0;
    double v1 = p1->dot_product(est_u);
    double v2 = p2->dot_product(est_u);
    if(v1 > v2)
    {
        estResult = 1;
    }
    else
    {
        estResult = -1;
    }
    double error = ground - estResult;
    for(int i = 0; i < est_u->dim; ++i)
    {
        est_u->attr[i] += learning_rate * error * (p1->attr[i] - p2->attr[i]);
    }

}










