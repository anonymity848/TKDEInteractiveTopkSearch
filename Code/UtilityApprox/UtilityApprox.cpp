#include "UtilityApprox.h"
#include <sys/time.h>



/**
 * @brief Calculated the regret ratio of utility range
 * @param L the lower bound of each dimension
 * @param U the upper bound of each dimension
 * @param D the dimension
 * @return the regret ratio
 */
double rr_bound(double *L, double *U, int D)
{
    double bound = 0;

    for (int i = 0; i < D; i++)
    {
        bound += U[i] - L[i];
    }
    return bound;
}


/**
 * @brief Set the value of the extreme points of the rectangle
 * @param place      Set high/low value for the points
 * @param high       The largest value for each dimension
 * @param low        The smallest value for each dimension
 * @param ext_pts    All the extreme points
 * @param dim        The number of dimensions
 * @param layer      The number of dimensions we have set for the points
 */
void set_ext_on_R(bool *place, double *high, double *low, std::vector<point_t *> &ext_pts, int dim, int layer, int istar)
{
    if (layer + 1 == dim && layer != istar)
    {
        point_t *p1 = new point_t(dim);
        point_t *p2 = new point_t(dim);
        for (int i = 0; i < dim - 1; i++)
        {
            if (place[i] == true)
            {
                p1->attr[i] = high[i];
                p2->attr[i] = high[i];
            }
            else
            {
                p1->attr[i] = low[i];
                p2->attr[i] = low[i];
            }
        }
        p1->attr[dim - 1] = high[dim - 1];
        p2->attr[dim - 1] = low[dim - 1];
        ext_pts.push_back(p1);
        ext_pts.push_back(p2);
    }
    else if (layer + 1 == dim && layer == istar)
    {
        point_t *p1 = new point_t(dim);
        for (int i = 0; i < dim - 1; i++)
        {
            if (place[i] == true)
            {
                p1->attr[i] = high[i];
            }
            else
            {
                p1->attr[i] = low[i];
            }
        }
        p1->attr[dim - 1] = high[dim - 1];
        ext_pts.push_back(p1);
    }
    else if(layer + 1 != dim && layer != istar)
    {
        place[layer] = true;
        set_ext_on_R(place, high, low, ext_pts, dim, layer + 1, istar);
        place[layer] = false;
        set_ext_on_R(place, high, low, ext_pts, dim, layer + 1, istar);
    }
    else
    {
        place[layer] = true;
        set_ext_on_R(place, high, low, ext_pts, dim, layer + 1, istar);
    }
}


/**
 * @brief Danupon's Fake-Points algorithm
 * @param P the dataset
 * @param u the real utility vector
 * @param s the number of points shown in each question
 * @param epsilon the threshold of regret ratio
 * @param maxRound the upper bound of number of questions
 */
int utilityapprox(point_set *pset, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int D = pset->points[0]->dim;
    double delta, gamma, beta, regret = 1.0, sum;
    //double bestUtility, currUtility;
    int l, qIndex, alpha, t;
    int istar = 0;
    double *U = new double[D];
    double *L = new double[D];
    double *chi = new double[3];
    point_t *v = new point_t(D);
    point_set *q_set = new point_set();
    Partition *R = new Partition(D);


    // determine the attribute with highest rating in prefs
    for (int i = 1; i < D; ++i)
    {
        if (u->attr[i] > u->attr[istar])
            istar = i;
    }
    int numOfQuestion = D - 1;


    // initialize U and L
    for (int i = 0; i < D; ++i)
    {
        U[i] = 1.0;
        if (i == istar)
            L[i] = 1.0;
        else
            L[i] = 0.0;
    }

    // create memory for the K points that will be shown
    for (int i = 0; i < 2; ++i)
        q_set->points.push_back(new point_t(D));

    int i = 0;
    t = 1;
    std::vector<int> C;
    while ((rr_bound(L, U, D) > EQN2))
    {
        //printf("%d\n", C.size());
        //if (VERBOSE) printf("i = %d, rounds = %d, t = %d\n", i, rounds, t);
        if (i != istar)
        {
            // compute v
            for (int j = 0; j < D; ++j)
                v->attr[j] = (U[j] + L[j]) / 2.0;
            qIndex = pset->findBest(v);

            // calculate chi
            for (int j = 0; j <= 2; ++j)
                chi[j] = L[i] + j * (U[i] - L[i]) / 2;

            // compute delta as in latest version
            sum = 0.0;
            for (int j = 1; j < 2; ++j)
                sum = sum + chi[j] - L[i];
            delta = 1.0 / (pow(2.0, (t / (D - 1)) * 1.0) * sum);

            if (delta > 1e-5)
                delta = 1e-5;

            // compute gamma
            gamma = pset->points[qIndex]->attr[istar] / pset->points[qIndex]->attr[i];

            // create the set of fake points
            for (l = 0; l < D; ++l)
                q_set->points[1]->attr[l] = pset->points[qIndex]->attr[l];


            for (l = 0; l < D; ++l)
            {
                if (l == istar)
                {
                    q_set->points[0]->attr[l] = q_set->points[1]->attr[l] +
                            chi[1] * delta * gamma * q_set->points[1]->attr[i];
                }
                else if (l == i)
                {
                    q_set->points[0]->attr[l] =
                            q_set->points[1]->attr[l] - delta * gamma * q_set->points[1]->attr[i];
                }
                else
                {
                    q_set->points[0]->attr[l] = pset->points[qIndex]->attr[l];
                }
            }


            // compute beta and rescale points by it
            sum = 0.0;
            for (int j = 0; j < 2; ++j)
                sum += chi[j + 1] - L[i];
            beta = 1.0 / (1.0 + sum * delta);  // multiply by delta at the end to try and avoid underflow
            for (int j = 0; j < 2; ++j)
                for (l = 0; l < D; ++l)
                    q_set->points[j]->attr[l] = q_set->points[j]->attr[l] * beta;

            // check which point the user would pick and update L and U accordingly
            alpha = q_set->findBest(u);
            L[i] = chi[alpha];
            U[i] = chi[alpha + 1];
            ++numOfQuestion;


            bool *place = new bool[D];
            R->ext_pts.clear();
            set_ext_on_R(place, U, L, R->ext_pts, D, 0, istar);
            //R->print();

            point_set *top_current = R->findTopk(pset, k, s);
            if(top_current != NULL)
            {
                top_current->printResult("UtilityApprox", numOfQuestion, s, t1, 0, 0);

                point_set *resultSet = new point_set();
                pset->findTopk(u, k, resultSet);
                bool resultcheck = true;
                for(int i = 0; i < top_current->points.size(); ++i)
                {
                    if(!resultSet->checkExist(top_current->points[i]))
                    {
                        resultcheck = false;
                        break;
                    }
                }
                std::cout << resultcheck <<"\n";


                return numOfQuestion;
            }


        }
        i = (i + 1) % D; // cycle through all the dimensions
        t += 1;
    }
    delete[] U;
    delete[] L;
    delete[] chi;



    point_set *resultSet = new point_set();
    pset->findTopk(v, k, resultSet);
    resultSet->printResult("UtilityApprox", numOfQuestion, s, t1, 0, 0);

    return numOfQuestion;
}