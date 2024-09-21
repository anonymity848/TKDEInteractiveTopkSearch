#include "Scale.h"


/**
 * @brief Asking user question and return one of the top-k points
 *        Find the top-1 point by sampling
 * @param p_set 		 The dataset
 * @param u 			 The linear function
 * @param k 			 The parameter
 */
int scale(point_set* pset, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    double Alpha = 1;
    std::vector<point_t*> expSet;
    int dim = pset->points[0]->dim;
    Partition *R = new Partition(dim);
    point_t *exp1, *exp2;
    point_t *p1, *p2;
    double v1, v2;
    int numOfQuestion = 0;

    while (1)
    {
        for (int i = 0; i < R->ext_pts.size(); ++i)
        {
            expSet.push_back(new point_t(dim));
            for (int j = 0; j < dim; ++j)
                expSet[i]->attr[j] = R->ext_pts[i]->attr[j] * Alpha;
        }

        R->print();
        double dis = -9999;
        for (int i = 0; i < expSet.size(); ++i)
        {
            for (int j = i + 1; j < expSet.size(); ++j)
            {
                point_set *topSet1 = new point_set(), *topSet2 = new point_set();
                pset->findTopk(expSet[i], k, topSet1);
                pset->findTopk(expSet[j], k, topSet2);

                for(int ii = 0; ii < k; ++ii)
                {
                    for(int jj = 0; jj < k; ++jj)
                    {
                        hyperplane *hhh = new hyperplane(topSet1->points[ii], topSet2->points[jj]);
                        hyperplane *hh2 = new hyperplane( topSet2->points[jj], topSet1->points[ii]);
                        bool isin = 0;
                        for(int t = 0; t < R->hyperplanes.size(); ++t)
                        {
                            if(R->hyperplanes[t]->is_same(hhh) || R->hyperplanes[t]->is_same(hh2))
                            {
                                isin = 1;
                                break;
                            }
                        }
                        if(R->check_relation(hhh) == 0 && !isin)
                        {
                            p1 = topSet1->points[ii];
                            p2 = topSet2->points[jj];
                            R->check_relation(hhh);


                            hyperplane *hhh = new hyperplane(p1, p2);
                            bool isin = 0;
                            for(int t = 0; t < R->hyperplanes.size(); ++t)
                            {
                                if(R->hyperplanes[t]->is_same(hhh))
                                {
                                    isin = 1;
                                    break;
                                }
                            }


                            goto label;;
                        }
                    }
                }

            }
        }


label:

        hyperplane *hy;
        numOfQuestion++;
        v1 = p1->dot_product(u);
        v2 = p2->dot_product(u);
        if (v1 >= v2)
        {
            hy = new hyperplane(p2, p1);
        }
        else
        {
            hy = new hyperplane(p1, p2);
        }
        hy->print();
        R->hyperplanes.push_back(hy);
        R->set_ext_pts();
        point_set *top_current;
        top_current = R->findTopk(pset, k, s);
        if (top_current != NULL)
        {
            top_current->printResult("Scale", numOfQuestion, s, t1, 0, 0);

            point_set *resultSet = new point_set();
            pset->findTopk(u, k, resultSet);
            bool resultcheck = true;
            for (int i = 0; i < top_current->points.size(); ++i)
            {
                if (!resultSet->checkExist(top_current->points[i]))
                {
                    resultcheck = false;
                    break;
                }
            }
            std::cout << resultcheck << "\n";
            return numOfQuestion;
        }
    }
}