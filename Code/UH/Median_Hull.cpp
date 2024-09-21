#include "Median_Hull.h"
#include "sys/time.h"

#ifdef WIN32
#ifdef __cplusplus
extern "C" {
#endif
#endif

//#include "data_utility.h"

#include "../Qhull/mem.h"
#include "../Qhull/qset.h"
#include "../Qhull/libqhull.h"
#include "../Qhull/qhull_a.h"

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#if __MWERKS__ && __POWERPC__
#include <SIOUX.h>
#include <Files.h>
#include <console.h>
#include <Desk.h>

#elif __cplusplus
extern "C" {
int isatty(int);
}

#elif _MSC_VER
#include <io.h>
#define isatty _isatty
int _isatty(int);

#else
int isatty(int);  /* returns 1 if stdin is a tty
                   if "Undefined symbol" this can be deleted along with call in main() */
#endif

#ifdef WIN32
#ifdef __cplusplus
}
#endif
#endif



// Algorithm Median
int Median(point_set* originalSet, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    point_set* pset = new point_set(originalSet);
    int dim = pset->points[0]->dim;
    double numOfQuestion = 0;
    Partition *R = new Partition(dim);

    int M = pset->points.size();
    int **is_considered = new int *[M];
    for (int i = 0; i < M; i++)
    {
        is_considered[i] = new int[M];
        pset->points[i]->place = i;
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
            is_considered[i][j] = 0;
    }

    while (pset->points.size() > k)
    {
        numOfQuestion++;
        sort(pset->points.begin(), pset->points.end(), angleCmp());
        int start = 0, end = pset->points.size() - 1, delete_id;
        int idx1 = end / 2, idx2 = idx1;
        bool choose = false;
        int testcount = 0;
        while(!choose && testcount < ((pset->points.size()) * (pset->points.size())/2))
        {
            testcount++;
            idx2++;
            idx2 = idx2 % (end + 1);
            if(!pset->points[idx1]->is_same(pset->points[idx2])
                && !is_considered[pset->points[idx1]->place][pset->points[idx2]->place])
            {
                hyperplane *h = new hyperplane(pset->points[idx1], pset->points[idx2]);
                if (R->check_relation(h) == 0)
                    choose = true;
                delete h;
            }

            is_considered[pset->points[idx1]->place][pset->points[idx2]->place] = 1;
            is_considered[pset->points[idx2]->place][pset->points[idx1]->place] = 1;

            if((idx2 == idx1 - 1) && (!choose))
            {
                idx1++;
                idx1 = idx1 % (end + 1);
                idx2 = idx1;
            }
        }
        if(testcount >= ((pset->points.size()) * (pset->points.size())/2))
            break;
        double v1 = pset->points[idx1]->dot_product(u);
        double v2 = pset->points[idx2]->dot_product(u);
        if (v1 > v2)
        {
            delete_id = pset->points[idx2]->id;
            hyperplane *h = new hyperplane(pset->points[idx2], pset->points[idx1]);
            R->hyperplanes.push_back(h);
        }
        else
        {
            delete_id = pset->points[idx1]->id;
            hyperplane *h = new hyperplane(pset->points[idx1], pset->points[idx2]);
            R->hyperplanes.push_back(h);
        }
        R->set_ext_pts();
        R->find_possible_topK(pset, k);

        /*
        point_set* top_current = R->find_Rdominate_topK(pset, k, s);
        if(top_current->points.size() >= s)
        {
            top_current->printResult("Median", numOfQuestion, s, t1, 0, 0);
            //check the results (whether they are top-k points)
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


        for (int i = 0; i < pset->points.size(); i++)
        {
            if (pset->points[i]->id == delete_id)
            {
                pset->points.erase(pset->points.begin() + i);
                break;
            }
        }
        */
    }

    pset->printResult("Median", numOfQuestion, s, t1, 0, 0);

    point_set *resultSet = new point_set();
    originalSet->findTopk(u, k, resultSet);
    bool resultcheck = true;
    for(int i = 0; i < s; ++i)
    {
        if(!resultSet->checkExist(pset->points[i]))
        {
            resultcheck = false;
            break;
        }
    }
    std::cout << resultcheck <<"\n";
    return numOfQuestion;
}

// Algorithm Hull
int Hull(point_set *originalSet, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    point_set *pset = new point_set(originalSet);
    int dim = pset->points[0]->dim;
    double numOfQuestion = 0;
    Partition *R = new Partition(dim);


    int M = pset->points.size();
    int **is_considered = new int *[M];
    for (int i = 0; i < M; i++)
    {
        is_considered[i] = new int[M];
        pset->points[i]->place = i;
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
            is_considered[i][j] = 0;
    }

    while (pset->points.size() > k)
    {
        numOfQuestion++;
        sort(pset->points.begin(), pset->points.end(), angleCmp());
        int start = 0, end = pset->points.size()-1, delete_id;
        int idx1 = end / 3, idx2 = 2*idx1;
        bool choose = false;
        int count_change = 0, test_cout = 0;
        while(!choose && test_cout < (pset->points.size()*pset->points.size()/2))
        {
            idx2++;
            idx2 = idx2 % (end + 1);
            count_change++;
            test_cout++;
            if(idx2 == idx1)
            {
                idx2++;
                idx2 = idx2 % (end + 1);
                count_change++;
                test_cout++;
            }
            if(!pset->points[idx1]->is_same(pset->points[idx2])
                && !is_considered[pset->points[idx1]->place][pset->points[idx2]->place])
            {
                hyperplane *h = new hyperplane(pset->points[idx1], pset->points[idx2]);
                if (R->check_relation(h) == 0)
                    choose = true;
                delete h;
            }
            is_considered[pset->points[idx1]->place][pset->points[idx2]->place] = 1;
            is_considered[pset->points[idx2]->place][pset->points[idx1]->place] = 1;
            if((count_change >= end + 1) && (!choose))
            {
                idx1++;
                idx1 = idx1 % (end + 1);
                idx2 = (2*idx1) % (end + 1);
                count_change = 0;
            }
        }
        if(test_cout >= (pset->points.size()*pset->points.size()/2))
            break;
        double v1 = pset->points[idx1]->dot_product(u);
        double v2 = pset->points[idx2]->dot_product(u);
        hyperplane *h;
        if (v1 > v2)
        {
            delete_id = pset->points[idx2]->id;
            h = new hyperplane(pset->points[idx2], pset->points[idx1]);

        }
        else
        {
            delete_id = pset->points[idx1]->id;
            h = new hyperplane(pset->points[idx1], pset->points[idx2]);
        }
        R->hyperplanes.push_back(h);
        R->set_ext_pts();
        R->find_possible_topK(pset, k);

        /*
        point_set* top_current = R->find_Rdominate_topK(pset, k, s);
        if(top_current->points.size() >= s)
        {
            top_current->printResult("Hull", numOfQuestion, s, t1, 0, 0);
            //check the results (whether they are top-k points)
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
        */
    }

    pset->printResult("Hull", numOfQuestion, s, t1, 0, 0);

    point_set *resultSet = new point_set();
    originalSet->findTopk(u, k, resultSet);
    bool resultcheck = true;
    for(int i = 0; i < s; ++i)
    {
        if(!resultSet->checkExist(pset->points[i]))
        {
            resultcheck = false;
            break;
        }
    }
    std::cout << resultcheck <<"\n";
    return numOfQuestion;
}