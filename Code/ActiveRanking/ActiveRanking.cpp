#include "ActiveRanking.h"


/**
 * @brief Ask user questions and give a ranking
 * @param original_set 		The original dataset
 * @param u 				The linear function
 * @param k 				The threshold top-k
 */
int ActiveRanking(point_set *pset, point_t *u, int k, int s, int error)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int dim = pset->points[0]->dim, numOfQuestion = 0, M = pset->points.size();;
    pset->random(0.5);

    //initialization
    Partition *R = new Partition(dim);
    point_set *current = new point_set();
    current->points.push_back(pset->points[0]);

    std::cout << M << "\n";
    //store all the points in order
    for (int i = 1; i < M; i++) //compare: p_set contains all the points
    {
        if(i % 1000 == 0)
            std::cout << i <<"\n";
        int num_point = current->points.size();
        int place = 0; //the place of the point inserted into the current_use
        //find the question asked user
        for (int j = 0; j < num_point; j++)
        {
            hyperplane *h = new hyperplane(pset->points[i], current->points[j]);
            int relation = R->check_relationlose(h);
            delete h;
            //if intersect, calculate the distance
            if (relation == 0)
            {
                numOfQuestion++;
                //double v1 = pset->points[i]->dot_product(u);
                //double v2 = current->points[j]->dot_product(u);
                int compareResult = u->compare(pset->points[i], current->points[j], error);
                if (compareResult == 1)
                {
                    hyperplane *h = new hyperplane(current->points[j], pset->points[i]);
                    R->hyperplanes.push_back(h);
                    if(!R->set_ext_pts())
                        R->hyperplanes.pop_back();
                    break;

                }
                else
                {
                    hyperplane *h = new hyperplane(pset->points[i], current->points[j]);
                    R->hyperplanes.push_back(h);
                    if(!R->set_ext_pts())
                        R->hyperplanes.pop_back();
                    place = j + 1;
                }
                //R->print();
            }
            else if (relation == -1)
            {
                place = j + 1;
            }
            else
            {
                break;
            }
        }
        current->points.insert(current->points.begin() + place, pset->points[i]);
    }

    point_set *resultSet = new point_set();
    pset->findTopk(u, k, resultSet);
    double groudtruthsum = 0;
    for(int i = 1; i <= s; ++i)
    {
        groudtruthsum += resultSet->points[resultSet->points.size() - i]->dot_product(u);
    }
    double testsum = 0;
    for (int i = 0; i < s; ++i)
    {
        testsum += current->points[i]->dot_product(u);
    }
    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

    current->printResult("ActiveRanking", numOfQuestion, s, t1, 0, accuracy);





    return numOfQuestion;

}










