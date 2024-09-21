#include "HDPI.h"
#include <sys/time.h>
#define Beta 0.1
/**
 * @brief Build all the partition(intersection of the halfspace), each partition corresponds to a top-1 point
 * @param p_set 				The set containing all the points which are possible to be the top-1 point
 *                              for some utility vectors
 * @param choose_item_points		The set containing points used to build choose_item
 * @param half_set_set			The returned partitions
 * @param considered_half_set	The set recorded all the partitions in the form of index
 */
void construct_partitions(point_set *pset, std::vector<point_set*> &topSet,
                             std::vector<Partition*> &partitionSet, std::vector<int> &considered_halfset)
{
    int M = pset->points.size(), dim = pset->points[0]->dim;
    Partition *halfset;
    //record we use i-th point as the pivot(p[i]>p)
    for(int i = 0; i < topSet.size(); ++i)
    {
        //std::cout << i << "\n";
        halfset = new Partition(dim);
        for (int j = 0; j < M; j++)
        {
            if (!topSet[i]->checkExist(pset->points[j]))
            {
                for(int l = 0; l < topSet[i]->points.size(); ++l)
                {
                    hyperplane *h = new hyperplane(pset->points[j], topSet[i]->points[l]);
                    halfset->hyperplanes.push_back(h);
                }
            }
        }
        if (halfset->set_ext_pts())
        {
            halfset->topSet = topSet[i];
            for(int k = 0; k < topSet[i]->points.size(); ++k)
                topSet[i]->points[k]->topk = 1;
            partitionSet.push_back(halfset);
            considered_halfset.push_back(partitionSet.size() - 1);
        }
    }
}




/**
 * @brief Build all the partition(intersection of the halfspace), each partition corresponds to a top-1 point
 * @param p_set 				The set containing all the points which are possible to be the top-1 point
 *                              for some utility vectors
 * @param choose_item_points		The set containing points used to build choose_item
 * @param half_set_set			The returned partitions
 * @param considered_half_set	The set recorded all the partitions in the form of index
 */
void construct_partitions(point_set *pset, Partition *R, std::vector<point_set*> &topSet, std::vector<hyperplane*> &hyperSet,
                          std::vector<Partition*> &partitionSet, std::vector<int> &considered_halfset)
{
    int M = pset->points.size(), dim = pset->points[0]->dim;
    Partition *halfset;
    //record we use i-th point as the pivot(p[i]>p)
    for(int i = 0; i < topSet.size(); ++i)
    {
        halfset = new Partition(dim);
        for (int j = 0; j < M; j++)
        {
            if (!topSet[i]->checkExist(pset->points[j]))
            {
                for(int l = 0; l < topSet[i]->points.size(); ++l)
                {
                    hyperplane *h = new hyperplane(pset->points[j], topSet[i]->points[l]);
                    halfset->hyperplanes.push_back(h);
                }
            }
        }
        for(int j = 0; j < R->hyperplanes.size(); ++j)
            halfset->hyperplanes.push_back(R->hyperplanes[j]);
        if (halfset->set_ext_pts())
        {
            halfset->topSet = topSet[i];
            halfset->insertHyperplane(hyperSet);
            partitionSet.push_back(halfset);
            halfset->ID = partitionSet.size() - 1;
            considered_halfset.push_back(partitionSet.size() - 1);
        }
    }
}




/**
 * @brief Build the choose_item table used for selecting the hyperplane(question) asked user
 * @param half_set_set 		All the partitions
 * @param p_set 			All the points which used to build choose_item
 * @param choose_item_set 	The returned choose_item
 * @return The index of the choose_item(hyperplane) which divide the half_set most evenly
 */
int build_choose_item_table(std::vector<Partition*> PartitionSet, std::vector<hyperplane*> hset,
                            std::vector<choose_item*> &choose_item_set)
{
    int HM = hset.size(); //The number of points used to construct hyperplane(questions)
    int PM = PartitionSet.size(); //The number of halfspace sets
    double ES_h = -1000; //Even score for finding hyperplane(question)
    int num_hyperplane = 0; //The index of the chosen hyperplane(question)
    for (int i = 0; i < HM; ++i)
    {
        //hset[i]->print();
        if(hset[i]->p_1 != NULL && hset[i]->p_2 != NULL)
        {
            choose_item *c_item = new choose_item();
            choose_item_set.push_back(c_item);
            for (int j = 0; j < PM; ++j)
            {
                c_item->h = hset[i];
                int which_side = PartitionSet[j]->check_relation(hset[i]);
                if (which_side == 1)
                {
                    c_item->positiveSide.push_back(PartitionSet[j]->ID);
                } else if (which_side == -1)
                {
                    c_item->negativeSide.push_back(PartitionSet[j]->ID);
                } else if (which_side == 0)
                {
                    c_item->intersectCase.push_back(PartitionSet[j]->ID);
                } else
                {
                    printf("Error: check side failed.\n");
                    return -1;
                }
            }

            //calculate the even score
            double ES_i;
            int p_s = c_item->positiveSide.size();
            int n_s = c_item->negativeSide.size();
            int i_s = c_item->intersectCase.size();
            if (p_s > n_s)
            {
                ES_i = n_s - (Beta * i_s);
            } else
            {
                ES_i = p_s - (Beta * i_s);
            }
            //renew even score if necessary
            if (ES_h < ES_i)
            {
                ES_h = ES_i;
                num_hyperplane = choose_item_set.size() - 1;
            }
            //std::cout << ES_i <<"\n";
        }
    }
    //std::cout << ES_h <<"\n";
    return num_hyperplane;
}





/**
 * @brief Based on the answer from the user, modify the choose_item table
 * @param choose_item_set 		The choose_item table which will be modified
 * @param half_set_set 			All the partitions
 * @param considered_half_set 	Based on the user's answer, the partitions still considered
 * @param H_num					The choose_item used for asking question
 * @param direction 			The user's answer(which side)
 *                              -true on the positive side
 *                              -false on the negative side
 * @return The choose_item(hyperplane) which divide the half_set most evenly
 */
int modify_choose_item_table(std::vector<choose_item*> &choose_item_set, std::vector<Partition*> PartitionSet,
                             std::vector<int> &considered_halfset, int H_num, bool direction)
{
    int M = choose_item_set.size();
    //refine the considered half_set
    if (direction == true)//delete all the negative side partitions
    {
        int consider_count = 0; //index for scanning the considered_half_set
        for (int i = 0; i < choose_item_set[H_num]->negativeSide.size(); i++)
        {
            while (consider_count < considered_halfset.size() &&
            considered_halfset[consider_count] <= choose_item_set[H_num]->negativeSide[i])
            {

                if (considered_halfset[consider_count] == choose_item_set[H_num]->negativeSide[i])
                    considered_halfset.erase(considered_halfset.begin() + consider_count);
                else
                    consider_count++;
            }
        }
    }
    else
    {
        int consider_count = 0;
        for (int i = 0; i < choose_item_set[H_num]->positiveSide.size(); i++)
        {
            while (consider_count < considered_halfset.size() &&
                   considered_halfset[consider_count] <= choose_item_set[H_num]->positiveSide[i])
            {
                if (considered_halfset[consider_count] == choose_item_set[H_num]->positiveSide[i])
                {
                    considered_halfset.erase(considered_halfset.begin() + consider_count);
                }
                else
                {
                    consider_count++;
                }
            }
        }
    }


    //build a halfspace based on choose_item[H_num]
    hyperplane *half;
    if (direction == true)
        half = new hyperplane(choose_item_set[H_num]->h->p_2, choose_item_set[H_num]->h->p_1);
    else
        half = new hyperplane(choose_item_set[H_num]->h->p_1, choose_item_set[H_num]->h->p_2);

    std::vector<int> intersectDelete;
    //refine all the half_set in the intersect.case
    for (int i = 0; i < choose_item_set[H_num]->intersectCase.size(); i++)
    {
        //num_set: index of the halfspace_set in half_set_set
        int num_set = choose_item_set[H_num]->intersectCase[i];
        PartitionSet[num_set]->hyperplanes.push_back(half);
        if(!PartitionSet[num_set]->set_ext_pts())
        {
            for(int j = 0; j < considered_halfset.size(); ++j)
            {
                if(considered_halfset[j] == num_set)
                    considered_halfset.erase(considered_halfset.begin() + j);
            }
            intersectDelete.push_back(num_set);
        }
    }


    //delete all the half_sets which are on the non-direction side
    if (direction)
    {
        //delete negative side
        for (int j = 0; j < M; j++)//scan all the choose_item, j-th choose_item
        {
            if (H_num != j)
            {
                //deal with one item, delete all the half_set on the negative side of the hyperplane
                int positive_count = 0, negative_count = 0, intersect_count = 0;
                for (int i = 0; i < choose_item_set[H_num]->negativeSide.size(); i++)
                {
                    //deal with vector positive_side
                    while (positive_count < choose_item_set[j]->positiveSide.size() &&
                           choose_item_set[j]->positiveSide[positive_count] <=
                           choose_item_set[H_num]->negativeSide[i])
                    {
                        if (choose_item_set[j]->positiveSide[positive_count] ==
                            choose_item_set[H_num]->negativeSide[i])
                        {
                            choose_item_set[j]->positiveSide.erase(
                                    choose_item_set[j]->positiveSide.begin() + positive_count);
                        }
                        else
                        {
                            positive_count++;
                        }
                    }

                    //deal with vector negative_side
                    while (negative_count < choose_item_set[j]->negativeSide.size() &&
                           choose_item_set[j]->negativeSide[negative_count] <=
                           choose_item_set[H_num]->negativeSide[i])
                    {
                        if (choose_item_set[j]->negativeSide[negative_count] ==
                            choose_item_set[H_num]->negativeSide[i])
                        {
                            choose_item_set[j]->negativeSide.erase(
                                    choose_item_set[j]->negativeSide.begin() + negative_count);
                        }
                        else
                        {
                            negative_count++;
                        }
                    }

                    //deal with vector intersect_side
                    while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                           choose_item_set[j]->intersectCase[intersect_count] <=
                           choose_item_set[H_num]->negativeSide[i])
                    {
                        if (choose_item_set[j]->intersectCase[intersect_count] ==
                            choose_item_set[H_num]->negativeSide[i])
                        {
                            choose_item_set[j]->intersectCase.erase(
                                    choose_item_set[j]->intersectCase.begin() + intersect_count);
                        }
                        else
                        {
                            intersect_count++;
                        }
                    }
                }
            }
        }
    }
    else
    {
        //delete positive side
        for (int j = 0; j < M; j++)//scan all the choose_item, j-th choose_item
        {
            if (H_num != j)
            {
                //deal with one item, delete all the half_set on the positive side of the hyperplane
                int positive_count = 0, negative_count = 0, intersect_count = 0;
                for (int i = 0; i < choose_item_set[H_num]->positiveSide.size(); i++)
                {
                    //deal with vector positive_side
                    while (positive_count < choose_item_set[j]->positiveSide.size() &&
                           choose_item_set[j]->positiveSide[positive_count] <=
                           choose_item_set[H_num]->positiveSide[i])
                    {
                        if (choose_item_set[j]->positiveSide[positive_count] ==
                            choose_item_set[H_num]->positiveSide[i])
                        {
                            choose_item_set[j]->positiveSide.erase(
                                    choose_item_set[j]->positiveSide.begin() + positive_count);
                        }
                        else
                        {
                            positive_count++;
                        }
                    }

                    //deal with vector negative_side
                    while (negative_count < choose_item_set[j]->negativeSide.size() &&
                           choose_item_set[j]->negativeSide[negative_count] <=
                           choose_item_set[H_num]->positiveSide[i])
                    {
                        if (choose_item_set[j]->negativeSide[negative_count] ==
                            choose_item_set[H_num]->positiveSide[i])
                        {
                            choose_item_set[j]->negativeSide.erase(
                                    choose_item_set[j]->negativeSide.begin() + negative_count);
                        }
                        else
                        {
                            negative_count++;
                        }
                    }

                    //deal with vector intersect_side
                    while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                           choose_item_set[j]->intersectCase[intersect_count] <=
                           choose_item_set[H_num]->positiveSide[i])
                    {
                        if (choose_item_set[j]->intersectCase[intersect_count] ==
                            choose_item_set[H_num]->positiveSide[i])
                        {
                            choose_item_set[j]->intersectCase.erase(
                                    choose_item_set[j]->intersectCase.begin() + intersect_count);
                        }
                        else
                        {
                            intersect_count++;
                        }
                    }
                }
            }
        }
    }

    for (int j = 0; j < M; j++)//scan all the choose_item, j-th choose_item
    {
        if (H_num != j)
        {
            //deal with one item, delete all the half_set on the negative side of the hyperplane
            int positive_count = 0, negative_count = 0, intersect_count = 0;
            for (int i = 0; i < intersectDelete.size(); i++)
            {
                //deal with vector positive_side
                while (positive_count < choose_item_set[j]->positiveSide.size() &&
                       choose_item_set[j]->positiveSide[positive_count] <= intersectDelete[i])
                {
                    if (choose_item_set[j]->positiveSide[positive_count] == intersectDelete[i])
                    {
                        choose_item_set[j]->positiveSide.erase(
                                choose_item_set[j]->positiveSide.begin() + positive_count);
                    }
                    else
                    {
                        positive_count++;
                    }
                }

                //deal with vector negative_side
                while (negative_count < choose_item_set[j]->negativeSide.size() &&
                       choose_item_set[j]->negativeSide[negative_count] <= intersectDelete[i])
                {
                    if (choose_item_set[j]->negativeSide[negative_count] == intersectDelete[i])
                    {
                        choose_item_set[j]->negativeSide.erase(
                                choose_item_set[j]->negativeSide.begin() + negative_count);
                    }
                    else
                    {
                        negative_count++;
                    }
                }

                //deal with vector intersect_side
                while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                       choose_item_set[j]->intersectCase[intersect_count] <= intersectDelete[i])
                {
                    if (choose_item_set[j]->intersectCase[intersect_count] == intersectDelete[i])
                    {
                        choose_item_set[j]->intersectCase.erase(
                                choose_item_set[j]->intersectCase.begin() + intersect_count);
                    }
                    else
                    {
                        intersect_count++;
                    }
                }
            }
        }
    }

    //deal with the refined half_set
    for (int i = 0; i < choose_item_set[H_num]->intersectCase.size(); i++)
    {
        int num_set = choose_item_set[H_num]->intersectCase[i];
        //scan all the choose_item, j-th choose_item
        for (int j = 0; j < M; j++)
        {
            if (H_num != j)
            {
                //deal with vector intersect_case
                int intersect_count = 0;
                while (intersect_count < choose_item_set[j]->intersectCase.size() &&
                       choose_item_set[j]->intersectCase[intersect_count] <= num_set)
                {
                    if (choose_item_set[j]->intersectCase[intersect_count] == num_set)
                    {
                        int which_side = PartitionSet[num_set]->check_relation(choose_item_set[j]->h);
                        if (which_side == -2)
                        {
                            printf("%s\n", "Error: check side failed.");
                            return -1;
                        }
                        if (which_side != 0)
                        {
                            choose_item_set[j]->intersectCase.erase(
                                    choose_item_set[j]->intersectCase.begin() + intersect_count);
                            if (which_side == 1)//the half_set is on the positive side
                            {
                                bool is_insert = false;
                                for (int q = 0; q < choose_item_set[j]->positiveSide.size(); q++)
                                {
                                    if (num_set < choose_item_set[j]->positiveSide[q])
                                    {
                                        choose_item_set[j]->positiveSide.insert(
                                                choose_item_set[j]->positiveSide.begin() + q, num_set);
                                        is_insert = true;
                                        break;
                                    }
                                }
                                if (is_insert == false)
                                {
                                    choose_item_set[j]->positiveSide.push_back(num_set);
                                }

                            }
                            else//the half_set is on the negative side
                            {
                                bool is_insert = false;
                                for (int q = 0; q < choose_item_set[j]->negativeSide.size(); q++)
                                {
                                    if (num_set < choose_item_set[j]->negativeSide[q])
                                    {
                                        choose_item_set[j]->negativeSide.insert(
                                                choose_item_set[j]->negativeSide.begin() + q, num_set);
                                        is_insert = true;
                                        break;
                                    }
                                }
                                if (is_insert == false)
                                {
                                    choose_item_set[j]->negativeSide.push_back(num_set);
                                }
                            }

                        }
                        break;
                    }
                    intersect_count++;
                }
            }
        }
    }

    //refine the choose_item_set and select the hyperplane asked user
    choose_item_set.erase(choose_item_set.begin() + H_num);

    int item_count = 0;                         //the index for scanning the choose_item_set
    int ES_h = -1000;                           //Even score
    int num_hyperplane = 0;                     //index of the chosen hyperplane(question)
    int choose_size = choose_item_set.size();   //the original size of the choose_item_set
    for (int i = 0; i < choose_size; i++)
    {
        choose_item *c_item = choose_item_set[item_count];
        int M_positive = choose_item_set[item_count]->positiveSide.size();
        int M_negative = choose_item_set[item_count]->negativeSide.size();
        int M_intersect = choose_item_set[item_count]->intersectCase.size();
        if (M_negative + M_intersect == 0 || M_positive + M_intersect == 0)// || M_negative + M_positive == 0)
        {
            choose_item_set.erase(choose_item_set.begin() + item_count);
        }
        else
        {
            if (M_positive < M_negative)
            {
                double ES_i = M_positive - (Beta * M_intersect);
                if (ES_h < ES_i)
                {
                    ES_h = ES_i;
                    num_hyperplane = item_count;
                }
            }
            else
            {
                double ES_i = M_negative - (Beta * M_intersect);
                if (ES_h < ES_i)
                {
                    ES_h = ES_i;
                    num_hyperplane = item_count;
                }
            }
            item_count++;
        }
    }
    if(choose_item_set[num_hyperplane]->positiveSide.size() == 0 || choose_item_set[num_hyperplane]->negativeSide.size() == 0)
        return -1;
    else
        return num_hyperplane;
}



/**
 * @brief Asking user question and return one of the top-k points
 *        Find the top-1 point by sampling
 * @param p_set 		 The dataset
 * @param u 			 The linear function
 * @param k 			 The parameter
 */
int HDPI_sampling(point_set* pset, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    std::vector<point_set*> topSet;
    int dim = pset->points[0]->dim;
    pset->findTopk_sampling(topSet, new point_t(dim), k, 0, 0);//use sampling method
    point_set *candPt = new point_set();
    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(pset->points[i]->topk == 1)
            candPt->points.push_back(pset->points[i]);
    }

    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";

    std::vector<Partition*> partitionSet;
    std::vector<int> considered_halfset;
    std::vector<hyperplane*> hyperSet;
    std::vector<choose_item*> chooseitemSet;
    construct_partitions(candPt, topSet, partitionSet, considered_halfset);

    gettimeofday(&t2, 0);
    time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";


    for(int i = 0; i < pset->points.size(); ++i)
    {
        if (pset->points[i]->topk == 1)
        {
            for (int j = i + 1; j < pset->points.size(); ++j)
            {
                if (pset->points[j]->topk == 1)
                    hyperSet.push_back(new hyperplane(pset->points[i], pset->points[j]));
            }
        }
    }
    int hyperIndex = build_choose_item_table(partitionSet, hyperSet, chooseitemSet);


    double v1 = chooseitemSet[hyperIndex]->h->p_1->dot_product(u);
    double v2 = chooseitemSet[hyperIndex]->h->p_2->dot_product(u);

    //initial
    Partition *R = new Partition(dim);
    hyperplane *hy;
    int numOfQuestion = 0;
    point_t* point_result = NULL;
    while (considered_halfset.size() > 1)
    {
        numOfQuestion++;
        if (v1 >= v2)
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_2, chooseitemSet[hyperIndex]->h->p_1);
            //index = modify_choose_item_table(choose_item_set, half_set_set, considered_half_set, index, true);
        }
        else
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_1, chooseitemSet[hyperIndex]->h->p_2);
            //index = modify_choose_item_table(choose_item_set, half_set_set, considered_half_set, index, false);
        }
        v1 = chooseitemSet[hyperIndex]->h->p_1->dot_product(u);
        v2 = chooseitemSet[hyperIndex]->h->p_2->dot_product(u);

        //Find whether there exist point which is the topk point w.r.t any u in R
        R->hyperplanes.push_back(hy);
        R->set_ext_pts();
        point_set* top_current;
        top_current = R->findTopk(pset, k, s);
        if(top_current != NULL)
        {
            top_current->printResult("RH", numOfQuestion, s, t1, 0, 0);

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

    return numOfQuestion;


}



/**
 * @brief Asking user question and return one of the top-k points
 *        Find the top-1 point by sampling
 * @param p_set 		 The dataset
 * @param u 			 The linear function
 * @param k 			 The parameter
 */
int HDPI_grid1(point_set* pset, point_t *u, int k, int s)
{
    for(int i = 0; i < pset->points.size(); ++i)
        pset->points[i]->value = 0;
    timeval t1;
    gettimeofday(&t1, 0);

    std::vector<point_set*> topSet;
    int dim = pset->points[0]->dim;
    pset->findTopk_sampling(topSet, new point_t(dim), k, 0, 0);//use sampling method
    point_set *candPt = new point_set();
    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(pset->points[i]->topk == 1)
            candPt->points.push_back(pset->points[i]);
    }

    std::vector<Partition*> partitionSet;
    std::vector<int> considered_halfset;
    std::vector<hyperplane*> hyperSet;
    std::vector<choose_item*> chooseitemSet;
    construct_partitions(candPt, topSet, partitionSet, considered_halfset);
    for(int i = 0; i < topSet.size(); i +=10)
        topSet[i]->points[0]->value = 1;


    /*
    point_set *hyperSetPoint = new point_set();
    if(pset->points.size() < 100)
    {
        for (int i = 0; i < pset->points.size(); ++i)
            hyperSetPoint->points.push_back(pset->points[i]);
    }
    while (hyperSetPoint->points.size() < 100)
    {
        int index = rand() % topSet.size();
        if(topSet[index]->points[0]->value != 1)
        {
            hyperSetPoint->points.push_back(topSet[index]->points[0]);
            topSet[index]->points[0]->value = 1;
        }
    }
    */

    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";


    for(int i = 0; i < pset->points.size(); ++i)
    {
        if (pset->points[i]->value == 1)
        {
            for (int j = i + 1; j < pset->points.size(); ++j)
            {
                if (pset->points[j]->value == 1)
                    hyperSet.push_back(new hyperplane(pset->points[i], pset->points[j]));
            }
        }
    }
    int hyperIndex = build_choose_item_table(partitionSet, hyperSet, chooseitemSet);


    t2; gettimeofday(&t2, 0);
    time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";

    double v1 = chooseitemSet[hyperIndex]->h->p_1->dot_product(u);
    double v2 = chooseitemSet[hyperIndex]->h->p_2->dot_product(u);

    //initial
    Partition *R = new Partition(dim);
    hyperplane *hy;
    int numOfQuestion = 0;
    point_t* point_result = NULL;
    while (considered_halfset.size() > 1)
    {
        numOfQuestion++;
        if (v1 >= v2)
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_2, chooseitemSet[hyperIndex]->h->p_1);
            //index = modify_choose_item_table(choose_item_set, half_set_set, considered_half_set, index, true);
        }
        else
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_1, chooseitemSet[hyperIndex]->h->p_2);
            //index = modify_choose_item_table(choose_item_set, half_set_set, considered_half_set, index, false);
        }
        v1 = chooseitemSet[hyperIndex]->h->p_1->dot_product(u);
        v2 = chooseitemSet[hyperIndex]->h->p_2->dot_product(u);

        //Find whether there exist point which is the topk point w.r.t any u in R
        R->hyperplanes.push_back(hy);
        R->set_ext_pts();
        point_set* top_current;
        top_current = R->findTopk(pset, k, s);
        if(top_current != NULL)
        {
            top_current->printResult("RH", numOfQuestion, s, t1, 0, 0);

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

    return numOfQuestion;


}


/**
 * @brief Asking user question and return one of the top-k points
 *        Find the top-1 point by sampling
 * @param p_set 		 The dataset
 * @param u 			 The linear function
 * @param k 			 The parameter
 */
int HDPI_grid(point_set* pset, point_t *u, int k, int s)
{
    for(int i = 0; i < pset->points.size(); ++i)
    {
        pset->points[i]->value = 0;
        pset->points[i]->result = 0;
        pset->points[i]->count = 0;
    }
    timeval t1;
    gettimeofday(&t1, 0);

    std::vector<point_set*> topSet;
    int dim = pset->points[0]->dim;
    Partition *R = new Partition(dim);
    pset->findTopk_sampling(topSet, new point_t(dim), k, 0, 0);//use sampling method
    point_set *candPt = new point_set();
    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(pset->points[i]->topk == 1)
            candPt->points.push_back(pset->points[i]);
    }

    std::vector<Partition*> partitionSet;
    std::vector<int> considered_halfset;
    std::vector<hyperplane*> hyperSet;
    std::vector<choose_item*> chooseitemSet;
    construct_partitions(candPt, R, topSet, hyperSet, partitionSet, considered_halfset);

    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";


    int hyperIndex = build_choose_item_table(partitionSet, hyperSet, chooseitemSet);

    t2; gettimeofday(&t2, 0);
    time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";

    double v1, v2;

    //initial
    hyperplane *hy;
    int numOfQuestion = 0;
    point_t* point_result = NULL;
    while (1)
    {

        numOfQuestion++;
        v1 = chooseitemSet[hyperIndex]->h->p_1->dot_product(u);
        v2 = chooseitemSet[hyperIndex]->h->p_2->dot_product(u);
        if (v1 >= v2)
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_2, chooseitemSet[hyperIndex]->h->p_1);
            hyperIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset, hyperIndex, true);
        }
        else
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_1, chooseitemSet[hyperIndex]->h->p_2);
            hyperIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset, hyperIndex, false);
        }


        //Find whether there exist point which is the topk point w.r.t any u in R
        R->hyperplanes.push_back(hy);
        R->set_ext_pts();
        //R->print();
        point_set* top_current;
        top_current = R->findTopk(pset, k, s);
        if(top_current != NULL)
        {
            top_current->printResult("HDPI", numOfQuestion, s, t1, 0, 0);

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

        if(considered_halfset.size() < 2)
        {
            topSet.clear();
            partitionSet.clear();
            hyperSet.clear();
            chooseitemSet.clear();
            considered_halfset.clear();
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);

            for(int i = 0; i < pset->points.size(); ++i)
            {
                pset->points[i]->topk = 0;
            }

            R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method

            for(int i = 0; i < pset->points.size(); ++i)
            {
                if(pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, hyperSet, partitionSet, considered_halfset);
            if(topSet.size() <= 1)
            {
                topSet[0]->printResult("HDPI", numOfQuestion, s, t1, 0, 0);

                point_set *resultSet = new point_set();
                pset->findTopk(u, k, resultSet);
                bool resultcheck = true;
                for(int i = 0; i < s; ++i)
                {
                    if(!resultSet->checkExist(topSet[0]->points[i]))
                    {
                        resultcheck = false;
                        break;
                    }
                }
                std::cout << resultcheck <<"\n";
                return numOfQuestion;
            }
            hyperIndex = build_choose_item_table(partitionSet, hyperSet, chooseitemSet);
        }
    }

    return numOfQuestion;


















}




/**
 * @brief Asking user question and return one of the top-k points
 *        Find the top-1 point by sampling
 * @param p_set 		 The dataset
 * @param u 			 The linear function
 * @param k 			 The parameter
 */
int HDPI_extreme(point_set* pset, point_t *u, int k, int s)
{
    for(int i = 0; i < pset->points.size(); ++i)
        pset->points[i]->value = 0;
    timeval t1;
    gettimeofday(&t1, 0);

    std::vector<point_set*> topSet;
    int dim = pset->points[0]->dim;
    Partition *R = new Partition(dim);
    R->findTopk_extreme(pset, topSet, k);
    point_set *candPt = new point_set();
    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(pset->points[i]->topk == 1)
            candPt->points.push_back(pset->points[i]);
    }

    std::vector<Partition*> partitionSet;
    std::vector<int> considered_halfset;
    std::vector<hyperplane*> hyperSet;
    std::vector<choose_item*> chooseitemSet;
    construct_partitions(candPt, R, topSet, hyperSet, partitionSet, considered_halfset);

    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";


    int hyperIndex = build_choose_item_table(partitionSet, hyperSet, chooseitemSet);


    t2; gettimeofday(&t2, 0);
    time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << time_cost << "\n";

    double v1, v2;

    //initial
    hyperplane *hy;
    int numOfQuestion = 0;
    point_t* point_result = NULL;
    while (1)
    {

        numOfQuestion++;
        v1 = chooseitemSet[hyperIndex]->h->p_1->dot_product(u);
        v2 = chooseitemSet[hyperIndex]->h->p_2->dot_product(u);
        if (v1 >= v2)
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_2, chooseitemSet[hyperIndex]->h->p_1);
            hyperIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset, hyperIndex, true);
        }
        else
        {
            hy = new hyperplane(chooseitemSet[hyperIndex]->h->p_1, chooseitemSet[hyperIndex]->h->p_2);
            hyperIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset, hyperIndex, false);
        }


        //Find whether there exist point which is the topk point w.r.t any u in R
        R->hyperplanes.push_back(hy);
        R->set_ext_pts();
        //R->print();
        point_set* top_current;
        top_current = R->findTopk(pset, k, s);
        if(top_current != NULL)
        {
            top_current->printResult("HDPI", numOfQuestion, s, t1, 0, 0);

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

        if(considered_halfset.size() < 2)
        {
            topSet.clear();
            partitionSet.clear();
            hyperSet.clear();
            chooseitemSet.clear();
            considered_halfset.clear();
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);


            for(int i = 0; i < pset->points.size(); ++i)
            {
                pset->points[i]->topk = 0;
            }

            R->findTopk_extreme(pset, topSet, k);
            //R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method

            for(int i = 0; i < pset->points.size(); ++i)
            {
                if(pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, hyperSet, partitionSet, considered_halfset);
            if(topSet.size() <= 1)
            {
                topSet[0]->printResult("HDPI", numOfQuestion, s, t1, 0, 0);

                point_set *resultSet = new point_set();
                pset->findTopk(u, k, resultSet);
                bool resultcheck = true;
                for(int i = 0; i < topSet[0]->points.size(); ++i)
                {
                    if(!resultSet->checkExist(topSet[0]->points[i]))
                    {
                        resultcheck = false;
                        break;
                    }
                }
                std::cout << resultcheck <<"\n";
                return numOfQuestion;
            }
            hyperIndex = build_choose_item_table(partitionSet, hyperSet, chooseitemSet);
        }
    }

    return numOfQuestion;


















}