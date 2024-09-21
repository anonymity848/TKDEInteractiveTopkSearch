#include "Unify.h"
#include "HD/HDPI/HDPI.h"

int unify(point_set* pset, point_t *u, int k, int s, int CH, int HS)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int numOfQuestion = 0;
    int dim = pset->points[0]->dim, M = pset->points.size();
    Partition *R = new Partition(dim);
    std::vector<point_t*> considered;
    hyperplane_set *candHy = new hyperplane_set();
    point_set *top_current = NULL;


    while(1)
    {
        std::vector<point_set *> topSet;
        point_set *candPt = new point_set();
        std::vector<Partition*> partitionSet;
        std::vector<int> considered_halfset;
        std::vector<choose_item*> chooseitemSet;

        //candidate hyper-plane generation
        if (CH == 1) //skyband
        {
            for (int i = 0; i < pset->points.size() - 1; ++i)
            {
                for (int j = i + 1; j < pset->points.size(); ++j)
                {
                    hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                    candHy->hyperplanes.push_back(h);
                }
            }
        }
        else if (CH == 2) //Randomized Grouping
        {
            if(considered.size() == 0)//initialization
            {
                int cut = 8 * sqrt(k) * sqrt(dim);
                if (cut > M) cut = M;
                for (int i = 0; i < cut; ++i)
                    considered.push_back(pset->points[i]);
                for (int i = 0; i < cut - 1; ++i)
                {
                    for (int j = i + 1; j < cut; ++j)
                    {
                        hyperplane *h = new hyperplane(considered[i], considered[j]);
                        int relation = R->check_relation(h);
                        if (relation == 0)
                            candHy->hyperplanes.push_back(h);
                    }
                }
            }
            else
            {
                int cIndex = considered.size();
                for(int i = 0; i < considered.size(); ++i)
                {
                    hyperplane *h = new hyperplane(considered[i], pset->points[cIndex]);
                    int relation = R->check_relation(h);
                    if (relation == 0)
                        candHy->hyperplanes.push_back(h);
                }
                considered.push_back(pset->points[cIndex]);
            }
        }
        else if (CH == 3) //Iterative Sampling
        {
            int dim = pset->points[0]->dim;
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);
            for(int i = 0; i < pset->points.size(); ++i)
                pset->points[i]->topk = 0;
            R->findTopk_extreme(pset, topSet, k);
            R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, candHy->hyperplanes, partitionSet, considered_halfset);
            /*
            for(int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                {
                    for (int j = i + 1; j < pset->points.size(); ++j)
                    {
                        if (pset->points[j]->topk == 1)
                            candHy->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[j]));
                    }
                }
            }
            */
            if(partitionSet.size() <= 1)
            {
                top_current = partitionSet[0]->topSet;
                if (top_current != NULL)
                {
                    char algName[100];
                    sprintf(algName, "Unify%d%d", CH, HS);
                    //correctness checking
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
                        testsum += top_current->points[i]->dot_product(u);
                    }
                    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                    top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                    return numOfQuestion;
                }
            }

        }
        else if(CH == 4 && s == 1) //top1
        {
            pset->find_convexHull(topSet);
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, candHy->hyperplanes, partitionSet, considered_halfset);
        }

        std::cout << candHy->hyperplanes.size() << " candidate hyperplanes. \n";
        std::cout << "CH finished. \n";

        int candSize = candHy->hyperplanes.size();
        int hIndex = -200;
        bool direction = true;
        while (candSize > 0)
        {

            //hyper-plane selection
            if (HS == 1)//Covering score
            {
                double maxScore = -1;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    //positive
                    int postive_num = 0;
                    Partition *testR_positive = new Partition(R);
                    testR_positive->hyperplanes.push_back(
                            new hyperplane(candHy->hyperplanes[i]->p_1, candHy->hyperplanes[i]->p_2));
                    testR_positive->set_ext_pts();
                    if(testR_positive->ext_pts.size() <= 0)
                        continue;
                    testR_positive->center = testR_positive->average_point();
                    for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                    {
                        if (i != j)
                        {
                            hyperplane *h = candHy->hyperplanes[j];
                            int relation = testR_positive->check_relation(h);
                            if (relation != 0)
                            {
                                postive_num++;
                            }
                        }
                    }

                    //negative
                    int negative_num = 0;
                    if (postive_num > maxScore)
                    {
                        Partition *testR_negative = new Partition(R);
                        testR_negative->hyperplanes.push_back(
                                new hyperplane(candHy->hyperplanes[i]->p_2, candHy->hyperplanes[i]->p_1));
                        testR_negative->set_ext_pts();
                        if(testR_negative->ext_pts.size() <= 0)
                            continue;
                        testR_negative->center = testR_negative->average_point();
                        for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                        {
                            if (i != j)
                            {
                                hyperplane *h = candHy->hyperplanes[j];
                                int relation = testR_negative->check_relation(h);
                                if (relation != 0)
                                {
                                    negative_num++;
                                }
                            }
                        }
                    }

                    if (std::min(postive_num, negative_num) > maxScore)
                    {
                        maxScore = std::min(postive_num, negative_num);
                        hIndex = i;
                    }
                }
            }
            else if (HS == 2)
            {
                if(partitionSet.size() <=0)
                {
                    int dim = pset->points[0]->dim;
                    double *min = new double[dim], *max = new double[dim];
                    R->findMinMax(min, max);
                    for(int i = 0; i < pset->points.size(); ++i)
                        pset->points[i]->topk = 0;
                    R->findTopk_extreme(pset, topSet, k);
                    R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
                    for (int i = 0; i < pset->points.size(); ++i)
                    {
                        if (pset->points[i]->topk == 1)
                            candPt->points.push_back(pset->points[i]);
                    }
                    hyperplane_set *candHy_test = new hyperplane_set();
                    construct_partitions(candPt, R, topSet, candHy_test->hyperplanes, partitionSet, considered_halfset);
                }
                if(hIndex == -200)
                    hIndex = build_choose_item_table(partitionSet, candHy->hyperplanes, chooseitemSet);
                else
                    hIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset,
                                                      hIndex, direction);
                if(hIndex == -1)
                    break;
                for(int i = 0; i < considered_halfset.size(); ++i)
                    std::cout << partitionSet[considered_halfset[i]]->ID << "  " << partitionSet[considered_halfset[i]]->topSet->points[0]->id << "\n";
            }
            else if (HS == 3)
            {
                double minDis = 100000000;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    hyperplane *h = candHy->hyperplanes[i];
                    double dst_h = h->check_distance(R->center);
                    if (dst_h < minDis)
                    {
                        minDis = dst_h;
                        hIndex = i;
                    }
                }
            }


            if(HS != 2 || chooseitemSet.size() >= 1 && considered_halfset.size() >= 2)
            {
                //interaction
                numOfQuestion++;
                point_t *p1, *p2;
                if(HS != 2)
                {
                    p1 = candHy->hyperplanes[hIndex]->p_1;
                    p2 = candHy->hyperplanes[hIndex]->p_2;
                } else
                {
                    p1 = chooseitemSet[hIndex]->h->p_1;
                    p2 = chooseitemSet[hIndex]->h->p_2;
                }
                double v1 = p1->dot_product(u);
                double v2 = p2->dot_product(u);
                if (v1 > v2)
                {
                    hyperplane *half = new hyperplane(p2, p1);
                    R->hyperplanes.push_back(half);
                    direction = true;
                } else
                {
                    hyperplane *half = new hyperplane(p1, p2);
                    R->hyperplanes.push_back(half);
                    direction = false;
                }
                if(HS != 2)
                    candHy->hyperplanes.erase(candHy->hyperplanes.begin() + hIndex);
                R->set_ext_pts();
                R->center = R->average_point();

                //stopping condition
                top_current = R->findTopk(pset, k, s);
                if (top_current != NULL)
                {
                    char algName[100];
                    sprintf(algName, "Unify%d%d", CH, HS);
                    //correctness checking
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
                        testsum += top_current->points[i]->dot_product(u);
                    }
                    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                    top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                    return numOfQuestion;
                }
                candHy->hyperplanes.erase(
                        std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                       [R](hyperplane *h)
                                       { return R->check_relation(h) != 0; }),
                        candHy->hyperplanes.end());
            }

            if(HS == 2)
            {
                candSize = std::min(chooseitemSet.size(), considered_halfset.size() - 1);
                std::cout << "cand partition: " << candSize + 1 << "   chooseitem: " << chooseitemSet.size() << std::endl;
            }
            else
            {
                candSize = candHy->hyperplanes.size();
                std::cout << "cand hyperplane: " << candSize << std::endl;
            }



        }

    }

}

/**
 * @brief Several questions accompanied by a single test question
 * @param pset
 * @param u
 * @param k
 * @param s
 * @param CH
 * @param HS
 * @param error
 * @return
 */
int unify1(point_set* pset, point_t *u, int k, int s, int CH, int HS, int error)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int numOfQuestion = 0, normalQuestion = 0;
    int dim = pset->points[0]->dim, M = pset->points.size();
    Partition *R = new Partition(dim), *RROrigin = new Partition(dim);
    std::vector<point_t*> considered;
    hyperplane_set *candHy = new hyperplane_set();
    std::vector<hyperplane*> potentialHyper, usedHyper;
    point_set *top_current = NULL;
    int testforeachQ = 1;

    while(1)
    {
        std::vector<point_set *> topSet;
        point_set *candPt = new point_set();
        std::vector<Partition*> partitionSet;
        std::vector<int> considered_halfset;
        std::vector<choose_item*> chooseitemSet;

        //candidate hyper-plane generation
        if (CH == 1) //skyband
        {
            for (int i = 0; i < pset->points.size() - 1; ++i)
            {
                for (int j = i + 1; j < pset->points.size(); ++j)
                {
                    hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                    candHy->hyperplanes.push_back(h);
                }
            }
        }
        else if (CH == 2) //Randomized Grouping
        {
            if(considered.size() == 0)//initialization
            {
                int cut = 8 * sqrt(k) * sqrt(dim);
                if (cut > M) cut = M;
                for (int i = 0; i < cut; ++i)
                    considered.push_back(pset->points[i]);
                for (int i = 0; i < cut - 1; ++i)
                {
                    for (int j = i + 1; j < cut; ++j)
                    {
                        hyperplane *h = new hyperplane(considered[i], considered[j]);
                        int relation = R->check_relation(h);
                        if (relation == 0)
                            candHy->hyperplanes.push_back(h);
                    }
                }
            }
            else
            {
                int cIndex = considered.size();
                for(int i = 0; i < considered.size(); ++i)
                {
                    hyperplane *h = new hyperplane(considered[i], pset->points[cIndex]);
                    int relation = R->check_relation(h);
                    if (relation == 0)
                        candHy->hyperplanes.push_back(h);
                }
                considered.push_back(pset->points[cIndex]);
            }
        }
        else if (CH == 3) //Iterative Sampling
        {
            int dim = pset->points[0]->dim;
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);
            for(int i = 0; i < pset->points.size(); ++i)
                pset->points[i]->topk = 0;
            R->findTopk_extreme(pset, topSet, k);
            R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, candHy->hyperplanes, partitionSet, considered_halfset);
            /*
            for(int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                {
                    for (int j = i + 1; j < pset->points.size(); ++j)
                    {
                        if (pset->points[j]->topk == 1)
                            candHy->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[j]));
                    }
                }
            }
            */
            if(partitionSet.size() <= 1)
            {
                top_current = partitionSet[0]->topSet;
                if (top_current != NULL)
                {
                    char algName[100];
                    sprintf(algName, "Unify%d%d", CH, HS);
                    //correctness checking
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
                        testsum += top_current->points[i]->dot_product(u);
                    }
                    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                    top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                    return numOfQuestion;
                }
            }

        }
        std::cout << "CH finished. \n";



        for(int i = 0; i < candHy->hyperplanes.size(); ++i)
            usedHyper.push_back(candHy->hyperplanes[i]);

        int candSize = candHy->hyperplanes.size();
        int hIndex = -200;
        bool direction = true;
        while (candSize > 0)
        {
            //hyper-plane selection
            if (HS == 1)//Covering score
            {
                double maxScore = -1;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    //positive
                    int postive_num = 0;
                    Partition *testR_positive = new Partition(R);
                    testR_positive->hyperplanes.push_back(
                            new hyperplane(candHy->hyperplanes[i]->p_1, candHy->hyperplanes[i]->p_2));
                    testR_positive->set_ext_pts();
                    if(testR_positive->ext_pts.size() <= 0)
                        continue;
                    testR_positive->center = testR_positive->average_point();
                    for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                    {
                        if (i != j)
                        {
                            hyperplane *h = candHy->hyperplanes[j];
                            int relation = testR_positive->check_relation(h);
                            if (relation != 0)
                            {
                                postive_num++;
                            }
                        }
                    }

                    //negative
                    int negative_num = 0;
                    if (postive_num > maxScore)
                    {
                        Partition *testR_negative = new Partition(R);
                        testR_negative->hyperplanes.push_back(
                                new hyperplane(candHy->hyperplanes[i]->p_2, candHy->hyperplanes[i]->p_1));
                        testR_negative->set_ext_pts();
                        if(testR_negative->ext_pts.size() <= 0)
                            continue;
                        testR_negative->center = testR_negative->average_point();
                        for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                        {
                            if (i != j)
                            {
                                hyperplane *h = candHy->hyperplanes[j];
                                int relation = testR_negative->check_relation(h);
                                if (relation != 0)
                                {
                                    negative_num++;
                                }
                            }
                        }
                    }

                    if (std::min(postive_num, negative_num) > maxScore)
                    {
                        maxScore = std::min(postive_num, negative_num);
                        hIndex = i;
                    }
                }
            }
            else if (HS == 2)
            {
                if(partitionSet.size() <=0)
                {
                    int dim = pset->points[0]->dim;
                    double *min = new double[dim], *max = new double[dim];
                    R->findMinMax(min, max);
                    for(int i = 0; i < pset->points.size(); ++i)
                        pset->points[i]->topk = 0;
                    R->findTopk_extreme(pset, topSet, k);
                    R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
                    for (int i = 0; i < pset->points.size(); ++i)
                    {
                        if (pset->points[i]->topk == 1)
                            candPt->points.push_back(pset->points[i]);
                    }
                    hyperplane_set *candHy_test = new hyperplane_set();
                    construct_partitions(candPt, R, topSet, candHy_test->hyperplanes, partitionSet, considered_halfset);
                }
                if(hIndex == -200)
                    hIndex = build_choose_item_table(partitionSet, candHy->hyperplanes, chooseitemSet);
                else
                    hIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset,
                                                      hIndex, direction);
                if(hIndex == -1)
                    break;
            }
            else if (HS == 3)
            {
                double minDis = 100000000;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    hyperplane *h = candHy->hyperplanes[i];
                    double dst_h = h->check_distance(R->center);
                    if (dst_h < minDis)
                    {
                        minDis = dst_h;
                        hIndex = i;
                    }
                }
            }


            //normal interaction
            if(HS != 2 || chooseitemSet.size() >= 1 && considered_halfset.size() >= 2)
            {
                numOfQuestion++;
                normalQuestion++;
                point_t *p1, *p2;
                if(HS != 2)
                {
                    p1 = candHy->hyperplanes[hIndex]->p_1;
                    p2 = candHy->hyperplanes[hIndex]->p_2;
                }
                else
                {
                    p1 = chooseitemSet[hIndex]->h->p_1;
                    p2 = chooseitemSet[hIndex]->h->p_2;
                }
                double v1 = p1->dot_product(u);
                double v2 = p2->dot_product(u);
                int compareResult = 0;
                for(int i = 0; i < testforeachQ; ++i)
                    compareResult += u->compare(p1, p2, error);
                if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                    std::cout << "correct" << std::endl;
                else
                    std::cout << "wrong" << std::endl;

                if (compareResult > 0)
                {
                    hyperplane *half = new hyperplane(p2, p1);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = true;
                }
                else
                {
                    hyperplane *half = new hyperplane(p1, p2);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = false;
                }
                if(HS != 2)
                    candHy->hyperplanes.erase(candHy->hyperplanes.begin() + hIndex);
                R->set_ext_pts_withoutHull();
                R->center = R->average_point();

                candHy->hyperplanes.erase(
                        std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                       [R](hyperplane *h)
                                       { return R->check_relation(h) != 0; }),
                        candHy->hyperplanes.end());
            }

            //print the current status
            if(HS == 2)
            {
                candSize = std::min(chooseitemSet.size(), considered_halfset.size() - 1);
                std::cout << "cand partition: " << candSize << "   chooseitem: " << chooseitemSet.size()
                          << std::endl;
            }
            else
            {
                candSize = candHy->hyperplanes.size();
                std::cout << "cand hyperplane: " << candSize << std::endl;
            }


            //stopping condition
            top_current = R->findTopk(pset, k, s);
            if (top_current != NULL)
            {
                char algName[100];
                sprintf(algName, "Unify%d%d", CH, HS);
                //correctness checking
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
                    testsum += top_current->points[i]->dot_product(u);
                }
                double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                return numOfQuestion;
            }


            //test question
            int Y = 2, testNum = 4;

            if ((normalQuestion) % Y == 0)
            {
                for(int t = 0; t < testNum; ++t)
                {
                    double maxDis = -10000;
                    int testhIndex = -1;
                    for (int i = 0; i < usedHyper.size(); ++i)
                    {
                        hyperplane *h = usedHyper[i];
                        double dst_h = h->check_distance(R->center);
                        if (R->check_relation(h) != 0 && dst_h > maxDis )
                        {
                            maxDis = dst_h;
                            testhIndex = i;
                        }
                    }

                    numOfQuestion++;
                    point_t *p1 = usedHyper[testhIndex]->p_1;
                    point_t *p2 = usedHyper[testhIndex]->p_2;
                    double v1 = p1->dot_product(u);
                    double v2 = p2->dot_product(u);
                    int compareResult = 0;
                    for(int i = 0; i < testforeachQ; ++i)
                        compareResult += u->compare(p1, p2, error);
                    if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                        std::cout << "correct" << std::endl;
                    else
                        std::cout << "wrong" << std::endl;

                    if (compareResult > 0)
                    {
                        hyperplane *half = new hyperplane(p2, p1);
                        R->hyperplanes.push_back(half);
                        potentialHyper.push_back(half);
                    } else
                    {
                        hyperplane *half = new hyperplane(p1, p2);
                        R->hyperplanes.push_back(half);
                        potentialHyper.push_back(half);
                    }
                }


                bool intersectionFound = R->set_ext_pts_withoutHull();
                if(intersectionFound && R->ext_pts.size() > 0)
                {
                    R->center = R->average_point();
                    RROrigin = new Partition(R);
                    RROrigin->center = RROrigin->average_point();
                }
                else
                {
                    //Partition *testR = new Partition(RROrigin);
                    for(int i = 0; i < potentialHyper.size(); ++i)
                    {
                        RROrigin->hyperplanes.push_back(potentialHyper[i]);
                    }
                    double epsilon = 0.85;
                    while (!intersectionFound || RROrigin->ext_pts.size() <= 0)
                    {
                        epsilon += 0.1;
                        for (int i = RROrigin->hyperplanes.size() - potentialHyper.size(); i < RROrigin->hyperplanes.size(); ++i)
                        {
                            /*
                            point_t *qminusp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qminusp->attr[j] = psmall->attr[j] - plarge->attr[j];
                            point_t *qaddp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qaddp->attr[j] = psmall->attr[j] + plarge->attr[j];
                            double off1 = testR->max_uq(RROrigin->hyperplanes[i], psmall);
                            off1 = qminusp->dot_product(qminusp) * (2 * epsilon - 1) * off1;
                            double off2 = qminusp->dot_product(psmall) - epsilon * qminusp->dot_product(qaddp);
                            RROrigin->hyperplanes[i]->offset = off1 / off2;
                             */
                            RROrigin->hyperplanes[i]->offset =  - log(epsilon / (1.0 - epsilon)); //- 1.0 / QuestionHyper[i]->p_1->distance(QuestionHyper[i]->p_2) *log(epsilon / (1.0 - epsilon));
                        }
                        intersectionFound = RROrigin->set_ext_pts_withoutHull();
                    }
                    R = new Partition(RROrigin);
                    R->center = R->average_point();
                    R->print();
                    std::cout << "IN: " << R->isIn(u) <<std::endl;

                    candHy->hyperplanes.clear();
                    for (int i = 0; i < usedHyper.size(); ++i)
                        candHy->hyperplanes.push_back(usedHyper[i]);
                    candHy->hyperplanes.erase(
                            std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                           [R](hyperplane *h)
                                           { return R->check_relation(h) != 0; }),
                            candHy->hyperplanes.end());

                    topSet.clear();
                    candPt = new point_set();
                    partitionSet.clear();
                    considered_halfset.clear();
                    chooseitemSet.clear();
                    hIndex = -200;
                }
                potentialHyper.clear();
            }


        }

    }


}


/**
 * @brief A few questions with several confirmation questions (the same questions)
 * @param pset
 * @param u
 * @param k
 * @param s
 * @param CH
 * @param HS
 * @param error
 * @return
 */
int unify2(point_set* pset, point_t *u, int k, int s, int CH, int HS, int error)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int numOfQuestion = 0, normalQuestion = 0;
    int dim = pset->points[0]->dim, M = pset->points.size();
    Partition *R = new Partition(dim), *RROrigin = new Partition(dim);
    std::vector<point_t*> considered;
    hyperplane_set *candHy = new hyperplane_set();
    std::vector<hyperplane*> potentialHyper, usedHyper;
    point_set *top_current = NULL;

    while(1)
    {
        std::vector<point_set *> topSet;
        point_set *candPt = new point_set();
        std::vector<Partition*> partitionSet;
        std::vector<int> considered_halfset;
        std::vector<choose_item*> chooseitemSet;

        //candidate hyper-plane generation
        if (CH == 1) //skyband
        {
            for (int i = 0; i < pset->points.size() - 1; ++i)
            {
                for (int j = i + 1; j < pset->points.size(); ++j)
                {
                    hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                    candHy->hyperplanes.push_back(h);
                }
            }
        }
        else if (CH == 2) //Randomized Grouping
        {
            if(considered.size() == 0)//initialization
            {
                int cut = 8 * sqrt(k) * sqrt(dim);
                if (cut > M) cut = M;
                for (int i = 0; i < cut; ++i)
                    considered.push_back(pset->points[i]);
                for (int i = 0; i < cut - 1; ++i)
                {
                    for (int j = i + 1; j < cut; ++j)
                    {
                        hyperplane *h = new hyperplane(considered[i], considered[j]);
                        int relation = R->check_relation(h);
                        if (relation == 0)
                            candHy->hyperplanes.push_back(h);
                    }
                }
            }
            else
            {
                int cIndex = considered.size();
                for(int i = 0; i < considered.size(); ++i)
                {
                    hyperplane *h = new hyperplane(considered[i], pset->points[cIndex]);
                    int relation = R->check_relation(h);
                    if (relation == 0)
                        candHy->hyperplanes.push_back(h);
                }
                considered.push_back(pset->points[cIndex]);
            }
        }
        else if (CH == 3) //Iterative Sampling
        {
            int dim = pset->points[0]->dim;
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);
            for(int i = 0; i < pset->points.size(); ++i)
                pset->points[i]->topk = 0;
            R->findTopk_extreme(pset, topSet, k);
            R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, candHy->hyperplanes, partitionSet, considered_halfset);
            /*
            for(int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                {
                    for (int j = i + 1; j < pset->points.size(); ++j)
                    {
                        if (pset->points[j]->topk == 1)
                            candHy->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[j]));
                    }
                }
            }
            */
            if(partitionSet.size() <= 1)
            {
                top_current = partitionSet[0]->topSet;
                if (top_current != NULL)
                {
                    char algName[100];
                    sprintf(algName, "Unify%d%d", CH, HS);
                    //correctness checking
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
                        testsum += top_current->points[i]->dot_product(u);
                    }
                    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                    top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                    return numOfQuestion;
                }
            }

        }
        std::cout << "CH finished. \n";



        for(int i = 0; i < candHy->hyperplanes.size(); ++i)
            usedHyper.push_back(candHy->hyperplanes[i]);

        int candSize = candHy->hyperplanes.size();
        int hIndex = -200;
        bool direction = true;
        while (candSize > 0)
        {
            //hyper-plane selection
            if (HS == 1)//Covering score
            {
                double maxScore = -1;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    //positive
                    int postive_num = 0;
                    Partition *testR_positive = new Partition(R);
                    testR_positive->hyperplanes.push_back(
                            new hyperplane(candHy->hyperplanes[i]->p_1, candHy->hyperplanes[i]->p_2));
                    testR_positive->set_ext_pts();
                    if(testR_positive->ext_pts.size() <= 0)
                        continue;
                    testR_positive->center = testR_positive->average_point();
                    for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                    {
                        if (i != j)
                        {
                            hyperplane *h = candHy->hyperplanes[j];
                            int relation = testR_positive->check_relation(h);
                            if (relation != 0)
                            {
                                postive_num++;
                            }
                        }
                    }

                    //negative
                    int negative_num = 0;
                    if (postive_num > maxScore)
                    {
                        Partition *testR_negative = new Partition(R);
                        testR_negative->hyperplanes.push_back(
                                new hyperplane(candHy->hyperplanes[i]->p_2, candHy->hyperplanes[i]->p_1));
                        testR_negative->set_ext_pts();
                        if(testR_negative->ext_pts.size() <= 0)
                            continue;
                        testR_negative->center = testR_negative->average_point();
                        for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                        {
                            if (i != j)
                            {
                                hyperplane *h = candHy->hyperplanes[j];
                                int relation = testR_negative->check_relation(h);
                                if (relation != 0)
                                {
                                    negative_num++;
                                }
                            }
                        }
                    }

                    if (std::min(postive_num, negative_num) > maxScore)
                    {
                        maxScore = std::min(postive_num, negative_num);
                        hIndex = i;
                    }
                }
            }
            else if (HS == 2)
            {
                if(partitionSet.size() <=0)
                {
                    int dim = pset->points[0]->dim;
                    double *min = new double[dim], *max = new double[dim];
                    R->findMinMax(min, max);
                    for(int i = 0; i < pset->points.size(); ++i)
                        pset->points[i]->topk = 0;
                    R->findTopk_extreme(pset, topSet, k);
                    R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
                    for (int i = 0; i < pset->points.size(); ++i)
                    {
                        if (pset->points[i]->topk == 1)
                            candPt->points.push_back(pset->points[i]);
                    }
                    hyperplane_set *candHy_test = new hyperplane_set();
                    construct_partitions(candPt, R, topSet, candHy_test->hyperplanes, partitionSet, considered_halfset);
                }
                if(hIndex == -200)
                    hIndex = build_choose_item_table(partitionSet, candHy->hyperplanes, chooseitemSet);
                else
                    hIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset,
                                                      hIndex, direction);
                if(hIndex == -1)
                    break;
            }
            else if (HS == 3)
            {
                double minDis = 100000000;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    hyperplane *h = candHy->hyperplanes[i];
                    double dst_h = h->check_distance(R->center);
                    if (dst_h < minDis)
                    {
                        minDis = dst_h;
                        hIndex = i;
                    }
                }
            }


            //normal interaction
            if(HS != 2 || chooseitemSet.size() >= 1 && considered_halfset.size() >= 2)
            {
                numOfQuestion++;
                normalQuestion++;
                point_t *p1, *p2;
                if(HS != 2)
                {
                    p1 = candHy->hyperplanes[hIndex]->p_1;
                    p2 = candHy->hyperplanes[hIndex]->p_2;
                }
                else
                {
                    p1 = chooseitemSet[hIndex]->h->p_1;
                    p2 = chooseitemSet[hIndex]->h->p_2;
                }
                double v1 = p1->dot_product(u);
                double v2 = p2->dot_product(u);
                int compareResult = u->compare(p1, p2, error);
                if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                    std::cout << "v1: " << v1 << "  v2: " << v2 << "  correct" << std::endl;
                else
                    std::cout << "v1: " << v1 << "  v2: " << v2 << "  wrong" << std::endl;
                if (compareResult > 0)
                {
                    hyperplane *half = new hyperplane(p2, p1);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = true;
                }
                else
                {
                    hyperplane *half = new hyperplane(p1, p2);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = false;
                }
                if(HS != 2)
                    candHy->hyperplanes.erase(candHy->hyperplanes.begin() + hIndex);
                R->set_ext_pts_withoutHull();
                R->center = R->average_point();

                candHy->hyperplanes.erase(
                        std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                       [R](hyperplane *h)
                                       { return R->check_relation(h) != 0; }),
                        candHy->hyperplanes.end());
            }

            if(HS == 2)
            {
                candSize = std::min(chooseitemSet.size(), considered_halfset.size() - 1);
                std::cout << "cand partition: " << candSize << "   chooseitem: " << chooseitemSet.size()
                          << std::endl;
            }
            else
            {
                candSize = candHy->hyperplanes.size();
                std::cout << "cand hyperplane: " << candSize << std::endl;
            }


            //stopping condition
            top_current = R->findTopk(pset, k, s);

            //R->print();
            //test question
            int Y = 4;
            if ((normalQuestion) % Y == 0 || top_current != NULL)
            {
                int potentialHyperSize = potentialHyper.size();
                for(int t = 0; t < potentialHyperSize; ++t)
                {
                    hyperplane *testh = potentialHyper[t];
                    numOfQuestion++;
                    point_t *p1 = testh->p_1;
                    point_t *p2 = testh->p_2;
                    double v1 = p1->dot_product(u);
                    double v2 = p2->dot_product(u);
                    int compareResult = u->compare(p1, p2, error);
                    if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                        std::cout << "v1: " << v1 << "  v2: " << v2 << "  correct" << std::endl;
                    else
                        std::cout << "v1: " << v1 << "  v2: " << v2 << "  wrong" << std::endl;
                    if (compareResult > 0)
                    {
                        hyperplane *half = new hyperplane(p2, p1);
                        R->hyperplanes.push_back(half);
                        //potentialHyper.push_back(half);
                    }
                    else
                    {
                        hyperplane *half = new hyperplane(p1, p2);
                        R->hyperplanes.push_back(half);
                        //potentialHyper.push_back(half);
                    }
                }


                bool intersectionFound = R->set_ext_pts_withoutHull();
                if(intersectionFound && R->ext_pts.size() > 0)
                {
                    R->center = R->average_point();
                    RROrigin = new Partition(R);
                    RROrigin->center = RROrigin->average_point();
                    //RROrigin->print();
                    std::cout << "IN: " << R->isIn(u) <<std::endl;

                }
                else
                {
                    for(int i = 0; i < potentialHyper.size(); ++i)
                    {
                        RROrigin->hyperplanes.push_back(new hyperplane(potentialHyper[i]));
                    }
                    double epsilon = 0.95;
                    while (!intersectionFound || RROrigin->ext_pts.size() <= 0)
                    {
                        epsilon += (1 - epsilon) / 2;
                        for (int i = RROrigin->hyperplanes.size() - potentialHyper.size(); i < RROrigin->hyperplanes.size(); ++i)
                        {
                            /*
                            point_t *qminusp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qminusp->attr[j] = psmall->attr[j] - plarge->attr[j];
                            point_t *qaddp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qaddp->attr[j] = psmall->attr[j] + plarge->attr[j];
                            double off1 = testR->max_uq(RROrigin->hyperplanes[i], psmall);
                            off1 = qminusp->dot_product(qminusp) * (2 * epsilon - 1) * off1;
                            double off2 = qminusp->dot_product(psmall) - epsilon * qminusp->dot_product(qaddp);
                            RROrigin->hyperplanes[i]->offset = off1 / off2;
                             */
                            RROrigin->hyperplanes[i]->offset =  - log(epsilon / (1.0 - epsilon)); //- 1.0 / QuestionHyper[i]->p_1->distance(QuestionHyper[i]->p_2) *log(epsilon / (1.0 - epsilon));
                        }
                        intersectionFound = RROrigin->set_ext_pts_withoutHull();
                    }
                    R = new Partition(RROrigin);
                    R->center = R->average_point();
                    R->print();
                    std::cout << "IN: " << R->isIn(u) <<std::endl;

                    /*
                    //re-ask the questions
                    int *answerCollection = new int[potentialHyper.size()];
                    for(int i = 0; i < potentialHyper.size(); ++i)
                        answerCollection[i] = -1;
                    Partition *testR;

                    int counttest = -1;
                    do{
                        counttest++;
                        if(counttest > 10)
                            testR = new Partition(dim);
                        else
                            testR = new Partition(RROrigin);
                        for(int i = 0; i < potentialHyper.size(); ++i)
                        {
                            numOfQuestion += 10;
                            point_t *p1 = potentialHyper[i]->p_1;
                            point_t *p2 = potentialHyper[i]->p_2;
                            double v1 = p1->dot_product(u);
                            double v2 = p2->dot_product(u);
                            int compareResult = u->compare(p1, p2, error);
                            if (compareResult > 0)
                                answerCollection[i]++;
                            else
                                answerCollection[i]--;
                            while(answerCollection[i] == 0)
                            {
                                numOfQuestion++;
                                p1 = potentialHyper[i]->p_1;
                                p2 = potentialHyper[i]->p_2;
                                v1 = p1->dot_product(u);
                                v2 = p2->dot_product(u);
                                compareResult = 0;
                                for(int i = 0; i < 10; ++i)
                                {
                                    compareResult = u->compare(p1, p2, error);
                                    if (compareResult > 0)
                                        answerCollection[i]++;
                                    else
                                        answerCollection[i]--;
                                }
                            }
                        }
                        for(int i = 0; i < potentialHyper.size(); ++i)
                        {
                            if(answerCollection[i] == -1)
                            {
                                hyperplane *half = new hyperplane(potentialHyper[i]);
                                testR->hyperplanes.push_back(half);
                            }
                            else
                            {
                                hyperplane *half = new hyperplane(potentialHyper[i]);
                                for(int j = 0; j < dim; ++j)
                                    half->norm[j] = -half->norm[j];
                                testR->hyperplanes.push_back(half);
                            }
                        }
                        intersectionFound = testR->set_ext_pts_withoutHull();
                    }while(!intersectionFound);
                    R = new Partition(testR); RROrigin = new Partition(testR);
                    R->center = R->average_point(); RROrigin->center = RROrigin->average_point();
                    R->print();
                    */
                    candHy->hyperplanes.clear();
                    for (int i = 0; i < usedHyper.size(); ++i)
                        candHy->hyperplanes.push_back(usedHyper[i]);
                    candHy->hyperplanes.erase(
                            std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                           [R](hyperplane *h)
                                           { return R->check_relation(h) != 0; }),
                            candHy->hyperplanes.end());

                    topSet.clear();
                    candPt = new point_set();
                    partitionSet.clear();
                    considered_halfset.clear();
                    chooseitemSet.clear();
                    hIndex = -200;
                }
                potentialHyper.clear();
            }

            top_current = R->findTopk(pset, k, s);
            if (top_current != NULL)
            {
                char algName[100];
                sprintf(algName, "Unify%d%d", CH, HS);
                //correctness checking
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
                    testsum += top_current->points[i]->dot_product(u);
                }
                double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                return numOfQuestion;
            }

        }

    }


}


/**
 * @brief A few questions with several confirmation questions (randomly generated)
 *        Use the confirmation questions to build the convex hull
 * @param pset
 * @param u
 * @param k
 * @param s
 * @param CH
 * @param HS
 * @param error
 * @return
 */
int unify3(point_set* pset, point_t *u, int k, int s, int CH, int HS, int error)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int numOfQuestion = 0, normalQuestion = 0;
    int dim = pset->points[0]->dim, M = pset->points.size();
    Partition *R = new Partition(dim), *RROrigin = new Partition(dim);
    std::vector<point_t*> considered;
    hyperplane_set *candHy = new hyperplane_set();
    std::vector<hyperplane*> potentialHyper, usedHyper;
    point_set *top_current = NULL;

    while(1)
    {
        std::vector<point_set *> topSet;
        point_set *candPt = new point_set();
        std::vector<Partition*> partitionSet;
        std::vector<int> considered_halfset;
        std::vector<choose_item*> chooseitemSet;

        //candidate hyper-plane generation
        if (CH == 1) //skyband
        {
            for (int i = 0; i < pset->points.size() - 1; ++i)
            {
                for (int j = i + 1; j < pset->points.size(); ++j)
                {
                    hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                    candHy->hyperplanes.push_back(h);
                }
            }
        }
        else if (CH == 2) //Randomized Grouping
        {
            if(considered.size() == 0)//initialization
            {
                int cut = 8 * sqrt(k) * sqrt(dim);
                if (cut > M) cut = M;
                for (int i = 0; i < cut; ++i)
                    considered.push_back(pset->points[i]);
                for (int i = 0; i < cut - 1; ++i)
                {
                    for (int j = i + 1; j < cut; ++j)
                    {
                        hyperplane *h = new hyperplane(considered[i], considered[j]);
                        int relation = R->check_relation(h);
                        if (relation == 0)
                            candHy->hyperplanes.push_back(h);
                    }
                }
            }
            else
            {
                int cIndex = considered.size();
                for(int i = 0; i < considered.size(); ++i)
                {
                    hyperplane *h = new hyperplane(considered[i], pset->points[cIndex]);
                    int relation = R->check_relation(h);
                    if (relation == 0)
                        candHy->hyperplanes.push_back(h);
                }
                considered.push_back(pset->points[cIndex]);
            }
        }
        else if (CH == 3) //Iterative Sampling
        {
            int dim = pset->points[0]->dim;
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);
            for(int i = 0; i < pset->points.size(); ++i)
                pset->points[i]->topk = 0;
            R->findTopk_extreme(pset, topSet, k);
            R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, candHy->hyperplanes, partitionSet, considered_halfset);
            /*
            for(int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                {
                    for (int j = i + 1; j < pset->points.size(); ++j)
                    {
                        if (pset->points[j]->topk == 1)
                            candHy->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[j]));
                    }
                }
            }
            */
            if(partitionSet.size() <= 1)
            {
                top_current = partitionSet[0]->topSet;
                if (top_current != NULL)
                {
                    char algName[100];
                    sprintf(algName, "Unify%d%d", CH, HS);
                    //correctness checking
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
                        testsum += top_current->points[i]->dot_product(u);
                    }
                    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                    top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                    return numOfQuestion;
                }
            }

        }
        std::cout << "CH finished. \n";



        for(int i = 0; i < candHy->hyperplanes.size(); ++i)
            usedHyper.push_back(candHy->hyperplanes[i]);

        int candSize = candHy->hyperplanes.size();
        int hIndex = -200;
        bool direction = true;
        while (candSize > 0)
        {
            //hyper-plane selection
            if (HS == 1)//Covering score
            {
                double maxScore = -1;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    //positive
                    int postive_num = 0;
                    Partition *testR_positive = new Partition(R);
                    testR_positive->hyperplanes.push_back(
                            new hyperplane(candHy->hyperplanes[i]->p_1, candHy->hyperplanes[i]->p_2));
                    testR_positive->set_ext_pts();
                    if(testR_positive->ext_pts.size() <= 0)
                        continue;
                    testR_positive->center = testR_positive->average_point();
                    for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                    {
                        if (i != j)
                        {
                            hyperplane *h = candHy->hyperplanes[j];
                            int relation = testR_positive->check_relation(h);
                            if (relation != 0)
                            {
                                postive_num++;
                            }
                        }
                    }

                    //negative
                    int negative_num = 0;
                    if (postive_num > maxScore)
                    {
                        Partition *testR_negative = new Partition(R);
                        testR_negative->hyperplanes.push_back(
                                new hyperplane(candHy->hyperplanes[i]->p_2, candHy->hyperplanes[i]->p_1));
                        testR_negative->set_ext_pts();
                        if(testR_negative->ext_pts.size() <= 0)
                            continue;
                        testR_negative->center = testR_negative->average_point();
                        for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                        {
                            if (i != j)
                            {
                                hyperplane *h = candHy->hyperplanes[j];
                                int relation = testR_negative->check_relation(h);
                                if (relation != 0)
                                {
                                    negative_num++;
                                }
                            }
                        }
                    }

                    if (std::min(postive_num, negative_num) > maxScore)
                    {
                        maxScore = std::min(postive_num, negative_num);
                        hIndex = i;
                    }
                }
            }
            else if (HS == 2)
            {
                if(partitionSet.size() <=0)
                {
                    int dim = pset->points[0]->dim;
                    double *min = new double[dim], *max = new double[dim];
                    R->findMinMax(min, max);
                    for(int i = 0; i < pset->points.size(); ++i)
                        pset->points[i]->topk = 0;
                    R->findTopk_extreme(pset, topSet, k);
                    R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
                    for (int i = 0; i < pset->points.size(); ++i)
                    {
                        if (pset->points[i]->topk == 1)
                            candPt->points.push_back(pset->points[i]);
                    }
                    hyperplane_set *candHy_test = new hyperplane_set();
                    construct_partitions(candPt, R, topSet, candHy_test->hyperplanes, partitionSet, considered_halfset);
                }
                if(hIndex == -200)
                    hIndex = build_choose_item_table(partitionSet, candHy->hyperplanes, chooseitemSet);
                else
                    hIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset,
                                                      hIndex, direction);
                if(hIndex == -1)
                    break;
            }
            else if (HS == 3)
            {
                double minDis = 100000000;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    hyperplane *h = candHy->hyperplanes[i];
                    double dst_h = h->check_distance(R->center);
                    if (dst_h < minDis)
                    {
                        minDis = dst_h;
                        hIndex = i;
                    }
                }
            }


            //normal interaction
            if(HS != 2 || chooseitemSet.size() >= 1 && considered_halfset.size() >= 2)
            {
                numOfQuestion++;
                normalQuestion++;
                point_t *p1, *p2;
                if(HS != 2)
                {
                    p1 = candHy->hyperplanes[hIndex]->p_1;
                    p2 = candHy->hyperplanes[hIndex]->p_2;
                }
                else
                {
                    p1 = chooseitemSet[hIndex]->h->p_1;
                    p2 = chooseitemSet[hIndex]->h->p_2;
                }
                double v1 = p1->dot_product(u);
                double v2 = p2->dot_product(u);
                int compareResult = u->compare(p1, p2, error);
                if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                    std::cout << "correct" << std::endl;
                else
                    std::cout << "wrong" << std::endl;
                if (compareResult > 0)
                {
                    hyperplane *half = new hyperplane(p2, p1);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = true;
                }
                else
                {
                    hyperplane *half = new hyperplane(p1, p2);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = false;
                }
                if(HS != 2)
                    candHy->hyperplanes.erase(candHy->hyperplanes.begin() + hIndex);
                R->set_ext_pts_withoutHull();
                R->center = R->average_point();

                candHy->hyperplanes.erase(
                        std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                       [R](hyperplane *h)
                                       { return R->check_relation(h) != 0; }),
                        candHy->hyperplanes.end());
            }

            if(HS == 2)
            {
                candSize = std::min(chooseitemSet.size(), considered_halfset.size() - 1);
                std::cout << "cand partition: " << candSize << "   chooseitem: " << chooseitemSet.size()
                          << std::endl;
            }
            else
            {
                candSize = candHy->hyperplanes.size();
                std::cout << "cand hyperplane: " << candSize << std::endl;
            }


            //stopping condition
            top_current = R->findTopk(pset, k, s);

            R->print();
            //test question
            int Y = 3, testNum = 3;
            if ((normalQuestion) % Y == 0 || top_current != NULL)
            {
                for (int t = 0; t < testNum; ++t)
                {
                    double minDis = 10000;
                    hyperplane *testh = NULL;
                    int countTest = 0;
                    while (countTest < 10000 || testh == NULL)
                    {
                        int i = rand() % pset->points.size();
                        int j = rand() % pset->points.size();
                        hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                        double dst_h = h->check_distance(R->center);
                        if (R->check_relation(h) != 0 && dst_h < minDis)
                        {
                            minDis = dst_h;
                            testh = h;
                        }
                        countTest++;
                    }

                    numOfQuestion++;
                    point_t *p1 = testh->p_1;
                    point_t *p2 = testh->p_2;
                    double v1 = p1->dot_product(u);
                    double v2 = p2->dot_product(u);
                    int compareResult = u->compare(p1, p2, error);
                    if (v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                        std::cout << "v1: " << v1 << "  v2: " << v2 << "  correct" << std::endl;
                    else
                        std::cout << "v1: " << v1 << "  v2: " << v2 << "  wrong" << std::endl;
                    if (compareResult > 0)
                    {
                        hyperplane *half = new hyperplane(p2, p1);
                        R->hyperplanes.push_back(half);
                        potentialHyper.push_back(half);
                    } else
                    {
                        hyperplane *half = new hyperplane(p1, p2);
                        R->hyperplanes.push_back(half);
                        potentialHyper.push_back(half);
                    }
                }


                bool intersectionFound = R->set_ext_pts_withoutHull();
                if (intersectionFound && R->ext_pts.size() > 0)
                {
                    double epsilon = 0.3;
                    for (int i = potentialHyper.size() - 3; i < potentialHyper.size(); ++i)
                    {
                        RROrigin->hyperplanes.push_back(new hyperplane(potentialHyper[i]));
                        RROrigin->hyperplanes[RROrigin->hyperplanes.size() - 1]->offset = 0; //-log(epsilon / (1.0 - epsilon));;
                    }
                    RROrigin->set_ext_pts_withoutHull();
                    R = new Partition(RROrigin);
                    R->center = R->average_point();
                    R->print();
                    std::cout << "IN: " << R->isIn(u) << std::endl;
                }
                else
                {
                    for (int i = potentialHyper.size() - 3; i < potentialHyper.size(); ++i)
                    {
                        RROrigin->hyperplanes.push_back(new hyperplane(potentialHyper[i]));
                    }
                    double epsilon = 0.98;
                    while (!intersectionFound || RROrigin->ext_pts.size() <= 0)
                    {
                        epsilon += (1 - epsilon) / 2;
                        for (int i = RROrigin->hyperplanes.size() - 3; i < RROrigin->hyperplanes.size(); ++i)
                        {
                            /*
                            point_t *qminusp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qminusp->attr[j] = psmall->attr[j] - plarge->attr[j];
                            point_t *qaddp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qaddp->attr[j] = psmall->attr[j] + plarge->attr[j];
                            double off1 = testR->max_uq(RROrigin->hyperplanes[i], psmall);
                            off1 = qminusp->dot_product(qminusp) * (2 * epsilon - 1) * off1;
                            double off2 = qminusp->dot_product(psmall) - epsilon * qminusp->dot_product(qaddp);
                            RROrigin->hyperplanes[i]->offset = off1 / off2;
                             */
                            RROrigin->hyperplanes[i]->offset = -log(epsilon / (1.0 - epsilon));
                            //- 1.0 / QuestionHyper[i]->p_1->distance(QuestionHyper[i]->p_2) * log(epsilon / (1.0 - epsilon));
                        }
                        RROrigin->print();
                        intersectionFound = RROrigin->set_ext_pts_withoutHull();
                    }
                    R = new Partition(RROrigin);
                    R->center = R->average_point();
                    R->print();
                    std::cout << "IN: " << R->isIn(u) << std::endl;
                }
                /*
                //re-ask the questions
                int *answerCollection = new int[potentialHyper.size()];
                for(int i = 0; i < potentialHyper.size(); ++i)
                    answerCollection[i] = -1;
                Partition *testR;

                int counttest = -1;
                do{
                    counttest++;
                    if(counttest > 10)
                        testR = new Partition(dim);
                    else
                        testR = new Partition(RROrigin);
                    for(int i = 0; i < potentialHyper.size(); ++i)
                    {
                        numOfQuestion += 10;
                        point_t *p1 = potentialHyper[i]->p_1;
                        point_t *p2 = potentialHyper[i]->p_2;
                        double v1 = p1->dot_product(u);
                        double v2 = p2->dot_product(u);
                        int compareResult = u->compare(p1, p2, error);
                        if (compareResult > 0)
                            answerCollection[i]++;
                        else
                            answerCollection[i]--;
                        while(answerCollection[i] == 0)
                        {
                            numOfQuestion++;
                            p1 = potentialHyper[i]->p_1;
                            p2 = potentialHyper[i]->p_2;
                            v1 = p1->dot_product(u);
                            v2 = p2->dot_product(u);
                            compareResult = 0;
                            for(int i = 0; i < 10; ++i)
                            {
                                compareResult = u->compare(p1, p2, error);
                                if (compareResult > 0)
                                    answerCollection[i]++;
                                else
                                    answerCollection[i]--;
                            }
                        }
                    }
                    for(int i = 0; i < potentialHyper.size(); ++i)
                    {
                        if(answerCollection[i] == -1)
                        {
                            hyperplane *half = new hyperplane(potentialHyper[i]);
                            testR->hyperplanes.push_back(half);
                        }
                        else
                        {
                            hyperplane *half = new hyperplane(potentialHyper[i]);
                            for(int j = 0; j < dim; ++j)
                                half->norm[j] = -half->norm[j];
                            testR->hyperplanes.push_back(half);
                        }
                    }
                    intersectionFound = testR->set_ext_pts_withoutHull();
                }while(!intersectionFound);
                R = new Partition(testR); RROrigin = new Partition(testR);
                R->center = R->average_point(); RROrigin->center = RROrigin->average_point();
                R->print();
                */
                candHy->hyperplanes.clear();
                for (int i = 0; i < usedHyper.size(); ++i)
                    candHy->hyperplanes.push_back(usedHyper[i]);
                candHy->hyperplanes.erase(
                        std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                       [R](hyperplane *h)
                                       { return R->check_relation(h) != 0; }),
                        candHy->hyperplanes.end());

                topSet.clear();
                candPt = new point_set();
                partitionSet.clear();
                considered_halfset.clear();
                chooseitemSet.clear();
                hIndex = -200;

                potentialHyper.clear();
            }

            top_current = R->findTopk(pset, k, s);
            if (top_current != NULL)
            {
                char algName[100];
                sprintf(algName, "Unify%d%d", CH, HS);
                //correctness checking
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
                    testsum += top_current->points[i]->dot_product(u);
                }
                double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                return numOfQuestion;
            }
        }

    }


}



/**
 * @brief A few questions with several confirmation questions (randomly generated)
 *        Use the questions generated by the algorithms to build the convex hull
 * @param pset
 * @param u
 * @param k
 * @param s
 * @param CH
 * @param HS
 * @param error
 * @return
 */
int unify(point_set* pset, point_t *u, int k, int s, int CH, int HS, int error)
{
    timeval t1;
    gettimeofday(&t1, 0);

    int numOfQuestion = 0, normalQuestion = 0;
    int dim = pset->points[0]->dim, M = pset->points.size();
    Partition *R = new Partition(dim), *RROrigin = new Partition(dim);
    std::vector<point_t*> considered;
    hyperplane_set *candHy = new hyperplane_set();
    std::vector<hyperplane*> potentialHyper, usedHyper;
    point_set *top_current = NULL;

    while(1)
    {
        std::vector<point_set *> topSet;
        point_set *candPt = new point_set();
        std::vector<Partition*> partitionSet;
        std::vector<int> considered_halfset;
        std::vector<choose_item*> chooseitemSet;

        //candidate hyper-plane generation
        if (CH == 1) //skyband
        {
            for (int i = 0; i < pset->points.size() - 1; ++i)
            {
                for (int j = i + 1; j < pset->points.size(); ++j)
                {
                    hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                    candHy->hyperplanes.push_back(h);
                }
            }
        }
        else if (CH == 2) //Randomized Grouping
        {
            if(considered.size() == 0)//initialization
            {
                int cut = 8 * sqrt(k) * sqrt(dim);
                if (cut > M) cut = M;
                for (int i = 0; i < cut; ++i)
                    considered.push_back(pset->points[i]);
                for (int i = 0; i < cut - 1; ++i)
                {
                    for (int j = i + 1; j < cut; ++j)
                    {
                        hyperplane *h = new hyperplane(considered[i], considered[j]);
                        int relation = R->check_relation(h);
                        if (relation == 0)
                            candHy->hyperplanes.push_back(h);
                    }
                }
            }
            else
            {
                int cIndex = considered.size();
                for(int i = 0; i < considered.size(); ++i)
                {
                    hyperplane *h = new hyperplane(considered[i], pset->points[cIndex]);
                    int relation = R->check_relation(h);
                    if (relation == 0)
                        candHy->hyperplanes.push_back(h);
                }
                considered.push_back(pset->points[cIndex]);
            }
        }
        else if (CH == 3) //Iterative Sampling
        {
            int dim = pset->points[0]->dim;
            double *min = new double[dim], *max = new double[dim];
            R->findMinMax(min, max);
            for(int i = 0; i < pset->points.size(); ++i)
                pset->points[i]->topk = 0;
            R->findTopk_extreme(pset, topSet, k);
            R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
            for (int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                    candPt->points.push_back(pset->points[i]);
            }
            construct_partitions(candPt, R, topSet, candHy->hyperplanes, partitionSet, considered_halfset);
            /*
            for(int i = 0; i < pset->points.size(); ++i)
            {
                if (pset->points[i]->topk == 1)
                {
                    for (int j = i + 1; j < pset->points.size(); ++j)
                    {
                        if (pset->points[j]->topk == 1)
                            candHy->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[j]));
                    }
                }
            }
            */
            if(partitionSet.size() <= 1)
            {
                top_current = partitionSet[0]->topSet;
                if (top_current != NULL)
                {
                    char algName[100];
                    sprintf(algName, "Unify%d%d", CH, HS);
                    //correctness checking
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
                        testsum += top_current->points[i]->dot_product(u);
                    }
                    double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                    top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                    return numOfQuestion;
                }
            }

        }
        std::cout << "CH finished. \n";



        for(int i = 0; i < candHy->hyperplanes.size(); ++i)
            usedHyper.push_back(candHy->hyperplanes[i]);

        int candSize = candHy->hyperplanes.size();
        int hIndex = -200;
        bool direction = true;
        while (candSize > 0)
        {
            //hyper-plane selection
            if (HS == 1)//Covering score
            {
                double maxScore = -1;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    //positive
                    int postive_num = 0;
                    Partition *testR_positive = new Partition(R);
                    testR_positive->hyperplanes.push_back(
                            new hyperplane(candHy->hyperplanes[i]->p_1, candHy->hyperplanes[i]->p_2));
                    testR_positive->set_ext_pts();
                    if(testR_positive->ext_pts.size() <= 0)
                        continue;
                    testR_positive->center = testR_positive->average_point();
                    for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                    {
                        if (i != j)
                        {
                            hyperplane *h = candHy->hyperplanes[j];
                            int relation = testR_positive->check_relation(h);
                            if (relation != 0)
                            {
                                postive_num++;
                            }
                        }
                    }

                    //negative
                    int negative_num = 0;
                    if (postive_num > maxScore)
                    {
                        Partition *testR_negative = new Partition(R);
                        testR_negative->hyperplanes.push_back(
                                new hyperplane(candHy->hyperplanes[i]->p_2, candHy->hyperplanes[i]->p_1));
                        testR_negative->set_ext_pts();
                        if(testR_negative->ext_pts.size() <= 0)
                            continue;
                        testR_negative->center = testR_negative->average_point();
                        for (int j = 0; j < candHy->hyperplanes.size(); ++j)
                        {
                            if (i != j)
                            {
                                hyperplane *h = candHy->hyperplanes[j];
                                int relation = testR_negative->check_relation(h);
                                if (relation != 0)
                                {
                                    negative_num++;
                                }
                            }
                        }
                    }

                    if (std::min(postive_num, negative_num) > maxScore)
                    {
                        maxScore = std::min(postive_num, negative_num);
                        hIndex = i;
                    }
                }
            }
            else if (HS == 2)
            {
                if(partitionSet.size() <=0)
                {
                    int dim = pset->points[0]->dim;
                    double *min = new double[dim], *max = new double[dim];
                    R->findMinMax(min, max);
                    for(int i = 0; i < pset->points.size(); ++i)
                        pset->points[i]->topk = 0;
                    R->findTopk_extreme(pset, topSet, k);
                    R->findTopk_sampling(pset, topSet, min, max,new point_t(dim), k, 0, 0);//use sampling method
                    for (int i = 0; i < pset->points.size(); ++i)
                    {
                        if (pset->points[i]->topk == 1)
                            candPt->points.push_back(pset->points[i]);
                    }
                    hyperplane_set *candHy_test = new hyperplane_set();
                    construct_partitions(candPt, R, topSet, candHy_test->hyperplanes, partitionSet, considered_halfset);
                }
                if(hIndex == -200)
                    hIndex = build_choose_item_table(partitionSet, candHy->hyperplanes, chooseitemSet);
                else
                    hIndex = modify_choose_item_table(chooseitemSet, partitionSet, considered_halfset,
                                                      hIndex, direction);
                if(hIndex == -1)
                    break;
            }
            else if (HS == 3)
            {
                double minDis = 100000000;
                for (int i = 0; i < candHy->hyperplanes.size(); ++i)
                {
                    hyperplane *h = candHy->hyperplanes[i];
                    double dst_h = h->check_distance(R->center);
                    if (dst_h < minDis)
                    {
                        minDis = dst_h;
                        hIndex = i;
                    }
                }
            }


            //normal interaction
            if(HS != 2 || chooseitemSet.size() >= 1 && considered_halfset.size() >= 2)
            {
                numOfQuestion++;
                normalQuestion++;
                point_t *p1, *p2;
                if(HS != 2)
                {
                    p1 = candHy->hyperplanes[hIndex]->p_1;
                    p2 = candHy->hyperplanes[hIndex]->p_2;
                }
                else
                {
                    p1 = chooseitemSet[hIndex]->h->p_1;
                    p2 = chooseitemSet[hIndex]->h->p_2;
                }
                double v1 = p1->dot_product(u);
                double v2 = p2->dot_product(u);
                int compareResult = u->compare(p1, p2, error);
                if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                    std::cout << "correct" << std::endl;
                else
                    std::cout << "wrong" << std::endl;
                if (compareResult > 0)
                {
                    hyperplane *half = new hyperplane(p2, p1);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = true;
                }
                else
                {
                    hyperplane *half = new hyperplane(p1, p2);
                    R->hyperplanes.push_back(half);
                    potentialHyper.push_back(half);
                    direction = false;
                }
                if(HS != 2)
                    candHy->hyperplanes.erase(candHy->hyperplanes.begin() + hIndex);
                R->set_ext_pts_withoutHull();
                R->center = R->average_point();

                candHy->hyperplanes.erase(
                        std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                       [R](hyperplane *h)
                                       { return R->check_relation(h) != 0; }),
                        candHy->hyperplanes.end());
            }

            if(HS == 2)
            {
                candSize = std::min(chooseitemSet.size(), considered_halfset.size() - 1);
                std::cout << "cand partition: " << candSize << "   chooseitem: " << chooseitemSet.size()
                          << std::endl;
            }
            else
            {
                candSize = candHy->hyperplanes.size();
                std::cout << "cand hyperplane: " << candSize << std::endl;
            }


            //stopping condition
            top_current = R->findTopk(pset, k, s);

            //R->print();
            //test question
            int Y = 3, testNum = 3;
            if ((normalQuestion) % Y == 0 || top_current != NULL)
            {
                for(int t = 0; t < testNum; ++t)
                {
                    double minDis = 10000;
                    hyperplane *testh = NULL;
                    int countTest = 0;
                    while(countTest < 10000)
                    {
                        int i = rand() % pset->points.size();
                        int j = rand() % pset->points.size();
                        hyperplane *h = new hyperplane(pset->points[i], pset->points[j]);
                        double dst_h = h->check_distance(R->center);
                        if (R->check_relation(h) != 0 && dst_h < minDis )
                        {
                            minDis = dst_h;
                            testh = h;
                        }
                        countTest++;
                    }

                    numOfQuestion++;
                    point_t *p1 = testh->p_1;
                    point_t *p2 = testh->p_2;
                    double v1 = p1->dot_product(u);
                    double v2 = p2->dot_product(u);
                    int compareResult = u->compare(p1, p2, error);
                    if(v1 > v2 && compareResult > 0 || v1 < v2 && compareResult < 0)
                        std::cout << "v1: " << v1 << "  v2: " << v2 << "  correct" << std::endl;
                    else
                        std::cout << "v1: " << v1 << "  v2: " << v2 << "  wrong" << std::endl;
                    if (compareResult > 0)
                    {
                        hyperplane *half = new hyperplane(p2, p1);
                        R->hyperplanes.push_back(half);
                        potentialHyper.push_back(half);
                    }
                    else
                    {
                        hyperplane *half = new hyperplane(p1, p2);
                        R->hyperplanes.push_back(half);
                        potentialHyper.push_back(half);
                    }
                }


                bool intersectionFound = R->set_ext_pts_withoutHull();
                if(intersectionFound && R->ext_pts.size() > 0)
                {
                    R->center = R->average_point();
                    RROrigin = new Partition(R);
                    RROrigin->center = RROrigin->average_point();
                    RROrigin->print();
                    std::cout << "IN: " << R->isIn(u) <<std::endl;

                }
                else
                {
                    for(int i = 0; i < potentialHyper.size(); ++i)
                    {
                        RROrigin->hyperplanes.push_back(new hyperplane(potentialHyper[i]));
                    }
                    double epsilon = 0.95;
                    while (!intersectionFound || RROrigin->ext_pts.size() <= 0)
                    {
                        epsilon += (1 - epsilon) / 2;
                        for (int i = RROrigin->hyperplanes.size() - potentialHyper.size(); i < RROrigin->hyperplanes.size(); ++i)
                        {
                            /*
                            point_t *qminusp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qminusp->attr[j] = psmall->attr[j] - plarge->attr[j];
                            point_t *qaddp = new point_t(dim);
                            for(int j = 0; j < dim; ++j)
                                qaddp->attr[j] = psmall->attr[j] + plarge->attr[j];
                            double off1 = testR->max_uq(RROrigin->hyperplanes[i], psmall);
                            off1 = qminusp->dot_product(qminusp) * (2 * epsilon - 1) * off1;
                            double off2 = qminusp->dot_product(psmall) - epsilon * qminusp->dot_product(qaddp);
                            RROrigin->hyperplanes[i]->offset = off1 / off2;
                             */
                            RROrigin->hyperplanes[i]->offset =  - log(epsilon / (1.0 - epsilon)); //- 1.0 / QuestionHyper[i]->p_1->distance(QuestionHyper[i]->p_2) *log(epsilon / (1.0 - epsilon));
                        }
                        intersectionFound = RROrigin->set_ext_pts_withoutHull();
                    }
                    R = new Partition(RROrigin);
                    R->center = R->average_point();
                    R->print();
                    std::cout << "IN: " << R->isIn(u) <<std::endl;

                    /*
                    //re-ask the questions
                    int *answerCollection = new int[potentialHyper.size()];
                    for(int i = 0; i < potentialHyper.size(); ++i)
                        answerCollection[i] = -1;
                    Partition *testR;

                    int counttest = -1;
                    do{
                        counttest++;
                        if(counttest > 10)
                            testR = new Partition(dim);
                        else
                            testR = new Partition(RROrigin);
                        for(int i = 0; i < potentialHyper.size(); ++i)
                        {
                            numOfQuestion += 10;
                            point_t *p1 = potentialHyper[i]->p_1;
                            point_t *p2 = potentialHyper[i]->p_2;
                            double v1 = p1->dot_product(u);
                            double v2 = p2->dot_product(u);
                            int compareResult = u->compare(p1, p2, error);
                            if (compareResult > 0)
                                answerCollection[i]++;
                            else
                                answerCollection[i]--;
                            while(answerCollection[i] == 0)
                            {
                                numOfQuestion++;
                                p1 = potentialHyper[i]->p_1;
                                p2 = potentialHyper[i]->p_2;
                                v1 = p1->dot_product(u);
                                v2 = p2->dot_product(u);
                                compareResult = 0;
                                for(int i = 0; i < 10; ++i)
                                {
                                    compareResult = u->compare(p1, p2, error);
                                    if (compareResult > 0)
                                        answerCollection[i]++;
                                    else
                                        answerCollection[i]--;
                                }
                            }
                        }
                        for(int i = 0; i < potentialHyper.size(); ++i)
                        {
                            if(answerCollection[i] == -1)
                            {
                                hyperplane *half = new hyperplane(potentialHyper[i]);
                                testR->hyperplanes.push_back(half);
                            }
                            else
                            {
                                hyperplane *half = new hyperplane(potentialHyper[i]);
                                for(int j = 0; j < dim; ++j)
                                    half->norm[j] = -half->norm[j];
                                testR->hyperplanes.push_back(half);
                            }
                        }
                        intersectionFound = testR->set_ext_pts_withoutHull();
                    }while(!intersectionFound);
                    R = new Partition(testR); RROrigin = new Partition(testR);
                    R->center = R->average_point(); RROrigin->center = RROrigin->average_point();
                    R->print();
                    */
                    candHy->hyperplanes.clear();
                    for (int i = 0; i < usedHyper.size(); ++i)
                        candHy->hyperplanes.push_back(usedHyper[i]);
                    candHy->hyperplanes.erase(
                            std::remove_if(candHy->hyperplanes.begin(), candHy->hyperplanes.end(),
                                           [R](hyperplane *h)
                                           { return R->check_relation(h) != 0; }),
                            candHy->hyperplanes.end());

                    topSet.clear();
                    candPt = new point_set();
                    partitionSet.clear();
                    considered_halfset.clear();
                    chooseitemSet.clear();
                    hIndex = -200;
                }
                potentialHyper.clear();
            }

            top_current = R->findTopk(pset, k, s);
            if (top_current != NULL)
            {
                char algName[100];
                sprintf(algName, "Unify%d%d", CH, HS);
                //correctness checking
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
                    testsum += top_current->points[i]->dot_product(u);
                }
                double accuracy = testsum / groudtruthsum > 1? 1: testsum / groudtruthsum;

                top_current->printResult(algName, numOfQuestion, s, t1, 0, accuracy);

                return numOfQuestion;
            }
        }

    }


}








