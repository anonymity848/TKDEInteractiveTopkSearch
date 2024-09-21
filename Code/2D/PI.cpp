#include "PI.h"


struct MinHeapCmp
{
    inline bool operator()(const u_vector* y, const u_vector* z)const
    {
        double v1=y->x;
        double v2=z->x;
        //printf("v1%lf v2%lf\n", v1, v2);
        if(v1 - v2<0.000000001&&v1 - v2>-0.000000001)
        {
            if(y->place == z->place)
            {
                return y->time < z->time;
            }
            else
            {
                return y->place > z->place;
            }

        }
        else if(v1 > v2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

/*
 * @brief Used for 2D dimensional datasets.
 * Algorithm Partition-based Binary Search
 * @param original_set       The original dataset
 * @param u                  User's utility function
 * @param k                  The threshold
 */
void PartitionBS(point_set *pset, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);

    //the sorted version
    std::vector<point_t *> Q;
    pset->sort_point(Q);
    int M = Q.size();
    //for(int i = 0; i < Q.size(); ++i)
    //    std::cout << Q[i]->attr[1] << "  " << Q[i]->dot_product(u0) << "\n";


    //initial the utility vectors for points order change
    std::vector<u_vector *> H;
    for (int i = 0; i < M - 1; i++)
    {
        //note that the point might be overlap
        double x0 = (Q[i]->attr[1] - Q[i + 1]->attr[1]) /
                    (Q[i + 1]->attr[0] - Q[i + 1]->attr[1] - Q[i]->attr[0] + Q[i]->attr[1]);

        //note that x0 in [0,1]
        if (EQN2 <= x0 && x0 <= 1)
        {
            u_vector *ut = new u_vector(x0, Q[i], Q[i + 1]);
            ut->inserted(H);
            //if (H[left]->x != x0 || H[left]->point_up->id != Q[i]->id || H[left]->point_down->id != Q[i + 1]->id)
            //H.insert(H.begin() + left, ut);
        }
    }

    std::vector<point_set *> representSet;
    int index = 1, count = 0;
    for (int i = 0; i < k; ++i)
        Q[i]->topk = index;
    std::vector<u_vector *> border;
    int change_num = 0;
    while (H.size() > 0) //&& H[0]->x < 1)
    {
        change_num++;
        /*
        printf("utility %d\n", change_num);
        {
            //for(int ttt=0; ttt < H.size(); ++ttt)
            int ttt = H.size() - 1;
            {
                printf("u %lf point1 %d %d point2 %d %d\n", H[ttt]->x, H[ttt]->point_up->id, H[ttt]->point_up->place,
                       H[ttt]->point_down->id, H[ttt]->point_down->place);
            }


            for(int ttt=0; ttt < H.size(); ++ttt)
            {
                if(H[ttt]->point_down->place == 10)
                    printf("uttt %lf point1 %d %d point2 %d %d\n", H[ttt]->x, H[ttt]->point_up->id, H[ttt]->point_up->place,
                       H[ttt]->point_down->id, H[ttt]->point_down->place);
            }

        }
        */

        //pop the utility vector
        while (H[H.size() - 1]->point_up->place == -1 || H[H.size() - 1]->point_down->place == -1)
        {
            H.pop_back();
            if(H.size() < 1)
                break;
        }
        if(H.size() < 1)
            break;
        u_vector *current_u = H[H.size() - 1];
        H.pop_back();


        //change the location of point_up, point_down
        point_t *p = current_u->point_down;
        std::vector<point_t*>::iterator iter = std::find(Q.begin(), Q.end(), p);
        int downPlace = std::distance(Q.begin(), iter);


        if (downPlace == k)//end the current partition and start a new one
        {
            //if the point is not top-k any more, delete it
            if (current_u->point_up->topk == index)
                count++;
            //if there is no more point in the current_top_k
            if (count >= k - s + 1)
            {
                representSet.push_back(new point_set());
                for (int j = 0; j < k; ++j)
                {
                    if (Q[j]->topk == index)
                        representSet[representSet.size() - 1]->points.push_back(Q[j]);
                }
                //representSet[representSet.size() - 1]->points.push_back(Q[k - 1]);
                border.push_back(current_u);
                count = 0;
                index++;

                for (int j = 0; j < k + 1; ++j)
                    Q[j]->topk = index;
            }
            current_u->point_up->topk = -1;
        }


        if(Q[downPlace - 1]->countPassed(current_u->x) >= k)
        {
            Q[downPlace - 1]->place = -1;
            Q.erase(Q.begin() + downPlace - 1);
            //insert new utility vector of point: 1. place_down - 2 and 2. place_down - 1
            if (downPlace > 1)
            {
                if (!Q[downPlace - 2]->is_changed(Q[downPlace - 1]))
                {
                    double x0 = (Q[downPlace - 2]->attr[1] - Q[downPlace - 1]->attr[1]) /
                                (Q[downPlace - 1]->attr[0] - Q[downPlace - 1]->attr[1] -
                                 Q[downPlace - 2]->attr[0] + Q[downPlace - 2]->attr[1]);

                    if (current_u->x < x0 && x0 <= 1)
                    {
                        u_vector *ut = new u_vector(x0, Q[downPlace - 2], Q[downPlace - 1]);
                        ut->inserted(H);
                    }
                }
            }
        }
        else
        {
            Q.erase(Q.begin() + downPlace);
            Q.insert(Q.begin() + downPlace - 1, p);

            //insert new utility vector of point: 1. place_down and 2. place_down + 1
            if (downPlace + 1 < Q.size())
            {
                if(!Q[downPlace + 1]->is_changed(Q[downPlace]))
                {
                    double x0 = (Q[downPlace]->attr[1] - Q[downPlace + 1]->attr[1]) /
                                (Q[downPlace + 1]->attr[0] - Q[downPlace + 1]->attr[1] -
                                 Q[downPlace]->attr[0] + Q[downPlace]->attr[1]);

                    if (current_u->x <= x0 && x0 <= 1)
                    {
                        u_vector *ut = new u_vector(x0, Q[downPlace], Q[downPlace + 1]);
                        ut->inserted(H);
                        //H.insert(H.begin() + left, ut);
                    }
                }
            }
            //insert new utility vector of point: 1. place_down - 2 and 2. place_down - 1
            if (downPlace > 1)
            {
                if (!Q[downPlace - 2]->is_changed(Q[downPlace - 1]))
                {
                    double x0 = (Q[downPlace - 2]->attr[1] - Q[downPlace - 1]->attr[1]) /
                                (Q[downPlace - 1]->attr[0] - Q[downPlace - 1]->attr[1] -
                                 Q[downPlace - 2]->attr[0] + Q[downPlace - 2]->attr[1]);

                    if (current_u->x < x0 && x0 <= 1)
                    {
                        u_vector *ut = new u_vector(x0, Q[downPlace - 2], Q[downPlace - 1]);
                        ut->inserted(H);
                    }
                }
            }
        }
    }

    representSet.push_back(new point_set());
    for (int i = 0; i < k && representSet[representSet.size() - 1]->points.size() < s; ++i)
    {
        if (Q[i]->topk == index)
        {
            representSet[representSet.size() - 1]->points.push_back(Q[i]);
        }
    }
    //representSet[representSet.size() - 1]->points.push_back(Q[k]);

    //now we have represent_point, border. Ask questions
    /*
    for(int i=0; i<border.size();i++)
    {
        printf("border %d u %lf\n", i, border[i]->x);
    }
    for(int i=0; i<represent_point.size();i++)
    {
        printf("point %d %lf %lf\n", represent_point[i]->id, represent_point[i]->coord[0], represent_point[i]->coord[1]);
    }
    */


    int numOfQuestion = 0;
    int b_left = 0, b_right = border.size() - 1;
    while (b_left <= b_right)
    {
        numOfQuestion++;
        int index = (b_left + b_right) / 2;
        double v1 = border[index]->point_up->dot_product(u);
        double v2 = border[index]->point_down->dot_product(u);
        if (v1 > v2)
        {
            b_right = index - 1;
        } else
        {
            b_left = index + 1;
        }
    }

    representSet[b_left]->printResult("PartitionBS", numOfQuestion, s, t1, 0, 0);

    point_set *resultSet = new point_set();
    pset->findTopk(u, k, resultSet);
    bool resultcheck = true;
    for(int i = 0; i < representSet[b_left]->points.size(); ++i)
    {
        if(!resultSet->checkExist(representSet[b_left]->points[i]))
        {
            resultcheck = false;
            break;
        }
    }
    std::cout << resultcheck <<"\n";


}


/*
 * @brief Used for 2D dimensional datasets.
 * Algorithm Segment-based Binary Search
 * @param original_set       The original dataset
 * @param u                  User's utility function
 * @param k                  The threshold
 */
void SegmentBS(point_set *pset, point_t *u, int k, int s)
{
    timeval t1;
    gettimeofday(&t1, 0);
    int numOfQuestion = 0;

    point_t* middle_u = new point_t(2);
    double lboundary = 0, rboundary = 1;

    while(1)
    {
        middle_u->attr[0] = (lboundary + rboundary) / 2;
        middle_u->attr[1] = 1 - middle_u->attr[0];
        middle_u->print();
        //find cell
        std::vector<point_t *> Q;
        pset->sort_point(Q, middle_u);
        int M = Q.size();


        double lr = lboundary, er = rboundary;
        int lr_index_l, lr_index_r, er_index_l, er_index_r;
        for (int i = 0; i < k; ++i)
        {
            for (int j = k; j < M; ++j)
            {
                double x0 = (Q[i]->attr[1] - Q[j]->attr[1]) /
                            (Q[j]->attr[0] - Q[j]->attr[1] - Q[i]->attr[0] + Q[i]->attr[1]);
                if (x0 > lr && x0 < middle_u->attr[0])
                {
                    lr = x0;
                    lr_index_l = j;
                    lr_index_r = i;
                } else if (x0 < er && x0 > middle_u->attr[0])
                {
                    er = x0;
                    er_index_l = i;
                    er_index_r = j;
                }
            }
        }


        //interaction
        if (lr != lboundary)
        {
            numOfQuestion++;
            double v1 = Q[lr_index_l]->dot_product(u);
            double v2 = Q[lr_index_r]->dot_product(u);
            //std::cout << v1 << "  " << v2 << "\n";
            if (v1 > v2)
            {
                rboundary = lr;
            } else
            {
                lboundary = lr;
            }
        }
        if (rboundary != lr && er != rboundary)
        {
            numOfQuestion++;
            double v1 = Q[er_index_l]->dot_product(u);
            double v2 = Q[er_index_r]->dot_product(u);
            //std::cout << v1 << "  " << v2 << "\n";
            if (v1 > v2)
            {
                rboundary = er;
            } else
            {
                lboundary = er;
            }
        }

        //check if there are s top-k points
        point_set *representSet = new point_set();
        if(lr == lboundary && er == rboundary)
        {
            for (int i = 0; i < s; ++i)
            {
                representSet->points.push_back(Q[i]);
            }
        }
        else
        {
            for (int i = 0; i < k; ++i)
            {
                int numLarger = 0;
                for (int j = i + 1; j < Q.size(); ++j)
                {
                    double x0 = (Q[i]->attr[1] - Q[j]->attr[1]) /
                                (Q[j]->attr[0] - Q[j]->attr[1] - Q[i]->attr[0] + Q[i]->attr[1]);
                    if (x0 > rboundary || x0 < lboundary)
                        numLarger++;
                }
                if (numLarger >= M - k)
                    representSet->points.push_back(Q[i]);
            }
        }
        if (representSet->points.size() >= s)
        {
            representSet->printResult("SegmentBS", numOfQuestion, s, t1, 0, 0);

            //check if the output is correct
            point_set *resultSet = new point_set();
            pset->findTopk(u, k, resultSet);
            bool resultcheck = true;
            for (int i = 0; i < representSet->points.size(); ++i)
            {
                if (!resultSet->checkExist(representSet->points[i]))
                {
                    resultcheck = false;
                    break;
                }
            }
            std::cout << "Accuracy: " << resultcheck << "\n";

            return;
        }
    }

}


