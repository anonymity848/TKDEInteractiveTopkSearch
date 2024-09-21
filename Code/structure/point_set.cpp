#include "point_set.h"
#include "../qhull/qhull_build.h"


/**
 * @brief Constructor
 */
point_set::point_set(){}

/**
 * @brief Constructor
 *        Create a point set the same as p_set, all the points are re-created
 * @param p_set     Point set
 */
point_set::point_set(point_set *p_set)
{
    int num = p_set->points.size();
    point_t *p;
    for (int i = 0; i < num; i++)
    {
        p = new point_t(p_set->points[i]);
        points.push_back(p);
    }

}

/**
 * @brief Constructor
 *        Record all the points in the input file to the point set
 * @param input     Name of the data file.
 */
point_set::point_set(const char* input)
{
    FILE* c_fp;
    char filename[MAX_FILENAME_LENG];
    sprintf(filename, "../input/%s", input);
    printf("%s\n", filename);
    if ((c_fp = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "Cannot open the data file %s.\n", filename);
        exit(0);
    }

    int num, dim;
    point_t* p;
    fscanf(c_fp, "%i%i", &num, &dim);

    // read points line by line
    for (int i = 0; i < num; i++)
    {
        p = new point_t(dim, i);
        for (int j = 0; j < dim; j++)
            fscanf(c_fp, "%lf", &p->attr[j]);
        points.push_back(p);
    }
    fclose(c_fp);
}

/**
 *@brief  Destructor
 *        Delete the points in the array
 */
point_set::~point_set()
{
    int i = points.size();
    point_t *p;
    while(i>0)
    {
        p = points[i-1];
        points.pop_back();
        delete p;
        i--;
    }
    points.clear();
}

/**
 *@brief   Print all the points in the set
 */
void point_set::print()
{
    for (int i = 0; i < points.size(); i++)
        points[i]->print();
    printf("\n");
}

/**
 * @brief           Reload the points Randomly
 * @param RandRate  The parameter used to control how many points reinserted
 */
void point_set::random(double RandRate)
{
    int size = points.size();
    //reinsert
    for (int i = 0; i < size * RandRate; i++)
    {
        int n = ((int) rand()) % size;
        point_t *p = points[n];
        points.erase(points.begin() + n);
        points.push_back(p);
    }
}

/**
 * @brief       Sort points based on their utilites w.r.t. u
 * @param u     The utility vector
 * @return      The point set which contains all the point in order
 */
point_set* point_set::sort(point_t *u)
{
    int size = points.size();
    if(size <= 0)
    {
        cout << "Warning: The point set is empty.";
        return NULL;
    }

    point_set *return_set = new point_set();
    return_set->points.push_back(points[0]);
    for (int i = 1; i < size; i++)
    {
        double v0 = points[i]->dot_product(u);
        int left = 0, right = return_set->points.size() - 1;
        //find the place for p_set[i] in return_point and record the place index in "left"
        while (left <= right)
        {
            int middle = (left + right) / 2;
            double v = return_set->points[middle]->dot_product(u);
            if (v0 < v)
            {
                right = middle - 1;
            }
            else
            {
                left = middle + 1;
            }
        }
        return_set->points.insert(return_set->points.begin() + left, points[i]);
    }
    /*
    for(int i=0; i<return_set->points.size();i++)
    {
        return_set->points[i]->print();
    }
    */
    return return_set;
}

/**
 * @brief Find the top-k points w.r.t. u in the dataset
 * @param u      The utility vector
 * @param k      The number of returned points
 * @param topSet The returned set
 */
void point_set::findTopk(point_t *u, int k, point_set *topSet)
{

    std::vector<double> valueSet;
    //set the initial k points
    topSet->points.push_back(points[0]);
    valueSet.push_back(points[0]->dot_product(u));
    for (int i = 1; i < k; i++)
    {
        double value = points[i]->dot_product(u);
        int j = 0;
        while (j < valueSet.size() && value < valueSet[j])
            j++;
        topSet->points.insert(topSet->points.begin() + j, points[i]);
        valueSet.insert(valueSet.begin() + j, value);
    }

    //insert the rest points
    for (int i = k; i < points.size(); i++)
    {
        double value = points[i]->dot_product(u);
        for (int j = 0; j < k; ++j)
        {
            if (value > valueSet[j])
            {
                topSet->points.insert(topSet->points.begin() + j, points[i]);
                valueSet.insert(valueSet.begin() + j, value);
                topSet->points.pop_back();
                valueSet.pop_back();
                break;
            }
        }
    }
}

/**
 * @brief Find all the convex hull points in the dataset
 * @param topSet The returned set
 */
void point_set::find_convexHull(std::vector<point_set *> &topSet)
{
    int dim = points[0]->dim;
    int M = points.size(), original_size = points.size();
    FILE *rPtr = NULL, *wPtr = NULL;
    rPtr = (FILE *) fopen("../output/point.txt", "w");
    // write the feasible point
    int nn= M+1;
    fprintf(rPtr, "%i\n%i\n", dim, nn);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < dim; j++)
            fprintf(rPtr, "%lf ", points[i]->attr[j]);
        fprintf(rPtr, "\n");
    }
    for(int j = 0; j < dim; j++)
    {
        fprintf(rPtr, "%lf ", 0.0);
    }
    fprintf(rPtr, "\n");
    fclose(rPtr);
    rPtr = (FILE *) fopen("../output/point.txt", "r");
    wPtr = (FILE *) fopen("../output/top.txt", "w");
    Convex_H(rPtr, wPtr);
    fclose(rPtr);
    fclose(wPtr);

    //obtain the points
    wPtr = (FILE *) fopen("../output/top.txt", "r");
    fscanf(wPtr, "%i", &M);
    int index;
    for (int i = 0; i < M; i++)
    {
        fscanf(wPtr, "%i", &index);
        if(index < original_size)
        {
            point_set *p_set = new point_set();
            points[index]->topk = 1;
            p_set->points.push_back(points[index]);
            topSet.push_back(p_set);
        }
    }
    fclose(wPtr);
}


/**
 * @brief Find the top-k points w.r.t. u in the dataset
 * @param u      The utility vector
 * @param k      The number of returned points
 * @param topSet The returned set
 */
void point_set::findRanking(point_t *u, point_set *rankingSet)
{
    int M = points.size();
    std::vector<double> valueSet;
    //set the initial k points
    rankingSet->points.push_back(points[0]);
    valueSet.push_back(points[0]->dot_product(u));

    //insert points
    for (int i = 1; i < points.size(); i++)
    {
        double value = points[i]->dot_product(u);
        for (int j = 0; j < valueSet.size(); ++j)
        {
            if (value > valueSet[j])
            {
                rankingSet->points.insert(rankingSet->points.begin() + j, points[i]);
                valueSet.insert(valueSet.begin() + j, value);
                break;
            }
        }
        if(value <= valueSet[valueSet.size() - 1])
        {
            rankingSet->points.push_back(points[i]);
            valueSet.push_back(value);
        }
    }
}


/**
 * @brief   Write the dataset to the txt file
 * @param   fileName  The name of the txt file
 */
void point_set::write(std::string fileName)
{
    ofstream wPtr;
    wPtr.open(fileName, std::ios::out);
    wPtr.setf(ios::fixed, ios::floatfield);  // set as fixed model
    wPtr.precision(6);  // set precision to 6


    int size = points.size(), dim = points[0]->dim;
    // write the points
    wPtr << size << "   " << dim << " \n";//record the offset as one dimension
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            wPtr << points[i]->attr[j] <<" ";
        }
        wPtr <<"\n";
    }
    wPtr.close();
}

/**
 * @param Delete point p in the point set
 * @param p The point
 */
void point_set::prunePt(point_t *p)
{
    for(int i = 0; i < points.size(); ++i)
    {
        if(points[i]->id == p->id)
        {
            points.erase(points.begin() + i);
            return;
        }
    }
}

/**
 * @brief Check the points in the two sets are the same
 * @param pset  Point set 2
 * @return      1 They are the same
 *              -1 They are different
 */
bool point_set::isSame(point_set *pset)
{
    for(int i = 0; i < points.size(); ++i)
    {
        bool exist = false;
        for(int j = 0; j < pset->points.size(); ++j)
        {
            if(points[i]->id == pset->points[j]->id)
            {
                exist = true;
                break;
            }
        }
        if(!exist)
            return false;
    }
    return true;
}


/**
 * @brief Check the points in the two sets are the same
 * @param pset  Point set 2
 * @return      1 They are the same
 *              -1 They are different
 */
bool point_set::isSame_exact(point_set *pset)
{
    for(int i = 0; i < points.size(); ++i)
    {
        if(points[i]->id != pset->points[i]->id)
            return false;
    }
    return true;
}

//@brief Sort points based on their utility w.r.t. linear function u
// 		 Does not change the original dataset
//@param p_set 			The original dataset
//@param u 				The utility vector
void point_set::sort_point(std::vector<point_t*> &return_point, point_t *u)
{
    int M = points.size();
    return_point.push_back(points[0]);
    for (int i = 1; i < M; ++i)
    {
        double v0 = points[i]->dot_product(u);
        int left = 0, right = return_point.size() - 1;
        //find the place for p_set[i] in return_point
        //record the place index in "left"
        while (left <= right)
        {
            int middle = (left + right) / 2;
            double v = return_point[middle]->dot_product(u);
            if (v0 > v)
            {
                right = middle - 1;
            }
            else
            {
                left = middle + 1;
            }
        }
        points[i]->topk = -1;
        return_point.insert(return_point.begin() + left, points[i]);
    }
    for (int i = 1; i < M; ++i)
        return_point[i]->place = i;
    /*
    for(int i=0; i<return_point.size();i++)
    {
        printf("point %d  %lf %lf\n", return_point[i]->id, return_point[i]->coord[0], return_point[i]->coord[1]);
    }
    */
}



//@brief Sort points based on their p[0], and then p[1]
// 		 Do not change the original dataset
//@param return_point	The returned dataset containing all the sorted points
void point_set::sort_point(std::vector<point_t*> &return_point)
{
    int M = points.size();
    return_point.push_back(points[0]);
    for (int i = 1; i < M; ++i)
    {
        int left = 0, right = return_point.size() - 1;
        //find the place for p_set[i] in return_point
        //record the place index in "left"
        while (left <= right)
        {
            int middle = (left + right) / 2;
            if (points[i]->attr[1] > return_point[middle]->attr[1] ||
                    (points[i]->attr[1] == return_point[middle]->attr[1] &&
                    points[i]->attr[0] - points[i]->attr[1] > return_point[middle]->attr[0] - return_point[middle]->attr[1]))
            {
                right = middle - 1;
            }
            else
            {
                left = middle + 1;
            }
        }
        points[i]->topk = -1;
        return_point.insert(return_point.begin() + left, points[i]);
    }
    for (int i = 0; i < M; ++i)
        return_point[i]->place = i;
    /*
    for(int i=0; i<return_point.size();i++)
    {
        printf("point %d  %lf %lf\n", return_point[i]->id, return_point[i]->coord[0], return_point[i]->coord[1]);
    }
    */
}


/**
 * @brief Print the result of the algorithm
 * @param out_cp    The name of the output file
 * @param name      The name of the algorithm
 * @param Qcount    The number of question asked
 * @param t1        The start time
 * @param preTime   The preprocessing time cost
 */
void point_set::printResult(char *name, int Qcount, int s, timeval t1, double preTime, double acc)
{
    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << "-----------------------------------------------------------------------------------\n";
    //printf("|%15s |%15d |%15lf |%15lf |%15d |%10d |\n", name, Qcount, preTime, time_cost - preTime, get_mem_usage() - mem_baseline, id);
    printf("|%15s |%15d |%15lf |%15lf |%10s |\n", name, Qcount, time_cost - preTime, acc, "Points");
    for(int i = 0; i < s; ++i)
        printf("|%15s |%15s |%15s |%15s |%10d |\n", "-", "-", "-", "-", points[i]->id);
    std::cout << "-----------------------------------------------------------------------------------\n";
    std::ofstream out_cp("../result.txt");
    out_cp << Qcount << "       " << time_cost - preTime << "       " << acc << "\n"; //<< get_mem_usage() - mem_baseline << "\n";
}



/**
 * @brief Find all the points which are not dominated by >=k points
 * @param returnSet 	The returned points which are not dominated by >=k points
 * @param k 			The threshold
 */
void point_set::skyband(point_set *returnSet, int k)
{
    int num = points.size();
    point_t *pt;
    int *dominated = new int[num + 1];
    for (int i = 0; i < num; i++)
    {
        pt = points[i];
        dominated[pt->id] = 0;
        //check if pt is dominated k times by the return_point so far
        for (int j = 0; j < returnSet->points.size() && dominated[pt->id] < k; j++)
        {
            if(returnSet->points[j]->is_same(pt))
            {
                dominated[pt->id] = k;
                break;
            }
            else if (returnSet->points[j]->dominate(pt))
                dominated[pt->id]++;
        }
        if (dominated[pt->id] < k)
        {
            //eliminate any points in return_point dominated k times
            int m = returnSet->points.size();
            int index = 0;
            for (int j = 0; j < m; j++)
            {
                if (pt->dominate(returnSet->points[index]))
                {
                    dominated[returnSet->points[index]->id]++;
                }
                if (dominated[returnSet->points[index]->id] >= k)
                {
                    returnSet->points.erase(returnSet->points.begin() + index);
                }
                else
                {
                    index++;
                }
            }
            returnSet->points.push_back(pt);
        }
    }

    delete[] dominated;
    //printf("size%d\n", return_point.size());
    /*
    for(int i=0; i<return_point.size();i++)
    {
        printf("point %d %lf %lf\n", return_point[i]->id, return_point[i]->coord[0], return_point[i]->coord[1]);
    }
    */
}

/**
 * @brief   Check whether the point is in the set
 * @return  true the point exist
 *          false the point does not exist
 */
bool point_set::checkExist(point_t* p)
{
    for(int i = 0; i < points.size(); ++i)
    {
        if(p->id == points[i]->id)
            return true;
    }
    return false;
}

/**
 * @brief Find the same points in the two sets
 * @param pset  The second point set
 * @return  The same points
 */
point_set *point_set::findsame(point_set *pset)
{
    point_set *resultSet = new point_set();
    for(int i = 0; i < points.size(); ++i)
    {
        for(int j = 0; j < pset->points.size(); ++j)
        {
            if(points[i]->id == pset->points[j]->id)
            {
                resultSet->points.push_back(points[i]);
                break;
            }
        }
    }
    return resultSet;
}

/**
 * @brief Find the point which has the highest utility w.r.t. u
 * @param u The utility vector
 * @return  The index of the point in the set
 */
int point_set::findBest(point_t *u)
{
    int index = 0;
    double value = points[0]->dot_product(u);
    for(int i = 1; i < points.size(); ++i)
    {
        double v = points[i]->dot_product(u);
        if(v > value)
        {
            value = v;
            index = i;
        }
    }
    return index;
}



//@brief Use sampling to find all the points which is able to be top-1 at some utility vector
//@param p_set 		The point set containing all the points
//@param top_set	The returned point set containing all the possible top-1 point
//@param u 			The utility vector. For user, point_t* u = alloc_point(dim)
//@param level		The number of dimensions we have set. For user, only need to set level=0
//@param used_seg	The range which has been assigned to the u[i]. For user, set rest_seg=0
void point_set::findTopk_sampling(std::vector<point_set*> &topSet, point_t *u, int k, int level, int used_seg)
{
    int dim = points[0]->dim, M = points.size();
    double segment = 3, length = 1 / segment;
    if (level >= dim - 2)
    {
        for (int j = 0; j <= segment - used_seg; j++)
        {
            u->attr[level] = j * length;
            u->attr[dim - 1] = 1;
            for (int i = 0; i < dim - 1; i++)
                u->attr[dim - 1] -= u->attr[i];

            //u->print();
            //Find the top-k set w.r.t u
            point_set *topk = new point_set();
            findTopk(u, k, topk);

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
                for(int i = 0; i < topk->points.size(); ++i)
                    topk->points[i]->topk = 1;
            }
        }
    }
    else
    {
        for (int i = 0; i <= segment - used_seg; i++)
        {
            u->attr[level] = i * length;
            findTopk_sampling(topSet, u, k, level + 1, used_seg + i);
        }
    }
}





/**
 * @brief Use sampling to find all the possible ranking of points
 * @param rankSet   The returned ranking
 * @param u         The utility vector. For user, point_t* u = alloc_point(dim)
 * @param level     The number of dimensions we have set. For user, only need to set level=0
 * @param used_seg  The range which has been assigned to the u[i]. For user, set rest_seg=0
 */
void point_set::findRanking_sampling(std::vector<point_set*> &rankSet, point_t *u, int level, int used_seg)
{
    int dim = points[0]->dim, M = points.size();
    double segment = 10, length = 1 / segment;
    if (level >= dim - 2)
    {
        for (int j = 0; j <= segment - used_seg; j++)
        {
            u->attr[level] = j * length;
            u->attr[dim - 1] = 1;
            for (int i = 0; i < dim - 1; i++)
                u->attr[dim - 1] -= u->attr[i];

            //u->print();
            //Find the ranking w.r.t u
            point_set *rank = new point_set();
            findRanking(u, rank);

            //Check if it is already in top_set
            bool is_inside = false;
            for (int i = 0; i < rankSet.size(); i++)
            {
                if (rankSet[i]->isSame_exact(rank))
                {
                    is_inside = true;
                    break;
                }
            }
            if (!is_inside)
                rankSet.push_back(rank);
        }
    }
    else
    {
        for (int i = 0; i <= segment - used_seg; i++)
        {
            u->attr[level] = i * length;
            findRanking_sampling(rankSet, u, level + 1, used_seg + i);
        }
    }
}



//@brief Use sampling to find all the points which is able to be top-1 at some utility vector
//@param p_set 		The point set containing all the points
//@param top_set	The returned point set containing all the possible top-1 point
//@param u 			The utility vector. For user, point_t* u = alloc_point(dim)
//@param level		The number of dimensions we have set. For user, only need to set level=0
//@param used_seg	The range which has been assigned to the u[i]. For user, set rest_seg=0
void point_set::findTopk_sampling(std::vector<point_set*> &topSet, double *max, double *min, point_t *u, int k, int level, int used_seg)
{
    int dim = points[0]->dim, M = points.size();
    double segment = 5;
    if (level >= dim - 2)
    {
        for (int j = 0; j <= segment; j++)
        {
            u->attr[level] = min[level] + j * (max[level] - min[level]) / segment;
            u->attr[dim - 1] = 1;
            for (int i = 0; i < dim - 1; i++)
                u->attr[dim - 1] -= u->attr[i];

            u->print();
            if(u->is_positive()) //Find the top-k set w.r.t u
            {
                point_set *topk = new point_set();
                findTopk(u, k, topk);

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
            findTopk_sampling(topSet, max, min, u, k, level + 1, used_seg + i);
        }
    }
}


















