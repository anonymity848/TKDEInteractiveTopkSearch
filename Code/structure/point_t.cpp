#include <operation.h>
#include "point_t.h"

/**
 * @brief Constructor
 * @param dim   The number of dimensions
 */
point_t::point_t(int d)
{
    this->dim = d;
    attr = new double[dim];
    this->id = -1;
    this->index = -1;
    this->value = 0;
}

/**
 * @brief Constructor
 * @param d_order   The number of ordered dimensions
 * @param d_unorder The number of unordered dimensions
 * @param id        ID
 */
point_t::point_t(int d, int id)
{
    this->dim = d;
    attr = new double[dim];
    this->id = id;
    this->index = -1;
    this->value = 0;
}

/**
 * @brief Constructor
 * @param p   The point
 */
point_t::point_t(point_t *p)
{
    dim = p->dim;
    id = p->id;
    attr = new double[dim];
    for(int i = 0;i < dim; ++i)
        attr[i] = p->attr[i];
    index = p->index;
    value = p->value;
}

/*
 * @brief Constructor
 * @param p1   point
 */
point_t::point_t(point_t *p1, point_t *p2)
{
    dim = p1->dim;
    attr = new double[dim];
    for(int i = 0; i < dim; ++i)
        attr[i] = p1->attr[i] - p2->attr[i];
    id = -1;
    index = -1;
    this->value = 0;
}

/*
 * @brief Destructor
 *        delete the memory of the array
 */
point_t::~point_t()
{
    delete []attr;
}

/*
 * @brief For debug purpose, print the coordinates for a given point
 */
void point_t::print()
{
    //std::cout<< id <<"  ";
    for (int i = 0; i < dim; i++)
        std::cout<< attr[i] << "  ";
    std::cout << "place:" << place <<"  ";
    std::cout << "isTopk:" << topk;
    std::cout << "\n";
}

/**
 * @brief       Check whether the point is the same as p
 * @param p     Point
 * @return      1 same
 *              0 different
 */
bool point_t::is_same(point_t *p)
{
    if(dim != p->dim)
        return false;
    for (int i = 0; i < dim; ++i)
    {
        if (attr[i] - p->attr[i] < -EQN2 || attr[i] - p->attr[i] > EQN2)
            return false;
    }
    return true;
}

/**
 * @brief   Check whether all the attribute values are 0
 * @return  -true   all attributes values are 0
 *          -false  there exists attribute value which is not 0
 */
bool point_t::is_zero()
{
    for(int i = 0; i < dim; ++i)
    {
        if(attr[i] < -EQN2 || attr[i] > EQN2)
            return false;
    }
    return true;
}

/**
 * @brief	    Calculate the dot product between two points
 * @param p     One point
 */
double point_t::dot_product(point_t *p)
{
    double result = 0;
    for(int i = 0; i < dim; i++)
    {
        result += attr[i] * p->attr[i];
    }
    return result;
}

/**
 * @brief	    Calculate the dot product between two points
 * @param v     One array
 */
double point_t::dot_product(double *v)
{
    double result = 0;
    for(int i = 0; i < dim; i++)
    {
        result += attr[i] * v[i];
    }
    return result;
}

/**
 * @brief	Calculate the subtraction between two points.
 *          Remember to release the returned point to save memory.
 * @param p The subtractor
 * @return  The subtraction(new points)
 */
point_t *point_t::sub(point_t *p)
{
    point_t* result = new point_t(dim);
    for(int i = 0; i < dim; i++)
    {
        result->attr[i] = attr[i] - p->attr[i];
    }
    return result;
}

/**
 * @brief	Calculate the addition between two points.
 *          Remember to release the returned point to save memory.
 * @param p The point
 * @return  The addition(new points)
 */
point_t *point_t::add(point_t *p)
{
    point_t* result = new point_t(dim);
    for(int i = 0; i < dim; i++)
    {
        result->attr[i] = attr[i] + p->attr[i];
    }
    return result;
}

/**
 * @brief	Scale the point
 *          Remember to release the returned point to save memory.
 * @param c The scaled coefficient
 * @return  The scaled point
 */
point_t *point_t::scale(double c)
{
    point_t* result = new point_t(dim);
    for(int i = 0; i < dim; i++)
    {
        result->attr[i] = attr[i] * c;
    }
    return result;
}

/**
 * @brief Calculate the length of a vector
 * @return The length
 */
double point_t::cal_len()
{
    double diff = 0;
    for(int i = 0; i < dim; i++)
    {
        diff += (double) pow(attr[i], 2);
    }
    return sqrt(diff);
}


/**
 * @brief       Calculate the distance between two points
 * @param p     The points
 * @return      The distance
 */
double point_t::distance(point_t *p)
{
    double diff = 0;
    for(int i = 0; i < dim; i++)
    {
        diff += (double) pow(attr[i] - p->attr[i], 2);
    }
    return sqrt(diff);
}

/**
 * @brief       Calculate the distance between a points and a vector
 * @param p     The vector
 * @return      The distance
 */
double point_t::distance(double *p)
{
    double diff = 0;
    for(int i = 0; i < dim; i++)
    {
        diff += (double) pow(attr[i] - p[i], 2);
    }
    return sqrt(diff);
}

/**
 * @brief Print the result of the algorithm
 * @param out_cp    The name of the output file
 * @param name      The name of the algorithm
 * @param Qcount    The number of question asked
 * @param t1        The start time
 */
void point_t::printResult(char *name, int Qcount, timeval t1)
{
    timeval t2;
    std::ofstream out_cp("../../result.txt");
    gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << "-----------------------------------------------------------------\n";
    printf("|%15s |%15d |%15lf |%10d |\n", name, Qcount, time_cost, id);
    std::cout << "-----------------------------------------------------------------\n";
    out_cp << Qcount << "       " << time_cost << "\n";
    out_cp.close();
}


/**
 * @brief Print the result of the algorithm
 * @param out_cp    The name of the output file
 * @param name      The name of the algorithm
 * @param Qcount    The number of question asked
 * @param t1        The start time
 * @param preTime   The preprocessing time cost
 */
void point_t::printResult(char *name, int Qcount, timeval t1, double preTime, long mem_baseline)
{
    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << "-----------------------------------------------------------------------------------\n";
    //printf("|%15s |%15d |%15lf |%15lf |%15d |%10d |\n", name, Qcount, preTime, time_cost - preTime, get_mem_usage() - mem_baseline, id);
    printf("|%15s |%15d |%15lf |%15lf |%10d |\n", name, Qcount, preTime, time_cost - preTime, id);
    std::cout << "-----------------------------------------------------------------------------------\n";
    //out_cp << Qcount << "       " << preTime << "       " << time_cost - preTime << "      " << get_mem_usage() - mem_baseline << "\n";
}


/**
 * @brief calculate the place where the order of the two points change
 * @param p  The second point
 * @return
 */
double point_t::bound(point_t *p, double x)
{
    double sum = 0;
    sum =  (this->attr[0] * this->attr[0]) - (p->attr[0] * p->attr[0])
             + (this->attr[1] * this->attr[1]) - (p->attr[1] * p->attr[1]);
    sum = sum / 2 - x * (this->attr[0] - p->attr[0]);
    sum = sum /  (this->attr[1] - p->attr[1]);
    return sum;
}

/**
 * @brief Check whether the point dominates p
 * @param p     The point p
 * @return
 */
bool point_t::dominate(point_t *p)
{
    int dnum = 0;
    for (int i = 0; i < dim; i++)
    {
        if (attr[i] < p->attr[i])
            return false;
        if (attr[i] == p->attr[i])
        {
            dnum++;
        }
    }
    if(dnum >= dim)
    {
        return false;
    }
    return true;
}

/**
 * @brief Check whether the location of two points has already changed
 * @param p Point 2
 * @return  true it changes
 *          false it does not change
 */
bool point_t::is_changed(point_t *p)
{
    for (int i = 0; i < surpass.size(); ++i)
    {
        if (surpass[i]->id == p->id)
            return true;
    }

    for (int i = 0; i < p->surpass.size(); ++i)
    {
        if (p->surpass[i]->id == id)
            return true;
    }
    return false;
}

/**
 * @brief Count the number of points which has surpassed the point
 * @param ut
 * @return
 */
int point_t::countPassed(double ut)
{
    int count = 0;
    for(int i = 0; i < surpassVector.size(); ++i)
    {
        if(surpassVector[i] <= ut)
            count++;
    }
    return count;
}

/**
 * @brief Check whether each p[i] >= 0
 * @return  1 Each dimension is positive
 *          -1 Exist dimension which is negative
 */
bool point_t::is_positive()
{
    for(int i = 0; i < dim; ++i)
    {
        if(attr[i] < 0)
            return false;
    }
    return true;
}



/**
 * @brief Check the utilities between p1 and p2
 * @param p1
 * @param p2
 * @param error  Indicate if there is error involved
 * @return
 */
int point_t::compare(point_t *p1, point_t *p2, int error)
{
    double v1 = this->dot_product(p1);
    double v2 = this->dot_product(p2);
    if(!error)
    {
        if(v1 > v2)
            return 1;
        else
            return -1;
    }
    else
    {
        double prob;
        if (v1 > v2)
            prob = v1 - v2;
        else
            prob = v2 - v1;
        prob = 1.0 / (1.0 + exp(-1.0 * prob));
        double sumr = (double)rand() / RAND_MAX;
        if(v1 > v2 && sumr < prob || v1 < v2 && sumr > prob)
            return 1;
        else
            return -1;
    }
}














