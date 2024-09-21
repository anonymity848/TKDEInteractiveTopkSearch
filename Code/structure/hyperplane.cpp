#include <iomanip>
#include "hyperplane.h"

/**
 * @brief Constructor
 */
hyperplane::hyperplane()
{
    p_1 = NULL;
    p_2 = NULL;
}

/**
 * @brief       Constructor
 * @param dim   The number of dimensions of the hyperplane
 */
hyperplane::hyperplane(int dim)
{
    this->dim = dim;
    norm = new double[dim];
    p_1 = NULL;
    p_2 = NULL;
}


/**
 * @brief Constructor
 * @param h The hyperplane
 */
hyperplane::hyperplane(hyperplane *h)
{
    dim = h->dim;
    norm = new double[dim];
    for (int j = 0; j < dim; j++)
        norm[j] = h->norm[j];
    offset = h->offset;
    p_1 = h->p_1;
    p_2 = h->p_2;
    is_seleted = h->is_seleted;
}

/**
 * @brief Constructor
 *        Guarantee that n has at least one element
 * @param n         The norm vector
 * @param offset    The offset
 */
hyperplane::hyperplane(int dim, double *n, double offset)
{
    this->dim = dim;
    norm = n;
    this->offset = offset;
    p_1 = NULL;
    p_2 = NULL;
    is_seleted = 0;
}

/**
 * @brief Constructor
 * @param p1        The first point consisting of the hyperplane
 * @param p2        The second point consisting of the hyperplane
 *
 */
hyperplane::hyperplane(point_t *p1, point_t *p2)
{
    dim = p1->dim;
    offset = 0;
    norm = new double[dim];
    for (int i = 0; i < dim; i++)
        norm[i] = p1->attr[i] - p2->attr[i];
    p_1 = p1;
    p_2 = p2;
    is_seleted = 0;
}


/**
 * @brief Constructor
 * @param d         The number of dimensions
 * @param p1        The first vector consisting of the hyperplane
 * @param p2        The second vector consisting of the hyperplane
 * @param offset    The offset
 * @param exp       The point which contains the value of the ordered dimension
 *
 */
hyperplane::hyperplane(int d, double *p1, double *p2)
{
    dim = d;
    offset = 0;
    norm = new double[dim];
    for (int i = 0; i < dim; i++)
        norm[i] = p1[i] - p2[i];
    p_1 = NULL;
    p_2 = NULL;
    is_seleted = 0;
}

/**
 * @brief Destructor
 *        Delete the array of norm
 */
hyperplane::~hyperplane()
{
    delete []norm;
}

/**
 * @brief  Print the information of the hyperplane
 */
void hyperplane::print()
{
    for (int i = 0; i < dim; i++)
    {
        std::cout << std::setprecision(10) << norm[i] << " ";
    }
    std::cout << offset << "\n";

}

/**
 * @brief Check the position of the point w.r.t. the hyperplane
 * @param p The point
 * @return   1 The point is in h+
 *          -1 The point is in h-
 */
int hyperplane::check_position(point_t *p)
{
    double sum = 0;
    for(int i = 0; i < dim; ++i)
        sum += p->attr[i] * norm[i];
    sum += offset;
    if(sum >= EQN2)
        return 1;
    else if (sum <= -EQN2)
        return -1;
    else
        return 0;
}


/**
 * @brief Check the position of the point w.r.t. the hyperplane
 * @param p The point
 * @return   1 The point is in h+
 *          -1 The point is in h-
 */
int hyperplane::check_positionlose(point_t *p)
{
    double sum = 0;
    for(int i = 0; i < dim; ++i)
        sum += p->attr[i] * norm[i];
    sum += offset;
    if(sum >= EQN_EPS)
        return 1;
    else if (sum <= -EQN_EPS)
        return -1;
    else
        return 0;
}


/**
 * @brief Check the position of the point w.r.t. the hyperplane
 * @param p The point
 * @return   1 The point is in h+
 *          -1 The point is in h-
 */
int hyperplane::check_positive(point_t *p)
{
    double sum = 0;
    for(int i = 0; i < dim; ++i)
        sum += p->attr[i] * norm[i];
    sum += offset;
    if(sum >= -EQN3)
        return 1;
    else
        return 0;
}

/**
 * @brief Calculate the distance from the point to the hyperplane
 * @param p     The point
 * @return      The distance
 */
double hyperplane::check_distance(point_t *p)
{
    double numerato = 0, dinominator = 0;
    for(int i = 0; i < dim; ++i)
    {
        numerato += p->attr[i] * norm[i];
        dinominator += norm[i] * norm[i];
    }
    numerato += offset;
    if(numerato < 0)
        numerato = -numerato;
    dinominator = sqrt(dinominator);
    return numerato/dinominator;

}

/**
 * @brief Calculate the priority of the hyperplane
 * @param pset The point set
 * @return The priority
 */
int hyperplane::priority(point_set *pset)
{
    if(p_1->value == 2 || p_2->value == 2)
        return -1;
    int d = pset->points[0]->dim;
    hyperplane *h = new hyperplane(d);
    h->norm = new double[d]; h->offset = 0;
    double weight = 0;
    for (int i = 0; i < d; i++)
    {
        h->norm[i] = p_1->attr[i] - p_2->attr[i];
        h->offset += (- p_1->attr[i] - p_2->attr[i]) * h->norm[i]/2;
        weight += h->norm[i];
    }

    weight = sqrt(weight);
    int up = 0, down = 0;
    for(int i = 0; i < pset->points.size(); ++i)
    {
        double sum = 0;
        for(int j = 0; j < dim; ++j)
            sum += pset->points[i]->attr[j] * norm[j];
        sum += offset;
        sum = abs(sum)/weight;
        //int pstion = h->check_position(pset->points[i]);
        if(sum >= 80)
            ++up;
        else if(sum <= -80)
            ++down;
        if(pset->points[i]->id == p_1->id)
            up += 3;
        else if (pset->points[i]->id == p_2->id)
            down += 3;
    }
    if(up < down)
        return up;
    else
        return down;
}


/**
 * @brief Calculate the priority of the hyperplane
 * @param pset The point set
 * @return The priority
 */
int hyperplane::priority2(point_set *pset)
{
    if(p_1->value == 2 || p_2->value == 2)
        return -1;
    int d = pset->points[0]->dim;

    int up = 0, down = 0;
    for(int i = 0; i < pset->points.size(); ++i)
    {
        int pstion = check_position(pset->points[i]);
        if(pstion == 1)
            ++up;
        else if(pstion == -1)
            ++down;
        if(pset->points[i]->id == p_1->id)
            up += 3;
        else if (pset->points[i]->id == p_2->id)
            down += 3;
    }
    if(up < down)
        return up;
    else
        return down;
}


/**
 * @brief Calculate the priority of the hyperplane
 * @param pset The point set
 * @return The priority
 */
int hyperplane::priority(point_set *pset, double dNN, double Beta)
{
    if(p_1->value == 2 || p_2->value == 2)
        return -1;
    int d = pset->points[0]->dim;

    int up = 0, down = 0;
    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(pset->points[i]->id == p_1->id)
            up += Beta;
        else if (pset->points[i]->id == p_2->id)
            down += Beta;
        else
        {
            double dis = check_distance(pset->points[i]);
            if (dis > dNN)
                up += Beta;
            else if (dis > 0)
                ++up;
            else if (dis < -dNN)
                down += Beta;
            else if (dis < 0)
                ++down;
        }
    }
    if(up < down)
        return up;
    else
        return down;
}

bool hyperplane::is_boundary(point_set *pset)
{
    int num = 0;
    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(p_1->id == pset->points[i]->id || p_2->id == pset->points[i]->id)
            ++num;
        if(num >=2 )
            return true;
    }
    return false;
}


/**
 * @brief Calculate the distance from a point to the hyperplane
 * @param p  The point
 * @return   The distance
 */
double hyperplane::distance(point_t *p)
{
    double sum = 0, value = 0;
    for (int i = 0; i < dim; i++)
    {
        value = value + norm[i] * p->attr[i];
        sum = sum + norm[i] * norm[i];
    }
    sum = sqrt(sum);
    if (value >= 0)
    {
        return value / sum;
    }
    else
    {
        return -value / sum;
    }
}



double hyperplane::bound(point_t *p1, point_t *p2)
{
    double x1 = p1->attr[0] - p2->attr[0];
    double y1 = p1->attr[1] - p2->attr[1];
    double off = (-p1->attr[0] - p2->attr[0]) * x1 / 2 + (-p1->attr[1] - p2->attr[1]) * y1 / 2;

    if(norm[0] != 0 && norm[1] != 0)
    {

        double x = (this->offset / this->norm[1] - off / y1) / (x1 / y1 - norm[0] / norm[1]);
        double y = (this->offset / this->norm[0] - off / x1) / (y1 / x1 - norm[1] / norm[0]);
        return x * norm[0] + y * norm[1];
    }
    else if (norm[0] == 0)
    {
        double y = -this->offset / norm[1];
        double x = -(off + y * y1) / x1;
        return x * norm[0] + y * norm[1];
    }
    else //if (norm[1] == 0)
    {
        double x= -this->offset / norm[0];
        double y = -(off + x * x1) / y1;
        return x * norm[0] + y * norm[1];
    }
}

/**
 * @brief Check whether the hyperplane is the same as h
 * @param h The hyperplane
 * @return  1 they are the same
 *          -1 they are different
 */
bool hyperplane::is_same(hyperplane *h)
{
    for(int i = 0; i < dim; ++i)
    {
        if (-EQN2 > (norm[i] - h->norm[i]) || (norm[i] - h->norm[i]) > EQN2)
            return false;
    }
    return true;
}










