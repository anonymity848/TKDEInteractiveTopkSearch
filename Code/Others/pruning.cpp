#include "pruning.h"

char hidden_options[] = " d n v Qbb QbB Qf Qg Qm Qr QR Qv Qx Qz TR E V Fa FA FC FD FS Ft FV Gt Q0 Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 ";

#ifdef WIN32
#ifdef __cplusplus
extern "C" {
#endif
#endif

//#include "data_utility.h"

#include "../qhull/mem.h"
#include "../qhull/qset.h"
#include "../qhull/libqhull.h"
#include "../qhull/qhull_a.h"
#include "../qhull/io.h"

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

/**
 * @brief Conduct halfspace intersection by invoking Qhull
 * @param rPtr   Contains the data of all the halfspace
 * @param wPtr   The extreme points of the intersection and the necessary halfspaces are shown in wPtr
 */
int halfspace(FILE *rPtr, FILE *wPtr)
{
    int curlong, totlong; /* used !qh_NOmem */
    int exitcode, numpoints, dim;
    coordT *points;
    boolT ismalloc;

    // the required parameters
    int argc = 3;
    char *argv[3];
    argv[0] = (char*)"qhalf";
    argv[1] = (char*)"Fp";
    argv[2] = (char*)"Fx";

    qh_init_A(rPtr, wPtr, stderr, argc, argv);  /* sets qh qhull_command */
    exitcode = setjmp(qh errexit); /* simple statement for CRAY J916 */
    if (!exitcode)
    {
        qh NOerrexit = False;
        qh_option("Halfspace", NULL, NULL);
        qh HALFspace = True;    /* 'H'   */
        qh_checkflags(qh qhull_command, hidden_options);
        qh_initflags(qh qhull_command);

        points = qh_readpoints(&numpoints, &dim, &ismalloc);

        if (dim >= 5)
        {
            qh_option("Qxact_merge", NULL, NULL);
            qh MERGEexact = True; /* 'Qx' always */
        }
        qh_init_B(points, numpoints, dim, ismalloc);
        qh_qhull();
        qh_check_output();
        qh_produce_output();
        //print_summary();

        if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
        {
            qh_check_points();
        }
        exitcode = qh_ERRnone;
    }
    qh NOerrexit = True;  /* no more setjmp */
#ifdef qh_NOmem
    qh_freeqhull(qh_ALL);
#else
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);
    if (curlong || totlong)
    {
        fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
                totlong, curlong);
    }
#endif
    return exitcode;
}

/**
 * @brief Use the branch-and-bound skyline (BBS) algorithm for maintaining the candidate set
 * @param pset      The point set
 * @param C_idx     The indexes of the current candidate favorite car in pset
 * @param ext_pts   The set of extreme points of R
 */
void rtree_pruning(point_set *pset, vector<int> &C_idx, Partition* R)
{
    if(C_idx.size() <= 1)
        return;

    vector<point_t *> hyperplanes;
    hyperplane *hp = NULL;

    // parameters for building the R-trees
    rtree_info *aInfo;
    aInfo = (rtree_info *) malloc(sizeof(rtree_info));
    memset(aInfo, 0, sizeof(rtree_info));
    aInfo->m = 18;
    aInfo->M = 36;  //number of children
    aInfo->dim = pset->points[0]->dim;
    aInfo->reinsert_p = 27;
    aInfo->no_histogram = C_idx.size();

    // construct R-tree
    node_type *root = contructRtree(pset, C_idx, aInfo);
    priority_queue<node_type *, vector<node_type *>, nodeCmp> heap;

    heap.push(root);
    int *sl = new int[C_idx.size()];
    int index = 0;
    int dim = aInfo->dim;

    // run the adapted BBS algorithm
    while (!heap.empty())
    {
        node_type *n = heap.top();
        heap.pop();
        if (n->attribute != LEAF)
        {
            int dominated = 0;
            vector<point_t*> TRpt; //find the bounding box of the set of points
            for(int i = 0; i < pow(2, dim); ++i)
                TRpt.push_back(new point_t(dim));
            for (int j = 0; j < dim; j++)
            {
                int groupSize = pow(2, dim)/pow(2, j + 1);
                for(int i = 0; i < pow(2, dim); ++i)
                {
                    if( (i / groupSize + 1) % 2 == 0)
                        TRpt[i]->attr[j] = n->a[j];
                    else
                        TRpt[i]->attr[j] = n->b[j];
                }
            }

            // check if TRpt is dominated by other points
            for (int j = 0; j < index && !dominated; ++j)
            {
                dominated = 1;
                for(int i = 0; i <  TRpt.size(); ++i)
                {
                    if(!R->R_dominate(pset->points[sl[j]], TRpt[i]))
                    {
                        dominated = 0;
                        break;
                    }
                }
            }

            if (!dominated)
            {
                for (int i = 0; i < aInfo->M - n->vacancy; i++)
                    heap.push(n->ptr[i]);
            }
        }
        else
        {
            int idx = n->id;
            //S = updateS(id, C, S, V);
            int dominated = 0;
            for (int j = 0; j < index && !dominated; ++j)
            {
                if(R->R_dominate(pset->points[sl[j]], pset->points[C_idx[idx]]))
                    dominated = 1;
            }
            if (dominated)
                continue;

            // eliminate any points in current Q' that it is dominated by the new inserted point
            int m = index;
            index = 0;
            for (int j = 0; j < m; ++j)
            {
                if(!R->R_dominate(pset->points[C_idx[idx]], pset->points[sl[j]]))
                    sl[index++] = sl[j];
            }

            // add this point as well
            sl[index++] = C_idx[idx];
        }
    }

    // clean up
    C_idx.clear();
    for (int i = 0; i < index; i++)
        C_idx.push_back(sl[i]);
    delete[] sl;
    free(aInfo);
}








































