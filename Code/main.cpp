#include "structure/data_utility.h"
#include "structure/data_struct.h"
#include "structure/point_set.h"
#include "structure/define.h"
#include <vector>
#include <ctime>
#include <sys/time.h>
#include "UH/UH.h"
#include "GroundTruth/Groundtruth.h"
#include "2D/PI.h"
#include "HD/RH/RH.h"
#include "ActiveRanking/ActiveRanking.h"
#include "PreferenceLearning/preferenceLearning.h"
#include "UH/UH.h"
#include "UH/Median_Hull.h"
#include "UtilityApprox/UtilityApprox.h"
#include "HD/HDPI/HDPI.h"
#include "HD/Unify/Unify.h"
#include "Scale/Scale.h"


int main(int argc, char *argv[])
{
    long mem_baseline = get_mem_usage();
    ifstream config("../config.txt");
    string alg_name;
    char data_name[MAX_FILENAME_LENG]; double k, s;
    config >> alg_name >> data_name >> k >> s;
    cout << alg_name << "   " << data_name << "   " << k << "   " << s << "\n";
    //initialization
    //point set
    point_set *originalSet = new point_set(data_name);
    point_set *pset = new point_set();
    originalSet->skyband(pset, k);
    point_t *u = new point_t(pset->points[0]->dim);
    for(int i = 0; i < pset->points[0]->dim; ++i)
        config >> u->attr[i];

    printf("-----------------------------------------------------------------------------------\n");
    printf("|%15s |%15s |%15s |%15s |%10s |\n", "Algorithm", "# of Questions", "Preprocessing", "Interaction", "Point #ID");
    printf("-----------------------------------------------------------------------------------\n");
    ground_truth(pset, u, k); //look for the ground truth maximum utility point

    int error = 1;

    if(alg_name == "partitionBS")
        PartitionBS(pset, u, k, s);
    else if(alg_name == "segmentBS")
        SegmentBS(pset, u, k, s);
    else if(alg_name == "Unify11")
        unify(pset, u, k, s, 1, 1);
    else if(alg_name == "Unify12")
        unify(pset, u, k, s, 1, 2);
    else if(alg_name == "Unify13")
        unify(pset, u, k, s, 1, 3);
    else if(alg_name == "Unify21")
        unify(pset, u, k, s, 2, 1);
    else if(alg_name == "Unify22")
        unify(pset, u, k, s, 2, 2);
    else if(alg_name == "Unify23")
        unify(pset, u, k, s, 2, 3);
    else if(alg_name == "Unify31")
        unify(pset, u, k, s, 3, 1);
    else if(alg_name == "Unify32")
        unify(pset, u, k, s, 3, 2);
    else if(alg_name == "Unify33")
        unify(pset, u, k, s, 3, 3);
    else if(alg_name == "Unify41")
        unify(pset, u, k, s, 4, 1);
    else if(alg_name == "Unify42")
        unify(pset, u, k, s, 4, 2);
    else if(alg_name == "Unify43")
        unify(pset, u, k, s, 4, 3);
    else if(alg_name == "Unify11Error")
        unify(pset, u, k, s, 1, 1, error);
    else if(alg_name == "Unify12Error")
        unify(pset, u, k, s, 1, 2, error);
    else if(alg_name == "Unify13Error")
        unify(pset, u, k, s, 1, 3, error);
    else if(alg_name == "Unify21Error")
        unify(pset, u, k, s, 2, 1, error);
    else if(alg_name == "Unify22Error")
        unify(pset, u, k, s, 2, 2, error);
    else if(alg_name == "Unify23Error")
        unify(pset, u, k, s, 2, 3, error);
    else if(alg_name == "Unify31Error")
        unify(pset, u, k, s, 3, 1, error);
    else if(alg_name == "Unify32Error")
        unify(pset, u, k, s, 3, 2, error);
    else if(alg_name == "Unify33Error")
        unify(pset, u, k, s, 3, 3, error);
    else if(alg_name == "Median")
        Median(pset, u, k, s);
    else if(alg_name == "Hull")
        Hull(pset, u, k, s);
    else if(alg_name == "ActiveRanking")
        ActiveRanking(pset, u, k, s, error);
    else if(alg_name == "PreferenceLearning")
        PreferenceLearning(pset, u, k, s, error);
    else if(alg_name == "UHRandom")
        UHRandom(pset, u, k, s, error);
    else if(alg_name == "UHSimplex")
        UHSimplex(pset, u, k, s, error);
    //else if(alg_name == "HDPI")
    //HDPI_grid(pset, u, k, s);
    //HDPI_extreme(pset, u, k, s);







    delete pset;
    return 0;
}
