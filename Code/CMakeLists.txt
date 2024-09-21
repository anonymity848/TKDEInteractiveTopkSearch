cmake_minimum_required(VERSION 3.15)
project(run)

set(CMAKE_CXX_STANDARD 14)
#set(INC_DIR /usr/local/Cellar/glpk/5.0/include)
#set(LINK_DIR /usr/local/Cellar/glpk/5.0/lib)
set(INC_DIR /opt/homebrew/opt/glpk/include)
set(LINK_DIR /opt/homebrew/opt/glpk/lib)

include_directories(.)
include_directories(Others)
include_directories(qhull)
include_directories(QuadProg)
include_directories(structure)
include_directories(UH)

include_directories(${INC_DIR})
link_directories(${LINK_DIR})
link_libraries(glpk)


add_executable(run
        Others/lp.cpp
        Others/lp.h
        Others/operation.cpp
        Others/operation.h
        Others/pruning.cpp
        Others/pruning.h
        output/ext_pt.txt
        output/hyperplane_data.txt
        output/point.txt
        output/top.txt
        qhull/geom.cpp
        qhull/geom.h
        qhull/geom2.cpp
        qhull/global.cpp
        qhull/io.cpp
        qhull/io.h
        qhull/libqhull.cpp
        qhull/libqhull.h
        qhull/mem.cpp
        qhull/mem.h
        qhull/merge.cpp
        qhull/merge.h
        qhull/poly.cpp
        qhull/poly.h
        qhull/poly2.cpp
        qhull/qhull_a.h
        qhull/qhull_build.cpp
        qhull/qhull_build.h
        qhull/qset.cpp
        qhull/qset.h
        qhull/random.cpp
        qhull/random.h
        qhull/rboxlib.cpp
        qhull/stat.cpp
        qhull/stat.h
        qhull/user.cpp
        qhull/user.h
        qhull/usermem.cpp
        qhull/userprintf.cpp
        qhull/userprintf_rbox.cpp
        QuadProg/Array.cc
        QuadProg/Array.hh
        QuadProg/QuadProg++.cc
        QuadProg/QuadProg++.hh
        structure/data_struct.h
        structure/data_utility.cpp
        structure/data_utility.h
        structure/define.h
        structure/hyperplane.cpp
        structure/hyperplane.h
        structure/Partition.cpp
        structure/Partition.h
        structure/point_set.cpp
        structure/point_set.h
        structure/point_t.cpp
        structure/point_t.h
        structure/rtree.cpp
        structure/rtree.h
        UH/frame.cpp
        UH/frame.h
        UH/UH.cpp
        UH/UH.h
        main.cpp
        GroundTruth/Groundtruth.cpp
        GroundTruth/Groundtruth.h
        structure/hyperplane_set.cpp
        structure/hyperplane_set.h
        2D/PI.cpp
        2D/PI.h
        structure/u_vector.cpp
        structure/u_vector.h
        HD/RH/RH.cpp
        HD/RH/RH.h
        ActiveRanking/ActiveRanking.cpp
        ActiveRanking/ActiveRanking.h
        PreferenceLearning/preferenceLearning.cpp
        PreferenceLearning/preferenceLearning.h
        structure/s_node.cpp
        structure/s_node.h
        structure/cluster_t.cpp
        structure/cluster_t.h
        UH/Median_Hull.cpp
        UH/Median_Hull.h
        UtilityApprox/UtilityApprox.cpp
        UtilityApprox/UtilityApprox.h
        HD/HDPI/HDPI.cpp
        HD/HDPI/HDPI.h
        structure/choose_item.cpp
        structure/choose_item.h
        Scale/Scale.cpp
        Scale/Scale.h
        HD/Unify/Unify.cpp
        HD/Unify/Unify.h)

target_link_libraries(run glpk)