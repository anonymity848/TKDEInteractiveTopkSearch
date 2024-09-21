#ifndef RUN_CHOOSE_ITEM_H
#define RUN_CHOOSE_ITEM_H
#include "hyperplane.h"

class choose_item
{
public:
    hyperplane *h;
    std::vector<int> positiveSide;
    std::vector<int> negativeSide;
    std::vector<int> intersectCase;
};


#endif //RUN_CHOOSE_ITEM_H
