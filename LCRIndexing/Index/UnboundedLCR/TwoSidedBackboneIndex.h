#include "../../Graph/Graph.h"
#include "Index.h"

#include <set>

using namespace indexns;
using namespace graphns;

#ifndef TWOSIDEDBACKBONEINDEX_H
#define TWOSIDEDBACKBONEINDEX_H

/**
* This class answers query on the fly. It uses breadth first search to answers
* the queries.
*/
class TwoSidedBackboneIndex : public Index
{
    public:

        TwoSidedBackboneIndex(Graph* mg);

        bool query(VertexID source, VertexID target, LabelSet ls);
        void buildIndex();
        bool queryShell(VertexID source, VertexID target, LabelSet ls, set< VertexID >& visitedSet);
        void queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach);

        int getVisitedSetSize();
        unsigned long getIndexSizeInBytes();

    private:
        int visitedSetSize;

};
#endif
