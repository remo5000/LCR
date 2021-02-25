#include "../../Graph/Graph.h"
#include "Index.h"

#include <set>

using namespace indexns;
using namespace graphns;

#ifndef TWOSIDEDBACKBONEINDEX_H
#define TWOSIDEDBACKBONEINDEX_H

namespace twosidedbackbonens {
    class LabelledReachabilityMap {
        public:
            void insert(VertexID source, VertexID destination, LabelSet ls);
            bool isPresent(VertexID source, VertexID destination, LabelSet ls);
    };
}

class TwoSidedBackboneIndex : public Index
{
    public:

        TwoSidedBackboneIndex(Graph* mg, unsigned int localSearchDistance);

        // To implement Index
        bool query(VertexID source, VertexID target, LabelSet ls);
        unsigned long getIndexSizeInBytes();
        void queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach);

    private:
        void buildIndex();
};
#endif
