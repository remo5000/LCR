#include "../../Graph/Graph.h"
#include "Index.h"

using namespace indexns;
using namespace graphns;

#ifndef TWOHOPINDEX_H
#define TWOHOPINDEX_H

class TwoHopIndex : public Index
{
    public:

        TwoHopIndex(Graph* mg);
        ~TwoHopIndex() = default;

        bool query(VertexID source, VertexID target, LabelSet ls);
        void queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach);
        unsigned long getIndexSizeInBytes();

    private:
        int visitedSetSize;
        void buildIndex();
        bool computeQuery(VertexID source, VertexID target, LabelSet ls);

	indexns::TuplesList inIndex;
	indexns::TuplesList outIndex;
	unsigned int N;
};
#endif
