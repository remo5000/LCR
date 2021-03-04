#include "../../Graph/Graph.h"
#include "Index.h"

#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>

using namespace indexns;
using namespace graphns;

#ifndef TWOSIDEDBACKBONEINDEX_H
#define TWOSIDEDBACKBONEINDEX_H

namespace twosidedbackbonens {
    class LabelledDistancedReachabilityMap {
        public:
            void insert(VertexID source, VertexID destination, LabelSet ls, unsigned int distance) {
                if (m[source][destination].count(ls)) {
                    this->m[source][destination][ls] = min(this->m[source][destination][ls], distance);
                    return;
                } else {
                    this->m[source][destination][ls] = distance;
                }
            }
            bool isPresent(VertexID source, VertexID destination, LabelSet ls) {
                if (!m.count(source)) return false;
                if (!m[source].count(destination)) return false;
                for (auto p : m[source][destination]) {
                    if (isLabelSubset(ls, p.first)) {
                        return true;
                    }
                }
                return false;
            }
            unsigned int getDistance(VertexID source, VertexID destination, LabelSet ls) {
                if (this->isPresent(source, destination, ls)) {
                    return this->m[source][destination][ls];
                } else {
                    return std::numeric_limits<unsigned int>::max();
                }
            }
            // source -> dest -> ls -> dist
            unordered_map<VertexID, unordered_map<VertexID, unordered_map<LabelSet, unsigned int> > > m;
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
        unsigned int localSearchDistance;
        void buildIndex();
};
#endif
