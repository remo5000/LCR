#include "TwoSidedBackboneIndex.h"

#include <vector>
#include <algorithm>
#include <utility>
#include <tuple>

using namespace twosidedbackbonens;
using namespace indexns;
using namespace graphns;


TwoSidedBackboneIndex::TwoSidedBackboneIndex(Graph* mg, unsigned int localSearchDistance) {
    this->graph = mg;
    this->localSearchDistance = localSearchDistance;
    this->indexType = IndexType::TwoSidedBackbone;

    this->indexDirection = BOTHINDEX;
    this->isBlockedMode = false; // TODO this has something to do with input

    // Construct index
    this->didComplete = false;
    this->buildIndex();
    this->didComplete = true;
}

unsigned long TwoSidedBackboneIndex::getIndexSizeInBytes()
{
    return getIndexSizeInBytesM();
};

bool TwoSidedBackboneIndex::query(VertexID source, VertexID target, LabelSet ls)
{
    cout << "TwoSidedBackboneIndex::query source=" << to_string(source) << ",target=" << to_string(target) << ",ls=" << labelSetToString(ls) << endl;
    queryStart = getCurrentTimeInMilliSec();
    bool b = true; // TODO
    queryEndTime = getCurrentTimeInMilliSec();
    cout << "TwoSidedBackboneIndex::query answer =" << b << endl;
    return b;
}

// For some reason this is blank for a lot of methods -- we shall leave it blank for now.
void TwoSidedBackboneIndex::queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach)
{

};

// Generate all vertex pairs that are distance epsilon+1 distance apart.
LabelledDistancedReachabilityMap generateGoundSet(Graph* graph, unsigned int localSearchDistance) {
    // Used to return the ground set
    LabelledDistancedReachabilityMap twoSidedReachability;
    // Used to track visited nodes
    LabelledDistancedReachabilityMap reachability;

    typedef unsigned int Distance;

    for (VertexID source = 0; source < graph->getNumberOfVertices(); source++) {
        vector< tuple<VertexID, LabelSet, Distance> > stack;
        assert(isEmptyLabelSet(0));
        stack.push_back(make_tuple(source,0,0));

        while (!stack.empty()) {
            VertexID vertex; LabelSet ls; Distance dist;
            std::tie(vertex, ls, dist) = stack.back();
            stack.pop_back();

            // Stop BFS-ing if dist > epsilon+1
            if (dist > localSearchDistance+1) continue;
            // Add this to the resultant ground set if the distance is exactly epsilon+1
            if (dist == localSearchDistance+1) {
                twoSidedReachability.insert(source, vertex, ls, dist);
                continue;
            }
            // else, continue the BFS.

            // If we have already visited this node, continue. Else, visit it.
            if (reachability.isPresent(source, vertex, ls) || dist >= reachability.getDistance(source, vertex, ls)) continue;
            reachability.insert(source, vertex, ls, dist);

            // Add all neighbours
            SmallEdgeSet ses;
            graph->getOutNeighbours(vertex, ses);
            for(const auto& p : ses)
            {
                VertexID neighbor = p.first;
                LabelSet ls2 = p.second;

                // Get the new LS
                LabelSet newLs = joinLabelSets(ls, ls2);
                // Get the (possibly) shorter distance
                unsigned int newDist = min(reachability.getDistance(source, neighbor, newLs), reachability.getDistance(source, vertex, newLs)+1);

                stack.push_back(make_tuple(neighbor, newLs, newDist));
            }
        }
    }

    return twoSidedReachability;
}

void TwoSidedBackboneIndex::buildIndex()
{
    Graph* graph = this->graph;
    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    // hasBeenIndexed = dynamic_bitset<>(N);

    LabelledReachabilityMap groundSet = generateGoundSet(graph, this->localSearchDistance);
};
