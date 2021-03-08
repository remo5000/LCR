#include "TwoSidedBackboneIndex.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <tuple>
#include <map>

#define watch(x) cout << (#x) << " is " << (x) << endl
#define watch_vector(vect) cout << (#vect) << ": [";for (int i = 0; i < vect.size(); i++) cout << vect[i] << " ";cout << "]" << endl;

#define log(x) cout << x << endl

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
            if (reachability.isPresent(source, vertex, ls) && dist >= reachability.getDistance(source, vertex, ls)) continue;
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

void generateCandidatesDfs(
    VertexID node,
    Graph* graph,
    vector<VertexID>& path,
    unordered_map<VertexID, LabelledDistancedReachabilityMap>& candidates,
    unsigned int localSearchDistance,
    VertexID dest,
    LabelSet ls
)
{
    // Dist too long
    if (path.size() > localSearchDistance+2) {
        return;
    }

    // Add middle nodes to candidate set
    if (path.size() == localSearchDistance+2 && node == dest) {
        for (int i = 1; i < path.size()-1; ++i) {
            VertexID source = path[0];
            VertexID candidateNode = path[i];
            candidates[candidateNode].insert(source, dest, ls, localSearchDistance+2);
        }
    }

    SmallEdgeSet ses;
    graph->getOutNeighbours(node, ses);
    for(const auto& p : ses) {
        VertexID neighbor = p.first;
        LabelSet ls2 = p.second;
        if (!isLabelSubset(ls, ls2)) {
            continue;
        }
        if (find(path.begin(), path.end(), neighbor) != path.end()) {
            continue;
        }

        path.push_back(neighbor);
        generateCandidatesDfs(
            neighbor,
            graph,
            path,
            candidates,
            localSearchDistance,
            dest,
            ls
        );
        path.pop_back();
    }
}


unordered_map<VertexID, LabelledDistancedReachabilityMap> generateCandidates(
    Graph* graph,
    unsigned int localSearchDistance,
    const LabelledDistancedReachabilityMap& twoSidedReachability
)
{
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidates;
    for (const auto& p1 :  twoSidedReachability.m) {
        VertexID source = p1.first;
        for (const auto& p2 : p1.second) {
            VertexID dest = p2.first;
            for (const auto& p3 : p2.second) {
                LabelSet ls = p3.first;
                assert(p3.second == localSearchDistance+1);
                cout << "start dfs " << source << " -> " << dest << " via " << labelSetToString(ls) << "." << endl;

                vector<VertexID> path = {source};
                generateCandidatesDfs(
                    source,
                    graph,
                    path,
                    candidates,
                    localSearchDistance,
                    dest,
                    ls
                );
            }
        }
    }

    return candidates;
}


// Used for set cover
typedef pair<VertexID, std::pair<VertexID, LabelSet>> Item;

// TODO refactor names
set<Item> reachabilityToSet(const LabelledDistancedReachabilityMap& groundSetMap) {
    set<Item> uncovered;
    for (const auto& p : groundSetMap.m) {
        VertexID source = p.first;
        for (const auto& p2 : p.second) {
            VertexID dest = p2.first;
            for (auto p3 : p2.second) {
                const LabelSet& ls = p3.first;
                uncovered.insert(make_pair(source, make_pair(dest, ls)));
            }
        }
    }
    return uncovered;
}

void TwoSidedBackboneIndex::buildIndex()
{
    Graph* graph = this->graph;
    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    // hasBeenIndexed = dynamic_bitset<>(N);

    LabelledDistancedReachabilityMap groundSetMap = generateGoundSet(graph, this->localSearchDistance);


    // Generate candidates
    // u -> reachability info that u covers
    const unordered_map<VertexID, LabelledDistancedReachabilityMap>& candidatesToReachabilityMap = generateCandidates(
        graph,
        localSearchDistance,
        groundSetMap
    );
    unordered_map<VertexID, set<Item>> candidates;
    for (const auto& p : candidatesToReachabilityMap) {
        VertexID candidate = p.first;
        candidates[candidate] = reachabilityToSet(p.second);
    }

    // Compute the uncovered vertices
    set<Item> uncovered = reachabilityToSet(groundSetMap);

    unordered_set<VertexID> backboneVertices;
    while (uncovered.size() > 0) {
        VertexID biggestCoverVertex;
        set<Item> biggestCover;
        for (const auto& p : candidates) {
            VertexID vertex = p.first;
            set<Item> coveredByVertex = p.second;

            set<Item> newCoveredItems;
            for (const Item& item : coveredByVertex) {
                if (!uncovered.count(item)) {
                    newCoveredItems.insert(item);
                }
            }

            if (newCoveredItems.size() > biggestCover.size()) {
                biggestCover = newCoveredItems; // TODO evaluate if we should use `move` here
                biggestCoverVertex = vertex;
            }
        }

        for (const Item& item : biggestCover) {
            uncovered.erase(item);
        }
        backboneVertices.insert(biggestCoverVertex);
    }
};
