#include "TwoSidedBackboneIndex.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <tuple>
#include <map>

#define DEBUG 0
#define watch(x) if (DEBUG) cout << (#x) << " is " << (x) << endl
#define watch_vector(vect) if (DEBUG) cout << (#vect) << ": [";for (int i = 0; i < vect.size(); i++) if (DEBUG) cout << vect[i] << " "; if (DEBUG) cout << "]" << endl;
#define log(x) if (DEBUG) cout << x << endl

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

// Generate all vertex pairs that are distance epsilon+1 distance apart,
// and generate candidate vertices for the backbone set in the process.
// This step is easier to not separate due to the fact that we would do almost-repeated computation otherwise.
pair<LabelledDistancedReachabilityMap, unordered_map<VertexID, LabelledDistancedReachabilityMap>>
generateGroundSetAndCandidates(Graph* graph, unsigned int localSearchDistance) {
    // Used to return the ground set
    LabelledDistancedReachabilityMap twoSidedReachability;
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidates;
    // Used to track visited nodes
    LabelledDistancedReachabilityMap reachability;

    typedef unsigned int Distance;
    typedef vector<VertexID> Path;

    for (VertexID source = 0; source < graph->getNumberOfVertices(); source++) {
        vector< tuple<VertexID, LabelSet, Distance, Path> > stack;
        Path startingPath = {source};
        LabelSet startingLabelSet = 0;
        assert(isEmptyLabelSet(0));
        stack.push_back(make_tuple(source, startingLabelSet, 0, startingPath));
        watch(source);

        while (!stack.empty()) {
            VertexID vertex; LabelSet ls; Distance dist; Path path;
            std::tie(vertex, ls, dist, path) = stack.back();
            stack.pop_back();

            // Stop BFS-ing if dist > epsilon+1
            if (dist > localSearchDistance+1) continue;

            // If we are in a loop, quit (there must be a more efficient reachability path)
            bool inPath = false;
            for (int i = 0; i < path.size()-1; i++) if (path[i] == vertex) inPath = true;
            if (inPath) break;

            // If we have already visited this node, continue. Else, visit it.
            if (reachability.isPresent(source, vertex, ls) && dist >= reachability.getDistance(source, vertex, ls)) continue;
            reachability.insert(source, vertex, ls, dist);

            // Add this to the resultant ground set if the distance is exactly epsilon+1
            if (dist == localSearchDistance+1) {
                twoSidedReachability.insert(source, vertex, ls, dist);
                for (int i = 1; i < path.size()-1; i++) {
                    VertexID intermediate = path[i];
                    candidates[intermediate].insert(source, vertex, ls, dist);
                }
                continue;
            }

            watch(vertex);
            watch(labelSetToString(ls));
            watch(dist);
            watch_vector(path);

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
                unsigned int newDist = min(reachability.getDistance(source, neighbor, newLs), dist+1);

                Path newPath = path;
                newPath.push_back(neighbor);

                stack.push_back(make_tuple(neighbor, newLs, newDist, newPath));
            }
        }
        log("\n\n");
    }

    return {twoSidedReachability, candidates};
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

    log("generating ground set");
    LabelledDistancedReachabilityMap groundSetMap;
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidatesToReachabilityMap;
    std::tie(groundSetMap, candidatesToReachabilityMap) = generateGroundSetAndCandidates(graph, this->localSearchDistance);
    log("generated ground set");


    // Generate candidates
    // u -> reachability info that u covers
    log("generating candidates");
    map<VertexID, set<Item>> candidates;
    for (const auto& p : candidatesToReachabilityMap) {
        VertexID candidate = p.first;
        candidates[candidate] = reachabilityToSet(p.second);
    }
    log("generated candidates");
    watch(candidates.size());
    watch(candidatesToReachabilityMap.size());

    // TODO remove
    if (DEBUG)
    for (const auto& p : candidates) {
        watch(p.first);
        for (const auto& item : p.second) {
            cout << "  ("  << item.first << "->" << item.second.first << ", " << labelSetToLetters(item.second.second) << ")\n";
        }
    }

    // Compute the uncovered vertices
    LabelledDistancedReachabilityMap& uncovered = groundSetMap;

    log("starting set cover");
    while (uncovered.size() > 0) {
        // TODO update to use LabelledDistancedReachabilityMap before using
        log("Uncovered:");
        watch(uncovered.size());
        // TODO FDSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSDJFIWHWEIUHSUDFL
        if (DEBUG)
        for (const auto& p1 : uncovered.m) {
            for (const auto& p2 : p1.second) {
                for (const auto& p3 : p2.second) {
                    cout << "  ("  << p1.first;
                    cout << "->" << p2.first << ", ";
                    cout << labelSetToString(p3.first) << ")\n";
                }
            }
        }
        log("----:");

        VertexID biggestCoverVertex;
        set<Item> biggestCover;
        for (const auto& p : candidates) {
            VertexID vertex = p.first;
            const set<Item>& coveredByVertex = p.second;

            set<Item> newCoveredItems;
            for (const Item& item : coveredByVertex) {
                if (uncovered.isPresent(item.first, item.second.first, item.second.second)) {
                    newCoveredItems.insert(item);
                }
            }

            if (newCoveredItems.size() > biggestCover.size()) {
                biggestCover = newCoveredItems; // TODO evaluate if we should use `move` here
                biggestCoverVertex = vertex;
            }
        }

        watch(biggestCover.size());
        watch(biggestCoverVertex);
        for (const Item& item : biggestCover) {
            uncovered.erase(item.first, item.second.first, item.second.second);
        }
        // Add the biggest vertex to the backbone vertices
        backboneVertices.insert(biggestCoverVertex);
    }
    log("backboneVertices:");
    for (auto v : backboneVertices) {
        log(v);
    }
};

const unordered_set<VertexID>& TwoSidedBackboneIndex::getBackBoneVertices() const {
    return backboneVertices;
};
