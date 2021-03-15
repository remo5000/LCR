#include "TwoSidedBackboneIndex.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <tuple>
#include <memory>
#include <map>
#include <deque>

#define DEBUG 0
#define watch(x) if (DEBUG) cout << (#x) << " is " << (x) << endl; cout.flush();
#define watch_vector(vect) if (DEBUG) cout << (#vect) << ": [";for (int i = 0; i < vect.size(); i++) if (DEBUG) cout << vect[i] << " "; if (DEBUG) cout << "]" << endl;
#define log(x) if (DEBUG) cout << x << endl; cout.flush();
#define print(x) cout << x << endl; cout.flush();

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

bool TwoSidedBackboneIndex::computeQuery(VertexID source, VertexID target, LabelSet ls) {
    log("Starting query");
    watch(source);
    watch(target);
    watch(labelSetToLetters(ls));

    typedef unsigned int Distance;

    deque<pair<VertexID, Distance>> sourceOut;
    sourceOut.emplace_back(source, 0);
    unordered_set<VertexID> outVisited;
    for (const auto& p : sourceOut) {
        log(p.first);
    }

    deque<pair<VertexID, Distance>> targetIn;
    targetIn.emplace_back(target, 0);
    unordered_set<VertexID> inVisited;

    bool bfsOutwards = true;

    // -- Local BFS for epsilon distance --
    // BFS locally until both queues are empty
    while(
        bfsOutwards = !bfsOutwards,
        sourceOut.size() > 0 || targetIn.size() > 0)
    {
        watch(bfsOutwards);

        VertexID vertex; Distance dist;
        deque<pair<VertexID, Distance>>& q = bfsOutwards ? sourceOut : targetIn;
        unordered_set<VertexID>& visitedSet = bfsOutwards ? outVisited : inVisited;
        const unordered_set<VertexID>& otherVisitedSet = bfsOutwards ? inVisited: outVisited;

        if (DEBUG) {
            cout << (bfsOutwards ? "out" : "in")
                 << " q: [";
            for (const auto& p : q) cout << "(" << p.first << "," << p.second << "), ";
            cout << "]\n";
        }

        if (q.empty()) continue;
        std::tie(vertex, dist) = q.front();
        q.pop_front();

        // Dont search more than epsilon for local search
        if (dist > this->localSearchDistance) continue;

        // Quick return if the vertices are locally reachable
        if (otherVisitedSet.count(vertex)) return true;

        if (visitedSet.count(vertex)) continue;
        else visitedSet.insert(vertex);


        SmallEdgeSet ses;
        if (bfsOutwards) {
            this->graph->getOutNeighbours(vertex, ses);
        } else {
            this->graph->getInNeighbours(vertex, ses);
        }

        for(const auto& p : ses) {
            VertexID neighbor = p.first;
            LabelSet ls2 = p.second;
            if (!isLabelSubset(ls2, ls)) continue;
            q.push_back({neighbor, dist+1});
        }
    }

    log("Source reached:");
    for (const auto& p : outVisited) {
        log(p);
    }
    log("Target reached:");
    for (const auto& p : inVisited) {
        log(p);
    }


    log("-- starting backbone bfs --:");


    // -- Backbone BFS --

    deque<VertexID> outgoingBackboneQueue;
    for (const VertexID& vertex : outVisited)
        if (this->backboneVertices.count(vertex))
            outgoingBackboneQueue.push_back(vertex);
    deque<VertexID> incomingBackboneQueue;
    for (const VertexID& vertex : inVisited)
        if (this->backboneVertices.count(vertex))
            incomingBackboneQueue.push_back(vertex);

    outVisited.clear();
    inVisited.clear();

    bfsOutwards = true;


    while(
        bfsOutwards = !bfsOutwards,
        outgoingBackboneQueue.size() > 0 || incomingBackboneQueue.size() > 0
    ) {
        VertexID vertex;
        deque<VertexID>& q = bfsOutwards ? outgoingBackboneQueue : incomingBackboneQueue;
        unordered_set<VertexID>& visitedSet = bfsOutwards ? outVisited : inVisited;
        const unordered_set<VertexID>& otherVisitedSet = bfsOutwards ? inVisited: outVisited;

        if (DEBUG) {
            cout << (bfsOutwards ? "backbone_out" : "backbone_in")
                 << " q: [";
            for (const auto& v : q) cout << v << ", ";
            cout << "]\n";
        }

        if (q.empty()) continue;
        vertex = q.front();
        q.pop_front();

        // Quick return if the vertices are locally reachable
        if (otherVisitedSet.count(vertex)) return true;

        if (visitedSet.count(vertex)) continue;
        else visitedSet.insert(vertex);


        SmallEdgeSet ses;
        if (bfsOutwards) {
            this->backbone->getOutNeighbours(vertex, ses);
            if (vertex == 1 && DEBUG) {
                cout << "ses: [";
                for (const auto& p : ses)
                    cout << "(" << p.first << "," << labelSetToLetters(p.second) << "), ";
                cout << "]\n";
            }
        } else {
            this->backbone->getInNeighbours(vertex, ses);
        }

        for(const auto& p : ses) {
            VertexID neighbor = p.first;
            LabelSet ls2 = p.second;
            if (!isLabelSubset(ls2, ls)) continue;

            q.push_back(neighbor);
        }
    }


    return false;
}

bool TwoSidedBackboneIndex::query(VertexID source, VertexID target, LabelSet ls)
{
    cout << "TwoSidedBackboneIndex::query source=" << to_string(source) << ",target=" << to_string(target) << ",ls=" << labelSetToString(ls) << endl;
    queryStart = getCurrentTimeInMilliSec();

    bool b = this->computeQuery(source, target, ls);
    watch(b);

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

    log("Starting BFS from every vertex");
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
            if (inPath) continue;

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

            // watch(vertex);
            // watch(labelSetToString(ls));
            // watch(dist);
            // watch_vector(path);

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
    }

    return {twoSidedReachability, candidates};
}


// Used for set cover
typedef pair<VertexID, std::pair<VertexID, LabelSet>> Item;

set<Item> reachabilityToSet(const LabelledDistancedReachabilityMap& reachabilityMap) {
    set<Item> result;
    for (const auto& p : reachabilityMap.m) {
        VertexID source = p.first;
        for (const auto& p2 : p.second) {
            VertexID dest = p2.first;
            for (auto p3 : p2.second) {
                const LabelSet& ls = p3.first;
                result.insert(make_pair(source, make_pair(dest, ls)));
            }
        }
    }
    return result;
}

void TwoSidedBackboneIndex::buildIndex()
{
    Graph* graph = this->graph;
    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    // hasBeenIndexed = dynamic_bitset<>(N);

    print("Generating ground set");
    LabelledDistancedReachabilityMap groundSetMap;
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidatesToReachabilityMap;
    std::tie(groundSetMap, candidatesToReachabilityMap) = generateGroundSetAndCandidates(graph, this->localSearchDistance);
    log("generated ground set");


    // Generate candidates
    // u -> reachability info that u covers
    log("generating candidates");
    unordered_map<VertexID, set<Item>> candidates;
    for (const auto& p : candidatesToReachabilityMap) {
        VertexID candidate = p.first;
        candidates[candidate] = reachabilityToSet(p.second);
    }
    log("generated candidates. Candidates size:");
    watch(candidates.size());

    // Compute the uncovered vertices
    LabelledDistancedReachabilityMap& uncovered = groundSetMap;

    // Compute the backbone vertices
    print("Starting set cover");
    while (uncovered.size() > 0) {
        log("Uncovered:");
        if (DEBUG) cout << uncovered.toString();
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

    // source -> dest -> {LS} in backbone
    // Minimal in the sense that if u -L-> w -L-> v, and u,w,v are all in backbone,
    // u-L->w is not included.
     LabelledDistancedReachabilityMap backboneReachability;


    // Compute the backbone reachability (to compute edges)
    print("Computing backbone");
    unsigned int DIST_NOT_USED = 0;
    for (const VertexID& source : backboneVertices) {
        // Keep a local reachability for the source vertex
        LabelledDistancedReachabilityMap dfsReachability;

        typedef vector<VertexID> Path;
        vector<tuple<VertexID, LabelSet, Path>> stack;

        Path startingPath = {source};
        stack.emplace_back(source, 0, startingPath);

        while(stack.size()) {
            VertexID vertex; LabelSet ls; Path path;
            std::tie(vertex, ls, path) = stack.back();
            stack.pop_back();

            // If we are in a loop, quit (there must be a more efficient reachability path)
            bool inPath = false;
            for (int i = 0; i < path.size()-1; i++) if (path[i] == vertex) inPath = true;
            if (inPath) continue;

            // If we have already visited this node, continue. Else, visit it.
            if (dfsReachability.isPresent(source, vertex, ls)) continue;
            dfsReachability.insert(source, vertex, ls, DIST_NOT_USED);

            // Track reachability between backbone vertices (this is a subset of reachability)
            if (backboneVertices.count(vertex)) backboneReachability.insert(source, vertex, ls, DIST_NOT_USED);

            // If the current vertex is in the backbone, we don't need to dfs further because we will
            // obtain that reachability information eventually.
            if (vertex != source && backboneVertices.count(vertex)) continue;


            SmallEdgeSet ses;
            graph->getOutNeighbours(vertex, ses);
            for(const auto& p : ses)
            {
                VertexID neighbor = p.first;
                LabelSet ls2 = p.second;

                // Get the new LS
                LabelSet newLs = joinLabelSets(ls, ls2);

                Path newPath = path;
                newPath.push_back(neighbor);

                stack.push_back(make_tuple(neighbor, newLs, newPath));
            }
        }
    }

    // Clean up self-edges
    for (const VertexID& source : backboneVertices) backboneReachability.erase(source, source);

    log("Computed backbone. Backbone:");
    log(backboneReachability.toString());

    // Generate edges
    EdgeSet emptyEdgeSet;
    DGraph* dg = new DGraph(&emptyEdgeSet, this->graph->getNumberOfVertices(), 0, true);

    for (const auto& p : backboneReachability.toEdgeMap()) {
        for (const SmallEdge& smallEdge : p.second) {
            VertexID u = p.first;
            VertexID v = smallEdge.first;
            LabelSet ls = smallEdge.second;
            dg->addMultiEdge(u,v,ls);
        }
    }

    // Set the bacbone
    backbone = std::unique_ptr<DGraph>(dg);
    print("Done building index");
};

const unordered_set<VertexID>& TwoSidedBackboneIndex::getBackBoneVertices() const {
    return backboneVertices;
};

const DGraph& TwoSidedBackboneIndex::getBackBone() const {
    return *backbone;
};
