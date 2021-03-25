#include "BackboneIndex.h"
#include "../../Index/UnboundedLCR/BFSIndex.cc"

#include <vector>
#include <algorithm>
#include <utility>
#include <tuple>
#include <memory>
#include <math.h>
#include <map>
#include <queue>

#define DEBUG 0
#define watch(x) if (DEBUG) cout << (#x) << " is " << (x) << endl; cout.flush();
#define watch_vector(vect) if (DEBUG) cout << (#vect) << ": [";for (int i = 0; i < vect.size(); i++) if (DEBUG) cout << vect[i] << " "; if (DEBUG) cout << "]" << endl;
#define log(x) if (DEBUG) cout << x << endl; cout.flush();
#define print(x) cout << x << endl; cout.flush();

#define DIST_NOT_USED 0

using namespace backbonens;
using namespace indexns;
using namespace graphns;

typedef unsigned int Distance;
typedef vector<VertexID> Path;

// Used for set cover
typedef pair<VertexID, std::pair<VertexID, LabelSet>> Item;



BackboneIndex::BackboneIndex(
    Graph* mg,
    unsigned int localSearchDistance,
    BackboneVertexSelectionMethod backboneVertexSelectionMethod,
    BackboneEdgeCreationMethod backboneEdgeCreationMethod,
    BackboneIndexingMethod backboneIndexingMethod,
    LocalSearchMethod localSearchMethod
) {
    this->graph = mg;
    this->localSearchDistance = localSearchDistance;
    this->indexType = IndexType::Backbone;

    this->indexDirection = BOTHINDEX;
    this->isBlockedMode = false; // TODO this has something to do with input

    // Set parameters
    this->backboneVertexSelectionMethod = backboneVertexSelectionMethod;
    this->backboneEdgeCreationMethod = backboneEdgeCreationMethod;
    this->backboneIndexingMethod = backboneIndexingMethod;
    this->localSearchMethod = localSearchMethod;

    // Construct index
    this->didComplete = false;
    this->buildIndex();
    this->didComplete = true;
}

unsigned long BackboneIndex::getIndexSizeInBytes()
{
    int N = graph->getNumberOfVertices();
    unsigned long size = getIndexSizeInBytesM();

    size += this->backboneVertices.size()*sizeof(VertexID);
    size += this->backbone->getGraphSizeInBytes();

    if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        size += this->backboneTransitiveClosure.getSizeInBytes();
    }

    size += this->graph->getGraphSizeInBytes();
    return size;
};

bool BackboneIndex::uniDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls) {
        queue<VertexID> sourceOut;
        sourceOut.push(source);
        for (int round = 0; round < this->localSearchDistance+1; round++) {
            for (int itemsInCurrentLevel = sourceOut.size(); itemsInCurrentLevel > 0; itemsInCurrentLevel--) {
                VertexID vertex = sourceOut.front();
                sourceOut.pop();

                if (vertex == target) return true;

                for(const auto& p : this->graph->getOutNeighbours(vertex)) {
                    VertexID neighbor = p.first;
                    LabelSet ls2 = p.second;

                    if (!isLabelSubset(ls2, ls)) continue;
                    if (this->backboneVertices.count(neighbor)) continue;

                    sourceOut.push(neighbor);
                }
            }
        }
        return false;
};

bool BackboneIndex::biDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls) {
        queue<VertexID> sourceOut;
        sourceOut.push(source);
        int outRounds = 0;
        int maxOutRounds = (this->localSearchDistance/2 + this->localSearchDistance%2);

        queue<VertexID> targetIn;
        targetIn.push(target);
        int inRounds = 0;
        int maxInRounds = this->localSearchDistance/2 ;
        dynamic_bitset<> outVisited = dynamic_bitset<>(this->graph->getNumberOfVertices());
        dynamic_bitset<> inVisited = dynamic_bitset<>(this->graph->getNumberOfVertices());

        // Increment once as adding the source should not count as 1 round
        maxOutRounds++;
        maxInRounds++;

        while (
            outRounds < maxOutRounds || inRounds < maxInRounds)
        {
            if (outRounds < maxOutRounds) {
                for (int itemsInCurrentLevel = sourceOut.size(); itemsInCurrentLevel > 0; itemsInCurrentLevel--) {
                    VertexID vertex = sourceOut.front();
                    sourceOut.pop();

                    if (outVisited[vertex]) continue;
                    else outVisited[vertex] = 1;

                    if (inVisited[vertex]) return true;

                    for(const auto& p : this->graph->getOutNeighbours(vertex)) {
                        VertexID neighbor = p.first;
                        LabelSet ls2 = p.second;

                        if (!isLabelSubset(ls2, ls)) continue;
                        if (this->backboneVertices.count(neighbor)) continue;

                        sourceOut.push(neighbor);
                    }
                }
                outRounds++;
            }

            if (inRounds < maxInRounds) {
                for (int itemsInCurrentLevel = targetIn.size(); itemsInCurrentLevel > 0; itemsInCurrentLevel--) {
                    VertexID vertex = targetIn.front();
                    targetIn.pop();

                    if (inVisited[vertex]) continue;
                    else inVisited[vertex] = 1;

                    if (outVisited[vertex]) return true;

                    for(const auto& p : this->graph->getInNeighbours(vertex)) {
                        VertexID neighbor = p.first;
                        LabelSet ls2 = p.second;

                        if (!isLabelSubset(ls2, ls)) continue;
                        if (this->backboneVertices.count(neighbor)) continue;

                        targetIn.push(neighbor);
                    }
                }
                inRounds++;
            }
        }

        return false;
};

bool BackboneIndex::bfsLocally(VertexID source, VertexID target, LabelSet ls) {
    if (this->localSearchMethod == LocalSearchMethod::UNIDIRECTIONAL_BFS) {
        return this->uniDirectionalLocalBfs(source, target, ls);
    } else if (this->localSearchMethod == LocalSearchMethod::BIDIRECTIONAL_BFS) {
        return this->biDirectionalLocalBfs(source, target, ls);
    } else {
        print("Unsupported localSearchMethod. Backtrace here to check how it happened.");
        exit(1);
    }
}

bool BackboneIndex::bfsBackbone(
    VertexID source,
    VertexID target,
    LabelSet ls
) {

    log("-- starting backbone BFS --:");

    unsigned int pos;
    queue<VertexID> outgoingBackboneQueue;
    if (this->backboneReachableOut.count(source)) {
        const SmallEdgeSet& outSes = this->backboneReachableOut.at(source);
        pos = 0;
        while(pos < outSes.size()) {
            VertexID v = outSes[pos].first;
            while(pos < outSes.size() && outSes[pos].first == v) {
                if (isLabelSubset(outSes[pos].second, ls)) {
                    outgoingBackboneQueue.push(v);
                    break;
                }
                pos++;
            }
            while(pos < outSes.size() && outSes[pos].first == v) pos++;
        }
    }

    queue<VertexID> incomingBackboneQueue;
    if (this->backboneReachableIn.count(target)) {
        const SmallEdgeSet& inSes = this->backboneReachableIn.at(target);
        pos = 0;
        while(pos < inSes.size()) {
            VertexID v = inSes[pos].first;
            while(pos < inSes.size() && inSes[pos].first == v) {
                if (isLabelSubset(inSes[pos].second, ls)) {
                    incomingBackboneQueue.push(v);
                    break;
                }
                pos++;
            }
            while(pos < inSes.size() && inSes[pos].first == v) pos++;
        }
    }

    // --
    for (const auto& u : outgoingBackboneQueue)
        for (const auto& v : incomingBackboneQueue)
            if (this->backboneIndex.get()->query(source, target, ls))
                return true;
    return false;
    // --

    dynamic_bitset<> outVisited = dynamic_bitset<>(this->graph->getNumberOfVertices());
    dynamic_bitset<> inVisited = dynamic_bitset<>(this->graph->getNumberOfVertices());

    VertexID vertex;

    while(
        outgoingBackboneQueue.size() > 0 || incomingBackboneQueue.size() > 0
    ) {
        if (outgoingBackboneQueue.size()) {
            vertex = outgoingBackboneQueue.front();
            outgoingBackboneQueue.pop();

            // Quick return if the vertices are locally reachable
            if (inVisited[vertex]) return true;

            outVisited[vertex] = 1;


            const SmallEdgeSet& ses = this->backbone->getOutNeighbours(vertex);

            for(const auto& p : ses) {
                VertexID neighbor = p.first;
                LabelSet ls2 = p.second;
                if (!isLabelSubset(ls2, ls)) continue;
                if (outVisited[neighbor]) continue;

                outgoingBackboneQueue.push(neighbor);
            }
        }

        // --- Inwards ---
        if (incomingBackboneQueue.size()) {
            if (incomingBackboneQueue.empty()) continue;
            vertex = incomingBackboneQueue.front();
            incomingBackboneQueue.pop();

            // Quick return if the vertices are locally reachable
            if (outVisited[vertex]) return true;

            inVisited[vertex] = 1;


            const SmallEdgeSet& ses = this->backbone->getInNeighbours(vertex);

            for(const auto& p : ses) {
                VertexID neighbor = p.first;
                LabelSet ls2 = p.second;
                if (!isLabelSubset(ls2, ls)) continue;
                if (inVisited[neighbor]) continue;

                incomingBackboneQueue.push(neighbor);
            }
        }
    }

    return false;
}

bool BackboneIndex::backboneQueryTransitiveClosure(
    VertexID source,
    VertexID target,
    LabelSet ls
) {

    log("-- starting backbone TC query --:");

    unsigned int pos;
    vector<VertexID> outgoingBackboneQueue;
    const SmallEdgeSet& outSes = this->backboneReachableOut[source];
    pos = 0;
    while(pos < outSes.size()) {
        VertexID v = outSes[pos].first;
        while(pos < outSes.size() && outSes[pos].first == v) {
            if (isLabelSubset(outSes[pos].second, ls)) {
                outgoingBackboneQueue.push_back(v);
                break;
            }
            pos++;
        }
        while(pos < outSes.size() && outSes[pos].first == v) pos++;
    }


    vector<VertexID> incomingBackboneQueue;
    const SmallEdgeSet& inSes = this->backboneReachableIn[target];
    pos = 0;
    while(pos < inSes.size()) {
        VertexID v = inSes[pos].first;
        while(pos < inSes.size() && inSes[pos].first == v) {
            if (isLabelSubset(inSes[pos].second, ls)) {
                incomingBackboneQueue.push_back(v);
                break;
            }
            pos++;
        }
        while(pos < inSes.size() && inSes[pos].first == v) pos++;
    }

    for (const auto& u : outgoingBackboneQueue)
        for (const auto& v : incomingBackboneQueue)
            if (u == v || this->backboneTransitiveClosure.isPresent(u, v, ls))
                return true;
    return false;

}

bool BackboneIndex::queryBackbone(VertexID source, VertexID target, LabelSet ls) {
    if (this->backboneIndexingMethod == BackboneIndexingMethod::BFS) {
        return this->bfsBackbone(source, target, ls);
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        return this->backboneQueryTransitiveClosure(source, target, ls);
    } else {
        print("Unsupported backboneIndexingMethod. Backtrace here to check how it happened.");
        exit(1);
    }
}

bool BackboneIndex::computeQuery(VertexID source, VertexID target, LabelSet ls) {
    log("Starting query");
    watch(source);
    watch(target);
    watch(labelSetToLetters(ls));

    // Quick checks
    if (source == target) return true;
    if (ls == 0) return false;

    if (this->bfsLocally(source, target, ls)) return true;

    return this->queryBackbone(source, target, ls);
}

bool BackboneIndex::query(VertexID source, VertexID target, LabelSet ls)
{
    // cout << "BackboneIndex::query source=" << to_string(source) << ",target=" << to_string(target) << ",ls=" << labelSetToString(ls) << endl;
    queryStart = getCurrentTimeInMilliSec();

    bool b = this->computeQuery(source, target, ls);

    queryEndTime = getCurrentTimeInMilliSec();
    // cout << "BackboneIndex::query answer =" << b << endl;
    return b;
}

// For some reason this is blank for a lot of methods -- we shall leave it blank for now.
void BackboneIndex::queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach)
{

};

// Generate all vertex pairs that are distance epsilon+1 distance apart,
// and generate candidate vertices for the backbone set in the process.
// This step is easier to not separate due to the fact that we would do almost-repeated computation otherwise.
pair<LabelledDistancedReachabilityMap, unordered_map<VertexID, LabelledDistancedReachabilityMap>>
generateGroundSetAndCandidates(Graph* graph, unsigned int localSearchDistance) {
    // Used to return the ground set
    LabelledDistancedReachabilityMap localMeetingReachability;
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidates;

    log("Starting BFS from every vertex");
    for (VertexID source = 0; source < graph->getNumberOfVertices(); source++) {
        // Used to track visited nodes from the source
        LabelledDistancedReachabilityMap reachability;

        vector< tuple<VertexID, LabelSet, Distance, Path> > stack;
        Path startingPath = {source};
        LabelSet startingLabelSet = 0;
        assert(isEmptyLabelSet(0));
        stack.push_back(make_tuple(source, startingLabelSet, 0, startingPath));

        // watch(source);
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
                localMeetingReachability.insert(source, vertex, ls, dist);
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
            const SmallEdgeSet& ses = graph->getOutNeighbours(vertex);
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

    return {localMeetingReachability, candidates};
}

void BackboneIndex::localMeetingCriteriaSetCover() {
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
        candidates[candidate] = p.second.toTuples();
    }
    log("generated candidates. Candidates size:");
    watch(candidates.size());

    // Compute the uncovered vertices
    LabelledDistancedReachabilityMap& uncovered = groundSetMap;

    // Compute the backbone vertices
    print("Starting set cover");
    while (uncovered.size() > 0) {
        // log("Uncovered:");
        // if (DEBUG) cout << uncovered.toString();
        // log("----:");

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

        // watch(biggestCover.size());
        // watch(biggestCoverVertex);
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

}

void BackboneIndex::oneSideConditionCover() {
    struct sort_pred
    {
        bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
        {
            return left.second > right.second;
        }
    };
    vector< pair< VertexID, int > > degreePerNode;
    for(int i = 0; i < this->graph->getNumberOfVertices(); i++)
    {
        degreePerNode.push_back(
            make_pair(
                i,
                // Use product of inDeg and outDeg
                this->graph->getOutNeighbours(i).size() *  this->graph->getInNeighbours(i).size()
            )
        );
    }
    sort(degreePerNode.begin(), degreePerNode.end(), sort_pred());

    int quotum = degreePerNode.size()/20;
    if (!quotum) quotum++;

    for (int i = 0; i < degreePerNode.size(); i++) {

        if( ((i+1)%quotum) == 0 && i > 0 )
        {
            double perc = i;
            perc /= degreePerNode.size();
            perc *= 100.0;
            double timePassed = getCurrentTimeInMilliSec()-constStartTime;
            cout << this->name << "::oneSideConditionCover " << perc << "%" << ", time(s)=" << (timePassed) << endl;
        }

        const auto& p = degreePerNode[i];
        const VertexID source = p.first;

        LabelledDistancedReachabilityMap depthMap;
        LabelledDistancedReachabilityMap distanceMap;
        LabelledDistancedReachabilityMap visitedMap;

        depthMap.insert(source, source, 0, 0);
        distanceMap.insert(source, source, 0, 0);

        deque<pair<VertexID, LabelSet>> q;
        q.push_front(make_pair(source, 0));

        while (!q.empty()) {
            VertexID vertex;
            LabelSet ls;
            std::tie(vertex, ls) = q.front();
            q.pop_front();

            if (visitedMap.isPresent(source, vertex, ls)) continue;
            else visitedMap.insert(source, vertex, ls, DIST_NOT_USED);

            if (
                depthMap.getDistance(source, vertex, ls) == this->localSearchDistance
                // TODO check if == works here
                && distanceMap.getDistance(source, vertex, ls) >= this->localSearchDistance) {
                this->backboneVertices.insert(source);
                break;
            }

            // Don't go beyond epsilon in depth
            if (depthMap.getDistance(source, vertex, ls) == this->localSearchDistance)
                continue;
            for (const auto& se : this->graph->getOutNeighbours(vertex)) {
                VertexID neighbor = se.first;
                LabelSet ls2 = se.second;
                LabelSet newLs = joinLabelSets(ls, ls2);

                unsigned int neighborDepth = min(
                    depthMap.getDistance(source, neighbor, newLs),
                    depthMap.getDistance(source, vertex, ls) + 1
                );
                depthMap.insert(source, neighbor, newLs, neighborDepth);

                unsigned int neighborDist = min(
                    distanceMap.getDistance(source, neighbor, newLs),
                    distanceMap.getDistance(source, vertex, ls) + 1
                );
                distanceMap.insert(
                    source,
                    neighbor,
                    newLs,
                    backboneVertices.count(neighbor)
                    ? 0
                    : neighborDist
                );

                q.push_back(make_pair(neighbor, newLs));
            }
        }
    }
};

void BackboneIndex::selectBackboneVertices() {
    if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA) {
        this->localMeetingCriteriaSetCover();
    } else if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION) {
        this->oneSideConditionCover();
    } else {
        print("Unsupported backboneVertexSelectionMethod. Backtrace here to check how it happened.");
        exit(1);
    }
};

void BackboneIndex::createBackboneEdges() {
    if (this->backboneEdgeCreationMethod == BackboneEdgeCreationMethod::BFS) {
        // source -> dest -> {LS} in backbone
        // Minimal in the sense that if u -L-> w -L-> v, and u,w,v are all in backbone,
        // u-L->w is not included.
        LabelledDistancedReachabilityMap backboneReachability;

        // Compute the backbone reachability (to compute edges)
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


                SmallEdgeSet ses = graph->getOutNeighbours(vertex);
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

        print("Computed backbone vertices. |V*|:");
        print(backboneVertices.size());

        log("Backbone:");
        log(backboneReachability.toString());

        // Generate edges
        this->backboneEdgeSet = std::unique_ptr<EdgeSet>(new EdgeSet());

        // Set the bacbone
        this->backbone = std::unique_ptr<DGraph>(new DGraph(this->backboneEdgeSet.get(), this->graph->getNumberOfVertices(), 0, true));

        for (const auto& p : backboneReachability.toEdgeMap()) {
            for (const SmallEdge& smallEdge : p.second) {
                VertexID u = p.first;
                VertexID v = smallEdge.first;
                LabelSet ls = smallEdge.second;
                this->backbone->addMultiEdge(u,v,ls);
            }
        }
    } else {
        print("Unsupported backboneEdgeCreationMethod. Backtrace here to check how it happened.");
        exit(1);
    }
};

void BackboneIndex::indexBackbone() {
    if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        print("Computing Backbone TC");
        for (const VertexID& source : backboneVertices) {
            // Keep a local reachability for the source vertex
            LabelledDistancedReachabilityMap dfsReachability;

            vector<tuple<VertexID, LabelSet>> stack;
            stack.emplace_back(source, 0);

            while(stack.size()) {
                VertexID vertex; LabelSet ls;
                std::tie(vertex, ls) = stack.back();
                stack.pop_back();

                dfsReachability.insert(source, vertex, ls, DIST_NOT_USED);

                // Track reachability between backbone vertices (this is a subset of reachability)
                if (backboneVertices.count(vertex))
                    backboneTransitiveClosure.insert(source, vertex, ls, DIST_NOT_USED);



                const SmallEdgeSet& ses = this->graph->getOutNeighbours(vertex);
                for(const auto& p : ses)
                {
                    VertexID neighbor = p.first;
                    LabelSet ls2 = p.second;

                    // Get the new LS
                    LabelSet newLs = joinLabelSets(ls, ls2);

                    if (dfsReachability.isPresent(source, neighbor, newLs)) continue;

                    stack.push_back(make_tuple(neighbor, newLs));
                }
            }
        }
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::BFS) {
        this->backboneIndex = std::unique_ptr<Index>(new BFSIndex(this->graph));
    } else {
        print("Unsupported backboneIndexingMethod. Backtrace here to check how it happened.");
        exit(1);
    }
};

void BackboneIndex::cacheVertexToBackboneReachability() {
    // Inwards
    LabelledDistancedReachabilityMap inReachability;
    for (VertexID source : this->backboneVertices)
    {

        queue<Triplet> q;
        Triplet t;
        t.x = source;
        t.ls = 0;
        t.dist = 0;
        q.push(t);


        while(q.size()) {
            const Triplet triplet = q.front();
            VertexID vertex = triplet.x;
            LabelSet ls = triplet.ls;
            Distance dist = triplet.dist;
            q.pop();


            inReachability.insert(vertex, source, ls, DIST_NOT_USED);

            // Dont add neighbours if its the final vertex in the frontier
            if (dist == this->localSearchDistance) continue;

            for (SmallEdge se : this->graph->getOutNeighbours(vertex)) {
                Triplet newTriplet;
                newTriplet.x = se.first;
                newTriplet.ls = joinLabelSets(ls, se.second);
                newTriplet.dist = dist+1;
                q.push(newTriplet);
            }
        }

    }
    this->backboneReachableIn = inReachability.toEdgeMap();

    // Outwards
    LabelledDistancedReachabilityMap outReachability;
    for (VertexID target : this->backboneVertices)
    {
        queue<Triplet> q;
        Triplet t;
        t.x = target;
        t.ls = 0;
        t.dist = 0;
        q.push(t);


        while(q.size()) {
            const Triplet triplet = q.front();
            VertexID vertex = triplet.x;
            LabelSet ls = triplet.ls;
            Distance dist = triplet.dist;
            q.pop();

            outReachability.insert(vertex, target, ls, DIST_NOT_USED);

            // Dont add neighbours if its the final vertex in the frontier
            if (dist == this->localSearchDistance) continue;

            for (SmallEdge se : this->graph->getInNeighbours(vertex)) {
                Triplet newTriplet;
                newTriplet.x = se.first;
                newTriplet.ls = joinLabelSets(ls, se.second);
                newTriplet.dist = dist+1;
                q.push(newTriplet);
            }
        }
    }
    this->backboneReachableOut = outReachability.toEdgeMap();
};


void BackboneIndex::buildIndex()
{
    Graph* graph = this->graph;
    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    watch(localSearchDistance);

    constStartTime = getCurrentTimeInMilliSec();

    print("Selecting backbone vertices");
    selectBackboneVertices();
    print("Computing backbone edges");
    createBackboneEdges();
    print("Indexing backbone");
    indexBackbone();


    print("Caching vertex->backbone reachability");
    cacheVertexToBackboneReachability();

    constEndTime = getCurrentTimeInMilliSec();
    totalConstTime = constEndTime - constStartTime;

    print("Done building index");
};

const unordered_set<VertexID>& BackboneIndex::getBackBoneVertices() const {
    return backboneVertices;
};

const DGraph& BackboneIndex::getBackBone() const {
    return *backbone;
};
