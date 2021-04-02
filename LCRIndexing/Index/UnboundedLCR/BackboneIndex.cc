#include "BackboneIndex.h"
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

    N = this->graph->getNumberOfVertices();

    // Backbone-specific
    this->backboneTransitiveClosure = LabelledDistancedReachabilityMap(N, N);

    // Set parameters
    this->backboneVertexSelectionMethod = backboneVertexSelectionMethod;
    this->backboneEdgeCreationMethod = backboneEdgeCreationMethod;
    this->backboneIndexingMethod = backboneIndexingMethod;
    this->localSearchMethod = localSearchMethod;

    // Construct index
    this->didComplete = false;
    this->buildIndex();
    this->didComplete = true;

    // Update name
    this->name = "Backbone(";
    if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA) {
        this->name += "LOCAL_MEETING_CRITERIA";
    } else if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER) {
        this->name += "ONE_SIDE_CONDITION_DEGREE_ORDER";
    } else if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER) {
        this->name += "ONE_SIDE_CONDITION_RANDOM_ORDER";
    } else {
        this->name += "UNKNOWN";
    }
    this->name += ", ";
    if (this->backboneIndexingMethod == BackboneIndexingMethod::BFS) {
        this->name += "BFS";
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        this->name += "TRANSITIVE_CLOSURE";
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS) {
        this->name += "LANDMARK_NO_EXTENSIONS";
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS) {
        this->name += "LANDMARK_ALL_EXTENSIONS";
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_FULL) {
        this->name += "LANDMARK_FULL";
    } else {
        this->name += "UNKNOWN";
    }
    this->name += ")";
}

unsigned long BackboneIndex::getIndexSizeInBytes()
{
    int N = graph->getNumberOfVertices();
    unsigned long size = getIndexSizeInBytesM();

    size += this->backboneVertices.size()*sizeof(VertexID);
    size += this->backbone->getGraphSizeInBytes();

    if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        size += this->backboneTransitiveClosure.getSizeInBytes();
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS 
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_FULL
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS) {
        size += this->backboneLi->getIndexSizeInBytes();
        size -= this->backbone->getGraphSizeInBytes();
        size -= this->graph->getGraphSizeInBytes();
    }

    size += sizeof(this->backboneReachableOut);
    for (const auto& p : this->backboneReachableOut) {
        size += sizeof(p.first);
        size += sizeof(p.second);
        for (const SmallEdge& se : p.second) {
            size += sizeof(se.first);
            size += sizeof(se.second);
        }
    }

    size += sizeof(this->backboneReachableIn);
    for (const auto& p : this->backboneReachableIn) {
        size += sizeof(p.first);
        size += sizeof(p.second);
        for (const SmallEdge& se : p.second) {
            size += sizeof(se.first);
            size += sizeof(se.second);
        }
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
        dynamic_bitset<> outVisited = dynamic_bitset<>(N);
        dynamic_bitset<> inVisited = dynamic_bitset<>(N);

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
    queue<VertexID> q;
    if (backboneVertices.count(source)) {
        q.push(source);
    } else {
        const SmallEdgeSet& outSes = this->backboneReachableOut[source];
        pos = 0;
        while(pos < outSes.size()) {
            VertexID v = outSes[pos].first;
            while(pos < outSes.size() && outSes[pos].first == v) {
                if (isLabelSubset(outSes[pos].second, ls)) {
                    q.push(v);
                    break;
                }
                pos++;
            }
            while(pos < outSes.size() && outSes[pos].first == v) pos++;
        }
    }

    unordered_set<VertexID> targets;
    if (this->backboneVertices.count(target)) {
        targets.insert(target);
    } else {
        const SmallEdgeSet& inSes = this->backboneReachableIn[target];
        pos = 0;
        while(pos < inSes.size()) {
            VertexID v = inSes[pos].first;
            while(pos < inSes.size() && inSes[pos].first == v) {
                if (isLabelSubset(inSes[pos].second, ls)) {
                    targets.insert(v);
                    break;
                }
                pos++;
            }
            while(pos < inSes.size() && inSes[pos].first == v) pos++;
        }
    }

    dynamic_bitset<> marked = dynamic_bitset<>(this->backbone->getNumberOfVertices());

    while (q.empty() == false) {
        VertexID x = q.front();
        q.pop();

        if( targets.count(x) )
            return true;

        if( marked[x] == 1 )
        {
            continue;
        }
        marked[x] = 1;

        for(const auto& se : this->backbone->getOutNeighbours(x)) {
            if( isLabelSubset(se.second, ls) == true )
            {
                q.push( se.first );
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

bool BackboneIndex::backboneQueryLandmarks(VertexID source, VertexID target, LabelSet ls) {
    log("-- starting backbone TC query --:");

    unsigned int pos;
    vector<VertexID> sources;
    const SmallEdgeSet& outSes = this->backboneReachableOut[source];
    pos = 0;
    while(pos < outSes.size()) {
        VertexID v = outSes[pos].first;
        while(pos < outSes.size() && outSes[pos].first == v) {
            if (isLabelSubset(outSes[pos].second, ls)) {
                sources.push_back(v);
                break;
            }
            pos++;
        }
        while(pos < outSes.size() && outSes[pos].first == v) pos++;
    }


    unordered_set<VertexID> targets;
    const SmallEdgeSet& inSes = this->backboneReachableIn[target];
    pos = 0;
    while(pos < inSes.size()) {
        VertexID v = inSes[pos].first;
        while(pos < inSes.size() && inSes[pos].first == v) {
            if (isLabelSubset(inSes[pos].second, ls)) {
                targets.insert(v);
                break;
            }
            pos++;
        }
        while(pos < inSes.size() && inSes[pos].first == v) pos++;
    }

    return this->backboneLi->query(sources, targets, ls);
};

bool BackboneIndex::queryBackbone(VertexID source, VertexID target, LabelSet ls) {
    if (this->backboneIndexingMethod == BackboneIndexingMethod::BFS) {
        return this->bfsBackbone(source, target, ls);
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        return this->backboneQueryTransitiveClosure(source, target, ls);
    } else if (
            this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_FULL) {
        return this->backboneQueryLandmarks(source, target, ls);
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
BackboneIndex::generateGroundSetAndCandidates() {
    // Used to return the ground set
    int N = graph->getNumberOfVertices();
    LabelledDistancedReachabilityMap localMeetingReachability(N);
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidates;

    auto constStartTime = getCurrentTimeInMilliSec();
    log("Starting BFS from every vertex");
    for (VertexID source = 0; source < graph->getNumberOfVertices(); source++) {
        if (source % (max(N/20, 1)) == 0) {
            double perc = source;
            perc /= graph->getNumberOfVertices();
            perc *= 100;
            double timePassed = getCurrentTimeInMilliSec()-constStartTime;
            cout << this->name << "::generateGroundSetAndCandidates " << perc << "%" << ", time(s)=" << (timePassed) << endl;
        }
        // Used to track visited nodes from the source
        LabelledDistancedReachabilityMap reachability(N);

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
    LabelledDistancedReachabilityMap groundSetMap(N);
    unordered_map<VertexID, LabelledDistancedReachabilityMap> candidatesToReachabilityMap;
    std::tie(groundSetMap, candidatesToReachabilityMap) = this->generateGroundSetAndCandidates();
    log("generated ground set");


    // Generate candidates
    // u -> reachability info that u covers
    log("generating candidates");
    unordered_map<VertexID, set<Item>> candidates;
    for (const auto& p : candidatesToReachabilityMap) {
        VertexID candidate = p.first;
        candidates[candidate] = p.second.toSetCoverItems();
    }
    log("generated candidates. Candidates size:");
    watch(candidates.size());

    // Compute the uncovered vertices
    set<Item> uncovered = groundSetMap.toSetCoverItems();

    print("Starting set cover");

    // For printing
    unsigned int originalSize = uncovered.size();
    unsigned int counts = 0;
    auto constStartTime = getCurrentTimeInMilliSec();

    while (uncovered.size() > 0) {
        if (counts++ % 10 == 0) {
            double timePassed = getCurrentTimeInMilliSec()-constStartTime;
            double perc = uncovered.size();
            perc /= originalSize;
            perc *= 100;
            perc = 100 - perc;
            cout << this->name 
                << "::localMeetingCriteriaSetCover " 
                << perc << "%" 
                << ", |Uncovered|=" << (uncovered.size()) 
                << ", time(s)=" << (timePassed) 
                << endl;
        }

        // log("Uncovered:");
        // if (DEBUG) cout << uncovered.toString();
        // log("----:");

        VertexID biggestCoverVertex;
        set<Item> biggestCover;

        vector<VertexID> toDelete;

        for (auto& p : candidates) {
            VertexID vertex = p.first;
            set<Item>& coveredByVertex = p.second;

            // Count intersection
            set<Item> newCoveredItems;
            for (const Item& item : coveredByVertex) {
                if (uncovered.count(item)) {
                    newCoveredItems.insert(item);
                }
            }

            // Delete useless candidates
            if (newCoveredItems.size() == 0) {
                toDelete.push_back(vertex);
                continue;
            }

            // Update the candidate's cover to a subset that is required
            coveredByVertex = newCoveredItems;

            // Update the best cover
            if (newCoveredItems.size() > biggestCover.size()) {
                biggestCover = newCoveredItems; // TODO evaluate if we should use `move` here
                biggestCoverVertex = vertex;
            }
        }

        // Prune empty candidates
        for (VertexID v : toDelete) candidates.erase(v);

        for (const Item& item : biggestCover) {
            uncovered.erase(item);
        }
        // Add the biggest vertex to the backbone vertices
        backboneVertices.insert(biggestCoverVertex);
        candidates.erase(biggestCoverVertex);
    }

    log("backboneVertices:");
    for (auto v : backboneVertices) {
        log(v);
    }

}


inline vector<VertexID> BackboneIndex::getVerticesInDegreeOrder() {
    struct sort_pred
    {
        bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
        {
            return left.second > right.second;
        }
    };
    vector< pair< VertexID, int > > degreePerNode;
    for(int i = 0; i < N; i++)
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

    vector<VertexID> res;
    for (const auto& p : degreePerNode) res.push_back(p.first);
    return res;
}

inline vector<VertexID> BackboneIndex::getVerticesInRandomOrder() {
    vector<VertexID> res;
    for(VertexID i = 0; i < N; i++) res.push_back(i);
    std::random_shuffle ( res.begin(), res.end() );
    return res;
}

void BackboneIndex::oneSideConditionCover() {

    print("    Ordering vertices...");
    vector<VertexID> vertices;
    if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER) {
        vertices = this->getVerticesInDegreeOrder();
    } else if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER) {
        vertices = this->getVerticesInRandomOrder();
    } else {
        print("Unsupported backboneVertexSelectionMethod. Backtrace here to check how it happened.");
        exit(1);
    }

    int quotum = vertices.size()/20;
    if (!quotum) quotum++;
    auto constStartTime = getCurrentTimeInMilliSec();

    for (int i = 0; i < vertices.size(); i++) {

        if( (i%quotum) == 0)
        {
            double perc = i;
            perc /= vertices.size();
            perc *= 100.0;
            double timePassed = getCurrentTimeInMilliSec()-constStartTime;
            cout << this->name << "::oneSideConditionCover " << perc << "%" << ", time(s)=" << (timePassed) << endl;
        }

        const VertexID& source = vertices[i];

        // Use unordered_map because 
	// dense graphs have big neighborhoods.
	// Even for small graphs, will be ~= 3^k <= 1000 most of the time.
        unordered_map<VertexID, vector<LabelSet>> m;

        queue<pair<VertexID, LabelSet>> q;
        q.push(make_pair(source, 0));

        for (int round = 0; round < this->localSearchDistance; round++) {
            for (int popRound = 0; popRound < q.size(); popRound++) {
                VertexID vertex;
                LabelSet ls;
                std::tie(vertex, ls) = q.front();
                q.pop();

                // Check ls
                if (!tryInsertLabelSet(ls, m[vertex])) continue;

                for (const auto& se : this->graph->getOutNeighbours(vertex)) {
                    const VertexID& neighbor = se.first;
                    LabelSet ls2 = se.second;

                    // Check if backbone
                    if (backboneVertices.count(neighbor)) continue;

                    q.push(make_pair(neighbor, joinLabelSets(ls, ls2)));
                }
            }
        }
        if (q.size())
            backboneVertices.insert(source);
    }
};

void BackboneIndex::selectBackboneVertices() {
    if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA) {
        this->localMeetingCriteriaSetCover();
    } else if (
            this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER
            || this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER) {
        this->oneSideConditionCover();
    } else {
        print("Unsupported backboneVertexSelectionMethod. Backtrace here to check how it happened.");
        exit(1);
    }
};

void BackboneIndex::createBackboneEdges() {
    if (this->backboneEdgeCreationMethod == BackboneEdgeCreationMethod::BFS) {

        // Set backbone
        EdgeSet* emptyEdgeSet = new EdgeSet();
        DGraph* dg = new DGraph(emptyEdgeSet, N, 0, true);
        backbone = std::unique_ptr<DGraph>(dg);

        // From each source in the backbone, do a BFS on the original graph until you hit a backbone.
        for (const VertexID& source : backboneVertices) {
            LabelledReachMap backboneReachability(this->graph->getNumberOfLabels());
            LabelledReachMap graphReachability(this->graph->getNumberOfLabels());

            queue<pair<VertexID, LabelSet>> q;
            q.push(make_pair(source, 0));

            while (q.empty() == false) {
                VertexID vertex;
                LabelSet ls;
                std::tie(vertex, ls) = q.front();
                q.pop();

                if (graphReachability.isPresent(vertex, ls))
                    continue;
                else
                    graphReachability.insert(vertex, ls);

                if (vertex != source && backboneVertices.count(vertex)) {
                    backboneReachability.insert(vertex, ls);
                    continue;
                }

                for(const auto& p : graph->getOutNeighbours(vertex)) {
                    VertexID neighbor = p.first;
                    LabelSet ls2 = p.second;
                    LabelSet newLs = joinLabelSets(ls, ls2);

                    q.push(make_pair(neighbor, newLs));
                }
            }

            // Add the (minimal) reachability information to the graph
            graphReachability.addToGraphWithSource(source, backbone.get());
        }

    } else {
        print("Unsupported backboneEdgeCreationMethod. Backtrace here to check how it happened.");
        exit(1);
    }
};

void BackboneIndex::indexBackbone() {
    // TODO use the backbone graph for this && rework.
    if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        print("Computing Backbone TC");
        for (const VertexID& source : backboneVertices) {
            // Keep a local reachability for the source vertex
            Tuples dfsReachability;

            vector<tuple<VertexID, LabelSet>> stack;
            stack.emplace_back(source, 0);

            while(stack.size()) {
                VertexID vertex; LabelSet ls;
                std::tie(vertex, ls) = stack.back();
                stack.pop_back();

                // Insert -ls->vertex to dfsReachability
                int pos;
                if (indexns::findTupleInTuples(vertex, dfsReachability, pos)) {
                    continue;
                } else {
                    indexns::LabelSets lss = indexns::LabelSets();
                    lss.reserve( this->graph->getNumberOfLabels() * 2 );
                    indexns::Tuple newTuple = make_pair(vertex, lss );
                    dfsReachability.insert( dfsReachability.begin() + pos, newTuple );
                }

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

                    stack.push_back(make_tuple(neighbor, newLs));
                }
            }
        }
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::BFS) {
        // Do nothing
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS) {
        int k = 1250 + sqrt(this->backboneVertices.size());
        if (k >= this->backboneVertices.size())
            k = sqrt(this->backboneVertices.size());
        int b = 20;
        this->backboneLi = unique_ptr<LandmarkedIndex>(new LandmarkedIndex(this->backbone.get(), false, false, k, b));
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS) {
        int k = 1250 + sqrt(this->backboneVertices.size());
        if (k >= this->backboneVertices.size())
            k = sqrt(this->backboneVertices.size());
        int b = 20;
        this->backboneLi = unique_ptr<LandmarkedIndex>(new LandmarkedIndex(this->backbone.get(), true, true, k, b));
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_FULL) {
        this->backboneLi = unique_ptr<LandmarkedIndex>(new LandmarkedIndex(this->backbone.get(), false, false, this->backbone->getNumberOfVertices(), 0));
    } else {
        print("Unsupported backboneIndexingMethod. Backtrace here to check how it happened.");
        exit(1);
    }
};

void BackboneIndex::cacheVertexToBackboneReachability() {

    // Inwards
    for (VertexID source : this->backboneVertices)
    {
        unordered_map<VertexID, vector<LabelSet>> inReachability;

        queue<SmallEdge> q;
        q.push(make_pair(source, 0));

        for (int round = 0; round < this->localSearchDistance+1; round++) {
            bool reachedLocalSearchDistance = round == this->localSearchDistance;
            for (int popRound = 0; popRound < q.size(); popRound++) {
                VertexID vertex;
                LabelSet ls;
                tie(vertex, ls) = q.front();
                q.pop();

                if (!tryInsertLabelSet(ls, inReachability[vertex]))
                    continue;

		// if (vertex != source && backboneVertices.count(vertex))
		// 	continue;

                if (reachedLocalSearchDistance) continue;
                for (SmallEdge se : this->graph->getOutNeighbours(vertex)) {
                    VertexID neighbor = se.first;
                    LabelSet newLs = joinLabelSets(ls, se.second);
                    q.push(make_pair(neighbor, newLs));
                }
            }
        }

        for (const auto& p : inReachability) {
            const VertexID& vertex = p.first;
            for (const LabelSet& ls : p.second) {
                this->backboneReachableIn[vertex].push_back(make_pair(source, ls));
            }
        }

    }

    // Outwards
    for (VertexID target : this->backboneVertices)
    {
        unordered_map<VertexID, vector<LabelSet>> outReachability;

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

	    if (!tryInsertLabelSet(ls, outReachability[vertex])) continue;

	//     if (vertex != target && backboneVertices.count(vertex))
	// 	    continue;

            if (dist == this->localSearchDistance) continue;

            for (SmallEdge se : this->graph->getInNeighbours(vertex)) {
                Triplet newTriplet;
                newTriplet.x = se.first;
                newTriplet.ls = joinLabelSets(ls, se.second);
                newTriplet.dist = dist+1;
                q.push(newTriplet);
            }
        }

	for (const auto& p : outReachability) {
	    const VertexID& vertex = p.first;
	    for (const LabelSet& ls : p.second) {
		this->backboneReachableOut[vertex].push_back(make_pair(target, ls));
	    }
	}
    }
};


void BackboneIndex::buildIndex()
{
    Graph* graph = this->graph;
    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();
    watch(localSearchDistance);

    constStartTime = getCurrentTimeInMilliSec();

    print("Using localDistance of ");
    print(this->localSearchDistance);

    print("Selecting backbone vertices");
    selectBackboneVertices();
    print("Computed backbone vertices. |V*|:");
    print(backboneVertices.size());
    print("|V*|/|V| :");
    print((double)backboneVertices.size() / (double)N);

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
