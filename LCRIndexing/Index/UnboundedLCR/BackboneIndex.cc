#include "BackboneIndex.h"
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <tuple>
#include <memory>
#include <math.h>
#include <map>
#include <queue>

#include <iostream>
#include <utility>
#include <thread>
#include <chrono>

#include <omp.h>


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
    this->isBlockedMode = false; 

    N = this->graph->getNumberOfVertices();

    // Backbone-specific
    this->backboneTransitiveClosure = LabelledDistancedReachabilityMap(N, N);

    this->bfsBackboneTargets = dynamic_bitset<>(N);

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


    // size += this->graph->getGraphSizeInBytes(); Exclude as this is not part of the index itself.
    size += this->backboneVertices.size()*sizeof(VertexID);
    size += this->backbone->getGraphSizeInBytes();

    print("Backbone vertices and graph size");
    print(size);

    if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
        size += this->backboneTransitiveClosure.getSizeInBytes();
    } else if (this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS 
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_FULL
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS) {
        size += this->backboneLi->getIndexSizeInBytes();
    }

    unsigned long beforeCache = size;
    for (const auto& p : this->backboneReachableOut) {
        for (const auto& se : p) {
            size += sizeof(se.first);
	    for (const auto& ls : se.second) {
		size += sizeof(ls);
	    }
        }
    }

    for (const auto& p : this->backboneReachableIn) {
        for (const auto& se : p) {
            size += sizeof(se.first);
	    for (const auto& ls : se.second) {
		size += sizeof(ls);
	    }
        }
    }

    print("Cache size:");
    print(size - beforeCache);

    return size;
};

inline bool BackboneIndex::uniDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls) {
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
                    if (this->isBackboneVertex[neighbor]) continue;

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
                        if (this->isBackboneVertex[neighbor]) continue;

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
                        if (this->isBackboneVertex[neighbor]) continue;

                        targetIn.push(neighbor);
                    }
                }
                inRounds++;
            }
        }

        return false;
};

inline bool BackboneIndex::bfsLocally(VertexID source, VertexID target, LabelSet ls) {
    if (this->localSearchMethod == LocalSearchMethod::UNIDIRECTIONAL_BFS) {
        return this->uniDirectionalLocalBfs(source, target, ls);
    } else if (this->localSearchMethod == LocalSearchMethod::BIDIRECTIONAL_BFS) {
        return this->biDirectionalLocalBfs(source, target, ls);
    } else {
        print("Unsupported localSearchMethod. Backtrace here to check how it happened.");
        exit(1);
    }
}

inline vector<VertexID> BackboneIndex::accessBackboneOut(VertexID source, LabelSet ls) {
    vector<VertexID> outgoingBackboneQueue;
	if (backboneVertices.count(source)) {
		outgoingBackboneQueue.push_back(source);
		return outgoingBackboneQueue;
	}
    for (const Tuple& tuple : this->backboneReachableOut[source]) {
	    const VertexID& v = tuple.first;

	    bool found = false;
	    for (const LabelSet& ls2 : tuple.second) {
		    if (isLabelSubset(ls2, ls)) {
			    found = true;
			    break;
		    }
	    }

	    if (found) {
		    outgoingBackboneQueue.push_back(v);
	    }
    }
    return outgoingBackboneQueue;
	
    // vector<VertexID> outgoingBackboneQueue;
    // const SmallEdgeSet& outSes = this->backboneReachableOut[source];
    // unsigned int pos = 0;
    // while(pos < outSes.size()) {
    //     VertexID v = outSes[pos].first;
    //     while(pos < outSes.size() && outSes[pos].first == v) {
    //         if (isLabelSubset(outSes[pos].second, ls)) {
    //             outgoingBackboneQueue.push_back(v);
    //             break;
    //         }
    //         pos++;
    //     }
    //     while(pos < outSes.size() && outSes[pos].first == v) pos++;
    // }
    // return outgoingBackboneQueue;
}

inline vector<VertexID> BackboneIndex::accessBackboneIn(VertexID target, LabelSet ls) {
    vector<VertexID> outgoingBackboneQueue;
	if (backboneVertices.count(target)) {
		outgoingBackboneQueue.push_back(target);
		return outgoingBackboneQueue;
	}
    for (const Tuple& tuple : this->backboneReachableIn[target]) {
	    const VertexID& v = tuple.first;

	    bool found = false;
	    for (const LabelSet& ls2 : tuple.second) {
		    if (isLabelSubset(ls2, ls)) {
			    found = true;
			    break;
		    }
	    }

	    if (found) {
		    outgoingBackboneQueue.push_back(v);
	    }
    }
    return outgoingBackboneQueue;
}


inline queue<VertexID> BackboneIndex::accessBackboneOutQueue(VertexID source, LabelSet ls) {
    queue<VertexID> outgoingBackboneQueue;
	if (backboneVertices.count(source)) {
		outgoingBackboneQueue.push(source);
		return outgoingBackboneQueue;
	}
    for (const Tuple& tuple : this->backboneReachableOut[source]) {
	    const VertexID& v = tuple.first;

	    bool found = false;
	    for (const LabelSet& ls2 : tuple.second) {
		    if (isLabelSubset(ls2, ls)) {
			    found = true;
			    break;
		    }
	    }

	    if (found) {
		    outgoingBackboneQueue.push(v);
	    }
    }
    return outgoingBackboneQueue;
}

inline unordered_set<VertexID> BackboneIndex::accessBackboneInSet(VertexID target, LabelSet ls) {
    unordered_set<VertexID> targets;
	if (backboneVertices.count(target)) {
		targets.insert(target);
		return targets;
	}
    for (const Tuple& tuple : this->backboneReachableIn[target]) {
	    const VertexID& v = tuple.first;

	    bool found = false;
	    for (const LabelSet& ls2 : tuple.second) {
		    if (isLabelSubset(ls2, ls)) {
			    found = true;
			    break;
		    }
	    }

	    if (found) {
		targets.insert(v);
	    }
    }
    return targets;
}

inline void BackboneIndex::markTargetsForBackboneBfs(VertexID target, LabelSet ls) {
    bfsBackboneTargets = dynamic_bitset<>(this->backbone->getNumberOfVertices());
	if (backboneVertices.count(target)) {
		bfsBackboneTargets[target] = 1;
		return;
	}
    for (const Tuple& tuple : this->backboneReachableIn[target]) {
	    const VertexID& v = tuple.first;

	    bool found = false;
	    for (const LabelSet& ls2 : tuple.second) {
		    if (isLabelSubset(ls2, ls)) {
			    found = true;
			    break;
		    }
	    }

	    if (found) {
		bfsBackboneTargets[v] = 1;
	    }
    }
}

inline bool BackboneIndex::queryBackbone(VertexID source, VertexID target, LabelSet ls) {

    markTargetsForBackboneBfs(target, ls);

    if (this->backboneIndexingMethod == BackboneIndexingMethod::BFS) {
		log("-- starting backbone BFS --:");

		unsigned int pos;
		queue<VertexID> q = accessBackboneOutQueue(source, ls);
		dynamic_bitset<>& targets = this->bfsBackboneTargets;
		dynamic_bitset<> marked = dynamic_bitset<>(this->backbone->getNumberOfVertices());

		while (q.empty() == false) {
			VertexID x = q.front();
			q.pop();

			if( targets[x] == 1) {
				return true;
			}

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
	} else if (this->backboneIndexingMethod == BackboneIndexingMethod::TRANSITIVE_CLOSURE) {
		vector<VertexID> outgoingBackboneQueue = accessBackboneOut(source, ls);
		vector<VertexID> incomingBackboneQueue = accessBackboneIn(target, ls);
		for (const auto& u : outgoingBackboneQueue)
			for (const auto& v : incomingBackboneQueue)
				if (u == v || this->backboneTransitiveClosure.isPresent(u, v, ls))
					return true;
		return false;
    } else if (
            this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS
            || this->backboneIndexingMethod == BackboneIndexingMethod::LANDMARK_FULL) {
		vector<VertexID> sources = accessBackboneOut(source, ls);
		unordered_set<VertexID> targets = accessBackboneInSet(target, ls);
		return this->backboneLi->query(sources, targets, bfsBackboneTargets, ls);
    } else {
        print("Unsupported backboneIndexingMethod. Backtrace here to check how it happened.");
        exit(1);
    }
}

inline bool BackboneIndex::computeQuery(VertexID source, VertexID target, LabelSet ls) {
    log("Starting query");

    // watch(source);
    // watch(target);
    // watch(labelSetToLetters(ls));

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
	isBackboneVertex[biggestCoverVertex] = 1;
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
    for(VertexID i = 0; i < N; i++) {
		res.push_back(i);
	}

	std::shuffle(std::begin(res), std::end(res), std::random_device());

    return res;
}

#define MIN_BACKBONE_RATIO 0.30
void BackboneIndex::oneSideConditionCover() {

    print("Ordering vertices...");
    vector<VertexID> vertices;
    if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER) {
        vertices = this->getVerticesInDegreeOrder();
    } else if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER) {
        vertices = this->getVerticesInRandomOrder();
    } else {
        print("Unsupported backboneVertexSelectionMethod. Backtrace here to check how it happened.");
        exit(1);
    }

    // Set bitset for fast checking
    this->isBackboneVertex = dynamic_bitset<>(N);

    int quotum = vertices.size()/100;
    if (!quotum) quotum++;
    auto constStartTime = getCurrentTimeInMilliSec();


    // Use unordered_map because 
    // dense graphs have big neighborhoods.
    // Even for small graphs, will be ~= 3^k <= 1000 most of the time.
    unordered_map<VertexID, vector<LabelSet>> m;
    // TODO initial perf analysis shows that this is good, but we need to measure more.
    size_t maxUsedSize = 4096;


    int numberOfConfirmedBackboneVertices = (int)((float)MIN_BACKBONE_RATIO * (float)vertices.size());
    // // TODO fix this for random ordering -- queries that should be true are getting false from query()
    // if (this->backboneVertexSelectionMethod == BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER) {
	// numberOfConfirmedBackboneVertices = 0;
    // }

    print("Resrving backbone vertices in advance...");
    print(numberOfConfirmedBackboneVertices);
    backboneVertices.reserve(numberOfConfirmedBackboneVertices);

    cout << this->name << "::oneSideConditionCover: begin" << endl;
    for (int i = 0; i < vertices.size(); i++) {

        if( ((i+1)%quotum) == 0)
        {
            double perc = i;
            perc /= vertices.size();
            perc *= 100.0;
            double timePassed = getCurrentTimeInMilliSec()-constStartTime;
            cout << this->name << "::oneSideConditionCover " << perc << "%" << ", time(s)=" << (timePassed) << endl;
        }

        const VertexID& source = vertices[i];

		// Add to backbone by default if required.
		if (i < numberOfConfirmedBackboneVertices) {
			backboneVertices.insert(source);
			isBackboneVertex[source] = 1;
			continue;
		}


        queue<pair<VertexID, LabelSet>> q;
        q.push(make_pair(source, 0));

	m.clear();
	m.reserve(maxUsedSize);

        for (int round = 0; round < this->localSearchDistance+1; round++) {

	    bool reachedLocalSearchDistance = round == this->localSearchDistance;

	    if (reachedLocalSearchDistance) {
		// At the last frontier, we only care about checking 
		// the vertices (as they are of depth localSearchDistance
		while (q.empty() == false) {
		    VertexID vertex;
		    LabelSet ls;
		    std::tie(vertex, ls) = q.front();
		    q.pop();

		    auto it = m.find(vertex);
		    // If we havent seen this vertex, add to backbone.
		    if (it == m.end()) {
			backboneVertices.insert(source);
			isBackboneVertex[source] = 1;
			break;
		    } 
		    // If we havent seen this vertex w.r.t this labelset, add to backbone.
		    if (tryInsertLabelSet(ls, it->second))  {
			backboneVertices.insert(source);
			isBackboneVertex[source] = 1;
			break;
		    }
		}
	    } else {
		// Before the last frontier, we simply add vertices to the queue
		// and process themn in order, the popRound forlooop ensures that 
		// we expand the frontier by one level.
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
			if (isBackboneVertex[neighbor]) continue;

			q.push(make_pair(neighbor, joinLabelSets(ls, ls2)));
		    }
		}
		maxUsedSize = max(maxUsedSize, q.size());
	    }

        }
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


void bfsFromBackboneVertex() {

}


void BackboneIndex::createBackboneEdges() {
    if (this->backboneEdgeCreationMethod != BackboneEdgeCreationMethod::BFS) {
        print("Unsupported backboneEdgeCreationMethod. Backtrace here to check how it happened.");
        exit(1);
    }

    // Set backbone
    EdgeSet* emptyEdgeSet = new EdgeSet();
    DGraph* dg = new DGraph(emptyEdgeSet, N, 0, true);
    backbone = std::unique_ptr<DGraph>(dg);

    auto constStartTime = getCurrentTimeInMilliSec();
    int quotum = backboneVertices.size()/20;
    if (!quotum) quotum++;
    int numProcessed = 0;

    // From each source in the backbone, do a BFS on the original graph until you hit a backbone.
    for (const VertexID& source : backboneVertices) {

        if( (numProcessed++%quotum) == 0  )
        {
            double perc = numProcessed;
            perc /= backboneVertices.size();
            perc *= 100.0;
            double timePassed = getCurrentTimeInMilliSec()-constStartTime;
            cout << this->name << "::createBackboneEdges " << perc << "%" << ", time(s)=" << (timePassed) << endl;
        }

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

	    if (vertex != source && isBackboneVertex[vertex]) {
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
	backboneReachability.addToGraphWithSource(source, backbone.get());
    }

    double timePassed = getCurrentTimeInMilliSec()-constStartTime;
    cout << this->name 
	<< "::createBackboneEdges" 
	<< ", time(s)=" << (timePassed) 
	<< endl;
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
                if (isBackboneVertex[vertex])
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

Graph* BackboneIndex::getGraph() {
    return  this->graph;
}

void BackboneIndex::cacheVertexToBackboneReachability() {
    auto constStartTime = getCurrentTimeInMilliSec();

    this->backboneReachableOut = TuplesList(N);
    this->backboneReachableIn = TuplesList(N);
	
	#pragma omp parallel for
    for (VertexID source = 0; source < N; source++) {
		if (backboneVertices.count(source)) 
			continue;

		for (int out = 0; out < 2; out++) {
			if (out) {
				// Outward
				queue<Triplet> q;
				Triplet t;
				t.x = source;
				t.ls = 0;
				t.dist = 0;
				q.push(t);

				unordered_map<VertexID, LabelSets> outReachability;

				while(q.size()) {
					const Triplet triplet = q.front();
					const VertexID vertex = triplet.x;
					const LabelSet ls = triplet.ls;
					const Distance dist = triplet.dist;
					q.pop();

					if (this->isBackboneVertex[vertex]) {
						auto it = outReachability.find(vertex);
						if (it == outReachability.end()) {

							// outReachability.reserve(outReachability.size()*2);

							indexns::LabelSets lss = indexns::LabelSets();
							lss.reserve( this->getGraph()->getNumberOfLabels() * 2 );
							lss.push_back(ls);

							outReachability.emplace(vertex, std::move(lss));
						} else {
							if (!tryInsertLabelSet(ls, it->second)) continue;
						}
					}


					if (dist == this->localSearchDistance) continue;

					for (const SmallEdge& se : this->getGraph()->getOutNeighbours(vertex)) {
						Triplet newTriplet;
						newTriplet.x = se.first;
						newTriplet.ls = joinLabelSets(ls, se.second);
						newTriplet.dist = dist+1;
						q.push(newTriplet);
					}
				}

				for (const auto& p : outReachability) {
					const VertexID& vertex = p.first;

					const LabelSets& lss = p.second;

					int pos = 0;
					Tuples& tuples = this->backboneReachableOut[source];
					tuples.emplace_back(vertex, std::move(lss) );
				}
			} else {
				// Inward
				queue<Triplet> q;
				Triplet t;
				t.x = source;
				t.ls = 0;
				t.dist = 0;
				q.push(t);

				unordered_map<VertexID, LabelSets> inReachability;

				while(q.size()) {
					const Triplet triplet = q.front();
					const VertexID vertex = triplet.x;
					const LabelSet ls = triplet.ls;
					const Distance dist = triplet.dist;
					q.pop();

					if (this->isBackboneVertex[vertex]) {
						auto it = inReachability.find(vertex);
						if (it == inReachability.end()) {

							// inReachability.reserve(inReachability.size()*2);

							indexns::LabelSets lss = indexns::LabelSets();
							lss.reserve( this->getGraph()->getNumberOfLabels() * 2 );
							lss.push_back(ls);

							inReachability.emplace(vertex, std::move(lss));
						} else {
							if (!tryInsertLabelSet(ls, it->second)) continue;
						}
					}

					if (dist == this->localSearchDistance) continue;

					for (const SmallEdge& se : this->getGraph()->getInNeighbours(vertex)) {
						Triplet newTriplet;
						newTriplet.x = se.first;
						newTriplet.ls = joinLabelSets(ls, se.second);
						newTriplet.dist = dist+1;
						q.push(newTriplet);
					}
				}

				for (const auto& p : inReachability) {
					const VertexID& vertex = p.first;

					const LabelSets& lss = p.second;

					int pos = 0;
					Tuples& tuples = this->backboneReachableIn[source];
					tuples.emplace_back(vertex, std::move(lss) );
				}
			}
		}
	}
	

    double timePassed = getCurrentTimeInMilliSec()-constStartTime;
    cout << this->name << "::cacheVertexToBackboneReachability " << "100%" << ", time(s)=" << (timePassed) << endl;
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
    print("Computed backbone vertices. |E*|:");
    print(this->backbone->getNumberOfEdges());
    print("|E*|/|E| :");
    print((double)this->backbone->getNumberOfEdges() / (double)this->graph->getNumberOfEdges());

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
