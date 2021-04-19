#include "../../Graph/Graph.h"
#include "../../Graph/DGraph.h"
#include "../../Index/UnboundedLCR/LandmarkedIndex.cc"
#include "../../Index/UnboundedLCR/TwoHopIndex.cc"
#include "Index.h"

#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <memory>
#include <limits>
#include <string>
#include <sstream>

using namespace indexns;
using namespace graphns;

#ifndef BACKBONEINDEX_H
#define BACKBONEINDEX_H

namespace backbonens {
    typedef unsigned int Distance;
    typedef pair<LabelSet, Distance> LabelSetAndDistance;
    typedef vector<pair<LabelSet, Distance>> LabelSetAndDistanceList;
    typedef unordered_map<VertexID, LabelSetAndDistanceList> DTuples;
    typedef vector<DTuples> DTuplesList;

    class LabelledReachMap {
        public:
            explicit LabelledReachMap(int numberOfLabels) { this->numberOfLabels = numberOfLabels;
            };
            inline bool insert(VertexID destination, LabelSet ls) {
                int pos = 0;
                bool b = indexns::findTupleInTuples(destination, tuples, pos);

                if( b == false )
                {
                    // entry was not found and should be inserted
                    indexns::LabelSets lss = indexns::LabelSets();
                    lss.reserve( this->numberOfLabels * 2 );
                    indexns::Tuple newTuple = make_pair(destination, lss );
                    tuples.insert( tuples.begin() + pos, newTuple );
                }

                b = indexns::tryInsertLabelSet(ls, tuples[pos].second);
                return b;
            }
            inline bool isPresent(VertexID destination, LabelSet ls) {
                int pos = 0;
                if (indexns::findTupleInTuples(destination, tuples, pos)) {
                    const LabelSets& lss = tuples[pos].second;
                    for (int i = 0; i < lss.size(); i++) {
                        const LabelSet& ls2 = lss[i];
                        if( isLabelSubset(ls2,ls)  ) return true;
                    }
                }
                return false;
            }
            inline void erase(VertexID destination) {
                int pos;
                if(!indexns::findTupleInTuples(destination, tuples, pos))
                    return;
                if( pos < 0 || pos >= tuples.size() )
                    return;
                tuples.erase(tuples.begin() + pos);
            }
            inline const Tuples& toTuples() {
                return tuples;
            }
            inline void addToGraphWithSource(VertexID source, DGraph* graph) {
                for (const Tuple& tuple : tuples) {
                    VertexID destination = tuple.first;
                    for (const LabelSet& ls : tuple.second) {
                        graph->addMultiEdge(source, destination, ls);
                    }
                }
            }
            inline SmallEdgeSet toSmallEdgeSet() const {
                SmallEdgeSet ses;
                for (const Tuple& tuple : tuples) {
                    for (const LabelSet& ls : tuple.second) {
                        ses.push_back(make_pair(tuple.first, ls));
                    }
                }
                return ses;
            }
        private:
            Tuples tuples;
            int numberOfLabels;
    };

    // TODO either use or deprecate
    class LabelledReachabilityMap {
        public:
            explicit LabelledReachabilityMap(int numberOfLabels, int N) {
                this->numberOfLabels = numberOfLabels;
                this->tIn = indexns::TuplesList(N);
            };
            inline bool insert(VertexID source, VertexID destination, LabelSet ls) {
                int pos = 0;
                bool b = indexns::findTupleInTuples(destination, tIn[source], pos);

                if( b == false )
                {
                    // entry was not found and should be inserted
                    indexns::LabelSets lss = indexns::LabelSets();
                    lss.reserve( this->numberOfLabels * 2 );
                    indexns::Tuple newTuple = make_pair(destination, lss );
                    tIn[source].insert( tIn[source].begin() + pos, newTuple );
                }

                b = indexns::tryInsertLabelSet(ls, tIn[source][pos].second);
                return b;
            }
            inline bool isPresent(VertexID source, VertexID destination, LabelSet ls) {
                int pos = 0;
                if (indexns::findTupleInTuples(destination, tIn[source], pos)) {
                    const LabelSets& lss = tIn[source][pos].second;
                    for (int i = 0; i < lss.size(); i++) {
                        const LabelSet& ls2 = lss[i];
                        if( isLabelSubset(ls2,ls)  ) return true;
                    }
                }
                return false;
            }
            inline void erase(VertexID source, VertexID destination) {
                int pos;
                if(!indexns::findTupleInTuples(destination, tIn[source], pos))
                    return;
                if( pos < 0 || pos >= tIn[source].size() )
                    return;

                Tuples& tuples = tIn[source];
                tuples.erase(tuples.begin() + pos);
            }
            inline vector<tuple<VertexID, VertexID, LabelSet>> toMultiEdges() {
                vector<tuple<VertexID, VertexID, LabelSet>> res;
                for (VertexID u = 0; u < tIn.size(); u++) {
                    for (int j = 0; j <  tIn[u].size(); j++) {
                        const VertexID& v = tIn[u][j].first;
                        const LabelSets& lss = tIn[u][j].second;
                        for (const LabelSet& ls : lss) {
                            res.push_back(make_tuple(u,v,ls));
                        }
                    }
                }
                return res;
            }
            inline map<VertexID, SmallEdgeSet> toEdgeMap() const {
                map<VertexID, SmallEdgeSet> res;
                for (VertexID u = 0; u < tIn.size(); u++) {
                    for (int j = 0; j <  tIn[u].size(); j++) {
                        const VertexID& v = tIn[u][j].first;
                        const LabelSets& lss = tIn[u][j].second;
                        for (const LabelSet& ls : lss) {
                            res[u].push_back(make_pair(v,ls));
                        }
                    }
                }
                return res;
            }
        private:
            TuplesList tIn;
            int numberOfLabels;
    };

    // TODO add subset/superset matching to speed up
    class LabelledDistancedReachabilityMap {
        public:
            LabelledDistancedReachabilityMap(VertexID N, VertexID M) {
                this->tuplesList.resize(N);
            };
            LabelledDistancedReachabilityMap(VertexID N) {
                this->tuplesList.resize(N);
            };
            LabelledDistancedReachabilityMap(): LabelledDistancedReachabilityMap(0) {};
            ~LabelledDistancedReachabilityMap() = default;

            inline void insert(VertexID source, VertexID destination, LabelSet ls, Distance distance) {
                ensureTuplesListHasSource(source);
                assert(source < tuplesList.size());

                int pos;
                DTuples& tuples = tuplesList.at(source);

                LabelSetAndDistanceList& lsdl = tuples[destination];

                if(!findLabelSetAndDistanceInList(ls, lsdl, pos)) {
                    lsdl.push_back(make_pair(ls, distance));
                    _size++;
                }
            }
            inline bool isPresent(VertexID source, VertexID destination, LabelSet ls) {
                ensureTuplesListHasSource(source);
                assert(source < tuplesList.size());
                DTuples& tuples = tuplesList[source];

                const auto it = tuples.find(destination);
                if (it == tuples.end())
                    return false;

                int pos2 = -1;
                return findLabelSetAndDistanceInList(ls, it->second, pos2);
            }
            inline void erase(VertexID source, VertexID destination, LabelSet ls) {
                ensureTuplesListHasSource(source);
                int pos;
                DTuples& tuples = tuplesList[source];

                auto it = tuples.find(destination);
                if (it == tuples.end())
                    return;

                LabelSetAndDistanceList& lsdl = it->second;
                while(findLabelSetAndDistanceInList(ls, lsdl, pos)) {
                    lsdl.erase(lsdl.begin() + pos);
                    _size--;
                }

            }
            inline void erase(VertexID source, VertexID destination) {
                ensureTuplesListHasSource(source);
                int pos;
                DTuples& tuples = tuplesList[source];
                if (!tuples.count(destination)) {
                    return;
                }
                _size -= tuples.find(destination)->second.size();
                tuples.erase(destination);
            }
            inline unsigned int getDistance(VertexID source, VertexID destination, LabelSet ls) {
                ensureTuplesListHasSource(source);
                int pos;
                DTuples& tuples = tuplesList[source];

                const auto it = tuples.find(destination);
                if (it == tuples.end())
                    return std::numeric_limits<unsigned int>::max();

                const LabelSetAndDistanceList& lsdl = it->second;

                if(findLabelSetAndDistanceInList(ls, lsdl, pos)) {
                    return lsdl[pos].second;
                }

                return std::numeric_limits<unsigned int>::max();
            }
            inline unsigned int size() {
                return _size;
            }
            inline string toString() {
                stringstream s;
                s << "[ LabelledDistancedReachabilityMap (" << this->size() << " item(s)\n";
                for(VertexID source = 0; source < tuplesList.size(); source++) {
                    const DTuples& dtls = tuplesList[source];
                    for (const auto& dt : dtls) {
                        VertexID destination = dt.first;
                        for (const LabelSetAndDistance& lsd : dt.second) {
                            LabelSet ls = lsd.first;
                            Distance dist = lsd.second;
                            s << "  ("  << source;
                            s << "->" << destination << ", ls: ";
                            s << labelSetToLetters(ls) << ", d: ";
                            s << dist << ")\n";
                        }
                    }
                }
                s << "]\n";
                return s.str();
            }
            inline map<VertexID, SmallEdgeSet> toEdgeMap() const {
                map<VertexID, SmallEdgeSet> result;
                for(VertexID source = 0; source < tuplesList.size(); source++) {
                    const DTuples& dtls = tuplesList[source];
                    for (const auto& dt : dtls) {
                        VertexID destination = dt.first;
                        for (const LabelSetAndDistance& lsd : dt.second) {
                            LabelSet ls = lsd.first;
                            Distance dist = lsd.second;

                            result[source].emplace_back(destination, ls);
                        }
                    }
                }
                return result;
            }
            typedef pair<VertexID, std::pair<VertexID, LabelSet>> Item;
            inline set<Item> toSetCoverItems() const {
                set<Item> result;
                for(VertexID source = 0; source < tuplesList.size(); source++) {
                    const DTuples& dtls = tuplesList[source];
                    for (const auto& dt : dtls) {
                        VertexID destination = dt.first;
                        for (const LabelSetAndDistance& lsd : dt.second) {
                            LabelSet ls = lsd.first;
                            Distance dist = lsd.second;

                            result.insert({source, {destination, dist}});
                        }
                    }
                }
                return result;
            }
            inline unsigned long getSizeInBytes() {
                unsigned long size = sizeof(tuplesList);
                for(VertexID source = 0; source < tuplesList.size(); source++) {
                    const DTuples& dtls = tuplesList[source];
                    size += sizeof(dtls);
                    for (const auto& dt : dtls) {
                        VertexID destination = dt.first;
                        // Don't add pair directly, as sizeof impl can vary.
                        size += sizeof(VertexID);
                        size += dt.second.size() * sizeof(LabelSetAndDistance);
                    }
                }
                return size;
            }
        private:
            // source -> dest -> (ls, dist)
            DTuplesList tuplesList;

            inline bool isPresent(VertexID destination, const DTuples& tuples, int& pos) {
                return tuples.count(destination);
            }

            inline bool findLabelSetAndDistanceInList(LabelSet ls, const LabelSetAndDistanceList& lsdl, int& pos)
            {
                LabelSetAndDistance res = {ls, std::numeric_limits<unsigned int>::max()};
                pos = -1;
                bool found = false;

                for(int i = 0; i < lsdl.size(); i++)
                {
                    const LabelSetAndDistance& lsd = lsdl[i];

                    if( isLabelSubset(lsd.first, ls) && lsd.second <= res.second)
                    {
                        res = lsd;
                        pos = i;
                        found = true;
                    }
                }

                return found;
            }

            inline void ensureTuplesListHasSource(VertexID source) {
                if (source >= tuplesList.size()) tuplesList.resize(source+1);
            }

            unsigned int _size;
    };
}

enum class BackboneVertexSelectionMethod {
LOCAL_MEETING_CRITERIA,
ONE_SIDE_CONDITION_DEGREE_ORDER,
ONE_SIDE_CONDITION_RANDOM_ORDER
};

enum class BackboneEdgeCreationMethod {
BFS
};

enum class BackboneIndexingMethod {
BFS,
TRANSITIVE_CLOSURE,
LANDMARK_NO_EXTENSIONS,
LANDMARK_ALL_EXTENSIONS,
LANDMARK_FULL,
TWOHOP,
};

enum class LocalSearchMethod {
UNIDIRECTIONAL_BFS,
BIDIRECTIONAL_BFS
};

class BackboneIndex : public Index
{
    public:
        BackboneIndex(
            Graph* mg,
            unsigned int localSearchDistance,
            BackboneVertexSelectionMethod backboneVertexSelectionMethod,
            BackboneEdgeCreationMethod backboneEdgeCreationMethod,
            BackboneIndexingMethod backboneIndexingMethod,
            LocalSearchMethod localSearchMethod
        );
        BackboneIndex(Graph* mg, unsigned int localSearchDistance)
            : BackboneIndex(
                mg,
                localSearchDistance,
                BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER,
                BackboneEdgeCreationMethod::BFS,
                BackboneIndexingMethod::BFS,
                LocalSearchMethod::UNIDIRECTIONAL_BFS
            ) {};
        ~BackboneIndex() = default;

        // To implement Index
        bool query(VertexID source, VertexID target, LabelSet ls);
        void queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach);

        // Backbone
        const unordered_set<VertexID>& getBackBoneVertices() const;
        const DGraph& getBackBone() const;

        unsigned long getIndexSizeInBytes();

    private:

        unsigned int N;

        BackboneVertexSelectionMethod backboneVertexSelectionMethod;
        BackboneEdgeCreationMethod backboneEdgeCreationMethod;
        BackboneIndexingMethod backboneIndexingMethod;
        LocalSearchMethod localSearchMethod;

        // The "epsilon" parameter of the backbone
        unsigned int localSearchDistance;
        // The backbone vertices
        unordered_set<VertexID> backboneVertices;
        dynamic_bitset<> isBackboneVertex;
        // The backbone itself
        // unique_ptr<const DGraph> backbone;
        unique_ptr<DGraph> backbone;

        // Indexing
        void selectBackboneVertices();
        void createBackboneEdges();
        void indexBackbone();
        void buildIndex();

        // General querying
        inline bool bfsLocally(VertexID source, VertexID target, LabelSet ls);
        inline bool queryBackbone(VertexID source, VertexID target, LabelSet ls);
        inline bool computeQuery(VertexID source, VertexID target, LabelSet ls);
        inline void markTargetsForBackboneBfs(VertexID target, LabelSet ls);
        dynamic_bitset<> bfsBackboneTargets;

        // -- Indexing method specific --
        // BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA
        void localMeetingCriteriaSetCover();
        pair<
            backbonens::LabelledDistancedReachabilityMap,
            unordered_map<VertexID, backbonens::LabelledDistancedReachabilityMap>> generateGroundSetAndCandidates();

        // BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_*
        void oneSideConditionCover();
        // BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER
        inline vector<VertexID> getVerticesInRandomOrder();
        // BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER
        inline vector<VertexID> getVerticesInDegreeOrder();


        // BackboneIndexingMethod::BFS
        // BackboneIndexingMethod::TRANSITIVE_CLOSURE
        backbonens::LabelledDistancedReachabilityMap backboneTransitiveClosure;
        // BackboneIndexingMethod::LANDMARK_*
        unique_ptr<LandmarkedIndex> backboneLi;
        // BackboneIndexingMethod::TWOHOP
        unique_ptr<TwoHopIndex> backboneTwoHopIndex;

        // LocalSearchMethod::UNIDIRECTIONAL_BFS
        inline bool uniDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls);
        // LocalSearchMethod::BIDIRECTIONAL_BFS
        bool biDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls);

        // -- Misc --
        // Speedup reachable backbone vertices discovery
        void cacheVertexToBackboneReachability();
        inline queue<VertexID> accessBackboneOutQueue(VertexID source, LabelSet ls);
        inline vector<VertexID> accessBackboneOut(VertexID source, LabelSet ls);
        TuplesList backboneReachableOut;
        inline unordered_set<VertexID> accessBackboneInSet(VertexID target, LabelSet ls);
        inline vector<VertexID> accessBackboneIn(VertexID target, LabelSet ls);
        TuplesList backboneReachableIn;
};
#endif
