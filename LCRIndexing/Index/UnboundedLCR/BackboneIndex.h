#include "../../Graph/Graph.h"
#include "../../Graph/DGraph.h"
#include "../../Index/UnboundedLCR/LandmarkedIndex.cc"
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
    typedef pair<VertexID, LabelSetAndDistanceList> DTuple;
    typedef vector<DTuple> DTuples;
    typedef vector<DTuples> DTuplesList;

    struct comp_distance_tuples
    {
        bool operator()( const DTuple& lhs,  const int& rhs ) const
        {
            return lhs.first < rhs;
        }

        bool operator()( const int& lhs, const DTuple& rhs ) const
        {
            return lhs<rhs.first;
        }
    };

    class LabelledDistancedReachabilityMap {
        public:
            LabelledDistancedReachabilityMap(VertexID N, VertexID M) {
                this->tuplesList.resize(N);
                for (auto& dtls : tuplesList)
                    for (VertexID i = 0; i < M; i++)
                        dtls.push_back(make_pair(i, LabelSetAndDistanceList()));
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

                if (!isPresent(destination, tuples, pos)) {
                    DTuple newTuple = {destination, {}};
                    tuples.insert(tuples.begin()+pos, newTuple);
                }

                DTuple& tuple = tuples[pos];
                LabelSetAndDistanceList& lsdl = tuple.second;

                if(!findLabelSetAndDistanceInList(ls, lsdl, pos)) {
                    lsdl.push_back(make_pair(ls, distance));
                    _size++;
                }
            }
            inline bool isPresent(VertexID source, VertexID destination, LabelSet ls) {
                ensureTuplesListHasSource(source);
                assert(source < tuplesList.size());
                const DTuples& tuples = tuplesList[source];

                auto it = std::lower_bound(tuples.begin(), tuples.end(), destination, comp_distance_tuples());
                auto pos =  it - tuples.begin();

                if( pos < 0 || pos >= tuples.size())
                {
                    return false;
                }

                const DTuple& dt = tuples[pos];
                if(dt.first != destination)
                {
                    return false;
                }

                const LabelSetAndDistanceList& lsdl = dt.second;

                int pos2 = -1;
                return findLabelSetAndDistanceInList(ls, lsdl, pos2);
            }
            inline void erase(VertexID source, VertexID destination, LabelSet ls) {
                ensureTuplesListHasSource(source);
                int pos;
                DTuples& tuples = tuplesList[source];
                if (!isPresent(destination, tuples, pos)) {
                    return;
                }
                DTuple& tuple = tuples[pos];
                LabelSetAndDistanceList& lsdl = tuple.second;

                // While source-ls->destination, erase those entries.
                while(findLabelSetAndDistanceInList(ls, lsdl, pos)) {
                    lsdl.erase(lsdl.begin() + pos);
                    _size--;
                }
            }
            inline void erase(VertexID source, VertexID destination) {
                ensureTuplesListHasSource(source);
                int pos;
                DTuples& tuples = tuplesList[source];
                if (!isPresent(destination, tuples, pos)) {
                    return;
                }
                _size -= tuples[pos].second.size();
                tuples.erase(tuples.begin() + pos);
            }
            inline unsigned int getDistance(VertexID source, VertexID destination, LabelSet ls) {
                ensureTuplesListHasSource(source);
                int pos;
                DTuples& tuples = tuplesList[source];
                if (isPresent(destination, tuples, pos)) {
                    DTuple& tuple = tuples[pos];
                    LabelSetAndDistanceList& lsdl = tuple.second;

                    if(findLabelSetAndDistanceInList(ls, lsdl, pos)) {
                        return lsdl[pos].second;
                    }
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
                    for (const DTuple& dt : dtls) {
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
                    for (const DTuple& dt : dtls) {
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
                    for (const DTuple& dt : dtls) {
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
                    for (const DTuple& dt : dtls) {
                        VertexID destination = dt.first;
                        // Don't add pair directly, as sizeof impl can vary.
                        size += sizeof(VertexID);
                        size += dt.second.size() * sizeof(LabelSetAndDistance);
                    }
                }
                return size;
            }
        private:
            // source -> dest -> ls -> dist
            // map<VertexID, map<VertexID, map<LabelSet, unsigned int> > > m;
            DTuplesList tuplesList;

            inline bool isPresent(VertexID destination, const DTuples& tuples, int& pos) {
                auto it = std::lower_bound(tuples.begin(), tuples.end(), destination, comp_distance_tuples());
                pos =  it - tuples.begin();
                return pos >= 0 && pos < tuples.size() && tuples[pos].first == destination;
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
    ONE_SIDE_CONDITION
};

enum class BackboneEdgeCreationMethod {
    BFS
};

enum class BackboneIndexingMethod {
    BFS,
    TRANSITIVE_CLOSURE,
    LANDMARK_NO_EXTENSIONS,
    LANDMARK_FULL,
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
                BackboneVertexSelectionMethod::ONE_SIDE_CONDITION,
                BackboneEdgeCreationMethod::BFS,
                BackboneIndexingMethod::BFS,
                LocalSearchMethod::UNIDIRECTIONAL_BFS
            ) {};

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
        // The backbone itself
        // unique_ptr<const DGraph> backbone;
        unique_ptr<DGraph> backbone;

        // Indexing
        void selectBackboneVertices();
        void createBackboneEdges();
        void indexBackbone();
        void buildIndex();

        // General querying
        bool bfsLocally(VertexID source, VertexID target, LabelSet ls);
        bool queryBackbone(VertexID source, VertexID target, LabelSet ls);
        bool computeQuery(VertexID source, VertexID target, LabelSet ls);



        // -- Indexing method specific --
        // BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA
        void localMeetingCriteriaSetCover();
        pair<
            backbonens::LabelledDistancedReachabilityMap,
            unordered_map<VertexID, backbonens::LabelledDistancedReachabilityMap>> generateGroundSetAndCandidates();

        // BackboneVertexSelectionMethod::ONE_SIDE_CONDITION
        void oneSideConditionCover();

        // BackboneIndexingMethod::BFS
        bool bfsBackbone(VertexID source, VertexID target, LabelSet ls);

        // BackboneIndexingMethod::TRANSITIVE_CLOSURE
        backbonens::LabelledDistancedReachabilityMap backboneTransitiveClosure;
        bool backboneQueryTransitiveClosure(VertexID source, VertexID target, LabelSet ls);

        // BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS
        unique_ptr<LandmarkedIndex> backboneLi;
        bool backboneQueryLandmarks(VertexID source, VertexID target, LabelSet ls);

        // LocalSearchMethod::UNIDIRECTIONAL_BFS
        bool uniDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls);
        // LocalSearchMethod::BIDIRECTIONAL_BFS
        bool biDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls);

        // -- Misc --
        // Speedup reachable backbone vertices discovery
        void cacheVertexToBackboneReachability();
        map<VertexID, SmallEdgeSet> backboneReachableOut;
        map<VertexID, SmallEdgeSet> backboneReachableIn;
};
#endif
