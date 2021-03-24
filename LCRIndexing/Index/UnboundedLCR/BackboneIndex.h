#include "../../Graph/Graph.h"
#include "../../Graph/DGraph.h"
#include "Index.h"

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
    class LabelledDistancedReachabilityMap {
        public:
            void insert(VertexID source, VertexID destination, LabelSet ls, unsigned int distance) {
                // Check if this LS is strictly dominated by any other LS
                for (const auto& p : this->m[source][destination]) {
                    LabelSet foundLs; unsigned int foundDist;
                    std::tie(foundLs, foundDist) = p;

                    if (isLabelSubset(foundLs, ls) && foundDist <= distance)
                        return;
                }

                // Remove all labelsets that it strictly dominates in terms of labels and distance;
                vector<LabelSet> toRemove;
                for (const auto& p : this->m[source][destination]) {
                    LabelSet foundLs; unsigned int foundDist;
                    std::tie(foundLs, foundDist) = p;

                    if (isLabelSubset(ls, foundLs) && distance <= foundDist) {
                        toRemove.push_back(foundLs);
                    }
                }
                for (LabelSet lsToRemove : toRemove)
                    m[source][destination].erase(lsToRemove);

                this->m[source][destination][ls] = distance;
            }
            bool isPresent(VertexID source, VertexID destination, LabelSet ls) {
                if (!m.count(source)) return false;
                if (!m[source].count(destination)) return false;
                for (auto p : m[source][destination]) {
                    if (isLabelSubset(p.first, ls)) {
                        return true;
                    }
                }
                return false;
            }
            void erase(VertexID source, VertexID destination, LabelSet ls) {
                m.at(source).at(destination).erase(ls);
                if (m.at(source).at(destination).size() == 0) m.at(source).erase(destination);
                if (m.at(source).size() == 0) m.erase(source);
            }
            void erase(VertexID source, VertexID destination) {
                if (!m.count(source)) return;
                if (!m[source].count(destination)) return;
                m.at(source).erase(destination);
                if (m.at(source).size() == 0) m.erase(source);
            }
            unsigned int getDistance(VertexID source, VertexID destination, LabelSet ls) {
                if (this->isPresent(source, destination, ls)) {
                    return this->m[source][destination][ls];
                }


                for (auto p : m[source][destination]) {
                    if (isLabelSubset(p.first, ls)) {
                        return p.second;
                    }
                }

                return std::numeric_limits<unsigned int>::max();
            }
            unsigned int size() {
                unsigned int ans = 0;
                for (const auto& p1 : m)
                    for (const auto& p2 : p1.second)
                        ans += p2.second.size();
                return ans;
            }
            string toString() {
                stringstream s;
                s << "[ LabelledDistancedReachabilityMap (" << this->size() << " item(s)\n";
                for (const auto& p1 : m) {
                    for (const auto& p2 : p1.second) {
                        for (const auto& p3 : p2.second) {
                            s << "  ("  << p1.first;
                            s << "->" << p2.first << ", ";
                            s << labelSetToLetters(p3.first) << ", ";
                            s << p3.second << ")\n";
                        }
                    }
                }
                s << "]\n";
                return s.str();
            }
            map<VertexID, SmallEdgeSet> toEdgeMap() const {
                map<VertexID, SmallEdgeSet> result;
                for (const auto& p1 : m) {
                    for (const auto& p2 : p1.second) {
                        for (const auto& p3 : p2.second) {
                            VertexID u = p1.first;
                            VertexID v = p2.first;
                            LabelSet ls = p3.first;
                            result[u].emplace_back(v, ls);
                        }
                    }
                }
                return result;
            }
            typedef pair<VertexID, std::pair<VertexID, LabelSet>> Item;
            // TODO change to triples
            set<Item> toTuples() const {
                set<Item> result;
                for (const auto& p1 : m) {
                    for (const auto& p2 : p1.second) {
                        for (const auto& p3 : p2.second) {
                            VertexID u = p1.first;
                            VertexID v = p2.first;
                            LabelSet ls = p3.first;
                            result.insert({u, {v, ls}});
                        }
                    }
                }
                return result;
            }
            unsigned long getSizeInBytes() {
                unsigned long size = sizeof(m);
                for (const auto& p1 : m) {
                    size += sizeof(p1.second);
                    for (const auto& p2 : p1.second) {
                        size += sizeof(p2.second);
                        for (const auto& p3 : p2.second) {
                            size += sizeof(p3);
                        }
                    }
                }
                return size;
            }
        private:
            // source -> dest -> ls -> dist
            map<VertexID, map<VertexID, map<LabelSet, unsigned int> > > m;
    };
}

enum class BackboneVertexSelectionMethod {
    LOCAL_MEETING_CRITERIA
};

enum class BackboneEdgeCreationMethod {
    BFS
};

enum class BackboneIndexingMethod {
    BFS,
    TRANSITIVE_CLOSURE
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
                BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA,
                BackboneEdgeCreationMethod::BFS,
                BackboneIndexingMethod::TRANSITIVE_CLOSURE,
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

        // BackboneIndexingMethod::BFS
        bool bfsBackbone(VertexID source, VertexID target, LabelSet ls);
        // BackboneIndexingMethod::TRANSITIVE_CLOSURE
        bool backboneQueryTransitiveClosure(VertexID source, VertexID target, LabelSet ls);

        // LocalSearchMethod::UNIDIRECTIONAL_BFS
        bool uniDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls);
        // LocalSearchMethod::BIDIRECTIONAL_BFS
        bool biDirectionalLocalBfs(VertexID source, VertexID target, LabelSet ls);

        // -- Misc --
        // Speedup reachable backbone vertices discovery
        void cacheVertexToBackboneReachability();
        map<VertexID, SmallEdgeSet> backboneReachableOut;
        map<VertexID, SmallEdgeSet> backboneReachableIn;


        backbonens::LabelledDistancedReachabilityMap backboneTransitiveClosure;
};
#endif
