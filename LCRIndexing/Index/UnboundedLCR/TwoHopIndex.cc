#include "TwoHopIndex.h"

#include <iostream>
#include <queue>
#include <set>
#include <unordered_set>

using namespace std;
using namespace indexns;

TwoHopIndex::TwoHopIndex(Graph* mg)
{
    constStartTime = getCurrentTimeInMilliSec();
    graph = mg;
    N = graph->getNumberOfVertices();
    indexType = IndexType::TwoHop;

    inIndex.resize(N);
    outIndex.resize(N);

    queryStart = 0.0;
    queryEndTime = 0.0;
    visitedSetSize = 0;
    indexDirection = BOTHINDEX;

    buildIndex();

    constEndTime = getCurrentTimeInMilliSec();
    totalConstTime = constEndTime - constStartTime;
    this->didComplete = true;
}

unsigned long TwoHopIndex::getIndexSizeInBytes()
{
    return getIndexSizeInBytesM();
};

bool TwoHopIndex::computeQuery(VertexID source, VertexID target, LabelSet ls)
{
    if (source == target) return true;
    if (ls == 0) return true;

    for (const Tuple& tuple : this->outIndex[source]) {
        const VertexID& outVertex = tuple.first;

        if (!labelSetInLabelSets(ls, tuple.second))
            continue;

        if (outVertex == target)
            return true;

        int pos;
        if (findTupleInTuples(outVertex, this->inIndex[target], pos)) {
            Tuple& tuple = this->inIndex[target][pos];
            if (labelSetInLabelSets(ls, tuple.second))
                return true;
        }
    }

    int pos;
    return findTupleInTuples(source, this->inIndex[target], pos)
        && labelSetInLabelSets(ls, this->inIndex[target][pos].second);
};

bool TwoHopIndex::query(VertexID source, VertexID target, LabelSet ls)
{
    queryStart = getCurrentTimeInMilliSec();
    bool b = this->computeQuery(source, target, ls);
    queryEndTime = getCurrentTimeInMilliSec();
    return b;
};

void TwoHopIndex::queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach)
{
};

void TwoHopIndex::buildIndex()
{
    int quorum = N/100;
    for (VertexID source = 0; source < this->N; source++) {
        if (source % quorum == 0) {
            double perc = source;
            perc /= N;
            perc *= 100;
            cout << "TwoHopIndex::buildIndex(): " << perc << "% done..." << endl;
        }

        for (int outwards = 0; outwards < 2; outwards++) {
            auto getNeighbors = [&](VertexID v) {
                return outwards 
                    ? this->graph->getOutNeighbours(v) 
                    : this->graph->getInNeighbours(v);
            };
            TuplesList& tuplesList =  outwards 
                ? this->inIndex
                : this->outIndex;

            queue<pair<VertexID, LabelSet>> q;
            q.push(make_pair(source, 0));
            while (q.empty() == false) {
                VertexID vertex;
                LabelSet ls;
                std::tie(vertex, ls) = q.front();
                q.pop();

                for (const auto& smallEdge : getNeighbors(vertex)) {
                    const VertexID& neighbor = smallEdge.first;
                    const LabelSet& ls2 = smallEdge.second;
                    const LabelSet& newLs = joinLabelSets(ls, ls2);

                    // Vertex has already been indexed (Rule 1)
                    if (neighbor < source) 
                        continue;

                    // Vertex already has a dominating entry for 'source'
                    Tuples& tuples = tuplesList[neighbor];
                    int pos;
                    if (!findTupleInTuples(source, tuples, pos)) {
                        indexns::LabelSets lss = indexns::LabelSets();
                        lss.reserve( this->graph->getNumberOfLabels() * 2 );
                        indexns::Tuple newTuple = make_pair(source, lss);
                        tuples.insert( tuples.begin() + pos, newTuple );
                    }
                    Tuple& tuple = tuples[pos];
                    LabelSets& lss = tuple.second;
                    if (!tryInsertLabelSet(newLs, lss)) {
                        continue;
                    }

                    q.push(make_pair(neighbor, newLs));
                }
            }
        }
    }
};
