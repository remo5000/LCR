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

bool TwoHopIndex::computeQueryBackbone(
		const vector<VertexID>& sources,
		const vector<VertexID>& targets,
		const dynamic_bitset<>& isInTargets,
		const LabelSet& ls)
{
    for (const VertexID& source : sources) {
        for (const Tuple& tuple : this->outIndex[source]) {
            const VertexID& outVertex = tuple.first;

            if (!labelSetInLabelSets(ls, tuple.second))
                continue;

            if (isInTargets[outVertex])
                return true;

            for (const VertexID& target : targets) {
                int pos;
                if (findTupleInTuples(outVertex, this->inIndex[target], pos)) {
                    Tuple& tuple = this->inIndex[target][pos];
                    if (labelSetInLabelSets(ls, tuple.second))
                        return true;
                }
            }

        }
    }

    return false;
};

bool TwoHopIndex::queryBackbone(
		const vector<VertexID>& sources,
		const vector<VertexID>& targets,
		const dynamic_bitset<>& isInTargets,
		const LabelSet& ls)
{
    queryStart = getCurrentTimeInMilliSec();
    bool b = this->computeQueryBackbone(sources, targets, isInTargets, ls);
    queryEndTime = getCurrentTimeInMilliSec();
    return b;
};

bool TwoHopIndex::computeQuery(VertexID source, VertexID target, LabelSet ls)
{
    if (source == target) return true;
    if (ls == 0) return false;

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

    return false;
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
    // Use degree ordering for vertices
    struct sort_pred
    {
        bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
        {
            return left.second > right.second;
        }
    };
    vector< pair< VertexID, int > > degreePerNode;
    for (VertexID i = 0; i < N; i++)
        degreePerNode.push_back(
            make_pair(
                i,
                // Use product of inDeg and outDeg
                this->graph->getOutNeighbours(i).size() *  this->graph->getInNeighbours(i).size()
            )
        );
    sort(degreePerNode.begin(), degreePerNode.end(), sort_pred());

    // Keep track of already-processed vertices
    dynamic_bitset<> done = dynamic_bitset<>(N);

    for (int i = 0; i < degreePerNode.size(); i++) {
    // for (const auto& p : degreePerNode) {
        const VertexID& source = degreePerNode[i].first;

        double perc = i;
        perc /= N;
        perc *= 100;
        cout << "TwoHopIndex::buildIndex(): " << perc << "% done..." << endl;

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


                // Vertex has already been indexed (Rule 1)
                if (done[vertex])
                    continue;

                for (const auto& smallEdge : getNeighbors(vertex)) {
                    const VertexID& neighbor = smallEdge.first;
                    const LabelSet& ls2 = smallEdge.second;
                    const LabelSet& newLs = joinLabelSets(ls, ls2);

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

        done[source] = 1;
    }
};
