#include "../../Graph/Graph.h"
#include "Index.h"

using namespace indexns;
using namespace graphns;

#ifndef TWOHOPINDEX_H
#define TWOHOPINDEX_H

/**
P2H Index, adapted from the paper, "Answering Billion-Scale Label-Constrained Reachability Queries within Microsecond".

ACM Ref:
You Peng, Ying Zhang, Xuemin Lin, Lu Qin, and Wenjie Zhang. 2020. Answering billion-scale label-constrained reachability queries within microsecond. Proc. VLDB Endow. 13, 6 (February 2020), 812–825. DOI:https://doi.org/10.14778/3380750.3380753
*/
class TwoHopIndex : public Index
{
    public:

	TwoHopIndex(Graph* mg);
	~TwoHopIndex() = default;

	bool query(VertexID source, VertexID target, LabelSet ls);
	bool queryBackbone(
		const vector<VertexID>& sources,
		const vector<VertexID>& targets,
		const dynamic_bitset<>& isInTargets,
		const LabelSet& ls);
	void queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach);
	unsigned long getIndexSizeInBytes();

    private:
	int visitedSetSize;
	void buildIndex();
	bool computeQuery(VertexID source, VertexID target, LabelSet ls);
	bool computeQueryBackbone(
		const vector<VertexID>& sources,
		const vector<VertexID>& targets,
		const dynamic_bitset<>& isInTargets,
		const LabelSet& ls);

	indexns::TuplesList inIndex;
	indexns::TuplesList outIndex;
	unsigned int N;
};
#endif
