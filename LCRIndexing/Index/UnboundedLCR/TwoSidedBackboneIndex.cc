#include "TwoSidedBackboneIndex.h"

using namespace twosidedbackbonens;
using namespace indexns;


TwoSidedBackboneIndex::TwoSidedBackboneIndex(Graph* mg) {
    this->graph = mg;
    this->indexType = IndexType::TwoSidedBackbone;
    int N = graph->getNumberOfVertices();
    int L = graph->getNumberOfLabels();

    this->indexDirection = BOTHINDEX;
    this->isBlockedMode = false; // TODO this has something to do with input

    // Construct index
    this->didComplete = false;
}

unsigned long TwoSidedBackboneIndex::getIndexSizeInBytes()
{
    return getIndexSizeInBytesM();
};

bool TwoSidedBackboneIndex::query(VertexID source, VertexID target, LabelSet ls)
{
    cout << "TwoSidedBackboneIndex::query source=" << to_string(source) << ",target=" << to_string(target) << ",ls=" << labelSetToString(ls) << endl;
    queryStart = getCurrentTimeInMilliSec();
    bool b = true; // TODO
    queryEndTime = getCurrentTimeInMilliSec();
    cout << "TwoSidedBackboneIndex::query answer =" << b << endl;
    return b;
}

// For some reason this is blank for a lot of methods -- we shall leave it blank for now.
void TwoSidedBackboneIndex::queryAll(VertexID source, LabelSet ls, dynamic_bitset<>& canReach)
{

};

void TwoSidedBackboneIndex::buildIndex()
{

};
