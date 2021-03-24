#include "../../Index/UnboundedLCR/BackboneIndex.cc"
#include "../../Graph/DGraph.cc"

#include <string>

using namespace std;
using namespace graphns;
using namespace indexns;

void test(const bool a, int& t)
{
    if( a == true )
    {
        cout << "testID=" << to_string(t) << " passed" << endl;
    }
    else
    {
        cout << "testID=" << to_string(t) << " failed" << endl;
    }

    t++;
}

void runQuery(Index* index, VertexID v, VertexID w, LabelSet ls, bool answer, int& t)
{
    test(index->query(v,w,ls) == answer, t);
    cout << "query(" << v << "," << w << "," << ls << ") took time(s)=" << index->getLastQueryTime() << endl;
}

void testBackboneVertices(
    const string& edgeFileName,
    const unsigned int& localDist,
    const unordered_set<VertexID>& expectedBackboneVertices,
    int& testId
)
{
    Graph* g1 = new DGraph(edgeFileName);
    const BackboneIndex index(g1, localDist);
    const unordered_set<VertexID>& actual = index.getBackBoneVertices();
    for (VertexID v : actual) cout << "vertex " << v << " is present in backbone\n";
    test(actual == expectedBackboneVertices, testId);
    free(g1);
}

int main(int argc, char *argv[])
{
    int t = 1;
    unordered_set<VertexID> expectedBackboneVertices = {};
    // -- Test for no possible backbone vertices (all reachability is local) --
    testBackboneVertices("tests/graphs/BackboneEmpty.edge", 2, expectedBackboneVertices, t);
    testBackboneVertices("tests/graphs/BackboneComplete4.edge", 2, expectedBackboneVertices, t);
    testBackboneVertices("tests/graphs/BackboneLine1.edge", 2, expectedBackboneVertices, t);

    // -- Test for >0 backbone vertices --
    expectedBackboneVertices = {1};
    testBackboneVertices("tests/graphs/BackboneBranch1.edge", 2, expectedBackboneVertices, t);
    expectedBackboneVertices = {2};
    testBackboneVertices("tests/graphs/BackboneLine2.edge", 2, expectedBackboneVertices, t);
    // Regression -- either 0 or 2 works.
    expectedBackboneVertices = {0};
    testBackboneVertices("tests/graphs/BackboneCycle1.edge", 2, expectedBackboneVertices, t);


    return 0;
}
