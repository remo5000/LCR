#include "Graph.h"

#include <set>
#include <stack>
#include <boost/dynamic_bitset.hpp>

using namespace graphns;
using namespace boost;

#ifndef DGRAPH_H
#define DGRAPH_H

/**
* This class contains a normal graph, just two lists with In and Out edges
* per node. It is a directed graph.
*/
class DGraph : public Graph
{
    public:
        DGraph(string fileName);
        DGraph(EdgeSet* edgeSet);
        DGraph(EdgeSet* edgeSet, int pN, int pL);
        DGraph(EdgeSet* edgeSet, int pN, int pL, bool allowMultipleEdges);
        ~DGraph() = default;
        void construct(EdgeSet* edgeSet, int pN, int pL, bool allowMultipleEdges);

        void buildGraph(graphns::EdgeSet* edges);
        EdgeSet* loadEdgeFile(string fileName);
        void loadEdgeStats(graphns::EdgeSet* edgeSet);

        int getGraphSizeInBytes();
        double getGraphConstructionTime();

        const SmallEdgeSet& getOutNeighbours(graphns::VertexID w);
        void getOutNeighbours(graphns::VertexID w, SmallEdgeSet& outNeighbours);
        const SmallEdgeSet& getInNeighbours(graphns::VertexID w);
        void getInNeighbours(graphns::VertexID w, SmallEdgeSet& inNeighbours);
        void getAllNeighbours(graphns::VertexID w, SmallEdgeSet& allNeighbours);

        // prints stats of the graph
        std::string getStats();
        std::string toString();

        int getNumberOfVertices();
        int getNumberOfLabels();
        int getNumberOfEdges() const;

        void addNode();
        void removeNode(graphns::VertexID w);
        void addEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel);
        void addMultiEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelSet newLabelSet);
        void removeEdge(graphns::VertexID v, graphns::VertexID w);
        void changeLabel(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel);
        bool hasMultiEdge(graphns::VertexID v , graphns::VertexID w, graphns::LabelSet ls);

        bool findInsertablePosition(graphns::VertexID w, graphns::SmallEdgeSet& ses, int& pos);
        void insertEdge(graphns::VertexID v, graphns::VertexID w, graphns::LabelID newLabel, SmallEdgeSet& ses);
        LabelID getLabelID(graphns::VertexID v , graphns::VertexID w);
        long getCountPerLabel(graphns::LabelID l);

        bool hasEdge(VertexID v, VertexID w);

        double computeAverageDiameter();
        int computerNumberOfTriangles();
        int computeNumberOfConnectedTriplets();
        double computeClusterCoefficient();
        void tarjan(vector< vector<VertexID> >& SCCs);
        void tarjanStrongConnect(int v, int& index, stack<VertexID>& q, vector< int >& indexPerNode,
            vector< int >& lowlinkPerNode, vector< bool >& onStack, vector< vector<VertexID> >& SCCs);
        int computeDiameter();

    private:
        SmallEdgeSets outE, inE;
        double constStartTime, constEndTime;
        bool allowMultipleEdges;

};
#endif
