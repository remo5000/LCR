#include <cstdlib>
#include <sstream>
#include <vector>
#include <random>
#include <queue>
#include <algorithm>
#include <string>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <unistd.h>
#include <sys/time.h>

#include "Graph/DGraph.h"
#include "Index/UnboundedLCR/Index.h"
#include "Index/UnboundedLCR/BFSIndex.h"
#include "Index/UnboundedLCR/BackboneIndex.h"
#include "Index/UnboundedLCR/LandmarkedIndex.h"
#include "Index/UnboundedLCR/TwoHopIndex.h"
#include "Index/UnboundedLCR/Zou.h"

using namespace std;
using namespace indexns;
using namespace graphns;

#ifdef _WIN32
#include <windows.h>
#define SYSERROR()  GetLastError()
#else
#include <errno.h>
#define SYSERROR()  errno
#endif

bool is_file_exist(string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

void loadQueryFile(string queryFileName, set<Query>& tmpSet)
{
    //cout << "queryFileName=" << queryFileName << endl;
    tmpSet.clear();
    string line = "";
    ifstream queryFile (queryFileName);

    VertexID f,t;
    LabelSet l2;

    if (queryFile.is_open())
    {
        while ( getline (queryFile,line) )
        {
            //cout << "line="<< line << endl;
            LabelID lID = 0;
            istringstream iss(line);
            string from, label, to;
            iss >> from >> to >> label;
            istringstream (from) >> f;
            istringstream (to) >> t;
            istringstream (label) >> lID;
            l2 = lID;

            Query q = make_pair( make_pair(f,t), l2);
            tmpSet.insert(q);
        }

        queryFile.close();
    }
}

/*
Runs all queries obtained by loading a query set for a given method. Index construction
time, index size and query times are logged.

Each query is executed 10 times and the performance is averaged.
*/
int runTestsPerIndex(Index* index, vector< vector< vector<double> > >& queryTimes,
        vector< double >& queryTimeSums,
        vector<unsigned long>& indexSizes, vector<double>& indexTimes, unsigned int noOfQuerySets,
        string edgeFile, unsigned long size, double indexingTime)
{
    string indexName = index->getIndexTypeAsString();

    // store the size and index time
    indexSizes.push_back( size );
    indexTimes.push_back( indexingTime );

    vector< vector< double > > s (0);

    // loop over all query sets
    for(int j = 0; j < noOfQuerySets; j++)
    {
        cout << "Method: " << indexName << " queryset: " << j << "\n";

        vector< double > t (0);
        double trueSum = 0.0;
        double falseSum = 0.0;

        string trueFileName = edgeFile + to_string(j) + ".true";
        string falseFileName = edgeFile + to_string(j) + ".false";

        set<Query> trueSet; //query set that returns true
        set<Query> falseSet; //qyery set that returns false

        loadQueryFile(trueFileName, trueSet); // load the two query files
        loadQueryFile(falseFileName, falseSet);

        int query_times = trueSet.size();
        int i = 0;

	int NUMBER_OF_RUNS = 5;
        for (auto p : trueSet)
        {
            VertexID from = p.first.first;
            VertexID to = p.first.second;
            LabelSet ls = p.second;

            // average over 5 runs
            double avg = 0.0;
            vector< double > avgs;

            for(int k = 0; k < NUMBER_OF_RUNS; k++)
            {
                if( index->query(from, to, ls) == false )
                {
                    cout << "Query " << i << ": (" << from << "," << to << "," << labelSetToString(ls) << ") should be true" << endl;
                    avg = -1.0;
                    return 1;
                }

                avgs.push_back( index->getLastQueryTime() );
            }

            sort( avgs.begin(), avgs.end() );
            for(int k = 0; k < NUMBER_OF_RUNS; k++)
            {
                avg += avgs[k];
            }
            avg /= NUMBER_OF_RUNS;

            cout << "*Query " << i << ": (" << from << "," << to << "," << labelSetToString(ls) << ", true), avg=" << print_digits(avg,11) << endl;
            i++;

            trueSum += avg;
            t.push_back( avg );
        }

        i = 0;
        for (auto p : falseSet)
        {
            VertexID from = p.first.first;
            VertexID to = p.first.second;
            LabelSet ls = p.second;

            // average over 5 runs
            double avg = 0.0;
            vector< double > avgs;

            for(int k = 0; k < NUMBER_OF_RUNS; k++)
            {
                if( index->query(from, to, ls) == true )
                {
                    cout << "Query " << i << ": (" << from << "," << to << "," << labelSetToString(ls) << ") should be false" << endl;
                    avg = -1.0;
                    return 1;
                }

                avgs.push_back( index->getLastQueryTime() );
            }

            sort( avgs.begin(), avgs.end() );
            for(int k = 0; k < NUMBER_OF_RUNS; k++)
            {
                avg += avgs[k];
            }

            avg /= NUMBER_OF_RUNS;
            cout << "*Query " << i << ": (" << from << "," << to << "," << labelSetToString(ls) << ", false), avg=" << print_digits(avg,11) << endl;
            i++;

            falseSum += avg;
            t.push_back( avg );
        }

        cout << "indexName=" << indexName << ",j=" << j << ",trueSum=" << trueSum << ",falseSum=" << falseSum << endl;

        queryTimeSums.push_back( trueSum );
        queryTimeSums.push_back( falseSum );

        s.push_back( t );
    }

    queryTimes.push_back( s );

    cout << "------" << endl;
    return 0;
}

int main(int argc, char *argv[]) {
    if(argc < 4)
    {
        cout << "usage: ./runExperiment <edge_file> <number of query sets> <output_file>" << endl;
        cout << "or: ./runExperiment <edge_file> <number of query sets> <output_file> <number of methods (BFS|Partial|Landmarked|DoubleBFS)>" << endl;
        exit(1);
    }

    string edge_file = argv[1]; // the file containing the graph
    int noOfQuerySets = atoi(argv[2]);
    string output_file = argv[3]; // the file where the output goes

    set<Query> tmpSet;
    loadQueryFile(edge_file + "0.true", tmpSet); // load the two query files
    int noOfQueries = tmpSet.size();

    int noOfMethods = 6; // the total number of methods
    if( argc == 5 )
    {
        noOfMethods = atoi(argv[4]);
    }
    int firstMethod = 0; // the method to start with: BFS, ExactIndex, ExactHopIndex, Joindex

    cout << print_digits(10000000.0, 0) << endl;
    cout << "number of methods: " << noOfMethods << endl;
    cout << "number of query sets: " << noOfQuerySets << endl;
    cout << "number of queries in each true or false: " << noOfQueries << endl;

    unique_ptr<DGraph> graph = unique_ptr<DGraph>(new DGraph(edge_file));
    int labelSetSize = graph->getNumberOfLabels();

    cout << "dataset: " << edge_file << " with |L|=" << labelSetSize << " and |V|=" << graph->getNumberOfVertices() << endl;

    // variables that store the output
    vector<unsigned long> indexSizes(0);
    vector<double> indexTimes(0);
    vector< vector< vector<double> > > queryTimes ( 0 ); // queryTimes[i][j][k] gives the time for method i, queryset j,
    // and query k where k >= 0 and k < 2*noOfQueries
    vector< string > methodNames(0);

    vector< double > queryTimeSums( 0 );

    int N = graph->getNumberOfVertices();
    int M = graph->getNumberOfEdges();

    // TODO formalize
    map<string, double> nameToBackboneVertexRatio;
    map<string, double> nameToBackboneEdgeRatio;


    // Here we loop over all methods
    for(int i = firstMethod; i < noOfMethods; i++)
    {
        cout << "method i=" << i << endl;

        Index* index;

        long altSize = -1;

        if( i == __COUNTER__ )
            index = new BFSIndex(graph.get());

        // P2H
        if( i == __COUNTER__ )
            index = new TwoHopIndex(graph.get());

        // // LI (no extensions)
        // if( i == __COUNTER__ )
	//       {
        //     int k = 1250 + sqrt(N);
        //     if (k >= N) k = sqrt(N);
        //     int b = 20;
        //     index = new LandmarkedIndex(graph.get(), false, false, k, b);
	//       }


        // LI+ (both extensions)
        if( i == __COUNTER__ )
	      {
            int k = 1250 + sqrt(N);
            if (k >= N) k = sqrt(N);
            int b = 20;
            index = new LandmarkedIndex(graph.get(), true, true, k, b);
	      }

        float logn = log2(N);
        float loglogn = log2(log2(N));
        // unsigned int localDist = max(2, (int)loglogn);
	// We switch to 2 as it is the only viable option for 
	// dense, large graphs.
        unsigned int localDist = 2;

        // // Backbone(LMC, BFS)
        // if (i == __COUNTER__) {
        //     index = new BackboneIndex(
        //             graph.get(),
        //             localDist,
        //             BackboneVertexSelectionMethod::LOCAL_MEETING_CRITERIA,
        //             BackboneEdgeCreationMethod::BFS,
        //             BackboneIndexingMethod::BFS,
        //             LocalSearchMethod::UNIDIRECTIONAL_BFS
        //     );
        // }
        // Backbone(1SCD, BFS)
        if (i == __COUNTER__) {
            index = new BackboneIndex(
                    graph.get(),
                    localDist,
                    BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER,
                    BackboneEdgeCreationMethod::BFS,
                    BackboneIndexingMethod::BFS,
                    LocalSearchMethod::UNIDIRECTIONAL_BFS
            );
        }
        // // Backbone(1SCR, BFS)
        // if (i == __COUNTER__) {
        //     index = new BackboneIndex(
        //             graph.get(),
        //             localDist,
        //             BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_RANDOM_ORDER,
        //             BackboneEdgeCreationMethod::BFS,
        //             BackboneIndexingMethod::BFS,
        //             LocalSearchMethod::UNIDIRECTIONAL_BFS
        //     );
        // }

        // // Backbone(1SCD, P2H)
        // if (i == __COUNTER__) {
        //     index = new BackboneIndex(
        //             graph.get(),
        //             localDist,
        //             BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER,
        //             BackboneEdgeCreationMethod::BFS,
        //             BackboneIndexingMethod::TWOHOP,
        //             LocalSearchMethod::UNIDIRECTIONAL_BFS
        //     );
        // }

        // // Backbone(1SCD, LINOEXT)
        // if (i == __COUNTER__) {
        //     index = new BackboneIndex(
        //             graph.get(),
        //             localDist,
        //             BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER,
        //             BackboneEdgeCreationMethod::BFS,
        //             BackboneIndexingMethod::LANDMARK_NO_EXTENSIONS,
        //             LocalSearchMethod::UNIDIRECTIONAL_BFS
        //     );
        // }

        // Backbone(1SCD, LIEXT)
        if (i == __COUNTER__) {
            index = new BackboneIndex(
                    graph.get(),
                    localDist,
                    BackboneVertexSelectionMethod::ONE_SIDE_CONDITION_DEGREE_ORDER,
                    BackboneEdgeCreationMethod::BFS,
                    BackboneIndexingMethod::LANDMARK_ALL_EXTENSIONS,
                    LocalSearchMethod::UNIDIRECTIONAL_BFS
            );
        }

        // // Zou
        // if( i == __COUNTER__ )
        // {
        //     index = new Zou(graph);
        // }


        if( index->didCompleteBuilding() == true )
        {
          string indexName = index->getIndexTypeAsString();
          cout << "runTestsPerIndex index=" << indexName;
          double indexingTime = index->getIndexConstructionTimeInSec();
          unsigned long size = index->getIndexSizeInBytes();
          cout << indexName << " has index size (byte): " << size << endl;
          cout << indexName << " required index construction time (s): " << indexingTime << endl;

          methodNames.push_back( indexName );
          int status = runTestsPerIndex(index, queryTimes, queryTimeSums, indexSizes, indexTimes, noOfQuerySets, edge_file , size, indexingTime );
          if( status == 1 )
          {
              cout << "Error with index=" << indexName << endl;
              return 1;
          }

	  if (index->getIndexType() == indexns::IndexType::Backbone) {
		BackboneIndex* bi = (BackboneIndex*)index;
		nameToBackboneVertexRatio[indexName] = ((double)bi->getBackBoneVertices().size() / (double)N);
		nameToBackboneEdgeRatio[indexName] = ((double)bi->getBackBone().getNumberOfEdges() / (double)M);
	  }

        }
        else
        {
          string indexName = index->getIndexTypeAsString();
          cout << "No experiments for index: " << indexName << " as it did not complete building the index successfully.";
        }

        delete index;
    }

    cout << "--- Writing out the data ---" << endl;

    // Write out the data
    fstream myfile;
    myfile.open (output_file, std::fstream::in | std::fstream::out | std::ofstream::trunc);

    if( !myfile.is_open() )
    {
        std::cerr<<"Failed to open file : "<<SYSERROR()<<std::endl;
        exit(1);
    }

    // index statistics
    myfile << "method name, index construction time (s), index size (MB)\n";
    cout << "method name & index construction time (s) & index size (MB) \\\\ \n";
    for(int i = 0; i < methodNames.size(); i++)
    {
        double sizeInMB = max(indexSizes.at(i) / 1000000.0, 0.01);
        cout << methodNames[i] << " & " << print_digits(indexTimes[i],2) << " & " << print_digits(sizeInMB,2) << "\\\\ \n";
        myfile << methodNames[i] << "," << print_digits(indexTimes[i],2) << "," << print_digits(sizeInMB,2) << "\n";
    }

    myfile << "\n\n";

    for(int i = 0; i < methodNames.size(); i++)
    {
        for(int j = 0; j < (noOfQuerySets); j++)
        {
            myfile << methodNames[i] << "_queryset-" << j << ",";
        }
    }

    myfile << "\n";

    // write out query times
    // queryTimes[i][j][k] gives the time for method i, queryset

    for(int k = 0; k < (2*noOfQueries); k++)
    {
        for(int i = 0; i < methodNames.size(); i++)
        {
            for(int j = 0; j < noOfQuerySets; j++)
            {
                double d = queryTimes[i][j][k];
                myfile << print_digits(d, 11);
                if( j < noOfQuerySets-1 )
                    myfile << ",";
            }

            if( i < noOfMethods-1 )
                myfile << ",";
        }

        myfile << endl;
    }

    myfile << endl;

    // compute mean and std dev query time for true and false
    myfile << ",";
    for(int i = 0; i < methodNames.size(); i++)
    {
        for(int j = 0; j < noOfQuerySets; j++)
        {
            myfile << methodNames[i] << "-" << j;

            if( j < noOfQuerySets-1 )
                myfile << ",";
        }

        if( i < noOfMethods-1 )
            myfile << ",";
    }

    myfile << endl;

    // compute mean and sd per method and query set
    vector< vector < double > > numbers = vector< vector < double > > ();

    for(int i = 0; i < methodNames.size(); i++)
    {
        for(int j = 0; j < noOfQuerySets; j++)
        {
            vector < double > A = vector < double >();

            double trueMean = 0.0;
            double trueSD = 0.0;
            double falseMean = 0.0;
            double falseSD = 0.0;

            for(int k = 0; k < noOfQueries; k++)
            {
                trueMean += queryTimes[i][j][k];
            }

            for(int k = 0; k < noOfQueries; k++)
            {
                falseMean += queryTimes[i][j][k + noOfQueries];
            }

            trueMean /= noOfQueries;
            falseMean /= noOfQueries;

            for(int k = 0; k < noOfQueries; k++)
            {
                trueSD += pow(queryTimes[i][j][k] - trueMean, 2);
            }

            trueSD /= noOfQueries;
            trueSD = sqrt(trueSD);

            for(int k = 0; k < noOfQueries; k++)
            {
                falseSD += pow(queryTimes[i][j][k + noOfQueries] - falseMean, 2);
            }

            falseSD /= noOfQueries;
            falseSD = sqrt(falseSD);

            A.push_back( trueMean );
            A.push_back( falseMean );
            A.push_back( trueSD );
            A.push_back( falseSD );

            numbers.push_back(A);
        }
    }

    string names[] = { "mean-true", "mean-false", "sd-true", "sd-false" };
    for(int k = 0; k < 4; k++)
    {
        myfile << names[k] << ",";
        cout << names[k] << " & ";

        for(int i = 0; i < methodNames.size(); i++)
        {
            for(int j = 0; j < noOfQuerySets; j++)
            {
                myfile << numbers[i * noOfQuerySets + j][k];
                cout << numbers[i * noOfQuerySets + j][k];

                if( j < noOfQuerySets-1 )
                {
                    myfile << ",";
                    cout << " & ";
                }
            }

            if( i < noOfMethods-1 )
            {
                myfile << ",";
                cout << " & ";
            }
        }

        myfile << endl;
        cout << endl;
    }

    myfile << endl;

    // Compute speed-ups compared BFS, the first method
    // The other methods should be faster
    for (int i = 0; i < methodNames.size(); i++) {
      for (int j = 0; j < noOfQuerySets; j++) {
        myfile << methodNames[i] << "-su-true-" << j << ",";
        myfile << methodNames[i] << "-su-false-" << j;

        if (j < noOfQuerySets - 1)
          myfile << ",";
      }

      if (i < noOfMethods - 1)
        myfile << ",";
    }

    myfile << endl;

    cout << "method name";
    for(int i = 0; i < noOfQuerySets; i++)
    {
        cout << " & queryset (true) " << i << " & queryset (false) " << i;
    }
    cout << "\\\\ \n";

    for(int i = 0; i < methodNames.size(); i++)
    {
        cout << methodNames[i];

        for(int j = 0; j < noOfQuerySets; j++)
        {
            double a = queryTimeSums[2*noOfQuerySets*i + j*2];
            double b = queryTimeSums[2*noOfQuerySets*i + j*2 + 1];

            double su1 = queryTimeSums[j*2] / a;
            double su2 = queryTimeSums[j*2 + 1] / b;

            cout << " & " << print_digits(su1,2) << " & " << print_digits(su2,2);
            //cout << "a=" << a << ",su1=" << su1 << endl;
            //cout << "b=" << b << ",su2=" << su2 << endl;

            myfile << print_digits(su1,2) << "," << print_digits(su2,2);

            if( j < noOfQuerySets-1 )
                myfile << ",";
        }
        cout << "\\\\ \n";

        if( i < noOfMethods-1 )
            myfile << ",";

    }

    // TODO make more robust
    cout << "\n";
    for(int i = 0; i < methodNames.size(); i++)
    {
	string methodName = methodNames[i];
	if (nameToBackboneVertexRatio.count(methodName)) {
	    cout << methodName << " |V*|/|V|: " << nameToBackboneVertexRatio[methodName] << endl;
	}
    }
    cout << "\n";
    for(int i = 0; i < methodNames.size(); i++)
    {
	string methodName = methodNames[i];
	if (nameToBackboneEdgeRatio.count(methodName)) {
	    cout << methodName << " |E*|/|E|: " << nameToBackboneEdgeRatio[methodName] << endl;
	}
    }


    myfile.flush();
    myfile.close();

    exit(EXIT_SUCCESS);
}
