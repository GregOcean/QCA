#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility> 
#include <algorithm> 
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
using namespace std;
using namespace boost;

typedef property<edge_weight_t, int> EdgeWeightProperty;
typedef boost::adjacency_list<listS, vecS, undirectedS, no_property, EdgeWeightProperty> WGraph;
typedef graph_traits<WGraph>::adjacency_iterator  AdjIterator;								//迭代器

struct WEdge{
	int u;
	int v;
	int w;
	bool operator < (const WEdge &rhs) const
	{
		if (u < rhs.u) return true;
		else if (u == rhs.u) return v < rhs.v;
		return false;
	}
};

class BGLL{
public:
	BGLL()
	{
		gTotalWeight = 0;
	}

private:
	WGraph m_wGraph;

public:
	map<int, int> m_comBGLL;				//节点->社团号

public:
	void Initialize( set<WEdge> edges);

public:
	int gTotalWeight;						//图的总度数
	map<int, int> selfLoopWeight;			//社团的自我度数（自己到自己边的权）
	map<int, int> nodeDegree;				//点的周围的所有边的权和
	map<int, int> weighted_degree;			//一个社团的周围所有边的权的和
	map<int, map<int, int>> neighbor;		//一个<社团,<社团邻居，到邻居的边权>>

public:
	void remove(int node, int comm);
	void insert(int node, int comm);

public:
	double modularity_gain(int node, int comm, int dnodecomm);
};
