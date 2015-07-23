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
typedef graph_traits<WGraph>::adjacency_iterator  AdjIterator;								//������

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
	map<int, int> m_comBGLL;				//�ڵ�->���ź�

public:
	void Initialize( set<WEdge> edges);

public:
	int gTotalWeight;						//ͼ���ܶ���
	map<int, int> selfLoopWeight;			//���ŵ����Ҷ������Լ����Լ��ߵ�Ȩ��
	map<int, int> nodeDegree;				//�����Χ�����бߵ�Ȩ��
	map<int, int> weighted_degree;			//һ�����ŵ���Χ���бߵ�Ȩ�ĺ�
	map<int, map<int, int>> neighbor;		//һ��<����,<�����ھӣ����ھӵı�Ȩ>>

public:
	void remove(int node, int comm);
	void insert(int node, int comm);

public:
	double modularity_gain(int node, int comm, int dnodecomm);
};
