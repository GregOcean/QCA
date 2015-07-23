/************
by gregocean
2015.7 BIT
copyright 
*************/
#include "cliques.h"

#include <utility> 
#include <algorithm> 
#include <iostream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
using namespace std;
using namespace boost;

struct Edge{
	int uID;
	int vID;
	Edge()
	{
		uID = 0;
		vID = 0;
	}

	bool operator < (const Edge &rhs) const
	{
		if (uID < rhs.uID) return true;
		else if (uID == rhs.uID) return vID < rhs.vID;
		return false;
	}
};																								//边结构

typedef adjacency_list<listS, vecS, undirectedS> Graph;											//图数据结构
typedef graph_traits<Graph>::adjacency_iterator  AdjacencyIterator;								//迭代器

class CQCA{
private:
	Graph m_graph;
	std::map<unsigned int , unsigned int> m_community;											//nID,communityID存每个节点对应的C
	std::map<unsigned int, std::set<int>> m_comToVertex;										//communityID,nID存每个C中的节点
	
public:
	/*使用前先初始化*/
	void InitGraph(std::list<Edge> & edges, std::map<unsigned int, unsigned int > &community);
	/*四个目标函数*/
	std::map<unsigned int, unsigned int> newNode(int Vertex, std::set<int> & Edges);
	std::map<unsigned int, unsigned int> newEdge(Edge newEdge);
	std::map<unsigned int, unsigned int> removeNode(int Vertex);
	std::map<unsigned int, unsigned int> removeEdge(Edge & edge);

private:
	TVec<TIntV> ConstructClique(std::set<int> & vertexes, std::set<Edge> & edges);				//初始化Clique子图并返回结果
	void GetCliqueEdge(std::set<int> vertexes, std::set<Edge> & edges);							//提取Clique步骤的输入子图的边
private:
	double Force_KeepInV(int CommunityID, int Vertex);
	double Force_BringOutV(int CommunityID, int Vertex);
	double Force_MaxOutV(int Vertex, int & comID);												//comID返回Vertex最终所属社团

private:
	int deltaQuCD(int u, int v);

private:
	int DegreeOfVertex(int Vertex);														//点的度数
	int DegreeOfCommunity(int CommunityID);												//社团中所有点的度数和
	int DegreeOutCommunity(int CommunityID);											//某社团之外点的度数和
	
	int EdgesOfVtoC(int Vertex, int CommunityID, bool numOrbool = true);				//点到社团的边数
	std::set<int> ComsAroundV(int Vertex);												//点所在社团之外相邻的社团

	int NumEdgeBetweenC(int C1,int C2);													//两社团间的边数

private:
	int GetGraphNodeNum();																//得到图的总节点数
	int GetGraphEdgeNum();																//得到图的总边数

public:
	void printGraph();																	//打印图的邻接表

public:
	void CombineCom(int Csrc , int Cdest);												//合并Csrc社团到Cdest
};

