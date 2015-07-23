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
};																								//�߽ṹ

typedef adjacency_list<listS, vecS, undirectedS> Graph;											//ͼ���ݽṹ
typedef graph_traits<Graph>::adjacency_iterator  AdjacencyIterator;								//������

class CQCA{
private:
	Graph m_graph;
	std::map<unsigned int , unsigned int> m_community;											//nID,communityID��ÿ���ڵ��Ӧ��C
	std::map<unsigned int, std::set<int>> m_comToVertex;										//communityID,nID��ÿ��C�еĽڵ�
	
public:
	/*ʹ��ǰ�ȳ�ʼ��*/
	void InitGraph(std::list<Edge> & edges, std::map<unsigned int, unsigned int > &community);
	/*�ĸ�Ŀ�꺯��*/
	std::map<unsigned int, unsigned int> newNode(int Vertex, std::set<int> & Edges);
	std::map<unsigned int, unsigned int> newEdge(Edge newEdge);
	std::map<unsigned int, unsigned int> removeNode(int Vertex);
	std::map<unsigned int, unsigned int> removeEdge(Edge & edge);

private:
	TVec<TIntV> ConstructClique(std::set<int> & vertexes, std::set<Edge> & edges);				//��ʼ��Clique��ͼ�����ؽ��
	void GetCliqueEdge(std::set<int> vertexes, std::set<Edge> & edges);							//��ȡClique�����������ͼ�ı�
private:
	double Force_KeepInV(int CommunityID, int Vertex);
	double Force_BringOutV(int CommunityID, int Vertex);
	double Force_MaxOutV(int Vertex, int & comID);												//comID����Vertex������������

private:
	int deltaQuCD(int u, int v);

private:
	int DegreeOfVertex(int Vertex);														//��Ķ���
	int DegreeOfCommunity(int CommunityID);												//���������е�Ķ�����
	int DegreeOutCommunity(int CommunityID);											//ĳ����֮���Ķ�����
	
	int EdgesOfVtoC(int Vertex, int CommunityID, bool numOrbool = true);				//�㵽���ŵı���
	std::set<int> ComsAroundV(int Vertex);												//����������֮�����ڵ�����

	int NumEdgeBetweenC(int C1,int C2);													//�����ż�ı���

private:
	int GetGraphNodeNum();																//�õ�ͼ���ܽڵ���
	int GetGraphEdgeNum();																//�õ�ͼ���ܱ���

public:
	void printGraph();																	//��ӡͼ���ڽӱ�

public:
	void CombineCom(int Csrc , int Cdest);												//�ϲ�Csrc���ŵ�Cdest
};

