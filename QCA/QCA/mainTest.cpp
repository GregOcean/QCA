#include "CQCA.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "cliques.h"
using namespace std;

int main(int, char*[])
{
	CQCA qcaObject;

	std::list<Edge>  edges;
	
	//����ͼ���ߣ������Ӧ�����ɣ������ж�������Ҫ��������
	ifstream inG("E:\\Jiang\\QCA\\GraphData_BGLL\\removeNode\\rmPaperGraphD.txt");

	if (!inG.is_open()){
		cout << "Error opening file"; exit(1);
	}
	struct Edge Etemp;
	while (!inG.eof())
	{
		inG >> Etemp.uID;
		inG >> Etemp.vID;
		edges.push_back(Etemp);
	}
	inG.close();

	//������������㣩
	std::map<unsigned int, unsigned int > community;
	ifstream inC("E:\\Jiang\\QCA\\GraphData_BGLL\\removeNode\\rmPaperGraphD_level1");
	if (!inC.is_open()){
		cout << "Error opening file"; exit(1);
	}
	int nID;
	int cID;
	while (!inC.eof())
	{
		inC >> nID;
		inC >> cID;
		community.insert(std::make_pair(nID, cID));
	}
	inC.close();
	//��ʼ��
	qcaObject.InitGraph(edges, community);
	
	//����
	std::map<unsigned int, unsigned int > resultCommunity;

	/*std::set<int> newNodeEdge;
	newNodeEdge.insert(15);
	newNodeEdge.insert(11);
	newNodeEdge.insert(13);
	newNodeEdge.insert(8);
	newNodeEdge.insert(10);
	newNodeEdge.insert(14);
	resultCommunity = qcaObject.newNode(16, newNodeEdge);*/

	
	/*Edge edge;
	edge.uID = 13;
	edge.vID = 6;
	resultCOmmunity = qcaObject.newEdge(edge);*/

	int Vertex = 6;
	resultCommunity = qcaObject.removeNode(Vertex);

	//Edge edge;
	//edge.uID = 1;
	//edge.vID = 4;
	//resultCommunity = qcaObject.removeEdge(edge);
	
	return 0;
}