#include "BGLL.h"

void BGLL::Initialize(set<WEdge> edges)
{
	for (set<WEdge>::iterator it = edges.begin(); it != edges.end(); it++)
	{
		add_edge(it->u, it->v, it->w, m_wGraph);
	}
}

double BGLL::modularity_gain(int node, int comm, int dnodecomm) {
	//assert(node >= 0 && node<size);

	double totc = (double)weighted_degree[comm];
	double degc = (double)nodeDegree[node];
	double m2 = gTotalWeight;
	double dnc = (double)dnodecomm;

	double gain = (dnc - totc*degc / m2);
	return gain;
}

void BGLL::remove(int node, int comm) {
	//assert(node >= 0 && node<size);

	//nodeDegree[comm] -= weighted_degree[node];

	//m_comBGLL[node] = -1;
}

void BGLL::insert(int node, int comm) {
	//assert(node >= 0 && node<size);

	/*nodeDegree[comm] += weighted_degree[node];
	selfLoopWeight[comm] += 2 * dnodecomm + selfLoopWeight[node];
	m_comBGLL[node] = comm;*/
}