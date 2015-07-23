#include "CQCA.h"
#include "BGLL.h"
void CQCA::InitGraph(std::list<Edge> & edges, std::map<unsigned int, unsigned int > &community)					//初始化m_graph、m_community、m_comToVertex
{
	for (std::list<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)									//插入边
	{	//cout << it->uID << " " << it->vID << endl;
		assert(it->uID >= 0 && it->vID > 0);
		add_edge(it->uID, it->vID, m_graph);
	}
	
	std::set<int> comID;																						//所有社团号，用作查重
	
	for (std::map<unsigned int, unsigned int >::iterator it = community.begin(); it != community.end(); ++it)
	{	//cout << it->first << " " << it->second << endl;
		m_community.insert(std::make_pair(it->first, it->second));
		
		if (comID.end() == comID.find(it->second))																	//构造 聚类->其中节点 的映射结构
		{
			comID.insert(it->second);
			std::set<int> vSet;
			m_comToVertex.insert(std::make_pair(it->second, vSet));
		}	
	}	

	for (std::map<unsigned int, unsigned int >::iterator it = community.begin(); it != community.end(); ++it)
	{
		m_comToVertex[it->second].insert(it->first);
	}
}

double CQCA::Force_KeepInV(int CommunityID, int Vertex)
{
	int Euc = EdgesOfVtoC(Vertex, CommunityID);
	int Du = DegreeOfVertex(Vertex);
	int Dc = DegreeOfCommunity(CommunityID);
	int M = GetGraphEdgeNum();

	return (double)(Euc - (Du * (Dc - Du)) / (double)(2 * M));
}

double CQCA::Force_BringOutV(int CommunityID, int Vertex)
{
	int Eus = EdgesOfVtoC(Vertex, CommunityID);
	int Du = DegreeOfVertex(Vertex);
	int Dos = DegreeOutCommunity(CommunityID);
	int M = GetGraphEdgeNum();

	return (double)(Eus - (Du * Dos) / (double)(2 * M));
}

double CQCA::Force_MaxOutV(int Vertex, int & comID)
{
	std::set<int> comSet = ComsAroundV(Vertex);
	double maxResult = -DBL_MAX;																					//先赋最小值
	double tempFbov;
	for (std::set<int>::iterator it = comSet.begin(); it != comSet.end(); ++it)
	{
		tempFbov = Force_BringOutV(*it, Vertex);
		if (tempFbov > maxResult)
		{
			maxResult = tempFbov;
			comID = *it;
		}
	}
	return maxResult;
}

std::map<unsigned int, unsigned int> CQCA::newNode(int Vertex, std::set<int> & Edges)
{
	m_graph.added_vertex(Vertex);
	m_community.insert(std::make_pair(Vertex,m_comToVertex.size()));
	std::set<int> newVSet;
	newVSet.insert(Vertex);
	m_comToVertex.insert(std::make_pair(m_comToVertex.size(), newVSet));
	if (Edges.empty())																								//没有边，直接返回
	{
		return m_community;
	}
	else
	{
		std::set<int> vertexJoinCv;
		for (std::set<int>::iterator it = Edges.begin(); it != Edges.end(); ++it)									//给图加边
		{
			add_edge(Vertex, *it, m_graph);
		}

		std::set<int> nearV = ComsAroundV(Vertex);

		for (std::set<int>::iterator it = Edges.begin(); it != Edges.end(); ++it)									//对u所有相邻的点判断F		
		{
			if (Force_BringOutV(m_community[Vertex], *it) > Force_KeepInV(m_community[*it], *it))					//若C(u)大，则对相邻的该社团的所有点判断吸引力		
			{
				if (nearV.find(m_community[*it]) != nearV.end())													//每个满足的社团判断一次
				{
					std::set<int> vertexesInNearCom = m_comToVertex[m_community[*it]];
					for (std::set<int>::iterator itC = vertexesInNearCom.begin(); itC != vertexesInNearCom.end(); )
					{
						if (Force_BringOutV(m_community[Vertex], *itC) <= Force_KeepInV(m_community[*itC], *itC))	//不满足 条件（>）则将点删除
						{
							vertexesInNearCom.erase(itC++);
						}
						else
							itC++;
					}
					nearV.erase(m_community[*it]);

					std::set_union(vertexJoinCv.begin(), vertexJoinCv.end(), vertexesInNearCom.begin(), vertexesInNearCom.end(), std::inserter(vertexJoinCv, vertexJoinCv.begin()));			//合并所有加入C(u)的点
				}
			}
		}

		if (vertexJoinCv.empty())
		{
			return m_community;
		}
		else 
		{																																					//更新m_community和m_comToVertex	
			for (std::set<int>::iterator it = vertexJoinCv.begin(); it != vertexJoinCv.end(); ++it)															//新点加入C(u)
			{
				m_comToVertex[m_community[*it]].erase(*it);
				m_comToVertex[m_community[Vertex]].insert(*it);
				m_community[*it] = m_community[Vertex];
			}

			int comID = -1;
			if (Force_MaxOutV(Vertex, comID) > Force_KeepInV(m_community[Vertex], Vertex))
			{
				for (std::set<int>::iterator it = m_comToVertex[m_community[Vertex]].begin(); it != m_comToVertex[m_community[Vertex]].end(); ++it)				//将C(u)中点移到社团comID中	
				{
					m_community[*it] = comID;
					m_comToVertex[comID].insert(*it);
				}
				m_comToVertex.erase(m_community[Vertex]);
			}
		}
	}
	return m_community;
}

int CQCA::deltaQuCD(int u, int v)
{
	int C = m_community[u];
	int D = m_community[v];;
	int M = GetGraphEdgeNum();
	
	int euD = EdgesOfVtoC(u, D);
	int euC = EdgesOfVtoC(u, C);
	
	int dC = DegreeOfCommunity(C);
	int dD = DegreeOfCommunity(D);
	int du = DegreeOfVertex(u);
	
	return 4 * (M + 1)*(euD + 1 - euC) + euC * (2 * dD - 2 * du - euC) - 2 * (du + 1)* (du + 1 + dD - dC);
}

std::map<unsigned int, unsigned int> CQCA::newEdge(Edge newEdge)
{
	add_edge(newEdge.uID, newEdge.vID, m_graph);
	if (m_community[newEdge.uID] == m_community[newEdge.vID])
	{
		return m_community;
	}
	else
	{
		int qu = deltaQuCD(newEdge.uID, newEdge.vID);
		int qv = deltaQuCD(newEdge.vID, newEdge.uID);
		if (qu < 0 && qv < 0)
		{
			return m_community;
		}
		else
		{
			int w = qu > qv ? newEdge.uID : newEdge.vID;
			m_comToVertex[m_community[w]].erase(w);
			if (w == newEdge.uID){																					//u加入v
				m_community[w] = m_community[newEdge.vID];
			}
			else{
				m_community[w] = m_community[newEdge.uID];
			}
			m_comToVertex[m_community[w]].insert(w);

			std::pair<AdjacencyIterator, AdjacencyIterator> ai = adjacent_vertices(w, m_graph);						//对点u邻居判断所属社团
			for (AdjacencyIterator ait = ai.first; ait != ai.second; ++ait) {
				if (Force_BringOutV(m_community[w], *ait) > Force_KeepInV(m_community[*ait], *ait))
				{
					m_comToVertex[m_community[*ait]].erase(*ait);
					m_comToVertex[m_community[w]].insert(*ait);
					m_community[*ait] = m_community[w];
				}
			}
		}
	}
	return m_community;
}

TVec<TIntV> CQCA::ConstructClique(std::set<int> & vertexes, std::set<Edge> & edges)							//形成临时图结构
{
	const int OverlapSz = Env.GetIfArgPrefixInt("-k:", 2, "Min clique overlap");							//OverlapSz						
	PUNGraph G = TUNGraph::New();

	for (std::set<int>::iterator it = vertexes.begin(); it != vertexes.end(); it++)
	{
		G->AddNode(*it);
	}
	for (std::set<Edge>::iterator it = edges.begin(); it != edges.end(); it++)
	{
		G->AddEdge(it->uID, it->vID);
	}
	TVec<TIntV> CmtyV;
	TCliqueOverlap::GetCPMCommunities(G, OverlapSz + 1, CmtyV);												//k = 2+1

	return CmtyV;
}

void CQCA::GetCliqueEdge(std::set<int> vertexes, std::set<Edge> & edges)													//点集不能为空
{
	for (std::set<int>::iterator its = vertexes.begin(); its != vertexes.end(); ++its)										//对每个社团内的点循环
	{
		std::pair<AdjacencyIterator, AdjacencyIterator> ai = adjacent_vertices(*its, m_graph);								//对该点所有边循环	
		for (AdjacencyIterator itd = ai.first; itd != ai.second; ++itd) 
		{																													//判断该边的另一个点是否为节点集中的点
			for (std::set<int>::iterator it = vertexes.begin(); it != vertexes.end(); ++it)
			{	std::cout << *its << " "<< *itd << "\t";
				if (*it == *itd && *its < *itd)
				{
					Edge edgeTmp;

					edgeTmp.uID = *its;
					edgeTmp.vID = *itd;
					
					if (edges.find(edgeTmp) == edges.end())
					{
						edges.insert(edgeTmp);
					}

					break;
				}
			}
		}
	}
}

std::map<unsigned int, unsigned int> CQCA::removeNode(int Vertex)
{
	int comID = m_community[Vertex];
	if (DegreeOfVertex(Vertex) == 0)
	{
		remove_vertex(Vertex, m_graph);
		m_community.erase(Vertex);
	}
	else
	{	// 共同构成临时图G的输入
		set<int> tempNodeSet;																				// 相邻节点集
		set<Edge> tempEdgeSet;																				// 相邻节点间构成的边集
		set<int> nearNodeSet;

		std::pair<AdjacencyIterator, AdjacencyIterator> ai = adjacent_vertices(Vertex, m_graph);				// 找到每个相邻节点
		for (AdjacencyIterator ait = ai.first; ait != ai.second; ++ait) {
			nearNodeSet.insert(*ait);
		}
		
		clear_vertex(Vertex, m_graph);																			//删除顶点相关边
		//remove_vertex(Vertex-1, m_graph);																			//删除顶点
		m_comToVertex[m_community[Vertex]].erase(Vertex);
		m_community.erase(Vertex);

		if (nearNodeSet.size() == 1){
			return m_community;
		}
		
		set<int> cliqueS;																				//待BGLL的新Clique社团号
		while (true)
		{
			bool nextClique = false;
			tempNodeSet = m_comToVertex[comID];
			GetCliqueEdge(tempNodeSet, tempEdgeSet);
			TVec<TIntV> CmtyV = ConstructClique(tempNodeSet, tempEdgeSet);									
			for (int i = 0; i < CmtyV.Len(); i++)
			{
				for (int j = 1; j < CmtyV[i].Len(); j++)
				{
					cout << (int)CmtyV[i][j].Val << " ";
				}
				cout << endl;
			}
			if (CmtyV.Len() > 0)																					
			{
				for (int i = 0; i < CmtyV.Len(); i++)
				{
					for (int j = 0; j < CmtyV[i].Len(); j++)														//从最大的clique开始找是否有C(u)	
					{
						if (nearNodeSet.find(CmtyV[i][j].Val) != nearNodeSet.end())									//如果clique中有C(u)中的点
						{
							nextClique = true;

							std::set<int> newCom;
							int cliqueV = 0;
							for (int k = 0; k < CmtyV[i].Len(); k++)
							{
								cliqueV = CmtyV[i][k].Val;
								
								nearNodeSet.erase(cliqueV);
								newCom.insert(cliqueV);
								m_comToVertex[m_community[cliqueV]].erase(cliqueV);
								m_community[cliqueV] = m_comToVertex.size();
							}
							cliqueS.insert(m_comToVertex.size());
							m_comToVertex[m_comToVertex.size()] = newCom;
							tempEdgeSet.clear();
							break;
						}																							//将最大的clique类的节点移出组成新社团

					}
					if (nextClique)	{
						nextClique = false;
						break;
					}
					else{																							//如果所有clique中都不包含N(u)中的点。
						break;
					}
				}
			}
			else
				break;
		}
		//tempNodeSet中都是单独节点
		for (std::set<int>::iterator it = tempNodeSet.begin(); it != tempNodeSet.end(); it++)							//其余节点各自判断所属社团
		{
			int communityID = -1;
			Force_MaxOutV(*it, communityID);

			m_comToVertex[m_community[*it]].erase(*it);
			m_community[*it] = communityID;
			m_comToVertex[communityID].insert(*it);
		}
		//初始化BGLL图
		//给出边集
		BGLL bgllGraph;

		set<WEdge> wEdges;
		for (map<unsigned int, std::set<int>>::iterator itc = m_comToVertex.begin(); itc != m_comToVertex.end(); ++itc)
		{
			for (map<unsigned int, std::set<int>>::iterator its = m_comToVertex.begin(); its != m_comToVertex.end(); ++its)
			{
				if (itc->first < its->first)
				{
					WEdge e;
					e.u = itc->first;
					e.v = its->first;
					e.w = NumEdgeBetweenC(e.u, e.v);
					if (e.w > 0)
					{
						wEdges.insert(e);

						bgllGraph.nodeDegree[e.u] += e.w;
						bgllGraph.nodeDegree[e.v] += e.w;

						bgllGraph.gTotalWeight += e.w;
						
						bgllGraph.neighbor[e.u].insert(make_pair(e.v, e.w));
						bgllGraph.neighbor[e.v].insert(make_pair(e.u, e.w));
					}
				}
			}
		}
		
		
		bgllGraph.Initialize(wEdges);
		//算in 和 community
		for (map<unsigned int, std::set<int>>::iterator it = m_comToVertex.begin(); it != m_comToVertex.end(); ++it)
		{
			bgllGraph.m_comBGLL[it->first] = it->first;

			set<int> SNodeSet = m_comToVertex[it->first];
			set<Edge> SEdgeSet;
			GetCliqueEdge(SNodeSet, SEdgeSet);
			bgllGraph.selfLoopWeight[it->first] = SEdgeSet.size();
		}
		
		//算tot
		for (map<unsigned int, std::set<int>>::iterator it = m_comToVertex.begin(); it != m_comToVertex.end(); ++it)
		{
			bgllGraph.weighted_degree[it->first] = 0;
			if (m_comToVertex[it->first].size() > 0)
			{
				set<int> nearComID = ComsAroundV(*m_comToVertex[it->first].begin());
				for (set<int>::iterator itc = nearComID.begin(); itc != nearComID.end(); itc++)
				{
					bgllGraph.weighted_degree[it->first] += NumEdgeBetweenC(it->first, *itc);
				}
			}
		}
		//按照BGLL方法对S进行社团归属
		map<int, int> comMove;
		for (set<int>::iterator its = cliqueS.begin(); its != cliqueS.end(); ++its)
		{
			map<int, int> ncomm = bgllGraph.neighbor[*its];
			int node_comm = bgllGraph.m_comBGLL[*its];
			//bgllGraph.remove(*its, node_comm);

			// compute the nearest community for node
			// default choice for future insertion is the former community
			int best_comm = node_comm;
			int best_nblinks = 0;														//ncomm.find(node_comm)->second;
			double best_increase = 0.;													//modularity_gain(node, best_comm, best_nblinks);

			for (map<int, int>::iterator it = ncomm.begin(); it != ncomm.end(); it++) {
				double increase = bgllGraph.modularity_gain(*its, it->first, it->second);
				if (increase > best_increase) {
					best_comm = it->first;
					best_nblinks = it->second;
					best_increase = increase;
				}
			}

			//返回值（set<消失的社团号+加入的社团号>）
			if (*its < best_comm)
			{
				comMove[*its] = best_comm;
			}
			else
			{
				comMove[best_comm] = *its;
			}
		}

		for (map<int, int>::iterator it = comMove.begin(); it != comMove.end(); ++it)
		{
			CombineCom(it->first, it->second);
		}
	}
	return m_community;
}

std::map<unsigned int, unsigned int> CQCA::removeEdge(Edge & edge)
{	
	remove_edge(edge.uID, edge.vID, m_graph);
	if (m_community[edge.uID] != m_community[edge.vID])															//两端属于不同社团
	{																											//在图中删除边
		return m_community;
	}
	int comID = m_community[edge.uID];

	if (DegreeOfVertex(edge.vID) == 0)																			//(u, v) is a single edge,删边后度就变为0
	{
		m_comToVertex[comID].erase(edge.vID);																	//将v移出C(u)，自立社团
		m_community[edge.vID] = m_comToVertex.size();
		std::set<int> vSet;
		vSet.insert(edge.vID);
		m_comToVertex[m_comToVertex.size()] = vSet;
	}
	else if (DegreeOfVertex(edge.uID) == 0)
	{	
		m_comToVertex[comID].erase(edge.uID);																	//将u移出C(v)，自立社团
		m_community[edge.uID] = m_comToVertex.size();
		std::set<int> uSet;
		uSet.insert(edge.uID);
		m_comToVertex[m_comToVertex.size()] = uSet;
	}
	else																																// 边的两个节点在同一个社团
	{
		std::set<Edge> tempEdgeSet;																						// 临时图G的输入
		
		GetCliqueEdge(m_comToVertex[comID], tempEdgeSet);
		TVec<TIntV> CmtyV = ConstructClique(m_comToVertex[comID], tempEdgeSet);
		
		if (CmtyV.Len() > 0)																											//将最大的伪类的节点移出组成新社团
		{
			std::set<int> newCom;
			for (int i = 0; i < CmtyV[0].Len(); i++) {
				cout << CmtyV[0][i].Val;
				newCom.insert(CmtyV[0][i].Val);
				m_comToVertex[comID].erase(CmtyV[0][i].Val);
				m_community[CmtyV[0][i].Val] = m_comToVertex.size();
				
			}
			m_comToVertex[m_comToVertex.size()] = newCom;
		}
		
		for (std::set<int>::iterator it = m_comToVertex[comID].begin(); it != m_comToVertex[comID].end(); it++)							//其余节点各自判断所属社团
		{
			int communityID = -1;
			Force_MaxOutV(*it, communityID);
			m_community[*it] = communityID;
			m_comToVertex[communityID].insert(*it);	
			
		}
		m_comToVertex[comID].clear();
	}
	
	return m_community;
}

int CQCA::DegreeOfVertex(int Vertex)
{
	return degree(Vertex, m_graph);
}

int CQCA::DegreeOfCommunity(int CommunityID)
{
	unsigned int comDegree = 0;
	std::set<int> vSet = m_comToVertex[CommunityID];
	for (std::set<int>::iterator it = vSet.begin(); it != vSet.end(); ++it)
	{
		comDegree += degree(*it, m_graph);
	}
	return comDegree;
}

int CQCA::DegreeOutCommunity(int CommunityID)
{
	return (2 * GetGraphEdgeNum()) - DegreeOfCommunity(CommunityID);
}

int CQCA::EdgesOfVtoC(int Vertex, int CommunityID, bool numOrbool)
{
	unsigned int edgeNumofVtoC = 0;
	std::pair<AdjacencyIterator, AdjacencyIterator> ai = adjacent_vertices(Vertex, m_graph);
	for (AdjacencyIterator ait = ai.first; ait!= ai.second; ++ait) {														//如果相邻节点在CommunityID中
		if (m_community[*ait] == CommunityID)
		{
			edgeNumofVtoC++;
			if (!numOrbool)	{
				return 1;
			}
		}
	}
	//boost::graph_traits<Graph>::edge_iterator  ei = edges(m_graph).first;	++ei;
	return edgeNumofVtoC;
}

std::set<int> CQCA::ComsAroundV(int Vertex)																					//NCu , 求V所在社团Cu周围的社团
{
	std::set<int> comsID;
	std::set<int> tmpVSet = m_comToVertex[m_community[Vertex]];
	for (std::map<unsigned int, std::set<int>>::iterator it = m_comToVertex.begin(); it != m_comToVertex.end(); ++it)		//循环所有社团
	{
		for (std::set<int>::iterator itor = tmpVSet.begin(); itor != tmpVSet.end(); ++itor)									//对C(u)中所有点
		{
			if (it->first != m_community[Vertex] && EdgesOfVtoC(*itor, it->first, false) > 0 )									//其邻接点在其他社团
			{
				comsID.insert(it->first);
				break;
			}
		}
	}
	return comsID;
}

int CQCA::GetGraphNodeNum()
{
	return num_vertices(m_graph);
}

int CQCA::GetGraphEdgeNum()
{
	return num_edges(m_graph);
}

int CQCA::NumEdgeBetweenC(int C1, int C2)
{
	int tot = 0;
	for (set<int>::iterator it = m_comToVertex[C1].begin(); it != m_comToVertex[C1].end(); ++it)
	{
		tot += EdgesOfVtoC(*it, C2);
	}
	return tot;
}

void CQCA::printGraph()
{
	Graph::vertex_iterator vertexIt, vertexEnd;
	Graph::adjacency_iterator neighbourIt, neighbourEnd;
	tie(vertexIt, vertexEnd) = vertices(m_graph);
	for (; vertexIt != vertexEnd; ++vertexIt)
	{
		cout << *vertexIt << " connect ";
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(*vertexIt, m_graph);
		for (; neighbourIt != neighbourEnd; ++neighbourIt)
			cout << *neighbourIt << " ";
		cout << "\n";
	}
	cout << endl;
}

void CQCA::CombineCom(int Csrc, int Cdest)
{
	set<int> nodesFromSrc = m_comToVertex[Csrc];
	for (set<int>::iterator it = nodesFromSrc.begin(); it != nodesFromSrc.end(); ++it)
	{
		m_comToVertex[Cdest].insert(*it);
		m_community[*it] = Cdest;
	}
}