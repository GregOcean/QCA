1、algo1的理解（对相邻点的条件判断，进而判断周围社团中的所有点？）
2、algo2的图示


CPM
BGLL

alg3和alg4
剩余节点的原社团完全消失？
从其余社团和新生成的S-clique中挑选？

极大完全子图识别

3和4中的剩余独立节点，用Force判断，需要foirce是正的么？

3中 clique不能自立社团一定要从属于某个社团么？（还是基于BGLL结果，deltaQ大于0就移动，否则自立）

	/*int qNb = deltaQuCD(*ait, w);																	
				int qu = deltaQuCD(w, *ait);
				if (qNb > 0 && qNb > qu)
				{
				m_comToVertex[m_community[*ait]].erase(*ait);
				m_comToVertex[m_community[newEdge.vID]].insert(*ait);
				m_community[*ait] = m_community[newEdge.vID];
				}*/

				/*int qNb = deltaQuCD(*ait, w);
				int qu = deltaQuCD(w, *ait);
				if (qNb > 0 && qNb > qu)
				{
				m_comToVertex[m_community[*ait]].erase(*ait);
				m_comToVertex[m_community[newEdge.uID]].insert(*ait);
				m_community[*ait] = m_community[newEdge.uID];
				}*/