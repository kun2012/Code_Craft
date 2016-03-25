
#include "route.h"
#include "lib_record.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <vector>
#include <queue>
using std::queue;
using std::set;
using std::vector;
using std::cout;
using std::endl;
using std::set;
using std::max;
using std::swap;

int edges[601][601][2];  //0 ==> weight, 1 ==> index
int edgesCnt[601];
int labels[601][601];
int N = 0;
int src = -1, dest = -1;
set<int> bridges;
int minDistance = INT_MIN;
bool is_time_out = false;
	
int tmpEdges[601][601];
void read_graph(char *topo[5000], int edge_num);
void read_demand(char *demand);
void getShortestPathSPFA(int start);
void getShortestPathBruteForce(int start);
void output_result();

int shortestPath[601];

set<int> excluded;
set<int> visitedB;
vector<int> path;
int bCnt;
int pathLen;
queue<int> que;
clock_t start_time;

void search_route(char *topo[5000], int edge_num, char *demand)
{
    start_time = clock();
    read_graph(topo, edge_num);
    read_demand(demand);


	path.push_back(0);  //distance
	path.push_back(src);  //source
	excluded.insert(src);
	
	if (bCnt <= 6 && N <= 20)
		getShortestPathBruteForce(src);
	else
		getShortestPathSPFA(src);
		
	output_result();
}

void output_result() {
    
    cout << minDistance << endl;

	for (int i = 0; i < pathLen - 1; i++)
        record_result(labels[shortestPath[i]][shortestPath[i + 1]]);
}

void read_graph(char *topo[5000], int edge_num)
{
    memset(tmpEdges, 0, sizeof(tmpEdges));
    for (int i = 0; i < edge_num; i++)
    {
        int linkID, srcID, destID, cost;
        sscanf(topo[i], "%d,%d,%d,%d", &linkID, &srcID, &destID, &cost);
        N = max(N, max(srcID, destID) + 1);
        if (tmpEdges[srcID][destID] == 0 || tmpEdges[srcID][destID] > cost)
        {
			tmpEdges[srcID][destID] = cost;
			labels[srcID][destID] = linkID;
		}
    }
    
    for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			if (tmpEdges[i][j] != 0) {
				edges[i][edgesCnt[i]][0] = tmpEdges[i][j];
				edges[i][edgesCnt[i]][1] = j;
				edgesCnt[i]++;
			}
}

void read_demand(char *demand)
{
    int i = 0;
    src = 0;
    while(demand[i] != ',')
    {
        src = src * 10 + demand[i] - '0';
        i++;
    }
    dest = 0;
    i++;
    while(demand[i] != ',')
    {
        dest = dest * 10 + demand[i] - '0';
        i++;
    }
    i++;
    while(demand[i] != '\n' && demand[i] != '\0')
    {
        int tmp = 0;
        while(demand[i] != '\0' && demand[i] != '|' && demand[i] != '\n')
        {
            tmp = tmp * 10 + demand[i] - '0';
            i++;
        }
        bridges.insert(tmp);
        i++;
    }
    printf("src = %d\n", src);
    printf("dest = %d\n", dest);
    
    /*
    set<int>::iterator it;
    for (it = bridges.begin(); it != bridges.end(); it++) {
        printf("%d ", *it);
    }
    printf("\n");
    */
}

bool color[601];

void getShortestPathSPFA(int start) {
    clock_t cur_time = clock();
	if (cur_time - start_time > 9900) {
		is_time_out = true;
		return ;
	}
	
	//distances[i] == 0 means: 1. start vertex, 2. unreabable
	int *distances = new int[N];
	int *parents = new int[N];
	memset(color, 0, sizeof(color));
	for (int i = 0; i < N; i++) {
	    distances[i] = 0;
	    parents[i] = 0;
	}
	
	if (!que.empty()) {
	    cout << "error que" << endl;
	    abort();
	}
	
	que.push(start);
	
	color[start] = true;
	parents[start] = start;
	
	//here, there may be a path like: start --> b1 --> b2
	while (!que.empty()) {
		int u = que.front();
		que.pop();
		color[u] = false;
		for (int c = 0; c < edgesCnt[u]; ++c) {
			int w = edges[u][c][0];
			int v = edges[u][c][1];
			if (excluded.find(v) != excluded.end()) continue;
			if (distances[v] == 0 || distances[v] > distances[u] + w) {
				distances[v] = distances[u] + w;
				parents[v] = u;
				if (color[v] == false) {
					que.push(v);
					color[v] = true;
				}
			}
		}
	}
	
	set<int>::iterator it;
	
	for (it = bridges.begin(); it != bridges.end(); ++it) {
	    int b = *it;
		if (visitedB.find(b) == visitedB.end() && distances[b] == 0) return;
	}
	if (distances[dest] == 0) return;
	if (distances[dest] + path[0] >= minDistance) return;
		
	//we can be sure that distances[dest] + path.get(0) < minDistance
	//and we update the minDistance and shortestPath
	if (bridges.size() == visitedB.size()) {
	    vector<int> res;
		int p = dest;
		while (parents[p] != p) {
			res.push_back(p);
			p = parents[p];
		}
		minDistance = path[0] + distances[dest];
		pathLen = path.size() + res.size() - 1;
		for (vector<int>::size_type i = 1; i < path.size(); ++i)
			shortestPath[i - 1] = path[i];
		for (vector<int>::size_type i = 0; i < res.size(); ++i)
			shortestPath[i + path.size() - 1] = res[res.size() - i - 1];
		return;
	}
		
	//here, we can use randomness
	int visitedBCnt = visitedB.size();
	int *tmpBridges = new int[bridges.size() - visitedBCnt];
	int tmpBCnt = 0;
	
	for (it = bridges.begin(); it != bridges.end(); ++it) {
	    int b = *it;
		if (visitedB.find(b) != visitedB.end() || (path[0] + distances[b] >= minDistance))
			continue;
		tmpBridges[tmpBCnt++] = b;
	}
	for (int i = 0; i < tmpBCnt; ++i) {
		int j = rand() % tmpBCnt;
		swap(tmpBridges[i], tmpBridges[j]);
	}
		
	vector<int> res;
	for (int bi = 0; bi < tmpBCnt; ++bi) {
		int b = tmpBridges[bi];
		
		if (4 < visitedBCnt && visitedBCnt + 5 < bCnt) {
			int p = rand() % visitedBCnt;
			if (p < visitedBCnt - 3)
				continue;
		}
		
		res.clear();
		int p = b;
		bool containsB = false;
		while (parents[p] != p) {
			if (p != b && bridges.find(p) != bridges.end()) {
				containsB = true;
				break;
			}
			res.push_back(p);
			p = parents[p];
		}
		//we use greedy tragedy. if start --> b1 --> b2, we won't select b2 (just b1)
		if (containsB) continue;
		
		//res does not contain start & path has already contained it
		path[0] += distances[b];
		for (int i = res.size() - 1; i >= 0; --i) {
			path.push_back(res[i]);
			excluded.insert(res[i]);
		}
		visitedB.insert(b);
		getShortestPathSPFA(b);
		
		if (is_time_out) {
		    //using return for speed but causing memory leak
		    return ;
		}
		
		visitedB.erase(b);
		for (int i = res.size() - 1; i >= 0; --i) {
			path.pop_back();
			excluded.erase(res[i]);
		}
		path[0] -= distances[b];
	}
	delete[] tmpBridges;
	delete[] parents;
	delete[] distances;
}


void getAllPaths(int start, int end, vector<int> &path, vector<vector<int> > &allPaths) {
	for (int c = 0; c < edgesCnt[start]; ++c) {
		int w = edges[start][c][0];
		int v = edges[start][c][1];
		if (v == end) {
			vector<int> r(path);
			r[0] += w;
			r.push_back(v);
			allPaths.push_back(r);
		} else if (excluded.find(v) == excluded.end() && bridges.find(v) == bridges.end()) {
			excluded.insert(v);
			path.push_back(v);
			path[0] += w;
			getAllPaths(v, end, path, allPaths);
			path[0] -= w;
			path.pop_back();
			excluded.erase(v);
		}
	}
}


void getNextBridges(int start, set<int> &nextB) {

    memset(color, 0, sizeof(color));
    
    if (!que.empty()) {
        cout << "que not empty error" << endl;
        abort();
    }
    
	que.push(start);
	color[start] = true;
	
	while (!que.empty()) {
		int u = que.front();
		que.pop();
		if (u != start && bridges.find(u) != bridges.end()) continue;
		for (int c = 0; c < edgesCnt[u]; ++c) {
			int v = edges[u][c][1];
			if (excluded.find(v) != excluded.end()) continue;
			if (color[v] == false) {
				que.push(v);
				color[v] = true;
				if (bridges.find(v) != bridges.end() && visitedB.find(v) == visitedB.end())
					nextB.insert(v);
			}
		}
	}
}


void getShortestPathBruteForce(int start) {

    vector<vector<int> > allPaths;

	//next is destination
    if (visitedB.size() == (unsigned int)bCnt) {
		vector<int> curPath;
		curPath.push_back(0);
		curPath.push_back(start);
		getAllPaths(start, dest, curPath, allPaths);
		
		if (allPaths.size() == 0) return;
		
		int mind = INT_MAX, ti = -1;
		//path[0] returns its distance
		for (size_t i = 0; i < allPaths.size(); i++) {
			if (mind > allPaths[i][0]) {
				mind = allPaths[i][0];
				ti = i;
			}
		}
		
	    if (minDistance > path[0] + mind) {
		    pathLen = path.size() + allPaths[ti].size() - 3;
		    for (size_t i = 1; i < path.size(); ++i)
			    shortestPath[i - 1] = path[i];
		    for (size_t i = 2; i < allPaths[ti].size(); ++i)
			    shortestPath[path.size() + i - 3] = allPaths[ti][i];
		    minDistance = path[0] + mind;
	    }
	    return;
	}
		
		
	//next vertex is bridge
	set<int> nextB;
	getNextBridges(start, nextB);
	
	set<int>::iterator it;
	for (it = nextB.begin(); it != nextB.end(); ++it) {
	    int b = *it;
		
		vector<int> curPath;
		curPath.push_back(0);
		curPath.push_back(start);
		getAllPaths(start, b, curPath, allPaths);
		
		if (allPaths.size() == 0) return;
		
		for (size_t j = 0; j < allPaths.size(); j++) {
			for (size_t i = 1; i < allPaths[j].size(); ++i) {
				excluded.insert(allPaths[j][i]);
				if (i >= 2)
					path.push_back(allPaths[j][i]);
			}
			path[0] += allPaths[j][0];
			
			visitedB.insert(b);
			
			getShortestPathBruteForce(b);
			
			visitedB.erase(b);
			
			path[0] -= allPaths[j][0];

			for (size_t i = 1; i < allPaths[j].size(); ++i) {
				excluded.erase(allPaths[j][i]);
				if (i >= 2)
					path.pop_back();
			}
		}
	}
}

