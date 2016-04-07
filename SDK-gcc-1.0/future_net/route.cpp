#include "route.h"
#include "lib_record.h"
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <setjmp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <queue>
#include <set>
using namespace std;

#define RANDOM_ON 0
#define REVERSE_ON 0
#define OUTPUT_ON 0

const int MAXN = 610;

jmp_buf jmpManiBuf;
//Global variables
int edges[MAXN][MAXN][2];  //0 ==> weight, 1 ==> index
int edgesCnt[MAXN];
int labels[MAXN][MAXN];

int edgesR[MAXN][MAXN][2];
int edgesRCnt[MAXN];
int distancesToDest[MAXN];

int N = 0;
int src = -1, dest = -1;
int minDistance = INT_MAX;
int tmpEdges[MAXN][MAXN];
int shortestPath[MAXN];
bool color[MAXN];
bool excluded[MAXN]; //indicate a node should not visited
bool visitedB[MAXN]; //indicate a bridge node has been visited
int visitedBCnt; // # of bridge nodes have been visited
//For bridge nodes
int bridges[MAXN];
bool isBridges[MAXN];
int bCnt;
int bEstimate[MAXN];

vector<int> path;
int pathLen = 0;
queue<int> que; //used in spfa

clock_t start_time; //start time of search_route function

//functions
void read_graph(char *topo[5000], int edge_num);
void read_demand(char *demand);
void getShortestPathSPFA(int start);

void getShortestPathBFS(int start);
bool reachableBFS(int start);

void getShortestPathBruteForce(int start);
void output_result();
void SPToDest();
void getBridgesEstimate();

bool first_try = true;
bool second_try = true;

void search_route(char *topo[5000], int edge_num, char *demand)
{
    srand(time(NULL));
    start_time = clock();
    read_graph(topo, edge_num);

    memset(isBridges, 0, sizeof(isBridges));
    bCnt = 0;
    read_demand(demand);

    memset(visitedB, 0, sizeof(visitedB));
    visitedBCnt = 0;

    memset(excluded, 0, sizeof(excluded));
    excluded[src] = true;

    path.clear();
    path.push_back(0);  //distance
    path.push_back(src);  //source


    int jstatus = setjmp(jmpManiBuf);
    if (jstatus == 0) {
        memset(distancesToDest, 0, sizeof(distancesToDest));
        SPToDest();
        if (N <= 120) {
            getShortestPathBruteForce(src);
        } else {
            //memset(distancesToDest, 0, sizeof(distancesToDest));
            //SPToDest();
            //getBridgesEstimate();
            getShortestPathSPFA(src);
        }
    } else if (jstatus == -1) {
        memset(visitedB, 0, sizeof(visitedB));
        visitedBCnt = 0;
        memset(excluded, 0, sizeof(excluded));
        excluded[src] = true;
        path.clear();
        path.push_back(0);  //distance
        path.push_back(src);  //source
        if (N <= 260) {
            getShortestPathSPFA(src);
        } else {
            getShortestPathBFS(src);
        }
    }

    if (pathLen == 0) {
        memset(visitedB, 0, sizeof(visitedB));
        visitedBCnt = 0;
        memset(excluded, 0, sizeof(excluded));
        excluded[src] = true;
        path.clear();
        path.push_back(0);  //distance
        path.push_back(src);  //source
        getShortestPathBruteForce(src);
    }

    output_result();
}

void output_result() {
#if OUTPUT_ON
    cout << minDistance << endl;
#endif
#if REVERSE_ON
    for (int i = pathLen - 1; i >= 1 ; i--)
        record_result(labels[shortestPath[i]][shortestPath[i - 1]]);
#else
    for (int i = 0; i < pathLen - 1; i++)
        record_result(labels[shortestPath[i]][shortestPath[i + 1]]);
#endif
}

void read_graph(char *topo[5000], int edge_num)
{
    memset(tmpEdges, 0, sizeof(tmpEdges));
    for (int i = 0; i < edge_num; i++) {
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
#if REVERSE_ON
                edgesR[i][edgesRCnt[i]][0] = tmpEdges[i][j];
                edgesR[i][edgesRCnt[i]][1] = j;
                edgesRCnt[i]++;

                edges[j][edgesCnt[j]][0] = tmpEdges[i][j];
				edges[j][edgesCnt[j]][1] = i;
				edgesCnt[j]++;
#else
                edges[i][edgesCnt[i]][0] = tmpEdges[i][j];
                edges[i][edgesCnt[i]][1] = j;
                edgesCnt[i]++;

                edgesR[j][edgesRCnt[j]][0] = tmpEdges[i][j];
				edgesR[j][edgesRCnt[j]][1] = i;
				edgesRCnt[j]++;
#endif
            }
}

void read_demand(char *demand)
{
    int i = 0;
    src = 0;
    while(demand[i] != ',') {
        src = src * 10 + demand[i] - '0';
        i++;
    }
    dest = 0;
    i++;
    while(demand[i] != ',') {
        dest = dest * 10 + demand[i] - '0';
        i++;
    }

#if REVERSE_ON
    swap(src, dest);
#endif

    i++;
    while(demand[i] != '\n' && demand[i] != '\0') {
        int tmp = 0;
        while(demand[i] != '\0' && demand[i] != '|' && demand[i] != '\n')
        {
            tmp = tmp * 10 + demand[i] - '0';
            i++;
        }
        bridges[bCnt++] = tmp;
        isBridges[tmp] = true;
        i++;
    }
/*
    printf("src = %d\n", src);
    printf("dest = %d\n", dest);
    set<int>::iterator it;
    for (it = bridges.begin(); it != bridges.end(); it++) {
        printf("%d ", *it);
    }
    printf("\n");
    */
}

void SPToDest() {
    memset(color, 0, sizeof(color));
	que.push(dest);
	color[dest] = true;
	while (!que.empty()) {
		int u = que.front();
		que.pop();
		color[u] = false;
		for (int c = 0; c < edgesRCnt[u]; ++c) {
			int w = edgesR[u][c][0];
			int v = edgesR[u][c][1];
			if (distancesToDest[v] == 0 || distancesToDest[v] > distancesToDest[u] + w) {
				distancesToDest[v] = distancesToDest[u] + w;
				if (color[v] == false) {
					que.push(v);
					color[v] = true;
				}
			}
		}
	}
}

void getBridgesEstimate() {
	for (int i = 0; i < bCnt; i++) {
	    int b = bridges[i];
		int minIn = INT_MAX;
		int minOut = INT_MAX;
		for (int c = 0; c < edgesCnt[b]; ++c) {
			if (minOut > edges[b][c][0])
				minOut = edges[b][c][0];
		}
		for (int c = 0; c < edgesRCnt[b]; ++c) {
			if (minIn > edgesR[b][c][0])
				minIn = edgesR[b][c][0];
		}
		bEstimate[b] = (minIn + minOut) / 2;
	}
}

bool cmp_dis_to_dest(int u, int v) {
    if (distancesToDest[u] < distancesToDest[v]) return true;
    return false;
}

struct BNode {
    int b, dis;
};

bool cmp_dis(const BNode &u, const BNode &v) {
    return u.dis < v.dis;
}

void getShortestPathSPFA(int start) {
    clock_t cur_time = clock();
//    if (first_try && (cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 3000) {
//        first_try = false;
//        longjmp(jmpManiBuf, -1);
//    } else if (second_try && (cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 6000) {
//        second_try = false;
//        longjmp(jmpManiBuf, -1);
//    } else
        if ((cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 9910) {
        longjmp(jmpManiBuf, -2);
    }

    if (!reachableBFS(start)) return;

    //distances[i] == 0 means: 1. start vertex, 2. unreabable
    int *distances = new int[N];
    int *parents = new int[N];
    for (int i = 0; i < N; i++) {
        distances[i] = 0;
        parents[i] = -1;
    }
    set<int> ss;
    ss.insert(start);
    parents[start] = start;
    int tBsize = visitedBCnt;

    while (!ss.empty()) {
        int mindis = INT_MAX, u = -1;
        for (set<int>::iterator it = ss.begin(); it != ss.end(); ++it) {
            if (mindis > distances[*it]) {
                mindis = distances[*it];
                u = *it;
            }
        }
        ss.erase(u);
        if (u == dest) {
            if (tBsize == bCnt) break;
            continue;
        }
        if (u != start && isBridges[u]) {
            if (++tBsize == bCnt) break;
            continue;
        }
        for (int c = 0; c < edgesCnt[u]; ++c) {
            int w = edges[u][c][0];
            int v = edges[u][c][1];
            if (excluded[v]) continue;
            if (distances[v] == 0) {
                distances[v] = distances[u] + w;
                parents[v] = u;
                ss.insert(v);
            } else if (distances[v] > distances[u] + w){
                distances[v] = distances[u] + w;
                parents[v] = u;
            }
        }
    }

    /*
    //here, there may be a path like: start --> b1 --> b2
    while (!que.empty()) {
        int u = que.front();
        que.pop();
        color[u] = false;
        for (int c = 0; c < edgesCnt[u]; ++c) {
            int w = edges[u][c][0];
            int v = edges[u][c][1];
            if (excluded[v]) continue;
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
    */

    //we can be sure that distances[dest] + path.get(0) < minDistance
    //and we update the minDistance and shortestPath
    if (bCnt == visitedBCnt) {
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
    BNode *tmpBridges = new BNode[bCnt - visitedBCnt];
    int tmpBCnt = 0;

    for (int i = 0; i < bCnt; i++) {
        int b = bridges[i];
        if (visitedB[b] || distancesToDest[b] == 0 || distances[b] == 0
            || (path[0] + distances[b] + distancesToDest[b] >= minDistance)
            )
            continue;
        tmpBridges[tmpBCnt].b = b;
        tmpBridges[tmpBCnt].dis = path[0] + distances[b];
        tmpBCnt++;
    }

#if RANDOM_ON
    for (int i = 0; i < tmpBCnt; ++i) {
        int j = rand() % tmpBCnt;
        swap(tmpBridges[i], tmpBridges[j]);
    }
#else
    sort(tmpBridges, tmpBridges + tmpBCnt, cmp_dis);
#endif
    vector<int> res;
    for (int bi = 0; bi < tmpBCnt; ++bi) {
        int b = tmpBridges[bi].b;
#if RANDOM_ON
        if (4 < visitedBCnt && visitedBCnt + 5 < bCnt) {
            int p = rand() % visitedBCnt;
            if (p < visitedBCnt - 3)
                continue;
        }
#endif
        res.clear();
        int p = b;
        while (parents[p] != p) {
            res.push_back(p);
            p = parents[p];
        }
        //res does not contain start & path has already contained it
        path[0] += distances[b];
        for (int i = res.size() - 1; i >= 0; --i) {
            path.push_back(res[i]);
            excluded[res[i]] = true;
        }
        visitedB[b] = true;
        visitedBCnt++;

        getShortestPathSPFA(b);

        visitedB[b] = false;
        visitedBCnt--;

        for (int i = res.size() - 1; i >= 0; --i) {
            path.pop_back();
            excluded[res[i]] = false;
        }
        path[0] -= distances[b];
    }
    delete[] tmpBridges;
    delete[] parents;
    delete[] distances;
}

bool reachableBFS(int start) {
	memset(color, 0, sizeof(color));

	que.push(start);
	color[start] = true;

	while (!que.empty()) {
		int u = que.front();
		que.pop();
		for (int c = 0; c < edgesCnt[u]; ++c) {
			int v = edges[u][c][1];
			if (excluded[v]) continue;
			if (color[v] == false) {
				que.push(v);
				color[v] = true;
			}
		}
	}
	for (int i = 0; i < bCnt; i++) {
	    int b = bridges[i];
		if (!visitedB[b] && color[b] == false)
			return false;
	}
	if (color[dest] == false) return false;
	return true;
}

void getShortestPathBFS(int start) {
    clock_t cur_time = clock();
    if (first_try && (cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 3000) {
        first_try = false;
        longjmp(jmpManiBuf, -1);
    } else if (second_try && (cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 6000) {
        second_try = false;
        longjmp(jmpManiBuf, -1);
    } else if ((cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 9900) {
        longjmp(jmpManiBuf, -2);
    }

    if (!reachableBFS(start)) return;

    //distances[i] == 0 means: 1. start vertex, 2. unreabable
    int *distances = new int[N];
    int *parents = new int[N];
    for (int i = 0; i < N; i++) {
        distances[i] = 0;
        parents[i] = 0;
    }
/*
    if (!que.empty()) {
        cout << "error que" << endl;
        abort();
    }
    */
    que.push(start);
    parents[start] = start;

    //here, there may be a path like: start --> b1 --> b2
    while (!que.empty()) {
        int u = que.front();
        que.pop();
        if (u != start && isBridges[u]) continue;
        for (int c = 0; c < edgesCnt[u]; ++c) {
            int w = edges[u][c][0];
            int v = edges[u][c][1];
            if (excluded[v]) continue;
            if (distances[v] == 0) {
                distances[v] = distances[u] + w;
                parents[v] = u;
                que.push(v);
            }
        }
    }

    //we update the minDistance and shortestPath
    if (bCnt == visitedBCnt) {
        if (distances[dest] + path[0] < minDistance) {
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
        }
        return;
    }

    int estimate = 0;
	for (int i = 0; i < bCnt; i++) {
	    int b = bridges[i];
		if (!visitedB[b])
			estimate += bEstimate[b];
	}

    //here, we can use randomness
    int *tmpBridges = new int[bCnt - visitedBCnt];
    int tmpBCnt = 0;

    for (int i = 0; i < bCnt; i++) {
        int b = bridges[i];
        if (visitedB[b] || distances[b] == 0 ||
                (path[0] + distances[b] + distancesToDest[b] + estimate - bEstimate[b] >= minDistance))
            continue;
        tmpBridges[tmpBCnt++] = b;
    }

#if RANDOM_ON
    for (int i = 0; i < tmpBCnt; ++i) {
        int j = rand() % tmpBCnt;
        swap(tmpBridges[i], tmpBridges[j]);
    }
#else
    sort(tmpBridges, tmpBridges + tmpBCnt, cmp_dis_to_dest);
#endif
    vector<int> res;
    for (int bi = 0; bi < tmpBCnt; ++bi) {
        int b = tmpBridges[bi];
//#if RANDOM_ON
        if (4 < visitedBCnt && visitedBCnt + 5 < bCnt) {
            int p = rand() % visitedBCnt;
            if (p < visitedBCnt - 3)
                continue;
        }
//#endif
        res.clear();
        int p = b;
        while (parents[p] != p) {
            res.push_back(p);
            p = parents[p];
        }
        //res does not contain start & path has already contained it
        path[0] += distances[b];
        for (int i = res.size() - 1; i >= 0; --i) {
            path.push_back(res[i]);
            excluded[res[i]] = true;
        }
        visitedB[b] = true;
        visitedBCnt++;

        getShortestPathBFS(b);

        visitedB[b] = false;
        visitedBCnt--;
        for (int i = res.size() - 1; i >= 0; --i) {
            path.pop_back();
            excluded[res[i]] = false;
        }
        path[0] -= distances[b];
    }
    delete[] tmpBridges;
    delete[] parents;
    delete[] distances;
}

void getAllPaths(int start, int end, vector<int> &curPath, vector<vector<int> > &allPaths) {

    clock_t cur_time = clock();
        if ((cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 8910) {
        longjmp(jmpManiBuf, -2);
    }

    if (path[0] + curPath[0] + distancesToDest[start] >= minDistance)
        return;

    for (int c = 0; c < edgesCnt[start]; ++c) {
        int w = edges[start][c][0];
        int v = edges[start][c][1];
        if (v == end) {
            vector<int> r(curPath);
            r[0] += w;
            r.push_back(v);
            allPaths.push_back(r);
        } else if (!excluded[v] && !isBridges[v]) {
            excluded[v] = true;
            curPath.push_back(v);
            curPath[0] += w;
            getAllPaths(v, end, curPath, allPaths);
            curPath[0] -= w;
            curPath.pop_back();
            excluded[v] = false;
        }
    }
}

void getNextBridges(int start, vector<int> &nextB) {
    memset(color, 0, sizeof(color));
    /*
    if (!que.empty()) {
        cout << "que not empty error" << endl;
        abort();
    }
    */
    que.push(start);
    color[start] = true;
    while (!que.empty()) {
        int u = que.front();
        que.pop();
        if (u != start && isBridges[u]) continue;
        for (int c = 0; c < edgesCnt[u]; ++c) {
            int v = edges[u][c][1];
            if (excluded[v]) continue;
            if (color[v] == false) {
                que.push(v);
                color[v] = true;
                if (isBridges[v] && !visitedB[v])
                    nextB.push_back(v);
            }
        }
    }
}

bool cmp_path_dis(const vector<int> &a, const vector<int> &b) {
    return a[0] < b[0];
}

void getShortestPathBruteForce(int start) {

    clock_t cur_time = clock();
        if ((cur_time - start_time) * 1.0 / CLOCKS_PER_SEC * 1000 > 8910) {
        longjmp(jmpManiBuf, -2);
    }

    if (!reachableBFS(start)) return;

    vector<vector<int> > allPaths;
    //next is destination
    if (visitedBCnt == bCnt) {
    //    cout << "in destination start = "  << start << endl;
        vector<int> curPath;
        curPath.push_back(0);
        curPath.push_back(start);
        getAllPaths(start, dest, curPath, allPaths);
     //   cout << "path_num = " << allPaths.size() << endl;
        if (allPaths.size() == 0) return;
        int mind = INT_MAX, ti = -1;
        //path[0] returns its distance
        for (size_t i = 0; i < allPaths.size(); i++) {
            if (mind > allPaths[i][0]) {
                mind = allPaths[i][0];
                ti = i;
            }
        }
      //  cout << "ti = "  << ti << endl;
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
    vector<int> nextB;
    getNextBridges(start, nextB);

    for (size_t i = 0; i < nextB.size(); i++) {
        int b = nextB[i];
        vector<int> curPath;
        curPath.push_back(0);
        curPath.push_back(start);
        getAllPaths(start, b, curPath, allPaths);
    }

    {
        sort(allPaths.begin(), allPaths.end(), cmp_path_dis);
        for (size_t j = 0; j < allPaths.size(); j++) {
            for (size_t i = 1; i < allPaths[j].size(); ++i) {
                excluded[allPaths[j][i]] = true;
                if (i >= 2)
                    path.push_back(allPaths[j][i]);
            }
            int b = path[path.size() - 1];

            path[0] += allPaths[j][0];
            visitedB[b] = true;
            visitedBCnt++;

            getShortestPathBruteForce(b);

            visitedB[b] = false;
            visitedBCnt--;
            path[0] -= allPaths[j][0];
            for (size_t i = 1; i < allPaths[j].size(); ++i) {
                excluded[allPaths[j][i]] = false;
                if (i >= 2)
                    path.pop_back();
            }
        }
    }
}

