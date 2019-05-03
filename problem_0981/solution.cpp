#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;


// ---------- DSU class: ---------- //


class DSU {
	private:
		vector<size_t> parent;

	public:
		DSU (const size_t &verticesCount) {
			parent.resize(verticesCount);

			for (size_t i = 0; i < verticesCount; i ++)
				parent[i] = i;
		}

		size_t findRoot (const size_t &vertex) {
			if (vertex == parent[vertex])
				return vertex;

			const size_t root = findRoot(parent[vertex]);
			parent[vertex] = root;
			return root;
		}

		void unite (const size_t &firstVertex, const size_t &secondVertex) {
			const size_t firstRoot = findRoot(firstVertex);
			const size_t secondRoot = findRoot(secondVertex);

			if (firstRoot != secondRoot)
				parent[secondRoot] = firstRoot;
		}
};


// ---------- Kruskal's algorithm implementation: ---------- //


template <typename T>
struct Edge {
	size_t from, to;
	T weight;

	Edge () {}
};


template <typename T>
struct Compare {
	bool operator () (const Edge<T> &firstEdge, const Edge<T> &secondEdge) {
		return firstEdge.weight < secondEdge.weight;
	}
};

template <typename T>
vector<Edge<T> > buildSpanningTree (const size_t &verticesCount, const vector<Edge<T> > &edges) {
	vector<Edge<T> > sortedEdges = edges; // O(E)

	// O(E*log(E)) = O(E*log(V))
	sort(sortedEdges.begin(), sortedEdges.end(), Compare<T>());

	// Create n(verticesCount) disjoint sets with one vertex in each.
	DSU dsu(verticesCount);


	// Main procedure:
	vector<Edge<T> > spanningTree;

	for (const Edge<T> &edge: sortedEdges) {
		if (dsu.findRoot(edge.from) != dsu.findRoot(edge.to)) {
			dsu.unite(edge.from, edge.to);
			spanningTree.push_back(edge);
		}
	}

	return spanningTree;
}

template <typename T>
T countSpanningTreeWeight (const vector<Edge<T> > &edges) {
	T spanningTreeWeight = 0;

	for (const Edge<T> &edge: edges) {
		spanningTreeWeight += edge.weight;
	}

	return spanningTreeWeight;
}


// ---------- Solving problem: ---------- //


int main () {
	// Reading data:
	typedef Edge<size_t> Edge;

	size_t verticesCount, edgesCount;
	cin >> verticesCount >> edgesCount;

	vector<Edge> edges(edgesCount);

	for (Edge &edge: edges) {
		cin >> edge.from >> edge.to >> edge.weight;
		edge.from --;
		edge.to --;
	}


	// Solving problem:
	vector<Edge> spanningTree = buildSpanningTree(verticesCount, edges);
	size_t spanningTreeWeight = countSpanningTreeWeight(spanningTree);

	cout << spanningTreeWeight << endl;

	return 0;
}
