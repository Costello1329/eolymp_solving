#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <vector>
#include <queue>
#include <vector>
#include <set>

using namespace std;


// Net class: --------------------------------------------------


template <class Edge>
struct Net {
	typedef typename Edge::IndexType IndexType;
	typedef Edge* EdgePointer;
	typedef pair<EdgePointer, EdgePointer> Arc;

	typename Edge::IndexType source;
	typename Edge::IndexType sink;
	vector<vector<Arc> > edges;

	explicit Net (const IndexType &verticesCount = 0) {
		edges.resize(verticesCount);
	}

	void addEdge (const Edge &edge) {
		EdgePointer edgePointer = new Edge(edge);
		EdgePointer reversedEdgePointer = new Edge(edge.reversed());
		edges[edgePointer->from].push_back(make_pair(edgePointer, reversedEdgePointer));
		edges[reversedEdgePointer->from].push_back(make_pair(reversedEdgePointer, edgePointer));
	}

	typename Edge::IndexType getVerticesCount () {
		return edges.size();
	}
};


// Min Cost Max Flow algorithm: --------------------------------------------------


template <class Edge>
vector<typename Edge::CostType> findMinCostPathsToAllVerticesFordBellman (const Net<Edge> &net, const typename Edge::IndexType &vertex) {
	vector<typename Edge::CostType> costs(net.edges.size(), numeric_limits<typename Edge::CostType>::max());
	costs[vertex] = 0;
	
	for (typename Edge::IndexType iteration = 0; iteration < net.edges.size() - 1; iteration ++) {
		bool changed = false;

		for (typename Edge::IndexType from = 0; from < net.edges.size(); from ++) {
			for (typename Net<Edge>::Arc arc: net.edges[from]) {
				Edge edge = *arc.first;

				if (edge.residualCapacity() > 0 && costs[edge.from] < numeric_limits<typename Edge::CostType>::max()) {
					changed = true;
					costs[edge.to] = min<typename Edge::CostType>(costs[edge.to], costs[edge.from] + edge.cost);
				}
			}
		}
		
		if (!changed)
			break;		
	}

	return costs;
}

template <class Edge>
vector<typename Edge::CostType> findMinCostPathsToAllVerticesDijkstraSET (const Net<Edge> &net, const vector<typename Edge::CostType> &potentials, const typename Edge::IndexType &vertex) {
	vector<typename Edge::CostType> costs = vector<typename Edge::CostType>(net.edges.size(), numeric_limits<typename Edge::CostType>::max());
	costs[vertex] = 0;

	set<pair<typename Edge::IndexType, typename Edge::CostType> > container;
	container.insert(make_pair(vertex, costs[vertex]));

	while (!container.empty()) {
		typename Edge::IndexType vertex = container.begin()->first;
		container.erase(container.begin());

		for (typename Net<Edge>::Arc arc: net.edges[vertex]) {	
			Edge edge = *arc.first;

			if (edge.residualCapacity() == 0)
				continue;

			if (potentials[edge.from] != numeric_limits<typename Edge::CostType>::max()) {
				edge.cost += potentials[edge.from];
				edge.cost -= potentials[edge.to];
			}
			
			else
				edge.cost = numeric_limits<typename Edge::CostType>::max();

			if (costs[edge.from] + edge.cost < costs[edge.to]) {
				container.erase(make_pair(edge.to, costs[edge.to]));
				costs[edge.to] = costs[edge.from] + edge.cost;
				container.insert(make_pair(edge.to, costs[edge.to]));
			}
		}
	}

	return costs;
}

template <class Edge>
vector<typename Edge::CostType> findMinCostPathsToAllVerticesDijkstra (const Net<Edge> &net, const vector<typename Edge::CostType> &potentials, const typename Edge::IndexType &start) {
	vector<typename Edge::CostType> minDistances = vector<typename Edge::CostType>(net.edges.size(), numeric_limits<typename Edge::CostType>::max());
	vector<typename Edge::IndexType> parents(net.edges.size());
	vector<char> used(net.edges.size(), false);
	minDistances[start] = 0;

	for (typename Edge::IndexType i = 0; i < net.edges.size(); i ++) {
		typename Edge::CostType currentMinDistance = numeric_limits<typename Edge::CostType>::max();
		typename Edge::IndexType minVertex;

		for (typename Edge::IndexType j = 0; j < net.edges.size(); j ++) {
			if (!used[j] && minDistances[j] < currentMinDistance) {
				currentMinDistance = minDistances[j];
				minVertex = j;
			}
		}

		used[minVertex] = true;

		if (currentMinDistance == numeric_limits<typename Edge::CostType>::max())
			break;

		for (typename Net<Edge>::Arc arc: net.edges[minVertex]) {
			Edge edge = *arc.first;

			if (edge.residualCapacity() == 0)
				continue;

			if (potentials[edge.from] != numeric_limits<typename Edge::CostType>::max()) {
				edge.cost += potentials[edge.from];
				edge.cost -= potentials[edge.to];
			}
			
			else
				edge.cost = numeric_limits<typename Edge::CostType>::max();

			if (minDistances[edge.from] + edge.cost < minDistances[edge.to]) {
				minDistances[edge.to] = minDistances[edge.from] + edge.cost;
				parents[edge.to] = edge.from;
			}
		}
	}

	return minDistances;
}

template <class Edge>
bool setFlowDFS (Net<Edge> &net, const vector<typename Edge::CostType> &potentials, const typename Edge::IndexType &vertex, const typename Edge::IndexType &destination, vector<char> &used, typename Edge::FlowType &augmentingFlow) {
	if (used[vertex])
		return false;

	used[vertex] = true;

	if (vertex == destination)
		return true;

	const typename Edge::FlowType temp = augmentingFlow;

	for (typename Net<Edge>::Arc arc: net.edges[vertex]) {
		Edge edge = *arc.first;

		if (edge.residualCapacity() > 0 && edge.cost == potentials[edge.to] - potentials[edge.from]) {
			augmentingFlow = min<typename Edge::FlowType>(augmentingFlow, edge.residualCapacity());

			if (setFlowDFS(net, potentials, edge.to, destination, used, augmentingFlow)) {
				arc.first->flow += augmentingFlow;
				arc.second->flow -= augmentingFlow;
				return true;
			}

			else
				augmentingFlow = temp;
		}
	}

	return false;
}

template <class Edge>
void setMaxFlowMinCost (Net<Edge> &net) {
	vector<typename Edge::CostType> potentials = findMinCostPathsToAllVerticesFordBellman(net, net.source);

	while (true) {
		// Building minCost path:
		vector<char> used(net.edges.size(), false);
		typename Edge::FlowType augmentingFlow = numeric_limits<typename Edge::FlowType>::max();
		if (!setFlowDFS(net, potentials, net.source, net.sink, used, augmentingFlow))
			break;

		// Counting costs:
		vector<typename Edge::CostType> costs = findMinCostPathsToAllVerticesDijkstraSET(net, potentials, net.source);

		// Recouting potentials:
		for (typename Edge::IndexType i = 0; i < net.edges.size(); i ++) {
			if (costs[i] != numeric_limits<typename Edge::CostType>::max())
				potentials[i] += costs[i];	

			else
				potentials[i] = costs[i];
		}
	}
}


// My non-lib code: --------------------------------------------------


// My custom Edge:
template <typename FirstType, typename SecondType, typename ThirdType>
struct CostedFlowedEdge {
	typedef FirstType IndexType;
	typedef SecondType FlowType;
	typedef ThirdType CostType;
	IndexType from;
	IndexType to;
	FlowType flow;
	FlowType capacity;
	CostType cost;

	explicit CostedFlowedEdge (const IndexType &from, const IndexType &to, const FlowType &capacity, const CostType &cost, const unsigned int index = 0) : from(from), to(to), flow(0), capacity(capacity), cost(cost) {}

	FlowType residualCapacity () const {
		return capacity - flow;
	}

	CostedFlowedEdge reversed () const {
		return CostedFlowedEdge(to, from, 0, - cost);
	}

	bool operator == (const CostedFlowedEdge &otherEdge) const {
		return from == otherEdge.from && to == otherEdge.to && flow == otherEdge.flow && capacity == otherEdge.capacity;
	}
};


// Matching:
template <typename Edge>
void inputMatching (Net<Edge> &net, istream &is) {
	typename Edge::IndexType verticesCount;
	is >> verticesCount;

	net = Net<Edge>(2 * verticesCount + 2);
	net.source = 0;
	net.sink = 2 * verticesCount + 1;

	for (typename Edge::IndexType from = 0; from < verticesCount; from ++) {
		net.addEdge(Edge(0, from + 1, 1, 0));
		net.addEdge(Edge(from + 1 + verticesCount, 2 * verticesCount + 1, 1, 0));

		for (typename Edge::IndexType to = 0; to < verticesCount; to ++) {
			typename Edge::CostType cost;
			is >> cost;
			net.addEdge(Edge(from + 1, to + 1 + verticesCount, 1, cost));
		}
	}
}

template <typename Edge>
void solveMatching (Net<Edge> &net, vector<Edge> &matching, typename Edge::CostType &minCost) {
 	setMaxFlowMinCost(net);

 	for (typename Edge::IndexType from = 0; from < net.edges.size(); from ++) {
 		for (typename Net<Edge>::Arc arc: net.edges[from]) {
			Edge edge = *arc.first;

			if (edge.flow > 0)
				minCost += edge.flow * edge.cost;

			if (edge.flow != 0 && edge.from != net.source && edge.to != net.source && edge.from != net.sink && edge.to != net.sink)
				if (edge.flow > 0)
					matching.push_back(edge);
		}
 	}		
}

template <typename Edge>
void outputMatching (const typename Edge::CostType &minCost, const vector<Edge> &matching, ostream &os) {
	os << minCost << endl;

	for (const Edge &edge: matching)
		os << edge.from << " " << edge.to - matching.size() << endl; 
}

template <typename Edge>
void outputMatchingE (const typename Edge::CostType &minCost, const vector<Edge> &matching, ostream &os) {
	os << minCost << endl;

	vector<typename Edge::IndexType> to(matching.size());

	for (const Edge &edge: matching)
		to[edge.from - 1] = edge.to - matching.size();
		
	for (const typename Edge::IndexType &edgeTo: to)
		os << edgeTo << " ";

	os << endl;
}



int main () {
	Net<CostedFlowedEdge<unsigned int, int, int> > net;
	inputMatching(net, cin);
	int minCost = 0;
	vector<CostedFlowedEdge<unsigned int, int, int > > matching;
	solveMatching(net, matching, minCost);
	outputMatchingE(minCost, matching, cout);
	return 0;
}
