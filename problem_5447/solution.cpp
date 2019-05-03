#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <vector>
#include <queue>
#include <list>
#include <set>

using namespace std;



// Graph class: --------------------------------------------------



template <class Edge>
class Graph {
	private:
		vector<vector<Edge> > edges;

	public:
		Graph () {}

		Graph (const typename Edge::IndexType &verticesCount) {
			edges.clear();
			edges.resize(verticesCount);
		}

		void addEdge (const Edge &edge) {
			edges[edge.from].push_back(edge);
		}

		const vector<Edge>& getConstOuterEdges (const typename Edge::IndexType &vertex) const {
			return edges[vertex];
		}

		typename Edge::IndexType getVerticesCount () const {
			return edges.size();
		}
};



// Net class: --------------------------------------------------



template <class Edge>
class Net {
	private:
		typedef Edge* EdgePointer;
		typedef pair<EdgePointer, EdgePointer> Arc;
		typename Edge::IndexType source;
		typename Edge::IndexType sink;
		vector<Arc>* edges;
		vector<vector<Arc>*> innerEdges;
		vector<vector<Arc>*> outerEdges;

		class UniversalEdgeIterator {
			friend class Net;

			private:
				vector<Arc>* edgesList;
				typename vector<Arc>::iterator edgeIterator;

				bool validated;
				enum UniversalIteratorType {inner, outer, all};
				UniversalIteratorType iteratorType;
	
				explicit UniversalEdgeIterator (vector<Arc> *edgesList, const UniversalIteratorType &iteratorType) : edgesList(edgesList), iteratorType(iteratorType), validated(false)  {};

				void validate () {
					if (!validated)
						begin();

					validated = true;
				}

				void jump () {
					while (edgeIterator != edgesList->end() && (edgeIterator->first->deleted || edgeIterator->second->deleted))
						++ edgeIterator;
				}

			public:
				typedef UniversalIteratorType IteratorType;

				explicit UniversalEdgeIterator () : edgesList(nullptr), validated(false) {}

				// Iterator methods:
				void begin () {
					edgeIterator = edgesList->begin();
					jump();
				}

				void move () {
					edgeIterator = next(edgeIterator);
					jump();
				}
	
				bool end () {
					jump();
					return edgesList->size() == 0 || edgeIterator == edgesList->end();
				}

				// Getters:
				Edge getEdge () const {
					return *(edgeIterator->first);
				}

				Edge getReversedEdge () const {
					return *(edgeIterator->second);
				}

				// Setter:
				void increaseFlow (const typename Edge::FlowType &delta) const {
					edgeIterator->first->flow += delta;
					edgeIterator->second->flow -= delta;
				}

				void deleteEdge () {
					edgeIterator->first->deleted = true;
					edgeIterator->second->deleted = true;
				}

				// Other:
				IteratorType getIteratorType () const {
					return iteratorType;
				}
		};

		UniversalEdgeIterator* allEdgesIterator;
		vector<UniversalEdgeIterator*> innerEdgesIterators;
		vector<UniversalEdgeIterator*> outerEdgesIterators;

		void expandIndexation (const typename Edge::IndexType &newVerticesCount) {
			const typename Edge::IndexType &oldVerticesCount = getVerticesCount();

			if (oldVerticesCount >= newVerticesCount)
				return;

			innerEdges.resize(newVerticesCount);
			outerEdges.resize(newVerticesCount);
			innerEdgesIterators.resize(newVerticesCount);
			outerEdgesIterators.resize(newVerticesCount);

			for (typename Edge::IndexType i = oldVerticesCount; i < newVerticesCount; i ++) {
				innerEdges[i] = new vector<Arc>();
				outerEdges[i] = new vector<Arc>();
				innerEdgesIterators[i] = new EdgeIterator(&(*innerEdges[i]), EdgeIterator::IteratorType::inner);
				outerEdgesIterators[i] = new EdgeIterator(&(*outerEdges[i]), EdgeIterator::IteratorType::outer);
			}
		}

	public:
		typedef UniversalEdgeIterator EdgeIterator;
		typedef typename Edge::IndexType IndexType;

		// Constructors, destructors:
		explicit Net (const IndexType &verticesCount = 0) {
			edges = new vector<Arc>();
			allEdgesIterator = new EdgeIterator(&(*edges), EdgeIterator::IteratorType::all);
		}

		~Net () {
			for (const Arc &arc: *edges)
				delete arc.first;

			delete edges;
			delete allEdgesIterator;

			for (typename Edge::IndexType vertex = 0; vertex < getVerticesCount(); vertex ++) {
				delete innerEdges[vertex];
				delete outerEdges[vertex];
				delete innerEdgesIterators[vertex];
				delete outerEdgesIterators[vertex];
			}
		}

		// Setters:
		void setSource (const typename Edge::IndexType &source) {
			expandIndexation(source + 1);
			this->source = source;
		}

		void setSink (const typename Edge::IndexType &sink) {
			expandIndexation(sink + 1);
			this->sink = sink;
		}

		void addEdge (const Edge &edge) {
			expandIndexation(max<typename Edge::IndexType>(edge.from, edge.to) + 1);

			EdgePointer edgePointer = new Edge(edge);
			EdgePointer reversedEdgePointer = new Edge(edge.reversed());

			edges->push_back(make_pair(edgePointer, reversedEdgePointer));
			outerEdges[edgePointer->from]->push_back(make_pair(edgePointer, reversedEdgePointer));
			innerEdges[edgePointer->to]->push_back(make_pair(edgePointer, reversedEdgePointer));

			edges->push_back(make_pair(reversedEdgePointer, edgePointer));
			outerEdges[reversedEdgePointer->from]->push_back(make_pair(reversedEdgePointer, edgePointer));
			innerEdges[reversedEdgePointer->to]->push_back(make_pair(reversedEdgePointer, edgePointer));

			allEdgesIterator->validate();
			outerEdgesIterators[edgePointer->from]->validate();
			outerEdgesIterators[edgePointer->to]->validate();
			innerEdgesIterators[edgePointer->from]->validate();
			innerEdgesIterators[edgePointer->to]->validate();
		}

		// Getters:
		const typename Edge::IndexType &getSource () const {
			return source;
		}

		const typename Edge::IndexType &getSink () const {
			return sink;
		}

		EdgeIterator& getAllEdgesIterator () {
			return *allEdgesIterator;
		}

		EdgeIterator getConstAllEdgesIterator () const {
			return *allEdgesIterator;
		}

		EdgeIterator& getInnerEdgeIterator (const IndexType &vertex) {
			return *innerEdgesIterators[vertex];
		}

		EdgeIterator getConstInnerEdgeIterator (const IndexType &vertex) const {
			return *innerEdgesIterators[vertex];
		}

		EdgeIterator& getOuterEdgeIterator (const IndexType &vertex) {
			return *outerEdgesIterators[vertex];
		}

		EdgeIterator getConstOuterEdgeIterator (const IndexType &vertex) const {
			return *outerEdgesIterators[vertex];
		}

		// Info:
		IndexType getVerticesCount () const {
			return innerEdges.size();
		}

		IndexType getEdgesCount () const {
			return edges->size();
		}
};


// Min Cost Max Flow algorithm: --------------------------------------------------


template <class Edge>
vector<typename Edge::CostType> findMinCostPathsToAllVerticesFordBellman (const Graph<Edge> &graph, const typename Edge::IndexType &vertex) {
	vector<typename Edge::CostType> costs(graph.getVerticesCount(), numeric_limits<typename Edge::CostType>::max());
	costs[vertex] = 0;
	
	for (typename Edge::IndexType iteration = 0; iteration < graph.getVerticesCount() - 1; iteration ++)
		for (typename Edge::IndexType from = 0; from < graph.getVerticesCount(); from ++)
			for (const Edge &edge: graph.getConstOuterEdges(from))
				if (costs[edge.from] < numeric_limits<typename Edge::CostType>::max())
					costs[edge.to] = min<typename Edge::CostType>(costs[edge.to], costs[edge.from] + edge.cost);

	return costs;
}

template <class Edge>
vector<typename Edge::CostType> findMinCostPathsToAllVerticesDijkstra1 (const Graph<Edge> &graph, const typename Edge::IndexType &vertex) {
	vector<typename Edge::CostType> costs = vector<typename Edge::CostType>(graph.getVerticesCount(), numeric_limits<typename Edge::CostType>::max());
	costs[vertex] = 0;

	set<pair<typename Edge::IndexType, typename Edge::CostType> > container;
	container.insert(make_pair(vertex, costs[vertex]));

	while (!container.empty()) {
		typename Edge::IndexType vertex = container.begin()->first;
		container.erase(container.begin());

		for (const Edge &edge: graph.getConstOuterEdges(vertex)) {	
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
vector<typename Edge::CostType> findMinCostPathsToAllVerticesDijkstra2 (const Graph<Edge> &graph, const typename Edge::IndexType &start) {
	vector<typename Edge::CostType> minDistances = vector<typename Edge::CostType>(graph.getVerticesCount(), numeric_limits<typename Edge::CostType>::max());
	vector<typename Edge::IndexType> parents(graph.getVerticesCount());
	vector<bool> used(graph.getVerticesCount(), false);
	minDistances[start] = 0;

	for (typename Edge::IndexType i = 0; i < graph.getVerticesCount(); i ++) {
		typename Edge::CostType currentMinDistance = numeric_limits<typename Edge::CostType>::max();
		typename Edge::IndexType minVertex;

		for (typename Edge::IndexType j = 0; j < graph.getVerticesCount(); j ++) {
			if (!used[j] && minDistances[j] < currentMinDistance) {
				currentMinDistance = minDistances[j];
				minVertex = j;
			}
		}

		used[minVertex] = true;

		if (currentMinDistance == numeric_limits<typename Edge::CostType>::max())
			break;

		for (const Edge &edge: graph.getConstOuterEdges(minVertex)) {	
			if (minDistances[edge.from] + edge.cost < minDistances[edge.to]) {
				minDistances[edge.to] = minDistances[edge.from] + edge.cost;
				parents[edge.to] = edge.from;
			}
		}
	}

	return minDistances;
}

template <class Edge>
bool buildMinCostPathDFS (const Graph<Edge> &potentialGraph, const typename Edge::IndexType &vertex, const typename Edge::IndexType &destination, vector<bool> &used, list<Edge> &path) {
	if (used[vertex])
		return false;

	used[vertex] = true;

	if (vertex == destination)
		return true;

	for (const Edge &edge: potentialGraph.getConstOuterEdges(vertex)) {
		if (edge.cost == 0) {
			path.push_back(edge);

			if (buildMinCostPathDFS(potentialGraph, edge.to, destination, used, path))
				return true;

			else
				path.pop_back();
		}
	}

	return false;
}

template <class Edge>
bool setFlowOnPathDFS (Net<Edge> &net, const typename Edge::IndexType &vertex, const typename Edge::FlowType &augmentingFlow, list<Edge> &path, vector<bool> &used) {
	if (used[vertex])
		return false;

	used[vertex] = true;

	if (path.empty())
		return true;

	typename Net<Edge>::EdgeIterator iter = net.getConstOuterEdgeIterator(vertex);
	for (iter.begin(); !iter.end(); iter.move()) {
		if (iter.getEdge().residualCapacity() <= 0)
			continue;

		if (iter.getEdge() == path.front()) {
			path.pop_front();
			iter.increaseFlow(augmentingFlow);

			if (setFlowOnPathDFS(net, iter.getEdge().to, augmentingFlow, path, used))
				return true;

			else {
				path.push_front(iter.getEdge());
				iter.increaseFlow(- augmentingFlow);
			}
		}
	}

	return false;
}


template <class Edge>
Graph<Edge> buildPotentialGraph (const Net<Edge> &net, const vector<typename Edge::CostType> &potentials) {
	Graph<Edge> potentialGraph(net.getVerticesCount());

	typename Net<Edge>::EdgeIterator iter = net.getConstAllEdgesIterator();
	for (iter.begin(); !iter.end(); iter.move()) {
		if (iter.getEdge().residualCapacity() == 0)
			continue;

		Edge newEdge = iter.getEdge();
		
		if (potentials[iter.getEdge().from] != numeric_limits<typename Edge::CostType>::max())
			newEdge.cost = iter.getEdge().cost + potentials[iter.getEdge().from] - potentials[iter.getEdge().to];
		
		else
			newEdge.cost = numeric_limits<typename Edge::CostType>::max();
		
		potentialGraph.addEdge(newEdge);
	}

	return potentialGraph;
}

template <class Edge>
void setMaxFlowMinCost (Net<Edge> &net) {
	typename Net<Edge>::EdgeIterator iter = net.getConstAllEdgesIterator();
	for (iter.begin(); !iter.end(); iter.move())
		if (iter.getEdge().flow != 0)
			iter.increaseFlow(- iter.getEdge().flow);


	Graph<Edge> potentialGraph = buildPotentialGraph(net, vector<typename Edge::CostType>(net.getVerticesCount(), 0));
	vector<typename Edge::CostType> potentials = findMinCostPathsToAllVerticesFordBellman(potentialGraph, net.getSource());

	while (true) {
		// Setting up potential graph to work with:
		potentialGraph = buildPotentialGraph(net, potentials);

		// Building minCost path:
		vector<bool> used(net.getVerticesCount(), false);
		list<Edge> path;
		
		if (!buildMinCostPathDFS(potentialGraph, net.getSource(), net.getSink(), used, path))
			break;

		// Setting up flow on minCost path:
		typename Edge::FlowType augmentingFlow = numeric_limits<typename Edge::FlowType>::max();
		
		for (const Edge &edge: path)
			augmentingFlow = min<typename Edge::FlowType>(augmentingFlow, edge.residualCapacity());

		used.assign(net.getVerticesCount(), false);
		setFlowOnPathDFS(net, net.getSource(), augmentingFlow, path, used);

		// Rebuilding graph:
		potentialGraph = buildPotentialGraph(net, potentials);

		// Counting costs:
		vector<typename Edge::CostType> costs = findMinCostPathsToAllVerticesDijkstra2(potentialGraph, net.getSource());

		// Recouting potentials:
		for (typename Edge::IndexType i = 0; i < net.getVerticesCount(); i ++) {
			if (costs[i] != numeric_limits<typename Edge::CostType>::max())
				potentials[i] += costs[i];	

			else
				potentials[i] = costs[i];
		}
	}
}




// My non-lib code: --------------------------------------------------




// My custom Edge:
template <typename FirstType = unsigned long long, typename SecondType = long long, typename ThirdType = long long>
struct CostedFlowedEdge {
	typedef FirstType IndexType;
	typedef SecondType FlowType;
	typedef ThirdType CostType;
	IndexType from;
	IndexType to;
	FlowType flow;
	FlowType capacity;
	CostType cost;
	bool deleted = false;

	explicit CostedFlowedEdge (const IndexType &from, const IndexType &to, const FlowType &capacity, const CostType &cost, const unsigned long long index = 0) : from(from), to(to), flow(0), capacity(capacity), cost(cost) {}

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


// MinCostMaxFlow:
template <typename Edge>
Net<Edge> inputMinCostMaxFlow (istream &is) {
	unsigned verticesCount, edgesCount;
	is >> verticesCount >> edgesCount;

	Net<Edge> net(verticesCount);
	net.setSource(0);
	net.setSink(verticesCount - 1);

	for (unsigned edge = 0; edge < edgesCount; edge ++) {
		typename Edge::IndexType from, to;
		typename Edge::FlowType capacity;
		typename Edge::CostType cost;
		is >> from >> to >> capacity >> cost;
		net.addEdge(Edge(from - 1, to - 1, capacity, cost));
	}

	return net;
}

template <typename Edge>
typename Edge::CostType solveMinCostMaxFlow (Net<Edge> &net) {
	setMaxFlowMinCost(net);
	typename Edge::CostType minCost = 0;

	typename Net<Edge>::EdgeIterator iter = net.getConstAllEdgesIterator();
	for (iter.begin(); !iter.end(); iter.move())
		if (iter.getEdge().flow > 0)
			minCost += iter.getEdge().flow * iter.getEdge().cost;

	return minCost;
}

void outputMinCostMaxFlow (const long long &minCost, ostream &os) {
	os << minCost << endl;
}

int main () {
	Net<CostedFlowedEdge<unsigned long long, long long, long long> > net = inputMinCostMaxFlow<CostedFlowedEdge<unsigned long long, long long, long long> >(cin);
	long long minCost = solveMinCostMaxFlow(net);
	outputMinCostMaxFlow(minCost, cout);
	
	return 0;
}