#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <queue>

using namespace std;


// ---------- Graph class: ---------- //


template <typename _size_type, typename _capacity_type, typename _flow_type>
class Edge {
	public:
		typedef _size_type size_type;
		typedef _capacity_type capacity_type;
		typedef _flow_type flow_type;

	private:
		_size_type from, to;
		_capacity_type capacity;
		_flow_type flow;

	public:
		pair<size_type, size_type> getEdge () const {
			return make_pair(from, to);
		}

		capacity_type getCapacity () const {
			return capacity;
		}

		flow_type getRecidualCapacity () const {
			return static_cast<flow_type>(capacity) - flow;
		}

		void increseFlow (const flow_type &increment) {
			flow += increment;
		}

		Edge (const size_type &from, const size_type &to, const capacity_type &capacity) :
			 from(from),
			 to(to),
			 capacity(capacity),
			 flow(static_cast<flow_type>(0)) {}

		template <typename size_type, typename capacity_type, typename flow_type>
		friend ostream& operator << (ostream &os, const Edge<size_type, capacity_type, flow_type> &edge);
};

template <typename size_type, typename capacity_type, typename flow_type>
ostream& operator << (ostream &os, const Edge<size_type, capacity_type, flow_type> &edge) {
	return os << edge.from << " -> " << edge.to << " " << edge.flow << "/" << edge.capacity;
}


template <class Edge>
class Net {
	public:
		typedef typename Edge::size_type size_type;
		typedef pair<Edge*, Edge*> Arc;

	private:
		size_type verticesCount;
		size_type source;
		size_type sink;
		vector<vector<Arc> > edges;

	public:
		Net (const size_type &verticesCount, const size_type &source, const size_type &sink) :
			verticesCount(verticesCount),
			source(source),
			sink(sink) {
			edges.resize(verticesCount);
		}

		// Setters:
		void addEdge (const Edge &edge) {
			Edge *straightEdge = new Edge(edge.getEdge().first, edge.getEdge().second, edge.getCapacity());
			Edge *reverseEdge = new Edge(edge.getEdge().second, edge.getEdge().first, 0);

			Arc straightArc = make_pair(straightEdge, reverseEdge);
			Arc reverseArc = make_pair(reverseEdge, straightEdge);

			edges[edge.getEdge().first].push_back(straightArc);
			edges[edge.getEdge().second].push_back(reverseArc);
		}

		// Getters:
		const vector<Arc>& getOuterEdges (const size_type &from) const {
			return edges[from];
		}

		size_type getVerticesCount () const {
			return verticesCount;
		}

		size_type getSource () const {
			return source;
		}

		size_type getSink () const {
			return sink;
		}

		// DEBUG:
		void print (ostream &os) {
			for (size_t i = 0; i < verticesCount; i ++) {
				for (size_t j = 0; j < edges[i].size(); j ++) {
					os << *(edges[i][j].first) << "; ";
				}

				os << endl;
			}
		}
};


// ---------- Edmonds-Karp algorythm implementation: ----------//


template <class Edge>
void setMaxFlow (Net<Edge> &net) {
	typename Edge::flow_type maxFlow = 0;

	while (true) {
		//cout << "Net :" << endl;
		//net.print(cout);

		vector<typename Edge::size_type> layers = buildLayeredNetwork(net);

		//cout << "Layers: ";

		//for (const typename Edge::size_type &layer: layers)
		//	cout << layer << " ";
	
		//cout << endl << endl;

		if (layers[net.getSink()] == static_cast<typename Edge::size_type>(-1))
			break;


		typename Edge::flow_type augmentation = 0;
		vector<typename Edge::size_type> edgeIterators(net.getVerticesCount(), 0);

		do {
			augmentation = dfs(net, layers, net.getSource(), edgeIterators);
			maxFlow += augmentation;
		}
	
		while (augmentation > 0);
	}

	cout << maxFlow << endl;
}

template <class Edge>
vector<typename Edge::size_type> buildLayeredNetwork (const Net<Edge> &net) {
	typedef typename Edge::size_type size_type;
	queue<size_type> que;
	vector<size_type> layers(net.getVerticesCount(), static_cast<size_type>(-1));
	layers[net.getSource()] = 0;
	vector<bool> used(net.getVerticesCount(), false);
	que.push(net.getSource());

	while (!que.empty()) {
		const size_type vertex = que.front();
		que.pop();

		if (used[vertex])
			continue;

		else {
			for (const typename Net<Edge>::Arc &arc: net.getOuterEdges(vertex)) {
				if (arc.first->getRecidualCapacity() <= 0)
					continue;

				const size_type nextVertex = arc.first->getEdge().second;

				if (layers[nextVertex] == static_cast<size_type>(-1)) {
					layers[nextVertex] = layers[vertex] + 1;
					que.push(nextVertex);
				}
			}
		}
		
		used[vertex] = true;
	}

	return layers;
}

template <class Edge>
typename Edge::flow_type dfs (Net<Edge> &net,
							  const vector<typename Edge::size_type> &layers,
							  const typename Edge::size_type &vertex,
							  vector<typename Edge::size_type> &edgeIterators,
							  const typename Edge::flow_type &pushingFlow = numeric_limits<typename Edge::flow_type>::max()) {
	typedef typename Edge::flow_type flow_type;

	if (vertex == net.getSink() || pushingFlow == 0)
		return pushingFlow;

	const vector<typename Net<Edge>::Arc> &arcs = net.getOuterEdges(vertex);

	for (size_t i = edgeIterators[vertex]; i < arcs.size(); i ++) {
		const Edge edge = *arcs[i].first;

		// if edge is not in residual network then continue:
		if (layers[edge.getEdge().first] + 1 != layers[edge.getEdge().second] || edge.getRecidualCapacity() <= 0)
			continue;

		const flow_type newPushingFlow = min<flow_type>(pushingFlow, edge.getRecidualCapacity());
		const flow_type pushedFlow = dfs(net, layers, edge.getEdge().second, edgeIterators, newPushingFlow);
		
		if (pushedFlow != 0) {
			arcs[i].first->increseFlow(+ pushedFlow);
			arcs[i].second->increseFlow(- pushedFlow);
			return pushedFlow;
		}
	}

	if (edgeIterators[vertex] != arcs.size())
		while (++ edgeIterators[vertex] != arcs.size()
			   && (
			   		arcs[edgeIterators[vertex]].first->getRecidualCapacity() <= 0
			   		|| (
			   			layers[arcs[edgeIterators[vertex]].first->getEdge().first] + 1
			   			!= layers[arcs[edgeIterators[vertex]].first->getEdge().second]
			   		)
			   	   )
			   );
	
	return 0;
}

/*
template <class Edge>
typename Edge::flow_type dfs (Net<Edge> &net, const typename Edge::IndexType &vertex, const typename Edge::FlowType &pushingFlow = numeric_limits<typename Edge::FlowType>::max()) {
	if (vertex == layeredNet.getSink() || pushingFlow == 0)
		return pushingFlow;

	for (typename Net<Edge>::EdgeIterator iter = layeredNet.getConstOuterEdgeIterator(vertex); !iter.end(); iter.move()) {
		Edge edge = iter.getEdge();

		if (edge.capacity == 0)
			continue;

		typename Edge::FlowType newPushingFlow = min<typename Edge::FlowType>(pushingFlow, edge.residualCapacity());
		typename Edge::FlowType pushedFlow = tryPushFlow(layeredNet, edge.to, newPushingFlow);

		if (pushedFlow != 0) {
			iter.increaseFlow(pushedFlow);
			return pushedFlow;
		}
	}

	if (!layeredNet.getOuterEdgeIterator(vertex).end())
		layeredNet.getOuterEdgeIterator(vertex).move();
	
	return 0;
}
*/


// ---------- Solving problem: ---------- //


int main () {
	typedef Edge<size_t, size_t, long long> MyEdge;
	size_t verticesCount;
	size_t edgesCount;

	cin >> verticesCount >> edgesCount;

	const MyEdge::size_type source = 0;
	const MyEdge::size_type sink = verticesCount - 1;
	Net<MyEdge> net(verticesCount, source, sink);

	for (size_t i = 0; i < edgesCount; i ++) {
		MyEdge::size_type from, to;
		MyEdge::capacity_type capacity;
		cin >> from >> to >> capacity;
		-- from;
		-- to;
		MyEdge edge(from, to, capacity);
		net.addEdge(edge);
		MyEdge reversedEdge(to, from, capacity);
		net.addEdge(reversedEdge);
	}


	setMaxFlow(net);


	return 0;
}
