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

		void setCapacity (const capacity_type &capacity) {
			this->capacity = capacity;
		}

		flow_type getFlow () const {
			return flow;
		}

		void setFlow (const flow_type &flow) {
			this->flow = flow;
		}

		flow_type getRecidualCapacity () const {
			return static_cast<flow_type>(capacity) - flow;
		}

		Edge (const size_type &from, const size_type &to, const capacity_type &capacity) :
			 from(from),
			 to(to),
			 capacity(capacity),
			 flow(static_cast<flow_type>(0)) {}

		template <typename size_type, typename capacity_type, typename flow_type>
		friend ostream& operator << (ostream &os, const Edge<size_type, capacity_type, flow_type> &edge);
};


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
};


// ---------- Dinic algorythm implementation: ----------//


template <class Edge>
typename Edge::flow_type setMaxFlow (Net<Edge> &net) {
	typename Edge::flow_type maxFlow = 0;

	while (true) {
		vector<typename Edge::size_type> layers = buildLayeredNetwork(net);

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

	return maxFlow;
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

		if (layers[edge.getEdge().first] + 1 != layers[edge.getEdge().second] || edge.getRecidualCapacity() <= 0)
			continue;

		const flow_type newPushingFlow = min<flow_type>(pushingFlow, edge.getRecidualCapacity());
		const flow_type pushedFlow = dfs(net, layers, edge.getEdge().second, edgeIterators, newPushingFlow);
		
		if (pushedFlow != 0) {
			arcs[i].first->setFlow(arcs[i].first->getFlow() + pushedFlow);
			arcs[i].second->setFlow(arcs[i].second->getFlow() - pushedFlow);
			return pushedFlow;
		}
	}

	if (edgeIterators[vertex] != arcs.size())
		edgeIterators[vertex] ++;

	return 0;
}


template <class Edge>
typename Edge::flow_type getMinFlowFromSource (const Net<Edge> &net) {
	typename Edge::flow_type ans = numeric_limits<typename Edge::flow_type>::max();

	for (const typename Net<Edge>::Arc &arc: net.getOuterEdges(net.getSource())) {
		if (arc.first->getCapacity() > 0) {
			ans = min<typename Edge::flow_type>(ans, arc.first->getFlow());
		}
	}

	return ans;
}


// ---------- Solving problem: ---------- //



template <class Edge>
void clearFlow (Net<Edge> &net) {
	for (typename Edge::size_type i = 0; i < net.getVerticesCount(); i ++) {
		for (const typename Net<Edge>::Arc &arc: net.getOuterEdges(i)) {
			arc.first->setFlow(0);
		}
	}
}


template <class Edge>
bool canDrink (Net<Edge> &net, const typename Edge::size_type &employesCount, const typename Edge::capacity_type &days) {
	for (const typename Net<Edge>::Arc &arc: net.getOuterEdges(net.getSource())) {
		arc.first->setCapacity(days);
	}

	const typename Edge::flow_type maxFlow = setMaxFlow(net);
	bool out = (maxFlow == employesCount * days);
	clearFlow(net);
	return out;
}


int main () {
	typedef Edge<size_t, size_t, long long> MyEdge;
	typename MyEdge::size_type employesCount, teaTypesCount;
	cin >> employesCount >> teaTypesCount;

	const typename MyEdge::size_type verticesCount = 2 + employesCount + teaTypesCount;
	const typename MyEdge::size_type source = 0;
	const typename MyEdge::size_type sink = verticesCount - 1;
	
	Net<MyEdge> net(verticesCount, source, sink);

	for (typename MyEdge::size_type i = 0; i < employesCount; i ++) {
		const typename MyEdge::capacity_type capacity = 100000000000;
		MyEdge edge(source, 1 + i, capacity);
		net.addEdge(edge);
	}

	for (typename MyEdge::size_type i = 0; i < teaTypesCount; i ++) {
		typename MyEdge::capacity_type capacity;
		cin >> capacity;
		MyEdge edge(1 + employesCount + i, sink, capacity);
		net.addEdge(edge);
	}

	for (typename MyEdge::size_type i = 0; i < employesCount; i ++) {
		typename MyEdge::size_type favouriteTeaTypesCount;
		cin >> favouriteTeaTypesCount;

		for (typename MyEdge::size_type j = 0; j < favouriteTeaTypesCount; j ++) {
			typename MyEdge::size_type favouriteTeaTypeId;
			cin >> favouriteTeaTypeId;
			-- favouriteTeaTypeId;

			const typename MyEdge::capacity_type capacity = 100000000000;
			MyEdge edge(1 + i, 1 + employesCount + favouriteTeaTypeId, capacity);
			net.addEdge(edge);
		}
	}

	typename MyEdge::capacity_type left = 0;
	typename MyEdge::capacity_type right = 100000000000;
	
	while (right - left > 1) {
		const typename MyEdge::capacity_type middle = (left + right) / 2;
		(canDrink(net, employesCount, middle) ? left : right) = middle;
	}

	cout << left << endl;

	return 0;
}
