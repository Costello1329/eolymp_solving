#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <limits>
#include <vector>
#include <queue>
#include <list>

using namespace std;


template <class Edge>
class Net {
	private:
		typedef Edge* EdgePointer;
		typedef pair<EdgePointer, EdgePointer> Arc;
		typename Edge::IndexType source;
		typename Edge::IndexType sink;
		list<Arc>* edges;
		vector<list<Arc>*> innerEdges;
		vector<list<Arc>*> outerEdges;

		class UniversalEdgeIterator {
			friend class Net;

			private:
				list<Arc>* edgesList;
				typename list<Arc>::iterator edgeIterator;

				bool validated;
				enum UniversalIteratorType {inner, outer, all};
				UniversalIteratorType iteratorType;
	
				UniversalEdgeIterator (list<Arc> *edgesList, const UniversalIteratorType &iteratorType) {
					this->edgesList = edgesList;
					this->iteratorType = iteratorType;
					validated = false;
				}

				void validate () {
					if (!validated)
						begin();

					validated = true;
				}

			public:
				typedef UniversalIteratorType IteratorType;

				UniversalEdgeIterator () {
					edgesList = nullptr;
					validated = false;
				}

			
				void begin () {
					edgeIterator = edgesList->begin();
				}

				void move () {
					edgeIterator = next(edgeIterator);
				}
	
				bool end () const {
					return edgesList->size() == 0 || edgeIterator == edgesList->end();
				}

			
				Edge getEdge () const {
					return *(edgeIterator->first);
				}

				Edge getReversedEdge () const {
					return *(edgeIterator->second);
				}

			
				void increaseFlow (const typename Edge::FlowType &delta) const {
					edgeIterator->first->flow += delta;
					edgeIterator->second->flow -= delta;
				}

			
				IteratorType getIteratorType () const {
					return iteratorType;
				}
		};

		UniversalEdgeIterator* allEdgesIterator;
		vector<UniversalEdgeIterator*> innerEdgesIterators;
		vector<UniversalEdgeIterator*> outerEdgesIterators;

		void expand (const typename Edge::IndexType &newVerticesCount) {
			const typename Edge::IndexType &oldVerticesCount = getVerticesCount();

			innerEdges.resize(newVerticesCount);
			outerEdges.resize(newVerticesCount);
			innerEdgesIterators.resize(newVerticesCount);
			outerEdgesIterators.resize(newVerticesCount);

			for (unsigned i = oldVerticesCount; i < newVerticesCount; i ++) {
				innerEdges[i] = new list<Arc>();
				outerEdges[i] = new list<Arc>();
				innerEdgesIterators[i] = new EdgeIterator(&(*innerEdges[i]), EdgeIterator::IteratorType::inner);
				outerEdgesIterators[i] = new EdgeIterator(&(*outerEdges[i]), EdgeIterator::IteratorType::outer);
			}
		}

	public:
		typedef UniversalEdgeIterator EdgeIterator;
		typedef typename Edge::IndexType IndexType;

	
		Net (const IndexType &verticesCount = 0) {
			edges = new list<Arc>();
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

	
		void setSource (const typename Edge::IndexType &source) {
			if (this->source >= getVerticesCount())
				expand(source + 1);

			this->source = source;
		}

		void setSink (const typename Edge::IndexType &sink) {
			if (this->sink >= getVerticesCount())
				expand(sink + 1);

			this->sink = sink;
		}

		void addEdge (const Edge &edge) {
			if (edge.from >= getVerticesCount() || edge.to >= getVerticesCount())
				expand(max<typename Edge::IndexType>(edge.from, edge.to) + 1);

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

	
		IndexType getVerticesCount () const {
			return innerEdges.size();
		}

		IndexType getEdgesCount () const {
			return edges->size();
		}
};


template <class Edge>
const typename Edge::FlowType getFlowFromOrToVertex (typename Net<Edge>::EdgeIterator iter) {
	typename Edge::FlowType flow = 0;

	for (iter.begin(); !iter.end(); iter.move())
		flow += iter.getEdge().flow;

	return flow;
}

template <class Edge>
bool isFlowValid (const Net<Edge> &net) {
	typename Net<Edge>::EdgeIterator iter = net.getConstAllEdgesIterator();

	for (iter.begin(); !iter.end(); iter.move())
		if (iter.getEdge().flow > iter.getEdge().capacity || iter.getEdge().flow != -iter.getReversedEdge().flow)
			return false;

	for (typename Edge::IndexType vertex = 0; vertex < net.getVerticesCount(); vertex ++) {
		if (vertex == net.getSource() || vertex == net.getSink())
			break;

		typename Edge::FlowType innerFlow = getFlowFromOrToVertex<Edge>(net.getConstInnerEdgeIterator(vertex));
		typename Edge::FlowType outerFlow = getFlowFromOrToVertex<Edge>(net.getConstOuterEdgeIterator(vertex));

		if (innerFlow != outerFlow)
			return false;
	}

	return true;
}

template <class Edge>
typename Edge::FlowType getFlowWithoutCheck (const Net<Edge> &net) {
	return getFlowFromOrToVertex<Edge>(net.getConstOuterEdgeIterator(net.getSource()));
}

template <class Edge>
typename Edge::FlowType getFlowWithCheck (const Net<Edge> &net) {
	if (isFlowValid(net))
		return getFlowWithoutCheck(net);

	else
		throw ("Net exeption: invalid flow!");
}


// Global random preflow-push operations O(V^2 E):
template <class Edge>
inline bool canPush (const Edge &edge, const vector<typename Edge::IndexType> &height) {
	if (height[edge.from] == height[edge.to] + 1 && edge.residualCapacity() > 0)
		return true;

	else
		return false;
}

template <class Edge>
inline bool canLift (const Net<Edge> &net, const typename Edge::IndexType &vertex, const vector<typename Edge::IndexType> &height) {
	typename Net<Edge>::EdgeIterator iter = net.getConstOuterEdgeIterator(vertex);

	for (iter.begin(); !iter.end(); iter.move()) {
		Edge edge = iter.getEdge();

		if (edge.residualCapacity() > 0 & height[edge.from] > height[edge.to])
			return true;
	}

	return false;
}

template <class Edge>
void randomPreflowPushOperations (Net<Edge> &net, vector<typename Edge::FlowType> &excess, vector<typename Edge::IndexType> &height) {
	while (true) {
		bool operated = false;
		typename Net<Edge>::EdgeIterator &iter = net.getAllEdgesIterator;

		for (iter.begin(); !iter.end(); iter.move()) {
			if (canPush(iter.getEdge(), height)) {
				cout << "ge" << endl;
				const typename Edge::FlowType pushingFlow = min<typename Edge::FlowType>(iter.getEdge().residualCapacity(), excess[iter.getEdge().from]);
				push(net, iter, pushingFlow, excess);
				operated = true;
			}
		}

		for (typename Edge::IndexType vertex = 0; vertex < net.getVerticesCount(); vertex ++) {
			if (canLift(net, vertex, height)) {
				lift(net, vertex, height);
				operated = true;
			}
		}

		if (!operated)
			break;
	}
}

// Global discharge preflow-push operations O(V^3):
template <class Edge>
inline void discharge (Net<Edge> &net, const typename Edge::IndexType &vertex, vector<typename Edge::FlowType> &excess, vector<typename Edge::IndexType> &height) {
	while (excess[vertex] > 0) {
		typename Net<Edge>::EdgeIterator &iter = net.getOuterEdgeIterator(vertex);

		if (iter.end()) {
			lift(net, vertex, height);
			iter.begin();
		}

		else {
			Edge edge = iter.getEdge();

			if (edge.residualCapacity() > 0 && height[vertex] == height[edge.to] + 1) {
				typename Edge::FlowType pushingFlow = min<typename Edge::FlowType>(excess[edge.from], edge.residualCapacity());
				push(net, iter, pushingFlow, excess);
			}

			else
				iter.move();
		}
	}
}

template <class Edge>
void globalDischargePreflowPushOperations (Net<Edge> &net, vector<typename Edge::FlowType> &excess, vector<typename Edge::IndexType> &height) {
	while (true) {
		bool discharged = false;

		for (typename Edge::IndexType vertex = 0; vertex < net.getVerticesCount(); vertex ++) {
			if (excess[vertex] && vertex != net.getSource() && vertex != net.getSink()) {
				discharge(net, vertex, excess, height);
				discharged = true;
			}
		}

		if (!discharged)
			break;
	}
}

// Global Preflow-Push algorithm:
template <class Edge>
void push (const Net<Edge> &net, const typename Net<Edge>::EdgeIterator &iter, const typename Edge::FlowType &pushingFlow, vector<typename Edge::FlowType> &excess) {
	iter.increaseFlow(pushingFlow);
	Edge edge = iter.getEdge();

	if (edge.from != net.getSource() && edge.from != net.getSink())
		excess[edge.from] -= pushingFlow;

	if (edge.to != net.getSource() && edge.to != net.getSink())
		excess[edge.to] += pushingFlow;
}

template <class Edge>
void lift (const Net<Edge> &net, const typename Edge::IndexType &vertex, vector<typename Edge::IndexType> &height) {
	typename Edge::IndexType minHeight = numeric_limits<typename Edge::IndexType>::max();

	typename Net<Edge>::EdgeIterator iter = net.getConstOuterEdgeIterator(vertex);

	for (iter.begin(); !iter.end(); iter.move()) {
		Edge edge = iter.getEdge();

		if (edge.residualCapacity() > 0)
			minHeight = min<typename Edge::IndexType>(minHeight, height[edge.to]);
	}

	height[vertex] = minHeight + 1;
}

template <class Edge>
void setMaxFlowWithPreflowPush (Net<Edge> &net) {
	vector<typename Edge::FlowType> excess(net.getVerticesCount(), 0);
	vector<typename Edge::IndexType> height(net.getVerticesCount(), 0);
	height[net.getSource()] = net.getVerticesCount();

	typename Net<Edge>::EdgeIterator iter = net.getConstOuterEdgeIterator(net.getSource());

	for (iter.begin(); !iter.end(); iter.move())
		push(net, iter, iter.getEdge().residualCapacity(), excess);

	globalDischargePreflowPushOperations(net, excess, height);
}


template <typename FirstType = unsigned, typename SecondType = long long>
struct SimpleEdge {
	typedef FirstType IndexType;
	typedef SecondType FlowType;
	IndexType from;
	IndexType to;
	FlowType flow;
	FlowType capacity;

	SimpleEdge (const IndexType &from, const IndexType &to, const FlowType &capacity) {
		this->from = from;
		this->to = to;
		this->flow = 0;
		this->capacity = capacity;
	}

	FlowType residualCapacity () const {
		return capacity - flow;
	}

	SimpleEdge reversed () const {
		return SimpleEdge(to, from, 0);
	}
};

template <class Edge>
void inputDirectedMaxFlow (istream &is, Net<Edge> &net) {
	typename Edge::IndexType verticesCount, edgesCount;
	is >> verticesCount >> edgesCount;

	net.setSource(0);
	net.setSink(verticesCount - 1);

	for (unsigned edgeIndex = 0; edgeIndex < edgesCount; edgeIndex ++) {
		typename Edge::IndexType from, to;
		typename Edge::FlowType capacity;
		cin >> from >> to >> capacity;
		-- from;
		-- to;

		net.addEdge(Edge(from, to, capacity));
	}
}

void outputDirectedMaxFlow (ostream &os, const long long &maxFlow) {
	os << maxFlow << endl;
}

void solveDirectedMaxFlow (istream &is, ostream &os) {
	Net<SimpleEdge<unsigned, long long> > net;
	
	inputDirectedMaxFlow(is, net);
		
	setMaxFlowWithPreflowPush(net);
	long long maxFlow = getFlowWithCheck(net);
	
	outputDirectedMaxFlow(os, maxFlow);
}


int main () {
	solveDirectedMaxFlow(cin, cout);
	return 0;
}
