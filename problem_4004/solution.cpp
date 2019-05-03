#include <iostream>
#include <vector>
#include <queue>
#include <list>

using namespace std;


// Класс (Ориентированного) статического (нельзя менять кол-во вершин и удалять ребра) графа:
template <class Edge>
class Graph {
	private:
		vector<list<Edge> > edges;

	public:
		// Конструкторы:
		Graph () {}

		Graph (const typename Edge::IndexType &verticesCount) {
			edges.clear();
			edges.resize(verticesCount);
		}

		// Добавить ребро:
		void addEdge (const Edge &edge) {
			edges[edge.from].push_back(edge);
		}

		// Получить константные изходящие ребра:
		const list<Edge>& getConstOuterEdges (const typename Edge::IndexType &vertex) const {
			return edges[vertex];
		}

		// Получить изходящие ребра (с возможностью что-то в них отредактировать и даже сломать):
		list<Edge>& getOuterEdges (const typename Edge::IndexType &vertex) const {
			return edges[vertex];
		}

		// Получить кол-во вершин в графе:
		typename Edge::IndexType getVerticesCount () const {
			return edges.size();
		}

		void print (ostream &os) const {
			for (typename Edge::IndexType vertex = 0; vertex < getVerticesCount(); vertex ++) {
				os << vertex << ": ";

				for (const Edge &edge: edges[vertex]) {
					edge.print(os);
					os << " ";
				}

				os << endl;
			}
		}
};

// Алгоритм (dfs с помечением цвета ребер):
enum vertexColor {white, grey, black};

template <typename Edge>
bool findCycleDFS (const Graph<Edge> &graph, const typename Edge::IndexType &vertex, vector<vertexColor> &colors, vector<typename Edge::IndexType> &parents, typename Edge::IndexType &cycleEnd) {
	if (colors[vertex] == grey) {
		cycleEnd = vertex;
		return true;
	}

	else if (colors[vertex] == black)
		return false;

	colors[vertex] = grey;

	for (const Edge &edge: graph.getConstOuterEdges(vertex)) {
		parents[edge.to] = edge.from;
		//edge.print(cout << endl);

		if (findCycleDFS(graph, edge.to, colors, parents, cycleEnd)) {
			return true;
		}
	}

	colors[vertex] = black;
	return false;
}

template <typename Edge>
bool findCycle (const Graph<Edge> &graph, list<typename Edge::IndexType> &cycle) {
	vector<vertexColor> colors(graph.getVerticesCount(), white);
	vector<typename Edge::IndexType> parents(graph.getVerticesCount());
	typename Edge::IndexType cycleEnd;
	bool cycleExists = false;

	for (typename Edge::IndexType vertex = 0; vertex < graph.getVerticesCount(); vertex ++) {
		cycleExists |= findCycleDFS(graph, vertex, colors, parents, cycleEnd);
		
		if (cycleExists)
			break;
	}

	//cout << cycleExists << endl;

	if (!cycleExists)
		return false;

	typename Edge::IndexType vertex = cycleEnd;

	while (true) {
		cycle.push_front(vertex);

		if (parents[vertex] == cycleEnd)
			break;

		else
			vertex = parents[vertex];
	}

	return true;
}


// ----------------------------------------------------------------------------------------------------


// Кастомный класс Edge (может меняться от задачи, к задаче:
// например, в одних задач ребра взвешены, в других нет,
// поэтому класс граф как раз написан шаблонным)

struct SimpleEdge {
	typedef unsigned IndexType;
	IndexType from, to;

	SimpleEdge (const IndexType &from, const IndexType &to) {
		this->from = from;
		this->to = to;
	}

	void print (ostream &os) const {
		os << "(" << from << " -> " << to << ")";
	}
};

void input (istream &is, Graph<SimpleEdge> &graph) {
	unsigned verticesCount, edgesCount;
	is >> verticesCount >> edgesCount;
	graph = Graph<SimpleEdge>(verticesCount);

	for (unsigned i = 0; i < edgesCount; i ++) {
		unsigned from, to;
		is >> from >> to;
		graph.addEdge(SimpleEdge(from - 1 ,to - 1));
	}
}

void inputMatrix (istream &is, Graph<SimpleEdge> &graph) {
	unsigned verticesCount;
	is >> verticesCount;
	graph = Graph<SimpleEdge>(verticesCount);

	for (unsigned i = 0; i < verticesCount; i ++) {
		for (unsigned j = 0; j < verticesCount; j ++) {
			bool edgeExists;
			is >> edgeExists;

			if (edgeExists)
				graph.addEdge(SimpleEdge(i, j));
		}
	}
}

void output (ostream &os, const bool &cycleExists, const list<unsigned> &path) {
	if (!cycleExists)
		cout << "NO" << endl;

	else {
		os << "YES" << endl;

		for (const unsigned &vertex: path)
			cout << vertex + 1 << " ";

		cout << endl;
	}
}

void outputMatrix (ostream &os, const bool &cycleExists, const list<unsigned> &path) {
	os << cycleExists << endl;
}

int main () {
	Graph<SimpleEdge> graph;
	inputMatrix(cin, graph);
	//graph.print(cout << endl);
	
	list<unsigned> cycle;
	
	bool cycleExists = findCycle(graph, cycle);
	outputMatrix(cout, cycleExists, cycle);
	
	return 0;
}
