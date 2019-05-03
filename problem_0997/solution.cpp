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

// Алгоритм (один bfs):
template <typename Edge>
bool getShortestPathBetweenTwoVertices (const Graph<Edge> &graph, const typename Edge::IndexType &from, const typename Edge::IndexType &to, list<typename Edge::IndexType> &path) {
	queue<typename Edge::IndexType> verticesQueue;
	vector<bool> used(graph.getVerticesCount(), false);
	vector<typename Edge::IndexType> parents(graph.getVerticesCount());
	
	verticesQueue.push(from);
	used[from] = true;
	parents[from] = from;

	
	while (!verticesQueue.empty()) {
		const typename Edge::IndexType vertex = verticesQueue.front();
		verticesQueue.pop();

		if (vertex == to)
			break;

		for (const Edge &edge: graph.getConstOuterEdges(vertex)) {
			if (!used[edge.to]) {
				used[edge.to] = true;
				verticesQueue.push(edge.to);
				parents[edge.to] = edge.from;
			}
		}
	}

	
	if (!used[to])
		return false;

	else {
		typename Edge::IndexType vertex = to;
		
		while (true) {
			path.push_front(vertex);

			if (vertex == parents[vertex])
				break;
			
			else
				vertex = parents[vertex];
		}
	
		return true;
	}
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

struct ChessfieldCell {
	unsigned x, y;
	unsigned chessfieldSize;

	ChessfieldCell () {}

	ChessfieldCell (const unsigned &index, const unsigned &chessfieldSize) {
		x = index % chessfieldSize;
		y = index / chessfieldSize;
		this->chessfieldSize = chessfieldSize;
	}

	ChessfieldCell (const unsigned &x, const unsigned &y, const unsigned &chessfieldSize) {
		this->x = x;
		this->y = y;
		this->chessfieldSize = chessfieldSize;
	}

	unsigned getIndex () const {
		return x + y * chessfieldSize;
	}
};

inline list<ChessfieldCell> getAdmissibleCells (const ChessfieldCell &cell) {
	list<ChessfieldCell> cells;
	const vector<short> deltaX = {1, 2, 2, 1, -1, -2, -2, -1};
	const vector<short> deltaY = {2, 1, -1, -2, -2, -1, 1, 2};

	for (ushort i = 0; i < 8; i ++)
		if ((int)cell.x + deltaX[i] >= 0 && (int)cell.y + deltaY[i] >= 0 && (int)cell.x + deltaX[i] < cell.chessfieldSize && (int)cell.y + deltaY[i] < cell.chessfieldSize)
			cells.push_back(ChessfieldCell(cell.x + deltaX[i], cell.y + deltaY[i], cell.chessfieldSize));

	return cells;
}

void buildChessfieldGraphDFS (Graph<SimpleEdge> &graph, const ChessfieldCell &cell, vector<bool> &used) {
	if (used[cell.getIndex()])
		return;

	used[cell.getIndex()] = true;

	list<ChessfieldCell> admissibleCells = getAdmissibleCells(cell);

	for (const ChessfieldCell &nextCell: admissibleCells) {
		graph.addEdge(SimpleEdge(cell.getIndex(), nextCell.getIndex()));
		buildChessfieldGraphDFS(graph, nextCell, used);
	}
}

void buildChessfieldGraph (Graph<SimpleEdge> &graph, const unsigned &chessfieldSize) {
	graph = Graph<SimpleEdge>(chessfieldSize * chessfieldSize);
	vector<bool> used(chessfieldSize * chessfieldSize, false);
	ChessfieldCell start = ChessfieldCell(0, 0, chessfieldSize);
	buildChessfieldGraphDFS(graph, start, used);
}

void output (ostream &os, const unsigned &chessfieldSize, const list<unsigned> &path) {
	os << path.size() - 1 << endl;

	for (const unsigned &index: path) {
		ChessfieldCell cell(index, chessfieldSize);
		//os << cell.x + 1 << " " << cell.y + 1 << endl;
	}
}

int main () {
	unsigned chessfieldSize;
	cin >> chessfieldSize;
	Graph<SimpleEdge> graph;
	buildChessfieldGraph(graph, chessfieldSize);
	//graph.print(cout);

	ChessfieldCell from;
	ChessfieldCell to;
	from.chessfieldSize = chessfieldSize;
	to.chessfieldSize = chessfieldSize;
	cin >> from.x >> from.y >> to.x >> to.y;
	from.x --;
	from.y --;
	to.x --;
	to.y --;

	list<unsigned> path;
	getShortestPathBetweenTwoVertices(graph, from.getIndex(), to.getIndex(), path);
	output(cout, chessfieldSize, path);
	
	return 0;
}
