#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;
const static size_t maxLength = 1000000000;



class SuffixTree {
	private:
		struct Node {
			map<char, Node*> edge;
			Node *parent;
			size_t edgeLeft;
			size_t edgeLength;
			size_t num;
			Node *suffixLink = nullptr;
			size_t occuranceCount = 0;

			Node () {};
			Node (Node *other) : edge(other->edge), parent(other->parent), edgeLeft(other->edgeLeft), edgeLength(other->edgeLength), suffixLink(other->suffixLink) {}

			size_t getEdgesCount () {
				size_t edgesCount = 0;

				for (auto it = edge.begin(); it != edge.end(); it ++)
					if (it->second != nullptr)
						edgesCount ++;

				return edgesCount;
			}
		};

		struct Position { 
			Node *node;
			size_t distance;

			Position () : node(nullptr), distance(0) {}
			Position (Node *node, const size_t &distance) : node(node), distance(distance) {}

			// RULE 2:
			void splitEdge (const string &str, size_t &size) {
				Node *middle = new Node();
				middle->edge[str[node->edgeLeft + node->edgeLength - distance]] = node;
				middle->parent = node->parent;
				middle->parent->edge[str[node->edgeLeft]] = middle;
				middle->edgeLeft = node->edgeLeft;
				middle->edgeLength = node->edgeLength - distance;
				middle->num = (size ++);
				middle->suffixLink = nullptr;

				node->parent = middle;
				node->edgeLeft += middle->edgeLength;
				node->edgeLength = distance;

				this->node = middle;
				distance = 0;
			}

			void appendEdge (const char &ch, const string &str, size_t &size) {
				Node *newLeaf = new Node();
				newLeaf->num = (size ++);
				newLeaf->parent = node;
				newLeaf->edgeLength = maxLength;
				newLeaf->edgeLeft = str.size();
				newLeaf->suffixLink = nullptr;

				node->edge[ch] = newLeaf;
			}

			// RULE 3:
			bool canMoveBySymbol (const char &ch, const string &str) {
				if (distance == 0)
					return node->edge[ch] != nullptr;

				else
					return str[node->edgeLeft + node->edgeLength - distance] == ch;
			}

			// This method launches only when we can move by symbol ch:
			void moveBySymbol (const char &ch) {
				if (distance == 0) {
					node = node->edge[ch];
					distance = node->edgeLength - 1;
				}

				else
					distance --;
			}
		};

		size_t size;
		string str;
		Node *root;
		Position posLastNotLeaf;
		size_t lastNotLeaf;

	public:
		SuffixTree (const string &str) {
			root = new Node();
			root->num = ((size = 0) ++);
			root->parent = nullptr;
			root->edgeLeft = 0;
			root->edgeLength = 0;

			lastNotLeaf = 0;
			posLastNotLeaf = Position(root, 0);
 
			for (size_t i = 0; i < str.size(); i ++) {
				const char ch = str[i];
				addSymbol(ch);
			}

			setTerminals();
		}

		void addSymbol (const char &ch) {
			while (true) {
				if (lastNotLeaf > str.size())
					break;

				// RULE 3:
				if (posLastNotLeaf.canMoveBySymbol(ch, str)) {
					posLastNotLeaf.moveBySymbol(ch);
					break;
				}

				// RULE 2:
				else {
					if (posLastNotLeaf.distance != 0) {
						posLastNotLeaf.splitEdge(str, size);
						posLastNotLeaf.node->suffixLink = nullptr;
					}

					posLastNotLeaf.appendEdge(ch, str, size);

					if (posLastNotLeaf.node != root)
						posLastNotLeaf = Position(countSuffixLink(posLastNotLeaf.node), 0);
					
					lastNotLeaf ++;
				}
			}

			str += ch;
		}

	private:
		Node* countSuffixLink (Node *node) {
			if (node->suffixLink != nullptr)
				return node->suffixLink;

			else {
				// Определяем откуда и по какой подстроке будем спускаться:
				Node *halfSuffixLink;
				size_t substringLeft;
				size_t substringLength;

				if (node->parent != root) {
					halfSuffixLink = countSuffixLink(node->parent);
					substringLeft = node->edgeLeft;
					substringLength = node->edgeLength;
				}

				else {
					halfSuffixLink = root;
					substringLeft = node->edgeLeft + 1;
					substringLength = node->edgeLength - 1;
				}

				return node->suffixLink = goDown(halfSuffixLink, substringLeft, substringLength);
			}
		}

		Node* goDown (Node *halfSuffixLink, const size_t &substringLeft, const size_t &substringLength) {
			// Объявляем счетчик, который равен пройденному пути от halfSuffixLink:
			size_t currentLength = 0;

			// Спускаемся много раз:
			while (currentLength < substringLength) {
				halfSuffixLink = halfSuffixLink->edge[str[substringLeft + currentLength]];
				currentLength += halfSuffixLink->edgeLength;
			}

			// Если пришли в вершину, тогда все хорошо:
			if (currentLength == substringLength)
				return halfSuffixLink;

			// Иначе:
			else {
				// Сдвигаемся по ребру:
				Position position(halfSuffixLink, halfSuffixLink->edgeLength);
				position.distance = currentLength - substringLength;

				// Разделяем ребро в этом месте и возвращаем вершину на ребре:
				position.splitEdge(str, size);
				return position.node;
			}	
		}

		void setTerminals () {
			Position endPos = Position(root, 0);

			for (const char &ch: str)
				endPos.moveBySymbol(ch);

			endPos.splitEdge(str, size);

			Node *cur = countSuffixLink(endPos.node);

			while (cur != root) {
				cur->occuranceCount = 1;
				cur = countSuffixLink(cur);
			}
		}


	public:
		// PRINT ALL SUBSTRINGS:
		void printAllSubstringsDFS (ostream &os, Node *node, const bool flag = true, const size_t &depth = 0) {
			for (size_t i = 0; i < depth; i ++)
				os << "  ";

			if (node != root) {
				os << str.substr(node->edgeLeft, min<size_t>(node->edgeLength, str.size() - node->edgeLeft));

				if (flag)
					os << "(" << node->num << ", ";

				if (flag && node->parent != nullptr)
					os << node->parent->num;

				else if (flag)
					os << "X";

				os << ", ";

				if (flag && node->suffixLink != nullptr)
					os << node->suffixLink->num;

				else if (flag)
					os << "X";

				if (flag)
					os << ")";

				os << endl;
			}
			
			else
				os << "X" << endl;

			for (auto it = node->edge.begin(); it != node->edge.end(); it ++)
				if (it->second != nullptr)
					printAllSubstringsDFS(os, it->second, flag, depth + 1);
		}

		void printAllSubstrings (ostream &os) {
			printAllSubstringsDFS(os, root);
		}


		// COUNT ALL SUBSTRINGS:
		size_t countAllSubstringsDFS (Node *node) {
			size_t count = min<size_t>(node->edgeLength, str.size() - node->edgeLeft);

			for (auto it = node->edge.begin(); it != node->edge.end(); it ++)
				if (it->second != nullptr)
					count += countAllSubstringsDFS(it->second);

			return count;
		}

		size_t countAllSubstrings () {
			return countAllSubstringsDFS(root);
		}
};



int main () {
	string str;
	cin >> str;

	SuffixTree suffixTree(str);
	cout << suffixTree.countAllSubstrings() << endl;
}
