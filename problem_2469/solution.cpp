#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdio.h>


using namespace std;


class SuffixAutomaton {
	private:
		struct Node {
			Node *suffixLink;
			Node *parent;
			vector<Node*> edges;
			size_t maxLength;
			size_t equivalenceClassSize;
			size_t num;
			bool used;

			Node () {
				equivalenceClassSize = 0;
				edges.resize(10, nullptr);
				parent = nullptr;
				used = false;
			}
		};

		Node *root;
		Node *last;
		Node *refrain;
		size_t size;


		void addSymbol (const char &ch) {
			Node *newLast = new Node;
			newLast->parent = last;
			newLast->maxLength = last->maxLength + 1;
			newLast->suffixLink = nullptr;
			newLast->num = (size ++);

			// Connect newLast to automaton:
			Node *current = last;
			
			while (current != nullptr && current->edges[ch - '0'] == nullptr) {
				current->edges[ch - '0'] = newLast;
				current = current->suffixLink;
			}

			// Calculate suffix link:
			if (current == nullptr)
				newLast->suffixLink = root;

			else if (current->maxLength + 1 == current->edges[ch - '0']->maxLength) 
				newLast->suffixLink = current->edges[ch - '0'];

			else {
				Node *nextCurrent = current->edges[ch - '0'];
				Node *clone = new Node;
				clone->parent = current;
				clone->num = (size ++);
				clone->suffixLink = nextCurrent->suffixLink;
				clone->maxLength = nextCurrent->maxLength;
				clone->edges = nextCurrent->edges;

				clone->maxLength = current->maxLength + 1;
				nextCurrent->suffixLink = newLast->suffixLink = clone;

				while (current != nullptr && current->edges[ch - '0'] == nextCurrent) {
					current->edges[ch - '0'] = clone;
					current = current->suffixLink;
				}
			}

			// Update last:
			last = newLast;
		}

		void setEquivalenceCountsOnTerminalVertices (const string &str) {
			Node *cur = root;

			for (const char &ch: str)
				cur = cur->edges[ch - '0'];

			while (cur != root) {
				cur->equivalenceClassSize = 1;
				cur = cur->suffixLink;
			}
		}

		void setEquivalenceCounts (Node *node) {
			node->used = true;

			for (char ch = '0'; ch <= '9'; ch ++)
				if (node->edges[ch - '0'] != nullptr && !node->edges[ch - '0']->used)
					setEquivalenceCounts(node->edges[ch - '0']);

			for (char ch = '0'; ch <= '9'; ch ++)
				if (node->edges[ch - '0'] != nullptr)
					node->equivalenceClassSize += node->edges[ch - '0']->equivalenceClassSize;
		}

		void getRefrain (Node *node) {
			node->used = false;

			if (refrain == nullptr || node->equivalenceClassSize*node->maxLength > refrain->equivalenceClassSize*refrain->maxLength)
				refrain = node;

			for (char ch = '0'; ch <= '9'; ch ++)
				if (node->edges[ch - '0'] != nullptr && node->edges[ch - '0']->used)
					getRefrain(node->edges[ch - '0']);
		}

		string getRefrainString () {
			string refrainString = "";

			vector<Node*> path;

			while (true) {
				path.push_back(refrain);

				if (refrain == root)
					break;

				refrain = refrain->parent;
			}

			reverse(path.begin(), path.end());

			for (size_t i = 0; i < path.size() - 1; i ++) {
				for (char ch = '0'; ch <= '9'; ch ++) {
					if (path[i]->edges[ch - '0'] == path[i + 1]) {
						refrainString.push_back(ch);
						break;
					}
				}
			}

			return refrainString;
		}

	public:
		explicit SuffixAutomaton (const string &str = "") {
			last = root = new Node();
			refrain = nullptr;
			root->suffixLink = nullptr;
			root->maxLength = 0;
			root->equivalenceClassSize = 0;
			root->num = ((size = 0) ++);

			for (const char &ch: str)
				addSymbol(ch);

			setEquivalenceCountsOnTerminalVertices(str);
			setEquivalenceCounts(root);
			getRefrain(root);

			const size_t refrainCharacteristic = refrain->equivalenceClassSize*refrain->maxLength;
			cout << refrainCharacteristic << endl;

			string refrainString = getRefrainString();
			cout << refrainString.size() << endl;

			for (const char &ch: refrainString)
				cout << ((size_t)(ch - '0') + 1) << " ";

			cout << endl;
		}

		~SuffixAutomaton () {
			// DELETE MEM:
		}

		void print (Node *node = nullptr) {
			if (node == nullptr)
				node = root;

			if (node->used)
				return;

			node->used = true;

			if (node != root)
				cout << node->num << "; suffixLink: " << node->suffixLink->equivalenceClassSize << "; len: " << node->maxLength << "; parent: " << node->parent->num << ": " << endl;

			else
				cout << "root: " << endl;

			for (char ch = '0'; ch <= '9'; ch ++)
				if (node->edges[ch - '0'] != nullptr)
					cout << "    " << ch << " -> " << node->edges[ch - '0']->num << endl;

			cout << endl;

			for (char ch = '0'; ch <= '9'; ch ++)
				if (node->edges[ch - '0'] != nullptr)
					print(node->edges[ch - '0']);
		}
};


int main () {
	//ifstream cin("input.txt");

	unsigned strSize;
	cin >> strSize;

	ushort c;
	cin >> c;

	string str;

	for (unsigned i = 0; i < strSize; i ++) {
		ushort ch;
		cin >> ch;
		ch --;
		str.push_back('0' + ch);
	}

	SuffixAutomaton suffixAutomaton(str);

	return 0;
}
