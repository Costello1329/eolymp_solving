#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
typedef long ld;


vector<ld> solveQuadraticEquasion (const ld &a, const ld &b, const ld &c) {
	vector<ld> roots;

	/*if (a == 0) {
		if (b == 0)
			roots = vector<ld>(3);

		else
			roots.push_back(- c/b);
	}

	else {*/
		ld D2 = b*b - 4*a*c;

		if (D2 > 0) {
			roots.push_back((-b + sqrt(D2)) / (2*a));
			roots.push_back((-b - sqrt(D2)) / (2*a));
		}

		else if (D2 == 0)
			roots.push_back(- b/(2*a));
	//}

	return roots;
}


int main () {
	ld a, b, c;
	cin >> a >> b >> c;

	vector<ld> roots = solveQuadraticEquasion(a, b, c);
	sort(roots.begin(), roots.end());

	if (roots.size() > 3)
		cout << "Infinite roots" << endl;

	else if (roots.size() == 2)
		cout << "Two roots: " << roots[0] << " " << roots[1] << endl;

	else if (roots.size() == 1)
		cout << "One root: " << roots[0] << endl;

	else
		cout << "No roots" << endl;

	return 0;
}
