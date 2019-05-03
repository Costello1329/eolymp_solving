#include <iostream>
#include <stack>
#include <vector>
#include <list>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <set>
#include <unordered_map>
 
using namespace std;

 
 
// Vector 3D: ----------------------------------------------------------------------------------------------------
 
template <typename Type>
struct Vector3D {
	Type x, y, z;
 
	Vector3D (const Type &x = 0, const Type &y = 0, const Type &z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	
	Vector3D (const Vector3D<Type> &other) {
		this->x = other.x;
		this->y = other.y;
		this->z = other.z;
	}
};
 
template <typename Type>
Vector3D<Type> operator + (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
	return Vector3D<Type>(firstVector.x + secondVector.x, firstVector.y + secondVector.y, firstVector.z + secondVector.z);
}
 
template <typename Type>
Vector3D<Type> operator * (const Type &number, const Vector3D<Type> &vec) {
	return Vector3D<Type>(number * vec.x, number * vec.y, number * vec.z);
}
 
template <typename Type>
Vector3D<Type> operator * (const Vector3D<Type> &vec, const Type &number) {
	return number * vec;
}
 
template <typename Type>
Vector3D<Type> operator / (const Type &number, const Vector3D<Type> &vec) {
	return vec / number;
}
 
template <typename Type>
Vector3D<Type> operator / (const Vector3D<Type> &vec, const Type &number) {
	return 1/number * vec;
}
 
template <typename Type>
Vector3D<Type> operator - (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
	return Vector3D<Type>(firstVector.x + -1 * secondVector.x, firstVector.y -1 * secondVector.y, firstVector.z + -1 * secondVector.z);
}
 
template <typename Type>
Type dotProduct (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
	return firstVector.x*secondVector.x + firstVector.y*secondVector.y + firstVector.z*secondVector.z;
}
 
template <typename Type>
Vector3D<Type> crossProduct (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
	return Vector3D<Type>(firstVector.y*secondVector.z - firstVector.z*secondVector.y, - firstVector.x*secondVector.z + firstVector.z*secondVector.x, firstVector.x*secondVector.y - firstVector.y*secondVector.x);
}
 
template <typename Type>
Type getLength (const Vector3D<Type> &vec) {
	return sqrt(dotProduct(vec, vec));
}
 
 
//Point3D<Type>: ----------------------------------------------------------------------------------------------------
 
template <typename Type>
struct Point3D {
	typedef Type CoordinateType;
	Type x, y, z;
 
	Point3D (const Type &x = 0, const Type &y = 0, const Type &z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
};
 
template <typename Type>
bool operator == (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return (firstPoint.x == secondPoint.x && firstPoint.y == secondPoint.y && firstPoint.z == secondPoint.z);
}
 
template <typename Type>
bool operator != (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return !(firstPoint == secondPoint);
}
 
template <typename Type>
bool operator < (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return firstPoint.x < secondPoint.x || (firstPoint.x == secondPoint.x && (firstPoint.y < secondPoint.y || (firstPoint.z == secondPoint.z && firstPoint.z < secondPoint.z)));
}
 
template <typename Type>
bool operator <= (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return firstPoint < secondPoint || firstPoint == secondPoint;
}
 
template <typename Type>
bool operator > (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return !(firstPoint < secondPoint) && firstPoint != secondPoint;
}
 
template <typename Type>
bool operator >= (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return firstPoint > secondPoint || firstPoint == secondPoint;
}
 
template <typename Type>
Vector3D<Type> operator - (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return Vector3D<Type>(firstPoint.x - secondPoint.x, firstPoint.y - secondPoint.y, firstPoint.z - secondPoint.z);
}
 
template <typename Type>
Point3D<Type> operator + (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
	return Point3D<Type>(firstPoint.x + secondPoint.x, firstPoint.y + secondPoint.y, firstPoint.z + secondPoint.z);
}
 
template <typename Type>
Point3D<Type> operator + (const Point3D<Type> &point, const Vector3D<Type> &vec) {
	return Point3D<Type>(point.x + vec.x, point.y + vec.y, point.z + vec.z);
}
 
template <typename Type>
ostream& operator << (ostream &os, const Point3D<Type> &point) {
	os << point.x << " " << point.y << " " << point.z;
	return os;
}
 
 
// Convex Hull 3D: ----------------------------------------------------------------------------------------------------
 
template<typename Type>
class ConvexHull3D {
	private:
		struct Point : public Point3D<Type> {
			unsigned number;

			Point (const Point3D<Type> &point = Point3D<Type>(), const unsigned &number = 0) {
				this->x = point.x;
				this->y = point.y;
				this->z = point.z;
				this->number = number;
			}
		};

		class Edge {
			private:
				pair<Point, Point> edge;

			public:
				Edge (const Point &firstPoint, const Point &secondPoint) {
					edge = make_pair(firstPoint, secondPoint);
				}

				bool operator == (const Edge &otherEdge) const {
					return edge.first == otherEdge.edge.first && edge.second == otherEdge.edge.second;
				}

				unsigned getEdgeNumber () const {
					unsigned first = edge.first.number;
					unsigned second = edge.second.number;
					return min<unsigned>(first, second) * 1000 ^ max<unsigned>(first, second);
				}

				const Point &from () const {
					return edge.first;
				}

				const Point &to () const {
					return edge.second;
				}
		};

		class Face {
			private:
				Point firstPoint;
				Point secondPoint;
				Point thirdPoint;
	   
			public:
				Face (const Point &firstPoint, const Point &secondPoint, const Point &thirdPoint) : firstPoint(firstPoint), secondPoint(secondPoint), thirdPoint(thirdPoint) {}

				bool isPointOnFace (const Point &point) const {
					return firstPoint == point || secondPoint == point || thirdPoint == point;
				}

				Vector3D<Type> getNormal () const {
					Vector3D<Type> firstVector = (Point3D<Type>)secondPoint - (Point3D<Type>)firstPoint;
					Vector3D<Type> secondVector = (Point3D<Type>)thirdPoint - (Point3D<Type>)firstPoint;
					Vector3D<Type> normal = crossProduct(firstVector, secondVector);
					normal = normal / getLength(normal);
					return normal;
				}

				void reorient () {
					swap(secondPoint, thirdPoint);
				}

				vector<Point> getPoints () const {
					vector<Point> points;
					points.push_back(firstPoint);
					points.push_back(secondPoint);
					points.push_back(thirdPoint);
					return points;
				}

				void print () const {
					const string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
					cout << alphabet[firstPoint.number] << alphabet[secondPoint.number] << alphabet[thirdPoint.number] << endl;
				}
		};
 
		vector<Face> faces;
 
	public:
		ConvexHull3D (const vector<Point3D<Type> > &points) {
			primarySetup(points);
			secondarySetup(points);
		}
 
 		// Первый полигон:
		void primarySetup (const vector<Point3D<Type> > &points) {
			// Ищем первую точку:
			size_t firstPointIndex = 0;

			for (size_t i = 0; i < points.size(); i ++)
				if (points[i].z < points[firstPointIndex].z)
					firstPointIndex = i;

			// Ищем вторую точку:
			size_t secondPointIndex = numeric_limits<unsigned>::max();

			for (size_t i = 0; i < points.size(); i ++) {
				if (i == firstPointIndex)
					continue;

				if (secondPointIndex == numeric_limits<unsigned>::max()) {
					secondPointIndex = i;
					continue;
				}

				Vector3D<Type> firstVector = points[i] - points[firstPointIndex];
				firstVector = firstVector / getLength(firstVector);
				Vector3D<Type> secondVector = points[secondPointIndex] - points[firstPointIndex];
				secondVector = secondVector / getLength(secondVector);

				if (dotProduct(firstVector, Vector3D<Type>(0, 0, 1)) < dotProduct(secondVector, Vector3D<Type>(0, 0, 1)))
					secondPointIndex = i;
			}

			// Создаем грань таким образом, чтобы все точки были с одной стороны от нее:
			const Edge edge(Point(points[firstPointIndex], firstPointIndex), Point(points[secondPointIndex], secondPointIndex));
			unsigned thirdPointIndex = numeric_limits<unsigned>::max();

			for (unsigned i = 0; i < points.size(); i ++) {
				if (firstPointIndex == i || secondPointIndex == i)
					continue;

				if (thirdPointIndex == numeric_limits<unsigned>::max())
					thirdPointIndex = i;

				else {
					Face oldFace(edge.from(), edge.to(), Point(points[thirdPointIndex], thirdPointIndex));
					Face newFace(edge.from(), edge.to(), Point(points[i], i));

					// Если новый полигон лучше старого, обновляем thirdPointIndex:
					if (dotProduct(crossProduct(oldFace.getNormal(), newFace.getNormal()), (Point3D<Type>)edge.to() - (Point3D<Type>)edge.from()) > 0)
						thirdPointIndex = i;
				}
			}

			// Создаем грань:
			Face face(Point(points[firstPointIndex], firstPointIndex), Point(points[secondPointIndex], secondPointIndex), Point(points[thirdPointIndex], thirdPointIndex));

			for (size_t i = 0; i < points.size(); i ++)
				if (i != firstPointIndex && i != secondPointIndex && i != thirdPointIndex)
					if (dotProduct(face.getNormal(), points[firstPointIndex] - points[i]) < 0)
						face.reorient();

			faces.push_back(face);
		}


		// Заворачивание:
		void addAllFaceEdgesToStack (const Face &face, stack<pair<Face, Edge> > &edgesStack, unordered_map<unsigned, ushort> &used) const {
			Edge firstEdge(face.getPoints()[0], face.getPoints()[1]);
			Edge secondEdge(face.getPoints()[1], face.getPoints()[2]);
			Edge thirdEdge(face.getPoints()[2], face.getPoints()[0]);

			if (used[firstEdge.getEdgeNumber()] != 2)
				edgesStack.push(make_pair(face, firstEdge));
			
			if (used[secondEdge.getEdgeNumber()] != 2)
				edgesStack.push(make_pair(face, secondEdge));
			
			if (used[thirdEdge.getEdgeNumber()] != 2)
				edgesStack.push(make_pair(face, thirdEdge));

			used[firstEdge.getEdgeNumber()] ++;
			used[secondEdge.getEdgeNumber()] ++;
			used[thirdEdge.getEdgeNumber()] ++;
		}

		void secondarySetup (const vector<Point3D<Type> > &points) {
			unordered_map<unsigned, ushort> used;

			stack<pair<Face, Edge> > edgesStack;
			addAllFaceEdgesToStack(faces.back(), edgesStack, used);
 
			while (!edgesStack.empty()) {
				const Face face = edgesStack.top().first;
				const Edge edge = edgesStack.top().second;
				const string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
				edgesStack.pop();

				if (used[edge.getEdgeNumber()] == 2)
					continue;

				// Если ребро еще не было 2 раза использовано в convex hull'e, то надо построить на нем новый полигон:
				unsigned bestPointIndex = numeric_limits<unsigned>::max();

				for (unsigned i = 0; i != points.size(); i ++) {
					if (face.isPointOnFace(Point(points[i], i)))
						continue;

					if (bestPointIndex == numeric_limits<unsigned>::max())
						bestPointIndex = i;

					else {
						Face oldFace(edge.to(), edge.from(), Point(points[bestPointIndex], bestPointIndex));
						Face newFace(edge.to(), edge.from(), Point(points[i], i));

						// Если новый полигон лучше старого, обновляем bestPointIndex:
						if (dotProduct(newFace.getNormal(), face.getNormal()) > dotProduct(oldFace.getNormal(), face.getNormal()))
							bestPointIndex = i;
					}
				}

				Face bestFace(edge.to(), edge.from(), Point(points[bestPointIndex], bestPointIndex));
				faces.push_back(bestFace);
				addAllFaceEdgesToStack(bestFace, edgesStack, used);
			}
		}

		// Получить вывод грани (относительно нужного циклического сдвига (в задаче надо вывести первую по номеру точку на гране первой))
		tuple<unsigned, unsigned, unsigned> rotateFace (const Face &face) {
			vector<unsigned> numbers;
			vector<Point> points = face.getPoints();

			for (const Point &point: points)
				numbers.push_back(point.number);

			unsigned minIndex = 0;

			for (unsigned i = 1; i < 3; i ++)
				if (numbers[i] < numbers[minIndex])
					minIndex = i;

			return make_tuple(numbers[minIndex % 3], numbers[(minIndex + 1) % 3], numbers[(minIndex + 2) % 3]);
		}

		// Просто печатаем в лексикографическом порядке Face'ы. (для этого нам нужен tuple, который сортится сам хорошо с std::sort, так как в STL написан компаратор less для tuple)
		void print (ostream &os) {
			set<tuple<unsigned, unsigned, unsigned> > rotatedFaces;
			
			for (const Face &face: faces)
				rotatedFaces.insert(rotateFace(face));
			
			os << rotatedFaces.size() << endl;
			
			for (auto i = rotatedFaces.begin(); i != rotatedFaces.end(); i ++)
				os << "3 " << get<0>(*i) << " " << get<1>(*i) << " " << get<2>(*i) << endl;
		}
};


int main () {
	unsigned T;
	cin >> T;

	for (unsigned i = 0; i < T; i ++) {
		unsigned pointsCount;
		cin >> pointsCount;

		vector<Point3D<double> > points(pointsCount);

		for (Point3D<double> &point: points)
			cin >> point.x >> point.y >> point.z;

		ConvexHull3D<double> convexHull(points);
		convexHull.print(cout);
	}

	return 0;
}
