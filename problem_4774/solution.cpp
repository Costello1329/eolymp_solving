#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const static long double epsilon = 0.00000000001f;
const static long double pi = 3.141592653589793238462643f;


struct Point2D {
	long long x;
	long long y;
	bool isNull;

	Point2D () {
		this->x = 0;
		this->y = 0;
	}

	Point2D (long long x, long long y) {
		this->x = x;
		this->y = y;
	}
};

struct Vector2D {
	long long x, y;

	Vector2D () {
		x = 0;
		y = 0;
	}

	Vector2D (const long long &x, const long long &y) {
		this->x = x;
		this->y = y;
	}

	double length () const {
		return sqrt(x*x + y*y);
	}
};

Vector2D operator - (const Point2D &point1, const Point2D &point2) {
	return Vector2D(point1.x - point2.x, point1.y - point2.y);
}

bool operator == (const Point2D &point1, const Point2D &point2) {
	return point1.x == point2.x && point1.y == point2.y;
}

bool operator != (const Point2D &point1, const Point2D &point2) {
	return !(point1 == point2);
}

ostream& operator << (ostream &os, const Point2D &point) {
	return os << point.x << " " << point.y;
}

ostream& operator << (ostream &os, const Vector2D &vector) {
	return os << vector.x << " " << vector.y;
}


double countPolarAngle (const Point2D &start, const Point2D &end) {
	long long x = end.x - start.x;
	long long y = end.y - start.y;
	long double angle = atan2(y, x);

	if (angle < 0)
		angle += 2 * pi;

	return angle;
}

bool isLeftTurn (const Point2D &point1, const Point2D &point2, const Point2D &point3) {
	return (countPolarAngle(point2, point3) - countPolarAngle(point1, point2) >= - epsilon);
}

double dot (const Vector2D &vector1, const Vector2D &vector2) {
	return vector1.x * vector2.x + vector1.y * vector2.y;
}

long double cross (const Vector2D &vector1, const Vector2D &vector2) {
	return vector1.x * vector2.y - vector1.y * vector2.x;
}

double countAngleBetweenVectors (const Vector2D &vector1, const Vector2D &vector2) {
	return acos(dot(vector1, vector2) / vector1.length() / vector2.length());
}

double distance (const Point2D &point1, const Point2D &point2) {
	return sqrt(pow(point2.x - point1.x, 2) + pow(point2.y - point1.y, 2));
}


void removeEqual (vector<Point2D> &points) {
	sort(points.begin(), points.end(), [] (const Point2D &point1, const Point2D &point2) {
		if (point1 == point2)
			return true;

		else
			return (point1.x <= point2.x && point1.y <= point2.y);
	});



	bool used = false;
	Point2D prev;
	vector<Point2D> goodPoints;

	for (const Point2D &point: points) {
		if (!used || point != prev) {
			used = true;
			goodPoints.push_back(prev = point);
		}
	}

	points = goodPoints;
}


vector<Point2D> buildConvexHull (const vector<Point2D> &points) {
	if (points.size() <= 3)
		return points;


	Point2D firstPoint = points[0];

	for (Point2D point: points)
		if (point.y < firstPoint.y || (point.y == firstPoint.y && point.x < firstPoint.x))
			firstPoint = point;


	vector<Point2D> sortedPoints = points;

	sort(sortedPoints.begin(), sortedPoints.end(), [firstPoint] (const Point2D &point1, const Point2D &point2) {
		long double angle1 = countPolarAngle(firstPoint, point1);
		long double angle2 = countPolarAngle(firstPoint, point2);

		if (abs(angle2 - angle1) < epsilon) {
			long double distance1 = distance(firstPoint, point1);
			long double distance2 = distance(firstPoint, point2);
			return distance1 < distance2;
		}

		else
			return angle1 < angle2;
	});

	vector<Point2D> convexHull;
	convexHull.push_back(sortedPoints[0]);
	convexHull.push_back(sortedPoints[1]);

	for (size_t i = 2; i < sortedPoints.size(); i ++) {
		while (convexHull.size() > 1 && !isLeftTurn(convexHull[convexHull.size() - 2], convexHull[convexHull.size() - 1], sortedPoints[i]))
			convexHull.pop_back();

		convexHull.push_back(sortedPoints[i]);
	}

	return convexHull;
}


void printArea (const vector<Point2D> &convexHull) {
	unsigned long long area = 0;

	for (size_t i = 0; i < convexHull.size() - 2; i ++) {
		Vector2D vec1 = convexHull[i + 1] - convexHull[0];
		Vector2D vec2 = convexHull[i + 2] - convexHull[0];
		area += abs(cross(vec1, vec2));
	}

	cout << area / 2;
	if (area % 2 == 1)
		cout << ".5";
	cout << endl;
}


int main () {
	size_t pointsNumber;
	cin >> pointsNumber;

	vector<Point2D> points(pointsNumber);

	for (size_t i = 0; i < pointsNumber; i ++)
		cin >> points[i].x >> points[i].y;


	removeEqual(points);

	vector<Point2D> convexHull = buildConvexHull(points);



	vector<Point2D> goodPoints;

	for (size_t i = 0; i < convexHull.size(); i ++) {
		Point2D prev = (i != 0 ? convexHull[i - 1] : convexHull[convexHull.size() - 1]);
		Point2D cur = convexHull[i];
		Point2D next = convexHull[(i + 1) % convexHull.size()];
		Vector2D vec1 = cur - prev;
		Vector2D vec2 = next - cur;

		//cout << vec1 << " " << vec2 << endl;

		if (cross(vec1, vec2) > epsilon)
			goodPoints.push_back(cur);
	}

	cout << goodPoints.size() << endl;


	for (const Point2D &point: goodPoints)
		cout << point << endl;

	//printArea(goodPoints);

	return 0;
}
