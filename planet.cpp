#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <utility>
#include <cassert>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace std::chrono;

#define point pair<double, double>
#define inter pair<bool, point>
const double doublePrecision = 1e-6;
const double delta = 1.0;
const double upperDelta = 1.0;

static void print(point u) {
	cerr << setprecision(15) << u.first << " " << u.second << endl;
}

// Read pair of doubles.
static vector<point> readBinary(char* fileName, int limitSize) {
	ifstream input(fileName, ios::binary);
    auto pos = input.tellg();
    input.seekg( 0, ios::end );
    auto size=input.tellg()-pos;
    input.seekg(0,ios::beg);
/*
	size >>= 1;
	cerr <<"size = " << size << endl;
	int sizeToRead = limitSize == -1 ? size : (limitSize < size ? limitSize : size);
	vector<point> coord(sizeToRead);
	cerr << "sizeToRead = " << sizeToRead << endl; 
	for (unsigned index = 0; index < sizeToRead; ++index) {
		double x, y;
		input.read(reinterpret_cast<char*>(&x), sizeof(double));
		input.read(reinterpret_cast<char*>(&y), sizeof(double));
		long double x0 = static_cast<long double>(x);
		long double y0 = static_cast<long double>(y);
		coord[index] = make_pair(x0, y0 + delta);
		if (index < 10)
			print(coord[index]);
	}
	return coord;
*/
    unsigned elements = size / sizeof(point);
    if (limitSize != -1)
    	elements = limitSize < elements ? limitSize : elements;
    vector<point> cdf;
    cdf.resize(elements);
    input.read(reinterpret_cast<char*>(cdf.data()), elements * sizeof(point));
    return cdf;
}

// Write the sample points in a binary file.
static void writeBinary(char* fileName, vector<point>& samplePoints) {
	ofstream output(fileName, ios::binary);
	for (point p: samplePoints) {
		double x = static_cast<double>(p.first);
		double y = static_cast<double>(p.second - delta);
		output.write(reinterpret_cast<char*>(&x), sizeof(double));
		output.write(reinterpret_cast<char*>(&y), sizeof(double));
	}
}

// Read from a .txt file.
static vector<point> readTxt(char* fileName) {
	ifstream input(fileName);
	vector<point> coord;
	double x, y;
	while (input >> x >> y) {
		long double x0 = static_cast<long double>(x);
		long double y0 = static_cast<long double>(y);
		coord.push_back(make_pair(x0, y0 + delta));
	}
	return coord;
}

// Write into a .txt file
static void writeTxt(char* fileName, vector<point>& samplePoints) {
	ofstream output(fileName);
	for (point p: samplePoints)
		output << p.first << " " << (p.second - delta) << "\n";
}

// Computes the position of the triple (u, v, w) in respect to pi.
static int angle(bool counterclockwise, point u, point v, point w) {
	double ret = ((v.first - u.first) * (w.second - u.second) - (v.second - u.second) * (w.first - u.first));
	return counterclockwise ? ret : -ret;
}

static double myabs(long double x) {
	return (x < -doublePrecision) ? -x : x;
}

// Intersection of the lines (u1, u2) and (v1, v2).
static inter intersection(point u1, point u2, point v1, point v2) {
	double a1 = u2.second - u1.second, b1 = u1.first - u2.first, c1 = a1 * u1.first + b1 * u1.second;
	double a2 = v2.second - v1.second, b2 = v1.first - v2.first, c2 = a2 * v1.first + b2 * v1.second;
	double det = a1 * b2 - a2 * b1;

	// compute with long double precision.
	if (abs(det) < doublePrecision) {
		cerr << "ba e determinantul zero!!! Ia sa vedem, totusi" << endl;
		print(u1);
		print(u2);
		print(v1);
		print(v2);
		return make_pair(false, u1);
	}
	return make_pair(true, make_pair((c1 * b2 - c2 * b1) / det, (a1 * c2 - a2 * c1) / det));
	/*if ((min(v1.first, v2.first) <= ret.first) && (ret.first <= max(v1.first, v2.first)) && (min(v1.second, v2.second) <= ret.second) && (ret.second <= max(v1.second, v2.second))) {
		return ret;
	} else {
		print(ret);
		print(u1);
		print(u2);
		print(v1);
		print(v2);
		//cerr << u1.firs << " "  << u2 << " " << v1 << " " << v2 << endl;
		assert(0);
	}*/
}

double globalEpsForChange;
double minusEpsForChange;

point change(bool sign, vector<point>& coord, unsigned index) {
    return sign ? make_pair(coord[index - 1].first, coord[index - 1].second * (1 + globalEpsForChange) + delta + upperDelta) 
                : make_pair(coord[index - 1].first, coord[index - 1].second + delta);
}

// Computes the samples points.
vector<point> samplePointsOfApproximateFunction(vector<point>& coord, double eps) { 
	unsigned n = coord.size();

    globalEpsForChange = eps;
    minusEpsForChange = 1e-2;
    
	cerr << "size of coord in funciton = " << n << endl;

	point *windowBound = new point[2], *leftBound = new point[2], *rightBound = new point[2];

	// lists of neighbours.
	unordered_map<point, point, boost::hash<point>> *rightList = new unordered_map<point, point, boost::hash<point>>[2];
	unordered_map<point, point, boost::hash<point>> *leftList = new unordered_map<point, point, boost::hash<point>>[2];
	
	//map<point, point> *rightList = new map<point, point>[2];
	//map<point, point> *leftList = new map<point, point>[2]; 

	cerr << "before alloc" << endl;

	vector<point> samplePoints;
	
	cerr << "wwas???" << endl;

	// Maybe we should reverse the order in points.
	/*for (unsigned index = 1; index <= n; ++index) {
		points[index][0] = make_pair(coord[index - 1].first, coord[index - 1].seconc * (1 - eps) + delta);
		points[index][1] = make_pair(coord[index - 1].first, coord[index - 1].second * (1 + eps) + delta + upperDelta);
	}*/

	cerr << "Umwandeln" << endl;

	for (int sign = 1; sign >= 0; --sign) {
        point p1 = change(sign, coord, 1), p2 = change(sign, coord, 2);
        
		windowBound[sign] = p1;//points[1][sign];
		leftBound[sign] = p1;//points[1][sign];
		rightBound[sign] = p1;//points[1][sign];

		rightList[sign][p1] = p2;//points[2][sign];
		leftList[sign][p2] = p1;//points[1][sign];
	}

	for (unsigned index = 3; index < n; ++index) {
		//cerr << "now at index = " << index << endl;
		bool nextWindow = false;
		// Updating convex hulls
		for (int sign = 1; sign >= 0; --sign) {
            point indexPoint = change(sign, coord, index);
			point currPoint = change(sign, coord, index - 1);
            
			while ((currPoint != windowBound[sign]) && (angle(sign, indexPoint, currPoint, leftList[sign][currPoint]) > doublePrecision))
				currPoint = leftList[sign][currPoint];
			rightList[sign][currPoint] = indexPoint; //[index][sign];
			leftList[sign][indexPoint] = currPoint;
		}

		for (int sign = 1; sign >= 0; --sign) {
			int star = sign, diamond = !sign;
            point indexPointStar = change(star, coord, index);
            point indexPointDiamond = change(diamond, coord, index);
            point lastIndexPoint = change(star, coord, index - 1);
            
			if ((!nextWindow) && (angle(star, indexPointStar, leftBound[star], rightBound[diamond]) < -doublePrecision)) {
				// Get new sample point.
				//cerr << "found sth!" << endl;
				inter value = intersection(leftBound[star], rightBound[diamond], windowBound[star], windowBound[diamond]);
				if (value.first == true) {
					samplePoints.push_back(value.second);
				} else {
					cerr << "For one sample point the determinant is 0!" << endl;
				}
				//cerr << "it's not him" << endl;

				// Update the window.
				windowBound[diamond] = rightBound[diamond];
				value = intersection(leftBound[star], rightBound[diamond], lastIndexPoint, indexPointStar);
				if (value.first == true)
					windowBound[star] = value.second;
				else 
					assert(0);

				// Update the lists.
				rightList[star][windowBound[star]] = indexPointStar;
				leftList[star][indexPointStar] = windowBound[star];

				// Update rightBound
				rightBound[star] = indexPointStar;
				rightBound[diamond] = indexPointDiamond;

				// Update leftBound
				leftBound[star] = windowBound[star];
				leftBound[diamond] = windowBound[diamond];
			
				while (angle(diamond, leftBound[diamond], rightBound[star], rightList[diamond][leftBound[diamond]]) < -doublePrecision)
					leftBound[diamond] = rightList[diamond][leftBound[diamond]];
				nextWindow = true;
			}
		}

		// Updating the supports and separating lines.
		if (!nextWindow) {
			for (int sign = 1; sign >= 0; --sign) {
				int star = sign, diamond = !sign;
                point indexPoint = change(star, coord, index);
                
				if (angle(star, indexPoint, leftBound[diamond], rightBound[star]) < -doublePrecision) {
					rightBound[star] = indexPoint;
					while (angle(star, indexPoint, leftBound[diamond], rightList[diamond][leftBound[diamond]]) < -doublePrecision) 
						leftBound[diamond] = rightList[diamond][leftBound[diamond]];
				}
			}
		}
	}

	// Compute the two last approximated points.
	// Quite weird here!
	cerr << "by now" << samplePoints.size() << endl;
	// I'm still confused if there is need for a condition here. 
	cerr << "n-ter Eintrag = " << endl;
	//print(points[n][1]);
	//print(points[n][0]);
	cerr << "leftBound + rightBound + windowBound" << endl;
	for (int sign = 1; sign >= 0; --sign) {
		print(leftBound[sign]);
		print(rightBound[sign]);
		print(windowBound[sign]);
	}

	inter value = intersection(leftBound[0], rightBound[1], windowBound[1], windowBound[0]);
	if (value.first == true)
		samplePoints.push_back(value.second);
	else
		cerr << "Penultimul punct nu exista!" << endl;

	value = intersection(leftBound[0], rightBound[1], change(1, coord, n), change(0, coord, n));
	if (value.first == true)
		samplePoints.push_back(value.second);
	else
		cerr << "Ultimul punct nu exista!" << endl;
	return samplePoints;
}

static double interpolate(vector<point>& spline, double pos, bool dump) {
   if (pos <= spline.front().first)
      return spline.front().second;
   if (pos >= spline.back().first)
      return spline.back().second;
  
	auto iter=lower_bound(spline.begin(),spline.end(),pos,[](const point& a,double b) { return a.first<b; });

   if (dump) cerr << "(" << (iter-1)->first << "," << (iter-1)->second << ") - " << pos << " - (" << (iter+0)->first << "," << (iter+0)->second << ")" << endl;
   
   if (iter->first == pos)
      return iter->second;

   double dx = iter->first - (iter - 1)->first;
   double dy = iter->second - (iter - 1)->second;

   double ofs = pos - (iter - 1)->first;
   return (iter - 1)->second + ofs * (dy / dx);
}

/*
static void convertBack(vector<point>& func) {
	for (auto& e: func)
		e.second -= delta;
}
*/

static vector<point> compressFunc(vector<point>& func, unsigned desiredSize)
// Compress to the desired size
{
	// Relax a bit to speed up compression
	unsigned maxSize = desiredSize + (desiredSize / 100), minSize = desiredSize - (desiredSize / 100);

	cerr << maxSize << " and min = " << minSize << endl;

	// Fits?
	if (func.size()<=maxSize)
		return func;

	// No, binary search
	long long capacity = 1e9, left = 1, right = capacity;
	double mul = 1.0 / capacity;
	//long double last = 0;
	while (left < right) {
		long long middle=(left+right)/2;
	

		double epsilon = mul * middle;
				cerr << "------------- Binary search for " << epsilon << endl;

		//if (abs(epsilon - last) < doublePrecision)return 

		if (myabs(epsilon - 1) < doublePrecision) {
			return samplePointsOfApproximateFunction(func, epsilon);
		}
		if (myabs(epsilon - capacity) < doublePrecision) {
			return samplePointsOfApproximateFunction(func, epsilon);
		}

		vector<point> candidate = samplePointsOfApproximateFunction(func, epsilon);
		
		cerr << "done sample" << endl;
		cerr << "------------- Binary search for " << epsilon << " has " << candidate.size() << endl;

		cerr << candidate.size() << " " << minSize << " " << maxSize << endl;

		if (candidate.size() < minSize) {
	 		right = middle;
		} else if (candidate.size()>maxSize) {
	 		left = middle;
		} else {
			cerr << "Found already in binary search = " << middle << " and epsilon " << mul * middle << endl; 
			//convertBack(candidate);
	 		return candidate;
		}
	}

	cerr << "At the end is left = " << left << " and so epsilon " << (mul * left) << endl;

   // Final call, this is the best we could get
   vector<point> candidate = samplePointsOfApproximateFunction(func, mul * left);
   //convertBack(candidate);
   return candidate;
}

int mainn(int argc, char** argv) {
	if (argc != 7) {
		cout << "Usage: ./a.out inputVariant(1 -> binary, 0 -> txt) inputFileName outputVariant(1 -> binary, 0 -> txt) outputFileName howManySupportingPoints howManyToRead\n";
		return 1;
	}
	int inputVariant = atoi(argv[1]);
	char* inputFileName = argv[2];
	int outputVariant = atoi(argv[3]);
	char* outputFileName = argv[4];
	int countSupportingPoints = atoi(argv[5]);
	int countToRead = atoi(argv[6]);

	cerr << "Begin of reading" << endl;
	auto readingStart = high_resolution_clock::now();

	vector<point> coord;
	if (inputVariant == 1)
		coord = readBinary(inputFileName, countToRead);
	else 
		coord = readTxt(inputFileName);

	auto readingStop = high_resolution_clock::now();
	cerr << "Reading done!" << endl;
	/*
	if (useSort) {
		cerr << "Begin of sort" << endl;
		cerr << "wieder" << coord.size() << endl;
		auto startSort = high_resolution_clock::now();

		sort(coord.begin(), coord.end());

		auto sortStop = high_resolution_clock::now();
		cerr << "Sort done!" << endl;

		cerr << "wieder" << coord.size() << endl;

		//writeBinary(outputFileName, coord);
		//return 0;
	}
	*/
	cerr << "test" << coord[0].first << " " << coord[0].second << endl;
	cerr << "test" << coord[1].first << " " << coord[1].second << endl;

	cerr << "Begin of the algorithm" << endl;
	auto startAlg = high_resolution_clock::now();

	cerr << "count = " << countSupportingPoints << endl;
/*
	unsigned n = coord.size();
	points = new point*[n + 1];
	for (unsigned index = 1; index <= n; ++index) {
		points[index] = new point[2];
	}	
*/
	//vector<point> samplePoints = compressFunc(coord, countSupportingPoints);
	vector<point> samplePoints = samplePointsOfApproximateFunction(coord, 1e-15);

	cerr << "End of the algorithm" << endl;
	auto algStop = high_resolution_clock::now();

	cerr << "Begin of writing" << endl;
	auto startWriting = high_resolution_clock::now();
/*
	if (outputVariant == 1)
		writeBinary(outputFileName, samplePoints);
	else
		writeTxt(outputFileName, samplePoints);
	
	auto writingDone = high_resolution_clock::now();
*/
	if (samplePoints.size() == 0) {
		cerr << "alles falsch" << endl;
		return 0;
	}
	  double sumError = 0, maxError = 0;
	  for (auto e: coord) {
          // TODO: I have put here a + delta! Don't forget that!
	     double estimate = interpolate(samplePoints, e.first, false);
	     double real = e.second + delta;
	     double diff = estimate-real;
	     if (diff<0) 
	     	diff=-diff;
	     if (diff>maxError)
	        maxError=diff;
	     sumError+=diff;
	  }
	  cerr << samplePoints.size() << " " << maxError << " " << (sumError / coord.size()) << endl;
	
	//auto duration = duration_cast<milliseconds>(writingDone - readingStart);
	//cerr << "Total time is = " << duration.count() << endl;
	return 0;
}
