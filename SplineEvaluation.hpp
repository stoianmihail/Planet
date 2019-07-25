#include <vector>
#include <boost/math/tools/polynomial.hpp>
#include <fstream>
#include <immintrin.h>
#include <map>
#include "../../Btree/btree_map.h"
//---------------------------------------------------------------------------
typedef std::pair<double,double> Coord;
typedef std::pair<double, unsigned> segCoord;
typedef btree::btree_map<double, double> btree_map;
typedef btree::btree_map<double, unsigned> btree_map_segments;
//---------------------------------------------------------------------------
class SplineEvaluation;
//---------------------------------------------------------------------------
/// A linear spline for evaluation
class SplineEvaluation 
{
    private:
    static const unsigned countOfBuckets = 2 + 1;
    static const unsigned limitSplineSize = 100000;
    int32_t offset;
    int32_t length;
    
    void transformSpline();
    double selfInterpolate(unsigned segment, double& x) const;
    unsigned searchLeft(int segment, double& x) const;
    unsigned searchRight(int segment, double& x) const;
    unsigned binarySearch(unsigned lower, unsigned upper, double& x) const;
    void changeIntoArray(long double** moveHere, boost::math::tools::polynomial<long double> poly);
    void convertValues();
    void splitSegmentSpline();
    public:
    static const unsigned polyGrade = 2;
    /// Constructor
    SplineEvaluation();
    /// Constructor with the function
    SplineEvaluation(const std::vector<Coord>& function, unsigned desiredSize, const unsigned useSplineMode, const unsigned searchMode, const unsigned testMode);
    /// Constructor with given data
    SplineEvaluation(unsigned size, const double* x, const double* y);

    void save(std::ofstream& file_out);
    void load(std::ifstream& file_in);
    
    unsigned size() const;
    unsigned findExactSegment(double segment, double& x) const;
    /// Return the approximate value of function(x)
    double evaluate(double& x) const;
    double hornerEvaluate(double& x) const;
    
    double splittedHornerEvaluate(double& x) const;
    double hornerEvaluateByBuckets(unsigned bucketIndex, double& x) const; 

    double chebyshevEvaluate(double& x) const;
    double splittedChebyshevEvaluate(double& x) const;
    
    double mapEvaluate(double& x) const;
    std::map<double, double> coordinates;
    
    double btreeEvaluate(double& x) const;
    btree_map* btreeMap;
    
    double btreeSegmentsEvaluate(double& x) const;
    btree_map_segments* btreeSegments;
    
    std::vector<Coord> spline;
    /// Spline which has to be fitted
    std::vector<segCoord> segmentSpline;
    /// Fitted polynomial to the spline
    long double* poly;
    double* values;
    unsigned splineSize;
    
    std::vector<segCoord>* splineBuckets;
    long double **polyBuckets;
    double splittedPositions[countOfBuckets];
    __m128d saveLimits;
    
    /// Fit a polynomial to spline
    boost::math::tools::polynomial<long double> fitSpline(const std::vector<segCoord>& spline);
};
//---------------------------------------------------------------------------

