#include <vector>
#include <boost/math/tools/polynomial.hpp>
#include <fstream>
//---------------------------------------------------------------------------
typedef std::pair<double,double> Coord;
typedef std::pair<double, unsigned> segCoord;
//---------------------------------------------------------------------------
class SplineEvaluation;
//---------------------------------------------------------------------------
/// A linear spline for evaluation
class SplineEvaluation 
{
    private:
    unsigned limitSplineSize = 100000;
    int32_t offset;
    int32_t length;
    
    void transformSpline();
    double selfInterpolate(unsigned segment, double& x) const;
    unsigned searchLeft(int segment, double& x) const;
    unsigned searchRight(int segment, double& x) const;
    unsigned binarySearch(unsigned lower, unsigned upper, double& x) const;
    void changeIntoArray(boost::math::tools::polynomial<long double> poly);
    void convertValues();
    public:
    unsigned polySize = 2;
    /// Constructor
    SplineEvaluation();
    /// Constructor with the function
    SplineEvaluation(const std::vector<Coord>& function, unsigned desiredSize, const unsigned mode);
    /// Constructor with given data
    SplineEvaluation(unsigned size, const double* x, const double* y);

    void save(std::ofstream& file_out);
    void load(std::ifstream& file_in);
    
    unsigned size() const;
    unsigned findExactSegment(double segment, double& x) const;
    /// Return the approximate value of function(x)
    double evaluate(double& x) const;
    double hornerEvaluate(double& x) const;
    double chebyshevEvaluate(double& x) const;
    
    std::vector<Coord> spline;
    /// Spline which has to be fitted
    std::vector<segCoord> segmentSpline;
    /// Fitted polynomial to the spline
    long double* poly;
    double* values;
    unsigned splineSize;
    
    /// Fit a polynomial to spline
    boost::math::tools::polynomial<long double> fitSpline(const std::vector<segCoord>& spline);
};
//---------------------------------------------------------------------------

