#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <unordered_map>

using namespace std;

/// for precision reasons
const double MIN_ANGLE = 1e-16; // 10^(-precision(double))
const double MAX_RUNS = 6;

typedef std::pair<double,double> Coord;
const double precision = std::numeric_limits<double>::epsilon();
const double pi = acos(0);

//---------------------------------------------------------------------------
static inline double dev(const Coord& a,const Coord& b)
   // Calculate the derivation
{
   double dx=b.first-a.first;
   double dy=b.second-a.second;
   return dy/dx;
}
//---------------------------------------------------------------------------
static int cmpDevs(const Coord& a,const Coord& b,const Coord& c)
   // Compare derivations. -1 if (a->b) is more swallow than (a->c), +1 if it is more steeper
{
   double dx1=b.first-a.first,dx2=c.first-a.first;
   double dy1=b.second-a.second,dy2=c.second-a.second;
   if ((dy1*dx2)<(dy2*dx1)) return -1;
   if ((dy1*dx2)>(dy2*dx1)) return 1;
   return 0;
}
//---------------------------------------------------------------------------
static double interpolate(const std::vector<Coord>& spline,double pos, bool dump=false);
//---------------------------------------------------------------------------
static double square(double x) {
    return x * x;
}
//---------------------------------------------------------------------------
static int cmp(double x, double y) {
// x greater than y
    return (x - y) > ((fabs(x) < fabs(y) ? fabs(y) : fabs(x)) * precision);
}
//---------------------------------------------------------------------------
static int equals(double x, double y) {
    return fabs(x - y) < precision;
}
//---------------------------------------------------------------------------
static double distance(Coord a, Coord b) 
// returns euclidian distance between both points
{
    return sqrt(square(a.first - b.first) + square(a.second - b.second));
}
//---------------------------------------------------------------------------
static double cosinus(Coord u, Coord v, Coord w) 
// computes cos(angle) between (u, v, w) with v in middle
{
    double a = distance(u, w), b = distance(u, v), c = distance(v, w);
    return (square(a) - square(b) - square(c)) / (2.0 * b * c);
}
//---------------------------------------------------------------------------
static double computeAngle(Coord u, Coord v, Coord w) 
// computes angle between (u, v, w) using cosinus
{
    // TODO: compare with arccos, maybe it's better
    return 1.0 - cosinus(u, v, w);
}
//---------------------------------------------------------------------------
static double computeDistance(Coord u, Coord v, Coord w) 
// compute the distance which we deviated by
{
    double ipot = distance(v, w);
    double cosAngle = cosinus(u, v, w);
    double c1 = ipot * cosAngle;
    double c2 = ipot * sqrt(1 - square(cosAngle)); // sin * ipot
    // the height on (v, w) = c1 * c2 / ipot
    return c1;
}
//---------------------------------------------------------------------------
static double determinant(Coord u, Coord v, Coord w) 
// computes the determinant, which represents the spanned area.
{
    return (v.first - u.first) * (w.second - u.second) - (v.second - u.second) * (w.first - u.first);
}
//---------------------------------------------------------------------------
static int orientation(Coord u, Coord v, Coord w) 
// computes in which direction w is placed in respect to (u, v)
{
    double det = determinant(u, v, w);
    if (det < -precision)
        return -1;
    else if (det > precision)
        return 1;
    else
        return 0;
}
//---------------------------------------------------------------------------
#if 0
static double computeMaxDistance(const std::vector<Coord>& data) 
// compute the max distance between any 2 consecutive points
{
    double maxDistance = 0, minDeviation = 1.0, currDist = 0, currSum = 0;
    if (data.size() == 1)
        return maxDistance;
    std::vector<Coord>::const_iterator iter = data.begin(), limit = data.end();
    Coord prev = *iter, curr = *(++iter);
    for (++iter; iter != limit; ++iter) {
        double deviation = 0.0;
        double distance = 0.0;
        int orient = orientation(prev, curr, *iter);
        if (orient) {
            deviation = computeAngle(prev, curr, *iter);
            distance = computeDistance(prev, curr, *iter);
            currSum += orient * deviation;
            currDist += orient * distance;
        }
        if (deviation < minDeviation)
            minDeviation = deviation;
        
        if (fabs(distance) > maxDistance)
            maxDistance = fabs(distance);
        
        prev = curr;
        curr = *iter;
    }
    //cerr << "minDeviation= " << minDeviation << " average = " << currDist << endl;
    return maxDistance;
}
#endif
//---------------------------------------------------------------------------
static std::vector<Coord> angleApproximation(const std::vector<Coord>& data, double maxAngle, double maxDistance = 0)
// Compress the functions by storing only the inflections with angle lower than maxAngle
{
   // Add the first point
   std::vector<Coord> result;
   if (data.empty())
      return result;
   result.push_back(data[0]);
   if (data.size()==1)
      return result;

   // Taut the string
   std::vector<Coord>::const_iterator iter = data.begin(), limit = data.end();
   Coord prev = *iter, curr = *(++iter);
   double currSum = 0, currDist = 0;
   
   // it depends if the angle is smaller -> accept larger distances
   //maxDistance = (maxDistance) * (2 - maxAngle);
    for (++iter; iter != limit; ++iter) {
        double deviation = 0.0;
        double distance = 0.0;
        int orient = orientation(prev, curr, *iter);
        if (orient) {
            deviation = computeAngle(prev, curr, *iter);
            currSum += orient * deviation;
#if 0
            distance = computeDistance(prev, curr, *iter);
            currDist += orient * distance;
#endif
        }
#if 0
        if (fabs(distance) > ppmaxDistance)
            ppmaxDistance = fabs(distance);
#endif
        //cerr << " now is deviation = " << deviation << " and orient = " << orient << endl;
#if 0
        if (deviation > maxDeviation) {
            maxDeviation = deviation;
        }
        if (distance(prev, curr) > maxDistance)
            maxDistance = distance(prev, curr);
        if (currSum > maxSum)
            maxSum = currSum;
#endif
        // Check for inflexion.
        if ((cmp(deviation, maxAngle)) || (cmp(fabs(currSum), maxAngle))) { 
         // || (cmp(fabs(currDist), maxDistance))) {
                //cerr << "choose (" << curr.first << ", " << curr.second << ")" << endl;
                // cerr << "choose and distance = " << distance << endl;
                result.push_back(curr);
                currSum = 0;
        }
        
        // update the indexes.
        prev = curr;
        curr = *iter;
    }
    // We should have at least 2 points in result.
    // For that, add the last point.
    if (result.back() != *(--data.end())) {
        result.push_back(*(--data.end()));
    }
     //cerr << "maxDistance = " << ppmaxDistance << endl;
#if 0
    cerr << "for angle = " << maxAngle << endl;
    cerr << "maxDeviation = " << maxDeviation << endl;
    cerr << "maxDistance = " << maxDistance << endl;
    cerr << "maxSum = " << maxSum << endl;
#endif
    return result;
}

//---------------------------------------------------------------------------
static std::vector<Coord> newAngleApproximation(const std::vector<Coord>& data,double maxAngle)
{
   // Add the first point
   std::vector<Coord> result;
   if (data.empty())
      return result;
   result.push_back(data[0]);
   if (data.size()==1)
      return result;

   // Taut the string
   std::vector<Coord>::const_iterator iter=data.begin(),limit=data.end();
   Coord curr = *(++iter);
    double currSum = 0;
   
   //double relEpsilon = epsilon;
   for (++iter;iter!=limit;++iter) {
       Coord back = result.back(); 
       double deviation = 0.0;
        double distance = 0.0;
        int orient = orientation(back, curr, *iter);
        if (orient) {
            deviation = computeAngle(back, curr, *iter);
            currSum += orient * deviation;
#if 0
            distance = computeDistance(prev, curr, *iter);
            currDist += orient * distance;
#endif
        }
#if 0
        if (fabs(distance) > ppmaxDistance)
            ppmaxDistance = fabs(distance);
#endif
        //cerr << " now is deviation = " << deviation << " and orient = " << orient << endl;
#if 0
        if (deviation > maxDeviation) {
            maxDeviation = deviation;
        }
        if (distance(prev, curr) > maxDistance)
            maxDistance = distance(prev, curr);
        if (currSum > maxSum)
            maxSum = currSum;
#endif
        // Check for inflexion.
        if ((cmp(deviation, maxAngle)) || (cmp(fabs(currSum), maxAngle))) { 
         // || (cmp(fabs(currDist), maxDistance))) {
                //cerr << "choose (" << curr.first << ", " << curr.second << ")" << endl;
                // cerr << "choose and distance = " << distance << endl;
                result.push_back(curr);
                back = curr;
                currSum = 0;
        }
        // And remember the last point
        curr = *iter;
   }

    // We should have at least 2 points in result.
    // For that, add the last point.
    if (result.back() != *(--data.end())) {
        result.push_back(*(--data.end()));
    }
   return result;
}
//-------------------------------------------------------------------------
#if 1
//-------------------------------------------------------------------------
static double formula(double angle) {
    return angle;
}
//-------------------------------------------------------------------------
static std::vector<Coord> packFunc(const std::vector<Coord>& func,unsigned desiredSize)
   // Compress to the desired size
{
   // Relax a bit to speed up compression
   unsigned maxSize=desiredSize+(desiredSize/100),minSize=desiredSize-(desiredSize/100);

   // Fits?
   if (func.size()<=maxSize)
      return func;

   //for (int chooseAlgorithm = 0; chooseAlgorithm <= 1; chooseAlgorithm++) {
        // define interval for cos(x)
        const double inverseCapacity = (double) 1 / 1e16;
        long long left = 0, right = 1e16;
        
        int count = 0;
        static int chooseAlgorithm = false;
        // double maxDistance = computeMaxDistance(func);
        while (left < right) {
            long long mid = (left + right) / 2;
            double angle = formula(inverseCapacity * mid);
            std::vector<Coord> candidate = chooseAlgorithm ? angleApproximation(func, angle) : newAngleApproximation(func, angle);
            
            //cerr << "size was " << candidate.size() << " with angle = " << angle << endl;
            
            //if (++count == 4) exit(0);
            
            if (candidate.size() > maxSize) {
                // Too many points, there was a small angle -> increase angle
                left = mid + 1;
            } else if (candidate.size() < minSize) {
                // Too few points, there was a big angle (small cosinus) -> decrease angle(bigger cosinus)
                right = mid;
            } else {
                return candidate;
            }
        }

        cerr << "couldn't get better " << left << endl;
        
        // The best we could get
        double angle = inverseCapacity * left;
        std::vector<Coord> candidate = chooseAlgorithm ? angleApproximation(func, angle) : newAngleApproximation(func, angle);
        if (false) {
            return packFunc(candidate, desiredSize);
        } else {
            return candidate;
        }
   //}
   //return func;
}
#endif
//---------------------------------------------------------------------------
static double interpolate(const vector<Coord>& spline,double pos, bool dump)
{
   if (pos<=spline.front().first)
      return spline.front().second;
   if (pos>=spline.back().first)
      return spline.back().second;
   auto iter=lower_bound(spline.begin(),spline.end(),pos,[](const Coord& a,double b) { return a.first<b; });
   if (dump) cerr << "(" << (iter-1)->first << "," << (iter-1)->second << ") - " << pos << " - (" << (iter+0)->first << "," << (iter+0)->second << ")" << endl;
   if (iter->first==pos)
      return iter->second;

   double dx=(iter+0)->first-(iter-1)->first;
   double dy=(iter+0)->second-(iter-1)->second;

   double ofs=pos-(iter-1)->first;
   return (iter-1)->second+ofs*(dy/dx);
}
//---------------------------------------------------------------------------

static void print(const char* msg, Coord& u) {
    //cerr << msg << " " << u.first << " " << u.second << endl;
}
//---------------------------------------------------------------------------
static std::vector<Coord> tautString(const std::vector<Coord>& data,double maxValue,double epsilon)
{
   // Add the first point
   std::vector<Coord> result;
   if (data.empty())
      return result;
   result.push_back(data[0]);
   if (data.size()==1)
      return result;

   // Taut the string
   std::vector<Coord>::const_iterator iter=data.begin(),limit=data.end();
   Coord upperLimit,lowerLimit,last=result.back();

   //double relEpsilon = epsilon;
   for (++iter;iter!=limit;++iter) {
      // Add the new bounds
      Coord u=*iter,l=*iter,b=result.back();
      //epsilon=1+((maxValue-(*iter).second)*relEpsilon);
      // epsilon=1+maxValue*relEpsilon;
      u.second+=epsilon; l.second-=epsilon;
#if 0
        print("upperLimit = ", upperLimit);
        print("lowerLimt = ", lowerLimit);

      print("U = ", u);
      print("L = ", l);
      
      print("b = ", b);
      print("last = ", last);
#endif
      // Check if we cut the error corridor
      if ((last!=b)&&((cmpDevs(b,upperLimit,*iter)<0)||(cmpDevs(b,lowerLimit,*iter)>0))) {
         result.push_back(last);
         b=last;
            
            //cerr << "found a new sample point" << endl;
        }

      // Update the error margins
      if ((last==b)||(cmpDevs(b,upperLimit,u)>0)) {
         //cerr << "first If" << endl;
         upperLimit=u;
        }
      if ((last==b)||(cmpDevs(b,lowerLimit,l)<0)) {
        lowerLimit=l;
        //cerr << "second If" << endl;
       }

      // And remember the current point
      last=*iter;
   }

   // Add the last point
   result.push_back(*(--data.end()));

   return result;
}
//---------------------------------------------------------------------------
static std::vector<Coord> compressFunc(const std::vector<Coord>& func,unsigned desiredSize)
   // Compress to the desired size
{
   // Relax a bit to speed up compression
   unsigned maxSize=desiredSize+(desiredSize/100),minSize=desiredSize-(desiredSize/100);

   // Fits?
   if (func.size()<=maxSize)
      return func;

   //cerr << "max value = " << func.back().second << endl;
   
   unsigned eps;
   // No, binary search
   unsigned left=0,right=func.size();
   while (left<right) {
      unsigned middle=(left+right)/2;
      std::vector<Coord> candidate=tautString(func,func.back().second,middle);
      if (candidate.size()<minSize) {
         right=middle;
      } else if (candidate.size()>maxSize) {
         left=middle+1;
      } else {
#if 0
    eps = middle;
    goto endCompressFunc;
#else
    return candidate;
#endif
    }
   }
#if 0    
    cerr << "not yet!" << endl;
    // find a lower epsilon, which generates the same number of sampling points.
    eps = left;
    endCompressFunc: {
        double step = 0.1;
        std::vector<Coord> last = tautString(func, func.back().second, eps);
        while (eps != 0) {
            eps -= step;
            std::vector<Coord> candidate = tautString(func, func.back().second, eps);
            if (candidate.size() != last.size())
                return last;
            last = candidate;
        }
        return last;
    }
#else
   // Final call, this is the best we could get
   return tautString(func,func.back().second,left);
#endif    
}
//---------------------------------------------------------------------------
static pair<double, double> computeOverallPrecision(const std::vector<Coord>& func, const std::vector<Coord>& reduced)
// returns (maxError, avgError)
{
    // compute errors
    double sumError=0, maxError=0, relError = 0;
    for (auto e: func) {
        double estimate=interpolate(reduced, e.first);
        double real=e.second;
        double diff=estimate-real;
        if (diff<0) 
            diff=-diff;
        if (diff>maxError)
            maxError=diff;
        sumError+=diff;
        relError += (!real) ? diff : diff / real;
    }
    relError /= func.size();
    sumError /= func.size();
    return make_pair(maxError, sumError);
}
//---------------------------------------------------------------------------
static std::vector<Coord> preciseCompressFunc(const std::vector<Coord>& function, unsigned desiredSize)
// first use a sieve for the smallest angles, then apply the previous algorithm
{
    // Apply the last algorithm to compare its results with the new ones
    std::vector<Coord> candidate = compressFunc(function, desiredSize);
    pair<double, double> previousAlg = computeOverallPrecision(function, candidate);
    double previousMaxError = previousAlg.first, previousAvgError = previousAlg.second;

    cerr << "previous has " << previousMaxError << " and " << previousAvgError << endl;
#if 0
    // TODO: put back
    return candidate;
#endif
    // Map the sizes of the splines to pairs of (maxError, avgError)
    unordered_map<int, pair<int, pair<double, double>>> answers;
    answers[candidate.size()] = make_pair(candidate.size(), previousAlg);
#if 0
    // TODO: get a good angle!
    double init = inflexionAngleOfFunction(function);
    double expInit = log10(init);
    unsigned expDiff = (unsigned)(expInit - log10(precision));
    assert(0 <= expDiff && expDiff <= 16);
    
    cerr << "expINit = " << expInit << " diff = " << expDiff << endl; 
    
    // reset init (this part could be changed, because the angle which we start with depends on how the function looks like)
    init = 1e-1;
    int count = (int)-expInit;
    while (count--)
        init /= 10;
    
    cerr << "after reset init = " << init << endl;
#else
    // it could be modified, for precision reasons
    double currAngle = MIN_ANGLE;
    bool foundBetter = true;
    int runs = MAX_RUNS;
#endif
    double smallestAvgError = previousAvgError, smallestMaxError = previousMaxError, cand = -1;
    for (int runs = MAX_RUNS; (runs--) && (foundBetter); currAngle *= 10) {
        foundBetter = false;
        
        cerr << "progress = " << runs << " with " << currAngle << endl;
        for (unsigned mul = 2; mul < 20; ++mul) {
            double eps = currAngle * mul / 2;
            
            // Sharp the function then apply the previous algorithm
            auto samplePoints = angleApproximation(function, eps);
            samplePoints = compressFunc(samplePoints, desiredSize);
            
            pair<double, double> errors = computeOverallPrecision(function, samplePoints);
            double maxError = errors.first, avgError = errors.second;
            
            cerr << "check with " << eps << " errors = " << maxError << " " << avgError << endl;
            
            if ((previousAvgError > avgError) && (previousMaxError > maxError)) {   
#if 1
                    cerr << "count Points : " << samplePoints.size() << endl;
                    cerr << "it's better at " << eps << " with " << maxError << " and " << avgError << endl;
                    cerr << "rate for avg is " << ((-avgError + previousAvgError) / previousAvgError) * 100 << endl;
                    cerr << "rate for max is " << ((-maxError + previousMaxError) / previousMaxError) * 100 << endl;
#endif
                
#if 0
                // Check the existance of this size
                if (answers.find(samplePoints.size()) == answers.end()) {
                    // not yet in map this size
                    candidate = compressFunction(function, samplePoints.size());
                    pair<double, double> curr = computeOverallPrecision(function, candidate);
                    
                    // And map it to its errors
                    answers[samplePoints.size()] = make_pair(candidate.size(), curr);
                }
                
                auto curr = answers.find(samplePoints.size())->second;
                
                cerr << "to " << samplePoints.size() << " found : " << curr.first << " with " << "(" << curr.second.first << "," << curr.second.second << ")" << endl;
#endif
                candidate = compressFunc(function, samplePoints.size());
                pair<double, double> curr = computeOverallPrecision(function, candidate);
            
                if (curr == previousAlg) {
                    cerr << "better " << endl;
                    foundBetter = true;
                    
                    // get the best avgError
                    if (smallestAvgError > avgError) {
                        smallestAvgError = avgError;
                        smallestMaxError = maxError;
                        cand = eps;
                    } else if (equals(smallestAvgError, avgError)) {
                        if (smallestMaxError > maxError) {
                            smallestMaxError = maxError;
                            cand = eps;
                        }
                    }
                } else {
                    if (candidate.size() <= samplePoints.size())
                        cerr << "Has been beaten at " << samplePoints.size() << endl;
                }
            }
        }
    }
    if (cand == -1) {
        cerr << "Not found better than that!" << endl;
        return candidate;
    }
    
    // Get the best through the sieve of angles
    auto ret = angleApproximation(function, cand);
    ret = compressFunc(ret, desiredSize);
    return ret;
}
#if 0
//---------------------------------------------------------------------------
static std::vector<Coord> preciseCompressFunc(const std::vector<Coord>& func,unsigned desiredSize) 
// first use a sieve for the smallest angles, then apply the previous algorithm
{
    std::vector<Coord> candidate = compressFunc(func, desiredSize);
    pair<double, double> previousAlg = computeOverallPrecision(func, candidate);
    double previousMaxError = previousAlg.first, previousAvgError = previousAlg.second;
    cerr << candidate.size() << "previous has max = " << previousMaxError << " and avg = " << previousAvgError << endl; 
    
    map<int, pair<double, double> > answers;
    
    answers[candidate.size()] = previousAlg;
    
    double init = 1e-14;
    double smallestAvgError = previousAvgError, smallestMaxError = previousMaxError, cand = -1;
    for (unsigned power = 1; power <= 4; power++) {
        init /= 10;
        cerr << "progress = " << power << endl;
        for (unsigned mul = 1; mul < 10; ++mul) {
            double eps = init * mul;
            auto ts = angleApproximation(func, eps);
            ts = compressFunc(ts, desiredSize);
            
            pair<double, double> errors = computeOverallPrecision(func, ts);
            double maxError = errors.first, avgError = errors.second;
            
            // TODO: get the best solution.
            if (cmp(previousAvgError, avgError) && cmp(previousMaxError, maxError)) {
#if 1
                cerr << "count Points : " << ts.size() << endl;
                cerr << "it's better at " << eps << " with " << maxError << " and " << avgError << endl;
                cerr << "rate for avg is " << ((-avgError + previousAvgError) / previousAvgError) * 100 << endl;
                cerr << "rate for max is " << ((-maxError + previousMaxError) / previousMaxError) * 100 << endl;
#endif      
                // Check the existance of this size
                if (answers.find(ts.size()) == answers.end()) {
                    // not yet in map this size
                    candidate = compressFunc(func, ts.size());
                    pair<double, double> curr = computeOverallPrecision(func, candidate);
                    
                    cerr << "Checked with " << candidate.size() << endl;
                    
                    cerr << "It has " << curr.second << " vs " << avgError << " and " << curr.first << " vs " << maxError << endl;
                    answers[ts.size()] = curr;
                }
                
                cerr << "Check me" << endl;
                
                auto curr = answers.find(ts.size())->second;
                if (cmp(curr.second, avgError) && cmp(curr.first, maxError)) {
                    // get the best avgError
                    if (cmp(smallestAvgError, avgError)) {
                        smallestAvgError = avgError;
                        smallestMaxError = maxError;
                        cand = eps;
                    } else if (equals(smallestAvgError, avgError)) {
                        if (cmp(smallestMaxError, maxError)) {
                            smallestMaxError = maxError;
                            cand = eps;
                        }
                    }
                } else {
                    cerr << "Has been beaten at " << ts.size() << endl;
                }
            }
        }
    }
    if (cand == -1) {
        cerr << "Not found better than that!" << endl;
        exit(0);
    }
    auto ts = angleApproximation(func, cand);
    ts = compressFunc(ts, desiredSize);
    return ts;
#if 0
    std::vector<Coord>::const_iterator ptr = candidate.begin(), iter = func.begin(), last = func.begin();
    ++ptr;
    for (++iter; iter != func.end(); ++iter) {
        if (*iter == *ptr && *last == *(ptr - 1)) {
            cerr << "found sth" << endl;
        }
        last = iter;
    }
    return candidate;

#endif
}
#endif
//---------------------------------------------------------------------------
int main(int argc,char* argv[])
{
   if (argc<4) {
      cerr << "usage: " << argv[0] << " file samplingPoints countToRead whichAlgoritm (0 -> worse, 1 -> better)" << endl;
      return 1;
   }
#if 0
   vector<double> file;
   {
      ifstream in(argv[1]);
      if (!in.is_open()) {
         cerr << "unable to open " << argv[1] << endl;
         return 1;
      }
#if 0
      double v;
      while (in>>v)
         file.push_back(v);
#else
      auto pos = in.tellg();
      in.seekg( 0, ios::end );
      auto size=in.tellg()-pos;
      in.seekg(0,ios::beg);

      unsigned elements=size/sizeof(double);
      file.resize(elements);
      in.read(reinterpret_cast<char*>(file.data()),elements*sizeof(double));
#endif
   }
   sort(file.begin(),file.end());
   vector<pair<double,double>> cdf;
   {
      unsigned pos=0;
      double last=file.front()+1;
      for (auto d:file) {
         if (d!=last)
            cdf.push_back({d,pos});
         pos++;
         last=d;
      }
      cdf.push_back({file.back()+1,pos});
   } vector<pair<double,double>> cdf;
   {
      ifstream in(argv[1]);
      auto pos = in.tellg();
      in.seekg( 0, ios::end );
      auto size=in.tellg()-pos;
      in.seekg(0,ios::beg);

      unsigned elements=size/sizeof(pair<double,double>);
      cdf.resize(elements);
      in.read(reinterpret_cast<char*>(cdf.data()),elements*sizeof(pair<double,double>));
   }

   if (false) {
      ofstream out("cfg");
      for (auto& e:cdf)
         out << e.first << " " << e.second << "\n";
   }
#else

#if 1
   unsigned countToRead = atoi(argv[3]);
   vector<pair<double,double>> cdf;
   {
      ifstream in(argv[1]);
      auto pos = in.tellg();
      in.seekg( 0, ios::end );  
      auto size=in.tellg()-pos;
      in.seekg(0,ios::beg);

      
      unsigned elements=size/sizeof(pair<double,double>);
      if (countToRead != -1 && countToRead < elements) {
        elements = countToRead;
      }
      
      cerr << "total " << elements << endl;
      
      cdf.resize(elements);
      in.read(reinterpret_cast<char*>(cdf.data()),elements*sizeof(pair<double,double>));
#if 0
    in.read(reinterpret_cast<char*>(cdf.data()),elements*sizeof(pair<double,double>));
    in.read(reinterpret_cast<char*>(cdf.data()),elements*sizeof(pair<double,double>));
#endif
   }
#else
    unsigned countToRead = atoi(argv[3]);
    ifstream in(argv[1]);
    
    unsigned n;
    in >> n;
    vector<pair<double, double>> cdf(n);
    
    for (unsigned index = 0; index < n; ++index) {
        double x, y;
        in >> x >> y;
        cdf[index] = make_pair(x, y);
    }
#endif

#endif


#if 0
   for (unsigned error=10;error<=600000;) {
      auto ts=tautString(cdf,cdf.back().second,error);
      auto total=ts.size();
      unsigned parts=0;
      double maxError=0,w=0,r=0;
      for (auto e:cdf) {
         double estimate=interpolate(ts,e.first);
         double real=e.second;
         double diff=estimate-real;
         if (diff<0) diff=-diff;
         if (diff>maxError) {
            maxError=diff;
            w=e.first; r=e.second;
         }
      }
      cerr << error << " " << total << " " << ((8*2*total)/1024) << " " << maxError << endl;
      cerr << "worst: at " << w << " estimate " << interpolate(ts,w) << " real " << r << endl;
      // interpolate(ts,w,true);
      if (error<100) error+=10; else if (error<1000) error+=100; else if (error<10000) error+=1000; else error+=10000;
   }
#else
 
   //for (unsigned size=100;size<=100;) {
      unsigned size = atoi(argv[2]);
   	//auto ts=compressFunc(cdf,size);
#if 0
      vector<Coord> samples;
      samples.push_back(make_pair(0, 0));
      samples.push_back(make_pair(2, 1));
      samples.push_back(make_pair(4.6, 4.9));
      samples.push_back(make_pair(9, 6));
      samples.push_back(make_pair(11.4, 8.4));
      samples.push_back(make_pair(23, 11));
      auto ts = samples;//tautString(cdf, cdf.back().second, 1);
#else
#if 0
   cerr << "precision = " << std::setprecision(12) << precision << endl;
   auto ts = angleApproximation(cdf, 3 * 1e-15);
   cerr << "my part done!" << endl;
   cerr << "remain with " << ts.size() << endl;
    ts = sieveCompressFunc(ts, size);
#endif

#endif
    int mode = atoi(argv[4]);
    cerr << mode << " 0 -> packFunc, 1 -> precise" << endl;  
    auto ts = mode ? preciseCompressFunc(cdf, size) : packFunc(cdf, size);
    pair<double, double> answer = computeOverallPrecision(cdf, ts);
        
    cerr << "mode = " << mode << "with comp = " << ts.size() << " maxError = " << answer.first << " and avgError = " << answer.second << endl; 
    
    
#if 0
        cerr << "Sample points of size " << ts.size() << endl;
        for (auto e: ts) {
            cerr << e.first << " " << e.second << endl;
        }
#endif
     // if (size<1000) size+=100; else if (size<10000) size+=1000; else size+=10000;
   //}
#endif
}
//---------------------------------------------------------------------------

