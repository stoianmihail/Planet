#ifndef H_cts_infra_SplineHistogram
#define H_cts_infra_SplineHistogram
//---------------------------------------------------------------------------
#include <vector>
//---------------------------------------------------------------------------
class SplineHistogramBuilder;
class Package;
//---------------------------------------------------------------------------
/// A linear spline histogram
class SplineHistogram
{
   public:
   /// The maximum number of slots
   static const unsigned maxSlotCount = 200;
   /// The number of slots
   static unsigned slotCount;

   private:
   /// Coordinates of the supporting points
   double x[maxSlotCount],y[maxSlotCount];

   /// Find the appropriate slot. Result in [0,slotCount]!
   unsigned findSlot(double score) const;
   /// Find the appropriate slot. Result in [0,slotCount]!
   unsigned findSlotCount(double count) const;

   /// Add a score interval
   void addInterval(double from,double to,double value);

   friend class SplineHistogramBuilder;
   friend class Package;

   public:
   /// Constructor
   SplineHistogram();
   /// Constructor with given data
   SplineHistogram(const double* x,const double* y);

   /// Get a known histogram
   static SplineHistogram getScoreDistribution(const Package& package,const char* name);

   /// Total number of documents
   unsigned totalCount() const { return static_cast<unsigned>(y[slotCount-1]); }
   /// Get the number of documents with a score >= x
   double documentsAboveOrEqual(double score) const;
   /// Get the score after n documents. This is a lower bound!
   double scoreAfter(unsigned n) const;

   /// Combine two histograms (join)
   SplineHistogram combineJoinSlow(const SplineHistogram& other,double resultSizeIfKnown = -1.0) const;
   /// Combine two histograms (join)
   SplineHistogram combineJoinTree(const SplineHistogram& other,double resultSizeIfKnown = -1.0) const;
   /// Combine two histograms (join). Not as exact as the other algorithms but faster
   SplineHistogram combineJoin(const SplineHistogram& other,double resultSizeIfKnown = -1.0) const;
   /// Remove all entries below a certain threshold
   SplineHistogram cutBelow(double score) const;

   /// Get the coordinate (debug only)
   double getXCoord(unsigned slot) const { return x[slot]; }
   /// Get the coordinate (debug only)
   double getYCoord(unsigned slot) const { return y[slot]; }
};
//---------------------------------------------------------------------------
/// Constructs a spline histogram from a function with ascending x and y coordinates
class SplineHistogramBuilder
{
   public:
   /// Build it
   static SplineHistogram build(const std::vector<std::pair<double,double> >& func);
   /// Build it
   static SplineHistogram buildDP(const std::vector<std::pair<double,double> >& func);
   /// Build it
   static SplineHistogram buildDPQuad(const std::vector<std::pair<double,double> >& func);
};
//---------------------------------------------------------------------------
#endif

