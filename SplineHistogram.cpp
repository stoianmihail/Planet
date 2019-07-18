#include "SplineHistogram.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <cmath>

using namespace std;
//---------------------------------------------------------------------------
// Number of slots
unsigned SplineHistogram::slotCount = 100;
//---------------------------------------------------------------------------
namespace {
//---------------------------------------------------------------------------
typedef std::pair<double,double> Coord;
//---------------------------------------------------------------------------
class AggregateIntervalTree
{
   private:
   /// A tree node
   struct Node
   {
      /// Possible colors
      enum Color { Red, Black };

      /// The tree structure
      Node*   parent,*left,*right;
      /// The payload
      double  from,to,value;
      /// The color
      Color color;

      /// Constructor
      Node(double from,double to,double value);
      /// Destructor
      ~Node();
   };

   /// The root
   Node* root;

   /// Insert a node
   void insert(Node* newOne);
   /// Find a node
   Node* searchNode (double to) const;
   /// Find a node
   Node* searchNodeOrNeighbour(double to) const;
   /// Rotate left
   void leftRotate(Node* node);
   /// Rotate right
   void rightRotate(Node* node);

   /// Find the minimum child
   static const Node* minimum(const Node* node);
   /// Find the successor
   static const Node* successor(const Node* node);
   /// Find the minimum child
   static Node* minimum(Node* node);
   /// Find the successor
   static Node* successor(Node* node);
   /// Find the maximum child
   static const Node* maximum(const Node* node);
   /// Find the predecessor
   static const Node* predecessor(const Node* node);
   /// Find the maximum child
   static Node* maximum(Node* node);
   /// Find the predecessor
   static Node* predecessor(Node* node);

   public:
   class const_iterator;
   friend class const_iterator;
   ///  An iterator
   class const_iterator {
      private:
      /// The node
      Node* node;

      public:
      /// Constructor
      const_iterator(Node* node) : node(node) {}

      /// Comparison
      bool operator==(const const_iterator& o) { return node==o.node; }
      /// Comparison
      bool operator!=(const const_iterator& o) { return node!=o.node; }
      /// Increment
      const_iterator& operator++() { node=successor(node); return *this; }

      /// From
      double from() const { return node->from; }
      /// To
      double to() const { return node->to; }
      /// Value
      double value() const { return node->value; }
   };

   /// Constructor
   AggregateIntervalTree();
   /// Destructor
   ~AggregateIntervalTree();

   /// Insert an interval
   void insert(double from,double to,double value);
   /// Lookup a single value
   bool lookup(double point,double& value) const;
   /// Clear the tree
   void clear() { delete root; root=0; }
   /// Is the tree empty?
   bool empty() const { return !root; };

   /// The first value
   const_iterator begin() const { return const_iterator(root?minimum(root):0); }
   /// Behind the last value
   const_iterator end() const { return const_iterator(0); }

   /// Get the function
   void getFunction(std::vector<Coord>& func) const;
};
//---------------------------------------------------------------------------
AggregateIntervalTree::Node::Node(double from,double to,double value)
   : parent(0),left(0),right(0),from(from),to(to),value(value),color(Black)
   // Constructor
{
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node::~Node()
   // Destructor
{
   delete left;
   delete right;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::AggregateIntervalTree()
   : root(0)
   // Constructor
{
}
//---------------------------------------------------------------------------
AggregateIntervalTree::~AggregateIntervalTree()
   // Destructor
{
   delete root;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node* AggregateIntervalTree::searchNode(double to) const
   // Find a node
{
   for (Node* iter=root;iter;) {
      if ((to<=iter->to)&&((to>iter->from)||((to==iter->from)&&(to==iter->to))))
         return iter;
      if (to<=iter->from)
         iter=iter->left; else
         iter=iter->right;
   }
   return 0;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node* AggregateIntervalTree::searchNodeOrNeighbour(double to) const
   // Find a node
{
   for (Node* iter=root;iter;) {
      if ((to<=iter->to)&&((to>iter->from)||((to==iter->from)&&(to==iter->to))))
         return iter;
      if (to<=iter->from) {
	 if (!iter->left)
	    return iter;
	 iter=iter->left;
      } else {
	 if (!iter->right)
	    return iter;
	 iter=iter->right;
      }
   }
   return 0;
}
//---------------------------------------------------------------------------
const AggregateIntervalTree::Node* AggregateIntervalTree::minimum(const Node* node)
   // Find the minimum
{
   while (node->left)
      node=node->left;
   return node;
}
//---------------------------------------------------------------------------
const AggregateIntervalTree::Node* AggregateIntervalTree::successor(const Node* node)
   // Find the successor
{
   if (node->right)
      return minimum(node->right);
   Node* iter=node->parent;
   while ((iter)&&(node==iter->right)) {
      node=iter;
      iter=iter->parent;
   }
   return iter;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node* AggregateIntervalTree::minimum(Node* node)
   // Find the minimum
{
   while (node->left)
      node=node->left;
   return node;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node* AggregateIntervalTree::successor(Node* node)
   // Find the successor
{
   if (node->right)
      return minimum(node->right);
   Node* iter=node->parent;
   while ((iter)&&(node==iter->right)) {
      node=iter;
      iter=iter->parent;
   }
   return iter;
}
//---------------------------------------------------------------------------
const AggregateIntervalTree::Node* AggregateIntervalTree::maximum(const Node* node)
   // Find the maximum
{
   while (node->right)
      node=node->right;
   return node;
}
//---------------------------------------------------------------------------
const AggregateIntervalTree::Node* AggregateIntervalTree::predecessor(const Node* node)
   // Find the predecessor
{
   if (node->left)
      return maximum(node->left);
   Node* iter=node->parent;
   while ((iter)&&(node==iter->left)) {
      node=iter;
      iter=iter->parent;
   }
   return iter;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node* AggregateIntervalTree::maximum(Node* node)
   // Find the maximum
{
   while (node->right)
      node=node->right;
   return node;
}
//---------------------------------------------------------------------------
AggregateIntervalTree::Node* AggregateIntervalTree::predecessor(Node* node)
   // Find the predecessor
{
   if (node->left)
      return maximum(node->left);
   Node* iter=node->parent;
   while ((iter)&&(node==iter->left)) {
      node=iter;
      iter=iter->parent;
   }
   return iter;
}
//---------------------------------------------------------------------------
void AggregateIntervalTree::insert(Node* n)
   // Insert a new node
{
   // Insert in the tree
   Node* iter1=root,*iter2=0;
   while (iter1) {
      iter2=iter1;
      if (n->to<=iter1->from)
         iter1=iter1->left; else
         iter1=iter1->right;
   }
   n->parent=iter2;
   if (!iter2) {
      root=n;
   } else {
      if (n->to<=iter2->from)
         iter2->left=n; else
         iter2->right=n;
   }

   // Rebalance
   iter1=n;
   iter1->color=Node::Red;
   while ((iter1!=root) && (iter1->parent->color==Node::Red)) {
      if (iter1->parent==root) break;
      if (iter1->parent==iter1->parent->parent->left) {
         iter2=iter1->parent->parent->right;
         if ((iter2) && (iter2->color==Node::Red)) {
            iter1->parent->color=Node::Black;
            iter2->color=Node::Black;
            iter1->parent->parent->color=Node::Red;
            iter1=iter1->parent->parent;
         } else  {
            if (iter1==iter1->parent->right) {
               iter1=iter1->parent; leftRotate(iter1);
            }
            iter1->parent->color=Node::Black;
            iter1->parent->parent->color=Node::Red;
            rightRotate(iter1->parent->parent);
         }
      } else {
         iter2=iter1->parent->parent->left;
         if ((iter2) && (iter2->color==Node::Red)) {
            iter1->parent->color=Node::Black;
            iter2->color=Node::Black;
            iter1->parent->parent->color=Node::Red;
            iter1=iter1->parent->parent;
         } else {
            if (iter1==iter1->parent->left) {
               iter1=iter1->parent; rightRotate(iter1);
            }
            iter1->parent->color=Node::Black;
            iter1->parent->parent->color=Node::Red;
            leftRotate(iter1->parent->parent);
         }
      }
   }
   root->color=Node::Black;
}
//---------------------------------------------------------------------------
void AggregateIntervalTree::insert(double from,double to,double value)
   // Inset an interval
{
   // Check bounds
   if (from>to) { double t=from; from=to; to=t; }
   if (value<=0.0) return;

   // A point?
   if (from==to) {
      Node* n=searchNode(to);
      // Check if a suitable range exists
      if (n) {
	 // A possible match?
	 if (n->to==to) {
	    if (n->from!=from) {
	       Node* n2=successor(n);
	       if (n2&&(n2->from==from)&&(n2->to==to))
	          n=n2;
	    }
	    if (n->from==from) {
	       n->value+=value;
	       return;
	    }
	 }
	 // Break the range down
	 if (n->to>to) {
	    double upper=n->to;
	    double remaining=n->value*(upper-to)/(upper-n->from);
	    n->value*=(to-n->from)/(upper-n->from);
	    n->to=to;
	    insert(new Node(to,upper,remaining));
	 }
	 if (n->from<to) {
	    double lower=n->from;
	    double remaining=n->value;
	    n->value=0;
	    n->from=to;
	    insert(new Node(lower,to,remaining));
	 }
	 n->value+=value;
      } else {
	 insert(new Node(from,to,value));
      }
      return;
   }

   // Lookup a node
   Node* n=searchNodeOrNeighbour(to);
   if (n&&(n->to<to)) {
      Node* n2=successor(n);
      if (n2&&(n2->from<to)) n=n2;
   }
   // And insert in chunks
   static const double minValue = 1.0E-200;
   while ((from!=to)&&(value>=minValue)) {
      if (!n) {
	 insert(new Node(from,to,value));
	 return;
      }
      if (n->to<to) {
	 double upper=to;
	 double remaining=value*(upper-n->to)/(upper-from);
	 value*=(n->to-from)/(upper-from);
	 to=n->to;
	 insert(new Node(n->to,upper,remaining));
      }
      if (from<n->from) {
	 if (n->to!=n->from) {
	    double upper=to;
	    double remaining=value*(upper-n->from)/(upper-from);
	    value*=(n->from-from)/(upper-from);
	    to=n->from;
	    n->value+=remaining;
	 }
	 n=predecessor(n);
      } else if (from>n->from) {
	 if (n->from==n->to) {
	    n=successor(n);
	    if (!n) break;
	 }
	 double lower=n->from;
	 double remaining=n->value*(from-lower)/(n->to-lower);
	 n->value*=(n->to-from)/(n->to-lower);
	 n->value+=value;
	 n->from=from;
	 insert(new Node(lower,from,remaining));
	 break;
      } else {
	 n->value+=value;
	 break;
      }
   }
}
//---------------------------------------------------------------------------
bool AggregateIntervalTree::lookup(double point,double& value) const
   // Lookup a value
{
   Node* iter=searchNode(point);
   if (!iter)
      return false;
   value=iter->value;
   return true;
}
//---------------------------------------------------------------------------
void AggregateIntervalTree::leftRotate(Node* node)
   // Rotate left
{
   Node* iter=node->right;
   if (!iter) return;
   node->right=iter->left;
   if (iter->left)
      iter->left->parent=node;
   iter->parent=node->parent;
   if (!node->parent) {
      root=iter;
   } else {
      if (node==node->parent->left)
         node->parent->left=iter; else
         node->parent->right=iter;
   }
   iter->left=node;
   node->parent=iter;
}
//---------------------------------------------------------------------------
void AggregateIntervalTree::rightRotate(Node* node)
   // Rotate right
{
   Node* iter=node->left;
   if (!iter) return;
   node->left=iter->right;
   if (iter->right)
      iter->right->parent=node;
   iter->parent=node->parent;
   if (!node->parent) {
      root=iter;
   } else {
      if (node==node->parent->left)
         node->parent->left=iter; else
         node->parent->right=iter;
   }
   iter->right=node;
   node->parent=iter;
}
//---------------------------------------------------------------------------
void AggregateIntervalTree::getFunction(std::vector<Coord>& func) const
   // Get the function
{
   func.clear();
   const_iterator iter=begin(),limit=end();
   func.push_back(std::pair<double,double>(iter.from(),0));

   double sum=0;
   for (;iter!=limit;++iter) {
      sum+=iter.value();
      func.push_back(std::pair<double,double>(iter.to(),sum));
   }
}
//---------------------------------------------------------------------------
}
//---------------------------------------------------------------------------
SplineHistogram::SplineHistogram()
  // Constructor
{
   memset(x,0,sizeof(x));
   memset(y,0,sizeof(y));
}
//---------------------------------------------------------------------------
SplineHistogram::SplineHistogram(const double* x,const double* y)
  // Constructor
{
   for (unsigned index=0;index<slotCount;index++) {
      this->x[index]=x[index];
      this->y[index]=y[index];
   }
}
//---------------------------------------------------------------------------
unsigned SplineHistogram::findSlot(double score) const
   // Find the appropriate slot. Result in [0,slotCount]!
{
   unsigned left=0,right=slotCount;
   while (left!=right) {
      unsigned middle=(left+right)/2;
      if (score>x[middle]) {
	 left=middle+1;
      } else if ((middle>0)&&(score<=x[middle-1])) {
	 right=middle;
      } else {
	 return middle;
      }
   }
   return left;
}
//---------------------------------------------------------------------------
unsigned SplineHistogram::findSlotCount(double count) const
   // Find the appropriate slot. Result in [0,slotCount]!
{
   unsigned left=0,right=slotCount;
   while (left!=right) {
      unsigned middle=(left+right)/2;
      if (count>y[middle]) {
	 left=middle+1;
      } else if ((middle>0)&&(count<=y[middle-1])) {
	 right=middle;
      } else {
	 return middle;
      }
   }
   return left;
}
//---------------------------------------------------------------------------
double SplineHistogram::documentsAboveOrEqual(double score) const
   // Get the number of documents with a score >= x
{
   unsigned slot=findSlot(score);
   if (slot>=slotCount) return 0;

   if (score==x[slot])
      return totalCount()-y[slot];

   if (!slot)
      return totalCount();

   double dx=x[slot]-x[slot-1];
   double dy=y[slot]-y[slot-1];
   double r=y[slot-1]+((score-x[slot-1])/dx)*dy;

   return totalCount()-r;
}
//---------------------------------------------------------------------------
double SplineHistogram::scoreAfter(unsigned n) const
   // Get the score after n documents
{
   if (n>=totalCount())
      return 0;
   unsigned count=totalCount()-n;

   unsigned slot=findSlotCount(count);
   if (slot>=slotCount)
      return 0;

   if ((n==y[slot])||(!slot))
      return x[slot];

   double dx=x[slot]-x[slot-1];
   double dy=y[slot]-y[slot-1];
   double r=x[slot-1]+((count-y[slot-1])/dy)*dx;

   return r;
}
//---------------------------------------------------------------------------
SplineHistogram SplineHistogram::combineJoinTree(const SplineHistogram& other,double resultSize) const
   // Combine two histograms (join)
{
   SplineHistogram result;

   // Calculate relative frequencies
   double a[slotCount-1];
   double b[slotCount-1];
   double div=totalCount();
   if (!div) div=1;
   for (unsigned index=0;index<slotCount-1;index++)
      a[index]=(y[index+1]-y[index])/div;
   a[0]=y[1]/div;
   div=other.totalCount();
   if (!div) div=1;
   for (unsigned index=0;index<slotCount-1;index++)
      b[index]=(other.y[index+1]-other.y[index])/div;
   b[0]=other.y[1]/div;

   // Estimate the matching tuples
   double intersectionPart,nonMatchingA,nonMatchingB;
   if (resultSize<0) {
      if (totalCount()>other.totalCount()) {
	 resultSize=totalCount();
	 intersectionPart=other.totalCount()/resultSize;
	 nonMatchingA=1-intersectionPart;
	 nonMatchingB=0;
      } else if (totalCount()<other.totalCount()) {
	 resultSize=other.totalCount();
	 intersectionPart=totalCount()/resultSize;
	 nonMatchingA=0;
	 nonMatchingB=1-intersectionPart;
      } else {
         resultSize=totalCount();
         intersectionPart=1;
         nonMatchingA=0;
         nonMatchingB=0;
      }
   } else {
      double intersection=totalCount()+other.totalCount()-resultSize;
      if (intersection<0) intersection=0;
      intersectionPart=intersection/resultSize;
      double div=totalCount()+other.totalCount(); if (!div) div=1;
      nonMatchingA=((resultSize-intersection)/resultSize)*totalCount()/div;
      nonMatchingB=((resultSize-intersection)/resultSize)*other.totalCount()/div;
   }

   // Construct the intervals
   AggregateIntervalTree tree;
   for (unsigned index=0;index<slotCount-1;index++)
      for (unsigned index2=0;index2<slotCount-1;index2++) {
	 double from=x[index]+other.x[index2];
	 double to=x[index+1]+other.x[index2+1];
	 tree.insert(from,to,a[index]*b[index2]*intersectionPart*resultSize);
      }
   for (unsigned index=0;index<slotCount-1;index++) {
      tree.insert(x[index],x[index+1],a[index]*nonMatchingA*resultSize);
      tree.insert(other.x[index],other.x[index+1],b[index]*nonMatchingB*resultSize);
   }

   // Construct a new spline
   std::vector<Coord> func;
   tree.getFunction(func);

   return SplineHistogramBuilder::build(func);
}
//---------------------------------------------------------------------------
namespace {
struct Interval {
   double from,to,value,dydx;

   Interval() {}
   Interval(double from,double to,double value) : from(from),to(to),value(value),dydx((to>from)?(value/(to-from)):0.0) {}

   /// Sort by to descending
   bool operator<(const Interval& i) const { return to>i.to; }
};
struct OrderIntervalFrom { bool operator()(const Interval& a,const Interval& b) const { return a.from<b.from; } };
}
//---------------------------------------------------------------------------
static inline void append(std::vector<Interval>& intervals,double from,double to,double value) { if (value>0.0) intervals.push_back(Interval(from,to,value)); }
//---------------------------------------------------------------------------
SplineHistogram SplineHistogram::combineJoinSlow(const SplineHistogram& other,double resultSize) const
   // Combine two histograms (join)
{
   SplineHistogram result;

   // Calculate relative frequencies
   double a[slotCount-1];
   double b[slotCount-1];
   double div=totalCount();
   if (!div) div=1;
   for (unsigned index=0;index<slotCount-1;index++)
      a[index]=(y[index+1]-y[index])/div;
   a[0]=y[1]/div;
   div=other.totalCount();
   if (!div) div=1;
   for (unsigned index=0;index<slotCount-1;index++)
      b[index]=(other.y[index+1]-other.y[index])/div;
   b[0]=other.y[1]/div;

   // Estimate the matching tuples
   double intersectionPart,nonMatchingA,nonMatchingB;
   if (resultSize<0) {
      if (totalCount()>other.totalCount()) {
	 resultSize=totalCount();
	 intersectionPart=other.totalCount()/resultSize;
	 nonMatchingA=1-intersectionPart;
	 nonMatchingB=0;
      } else if (totalCount()<other.totalCount()) {
	 resultSize=other.totalCount();
	 intersectionPart=totalCount()/resultSize;
	 nonMatchingA=0;
	 nonMatchingB=1-intersectionPart;
      } else {
         resultSize=totalCount();
         intersectionPart=1;
         nonMatchingA=0;
         nonMatchingB=0;
      }
   } else {
      double intersection=totalCount()+other.totalCount()-resultSize;
      if (intersection<0) intersection=0;
      intersectionPart=intersection/resultSize;
      double div=totalCount()+other.totalCount(); if (!div) div=1;
      nonMatchingA=((resultSize-intersection)/resultSize)*totalCount()/div;
      nonMatchingB=((resultSize-intersection)/resultSize)*other.totalCount()/div;
   }

   // Construct the intervals
   std::vector<Interval> intervals;
   for (unsigned index=0;index<slotCount-1;index++)
      for (unsigned index2=0;index2<slotCount-1;index2++) {
	 double from=x[index]+other.x[index2];
	 double to=x[index+1]+other.x[index2+1];
         append(intervals,from,to,a[index]*b[index2]*intersectionPart*resultSize);
      }
   for (unsigned index=0;index<slotCount-1;index++) {
      append(intervals,x[index],x[index+1],a[index]*nonMatchingA*resultSize);
      append(intervals,other.x[index],other.x[index+1],b[index]*nonMatchingB*resultSize);
   }
   std::sort(intervals.begin(),intervals.end(),OrderIntervalFrom());

   // Construct the function
   std::vector<Coord> func;
   std::vector<Interval> active;
   double scanPos=0;
   double pastValues=0;
   if (!intervals.empty()) {
      scanPos=intervals.front().from;
      func.push_back(Coord(scanPos,0));
   }
   for (std::vector<Interval>::const_iterator iter=intervals.begin(),limit=intervals.end();;++iter) {
      // Belongs to the current block?
      if ((iter!=limit)&&((*iter).from==scanPos)) {
         active.push_back(*iter);
         push_heap(active.begin(),active.end());
         continue;
      }
      // Find the next scan boundary
      double nextScan;
      if (iter==limit) {
         nextScan=scanPos;
         for (std::vector<Interval>::const_iterator iter2=active.begin(),limit2=active.end();iter2!=limit2;++iter2)
            if ((*iter2).to>nextScan)
               nextScan=(*iter2).to;
      } else nextScan=(*iter).from;

      // Flush redundant entries
      while (!active.empty()) {
         // Consider the next small interval
         double step=(active.front().to<nextScan)?active.front().to:nextScan;
         double sum=pastValues;
         for (std::vector<Interval>::const_iterator iter2=active.begin(),limit2=active.end();iter2!=limit2;++iter2)
            sum+=(*iter2).dydx*(step-(*iter2).from);
         func.push_back(Coord(step,sum));
         // Remove all entries inside this interval
         while ((!active.empty())&&(active.front().to<=step)) {
            pastValues+=active.front().value;
            pop_heap(active.begin(),active.end());
            active.pop_back();
         }
         if (!(step<nextScan)) break;
      }

      // Terminate if done
      if (iter==limit)
         break;
      scanPos=nextScan;
      active.push_back(*iter);
      push_heap(active.begin(),active.end());
   }

   return SplineHistogramBuilder::build(func);
}
//---------------------------------------------------------------------------
void SplineHistogram::addInterval(double from,double to,double value)
   // Add a score interval
{
   unsigned fromSlot=findSlot(from),toSlot=findSlot(to);
   while (fromSlot<toSlot) {
      double frac=(to-x[fromSlot])/(to-from);
      double v=value*frac;
      y[fromSlot]+=v;
      value-=v;
      from=x[fromSlot];
      fromSlot++;
   }
   y[fromSlot]+=value;
}
//---------------------------------------------------------------------------
SplineHistogram SplineHistogram::combineJoin(const SplineHistogram& other,double resultSize) const
   // Combine two histograms (join)
{
   // Caclculate some bounds
   SplineHistogram result;
   for (unsigned index=0;index<slotCount;index++) {
      result.x[index]=x[index]+other.x[index];
      result.y[index]=0;
   }

   // Calculate relative frequencies
   double a[slotCount-1];
   double b[slotCount-1];
   double div=totalCount();
   if (!div) div=1;
   for (unsigned index=0;index<slotCount-1;index++)
      a[index]=(y[index+1]-y[index])/div;
   a[0]=y[1]/div;
   div=other.totalCount();
   if (!div) div=1;
   for (unsigned index=0;index<slotCount-1;index++)
      b[index]=(other.y[index+1]-other.y[index])/div;
   b[0]=other.y[1]/div;

   // Estimate the matching tuples
   double intersectionPart,nonMatchingA,nonMatchingB;
   if (resultSize<0) {
      if (totalCount()>other.totalCount()) {
	 resultSize=totalCount();
	 intersectionPart=other.totalCount()/resultSize;
	 nonMatchingA=1-intersectionPart;
	 nonMatchingB=0;
      } else if (totalCount()<other.totalCount()) {
	 resultSize=other.totalCount();
	 intersectionPart=totalCount()/resultSize;
	 nonMatchingA=0;
	 nonMatchingB=1-intersectionPart;
      } else {
         resultSize=totalCount();
         intersectionPart=1;
         nonMatchingA=0;
         nonMatchingB=0;
      }
   } else {
      double intersection=totalCount()+other.totalCount()-resultSize;
      if (intersection<0) intersection=0;
      intersectionPart=intersection/resultSize;
      double div=totalCount()+other.totalCount(); if (!div) div=1;
      nonMatchingA=((resultSize-intersection)/resultSize)*totalCount()/div;
      nonMatchingB=((resultSize-intersection)/resultSize)*other.totalCount()/div;
   }

   // Construct the intervals
   for (unsigned index=0;index<slotCount-1;index++)
      for (unsigned index2=0;index2<slotCount-1;index2++) {
	 double from=x[index]+other.x[index2];
	 double to=x[index+1]+other.x[index2+1];
         result.addInterval(from,to,a[index]*b[index2]*intersectionPart*resultSize);
      }
   for (unsigned index=0;index<slotCount-1;index++) {
      result.addInterval(x[index],x[index+1],a[index]*nonMatchingA*resultSize);
      result.addInterval(other.x[index],other.x[index+1],b[index]*nonMatchingB*resultSize);
   }

   // Aggregate
   for (unsigned index=1;index<slotCount;index++)
      result.y[index]+=result.y[index-1];

   return result;
}
//---------------------------------------------------------------------------
SplineHistogram SplineHistogram::cutBelow(double score) const
   // Remove all entries below a certain threshold
{
   double below=totalCount()-documentsAboveOrEqual(score);

   SplineHistogram result;
   for (unsigned index=0;index<slotCount;index++) {
      result.x[index]=x[index];
      if (y[index]<=below)
         result.y[index]=0; else
         result.y[index]=y[index]-below;
   }

   return result;
}
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
SplineHistogram SplineHistogramBuilder::build(const std::vector<std::pair<double,double> >& func)
   // Build it
{
   // Construct a taut string
   std::vector<Coord> string;
   SplineHistogram result;
   unsigned bestError=~static_cast<unsigned>(0),left=1,right=100,expanding=10;
   double maxValue=func.size()?func[func.size()-1].second:0;
   while (left!=right) {
      unsigned middle=expanding?expanding:((left+right)/2);
      std::vector<Coord> string=tautString(func,maxValue,middle/1000.0);
      if (string.size()<=SplineHistogram::slotCount) {
         if (middle<bestError) {
            if (!string.size()) string.push_back(Coord(0,0));
            while (string.size()<SplineHistogram::slotCount)
               string.push_back(string.back());
            for (unsigned index=0;index<SplineHistogram::slotCount;index++) {
               result.x[index]=string[index].first;
               result.y[index]=string[index].second;
            }
            bestError=middle;
         }
         if (expanding) expanding=0;
         right=middle;
      } else if (expanding) {
         left=expanding+1;
         expanding*=2;
         if (expanding>right) expanding=0;
      } else {
         left=middle+1;
      }
   }
   // Try one more time
   if (left<bestError) {
      std::vector<Coord> string=tautString(func,maxValue,left/1000.0);
      if (string.size()<=SplineHistogram::slotCount) {
         if (!string.size()) string.push_back(Coord(0,0));
         while (string.size()<SplineHistogram::slotCount)
            string.push_back(string.back());
         for (unsigned index=0;index<SplineHistogram::slotCount;index++) {
            result.x[index]=string[index].first;
            result.y[index]=string[index].second;
         }
      }
   }

   return result;
}
//---------------------------------------------------------------------------
static double lineError(const std::vector<Coord>& data,unsigned from,unsigned to)
   // Compute the error by induced by a straight line
{
   // Check bounds
   if (from>to)
      std::swap(from,to);
   if ((from==to)||((from+1)==to))
      return 0;

   // Determine the function parameters
   double a=(data[to].second-data[from].second)/(data[to].first-data[from].first);
   double b=data[from].second-a*data[from].first;
   double maxValue=data.back().second;

   // And check all points in between
   double maxError=0;
   for (std::vector<Coord>::const_iterator iter=data.begin()+(from+1),limit=data.begin()+to;iter!=limit;++iter) {
      double y=a*((*iter).first)+b;
      double error=(*iter).second-y;
      if (error<0)
         error=-error;

      if (error>1.0)
         error-=1.0; else
         error=0;
      error/=maxValue-(*iter).second;

      if (error>maxError)
         maxError=error;
   }

   return maxError;
}
//---------------------------------------------------------------------------
namespace {
struct DPEntry {
   struct Step {
      double maxError;
      unsigned previous;
   };
   Step steps[SplineHistogram::maxSlotCount];
};
//---------------------------------------------------------------------------
static std::vector<Coord> compressFunc(const std::vector<Coord>& func,unsigned desiredSize)
   // Compress to the desired size
{
   // Relax a bit to speed up compression
   unsigned maxSize=desiredSize+(desiredSize/100),minSize=desiredSize-(desiredSize/100);

   // Fits?
   if (func.size()<=maxSize)
      return func;

   cerr << "max value = " << func.back().second << endl;
   
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
}
//---------------------------------------------------------------------------
#if 0
static double estimate(const std::vector<std::pair<double,double> >& func,double score)
{
   if (score<=func.front().first) return func.front().second;
   if (score>=func.back().second) return func.back().second;

   unsigned left=0,right=func.size();
   while (left<right) {
      unsigned middle=(left+right)/2;
      if (score>func[middle].first) left=middle+1; else
      if ((middle>0)&&(score<=func[middle-1].first)) right=middle; else
         { left=middle; break; }
   }

   if (left>=func.size()) return func.back().second;
   if ((left)&&(score<=func[left-1].first)) --left;

   if (score==func[left].first)
      return func[left].second;

   double dx=func[left].first-func[left-1].first;
   double dy=func[left].second-func[left-1].second;
   double r=func[left-1].second+((score-func[left-1].first)/dx)*dy;

   return r;
}
#endif
//---------------------------------------------------------------------------
SplineHistogram SplineHistogramBuilder::buildDP(const std::vector<std::pair<double,double> >& originalFunc)
   // Build it
{
   if (!originalFunc.size())
      return SplineHistogram();

   // Compress the input
   std::vector<Coord> func=compressFunc(originalFunc,5000);

#if 0
double maxRel=0;
   for (std::vector<Coord>::const_iterator iter=originalFunc.begin(),limit=originalFunc.end();iter!=limit;++iter) {
      double real = func.back().second - (*iter).second;
      double estimated = func.back().second - estimate(func,(*iter).first);

      double absError=(real<estimated)?(estimated-real):(real-estimated);
      if (absError>1) absError=absError-1; else absError=0;
      double relError=absError/real;
      if (relError>maxRel) {
         maxRel=relError;
      }
   }
std::cout << maxRel << std::endl;
#endif

   // Prepare the DP table
   std::vector<DPEntry> table;
   table.resize(func.size());
   for (unsigned index=0;index<SplineHistogram::slotCount;index++) {
      table[0].steps[index].maxError=0;
      table[0].steps[index].previous=0;
   }
   // Fill the remaining entries
   for (unsigned index=1;index<table.size();index++) {
      // Consider previous points
      for (unsigned index2=0;index2<index;index2++) {
         // Calculate the line error
         double error=lineError(func,index2,index);
         // First entry?
         if (index2==0) {
            for (unsigned index3=0;index3<SplineHistogram::slotCount;index3++) {
               table[index].steps[index3].maxError=error;
               table[index].steps[index3].previous=index2;
            }
         } else {
            // Fill the DP slots
            for (unsigned index3=1;index3<SplineHistogram::slotCount;index3++) {
               double maxError=error;
               if (table[index2].steps[index3-1].maxError>maxError)
                  maxError=table[index2].steps[index3-1].maxError;
               if (maxError<table[index].steps[index3].maxError) {
                  table[index].steps[index3].maxError=maxError;
                  table[index].steps[index3].previous=index2;
               }
            }
         }
      }
   }

   // Reconstruct the best path
   std::vector<unsigned> steps;
   for (unsigned index=SplineHistogram::slotCount,pos=table.size()-1;index>0;pos=table[pos].steps[index?index-1:0].previous) {
      --index;
      steps.push_back(pos);
   }

   // And build the histogram
   SplineHistogram result;
   unsigned slot=0;
   for (std::vector<unsigned>::const_reverse_iterator iter=steps.rbegin(),limit=steps.rend();iter!=limit;++iter,++slot) {
      result.x[slot]=func[*iter].first;
      result.y[slot]=func[*iter].second;
   }

   return result;
}
//---------------------------------------------------------------------------
namespace {
//---------------------------------------------------------------------------
/// Maintains an error front
class ErrorFront
{
   private:
   /// A front entry
   struct FrontEntry {
      /// Dominates at
      double dominates;
      /// Start
      double start;
      /// Weight
      double weight;

      /// Constructor
      FrontEntry(double dominates,double start,double weight) : dominates(dominates),start(start),weight(weight) {}
   };
   /// The front
   std::vector<FrontEntry> front;
   /// The last operations
   unsigned lastInsert,lastLookup;

   /// Prune entries in the front
   void pruneFront(unsigned index);
   /// Prune entries in the front
   void pruneFrontBackwards(unsigned index);

   public:
   /// Constructor
   ErrorFront();
   /// Destructor
   ~ErrorFront();

   /// Insert an error function
   void insert(double a,double weight);
   /// Lookup an error entry
   double error(double a);
};
//---------------------------------------------------------------------------
ErrorFront::ErrorFront()
   : lastInsert(0),lastLookup(0)
   // Constructor
{
}
//---------------------------------------------------------------------------
ErrorFront::~ErrorFront()
   // Destructor
{
}
//---------------------------------------------------------------------------
static inline double breakEven(double s1,double w1,double s2,double w2)
   // The break even point
{
   return -(s1*w1-s2*w2)/(w2-w1);
}
//---------------------------------------------------------------------------
void ErrorFront::pruneFront(unsigned index)
   // Prune entries in the front
{
   while (index+1<front.size()) {
      double even=breakEven(front[index].start,front[index].weight,front[index+1].start,front[index+1].weight);
      if (even<=front[index].dominates) {
         front.erase(front.begin()+(index+1));
      } else break;
   }
}
//---------------------------------------------------------------------------
void ErrorFront::pruneFrontBackwards(unsigned index)
   // Prune entries in the front
{
   while (index) {
      double even=breakEven(front[index].start,front[index].weight,front[index-1].start,front[index-1].weight);
      if (even<=front[index-1].dominates) {
         front.erase(front.begin()+(--index));
      } else break;
   }
   pruneFront(index);
}
//---------------------------------------------------------------------------
void ErrorFront::insert(double a,double weight)
   // Insert an error functio
{
   // First entry?
   if (!front.size()) {
      front.push_back(FrontEntry(a,a,weight));
      return;
   }
   // No, try to find the first possible position
   if (lastInsert>=front.size())
      lastInsert=front.size()-1;
   while (lastInsert&&(a<front[lastInsert].start))
      lastInsert--;
   if ((!lastInsert)&&(a<front[lastInsert].start)) {
      front.insert(front.begin(),FrontEntry(a,a,weight));
      pruneFront(lastInsert);
      return;
   }

   // Identical start?
   if (a==front[lastInsert].start) {
      front[lastInsert].weight=weight;
      pruneFrontBackwards(lastInsert);
      return;
   }

   // Check for dominance
   double even = breakEven(front[lastInsert].start,front[lastInsert].weight,a,weight);
   if (even<=front[lastInsert].dominates) {
      while (true) {
         front.erase(front.begin()+lastInsert);
         if (!lastInsert) {
            even=a;
            break;
         }
         even=breakEven(front[lastInsert-1].start,front[lastInsert-1].weight,a,weight);
         if (even<=front[lastInsert-1].dominates)
            --lastInsert; else
            break;
      }
      if (even<a) even=a;
      front.insert(front.begin()+lastInsert,FrontEntry(even,a,weight));
      pruneFront(lastInsert);
      return;
   }
   while (true) {
      // Last entry? Then append
      if (lastInsert==front.size()-1) {
         if (even<a) even=a;
         front.push_back(FrontEntry(even,a,weight));
         ++lastInsert;
         return;
      }
      // No, check with the next entry
      if (even<=front[lastInsert+1].dominates) {
         if (even<a) even=a;
         front.insert(front.begin()+(++lastInsert),FrontEntry(even,a,weight));
         pruneFront(lastInsert);
         return;
      }
      // Move forward
      ++lastInsert;
      even = breakEven(front[lastInsert].start,front[lastInsert].weight,a,weight);
   }
}
//---------------------------------------------------------------------------
double ErrorFront::error(double a)
   // Lookup an error entry
{
   // Empty front?
   if (!front.size())
      return 0;
   if (lastLookup>=front.size())
      lastLookup=front.size();

   // Find the first possible pos
   while (lastLookup&&(a<front[lastLookup].dominates))
      --lastLookup;
   if (a<front[lastLookup].dominates)
      return 0;

   // Find the correct pos
   while ((lastLookup+1<front.size())&&(a>=front[lastLookup+1].dominates))
      lastLookup++;
   return (a-front[lastLookup].start)*front[lastLookup].weight;
}
//---------------------------------------------------------------------------
}
//---------------------------------------------------------------------------
SplineHistogram SplineHistogramBuilder::buildDPQuad(const std::vector<std::pair<double,double> >& originalFunc)
   // Build it
{
   if (!originalFunc.size())
      return SplineHistogram();

#if 0
{
   ErrorFront upperFront,lowerFront;
   Coord base=originalFunc[0];
   double maxValue=originalFunc.back().second;
   for (unsigned to=1;to<originalFunc.size();to++) {
      Coord current=originalFunc[to];
      current.first-=base.first;
      current.second-=base.second;
      double a=current.second/current.first;

      // Determine the error
      double realError=lineError(originalFunc,0,to);
      double estimatedError=std::max(upperFront.error(a),lowerFront.error(-a));
      std::cout << to << " " << realError << " " << estimatedError << std::endl;

      // Update the front
      double ua=(current.second+1)/current.first,la=-((current.second-1)/current.first);
      double weight=current.first/(maxValue-originalFunc[to].second);

      upperFront.insert(ua,weight);
      lowerFront.insert(la,weight);
   }
}
#endif

#if 1
   // Compress the input
   std::vector<Coord> func=compressFunc(originalFunc,5000);
#else
   const std::vector<Coord>& func = originalFunc;
#endif

   // Prepare the DP table
   std::vector<DPEntry> table;
   table.resize(func.size());
   for (unsigned index=0;index<SplineHistogram::slotCount;index++) {
      table[0].steps[index].maxError=0;
      table[0].steps[index].previous=0;
   }
   // Propagate the info
   double maxValue=func.back().second;
   for (unsigned from=0;from<table.size();from++) {
      ErrorFront upperFront,lowerFront;
      Coord base=originalFunc[from];

      // Consider all possible targets
      for (unsigned to=from+1;to<table.size();to++) {
         // Calculate the error
         Coord current=func[to];
         current.first-=base.first;
         current.second-=base.second;
         double a=current.second/current.first;
         double error=std::max(upperFront.error(a),lowerFront.error(-a));

         // First entry?
         if (from==0) {
            for (unsigned index=0;index<SplineHistogram::slotCount;index++) {
               table[to].steps[index].maxError=error;
               table[to].steps[index].previous=from;
            }
         } else {
            // No, scan all slots
            for (unsigned index=1;index<SplineHistogram::slotCount;index++) {
               double maxError=std::max(error,table[from].steps[index-1].maxError);
               if (maxError<table[to].steps[index].maxError) {
                  table[to].steps[index].maxError=maxError;
                  table[to].steps[index].previous=from;
               }
            }
         }

         // Update the error front
         double ua=(current.second+1)/current.first,la=-((current.second-1)/current.first);
         double weight=current.first/(maxValue-func[to].second);

         upperFront.insert(ua,weight);
         lowerFront.insert(la,weight);
      }
   }

   // Reconstruct the best path
   std::vector<unsigned> steps;
   for (unsigned index=SplineHistogram::slotCount,pos=table.size()-1;index>0;pos=table[pos].steps[index].previous) {
      --index;
      steps.push_back(pos);
   }

   // And build the histogram
   SplineHistogram result;
   unsigned slot=0;
   for (std::vector<unsigned>::const_reverse_iterator iter=steps.rbegin(),limit=steps.rend();iter!=limit;++iter,++slot) {
      result.x[slot]=func[*iter].first;
      result.y[slot]=func[*iter].second;
   }

   return result;
}
//---------------------------------------------------------------------------
using namespace std;
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
int main(int argc,char* argv[])
{
   if (argc<4) {
      cerr << "usage: " << argv[0] << " file samplingPoints countToRead" << endl;
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
    auto ts = compressFunc(cdf, size);
    //auto ts = tautString(cdf, cdf.back().second, 0.60001);
#endif
      
      double sumError=0, maxError=0, relError = 0;
      for (auto e:cdf) {
         double estimate=interpolate(ts,e.first);
         double real=e.second;
         double diff=estimate-real;
         if (diff<0) diff=-diff;
         if (diff>maxError)
            maxError=diff;
         sumError+=diff;
         relError += (!real) ? diff : (diff / real);
      }
      cerr << ts.size() << " " << maxError << " " << (sumError/cdf.size()) << " rel :" << (relError / cdf.size()) << endl;
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

