#ifndef _HashLink_h
#define _HashLink_h

#include "HashResult.h"
#include "macros.h"

#include <math.h>
#include <cstdlib>
#include <iostream>
#include <cstddef>

struct GeomHashParams {
  GeomHashParams(int dim, std::vector<float> cube_sizes, float resize_coefficient)
    : dimension(dim), cubeSizes(cube_sizes), resizeCoefficient(resize_coefficient)
  {
    if (cubeSizes.size() != (unsigned int) dimension) {
      std::cerr << "Invalid Parameter: the size of the cube_sizes vector has to be equal to the dimension of the hash" << std::endl;
      exit(1);
    }
  }

  GeomHashParams(int dim, float cube_size, float resize_coefficient)
    : dimension(dim), cubeSizes(dim, cube_size), resizeCoefficient(resize_coefficient)
  {}

  int dimension;
  std::vector<float> cubeSizes;
  float resizeCoefficient;
};

template<class KeyT, class DataT, class Bucket>
/*
CLASS
  HashLink

  Template serves as the main engine of the geometric hash table. Classes
  defined are used by GeomHash to build, traverse and query a multi-dimensional
  hash-table.

KEYWORDS
  geometric, hash, iterator, query, container

AUTHORS
  Meir Fuchs (meirfux@math.tau.ac.il).
  Zipi Fligelman (zipo@math.tau.ac.il).

CHANGES LOG
<UL>
  <LI>5/12/13 - change short to int
  </LI>
  <LI>27/02/02 - Dina Duhovny & Oranit Shem-Tov:
    - GeomHashParams - We refer to the dimension attribute as the dimension
      of the overall geometric hash. Previously, this attribute was reffered
      to as the index of the last axis (i.e. dimension-1).
      Axis explanation by an example:
      Assume a 2D geometric hash. It's axes are numbered as 0 (the x-axis) & 1 (the y-axis).
    - Defining cube size for each axis. Previously, there was only one value
      for all the axis. It is usefull for storing data and a key, composed of
      several components, each one of them has a different scale. For example,
      a 4D key, composed of a 3D coordinates and a type.
    - Defining query ranges - for querying, within a different radius for axis.
      the previous get function turned to getInfinite. New get function returns
      the results of query by cutting each axis with different radii
  </LI>
  <LI>8/8/99 - Ram Nathaniel:<BR>
    Add support for GeomHash bucket iteration - method
    getAll(int dim, vector<Bucket *>* result) const;
   </LI>
   <LI>7/11/99- Zipi Fligelman <br>
   getAll added by Ram does not exist in that format anymore.
   There are 4 new functions:
   1. getAllNotEmpty - that recieves either a hash result or a pointer to a vector of
      pointers to buckets, this will return all the non empty backets in the ALL
      the table
   2. getAll - that really gets it all including the empty buckets in the last level
      it can recieve the same parameters as getAllNotEmpty, but the user must take
      extra care when processing it later not to fall into a segmentation fault.
   3. There is a defintion of vector<const Bucket*> as BucketsPointerList.
      </LI>
  <LI> 15/5/03 - M. Shatsky
  Copy constructor is added
  </LI>
</UL>

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

GOALS
  HashLink is the main engine of the geometric hash table. It is responisble
  for construction, traversal and querying of a multi-dimensional hash table.
  The class is wrapped with the GeomHash class which packages the complex
  methods and intricate implementation of HashLink with simple methods
  that can be used immediately. No methods of the class are public because
  this class was not meant to be used directly.

  A template instance is declared using a key class KeyT which defines the
  means by which the hash table is referenced and queried and a data class
  DataT which defines the data the hash table holds. It is assumed that DataT
  has an operator[] which will return floating point coordinates for each
  dimension. An optional template parameter, allows the user to define
  different bucket types. The default bucket type is a vector<DataT>.

  The HashLink holds a dynamic array of HashLinks. This array is dynamic in
  the sense that it may grow but it is also dynamic in the sense that it
  has a virtual starting point. It is a "ranged" array. Although
  there may be a slight performance cost during data insertion, this
  implementation saves memory and rids the user of the need to pre-define
  sizes and ranges.

USAGE
  Since the class is not meant for direct use usage will not be described.
  See the GeomHash implementation for usage examples. The GeomHash methods
  pretty much coincide with the HashLink methods.
*/
class HashLink
{
public:
  typedef std::vector<const Bucket* > BucketsPointerList;

private:
  HashLink(const HashLink<KeyT, DataT, Bucket>&){
    std::cerr<<"HashLinkCC" << std::endl; exit(0);}
protected:
  //// initialized as NULL. empty HashLink.
  HashLink();

  ////
  HashLink(const HashLink<KeyT, DataT, Bucket>& hl,const int axis);

  //// Returns the cube index of coord using cube size cs. The
  // cube index is returned without normalizing to range.
  static int cube(const float coord, const float cs);

  //// resize the array of HashLinks using a new range defined by newStart
  // and newStop. It is assumed that the new range is bigger than the old one.
  void resize(const int newStart, const int newStop);

  //// Resize for the bottom axis. This resize is used for the last
  // axis which only contains pointers to buckets.
  void buckResize(const int newStart, const int newStop);

  //// put will insert a data element into the hash table. put is recursive
  // it calls itself until at axis=0. it adds the data element to a bucket
  // object it either finds or creates. The bucket objects are by default of
  // vector<DataT> type.
  //
  // put must make sure that it can continue to recurse until axis = 0.
  // this means that the array of HashLinks must be not null, and the the
  // cube index must fall in the array's range. If this is not the case put
  // must resize the array. Arrays are allocated with extra elements at
  // their head and tail in an order to reduce the number of resize
  // operations.
  void put(const KeyT& key, const DataT& data, const GeomHashParams& params,
           int axis);

  //// getInfinite returns a list of elements that may be within a given radius to a
  // given key. The elements returned are all the elements found in the
  // neighboring cubes to the given key's cube. The data elements from all
  // these cubes are gathered in a HashResult container class passed on
  // to the method as an argument. get is a recursive method that recurses
  // until axis=0 at which point a bucket pointer is added to the HashResult.
  void getInfinite(const KeyT &key, const float radius, const GeomHashParams& params, int axis,
		   HashResult<DataT, Bucket>& result) const;

  //// getEuclidean returns a list of elements according to Euclidean norm
  void getEuclidean(const KeyT &key, const float radius, const GeomHashParams& params,
                    int axis,
                    HashResult<DataT, Bucket>& result) const;

  //// getManhattan returns a list of elements according to Manhattan norm
  void getManhattan(const KeyT &key, const float radius, const GeomHashParams& params,
                    int axis,
                    HashResult<DataT, Bucket>& result) const;


  //// get returns a list of elements that may be within a given set of radii to a
  // given key. The elements returned are all the elements found in the
  // neighboring cubes to the given key's cube. The data elements from all
  // these cubes are gathered in a HashResult container class passed on
  // to the method as an argument. get is a recursive method that recurses
  // until axis=0 at which point a bucket pointer is added to the HashResult.
  //
  // params:
  // ranges - a pointer to an array of size equal to the axis of the hash.
  // The user is the responsible for a valid array.
  // axis - is an index of the current checked axis.
  void get(const KeyT &key, const std::vector<float>& radii,
	   const GeomHashParams& params, int axis,
	   HashResult<DataT, Bucket>& result) const;

  //// This version creates a vector of all the non-empty  buckets.
  // The buckets are never changed during the use of the GeomHash, and therefore
  // this function create a snapshot of total buckets currently in use in the hash.
  // The function is also usefull for collecting information about the buckets
  // which may be traversed in a Bucket iterator now.
  void getAllNotEmpty(int axis, BucketsPointerList *result) const;

  //// This funciton will return all the buckets including the empty (null) ones
  // that exists in the hash table. This function enable to check the saprsity
  // of the hash table at least in the last level of the recursion
  void getAll(int axis, BucketsPointerList *result) const;

  //// Using get with radius = 0 makes things simpler in terms of code and
  // method interface. This method returns a bucket pointer to the bucket
  // found in the position pointed to by the given key. If no bucket is found
  // a NULL pointer is returned. This method should be faster then the general
  // get method.
  const Bucket* getBucket(const KeyT &key, const GeomHashParams& params, int axis) const;

  //// erase recursively erases all the data elements in the array.
  // This method is not called by the destructor since the HashLink elements
  // are destructed during the insertion operations in the put function.
  /// need to add a possibilty win axis 0 when the data is a pointer
  void erase(const int axis);

private:
  HashLink* next;
  int start, stop;
};

/*******************************************
 * Template class HashLink implementation. *
 *******************************************/

template<class KeyT, class DataT, class Bucket>
HashLink<KeyT, DataT, Bucket>::HashLink(): next(NULL) {}

template<class KeyT, class DataT, class Bucket>
HashLink<KeyT, DataT, Bucket>::HashLink(const HashLink<KeyT, DataT, Bucket> &hl,const int axis)
{
  if(hl.next == NULL){
    next=NULL;
    return;
  }

  start=hl.start;
  stop=hl.stop;
  unsigned int size = stop - start + 1;

  if (axis) {
    // allocate memory for new array.
    next = (HashLink*)::operator new(size*sizeof(HashLink));

    // copy elements
    HashLink* copyEnd = next+(stop -start);
    for (HashLink* p=next, *q=hl.next; p<=copyEnd; ++p, ++q)
      *p = HashLink(*q,axis-1);

  }else{
    Bucket** newnext = (Bucket**)::operator new(size*sizeof(Bucket*));

    // copy elements and instantiate empty ones.
    Bucket** copyEnd = newnext + (stop -start);
    for (Bucket** p=newnext, **q=(Bucket**)hl.next; p<=copyEnd; ++p, ++q){
      if(*q)
      	*p = new Bucket(**q);
      else
        *p=NULL;
    }

    next=(HashLink*)newnext;
  }
}

template<class KeyT, class DataT, class Bucket>
int HashLink<KeyT, DataT, Bucket>::cube(const float coord,
                                        const float cs)
{
  // There is a problem with truncation of negative numbers. Negative numbers
  // are truncated and one recieves the number rounded up while positive
  // numbers are rounded down.
  // There is a slight unimportant bug here. But I've ignored it for the sake
  // of speed. Assuming cs=1.0 then coord=0.0 will return 0 while coord=-1.0
  // will return -2. Cube -1 is "infinitesimally smaller". Since all HashLink
  // methods use the cube method this should have no bearings on the actual
  // output of the hash table.
  if (coord < 0)  return (int)(coord/cs) - 1;
  else            return (int)(coord/cs);
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::resize(const int newStart,
                                           const int newStop)
{
  // allocate memory for new array.
  unsigned int newSize = newStop - newStart + 1;
  HashLink* newNext = (HashLink*)::operator new(newSize*sizeof(HashLink));

  // copy elements and instantiate empty ones.
  HashLink* copyStart = newNext + (start-newStart);
  HashLink* copyEnd = copyStart + (stop -start + 1);
  HashLink* newEnd = newNext + newSize;
  for (HashLink* p=newNext; p<copyStart; ++p)               *p = HashLink();
  for (HashLink* p=copyStart, *q=next; p<copyEnd; ++p, ++q) *p = *q;
  for (HashLink* p=copyEnd; p<newEnd; ++p)                  *p = HashLink();

  // free space of old array
  ::operator delete(next);

  // update the HashLink parameters
  next = newNext;
  start = newStart;
  stop = newStop;
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::buckResize(const int newStart,
                                               const int newStop)
{
  unsigned int newSize = newStop - newStart + 1;
  Bucket** newNext = (Bucket**)::operator new(newSize*sizeof(Bucket*));

  // copy elements and instantiate empty ones.
  Bucket** copyStart = newNext + (start-newStart);
  Bucket** copyEnd = copyStart + (stop -start + 1);
  Bucket** newEnd = newNext + newSize;
  for (Bucket** p=newNext; p<copyStart; ++p)
    *p = NULL;
  for (Bucket** p=copyStart, **q=(Bucket**)next; p<copyEnd; ++p, ++q)
    *p = *q;
  for (Bucket** p=copyEnd; p<newEnd; ++p)
    *p = NULL;

  // free space of old array
  ::operator delete(next);

  // update the HashLink parameters
  next = (HashLink*)newNext;
  start = newStart;
  stop = newStop;
}


template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::put(const KeyT& key, const DataT& data,
                                        const GeomHashParams& params,
                                        const int axis)
{
  int pos = cube(key[axis], params.cubeSizes[axis]);

  if (axis) {
    // axis is not 0. have to go deeper.
    // check that new element fits into array. if not resize with a
    // few extra elements (to prevent constant resizing). number of
    // extra elements depends on axis.

    if (next==NULL) {    // if array not started create new array
      start = pos; stop = pos-1;
      resize(pos-axis, pos+axis);
    }
    else if (pos<start)
      resize(pos-(int)(params.resizeCoefficient*axis), stop);
    else if (pos>stop)
      resize(start, pos+(int)(params.resizeCoefficient*axis));
    // activate the put method on the correct element within the array.
    next[pos-start].put(key, data, params, axis-1);
  } else {
    // axis=0. time to insert the data element into a bucket. at axis=0 next
    // points to a bucket. make sure one exists and add the data to the
    // bucket.
    if (next==NULL) {    // if array not started create new array
      start = pos; stop = pos-1;
      buckResize(pos-1, pos+1);
    }
    else if (pos<start)
      buckResize(pos-(int)(params.resizeCoefficient+1), stop);
    else if (pos>stop)
      buckResize(start, pos+(int)(params.resizeCoefficient+1));
    // activate the put method on the correct element within the array.
    Bucket** bucket = ((Bucket**)next)+(pos-start);
    if (*bucket==NULL)   *bucket = new Bucket;
    (*bucket)->push_back(data);
  }
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::getInfinite(const KeyT &key, const float radius,
						const GeomHashParams& params, int axis,
						HashResult<DataT, Bucket>& result) const
{
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int i = cube(key[axis]-radius, params.cubeSizes[axis]);
    i = ((start > i) ? start : i) - start;
    int to = cube(key[axis]+radius, params.cubeSizes[axis]);
    to = ((stop < to) ? stop : to) - start;
    if (axis)
      for (; i<=to; ++i)
        next[i].getInfinite(key, radius, params, axis-1, result);
    else
      for (Bucket** buckets = (Bucket**)next; i<=to; ++i)
        if (buckets[i] != NULL)
          result.push_back(buckets[i]);
  }
  return;
}


template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::get(const KeyT &key, const std::vector<float>& radii,
                                        const GeomHashParams& params, int axis,
                                        HashResult<DataT, Bucket>& result) const
{
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int i = cube(key[axis]-radii[axis], params.cubeSizes[axis]);
    i = ((start > i) ? start : i) - start;
    int to = cube(key[axis] + radii[axis], params.cubeSizes[axis]);
    to = ((stop < to) ? stop : to) - start;
    if (axis)
      for (; i<=to; ++i)
        next[i].get(key, radii, params, axis-1, result);
    else
      for (Bucket** buckets = (Bucket**)next; i<=to; ++i)
        if (buckets[i] != NULL)
          result.push_back(buckets[i]);
  }
  return;
}


template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::getManhattan(const KeyT &key, const float radius,
                                            const GeomHashParams& params, int axis,
                                            HashResult<DataT, Bucket>& result) const
{
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int i = cube(key[axis]-radius, params.cubeSizes[axis]);
    i = ((start > i) ? start : i);
    int to = cube(key[axis]+radius, params.cubeSizes[axis]);
    to = ((stop < to) ? stop : to);
    if (axis) {
      int center = cube(key[axis], params.cubeSizes[axis]);
      for (; i<=to; ++i) {
        if (i==center)
          next[i-start].getManhattan(key, radius, params, axis-1, result);
        else {
          float rdiff;
          if (i<center)
            rdiff = key[axis]-(i+1) * params.cubeSizes[axis];
          else
            rdiff = i * params.cubeSizes[axis] - key[axis];
          if (rdiff <= radius)
            next[i-start].getManhattan(key, radius-rdiff,
                                       params, axis-1, result);
        }
      }
    }
    else {
      Bucket* buck;
      for (Bucket** buckets = (Bucket**)next; i<=to; ++i)
        if ((buck = buckets[i-start]) != NULL)
          result.push_back(buck);
    }
    return;
  }
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::getEuclidean(const KeyT &key, const float radius,
                                            const GeomHashParams& params, int axis,
                                            HashResult<DataT, Bucket>& result) const
{
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int i = cube(key[axis]-radius, params.cubeSizes[axis]);
    i = ((start > i) ? start : i);
    int to = cube(key[axis]+radius, params.cubeSizes[axis]);
    to = ((stop < to) ? stop : to);
    if (axis) {
      int center = cube(key[axis], params.cubeSizes[axis]);
      for (; i<=to; ++i) {
        if (i==center)
          next[i-start].getEuclidean(key, radius, params, axis-1, result);
        else {
          float rdiff;
          if (i<center)
            rdiff = key[axis]-(i+1) * params.cubeSizes[axis];
          else
            rdiff = i* params.cubeSizes[axis] - key[axis];
          if (rdiff <= radius)
            next[i-start].getEuclidean(key, sqrt(sqr(radius)-sqr(rdiff)),
                                       params, axis-1, result);
        }
      }
    }
    else {
      Bucket* buck;
      for (Bucket** buckets = (Bucket**)next; i<=to; ++i)
        if ((buck = buckets[i-start]) != NULL)
          result.push_back(buck);
    }
    return;
  }
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::getAllNotEmpty(int axis, BucketsPointerList* result) const
{
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int i = 0;
    int to = stop - start;
    if (axis)
      for (; i<=to; ++i)
	next[i].getAllNotEmpty(axis-1, result);
    else
      for (Bucket** buckets = (Bucket**)next; i<=to; ++i)
	if (buckets[i] != NULL)
	  result->push_back(buckets[i]);
  }
  return;
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::getAll(int axis, BucketsPointerList* result) const
{
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int i = 0;
    int to = stop - start;
    if (axis)
      for (; i<=to; ++i)
	next[i].getAll(axis-1, result);
    else
      for (Bucket** buckets = (Bucket**)next; i<=to; ++i)
	result->push_back(buckets[i]);
  }
  return;
}

template<class KeyT, class DataT, class Bucket>
const Bucket* HashLink<KeyT, DataT, Bucket>::getBucket(const KeyT &key, const GeomHashParams& params,
						       int axis) const {
  if (next) {
    // go over all cubes in the array that may contain members that are in
    // radius proximity to key and recursively call get.
    int pos = cube(key[axis], params.cubeSizes[axis]);
    if (pos >= start && pos <= stop) {
      if (axis)
        return next[pos-start].getBucket(key, params, axis-1);
      else
        return *(((Bucket**)next)+pos-start);
    }
  }
  return NULL;
}

template<class KeyT, class DataT, class Bucket>
void HashLink<KeyT, DataT, Bucket>::erase(const int axis)
{
  // recursively call for each element in the array and then free the array's
  // memory.
  if (next){
    if (axis) {
      HashLink* from = next;
      HashLink* to = next+(stop-start);
      for (; from<=to; ++from)
        from->erase(axis-1);
    }
    else {
      Bucket** from = (Bucket**)next;
      Bucket** to = from+(stop-start);
      for (; from<=to; ++from)
        if (*from != NULL) delete(*from);
    }
    ::operator delete(next);

    next=NULL;
  }
}

#endif
