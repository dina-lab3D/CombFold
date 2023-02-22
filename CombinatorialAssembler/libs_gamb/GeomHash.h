#ifndef _GeomHash_h
#define _GeomHash_h

#include "HashLink.h"

template<class KeyT, class DataT, class BucketT = std::vector<DataT> >
/*
CLASS
  GeomHash

  The template defines a geometric hash table class. The class is declared
  using a key type with which the data type instances are inserted and
  accessed. The geometric hash table may be instantiated using any dimension
  and grid-size.

KEYWORDS
  geometric, hash, iterator, query, container

AUTHORS
  Meir Fuchs (meirfux@math.tau.ac.il).
  Zipi Fligelman (zipo@math.tau.ac.il).

CHANGES LOG
<UL>
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
  - Defining query ranges - for querying, within a different radius for each axis.
  </LI>
  <LI>8/8/99 - Ram Nathaniel:<BR>
  Add bucket iteration - method vector<Bucket *> *getBuckets() const;
  </LI>
  <LI> 7/11/99 - Zipi Fligelman
  Changed the method by Ram and added parameter at the whole table query
  In the full table query there is a flag that by default is set to false, that
  means you collect only non empty buckets. But when it is set to true non-empty
  buckets are also collected. The same is true for the method getBuckets.
  That returns a vector of <b> const Bucket * </b> with all the buckets
  The empty backets are returned only if the flag is set to true.
  </LI>
  <LI> 15/5/03 - M. Shatsky
  Copy constructor is added
  </LI>
</UL>
  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

GOALS
  GeomHash defines a high performance geometric hash table with an intuitive
  interface. The class can be defined for any key type that supports the []
  operator for any data type. The hash table can be instantiated using any
  dimension and any cube or grid size. It is possible to specify different
  cube size for each axis. The query function also may be called with a set of
  radii, one for each axis. User may also decide how the data
  elements are collected in each bucket by determining the container used in
  each bucket.

  GeomHash actually is a wrap for the HashLink class which does most of
  the work of hash table traversal, insertion and querying. The user need
  not know anything about HashLink to use the hash table.

  Geometric hash table results are returned using an STL-like container class
  called HashResult. This class and its accompanying iterator class
  HashIterator are used to collect hash table query results and to iterate
  through the result elements.

USAGE
  Instantiation and declaration are done using a key type and a data type.
  The key type must support the [] operator. For example an array of floats
  or a Vector3 for 3 dimensional hash tables. User may also choose to determine
  the bucket type. By default the bucket type is a vector<DataT>. The user
  may decide on other STL containers or write one of his own. The bucket type
  must support the push_back, begin, end, size methods and the iterator and
  const_iterator typedefs.

  Instantiation may include a grid size and a dimension.
  Here is an example using Vector3's as the keys,
  integers as the data stored in the hash table, 3 is the dimension and 1.3 is
  the cube size of the grid.
  EXAMPLE
    GeomHash<Vector3, int> gHash(3, 1.3);
  END

  Members may be inserted using the insert command. The insert command
  recieves a key and a data element to be inserted into the proper bucket
  in the hash table.
  EXAMPLE
    gHash.insert(Vector3(1.2, 2.3, 3.4), -1);
  END

  When querying a geometric hash table one searches for data that was
  inserted into the hash table using keys that were close to the query
  key. A radius may be given to the query defining how far the neighboring
  points should be. The query returns all points within that distance but
  will return other points that lie in neighboring grid cubes. If you wish
  to extract points that are within a given euclidean radius you'll have
  to sift through the query results.

  Before running a query one must declare a HashResult container in which the
  results will be collected. The container is sent with the query request
  and results are found in the container after the query is done.

  An example of querying the hash table with a key and sifting through the
  data to find the data's that are within a euclidean radius of the key is
  given. The example assumes 10000 points are given inside a Vector3 array
  called points and a query point is given in Vector3 q.
  EXAMPLE
  float radius=1.0;
  GeomHash<Vector3, int> gHash(3, 2*radius);  // grid cube size is 2*radius
  for (int i=0; i<10000; i++)   // insert 10000 points
    gHash.insert(points[i], i);

  Vector3 q(1,1,1);
  HashResult<int> result;
  gHash.query(q, radius, result);  // radius = half of cube size.

  // traverse list STL iterator style.
  HashResult<int>::iterator it;
  for (it = result.begin(); it != result.end(); it++)
    if (q | points[*it]) <= radius
      cout << *it << " : " << points[*it] << '\n';
  END

  In case of more complex key, the user can specify the cube size of every axis.
  For example, in case of transformation: each transformation can be represented
  as 3 rotation parameters and 3 translation parameters. The rotation parameters
  are values between -PI and PI, while the translational parameters are not limited.
  The user may want smaller bin size for rotation then for translation.
  EXAMPLE
  vector<float> radii(6);
  radii[0]=radii[1]=radii[2]=0.1;
  radii[3]=radii[4]=radii[5]=2.0;
  GeomHash<Transformation,int> gHash(6, radii);

  The same idea works in query:
  HashResult<int> result;
  gHash.query(t, radii, result);
  END
*/

class GeomHash :
    private HashLink<KeyT, DataT, BucketT>
{
  public:
  typedef HashResult<DataT, BucketT> Result;
  typedef typename Result::iterator Iterator;
  typedef BucketT Bucket;
  typedef typename HashLink<KeyT,DataT, BucketT>::BucketsPointerList BucketsPointerList;

  typedef enum {
    manhattanNorm,
    euclideanNorm,
    infiniteNorm
  } NormType;

  //// The hash-table is instantiated by a given dimension (defualt 1)
  // and a given cube-size (default 1.0).
  explicit GeomHash(const unsigned short d=1, const float cs=1.0);

  //// The hash-table is instantiated by a given dimension and
  //// a given cube size for each dimension.
  explicit GeomHash(const unsigned short d, const std::vector<float>& cs);

  //// Destructor. Frees all hash table memory.
  ~GeomHash()  {
    erase();
  }

  //// Copy constructor
  GeomHash(const GeomHash<KeyT,DataT,BucketT>&);

  //// Data elements can be inserted according to keys using the insert method.
  // The key determines in which d-dimensional cube the data is inserted.
  void insert(const KeyT& key, const DataT& data);

  //// Searches for all data elements lying in radius neighborhood to key.
  // Note that query will return all points lying within a given radius
  // but may also return extra points that were not inserted with keys
  // that are within a radius to the query key.
  // Moreover, the default norm is the infinite norm
  // All results are collected using a results container of type HashResult<DataT>.
  void query(const KeyT &key, const float radius, Result& result) const;

  //// Searches for all data elements according to a given key.
  // The user may specify the type of the norm: Euclidean, Manhattan or Infinite.
  // Note that query will return all points lying within a given radius
  // but may also return extra points that were not inserted with keys
  // that are within a radius to the query key.
  // All results are collected using a results container of type HashResult<DataT>.
  void query(const KeyT &key, const float radius, const NormType ballType,
             Result& result) const;


  //// Searches for all data elements according to a given key and a set of radii.
  // the size of the radii vector must be equal to dimension. Each axis is searches
  // with a radius specified for it.
  // This query may return extra points as well.
  void query(const KeyT &key, const std::vector<float>& radii, Result& result) const;

  //// Searches for all data elements lying in radius=cube size neighborhood
  // to key. The query will return all data's that were inserted with keys
  // whose distance to the query key is less then radius. More data
  // elements may be returned due to the grid structure of the hash table.
  // All results are
  // collected using a results container of type HashResult<DataT>.
  void query(const KeyT &key, Result& result) const;

  //// Query the hash table with an infinite radius. Calling this method
  // results in a hash table results container which points to all hash table
  // elements. There is a choice whether to collect only non empty
  // buckets (ehich is also the default), but if one wished to check for
  // sparsity of the table you can set the includeEmpty to true and receive all
  // the empty buckets. Notice that travesing the result must be done with caution
  // mind the NULLs ....
  void query(Result& result, bool includeEmpty=false) const;

  //// Return a vector of all the buckets currently used in the Hash. This is
  // usefull for traversing the hash one bucket at a time, and also usefull for
  // creating a set of all the buckets now in use for follow up checks. There
  // is a flag that can be set if you want to add all the empty buckets in the
  // last dimesion and thus check the change in the sparsity of the table.(The
  // Hash class never changes the buckets but only copy a constant version of the
  // pointers so that even after insertions and deletions in the Hash the vector
  // will indicate the same set of buckets).
  BucketsPointerList *getBuckets(bool includeEmpty=false) const;

  //// Returns a pointer to a bucket at the location given by key. If no such
  // bucket exists a NULL pointer is returned.
  const Bucket* bucket(const KeyT& key) const;

  //// frees all memory used by hash table elements.
  void erase();

  //// Sets the hash table insertion speed. This parameter determines how
  // large the hash table arrays are made when an array resize is forced.
  // When this parameter is set to 0, arrays will be resized to be slightly
  // larger then what is currently needed by the hash table. If the parameter
  // is set to > 0, resizes will be more generous. Array sizes will grow
  // causing a decrease in the number of resize operations and hence an
  // increase in insertion speed.
  //
  // This feature should be useful for large but sparsely populated hash tables
  // In densely populated hash tables, a gain in speed should be noticed only
  // with the first insertions. Once the table is populated resize operations
  // are rare and hence the speed gain is limited.
  //
  // Note that raising the parameter too high could actually cause the
  // insertion operations to slow down because of the overhead of managing
  // a large hash table. This parameter should probably not exceed 4.0.
  void setInsertSpeed(const float speed);

  //// returns hash table dimension
  unsigned short dimension() const;

  //// returns hash table resolution or cube size.
  float resolution(unsigned short axis=0) const;

  //// Returns hash table insertion speed. See setInsertSpeed.
  float insertSpeed() const;

private:
  GeomHashParams params;
};

/*******************************************
 * Template class GeomHash implementation. *
 *******************************************/

template<class KeyT, class DataT, class BucketT>
GeomHash<KeyT, DataT, BucketT>::GeomHash(const unsigned short d,
                                         const float cs)
  : HashLink<KeyT, DataT, BucketT>(), params(d, cs, 1.0) {
}


template<class KeyT, class DataT, class BucketT>
GeomHash<KeyT, DataT, BucketT>::GeomHash(const unsigned short d,
                                         const std::vector<float>& cs)
  : HashLink<KeyT, DataT, BucketT>(), params(d, cs, 1.0) {
}

template<class KeyT, class DataT, class BucketT>
GeomHash<KeyT, DataT, BucketT>::GeomHash(const GeomHash<KeyT,DataT,BucketT>& gh)
  : HashLink<KeyT, DataT, BucketT>(gh,gh.params.dimension-1),params(gh.params){}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::insert(const KeyT& key, const DataT& data)
{
  // dim-1 for starting with the last axis
  HashLink<KeyT, DataT, BucketT>::put(key, data, params, params.dimension-1);
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::query(const KeyT &key, const float radius,
                                           const NormType ballType, Result& result) const
{
  switch (ballType) {
  case euclideanNorm:
    this->getEuclidean(key, radius, params, params.dimension-1, result);
    break;
  case manhattanNorm:
    this->getManhattan(key, radius, params, params.dimension-1, result);
    break;
  default:
    this->getInfinite(key, radius, params, params.dimension-1, result);
    break;
  }
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::query(const KeyT& key, const float radius,
                                           HashResult<DataT, BucketT>& result) const
{
  this->getInfinite(key, radius, params, params.dimension-1, result);
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::query(const KeyT& key, const std::vector<float>& radii,
                                           HashResult<DataT, BucketT>& result) const
{
  HashLink<KeyT, DataT, BucketT>::get(key, radii, params, params.dimension-1, result);
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::query(const KeyT& key, HashResult<DataT, BucketT>& result) const
{
  HashLink<KeyT, DataT, BucketT>::get(key, params.cubeSizes, params, params.dimension-1, result);
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::query(Result& result,bool includeEmpty) const
{
  BucketsPointerList bresult;

  if (!includeEmpty)
    this->getAllNotEmpty(params.dimension-1, &bresult);
  else
    this->getAll(params.dimension-1, &bresult);

  for(int i=0;i<(int)bresult.size();++i)
    result.push_back((Bucket*)bresult[i]);
}


template<class KeyT, class DataT, class BucketT>
typename GeomHash<KeyT, DataT, BucketT>::BucketsPointerList*
GeomHash<KeyT, DataT, BucketT>::getBuckets(bool includeEmpty) const
{
  BucketsPointerList *result=new BucketsPointerList();
  if (!includeEmpty)
    this->getAllNotEmpty(params.dimension-1, result);
  else
    this->getAll(params.dimension-1 , result);
  return result;
}

template<class KeyT, class DataT, class BucketT>
const BucketT* GeomHash<KeyT, DataT, BucketT>::bucket(const KeyT& key) const
{
  return getBucket(key, params, params.dimension-1);
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::erase()
{
  HashLink<KeyT, DataT, BucketT>::erase(params.dimension-1);
}

template<class KeyT, class DataT, class BucketT>
void GeomHash<KeyT, DataT, BucketT>::setInsertSpeed(const float speed)
{
  params.resizeCoefficient = max2(0.0f, speed);
}

template<class KeyT, class DataT, class BucketT>
unsigned short GeomHash<KeyT, DataT, BucketT>::dimension() const
{
  return params.dimension;
}

template<class KeyT, class DataT, class BucketT>
float GeomHash<KeyT, DataT, BucketT>::resolution(unsigned short axis) const
{
  return params.cubeSizes[axis];
}

template<class KeyT, class DataT, class BucketT>
float GeomHash<KeyT, DataT, BucketT>::insertSpeed() const
{
  return params.resizeCoefficient;
}

#endif
