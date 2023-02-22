#ifndef _HashResult_h
#define _HashResult_h

#include <vector>
#include <cstddef>

template<class DataT, class Bucket = std::vector<DataT> >
/*
CLASS
  HashIterator

  Template defines an STL-like iterator over geometric hash table query
  results. (See GeomHash). The iterator is used within the HashResult
  container.

KEYWORDS
  geometric, hash, iterator, query, container

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL

GOALS
  The HashIterator and HashResult families are used together. The HashIterator
  is an STL-like iterator class serving the HashResult class which is an
  STL-like container class.

  A HashResult is formed when a query is posted to a geometric hash table
  (GeomHash). A HashResult holds a set of pointers into hash buckets within
  the hash table that are the query's results. The HashIterator is
  of a const iterator type and it does not allow changes to its members
  because those changes will show up in the hash table itself. The HashIterator
  allows the user to iterate foreward and backward over the set of values
  within the hash table. This set is a query result.

USAGE
  HashResult is used like an STL iterator. This means iterating through the
  list of results can be done using the ++, -- operators and accessing the
  results pointed to by the iterators is done using the * operator. The
  HashResult supplies the begin() and end() methods.

  For an example of using both HashIterator and HashResult see GeomHash.
*/
class HashIterator
{
public:

  // GROUP: STL compliance defintions
  // Some more may be added.
  typedef HashIterator<DataT> const_iterator;

  typedef DataT value_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;

  // GROUP: typedefs & contructors.

  //// Default constructor. Seldom used.
  HashIterator();

  //// Constructor recieves pointer to Bucket and an index within the
  // Bucket. A Bucket is an STL vector and idx is the index within that
  // vector.
  HashIterator(Bucket* const* ptr);

  //// Copy constructor.
  HashIterator(const HashIterator& hi);

  //// Destructor.
  ~HashIterator();

  // GROUP: Equality operators.

  //// Equality operator. Hash iterators from different queries but pointing
  // to the same element in the hash table are equal. If the hash iterators
  // are at the end point of the HashResult container then the equality
  // operator will return true only if both hash iterators were created
  // from the same HashResult object.
  bool operator== (const HashIterator &hi2) {
    if ((*data_ptr != NULL) && (*hi2.data_ptr != NULL))
      return (index == hi2.index);
    else
      return (data_ptr == hi2.data_ptr);
  }

  bool operator!= (const HashIterator &hi2) { return !((*this)==hi2);}

  // GROUP: Iterating operators.

  //// Contents of. Returns data pointed to by const_iterator.
  DataT& operator*();

  //// Contents of. Returns data pointed to by const_iterator.
  const DataT& operator*() const;

  //// Iterator increment operator. prefix.
  HashIterator& operator++();

  //// Iterator decrement operator. prefix.
  HashIterator& operator--();

  //// Iterator increment operator. postfix.
  HashIterator operator++(int);

  //// Iterator decrement operator. postfix.
  HashIterator operator--(int);

  HashIterator nextBucket() const;
  HashIterator prevBucket() const;

private:
  Bucket* const* data_ptr;  // pointer to Bucket - a vector of data elements.
  typename Bucket::iterator index;
};


template<class DataT, class Bucket = std::vector<DataT> >
/*
CLASS
  HashResult

  Template defines an STL-like container class used to hold geometric hash
  table query results. (See GeomHash). The results are iterated through using
  the HashIterator STL-like iterator class.

KEYWORDS
  geometric, hash, iterator, query, container

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

GOALS
  The class defines an STL-like container that holds the geometric hash table
  query results. A HashResult is formed when a query is posted to a geometric
  hash table (GeomHash). A HashResult holds a set of pointers into hash
  Buckets within the hash table that are the query's results. The container
  elemenets are constant and no attempt to change them should be made because
  changes to the container elements are changes within the hash table itself.

USAGE
  The begin() and end() methods are used to generate head and tail iterators
  to the result list. The HashIterator-s can then be used to iterate through
  the result list using the ++, -- operators and using the * operator to
  reference the data itself.

  For an example of using the query result container and iterator class see
  GeomHash.
*/
class HashResult
{
public:
  // STL like typedef iterators found in containers.
  typedef HashIterator<DataT, Bucket> iterator;
  typedef DataT value_type;
  typedef value_type* pointer;
  typedef const value_type& const_reference;
  typedef const value_type& reference;

  //// The default size of the Bucket pointers array. The array grows
  // dynamically as more results are collected. The default size and
  // increment size are defined by defSize but could be changed using the
  // constructor. Note that the pointers are to Buckets of elements and not
  // to the elements themselves so the pointer array size should not be as
  // large as the number of results.
  static const int defSize = 32;

  //// Instantiate a HashResult. Allocates space for defSize bucket pointers.
  HashResult();

  //// Instantiate a HashResult. Allocates space for defAlloc bucket pointers.
  HashResult(const unsigned int defAlloc);

  //// Destructor.
  ~HashResult();

  //// Add a Bucket pointer to the dynamic pointer array held by the container.
  // Generally, for normal use of the geometric hash tables this method
  // should not be called. It is used by the HashLink template class to collect
  // query results.
  void push_back(Bucket* bkt);

  //// Returns number of elements in the query result container. It should be
  // stressed that the method does not return the number of Bucket pointers
  // or the size of the pointer array but the number of actual data elements
  // in the query results which could be much larger.
  size_t size() const;

  //// Returns an iterator to the beginning of the query results list.
  iterator begin();

  //// Returns an iterator to the end of the query results list.
  iterator end();

  //// Returns true if container is empty.
  bool empty() const;

  //// Same as size(). Here for STL compliance.
  size_t max_size() const;

  //// Returns first element of the container (the query results-
  // assuming its a Bucket)
  reference front() const;

  //// Returns last element of the container (the query results -
  // assuming its a Bucket).
  reference back() const;

  //// Returns the number of Buckets referenced by the results container.
  unsigned int bucketCount() const;

  //// Clears hash result container. Does not affect the geometric hash table
  // itself. After a clear call a HashResult object may be used for another
  // query. Memory used is not deallocated until object destruction.
  void clear();

private:
  Bucket** buckets;
  unsigned int count;
  unsigned int alloc;
};

/***********************************************
 * Template class HashIterator implementation. *
 ***********************************************/

template<class DataT, class Bucket>
HashIterator<DataT, Bucket>::HashIterator() {}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket>::HashIterator(Bucket* const* ptr)
  : data_ptr(ptr)
{
  if (*data_ptr != NULL)
    index = (*data_ptr)->begin();
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket>::HashIterator(const HashIterator<DataT, Bucket>& hi)
  : data_ptr(hi.data_ptr), index(hi.index) {}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket>::~HashIterator() {}

template<class DataT, class Bucket>
DataT& HashIterator<DataT, Bucket>::operator*()
{
  return *index;
}

template<class DataT, class Bucket>
const DataT& HashIterator<DataT, Bucket>::operator*() const
{
  return *index;
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket>& HashIterator<DataT, Bucket>::operator++()
{
  if (++index == (*data_ptr)->end())
    if (*(++data_ptr) != NULL)
      index = (*data_ptr)->begin();
  return *this;
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket>& HashIterator<DataT, Bucket>::operator--()
{
  if (*data_ptr == NULL || index == (*data_ptr)->begin())
    index = (*(--data_ptr))->end();
  --index;
  return *this;
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket> HashIterator<DataT, Bucket>::operator++(int)
{
  HashIterator temp = *this;
  ++*this;
  return temp;
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket> HashIterator<DataT, Bucket>::operator--(int) {
  HashIterator temp = *this;
  --*this;
  return temp;
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket> HashIterator<DataT, Bucket>::nextBucket() const
{
  return HashIterator(data_ptr+1);
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket> HashIterator<DataT, Bucket>::prevBucket() const
{
  if (*data_ptr != NULL && index != (*data_ptr)->begin())
    return *this;
  else
    return HashIterator(data_ptr-1);
}

// template<class DataT, class Bucket=vector<DataT> >
// inline bool operator==(const HashIterator<DataT, Bucket>& hi1,
// 		       const HashIterator<DataT, Bucket>& hi2)
// {
//   if ((*hi1.data_ptr != NULL) && (*hi2.data_ptr != NULL))
//     return (hi1.index == hi2.index);
//   else
//     return (hi1.data_ptr == hi2.data_ptr);
// }

// template<class DataT, class Bucket=vector<DataT> >
// inline bool operator!=(const HashIterator<DataT, Bucket>& hi1,
// 		       const HashIterator<DataT, Bucket>& hi2)
// {
//   return !(hi1==hi2);
// }

/*********************************************
 * Template class HashResult implementation. *
 *********************************************/

template<class DataT, class Bucket>
HashResult<DataT, Bucket>::HashResult()
  : count(0), alloc(defSize)
{
  buckets = (Bucket**)operator new(alloc*sizeof(Bucket*));
  *buckets = NULL;
}

template<class DataT, class Bucket>
HashResult<DataT, Bucket>::HashResult(const unsigned int defAlloc)
  : count(0), alloc(defAlloc+1)
{
  buckets = (Bucket**)operator new(alloc*sizeof(Bucket*));
  *buckets = NULL;
}

template<class DataT, class Bucket>
HashResult<DataT, Bucket>::~HashResult()
{
  if (buckets) operator delete(buckets);
}

template<class DataT, class Bucket>
void HashResult<DataT, Bucket>::push_back(Bucket* bkt)
{
  if (count == alloc-1) {
    if (count) alloc *= 2;
    else       alloc = defSize;
    Bucket** newBuckets = (Bucket**)operator new(alloc*sizeof(Bucket*));
    for (unsigned int i=0; i<count; ++i) newBuckets[i] = buckets[i];
    if (buckets) operator delete(buckets);
    buckets = newBuckets;
  }
  buckets[count++] = bkt;
  buckets[count] = NULL;
}

template<class DataT, class Bucket>
size_t HashResult<DataT, Bucket>::size() const
{
  size_t sz = 0;
  for (unsigned int i=0; i<count; i++)  sz+= buckets[i]->size();
  return sz;
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket> HashResult<DataT, Bucket>::begin()
{
  return iterator(buckets);
}

template<class DataT, class Bucket>
HashIterator<DataT, Bucket> HashResult<DataT, Bucket>::end()
{
  return iterator(buckets+count);
}

template<class DataT, class Bucket>
bool HashResult<DataT, Bucket>::empty() const
{
  return (count==0);
}

template<class DataT, class Bucket>
size_t HashResult<DataT, Bucket>::max_size() const
{
  return size();
}

template<class DataT, class Bucket>
const DataT& HashResult<DataT, Bucket>::front() const
{
  return *((*buckets)->begin());
}

template<class DataT, class Bucket>
const DataT& HashResult<DataT, Bucket>::back() const
{
  return *(--end());
}

template<class DataT, class Bucket>
unsigned int HashResult<DataT, Bucket>::bucketCount() const
{
  return count;
}

template<class DataT, class Bucket>
void HashResult<DataT, Bucket>::clear()
{
  count=0;
}
#endif
