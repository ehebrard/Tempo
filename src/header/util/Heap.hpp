
#ifndef __HEAP_HPP__
#define __HEAP_HPP__

namespace heap {

constexpr int left(int i) { return 2 * i + 1; }

constexpr int right(int i) { return 2 * i + 2; }

template <class RandomIt, class Compare>
int percolate_down(RandomIt begin, RandomIt end, int idx, Compare comp) {

  int size = std::distance(begin, end);
  int current = idx;
  int child;

  // limits
  while ((child = 2 * current + 1) < size) {
    RandomIt itr_child = begin + child;
    RandomIt itr_current = begin + current;

    // check if there is a right child and it is lower than the left child
    if (child < size - 1 && comp(*(itr_child + 1), *itr_child)) {
      ++itr_child;
      ++child;
    }

    // stop when the parent is lower than the child
    //    if (comp(*itr_current, *itr_child)) {

    // stop when the child is not lower than the parent
    if (not comp(*itr_child, *itr_current)) {
      break;
    }

    // or advance
    std::swap(*itr_current, *itr_child);
    current = child;
  }

  return current;
}

template <class RandomIt, class Array, class Compare>
int percolate_down(RandomIt begin, RandomIt end, int idx, Array &indexing,
                   Compare comp) {

  int size = std::distance(begin, end);
  int current = idx;
  int child;

  // limits
  while ((child = 2 * current + 1) < size) {
    RandomIt itr_child = begin + child;
    RandomIt itr_current = begin + current;

    // check if there is a right child and it is lower than the left child
    if (child < size - 1 && comp(*(itr_child + 1), *itr_child)) {
      ++itr_child;
      ++child;
    }

    // stop when the parent is lower than the child
    //    if (comp(*itr_current, *itr_child)) {

    // stop when the child is not lower than the parent
    if (not comp(*itr_child, *itr_current)) {
      break;
    }

    // or advance

    indexing[*itr_current] = child;
    indexing[*itr_child] = current;

    std::swap(*itr_current, *itr_child);
    current = child;
  }

  return current;
}

template <class RandomIt, class Compare>
int percolate_up(RandomIt begin, int idx, Compare comp) {
  int current = idx;
  int parent;

  while (current > 0) {
    parent = (current - 1) / 2;
    RandomIt itr_parent = begin + parent;
    RandomIt itr_current = begin + current;

    // stop when the parent is lower than the child
    if (comp(*itr_parent, *itr_current)) {
      break;
    }

    // or back off
    std::swap(*itr_current, *itr_parent);
    current = parent;
  }

  return current;
}

template <class RandomIt, class Array, class Compare>
int percolate_up(RandomIt begin, int idx, Array &indexing, Compare comp) {
  int current = idx;
  int parent;

  while (current > 0) {
    parent = (current - 1) / 2;
    RandomIt itr_parent = begin + parent;
    RandomIt itr_current = begin + current;

    // stop when the parent is lower than the child
    if (comp(*itr_parent, *itr_current)) {
      break;
    }

    // or back off
    indexing[*itr_current] = parent;
    indexing[*itr_parent] = current;
    std::swap(*itr_current, *itr_parent);
    current = parent;
  }

  return current;
}

template <class RandomIt, class Compare>
void heapify(RandomIt begin, RandomIt end, Compare comp) {
  int size = std::distance(begin, end);
  for (int i = size / 2 - 1; i > -1; --i) {
    percolate_down(begin, end, i, comp);
  }
}

template <class RandomIt, class Compare>
void sort(RandomIt begin, RandomIt end, Compare comp) {
  heapify(begin, end, comp);
  while (begin != --end) {
    std::swap(*begin, *end);
    percolate_down(begin, end, 0, comp);
  }
}

template <class RandomIt> void sort(RandomIt begin, RandomIt end) {
  std::less<decltype(*begin)> less;
  sort(begin, end, less);
}

template <class RandomIt, class Compare>
void remove_min(RandomIt begin, RandomIt end, Compare comp) {
  std::swap(*begin, *(end - 1));
  percolate_down(begin, end - 1, 0, comp);
}

template <class RandomIt, class Array, class Compare>
void remove_min(RandomIt begin, RandomIt end, Array &indexing, Compare comp) {
  indexing[*(end - 1)] = 0;
  indexing[*begin] = (end - begin - 1);
  std::swap(*begin, *(end - 1));
  percolate_down(begin, end - 1, 0, indexing, comp);
}

} // namespace heap

template <class T, class Comparator> class Heap {

  std::vector<T> heap;
  Comparator comp;

public:
  template <class RandomIt>
  void initialise(RandomIt begin, RandomIt end, Comparator c) {
    heap.reserve(end - begin);
    for (auto x{begin}; x != end; ++x) {
      heap.push_back(*x);
    }
    comp = c;
    heap::heapify(heap.begin(), heap.end(), comp);
  }

  Heap() {}

  template <class RandomIt> Heap(RandomIt begin, RandomIt end) {
    std::less<decltype(*begin)> less;
    initialise(begin, end, less);
  }

  template <class RandomIt> Heap(RandomIt begin, RandomIt end, Comparator c) {
    initialise(begin, end, c);
  }

  size_t size() { return heap.size(); }

  void add(const T x) {
    heap.push_back(x);
    heap::percolate_up(heap.begin(), heap.size() - 1, comp);
  }

  T pick() { return *(heap.begin()); }

  T extractMin() {
    auto x{*heap.begin()};
    heap.pop_back();
    std::swap(*(heap.begin()), *(heap.end()));
    heap::percolate_down(heap.begin(), heap.end(), 0, comp);
    return x;
  }

  bool empty() const { return heap.empty(); }

  typename std::vector<T>::iterator begin() { return heap.begin(); }
  typename std::vector<T>::iterator end() { return heap.end(); }
};

#endif
