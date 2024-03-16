
#ifndef _TEMPO_THETATREE_HPP
#define _TEMPO_THETATREE_HPP

#include <vector>

namespace tempo {

int parent(const int x);

class ThetaNode {
public:
  static const int infbound;

  int earliest_start;
  int duration;
  int bound;
  int gray_duration;
  int gray_bound;
  int responsibleDuration;
  int responsibleBound;
  ThetaNode();
  ~ThetaNode();
  void clear();
};

class ThetaTree {
private:
  std::vector<ThetaNode> node;
  size_t N;

public:
  ThetaTree();
  ThetaTree(const size_t ntask);
  ~ThetaTree();
  void insert(const int i, const int est, const int dur);
  void paint_gray(const int i, const int a);
  void remove(int i);
  void update(int i);
  void update_gray(int i);
  int leftChild(int x);
  int rightChild(int x);
  void clear();
  int getBound();
  int grayBound() const;
  int getResponsible() const;

  void printNodeDuration(std::ostream &os, const int i) const;
  void printNodeBound(std::ostream &os, const int i) const;
  std::ostream &display(std::ostream &os) const;
};

std::ostream &operator<<(std::ostream &os, const ThetaTree &x);
std::ostream &operator<<(std::ostream &os, const ThetaTree *x);

} // namespace tempo

#endif // _TEMPO_THETATREE_HPP
