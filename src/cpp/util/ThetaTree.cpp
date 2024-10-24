

#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>

#include "util/ThetaTree.hpp"

//string pretty(const int x) {
//    std::stringstream ss;
//    if(x == Constant)
//}


std::ostream &tempo::operator<<(std::ostream &os, const tempo::ThetaTree &x) {
  return x.display(os);
}

std::ostream &tempo::operator<<(std::ostream &os, const tempo::ThetaTree *x) {
  return (x ? x->display(os) : os);
}

int tempo::parent(const int x) { return x / 2; }

using namespace tempo;

const int ThetaNode::infbound = std::numeric_limits<int>::min();

ThetaNode::ThetaNode() { ThetaNode::clear(); }

void ThetaNode::clear() {
  earliest_start = 0;
  duration = 0;
  bound = ThetaNode::infbound;
  gray_est = 0;
  gray_duration = 0;
  gray_bound = ThetaNode::infbound;
  responsibleDuration = -1;
  responsibleBound = -1;
}

ThetaNode::~ThetaNode() {}

ThetaTree::ThetaTree(const size_t ntasks) {
  N = 1;
  while (N < ntasks) {
    N *= 2;
  }
  auto number_of_nodes = N + ntasks;
  node = std::vector<ThetaNode>(number_of_nodes, ThetaNode());
}

ThetaTree::ThetaTree() {}

void ThetaTree::remove(int i) {
  node[N + i].clear();
  update(N + i);
  update_gray(N + i);
}

int ThetaTree::getEst() const { return node[1].earliest_start; }

int ThetaTree::grayEst() const { return node[1].gray_est; }

int ThetaTree::grayBound() const { return node[1].gray_bound; }

int ThetaTree::getResponsible() const { return node[1].responsibleBound; }

void ThetaTree::update(int i) {
  while (i > 1) {
    i = parent(i);
    int l = leftChild(i);
    int r = rightChild(i);
    node[i].duration = node[l].duration + node[r].duration;
    node[i].bound = node[r].bound;
    if (node[i].bound < node[l].bound + node[r].duration) {
      node[i].bound = node[l].bound + node[r].duration;
      node[i].earliest_start = node[l].earliest_start;
    } else {
      node[i].earliest_start = node[r].earliest_start;
    }
  }
}

void ThetaTree::update_gray(int i) {
  while (i > 1) {
    i = parent(i);
    int l = leftChild(i);
    int r = rightChild(i);

    // duration: case where we take the gray task from the right
    auto ldrgd{node[l].duration + node[r].gray_duration};

    // duration: case where we take the gray task from the left
    auto rdlgd{node[l].gray_duration + node[r].duration};

    // bound in case of a gap
    auto rge{node[r].gray_bound};

    // carry-over bound: case where we take the gray task from the right
    auto lergd{node[l].bound + node[r].gray_duration};

    // carry-over bound: case where we take the gray task from the left
    auto lgerd{node[l].gray_bound + node[r].duration};

    node[i].gray_duration = std::max(ldrgd, rdlgd);

    if (rge >= lergd) {
      if (rge >= lgerd) {
        //              assert(node[r].responsibleBound >= 0);
        node[i].responsibleBound = node[r].responsibleBound;
        node[i].gray_bound = rge;
        node[i].gray_est = node[r].gray_est;
      } else {
        //              assert(node[l].responsibleBound >= 0);
        node[i].responsibleBound = node[l].responsibleBound;
        node[i].gray_bound = lgerd;
        node[i].gray_est = node[l].gray_est;
      }
    } else {
      if (lgerd <= lergd) {
        //              assert(node[r].responsibleDuration >= 0);
        node[i].responsibleBound = node[r].responsibleDuration;
        node[i].gray_bound = lergd;
        node[i].gray_est = node[l].earliest_start;
      } else {
        //              assert(node[l].responsibleBound >= 0);
        node[i].responsibleBound = node[l].responsibleBound;
        node[i].gray_bound = lgerd;
        node[i].gray_est = node[l].gray_est;
      }
    }

    //    node[i].gray_bound = std::max(std::max(rge, lgerd), lergd);
    assert(node[i].gray_bound == std::max(std::max(rge, lgerd), lergd));

    if (node[i].gray_duration == ldrgd) {
      node[i].responsibleDuration = node[r].responsibleDuration;
    } else if (node[i].gray_duration == rdlgd) {
      node[i].responsibleDuration = node[l].responsibleDuration;
    }

    //    if (node[r].responsibleBound >= 0 and node[i].gray_bound == rge) {
    //      node[i].responsibleBound = node[r].responsibleBound;
    //    } else if (node[l].responsibleBound >= 0 and node[i].gray_bound ==
    //    lgerd) {
    //      node[i].responsibleBound = node[l].responsibleBound;
    //    } else if (node[r].responsibleDuration >= 0 and
    //               node[i].gray_bound == lergd) {
    //      node[i].responsibleBound = node[r].responsibleDuration;
    //    }
  }
}

void ThetaTree::insert(const int i, const int est, const int dur) {
  auto l{i + N};
  node[l].gray_est = node[l].earliest_start = est;
  node[l].gray_duration = node[l].duration = dur;
  node[l].gray_bound = node[l].bound = est + dur;
  update(l);
}

void ThetaTree::paint_gray(const int i, const int a) {
  auto l{i + N};
  assert(node[l].gray_duration == node[l].duration);
  assert(node[l].gray_bound == node[l].bound);
  assert(node[l].gray_est == node[l].earliest_start);

  node[l].earliest_start = 0;
  node[l].duration = 0;
  node[l].bound = ThetaNode::infbound;

  node[l].responsibleDuration = a;
  node[l].responsibleBound = a;
  update(l);
  update_gray(l);
}

void ThetaTree::clear() {
  for (ThetaNode &n : node) {
    n.clear();
  }
}

int ThetaTree::leftChild(int x) { return 2 * x; }

int ThetaTree::rightChild(int x) {
  unsigned r = 2 * x + 1;
  if (r >= node.size()) {
    return 0;
  }
  return r;
}

int ThetaTree::getBound() { return node[1].bound; }

ThetaTree::~ThetaTree() {}

void ThetaTree::printNodeDuration(std::ostream &os, const int i) const {
  if (node[i].duration == 0) {
    if (node[i].responsibleDuration >= 0) {
      os << "gdu=" << std::left << std::setw(4) << node[i].gray_duration;
    } else {
      os << "   --   ";
    }
  } else if (node[i].responsibleDuration >= 0) {
    os << std::setfill('.') << std::left << std::setw(4) << node[i].duration
       << std::right << std::setw(4) << node[i].gray_duration
       << std::setfill(' ');
  } else {
    os << "dur=" << std::left << std::setw(4) << node[i].duration;
  }
}

void ThetaTree::printNodeStart(std::ostream &os, const int i) const {
  if (node[i].duration == 0) {
    if (node[i].responsibleBound >= 0) {
      os << "gst=" << std::left << std::setw(4) << node[i].gray_est;
    } else {
      os << "   --   ";
    }
  } else {
    if (node[i].responsibleBound >= 0) {
      os << std::setfill('.') << std::left << std::setw(4)
         << node[i].earliest_start << std::right << std::setw(4)
         << node[i].gray_est;
      os << std::setfill(' ');
      //          os << "gst=" << std::left << std::setw(4) << node[i].gray_est;
    } else {
      os << "est=" << std::left << std::setw(4) << node[i].earliest_start;
    }
  }
}

void ThetaTree::printNodeBound(std::ostream &os, const int i) const {
  if (node[i].duration == 0) {
    if (node[i].responsibleBound >= 0) {
      os << "get=" << std::left << std::setw(4);
      if (node[i].gray_bound <= ThetaNode::infbound / 2)
        os << "-oo";
      // os << ThetaNode::infbound;
      else
        os << node[i].gray_bound;
    } else {
      os << "   --   ";
    }
  } else if (node[i].responsibleBound >= 0) {
    os << std::setfill('.') << std::left << std::setw(4) << node[i].bound
       << std::right << std::setw(4);
    if (node[i].gray_bound <= ThetaNode::infbound / 2)
      os << "-oo";
    else
      os << node[i].gray_bound;
    os << std::setfill(' ');
  } else {
    os << "ect=" << std::left << std::setw(4) << node[i].bound;
  }
}

std::ostream &ThetaTree::display(std::ostream &os) const {
  int gap{2};
  int skip{1};
  int first = N;
  int last = node.size();
  int sz = 12;
  do {
    os << std::setw((gap + sz * (skip - 1)) / 2) << " ";
    printNodeStart(os, first);
    std::cout << std::setw(4) << "";
    for (auto i{first + 1}; i < last; ++i) {
      os << std::setw(gap + sz * (skip - 1)) << " ";
      printNodeStart(os, i);
      std::cout << std::setw(4) << "";
    }
    os << std::endl;
    os << std::setw((gap + sz * (skip - 1)) / 2) << " ";
    printNodeDuration(os, first);
    std::cout << "(" << std::setw(2) << node[first].responsibleDuration << ")";
    for (auto i{first + 1}; i < last; ++i) {
      os << std::setw(gap + sz * (skip - 1)) << " ";
      printNodeDuration(os, i);
      std::cout << "(" << std::setw(2) << node[i].responsibleDuration << ")";
    }
    os << std::endl << std::setw((gap + sz * (skip - 1)) / 2) << " ";
    printNodeBound(os, first);
    std::cout << "(" << std::setw(2) << node[first].responsibleBound << ")";
    for (auto i{first + 1}; i < last; ++i) {
      os << std::setw(gap + sz * (skip - 1)) << " ";
      printNodeBound(os, i);
      std::cout << "(" << std::setw(2) << node[i].responsibleBound << ")";
    }
    os << std::endl; //<< first << "-" << last << std::endl;
    // if (first + 1 >= last)
    //   break;
    if (first == 1)
      break;

    first /= 2;
    last /= 2;
    gap *= 2;
    skip *= 2;
    // std::cout << first << ".." << last << std::endl;
  } while (true);

  return os;
}
