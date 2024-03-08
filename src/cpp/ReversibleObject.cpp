
#include "ReversibleObject.hpp"

using namespace tempo;

BacktrackEnvironment* ReversibleObject::env = new BacktrackEnvironment();

BacktrackEnvironment::BacktrackEnvironment() { stamps.push_back(0); }

int BacktrackEnvironment::level() const {
  return static_cast<int>(stamps.size()) - 1;
}

void BacktrackEnvironment::subscribe(ReversibleObject *o) {
  subscribers.push_back(o);
}

void BacktrackEnvironment::save(ReversibleObject *o) { trail.push_back(o); }
void BacktrackEnvironment::save() {

  // std::cout << "SAVE @" << level() << std::endl;

  *(stamps.rbegin()) = trail.size();
  stamps.push_back(trail.size());
  for (auto o : subscribers)
    o->checkpoint();
}
// void BacktrackEnvironment::override() { *(stamps.rbegin()) = trail.size(); }
void BacktrackEnvironment::restore(int lvl) {

  if (lvl < level())
    for (auto o : subscribers)
      o->undo();

  size_t stamp = stamps[lvl];
  while (trail.size() > stamp) {
    trail.back()->undo();
    trail.pop_back();
  }
  stamps.resize(lvl + 1);

  // std::cout << "RESTORE TO @" << level() << std::endl;

}

void BacktrackEnvironment::print() const {
  std::cout << "levels:";
  for (auto i{0}; i < level(); ++i)
    std::cout << std::setw(5) << i;
  std::cout << "  cur" << std::endl;
  std::cout << "stamps:";
  for (auto s : stamps)
    std::cout << std::setw(5) << s;
  std::cout << " | " << trail.size() << std::endl;
}

void BacktrackEnvironment::reset() {
    trail.clear();
    stamps.clear();
    subscribers.clear();
    stamps.push_back(0);
}


// void ReversibleObject::undo() = 0;
void ReversibleObject::save() { local_env->save(this); }

void ReversibleObject::checkpoint() {}




ReversibleBool::ReversibleBool(BacktrackEnvironment *e) : ReversibleObject(e) {}

void ReversibleBool::undo() {
    value = Unknown;
}

ReversibleBool &ReversibleBool::operator=(const boolean_state &rhs) {
    assert(value == Unknown);
    
    ReversibleObject::save();
    value = rhs;

    return *this;
}

boolean_state ReversibleBool::val() const { return value; }

std::ostream &ReversibleBool::display(std::ostream &os) const {
    os << value;
    return os;
}
