
#ifndef __TEMPO_TASK_HPP
#define __TEMPO_TASK_HPP


#include "Global.hpp"


namespace tempo {

template<typename T> class Scheduler;

/**********************************************
 * Task
 **********************************************/

template<typename T>
class Event {

public:
    Event(const int i, const T o=0) : _id_(i), _offset(o) {}
    
    T earliest(Scheduler<T>& sc) const;
    T latest(Scheduler<T>& sc) const;
    
    BoundConstraint<T> after(const T t) const;
    BoundConstraint<T> before(const T t) const;
    
    DistanceConstraint<T> after(const Event<T>& e, const T t=0) const;
    DistanceConstraint<T> before(const Event<T>& e, const T t=0) const;
    
    
    int id() const {return _id_;}
    T offset() const {return _offset;}
    
private:
    
    const int _id_;
    const T _offset;
    
};

template<typename T>
class Task {
public:
  Task(Scheduler<T> &sc);
  Task(Scheduler<T> &sc, const T mindur, const T maxdur);

  T getEarliestStart() const;
  T getLatestStart() const;
  T getEarliestEnd() const;
  T getLatestEnd() const;

  bool mustExist() const;
  bool cannotExist() const;

  T minDuration() const;
  T maxDuration() const;

  event getStart() const;
  event getEnd() const;

  //    BoundConstraint<T> startAfter(const T t) const;
  //    BoundConstraint<T> endAfter(const T t) const;
  //    BoundConstraint<T> startBefore(const T t) const;
  //    BoundConstraint<T> endBefore(const T t) const;

  int id() const;
  bool operator==(const Task<T> &t) const;
    
    std::ostream &display(std::ostream &os) const;
//    std::string& print(const event e) const;

private:
    int _id_{-1};
    
    Scheduler<T>& sched;
    
public:
    Event<T> start;
    Event<T> end;
    
private:
    T min_duration;
    T max_duration;
    var optional{NoVar};
};

template <typename T>
Task<T>::Task(Scheduler<T> &sc)
    : _id_(-1), sched(sc), start(START(_id_), 0), end(END(_id_), 0) {}

template <typename T>
Task<T>::Task(Scheduler<T> &sc, const T mindur, const T maxdur)
    : _id_(static_cast<int>(sc.numTask())), sched(sc), 
start(static_cast<event>(sched.numEvent()), 0),
end(static_cast<event>(sched.numEvent() + (mindur != maxdur)), mindur)
//end(static_cast<event>(sched.numEvent() + 1), 0)
// end((mindur == maxdur ? START(_id_) : END(_id_)), (mindur == maxdur ? mindur : 0))
{

  //        std::cout << _id_ << std::endl;
  //
  //        std::cout << (getEnd() + 1) << std::endl;
  //
  //
  min_duration = mindur;
  max_duration = maxdur;

  //    std::cout << (getEnd() + 1) << std::endl;

  sched.resize(getEnd() + 1);
     
//     std::cout << "resize("<< (getEnd() + 1) << ")\n";

  //        std::cout << sched << std::endl;

  if (start.id() != end.id()) {
      
//      std::cout << "variable duration task: " << *this << std::endl;
      
//    if (max_duration != INFTY)
      sched.newMaximumLag(getStart(), getEnd(), max_duration);
    sched.newMaximumLag(getEnd(), getStart(), -min_duration);
  } 
//  else {
//      std::cout << "fixed duration task: " << *this << std::endl;
//  }
}

//
//template <typename T>
//Task<T>::Task(Scheduler<T> &sc, const T mindur, const T maxdur)
//    : _id_(static_cast<int>(sc.numTask())), sched(sc),
//start(START(_id_), 0),
////      end(END(_id_), 0)
// end((mindur == maxdur ? START(_id_) : END(_id_)), (mindur == maxdur ? mindur : 0))
//{
//
//  //        std::cout << _id_ << std::endl;
//  //
//  //        std::cout << (getEnd() + 1) << std::endl;
//  //
//  //
//  min_duration = mindur;
//  max_duration = maxdur;
//
//  //    std::cout << (getEnd() + 1) << std::endl;
//
//  sched.resize(getEnd() + 1);
//     
//     std::cout << "resize("<< (getEnd() + 1) << ")\n";
//
//  //        std::cout << sched << std::endl;
//
//  if (start.id() != end.id()) {
//      
//      std::cout << "variable duration task: " << *this << std::endl;
//      
////    if (max_duration != INFTY)
//      sched.newMaximumLag(getStart(), getEnd(), max_duration);
//    sched.newMaximumLag(getEnd(), getStart(), -min_duration);
//  } else {
//      std::cout << "fixed duration task: " << *this << std::endl;
//  }
//}

template<typename T>
T Event<T>::earliest(Scheduler<T>& sc) const {
    return sc.lower(_id_) + _offset;
}

template<typename T>
T Event<T>::latest(Scheduler<T>& sc) const {
    return sc.upper(_id_) + _offset;
}

template<typename T>
BoundConstraint<T> Event<T>::after(const T t) const {
    return {LOWERBOUND(_id_), _offset-t};
}

template<typename T>
BoundConstraint<T> Event<T>::before(const T t) const {
    return {UPPERBOUND(_id_), t-_offset};
}




//{START(job[i]), END(job[j]), -job.getTransitionTime(j, i)},
// s_i -> e_j
// e_j - s_i <= - t


template<typename T>
DistanceConstraint<T> Event<T>::after(const Event<T>& e, const T t) const {
//    return {_id_, e.id(), _offset - e.offset() - t};
    return e.before(*this, t);
}

template<typename T>
DistanceConstraint<T> Event<T>::before(const Event<T>& e, const T t) const {
    return {e.id(), _id_, e.offset() - _offset - t};
}

template<typename T>
int Task<T>::id() const {
    return _id_;
}

template<typename T>
bool Task<T>::operator==(const Task<T>& t) const {
    return id() == t.id();
}

template<typename T>
T Task<T>::getEarliestStart() const {
    return start.earliest(sched);
}

template<typename T>
T Task<T>::getLatestStart() const {
    return start.latest(sched);
}

template<typename T>
T Task<T>::getEarliestEnd() const {
    return end.earliest(sched);
//    return sched.lower(end_event);
}

template<typename T>
T Task<T>::getLatestEnd() const {
    return end.latest(sched);
//    return sched.upper(end_event);
}

//template<typename T>
//void Task<T>::setEarliestStart(const T t)  {
//    sched.newPrecedence(ORIGIN, start.id(), t-start.offset());
//}
//
//template<typename T>
//void Task<T>::setLatestStart(const T t) {
//    sched.newMaximumLag(ORIGIN, start.id(), t);
//}
//
//template<typename T>
//void Task<T>::setEarliestEnd(const T t) {
//    sched.newPrecedence(ORIGIN, end_event, t);
//}
//
//template<typename T>
//void Task<T>::setLatestEnd(const T t) {
//    sched.newMaximumLag(ORIGIN, end_event, t);
//}

template<typename T>
bool Task<T>::mustExist() const {
    return true;
}

template<typename T>
bool Task<T>::cannotExist() const {
    return false;
}

template<typename T>
T Task<T>::minDuration() const {
    return min_duration;
}

template<typename T>
T Task<T>::maxDuration() const {
    return max_duration;
}

//template<typename T>
//BoundConstraint<T> Task<T>::startAfter(const T t) const {
////    return {LOWERBOUND(start_event), -t};
//    return start.after(t); // {LOWERBOUND(start.id()), -t};
//}
//
//template<typename T>
//BoundConstraint<T> Task<T>::endAfter(const T t) const {
//    return end.after(t); //{LOWERBOUND(end_event), -t};
//}
//
//template<typename T>
//BoundConstraint<T> Task<T>::startBefore(const T t) const {
//    return start.before(t); //{UPPERBOUND(start_event), t};
//}
//
//template<typename T>
//BoundConstraint<T> Task<T>::endBefore(const T t) const {
//    return end.before(t); //{UPPERBOUND(end_event), t};
//}

template<typename T>
event Task<T>::getStart() const {
    return start.id();
}

template<typename T>
event Task<T>::getEnd() const {
    return end.id();
}


template<typename T>
std::ostream &Task<T>::display(std::ostream &os) const {
    os << "t" << id() << ": [" << start.earliest(sched) << ".." << end.latest(sched) << "]";
    return os;
}

//template<typename T>
//std::string &Task<T>::to_string() const {
//    std::stringstream ss;
//    display(ss);
//    return ss;
//}

//template<typename T>
//std::string& Task<T>::print(const event e) const {
//    std::stringstream ss;
//    if(start.id() == e) {
//        ss << "s" << id();
//        if(start.offset() != 0) {
//            ss << " + " << start.offset();
//        }
//    } else {
//        ss << "s" << id();
//        if(end.offset() != 0) {
//            ss << " + " << end.offset();
//        }
//    }
//    return ss;
//}


template<typename T>
std::ostream &operator<<(std::ostream &os, const Task<T> &x) {
  return x.display(os);
}







//
//
//
//
//template<typename T>
//class TaskImpl {
//public:
//    
////    Task(Scheduler<T>& sc, const event s, const event e, const T mindur, const T maxdur);
////    TaskImpl(const int i);
////    virtual ~TaskImpl() = default;
//    //
//    virtual T getEarliestStart() const = 0;
//    virtual T getLatestStart() const = 0;
//    virtual T getEarliestEnd() const = 0;
//    virtual T getLatestEnd() const = 0;
//    
//    virtual bool mustExist() const = 0;
//    virtual bool cannotExist() const = 0;
//    
//    virtual T minDuration() const = 0;
//    virtual T maxDuration() const = 0;
//    
////    virtual void setEarliestStart(const T t) = 0;
////    void setLatestStart(const T t) ;
////    void setEarliestEnd(const T t) ;
////    void setLatestEnd(const T t) ;
//    
//    virtual event getStart() const = 0;
//    virtual event getEnd() const = 0;
//    
//    int id() const;
//    
//    virtual BoundConstraint<T> startAfter(const T t) const = 0;
//    virtual BoundConstraint<T> endAfter(const T t) const = 0;
//    virtual BoundConstraint<T> startBefore(const T t) const = 0;
//    virtual BoundConstraint<T> endBefore(const T t) const = 0;
//    
////    bool operator==(const Task<T>& t) const;
//    
//};
//
//template<typename T>
//class Task {
//public:
//    
////    Task(Scheduler<T>& sc, const event s, const event e, const T mindur, const T maxdur);
//    Task(Scheduler<T>& sc, const event s, const event e, const T mindur, const T maxdur);
////    ~Task() {delete impl;}
//    //
//    T getEarliestStart() const {return impl->getEarliestStart();}
//    T getLatestStart() const {return impl->getLatestStart();}
//    T getEarliestEnd() const {return impl->getEarliestEnd();}
//    T getLatestEnd() const {return impl->getLatestEnd();}
//    
//    bool mustExist() const {return impl->mustExist();}
//    bool cannotExist() const {return impl->cannotExist();}
//    
//    T minDuration() const {return impl->minDuration();}
//    T maxDuration() const {return impl->maxDuration();}
//    
////    virtual void setEarliestStart(const T t) = 0;
////    void setLatestStart(const T t) ;
////    void setEarliestEnd(const T t) ;
////    void setLatestEnd(const T t) ;
//    
//    event getStart() const {return impl->getStart();}
//    event getEnd() const {return impl->getEnd();}
//    
//    BoundConstraint<T> startAfter(const T t) const {return impl->startAfter(t);}
//    BoundConstraint<T> endAfter(const T t) const {return impl->endAfter(t);}
//    BoundConstraint<T> startBefore(const T t) const {return impl->startBefore(t);}
//    BoundConstraint<T> endBefore(const T t) const {return impl->endBefore(t);}
//    
//    int id() const;
//    bool operator==(const Task<T>& t) const;
//    
//protected:
//    int _id_{-1};
//    TaskImpl<T> *impl;
//};
//
//template<typename T>
//class ClassicTask : public TaskImpl<T> {
//    
//public:
//        ClassicTask(Scheduler<T>& sc, const event s, const event e, const T mindur, const T maxdur);
////    virtual ~ClassicTask() = default;
//        
//         T getEarliestStart() const override;
//         T getLatestStart() const override;
//         T getEarliestEnd() const override;
//         T getLatestEnd() const override;
//        
//         bool mustExist() const override;
//         bool cannotExist() const override;
//        
//         T minDuration() const override;
//         T maxDuration() const override;
//        
//         event getStart() const override;
//         event getEnd() const override;
//        
////        int id() const;
//        
//         BoundConstraint<T> startAfter(const T t) const override;
//         BoundConstraint<T> endAfter(const T t) const override;
//         BoundConstraint<T> startBefore(const T t) const override;
//         BoundConstraint<T> endBefore(const T t) const override;
//    
//private:
//    Scheduler<T>& sched;
//    
//    event start_event;
//    event end_event; // equal to start_evt for fixed duration tasks
//    T offset{0}; // used for fixed duration tasks otherwise 0
//    T min_duration;
//    T max_duration;
//    var optional{NoVar};
//};
//
//template<typename T>
//Task<T>::Task(Scheduler<T>& sc, const event s, const event e, const T mindur, const T maxdur) : _id_(static_cast<int>(sc.numTask())) {
//    impl = new ClassicTask(sc,s,e,mindur,maxdur);
//}
//
//template<typename T>
//int Task<T>::id() const {
//    return _id_;
//}
//
//template<typename T>
//bool Task<T>::operator==(const Task<T>& t) const {
//    return id() == t.id();
//}
//
//
//
//template<typename T>
//ClassicTask<T>::ClassicTask(Scheduler<T>& sc, const event s, const event e, const T mindur, const T maxdur) : sched(sc) {
////    assert(i = static_cast<int>(sched.numTask());
//    start_event = s;
//    end_event = e;
//    min_duration = mindur;
//    max_duration = maxdur;
//}
//
//template<typename T>
//T ClassicTask<T>::getEarliestStart() const {
//    return sched.lower(start_event);
//}
//
//template<typename T>
//T ClassicTask<T>::getLatestStart() const {
//    return sched.upper(start_event);
//}
//
//template<typename T>
//T ClassicTask<T>::getEarliestEnd() const {
//    return sched.lower(end_event);
//}
//
//template<typename T>
//T ClassicTask<T>::getLatestEnd() const {
//    return sched.upper(end_event);
//}
//
////template<typename T>
////void ClassicTask<T>::setEarliestStart(const T t)  {
////    sched.newPrecedence(ORIGIN, start_event, t);
////}
////
////template<typename T>
////void ClassicTask<T>::setLatestStart(const T t) {
////    sched.newMaximumLag(ORIGIN, start_event, t);
////}
////
////template<typename T>
////void ClassicTask<T>::setEarliestEnd(const T t) {
////    sched.newPrecedence(ORIGIN, end_event, t);
////}
////
////template<typename T>
////void ClassicTask<T>::setLatestEnd(const T t) {
////    sched.newMaximumLag(ORIGIN, end_event, t);
////}
//
//template<typename T>
//bool ClassicTask<T>::mustExist() const {
//    return true;
//}
//
//template<typename T>
//bool ClassicTask<T>::cannotExist() const {
//    return false;
//}
//
//template<typename T>
//T ClassicTask<T>::minDuration() const {
//    return min_duration;
//}
//
//template<typename T>
//T ClassicTask<T>::maxDuration() const {
//    return max_duration;
//}
//
//template<typename T>
//BoundConstraint<T> ClassicTask<T>::startAfter(const T t) const {
//    return {LOWERBOUND(start_event), -t};
//}
//
//template<typename T>
//BoundConstraint<T> ClassicTask<T>::endAfter(const T t) const {
//    return {LOWERBOUND(end_event), -t};
//}
//
//template<typename T>
//BoundConstraint<T> ClassicTask<T>::startBefore(const T t) const {
//    return {UPPERBOUND(start_event), t};
//}
//
//template<typename T>
//BoundConstraint<T> ClassicTask<T>::endBefore(const T t) const {
//    return {UPPERBOUND(end_event), t};
//}
//
//template<typename T>
//event ClassicTask<T>::getStart() const {
//    return start_event;
//}
//
//template<typename T>
//event ClassicTask<T>::getEnd() const {
//    return end_event;
//}

}

#endif // __Task_HPP
