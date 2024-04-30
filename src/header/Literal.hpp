
#ifndef __TEMPO_LITERAL_HPP
#define __TEMPO_LITERAL_HPP



namespace tempo {


/**********************************************
* Literal
**********************************************/

using var_t = std::uint32_t;
using info_t = std::uint32_t;

//// this
// struct LitInfo {
//     LitInfo(uint32_t d);
//
//     operator info_t() const;
//     uint32_t _data_;
// };

// type_t event_type(const bool sign, const var_t x);

enum bound { lower = 0, upper = 1 };

template<typename T>
struct Literal {
    Literal() = default;
    Literal(const info_t d);
    Literal(const bool sign, const var_t x, const T v=0);

    operator info_t() const;

    bool sign() const;
    
    var_t variable() const;
    T value() const;

    T _value_{0};
    uint32_t _data_{0};
};

template <typename T> Literal<T>::Literal(const info_t d) : _data_(d) {}

template <typename T>
Literal<T>::Literal(const bool sign, const var_t x, const T v)
    : _data_(sign + 2 * x) {
  _value_ = v;
}

template <typename T> Literal<T>::operator info_t() const { return _data_; }

template <typename T> bool Literal<T>::sign() const { return (_data_ & 1); }

template<typename T>
var_t Literal<T>::variable() const {
    return _data_ / 2;
}

template<typename T>
T Literal<T>::value() const {
    return _value_;
}


template <typename T> Literal<T> operator~(const Literal<T> l) {
  return Literal<T>(l._data_ ^ 1);
}

}

#endif // __Literal_HPP
