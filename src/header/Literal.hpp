
#ifndef __TEMPO_LITERAL_HPP
#define __TEMPO_LITERAL_HPP



namespace tempo {


/**********************************************
* Literal
**********************************************/

using var_t = uint32_t;

template<typename T>
struct Literal {
    Literal() = default;
    Literal(const bool sign, const var_t x, const T v=0);
    
    bool sign() const;
    
    var_t variable() const;
    T value() const;
    
    T _value_;
    uint32_t _data_;
};


template<typename T>
Literal<T>::Literal(const bool sign, const var_t x, const T v) {
    _value_ = v;
    _data_ = sign + 2*x;
}

template<typename T>
bool Literal<T>::sign() const {
    return (_data_ & 1);
}

template<typename T>
var_t Literal<T>::variable() const {
    return _data_ / 2;
}

template<typename T>
T Literal<T>::value() const {
    return _value_;
}


}

#endif // __Literal_HPP
