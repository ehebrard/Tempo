#include "Global.hpp"


std::string tempo::prettyEvent(const event e) {
    if(e == ORIGIN)
        return "ORIGIN";
    if(e == HORIZON)
        return "C_max";
    return std::string(1,etype(e)) + std::to_string(TASK(e));
}

std::string tempo::prettyEventLit(const lit el) {
    event e{EVENT(el)};
    if(e == ORIGIN)
        return "ORIGIN";
    return (SIGN(el) ? "ub(" : "lb(") + prettyEvent(e) + ")";
}

//std::string prettyLiteral(const genlit el);

char tempo::etype(const event x) {
    return x%2 ? 'e' : 's';
}

double tempo::cpu_time(void) {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static unsigned long x = 123456789, y = 362436069, z = 521288629;

void tempo::seed(const unsigned long s) {
  x = s;
  y = 362436069;
  z = 521288629;
}

unsigned long tempo::random(void) { // period 2^96-1
  unsigned long t;
  x ^= x << 16;
  x ^= x >> 5;
  x ^= x << 1;

  t = x;
  x = y;
  y = z;
  z = t ^ x ^ y;

  return z;
}

// 0 is for the origin event, 1 is for the completion event
tempo::task tempo::TASK( event x ) {
		return x/2-1;
}

tempo::event tempo::START( task x ) {
		return 2*x+2;
}

tempo::event tempo::END( task x ) {
		return 2*x+3;
}

tempo::lit tempo::POS(const var x) {
    return 2*x+1;
}

tempo::lit tempo::LIT(const var x, const boolean_state v) {
    return 2*x+v;
}

tempo::lit tempo::NEG(const var x) {
    return 2*x;
}

tempo::var tempo::VAR(const lit l) {
    return l/2;
}

tempo::boolean_state tempo::SIGN(const lit l) {
    return l&1;
}

tempo::lit tempo::NOT(const lit l) {
    return l^1;
}

tempo::genlit tempo::GNOT(const genlit l) {
    return l^2;
}

int tempo::TODIMACS(const lit l) {
    return (SIGN(l) ? VAR(l)+1 : ~VAR(l));
}

tempo::lit tempo::LOWERBOUND(const event x) {
    return 2*x;
}

tempo::lit tempo::UPPERBOUND(const event x) {
    return 2*x+1;
}

tempo::event tempo::EVENT(const lit l) {
    return l/2;
}


tempo::lit_type tempo::LTYPE(const hint h) {
    return h&1;
}
tempo::genlit tempo::BOUND(const lit l) {
    return 2*l;
}
tempo::genlit tempo::EDGE(const lit l) {
    return 2*l+1;
}
tempo::lit tempo::FROM_GEN(const genlit h) {
    return h/2;
}
