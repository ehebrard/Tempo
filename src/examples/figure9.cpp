/************************************************
 * Tempodisjunctive scheduling solver
 *
 * Copyright 2024 Emmanuel Hebrard
 *
 * Tempo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 * Tempo is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tempo.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************/

#include <iostream>
#include <vector>

#include "Solver.hpp"

using namespace tempo;

int main() {

  Solver<> S;

 
    auto schedule = S.newInterval(0, 12, 0, 0, 0, 12);

    auto sx{S.newNumeric()};
    auto x{S.between(sx, sx + 3)};
    auto dx{S.newNumeric(1, 1)};
    
    auto sy{S.newNumeric()};
    auto y{S.between(sy, sy + 1)};
    auto dy{S.newNumeric(1, 1)};
    
    auto sz{S.newNumeric()};
    auto z{S.between(sz, sz + 3)};
    auto dz{S.newNumeric(1, 1)};

    auto sw{S.newNumeric()};
    auto w{S.between(sw, sw + 3)};
    auto dw{S.newNumeric(1, 1)};

    auto su{S.newNumeric()};
    auto u{S.between(su, su + 3)};
    auto du{S.newNumeric(1, 1)};

    auto sv{S.newNumeric()};
    auto v{S.between(sv, sv + 3)};
    auto dv{S.newNumeric(1, 1)};

    S.post(su >= 0);
    S.post(sv >= 0);
    S.post(sw >= 0);
    S.post(sz >= 0);
    
    S.post(sx >= 4);
    S.post(sy >= 4);

    
    S.post(x.end <= 8);
    S.post(y.end <= 8);
    
    
    S.post(u.end <= 6);
    S.post(v.end <= 6);
    S.post(w.end <= 6);
    
    S.post(z.end <= 12);
    
    

  S.post(
      Cumulative<int>(schedule, S.newConstant(2), {u, v, w, x, y, z}, {du, dv, dw, dx, dy, dz}));

  S.initializeSearch();
  S.propagate();

  std::cout << u << " (" << du.min(S) << "): " << u.start.min(S) << ".."
            << u.end.max(S) << std::endl;
    
    std::cout << v << " (" << dv.min(S) << "): " << v.start.min(S) << ".."
              << v.end.max(S) << std::endl;
    
    std::cout << w << " (" << dw.min(S) << "): " << w.start.min(S) << ".."
              << w.end.max(S) << std::endl;
    
    std::cout << x << " (" << dx.min(S) << "): " << x.start.min(S) << ".."
              << x.end.max(S) << std::endl;
    
    std::cout << y << " (" << dy.min(S) << "): " << y.start.min(S) << ".."
              << y.end.max(S) << std::endl;
    
    std::cout << z << " (" << dz.min(S) << "): " << z.start.min(S) << ".."
              << z.end.max(S) << std::endl;
     
}
