#pragma once

#include <vector>
#include <string>
#include "Types.h"

#ifndef COORDLISTMATRIX_H
#define COORDLISTMATRIX_H

class CoordListMatrix{
  public:
    CoordListMatrix(const std::string &filename, int &M, int &N);
    std::vector<Coord>& getCoords();

  private:
    std::vector<Coord> coords;

};

#endif //COORDLISTMATRIX_H
