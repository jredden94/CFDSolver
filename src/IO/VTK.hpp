#pragma once

#include <fstream>
#include <iostream>
#include "../DualGrid/DualGrid.hpp"
#include "../Geometry/Cell.hpp"
#include "../Solver/Solver.hpp"
#include "../Common/Config.hpp"

namespace VTK {

    using namespace std;

    void WriteVTK(const string &filename, const Grid &grid, const Solver &solver);

    int GetVTKCellType(const Cell::Type& cellType);
    int GetVTKOffset(const Cell::Type& cellType);
    void InsertIndent(ofstream &file, size_t nIndents = 1);
    string Indent(size_t nIndents = 1);
}
