#pragma once

#include <vector>
#include <set>
#include <cmath>
#include <fstream>
#include "Node.hpp"
#include "../Geometry/Cell.hpp"
#include "Edge.hpp"
#include "Boundary.hpp"
#include "../Common/Config.hpp"

using namespace std;

class Grid {
    public:
        Grid();
        ~Grid();
        Grid(const Grid&) = delete;
        Grid& operator=(const Grid&) = delete;

        const vector<Node>& Nodes() const;
        const vector<unique_ptr<Cell>>& Cells() const;
        const vector<Edge>& Edges() const;
        const vector<unique_ptr<Boundary>>& Boundaries() const;

        void BuildDualGrid(void);

        void FlipBoundaryNorms(void);

        friend class GridReader;
        friend void ReadGrid2D(Grid&, string);
        friend void ReadGrid3D(Grid&, string);

    protected:
        void CellData(void);
        void EdgeData(void);
        void BoundaryData(void);

        Config *config;
        unsigned short nDim;
        vector<Node> nodes;
        vector<unique_ptr<Cell>> cells;
        vector<Edge> edges;
        vector<unique_ptr<Boundary>> bndrys;
};
