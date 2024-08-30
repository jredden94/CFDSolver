#include <iostream>
#include <cgnslib.h>
#include <cstring>
#include "DualGrid/Node.hpp"
#include "IO/GridReader.hpp"
#include "Common/Config.hpp"
#include "Numerics/Viscous/Viscous.hpp"
#include "Solver/Solver.hpp"
#include "LinAlg/BLAS.hpp"
#include "IO/VTK.hpp"

using namespace std;

void ReadGrid2D(Grid&grid, string filename) {
    std::ifstream file(filename);

    vector<Node> nodes;
    vector<unique_ptr<Cell>> cells;
    size_t nNodes, nTri, nQuad;
    file >> nNodes >> nTri >> nQuad;

    cout << "nNodes: " << nNodes << "\tnTri: " << nTri << "\tnQuad: " << nQuad << endl;

    for (size_t i = 0; i < nNodes; i++) {
        zdouble x = 0, y = 0, z = 0;
        file >> x >> y;
        nodes.push_back(Node(x, y, z));
    };

    for (size_t i = 0; i < nTri; i++) {
        unsigned long n1, n2, n3;
        file >> n1 >> n2 >> n3;
        vector<unsigned long> conn {n1-1, n2-1, n3-1};
        vector<Node*> cNodes {&nodes[n1-1], &nodes[n2-1], &nodes[n3-1]};
        cells.push_back(Cell::CreateCell(Cell::Type::TRI_3, conn, cNodes));
    };
    for (size_t i = 0; i < nQuad; i++) {
        unsigned long n1, n2, n3, n4;
        file >> n1 >> n2 >> n3 >> n4;
        vector<unsigned long> conn {n1-1, n2-1, n3-1, n4-1};
        vector<Node*> cNodes {&nodes[n1-1], &nodes[n2-1], &nodes[n3-1], &nodes[n4-1]};
        cells.push_back(Cell::CreateCell(Cell::Type::QUAD_4, conn, cNodes));
    };
    for (size_t i = 0; i < nQuad; i++) {
        unsigned long n1, n2, n3, n4;
        file >> n1 >> n2 >> n3 >> n4;
        vector<unsigned long> conn {n1-1, n2-1, n3-1, n4-1};
        vector<Node*> cNodes {&nodes[n1-1], &nodes[n2-1], &nodes[n3-1], &nodes[n4-1]};
        cells.push_back(Cell::CreateCell(Cell::Type::QUAD_4, conn, cNodes));
    };


    vector<vector<unsigned long>> bnodes;
    vector<Boundary> bndrys;

    size_t nbndry;
    file >> nbndry;
    bnodes.resize(nbndry);
    for (size_t i = 0; i < nbndry; i++) {
        size_t nbnodes;
        file >> nbnodes;

        bnodes[i].resize(nbnodes);
    }

    for (size_t i = 0; i < nbndry; i++) {
        size_t nbnodes = bnodes[i].size();
        for (size_t j = 0; j < nbnodes; j++) {
            int n1, n2;
            file >> n1;// >> n2;
            bnodes[i][j] = n1 - 1;
        }
    }

    vector<vector<unique_ptr<Cell>>> bCells(nbndry);

    for (size_t i = 0; i < nbndry; i++) {
        vector<unique_ptr<Cell>> &bc = bCells[i];
        const vector<unsigned long> &conn = bnodes[i];
        size_t nElms = conn.size() - 1;

        size_t counter = 0;
        for (size_t j = 0; j < nElms; j++) {
            vector<unsigned long> line(2,0);
            vector<Node*> bnptrs;
            line[0] = conn[counter++];
            line[1] = conn[counter];
            
            bnptrs.push_back(&nodes[line[0]]);
            bnptrs.push_back(&nodes[line[1]]);
            bc.push_back(Cell::CreateCell(Cell::Type::BAR_2, line, bnptrs));
        }
    }

    //vector<Boundary> bounds;
    vector<unique_ptr<Boundary>> bounds;
    string name {"name"};
    for (vector<unique_ptr<Cell>> &bc : bCells) {
        bounds.push_back(make_unique<Boundary>(name, bc));
    }

    //This line is fine
    Boundary b{name, bCells[0]};

    cout << "NodeCount: " << nodes.size() << endl;
    cout << "CellCount: " << cells.size() << endl;

    grid.nodes.swap(nodes);
    grid.cells.swap(cells);
    grid.bndrys.swap(bounds);
}

void ReadGrid3D(Grid &grid, string filename) {
    std::ifstream file(filename);

    vector<Node> nodes;
    vector<unique_ptr<Cell>> cells;
    size_t nNodes, nTri, nQuad, nTet, nHex, nPyr, nPrism;
    file >> nNodes >> nTri >> nQuad >> nTet >> nPyr >> nPrism >> nHex;

    cout << "nNodes: " << nNodes << "\tnTri: " << nTri << "\tnTet: " << nTet << endl;

    for (size_t i = 0; i < nNodes; i++) {
        zdouble x = 0, y = 0, z = 0;
        file >> x >> y >> z;
        nodes.push_back(Node(x, y, z));
    }

    vector<unique_ptr<Cell>> bCells;

    if (nTri > 0) {
        for (size_t i = 0; i < nTri; i++) {
            unsigned long n1, n2, n3;
            file >> n1 >> n2 >> n3;
            vector<unsigned long> conn {n1-1, n2-1, n3-1};
            vector<Node*> cNodes {&nodes[n1-1], &nodes[n2-1], &nodes[n3-1]};
            bCells.push_back(Cell::CreateCell(Cell::Type::TRI_3, conn, cNodes));
        }
    }
    if (nQuad > 0) {
        for (size_t i = 0; i < nQuad; i++) {
            unsigned long n1, n2, n3, n4;
            file >> n1 >> n2 >> n3 >> n4;
            vector<unsigned long> conn {n1-1, n2-1, n3-1, n4-1};
            vector<Node*> cNodes {&nodes[n1-1], &nodes[n2-1], &nodes[n3-1], &nodes[n4-1]};
            bCells.push_back(Cell::CreateCell(Cell::Type::QUAD_4, conn, cNodes));
        }
    }

    vector<unsigned long> bcMap;
    for (size_t i = 0; i < nTri + nQuad; i++) {
        unsigned long n1;
        file >> n1;
        bcMap.push_back(n1);
    }

    if (nTet > 0) {
        for (size_t i = 0; i < nTet; i++) {
            unsigned long n1, n2, n3, n4;
            file >> n1 >> n2 >> n3 >> n4;
            vector<unsigned long> conn {n1-1, n2-1, n3-1, n4-1};
            vector<Node*> cNodes {&nodes[n1-1], &nodes[n2-1], &nodes[n3-1], &nodes[n4-1]};
            cells.push_back(Cell::CreateCell(Cell::Type::TETRA_4, conn, cNodes));
        }
    }

    size_t nBound = 3;
    vector<vector<unique_ptr<Cell>>> bbcells(6);

    for (size_t i = 0; i < bCells.size(); i++) {
        bbcells[bcMap[i]-1].push_back(std::move(bCells[i]));
    }

    vector<unique_ptr<Boundary>> bounds;;
    for (vector<unique_ptr<Cell>> &bbcell : bbcells) {
        //unique_ptr<Boundary> bndry = make_unique<Boundary>("f", bbcell);
        bounds.push_back(make_unique<Boundary>("f", bbcell));
    }

    grid.nodes.swap(nodes);
    grid.cells.swap(cells);
    grid.bndrys.swap(bounds);
}

/* TODO - Check geometry data, make sure everything is 2D/3D ready. */
int main(void) {
    Config &config = Config::GetConfig();
    config.SetDefault();
    Grid grid;
    unique_ptr<GridReader> reader = GridReader::GetReader();
    reader->ReadFile();
    reader->TransferToGrid(grid);

    grid.BuildDualGrid();
    //grid.FlipBoundaryNorms();

    config.SetNumEqn(grid.Nodes().size() );

    unique_ptr<Solver> solver = Solver::CreateSolver(&grid);
    solver->Solve();

    VTK::WriteVTK("Testing.vtu", grid, *solver);

    return 0;
}
