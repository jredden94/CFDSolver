#pragma once

#include <fstream>
#include <sstream>
#include <istream>
#include "../DualGrid/Grid.hpp"
#include "../Common/Config.hpp"
#include "../Common/AD.hpp"

class GridReader {
    public:
        GridReader();
        virtual ~GridReader();

        static unique_ptr<GridReader> GetReader(void);

        void TransferToGrid(Grid &grid);
        virtual void ReadFile(void) = 0;

    protected:
        Config *config;
        string filename;

        unsigned short nDim, nVar;
        unsigned long nNodes, nCells, nBoundaries;
        vector<Node> nodes;
        vector<unique_ptr<Cell>> cells;
        vector<unique_ptr<Boundary>> boundaries;
};
