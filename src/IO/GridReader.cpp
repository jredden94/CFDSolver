#include "GridReader.hpp"
#include "SU2Reader.hpp"

GridReader::GridReader() { 
    config = &Config::GetConfig();
    filename = config->GetGridFilename();
}
GridReader::~GridReader() { }

unique_ptr<GridReader> GridReader::GetReader() { 
    Config *config = &Config::GetConfig();
    string filename = config->GetGridFilename();
    size_t ext_pos = filename.find('.');
    string extension = filename.substr(ext_pos);

    if (extension == ".su2") return make_unique<SU2Reader>();
    else return nullptr;
}

void GridReader::ReadFile(void) { }
void GridReader::TransferToGrid(Grid &grid) {
    nodes.swap(grid.nodes);
    cells.swap(grid.cells);
    boundaries.swap(grid.bndrys);
}
