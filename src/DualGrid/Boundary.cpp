#include "Boundary.hpp"

string Boundary::BCTypeName(const Boundary::BCType bcType) {
    switch (bcType) {
        case(Boundary::BCType::Farfield) : return "Farfield"; break;
        case(Boundary::BCType::WallInviscid) : return "InviscidWall"; break;
        default : return "Empty"; break;
    }
}

Boundary::Boundary() {}
Boundary::Boundary(string name, vector<unique_ptr<Cell>> &cells) 
    : name(name) {
        this->cells.swap(cells);
    }
Boundary::~Boundary() {}

const vector<unique_ptr<Cell>>& Boundary::Cells(void) const { return cells; }
const Boundary::BCType& Boundary::Type(void) const { return type; }
const Cell::Type& Boundary::CellType(void) const { return cellType; }
const string& Boundary::Name(void) const { return name; }

const vector<zdouble>& Boundary::NormX(void) const { return norm_x; }
const vector<zdouble>& Boundary::NormY(void) const { return norm_y; }
const vector<zdouble>& Boundary::Area(void) const { return faceArea; }
const vector<vector<zdouble>>& Boundary::NodeNorms(void) const { return nodeNorms; }
const vector<zdouble>& Boundary::DualAreas(void) const { return nodeAreas; }
const vector<unsigned long>& Boundary::Nodes(void) const { return iNodes; }

void Boundary::SetType(Boundary::BCType type) {
    this->type = type;
}
