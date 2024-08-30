#include "Node.hpp"

Node::Node() : x(0), y(0), z(0) { }
Node::Node(zdouble x, zdouble y, zdouble z) : x(x), y(y), z(z), dualVol(0){ }
Node::Node(const Node &node) : x(node.x), y(node.y), z(node.z), dualVol(0) { }
Node::~Node() { }

const zdouble& Node::X() const { return x; }
const zdouble& Node::Y() const { return y; }
const zdouble& Node::Z() const { return z; }
void Node::SetX(zdouble x) { this->x = x; } 
void Node::SetY(zdouble y) { this->y = y; } 
void Node::SetZ(zdouble z) { this->z = z; } 
const vector<zdouble> Node::Coordinates() const { return vector<zdouble>{x, y, z}; }
const zdouble& Node::DualVolume() const { return dualVol; }
const vector<unsigned long> Node::Cell() const { return iCells; }
const vector<unsigned long>& Node::Neighbors() const { return iNbrs; }

void Node::AddNeighbor(const unsigned long &iNode) { iNbrs.push_back(iNode); }

bool Node::HasNeighbor(unsigned long iNbr) const {
    for (const auto &nbr : iNbrs) {
        if (nbr == iNbr) return true;
    }
    return false;
}
