#include "Node.hpp"

Node::Node() : x(0), y(0), z(0) { }
Node::Node(double x, double y, double z) : x(x), y(y), z(z), dualVol(0){ }
Node::Node(const Node &node) : x(node.x), y(node.y), z(node.z), dualVol(0) { }
Node::~Node() { }

const double& Node::X() const { return x; }
const double& Node::Y() const { return y; }
const double& Node::Z() const { return z; }
void Node::SetX(double x) { this->x = x; } 
void Node::SetY(double y) { this->y = y; } 
void Node::SetZ(double z) { this->z = z; } 
const vector<double> Node::Coordinates() const { return vector<double>{x, y, z}; }
const double& Node::DualVolume() const { return dualVol; }
const vector<unsigned long> Node::Cell() const { return iCells; }
const vector<unsigned long>& Node::Neighbors() const { return iNbrs; }

void Node::AddNeighbor(const unsigned long &iNode) { iNbrs.push_back(iNode); }

bool Node::HasNeighbor(unsigned long iNbr) const {
    for (const auto &nbr : iNbrs) {
        if (nbr == iNbr) return true;
    }
    return false;
}
