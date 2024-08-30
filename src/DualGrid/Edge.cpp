#include "Edge.hpp"

Edge::Edge() : mid_vec(3,0), edge_vec(3,0), area_vec(3,0)  { }
Edge::~Edge() { }

const unsigned long Edge::Node1() const { return iNode1; }
const unsigned long Edge::Node2() const { return iNode2; }
const vector<unsigned long> Edge::Cells() const { return iCells; }
const vector<zdouble>& Edge::EdgeVector() const { return edge_vec; } 
const vector<zdouble>& Edge::AreaVector() const { return area_vec; }
const vector<zdouble>& Edge::Midpoint() const { return mid_vec; }
const zdouble Edge::Area() const { return area_mag; }
const zdouble Edge::Length() const { return length; }

void Edge::ComputeMetrics(const Node &n1, const Node &n2) {
    mid_vec[0] = 0.5 * (n1.X() + n2.X() );
    mid_vec[1] = 0.5 * (n1.Y() + n2.Y() );
    mid_vec[2] = 0.5 * (n1.Z() + n2.Z() );

    edge_vec[0] = n2.X() - n1.X();
    edge_vec[1] = n2.Y() - n1.Y();
    edge_vec[2] = n2.Z() - n1.Z();
    length = std::sqrt(edge_vec[0] * edge_vec[0] + 
            edge_vec[1] * edge_vec[1] + edge_vec[2] * edge_vec[2]);
    edge_vec[0] = edge_vec[0] / length;
    edge_vec[1] = edge_vec[1] / length;
    edge_vec[2] = edge_vec[2] / length;

    // Edge Cells
    const auto &n1Cells = n1.Cell();
    const auto &n2Cells = n2.Cell();
    for (size_t ii = 0; ii < n1Cells.size(); ii++) {
        const unsigned long &n1c = n1Cells[ii];
        for (size_t jj = 0; jj < n2Cells.size(); jj++) {
            const unsigned long &n2c = n2Cells[jj];
            if (n1c == n2c) iCells.push_back(n1c);
        }
    }
}

void Edge::AddDirArea(const vector<zdouble> &dirArea) {
    area_vec[0] += dirArea[0];
    area_vec[1] += dirArea[1];
    area_vec[2] += dirArea[2];
}

void Edge::NormalizeDirAreaVec(void) {
    area_mag = sqrt(area_vec[0] * area_vec[0] 
            + area_vec[1] * area_vec[1] + area_vec[2] * area_vec[2]);
    area_vec[0] /= area_mag;
    area_vec[1] /= area_mag;
    area_vec[2] /= area_mag;
}
