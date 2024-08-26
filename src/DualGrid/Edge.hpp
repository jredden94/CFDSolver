#pragma once

#include <vector>
#include "../DualGrid/Node.hpp"

using namespace std;

class Edge {
    public:
        Edge();
        ~Edge();

        const unsigned long Node1(void) const;
        const unsigned long Node2(void) const;
        const vector<unsigned long> Cells(void) const;
        const vector<double>& EdgeVector(void) const;
        const vector<double>& AreaVector(void) const;
        const vector<double>& Midpoint(void) const;
        const double Area(void) const;
        const double Length(void) const;
        void ComputeMetrics(const Node &n1, const Node &n2);
        void AddDirArea(const vector<double> &dirArea);
        void NormalizeDirAreaVec(void);

        friend class Grid;

    private:
        unsigned long iNode1, iNode2;
        vector<unsigned long> iCells;

        double length, area_mag;
        vector<double> mid_vec;
        vector<double> edge_vec;
        vector<double> area_vec;
};
