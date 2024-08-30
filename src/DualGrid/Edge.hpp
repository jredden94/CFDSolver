#pragma once

#include <vector>
#include "../DualGrid/Node.hpp"
#include "../Common/AD.hpp"

using namespace std;

class Edge {
    public:
        Edge();
        ~Edge();

        const unsigned long Node1(void) const;
        const unsigned long Node2(void) const;
        const vector<unsigned long> Cells(void) const;
        const vector<zdouble>& EdgeVector(void) const;
        const vector<zdouble>& AreaVector(void) const;
        const vector<zdouble>& Midpoint(void) const;
        const zdouble Area(void) const;
        const zdouble Length(void) const;
        void ComputeMetrics(const Node &n1, const Node &n2);
        void AddDirArea(const vector<zdouble> &dirArea);
        void NormalizeDirAreaVec(void);

        friend class Grid;

    private:
        unsigned long iNode1, iNode2;
        vector<unsigned long> iCells;

        zdouble length, area_mag;
        vector<zdouble> mid_vec;
        vector<zdouble> edge_vec;
        vector<zdouble> area_vec;
};
