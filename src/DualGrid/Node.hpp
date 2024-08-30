#pragma once

#include <vector>
#include <cmath>
#include "../Common/AD.hpp"

using namespace std;

class Node {
    public:
        Node();
        Node(zdouble x, zdouble y, zdouble z);
        Node(const Node &node);
        ~Node();

        void operator=(const Node &node) = delete;

        const zdouble& X(void) const;
        const zdouble& Y(void) const;
        const zdouble& Z(void) const;
        void SetX(zdouble x);
        void SetY(zdouble y);
        void SetZ(zdouble z);
        const vector<zdouble> Coordinates(void) const;
        const vector<unsigned long> Cell(void) const;
        const zdouble& DualVolume(void) const;
        const vector<unsigned long>& Neighbors(void) const;
        bool HasNeighbor(unsigned long iNbr) const;
        void AddNeighbor(const unsigned long&);

        friend class Grid;

    private:
        zdouble x;
        zdouble y;
        zdouble z;
        zdouble dualVol;
        vector<unsigned long> iCells;
        vector<unsigned long> iNbrs;
};
