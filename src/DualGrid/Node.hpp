#pragma once

#include <vector>
#include <cmath>

using namespace std;

class Node {
    public:
        Node();
        Node(double x, double y, double z);
        Node(const Node &node);
        ~Node();

        void operator=(const Node &node) = delete;

        const double& X(void) const;
        const double& Y(void) const;
        const double& Z(void) const;
        void SetX(double x);
        void SetY(double y);
        void SetZ(double z);
        const vector<double> Coordinates(void) const;
        const vector<unsigned long> Cell(void) const;
        const double& DualVolume(void) const;
        const vector<unsigned long>& Neighbors(void) const;
        bool HasNeighbor(unsigned long iNbr) const;
        void AddNeighbor(const unsigned long&);

        friend class Grid;

    private:
        double x;
        double y;
        double z;
        double dualVol;
        vector<unsigned long> iCells;
        vector<unsigned long> iNbrs;
};
