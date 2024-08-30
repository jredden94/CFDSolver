#pragma once

#include <vector>
#include <string>
#include "../Geometry/Cell.hpp"
#include "../Common/AD.hpp"

using namespace std;

class Boundary {
    public:
        enum class BCType{ Farfield, WallInviscid, WallViscous, SupersonicInlet, SupersonicOutlet, Inlet, Outlet, Empty };
        static string BCTypeName(const Boundary::BCType);

        Boundary();
        Boundary(string, vector<unique_ptr<Cell>>&);
        ~Boundary();

        const BCType& Type(void) const;
        const Cell::Type& CellType(void) const;
        const vector<unique_ptr<Cell>>& Cells(void) const;
        const string& Name(void) const;

        const vector<zdouble>& NormX(void) const;
        const vector<zdouble>& NormY(void) const;
        const vector<zdouble>& Area(void) const;

        const vector<vector<zdouble>>& NodeNorms(void) const;
        const vector<zdouble>& DualAreas(void) const;
        const vector<unsigned long>& Nodes(void) const;
        void SetType(Boundary::BCType);

        friend class Grid;

    private:
        string name;
        vector<unique_ptr<Cell>> cells;
        vector<unsigned long> iNodes;
        BCType type;
        Cell::Type cellType;

        // Face normals. Make these cell-specific later
        vector<zdouble> faceArea;
        vector<zdouble> norm_x;
        vector<zdouble> norm_y;

        // Node normals and dual areas
        vector<vector<zdouble>> nodeNorms;
        vector<zdouble> nodeAreas;
};
