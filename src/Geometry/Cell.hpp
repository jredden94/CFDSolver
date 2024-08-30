#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <memory>
#include "../DualGrid/Node.hpp"
#include "../DualGrid/Edge.hpp"
#include "../Common/AD.hpp"

using namespace std;

class Cell {
    public:
        enum class Type { BAR_2, TRI_3, QUAD_4, TETRA_4, PYRA_5, PENTA_6, HEXA_8};
        static string CellTypeName(const Cell::Type &type);

        Cell();
        Cell(const Cell &cell) = delete;
        Cell(const Type& type, vector<unsigned long>& connectivity, vector<Node*>& nodePtr);
        virtual ~Cell();

        Cell& operator=(const Cell &cell) =delete;;

        const vector<unsigned long>& GetConnectivity(void) const;
        const Type& GetType(void) const;
        const size_t GetNodeCount(void) const;
        void AddEdge(unsigned long iEdge, Edge *pEdge);

        virtual vector<vector<size_t>> GetEdges(void) const = 0;
        virtual zdouble ComputeVolume(void) const = 0;
        virtual vector<zdouble> AreaVector(void) const = 0;
        virtual void ComputeCentroid(void) = 0;
        virtual void ComputeDirectedArea(void) = 0;
        virtual void FindNeighbors(unsigned long iCell) = 0;
        virtual void NodeNeighbors(void) = 0;

        size_t FindEdges(void) const;

        static unique_ptr<Cell> CreateCell(const Type &type, 
                vector<unsigned long> &connectivity, vector<Node*> &nodes);
        static unique_ptr<Cell> CreateCell(const int vtk, 
                vector<unsigned long> &connectivity, vector<Node*> &nodes);

        friend class Grid;

    protected:
        vector<zdouble> TriAreaVector(const vector<zdouble>&, const vector<zdouble>&, const vector<zdouble>&, const zdouble &coeff=1) const;
        vector<zdouble> TriCentroid(const Node* n1, const Node* n2, const Node* n3) const;

        Type type;
        vector<unsigned long> iNodes;
        vector<unsigned long> iNbrs;
        vector<unsigned long> iEdges;
        vector<Node*> pNodes;
        vector<Edge*> pEdges;

        zdouble cx, cy, cz; // centroid
};

class Bar : public Cell {
    public:
        Bar() = delete;
        Bar(const Bar&) = delete;
        Bar(vector<unsigned long>& connectivity, vector<Node*>& nodePtr)
            : Cell(Type::BAR_2, connectivity, nodePtr) { };
        ~Bar();

        zdouble ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<zdouble> AreaVector(void) const override;
        void ComputeCentroid(void) override;
        void ComputeDirectedArea(void) override;
        void FindNeighbors(unsigned long iCell) override;
        void NodeNeighbors(void) override;
};

class Tri : public Cell {
    public:
        Tri() = delete;
        Tri(const Tri&) = delete;
        Tri(vector<unsigned long>& connectivity, vector<Node*>& nodePtr)
            : Cell(Type::TRI_3, connectivity, nodePtr) { };
        ~Tri();

        zdouble ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<zdouble> AreaVector(void) const override;
        void ComputeCentroid(void) override;
        void FindNeighbors(unsigned long iCell) override;
        void ComputeDirectedArea(void) override;
        void NodeNeighbors(void) override;
};

class Quad : public Cell {
    public:
        Quad() = delete;
        Quad(const Quad&) = delete;
        Quad(vector<unsigned long>& connectivity, vector<Node*>& nodePtr)
            : Cell(Type::QUAD_4, connectivity, nodePtr) { };
        ~Quad();

        zdouble ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<zdouble> AreaVector(void) const override;
        void ComputeCentroid(void) override;
        void FindNeighbors(unsigned long iCell) override;
        void ComputeDirectedArea(void) override;
        void NodeNeighbors(void) override;
};

class Tetra : public Cell {
    public:
        Tetra() = delete;
        Tetra(const Tetra&) = delete;
        Tetra(vector<unsigned long>& connectivity, vector<Node*>& nodePtr)
            : Cell(Type::TETRA_4, connectivity, nodePtr) { };
        ~Tetra();

        zdouble ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<zdouble> AreaVector(void) const override;
        void ComputeCentroid(void) override;
        void FindNeighbors(unsigned long iCell) override;
        void ComputeDirectedArea(void) override;
        void NodeNeighbors(void) override;
};

/*
   class Pyra : public Cell {

   private:
   void CalculateVolume(void) override;
   };

   class Penta : public Cell {

   private:
   void CalculateVolume(void) override;
   };

   class Hexa : public Cell {

   private:
   void CalculateVolume(void) override;
   };
   */
