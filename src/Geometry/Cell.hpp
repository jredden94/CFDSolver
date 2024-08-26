#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <memory>
#include "../DualGrid/Node.hpp"
#include "../DualGrid/Edge.hpp"

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
        virtual double ComputeVolume(void) const = 0;
        virtual vector<double> AreaVector(void) const = 0;
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
        vector<double> TriAreaVector(const vector<double>&, const vector<double>&, const vector<double>&, const double &coeff=1) const;
        vector<double> TriCentroid(const Node* n1, const Node* n2, const Node* n3) const;

        Type type;
        vector<unsigned long> iNodes;
        vector<unsigned long> iNbrs;
        vector<unsigned long> iEdges;
        vector<Node*> pNodes;
        vector<Edge*> pEdges;

        double cx, cy, cz; // centroid
};

class Bar : public Cell {
    public:
        Bar() = delete;
        Bar(const Bar&) = delete;
        Bar(vector<unsigned long>& connectivity, vector<Node*>& nodePtr)
            : Cell(Type::BAR_2, connectivity, nodePtr) { };
        ~Bar();

        double ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<double> AreaVector(void) const override;
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

        double ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<double> AreaVector(void) const override;
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

        double ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<double> AreaVector(void) const override;
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

        double ComputeVolume(void) const override;
        vector<vector<size_t>> GetEdges(void) const override;
        vector<double> AreaVector(void) const override;
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
