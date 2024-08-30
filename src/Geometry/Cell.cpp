#include "Cell.hpp"

string Cell::CellTypeName(const Cell::Type &type) {
    switch (type) {
        case (Cell::Type::TETRA_4) : return "TETRA_4"; break;
        case (Cell::Type::HEXA_8) : return "HEXA_8"; break;
        case (Cell::Type::PYRA_5) : return "PYRA_5"; break;
        case (Cell::Type::PENTA_6) : return "PENTA_6"; break;
        case (Cell::Type::TRI_3) : return "TRI_3"; break;
        case (Cell::Type::QUAD_4) : return "QUAD_4"; break;
        case (Cell::Type::BAR_2) : return "BAR_2"; break;
        default : return ""; 
    }
}

Cell::Cell() {};
Cell::Cell(const Type &type, vector<unsigned long> &connectivity, vector<Node*>& pNodes) 
: type(type) { 
    this->iNodes.swap(connectivity); 
    this->pNodes.swap(pNodes);
}

Cell::~Cell() {};

const vector<unsigned long>& Cell::GetConnectivity() const {
    return iNodes;
}

const Cell::Type& Cell::GetType() const {
    return type;
}

const size_t Cell::GetNodeCount() const {
    return iNodes.size();
}

void Cell::AddEdge(unsigned long iEdge, Edge *pEdge) {
    bool hasEdge = false;
    for (const unsigned long &e : iEdges) {
        if (e == iEdge) {
            hasEdge = true;
            break;
        }
    }

    if (!hasEdge) {
        iEdges.push_back(iEdge);
        pEdges.push_back(pEdge);
    }
}

unique_ptr<Cell> Cell::CreateCell(const Type &type, 
        vector<unsigned long> &connectivity, vector<Node*> &nodes) {
    switch (type) {
        case(Type::BAR_2) : return make_unique<Bar>(connectivity, nodes); break;
        case(Type::TRI_3) : return make_unique<Tri>(connectivity, nodes); break;
        case(Type::QUAD_4) : return make_unique<Quad>(connectivity, nodes); break;
        case(Type::TETRA_4) : return make_unique<Tetra>(connectivity, nodes); break;
        default : return nullptr;
    }
}

unique_ptr<Cell> Cell::CreateCell(const int vtk, 
        vector<unsigned long> &connectivity, vector<Node*> &nodes) {
    switch (vtk) {
        case(3) : return make_unique<Bar>(connectivity, nodes); break;
        case(5) : return make_unique<Tri>(connectivity, nodes); break;
        case(9) : return make_unique<Quad>(connectivity, nodes); break;
        case(10) : return make_unique<Tetra>(connectivity, nodes); break;
        default : return nullptr;
    }
}

size_t Cell::FindEdges(void) const {
    size_t nEdges=0;
    size_t nNodes = iNodes.size();
    for (size_t i = 0; i < nNodes; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        unsigned long n2Ind = i == nNodes-1 ? 0 : i+1;
        Node *n2 = pNodes[n2Ind]; // = i == nNodes-1 ? nodePtr[0] : nodePtr[i+1];
        unsigned long n2_I = iNodes[n2Ind]; // = i == nNodes-1 ? nodeI[0] : nodeI[i+1];


        const vector<unsigned long> &n1Cells = n1->Cell();
        const vector<unsigned long> &n2Cells = n2->Cell();

        size_t nSharedCells = 0;

        for (const auto ii : n1Cells) {
            for (const auto jj : n2Cells) {
                if (ii == jj) nSharedCells++;
            }
        }

        if (nSharedCells == 1 || n2_I > n1_I) nEdges++;
    }

    return nEdges;
}

//normal1 = triangle_area_vector(xm,xc,xcl, ym,yc,ycl, zm,zc,zcl)
vector<zdouble> Cell::TriAreaVector(const vector<zdouble> &p1, 
        const vector<zdouble> &p2, const vector<zdouble> &p3, const zdouble &coeff) const {
    vector<zdouble> aVec(3,0);
    aVec[0] = coeff * 0.5 * ((p1[1]-p3[1]) * (p2[2]-p3[2]) - (p1[2]-p3[2]) * (p2[1]-p3[1]));
    aVec[1] = coeff * 0.5 * ((p1[2]-p3[2]) * (p2[0]-p3[0]) - (p1[0]-p3[0]) * (p2[2]-p3[2]));
    aVec[2] = coeff * 0.5 * ((p1[0]-p3[0]) * (p2[1]-p3[1]) - (p1[1]-p3[1]) * (p2[0]-p3[0]));
    return aVec;
}
vector<zdouble> Cell::TriCentroid(const Node* n1, 
        const Node* n2, const Node* n3) const {
    vector<zdouble> ctr(3,0);
    ctr[0] = (n1->X() + n2->X() + n3->X() ) / 3;
    ctr[1] = (n1->Y() + n2->Y() + n3->Y() ) / 3;
    ctr[2] = (n1->Z() + n2->Z() + n3->Z() ) / 3;

    return ctr;
}
