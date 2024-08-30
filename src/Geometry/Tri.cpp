#include "Cell.hpp"

Tri::~Tri() { }

zdouble Tri::ComputeVolume() const {
    const zdouble &x1 = pNodes[0]->X();
    const zdouble &x2 = pNodes[1]->X();
    const zdouble &x3 = pNodes[2]->X();

    const zdouble &y1 = pNodes[0]->Y();
    const zdouble &y2 = pNodes[1]->Y();
    const zdouble &y3 = pNodes[2]->Y();

    zdouble volume = 0.5 * ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
    volume *= 0.5;

    return volume;
}

void Tri::ComputeCentroid() {
    cx = cy = cz = 0;
    for (const auto *n : pNodes) {
        cx += n->X();
        cy += n->Y();
    }

    cx /= 3;
    cy /= 3;
}

vector<vector<size_t>> Tri::GetEdges() const {
    vector<vector<size_t>> edgeI;
    size_t nNodes = iNodes.size();
    for (size_t i = 0; i < nNodes; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        unsigned long n2Ind = i == nNodes-1 ? 0 : i+1;
        Node *n2 = pNodes[n2Ind]; // = i == nNodes-1 ? pNodes[0] : pNodes[i+1];
        unsigned long n2_I = iNodes[n2Ind]; // = i == nNodes-1 ? iNodes[0] : iNodes[i+1];


        const vector<unsigned long> &n1Cells = n1->Cell();
        const vector<unsigned long> &n2Cells = n2->Cell();

        size_t nSharedCells = 0;

        for (const auto ii : n1Cells) {
            for (const auto jj : n2Cells) {
                if (ii == jj) nSharedCells++;
            }
        }

        if (nSharedCells == 1 || n2_I > n1_I) edgeI.emplace_back(vector<size_t>{n1_I, n2_I});
    }

    return edgeI;
}

vector<zdouble> Tri::AreaVector() const {
    Node *n1 = pNodes[0];
    Node *n2 = pNodes[1];
    Node *n3 = pNodes[2];
    return TriAreaVector(n1->Coordinates(), n2->Coordinates(), n3->Coordinates(), -1);
}

void Tri::FindNeighbors(unsigned long iCell) {
    iNbrs.clear();
    size_t nNodes = iNodes.size();
    for (size_t j = 0; j < nNodes; j++) {
        const unsigned long n1_i = j;
        const unsigned long n2_i = j != nNodes - 1 ? j + 1 : 0;
        const vector<unsigned long> &n1Cells = pNodes[n1_i]->Cell();
        const vector<unsigned long> &n2Cells = pNodes[n2_i]->Cell();

        for (size_t ii = 0; ii < n1Cells.size(); ii++) {
            for (size_t jj = 0; jj < n2Cells.size(); jj++) {
                if (n1Cells[ii] == n2Cells[jj] && n1Cells[ii] != iCell) 
                    iNbrs.push_back(n1Cells[ii]);
            }
        }
    }
}

void Tri::ComputeDirectedArea(void) {
    for (Edge *e : pEdges) {
        const vector<zdouble> &midpoint = e->Midpoint();
        const vector<zdouble> &edgevec = e->EdgeVector();
        vector<zdouble> dirArea {0, 0, 0};
        dirArea[0] = -(midpoint[1] - cy);
        dirArea[1] = midpoint[0] - cx;
        zdouble dot = dirArea[0] * edgevec[0] + dirArea[1] * edgevec[1];
        if (dot < 0) {
            dirArea[0] *= -1;
            dirArea[1] *= -1;
        }

        e->AddDirArea(dirArea);
    }
}

void Tri::NodeNeighbors(void) {
    for (size_t i = 0; i < 3; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        unsigned long n2_ind = i == 2 ? 0 : i+1;
        Node *n2 = pNodes[n2_ind];

        unsigned long n2_I = iNodes[n2_ind];
        if (!n1->HasNeighbor(n2_I)) n1->AddNeighbor(n2_I);
        if (!n2->HasNeighbor(n1_I))n2->AddNeighbor(n1_I);
    }
}
