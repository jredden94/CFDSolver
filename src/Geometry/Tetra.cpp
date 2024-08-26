#include "Cell.hpp"

Tetra::~Tetra() { }

//  !x-component of the area vector
//   area_vector(ix) = 0.5_dp*( (y1-y3)*(z2-z3)-(z1-z3)*(y2-y3) )
//
//  !y-component of the area vector
//   area_vector(iy) = 0.5_dp*( (z1-z3)*(x2-x3)-(x1-x3)*(z2-z3) )
//
//  !z-component of the area vector
//   area_vector(iz) = 0.5_dp*( (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3) )

double Tetra::ComputeVolume() const {
    double volume = 0;
    double xc, yc, zc;
    vector<double> aVec(3,0);
    Node *n1, *n2, *n3;
    // Triangle 1,2,3
    n1 = pNodes[0];
    n2 = pNodes[1];
    n3 = pNodes[2];
    xc = (n1->X() + n2->X() + n3->X()) / 3;
    yc = (n1->Y() + n2->Y() + n3->Y()) / 3;
    zc = (n1->Z() + n2->Z() + n3->Z()) / 3;
    aVec[0] = 0.5 * ( (n1->Y() - n3->Y() ) * (n2->Z() - n3->Z() ) 
            - (n1->Z() - n3->Z() ) * (n2->Y() - n3->Y() ));
    aVec[1] = 0.5 * ( (n1->Z() - n3->Z() ) * (n2->X() - n3->X() ) 
            - (n1->X() - n3->X() ) * (n2->Z() - n3->Z() ));
    aVec[2] = 0.5 * ( (n1->X() - n3->X() ) * (n2->Y() - n3->Y() ) 
            - (n1->Y() - n3->Y() ) * (n2->X() - n3->X() ));
    volume += xc * aVec[0] + yc * aVec[1] + zc * aVec[2];

    // Triangle 1,4,3
    n1 = pNodes[0];
    n2 = pNodes[3];
    n3 = pNodes[2];
    xc = (n1->X() + n2->X() + n3->X()) / 3;
    yc = (n1->Y() + n2->Y() + n3->Y()) / 3;
    zc = (n1->Z() + n2->Z() + n3->Z()) / 3;
    aVec[0] = 0.5 * ( (n1->Y() - n3->Y() ) * (n2->Z() - n3->Z() ) 
            - (n1->Z() - n3->Z() ) * (n2->Y() - n3->Y() ));
    aVec[1] = 0.5 * ( (n1->Z() - n3->Z() ) * (n2->X() - n3->X() ) 
            - (n1->X() - n3->X() ) * (n2->Z() - n3->Z() ));
    aVec[2] = 0.5 * ( (n1->X() - n3->X() ) * (n2->Y() - n3->Y() ) 
            - (n1->Y() - n3->Y() ) * (n2->X() - n3->X() ));
    volume += xc * aVec[0] + yc * aVec[1] + zc * aVec[2];

    // Triangle 1,2,4
    n1 = pNodes[0];
    n2 = pNodes[1];
    n3 = pNodes[3];
    xc = (n1->X() + n2->X() + n3->X()) / 3;
    yc = (n1->Y() + n2->Y() + n3->Y()) / 3;
    zc = (n1->Z() + n2->Z() + n3->Z()) / 3;
    aVec[0] = 0.5 * ( (n1->Y() - n3->Y() ) * (n2->Z() - n3->Z() ) 
            - (n1->Z() - n3->Z() ) * (n2->Y() - n3->Y() ));
    aVec[1] = 0.5 * ( (n1->Z() - n3->Z() ) * (n2->X() - n3->X() ) 
            - (n1->X() - n3->X() ) * (n2->Z() - n3->Z() ));
    aVec[2] = 0.5 * ( (n1->X() - n3->X() ) * (n2->Y() - n3->Y() ) 
            - (n1->Y() - n3->Y() ) * (n2->X() - n3->X() ));
    volume += xc * aVec[0] + yc * aVec[1] + zc * aVec[2];

    // Triangle 2,3,4
    n1 = pNodes[1];
    n2 = pNodes[2];
    n3 = pNodes[3];
    xc = (n1->X() + n2->X() + n3->X()) / 3;
    yc = (n1->Y() + n2->Y() + n3->Y()) / 3;
    zc = (n1->Z() + n2->Z() + n3->Z()) / 3;
    aVec[0] = 0.5 * ( (n1->Y() - n3->Y() ) * (n2->Z() - n3->Z() ) 
            - (n1->Z() - n3->Z() ) * (n2->Y() - n3->Y() ));
    aVec[1] = 0.5 * ( (n1->Z() - n3->Z() ) * (n2->X() - n3->X() ) 
            - (n1->X() - n3->X() ) * (n2->Z() - n3->Z() ));
    aVec[2] = 0.5 * ( (n1->X() - n3->X() ) * (n2->Y() - n3->Y() ) 
            - (n1->Y() - n3->Y() ) * (n2->X() - n3->X() ));
    volume += xc * aVec[0] + yc * aVec[1] + zc * aVec[2];


    return volume;
}

void Tetra::ComputeCentroid() {
    cx = 0, cy = 0, cz = 0;
    for (Node* n : pNodes) {
        cx += n->X();
        cy += n->Y();
        cz += n->Z();
    }

    size_t nNodes = pNodes.size();
    cx /= nNodes;
    cy /= nNodes;
    cz /= nNodes;
}

vector<vector<size_t>> Tetra::GetEdges() const {
    vector<vector<size_t>> edgeI;

    for (size_t i = 0; i < 3; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        unsigned long n2Ind = i == 2 ? 0 : i+1;
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
    for (size_t i = 0; i < 3; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        unsigned long n2Ind = 3;
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

vector<double> Tetra::AreaVector() const {
    return vector<double>();
}

void Tetra::FindNeighbors(unsigned long iCell) {
    iNbrs.clear();
    size_t nNodes = iNodes.size();
    const vector<unsigned long> &n4Cells = pNodes[3]->Cell();
    for (size_t j = 0; j < 3; j++) {
        const unsigned long n1_i = j;
        const unsigned long n2_i = j + 1;
        const vector<unsigned long> &n1Cells = pNodes[n1_i]->Cell();
        const vector<unsigned long> &n2Cells = pNodes[n2_i]->Cell();

        for (size_t ii = 0; ii < n1Cells.size(); ii++) {
            for (size_t jj = 0; jj < n2Cells.size(); jj++) {
                if (n1Cells[ii] == n2Cells[jj] && n1Cells[ii] != iCell) 
                    iNbrs.push_back(n1Cells[ii]);
            }
            for (size_t kk = 0; kk < n4Cells.size(); kk++) {
                if (n1Cells[ii] == n4Cells[kk] && n1Cells[ii] != iCell)
                    iNbrs.push_back(n1Cells[ii]);
            }
        }
    }
}

void Tetra::ComputeDirectedArea(void) {
    Edge *e;
    Node *n1, *n2, *n3, *n4;
    vector<double> em, f1m, f2m, aVec;
    vector<double> cellCtr{cx, cy, cz};
    double coeff = 0;

    n1 = pNodes[0];
    n2 = pNodes[1];
    n3 = pNodes[2];
    n4 = pNodes[3];

    // Find edge 1
    for (Edge *pEdge : pEdges) {
        unsigned long iNode1 = pEdge->Node1();
        unsigned long iNode2 = pEdge->Node2();
        if (iNodes[0] == iNode1 && iNodes[1] == iNode2) {
            coeff = 1;
            e = pEdge;
            break;
        }
        else if (iNodes[1] == iNode1 && iNodes[0] == iNode2) {
            coeff = -1;
            e = pEdge;
            break;
        }
    }

    em = e->Midpoint();
    f1m = TriCentroid(n1, n2, n3);
    f2m = TriCentroid(n1, n2, n4);

    aVec = TriAreaVector(em, f1m, cellCtr, coeff);
    //aVec = TriAreaVector(em, cellCtr, f1m, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << f1m[0] << "\t" << f1m[1] << "\t" << f1m[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";

    aVec = TriAreaVector(em, cellCtr, f2m, coeff);
    //aVec = TriAreaVector(em, f2m, cellCtr, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << endl;
    //cout << f2m[0] << "\t" << f2m[1] << "\t" << f2m[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";

    // Find edge 2
    for (Edge *pEdge : pEdges) {
        unsigned long iNode1 = pEdge->Node1();
        unsigned long iNode2 = pEdge->Node2();
        if (iNodes[1] == iNode1 && iNodes[2] == iNode2) {
            coeff = 1;
            e = pEdge;
            break;
        }
        else if (iNodes[2] == iNode1 && iNodes[1] == iNode2) {
            coeff = -1;
            e = pEdge;
            break;
        }
    }
    em = e->Midpoint();
    f1m = TriCentroid(n1, n2, n3);
    f2m = TriCentroid(n2, n3, n4);

    aVec = TriAreaVector(em, f1m, cellCtr, coeff);
    //aVec = TriAreaVector(em, cellCtr, f1m, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << f1m[0] << "\t" << f1m[1] << "\t" << f1m[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";
    aVec = TriAreaVector(em, cellCtr, f2m, coeff);
    //aVec = TriAreaVector(em, f2m, cellCtr, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << endl;
    //cout << f2m[0] << "\t" << f2m[1] << "\t" << f2m[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";

    // Find edge 3
    for (Edge *pEdge : pEdges) {
        unsigned long iNode1 = pEdge->Node1();
        unsigned long iNode2 = pEdge->Node2();
        if (iNodes[2] == iNode1 && iNodes[0] == iNode2) {
            coeff = 1;
            e = pEdge;
            break;
        }
        else if (iNodes[0] == iNode1 && iNodes[2] == iNode2) {
            coeff = -1;
            e = pEdge;
            break;
        }
    }
    em = e->Midpoint();
    f1m = TriCentroid(n1, n2, n3);
    f2m = TriCentroid(n1, n3, n4);
    aVec = TriAreaVector(em, f1m, cellCtr, coeff);
    //aVec = TriAreaVector(em, cellCtr, f1m, coeff);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << f1m[0] << "\t" << f1m[1] << "\t" << f1m[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";
    e->AddDirArea(aVec);
    aVec = TriAreaVector(em, cellCtr, f2m, coeff);
    //aVec = TriAreaVector(em,  f2m, cellCtr, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << endl;
    //cout << f2m[0] << "\t" << f2m[1] << "\t" << f2m[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";

    // Find edge 4
    for (Edge *pEdge : pEdges) {
        unsigned long iNode1 = pEdge->Node1();
        unsigned long iNode2 = pEdge->Node2();
        if (iNodes[0] == iNode1 && iNodes[3] == iNode2) {
            coeff = 1;
            e = pEdge;
            break;
        }
        else if (iNodes[3] == iNode1 && iNodes[0] == iNode2) {
            coeff = -1;
            e = pEdge;
            break;
        }
    }
    em = e->Midpoint();
    f1m = TriCentroid(n1, n3, n4);
    f2m = TriCentroid(n1, n2, n4);
    //aVec = TriAreaVector(em, f1m, cellCtr, coeff);
    aVec = TriAreaVector(em, cellCtr, f1m, coeff);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << f1m[0] << "\t" << f1m[1] << "\t" << f1m[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";
    e->AddDirArea(aVec);
    //aVec = TriAreaVector(em, cellCtr, f2m, coeff);
    aVec = TriAreaVector(em, f2m, cellCtr, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << endl;
    //cout << f2m[0] << "\t" << f2m[1] << "\t" << f2m[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";

    // Find edge 5
    for (Edge *pEdge : pEdges) {
        unsigned long iNode1 = pEdge->Node1();
        unsigned long iNode2 = pEdge->Node2();
        if (iNodes[1] == iNode1 && iNodes[3] == iNode2) {
            coeff = 1;
            e = pEdge;
            break;
        }
        else if (iNodes[3] == iNode1 && iNodes[1] == iNode2) {
            coeff = -1;
            e = pEdge;
            break;
        }
    }
    em = e->Midpoint();
    f1m = TriCentroid(n2, n3, n4);
    f2m = TriCentroid(n1, n2, n4);
    aVec = TriAreaVector(em, f1m, cellCtr, coeff);
    //aVec = TriAreaVector(em, cellCtr, f1m, coeff);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << f1m[0] << "\t" << f1m[1] << "\t" << f1m[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";
    e->AddDirArea(aVec);
    aVec = TriAreaVector(em, cellCtr, f2m, coeff);
    //aVec = TriAreaVector(em, f2m, cellCtr, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << endl;
    //cout << f2m[0] << "\t" << f2m[1] << "\t" << f2m[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";

    // Find edge 6
    for (Edge *pEdge : pEdges) {
        unsigned long iNode1 = pEdge->Node1();
        unsigned long iNode2 = pEdge->Node2();
        if (iNodes[2] == iNode1 && iNodes[3] == iNode2) {
            coeff = 1;
            e = pEdge;
            break;
        }
        else if (iNodes[3] == iNode1 && iNodes[2] == iNode2) {
            coeff = -1;
            e = pEdge;
            break;
        }
    }
    em = e->Midpoint();
    f1m = TriCentroid(n2, n3, n4);
    f2m = TriCentroid(n1, n3, n4);
    //aVec = TriAreaVector(em, f1m, cellCtr, coeff);
    aVec = TriAreaVector(em, cellCtr, f1m, coeff);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << f1m[0] << "\t" << f1m[1] << "\t" << f1m[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";
    e->AddDirArea(aVec);
    //aVec = TriAreaVector(em, cellCtr, f2m, coeff);
    aVec = TriAreaVector(em, f2m, cellCtr, coeff);
    e->AddDirArea(aVec);
    //cout << em[0] << "\t" << em[1] << "\t" << em[2] << endl;
    //cout << cellCtr[0] << "\t" << cellCtr[1] << "\t" << cellCtr[2] << endl;
    //cout << f2m[0] << "\t" << f2m[1] << "\t" << f2m[2] << "\n\n";
    //cout << aVec[0] << "\t" << aVec[1] << "\t" << aVec[2] << "\n\n\n\n";
}

void Tetra::NodeNeighbors(void) {
    for (size_t i = 0; i < 3; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        unsigned long n2_ind = i == 2 ? 0 : i+1;
        Node *n2 = pNodes[n2_ind];

        unsigned long n2_I = iNodes[n2_ind];
        if (!n1->HasNeighbor(n2_I)) n1->AddNeighbor(n2_I);
        if (!n2->HasNeighbor(n1_I))n2->AddNeighbor(n1_I);
    }
    for (size_t i = 0; i < 3; i++) {
        Node *n1 = pNodes[i];
        unsigned long n1_I = iNodes[i];
        Node *n2 = pNodes[3];

        unsigned long n2_I = iNodes[3];
        if (!n1->HasNeighbor(n2_I)) n1->AddNeighbor(n2_I);
        if (!n2->HasNeighbor(n1_I))n2->AddNeighbor(n1_I);
    }
}
