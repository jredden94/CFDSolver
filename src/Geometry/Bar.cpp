#include "Cell.hpp"

/*
Bar::Bar(vector<unsigned long>& connectivity, vector<Node*>& nodePtr)
            : Cell(Type::BAR_2, connectivity, nodePtr) { }
            */
Bar::~Bar() { }

zdouble Bar::ComputeVolume(void) const { 
    cout << "Tried to calculate volume of a 'Bar' object. Returning negative value.\n";
    return -1; 
}

vector<vector<size_t>> Bar::GetEdges(void) const {
    cout << "Tried to call 'GetEdges' from a Bar object.\n";
    return vector<vector<size_t>>(); 
}

vector<zdouble> Bar::AreaVector() const {
    vector<zdouble> norm(3,0);
    const Node *n1 = pNodes[0];
    const Node *n2 = pNodes[1];
    //norm[0] = n1->Y() - n2->Y();
    //norm[1] = n2->X() - n1->X();
    norm[0] = n2->Y() - n1->Y();
    norm[1] = n1->X() - n2->X();

    return norm;
}

void Bar::ComputeCentroid(void) {
}

void Bar::FindNeighbors(unsigned long iCell) {
}

void Bar::ComputeDirectedArea(void) {
}

void Bar::NodeNeighbors(void) {
}
