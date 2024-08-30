#include "Grid.hpp"
#include <cstdlib>

Grid::Grid() { 
    config = &Config::GetConfig(); 
    nDim = config->GetNumDims();
}
Grid::~Grid() { }

const vector<Node>& Grid::Nodes() const { return nodes; }
const vector<unique_ptr<Cell>>& Grid::Cells() const { return cells; }
const vector<Edge>& Grid::Edges() const { return edges; }
const vector<unique_ptr<Boundary>>& Grid::Boundaries() const { return bndrys; }

void Grid::BuildDualGrid() {
    CellData();
    EdgeData();
    BoundaryData();
}

void Grid::CellData() {
    // Compute dual volume, centroids, and add cells to nodes 
    unsigned long nCells = cells.size();
    for (size_t i = 0; i < nCells; i++) {
        auto &cell = cells[i];
        zdouble vol = cell->ComputeVolume();
        cell->ComputeCentroid();
        vector<Node*> &cNodes = cell->pNodes;
        zdouble nodeAvgVol = vol / cNodes.size();

        for (Node *n : cNodes) {
            n->iCells.push_back(i);
            n->dualVol += nodeAvgVol;
        }
    }

    // Find neighbor cells
    for (size_t i = 0; i < nCells; i++) {
        auto &cell = cells[i];
        cell->FindNeighbors(i);
    }
}

void Grid::EdgeData(void) {
    size_t nEdges = 0;
    unsigned long nNodes = nodes.size();
    for (const auto &cell : cells) cell->NodeNeighbors();

    size_t tot_nbrs = 0;
    for (auto i = 0ul; i < nNodes; i++) {
        const Node &n = nodes[i];
        vector<unsigned long> iNodes = n.iNbrs;
        tot_nbrs += iNodes.size();
        for (const auto iNode : iNodes) 
            if (iNode > i) nEdges++;
    }

    edges.resize(nEdges);

    cout << "nEdges: " << nEdges << endl;
    cout << "totNbrs: " << tot_nbrs << endl;

    unsigned long eCnt = 0;
    for (auto i = 0ul; i < nNodes; i++) {
        const Node &n = nodes[i];
        vector<unsigned long> iNodes = n.iNbrs;
        tot_nbrs += iNodes.size();
        for (const auto iNode : iNodes) 
            if (iNode > i) {
                Edge &e = edges[eCnt++];
                e.iNode1 = i;
                e.iNode2 = iNode;
            }
    }

    for (auto i = 0ul; i < nEdges; i++) {
        Edge &e = edges[i];
        const Node &n1 = nodes[e.iNode1];
        const Node &n2 = nodes[e.iNode2];
        e.ComputeMetrics(n1, n2);
        const auto &n1Cells = n1.iCells;
        const auto &n2Cells = n2.iCells;
        for (size_t ii = 0; ii < n1Cells.size(); ii++) {
            const unsigned long &n1Cell = n1Cells[ii];
            for (size_t jj = 0; jj < n2Cells.size(); jj++) {
                const unsigned long &n2Cell = n2Cells[jj];
                if (n1Cell == n2Cell) cells[n1Cell]->AddEdge(i, &e);
            }
        }
    }

    for (const auto &cell : cells) cell->ComputeDirectedArea();
    for (Edge &e : edges) e.NormalizeDirAreaVec();
}

/* Needs cell-specific implementations */
void Grid::BoundaryData() {

    // Init local vector of dual areas
    vector<vector<zdouble>> areaVectors(nodes.size() );
    for (size_t i = 0; i < nodes.size(); i++) {
        areaVectors[i].resize(3);
        areaVectors[i][0] = 0;
        areaVectors[i][1] = 0;
        areaVectors[i][2] = 0;
    }

    for (auto &b : bndrys) {

        size_t nbfaces = b->cells.size();
        set<unsigned long> uniqueNodes;
        for (size_t i = 0; i < nbfaces; i++) {
            auto &c = b->cells[i];
            vector<zdouble> areaVec = c->AreaVector();
            size_t nCellNodes = c->iNodes.size();
            for (const auto &iNode : c->iNodes) { 
                uniqueNodes.insert(iNode);
                areaVectors[iNode][0] += areaVec[0] / nCellNodes;
                areaVectors[iNode][1] += areaVec[1] / nCellNodes;
                areaVectors[iNode][2] += areaVec[2] / nCellNodes;
            }
        }

        vector<unsigned long> iNodes;
        iNodes.assign(uniqueNodes.begin(), uniqueNodes.end() );
        b->iNodes.swap(iNodes);

        size_t nbNodes = b->iNodes.size();
        b->nodeNorms.resize(nbNodes);
        b->nodeAreas.resize(nbNodes);
        for (size_t i = 0; i < nbNodes; i++) {
            const unsigned long &iNode = b->iNodes[i];
            const vector<zdouble> &area = areaVectors[iNode];
            zdouble area_mag = sqrt(area[0] * area[0] + area[1] * area[1] + area[2] * area[2]);
            b->nodeAreas[i] = area_mag;

            b->nodeNorms[i].resize(3);
            b->nodeNorms[i][0] = area[0] / area_mag;
            b->nodeNorms[i][1] = area[1] / area_mag;
            b->nodeNorms[i][2] = area[2] / area_mag;

            areaVectors[iNode][0] = 0;
            areaVectors[iNode][1] = 0;
            areaVectors[iNode][2] = 0;
        }
    }
}

void Grid::FlipBoundaryNorms() {
    auto &bound = bndrys[0];
    auto &norms = bound->nodeNorms;
    for (auto &norm : norms) {
        for (auto i = 0ul; i < norm.size(); i++) norm[i] *= -1.0;
    }
    /*
    for (auto &bound : bndrys) {
        auto &norms = bound->nodeNorms;
        for (auto &norm : norms) {
            for (auto i = 0ul; i < norm.size(); i++) norm[i] *= -1.0;
        }
    }
    */
}
