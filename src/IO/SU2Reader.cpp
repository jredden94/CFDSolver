#include "SU2Reader.hpp"
#include <string>

SU2Reader::SU2Reader() { }
SU2Reader::~SU2Reader() { }

string SU2Reader::FindTag(const string& tag, ifstream &file) const {
    string line;
    while (true) {
        getline(file, line);
        istringstream iss(line);
        if (line.empty() || line[0] == '%') continue;
        if (line.substr(0, tag.length()) == tag) {
            size_t pos = line.find('=');
            return (line.substr(pos+1));
            break;
        }
    }
}

void SU2Reader::ReadFile() {
    ifstream file(filename);
    string line;

    nDim = stoi(FindTag("NDIME", file));
    nCells = stoi(FindTag("NELEM", file));
    cells.resize(nCells);

    vector<vector<unsigned long>> connectivity;
    connectivity.resize(nCells);
    for (auto i = 0ul; i < nCells; i++) connectivity[i].resize(15);
    for (auto iLine = 0ul; iLine < nCells; iLine++) {
        getline(file, line);
        istringstream iss(line);
        unsigned long value;
        size_t counter = 0;

        while (iss >> value) connectivity[iLine][counter++] = value;
    }

    nNodes = stoi(FindTag("NPOIN", file));
    nodes.resize(nNodes);

    double coord[15];
    if (nDim == 3) {
        for (auto iNode = 0ul; iNode < nNodes; iNode++) {
            getline(file, line);
            istringstream iss(line);
            size_t count = 0;
            double value;
            while (iss >> value) coord[count++] = value;

            Node &node = nodes[iNode];
            node.SetX(coord[0]);
            node.SetY(coord[1]);
            node.SetZ(coord[2]);
        }
    }
    else {
        for (auto iNode = 0ul; iNode < nNodes; iNode++) {
            getline(file, line);
            istringstream iss(line);
            size_t count = 0;
            double value;
            while (iss >> value) coord[count++] = value;

            Node &node = nodes[iNode];
            node.SetX(coord[0]);
            node.SetY(coord[1]);
        }
    }

    nBoundaries = stoi(FindTag("NMARK", file));

    unsigned long row[15];
    for (auto i = 0ul; i < nBoundaries; i++) {
        string bName = FindTag("MARKER_TAG", file);
        unsigned long nbCells = stoi(FindTag("MARKER_ELEMS", file));
        vector<unique_ptr<Cell>> bCells;
        bCells.resize(nbCells);

        for (auto iLine = 0ul; iLine < nbCells; iLine++) {
            getline(file, line);
            istringstream iss(line);
            size_t count = 0;
            unsigned long value;

            while (iss >> value) row[count++] = value;

            int vtk = row[0];
            vector<unsigned long> conn; 
            vector<Node*> node_ptr;
            if (vtk == 3) {
                conn.resize(2);
                node_ptr.resize(2);
                conn[0] = row[1];
                conn[1] = row[2];
                node_ptr[0] = &nodes[conn[0]];
                node_ptr[1] = &nodes[conn[1]];

                bCells[iLine] = Cell::CreateCell(vtk, conn, node_ptr);
            }
        }

        boundaries.push_back(make_unique<Boundary>(bName, bCells));
    }

    for (auto i = 0ul; i < nCells; i++) {
        vector<unsigned long> &conn = connectivity[i];
        vector<unsigned long> cell_conn;
        vector<Node*> node_ptr;
        int vtk = conn[0];
        if (vtk == 5) {
            cell_conn.resize(3);
            node_ptr.resize(3);
            cell_conn[0] = conn[1];
            cell_conn[1] = conn[2];
            cell_conn[2] = conn[3];
            node_ptr[0] = &nodes[cell_conn[0]];
            node_ptr[1] = &nodes[cell_conn[1]];
            node_ptr[2] = &nodes[cell_conn[2]];

            cells[i] = Cell::CreateCell(vtk, cell_conn, node_ptr);
        }
        else if (vtk == 9) {
            cell_conn.resize(4);
            node_ptr.resize(4);
            cell_conn[0] = conn[1];
            cell_conn[1] = conn[2];
            cell_conn[2] = conn[3];
            cell_conn[3] = conn[4];
            node_ptr[0] = &nodes[cell_conn[0]];
            node_ptr[1] = &nodes[cell_conn[1]];
            node_ptr[2] = &nodes[cell_conn[2]];
            node_ptr[3] = &nodes[cell_conn[3]];

            cells[i] = Cell::CreateCell(vtk, cell_conn, node_ptr);
        }
    }

}
