#include "VTK.hpp"

void VTK::WriteVTK(const string &filename, const Grid &grid, const Solver &solver) {
    ofstream f(filename);
    const vector<unique_ptr<Cell>> &cells = grid.Cells();
    const vector<Node> &nodes = grid.Nodes();

    Config &config = Config::GetConfig();
    unsigned short nVar = config.GetNumVars();
    unsigned short nDim = config.GetNumDims();

    unsigned long nCells = cells.size();
    unsigned long nNodes = nodes.size();

    // Dummy data for boundary testing. Remove this and the body data below
    /*
    vector<zdouble> sol(nNodes, 1);
    const vector<unique_ptr<Boundary>> &boundaries = grid.Boundaries();
    cout << "Hello :)\n";
    cout << "nBound: " << boundaries.size() << endl;
    int counter = 0;
    for (const auto &b : boundaries) {
        zdouble b_val = 0;

        cout << b->Name() << "\t" << counter << endl;

        if (counter == 0) b_val = 5000.0;
        else if (counter == 1) b_val = 0.0;
        else if (counter == 2) b_val = -5000.0;
        else b_val = 0.0;

        counter++;
        const vector<unique_ptr<Cell>> &cells = b->Cells();
        for (const auto &c : cells) {
            const vector<unsigned long> &conn = c->GetConnectivity();
            for (const auto &i : conn) {
                sol[i] = b_val;
            }
        }
        //ofstream f(b.Name() );
        //f.close();
    }
*/
    
    const SolVector &w = solver.GetPrimitives();

    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
    f << Indent() << "<UnstructuredGrid>\n"; 
    f << Indent(2) << "<Piece NumberOfPoints=" << "\"" << nNodes << "\"" << " NumberOfCells=" 
        << "\"" << nCells << "\">\n"; 

    f << Indent(3) << "<Points>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Node &n : nodes) {
        f << Indent(5) << n.X() << " " << n.Y() << " " << n.Z() << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</Points>\n";

    // Dummy data to check boundaries
    /*
    f << Indent(3) << "<PointData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"test\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        //f << Indent(5) << sol[i] << "\n";
        f << Indent(5) << sol[i] << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</PointData>\n";
    */

    f << Indent(3) << "<PointData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        f << Indent(5) << w[i][0] << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
            for (unsigned short j = 0; j < nDim; j++) {
                f << Indent(5) << w[i][j+1] << " ";
            }
            if (nDim == 2) f << Indent(5) << 0 << endl;
            else f << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        f << Indent(5) << w[i][nVar-1] << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(3) << "</PointData>\n";

    f << Indent(3) <<  "<Cells>\n";
    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const unique_ptr<Cell> &cell : cells) {
        const vector<unsigned long> &conn = cell->GetConnectivity();
        f << Indent(5);
        for (size_t i = 0; i < conn.size(); i++) f << conn[i] << " ";
        f << endl;
    }

    f << Indent(4) << "</DataArray>\n";
    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    f << Indent(5);
    for (const unique_ptr<Cell> &cell : cells) {
        offset += cell->GetConnectivity().size();
        f << offset << " ";
    }
    f << endl;

    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    f << Indent(5);
    for (const unique_ptr<Cell> &cell : cells) {
        int vtkType = GetVTKCellType(cell->GetType());
        f << vtkType << " ";
    }
    f << endl;

    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</Cells>\n";
    f << Indent(3) << "<CellData>\n";
    f << Indent(3) << "</CellData>\n";
    f << Indent(2) << "</Piece>\n";
    f << Indent() << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    f.close();
}

void VTK::InsertIndent(ofstream &file, size_t nIndent) {
    string indent(nIndent * 2, ' ');
    file << indent;
}

string VTK::Indent(size_t nIndent) {
    return string(nIndent * 2, ' ');
}

int VTK::GetVTKOffset(const Cell::Type &cellType) {
    switch (cellType) {
        case (Cell::Type::BAR_2) : return 2; break;
        case (Cell::Type::TRI_3) : return 3; break;
        case (Cell::Type::QUAD_4) : return 4; break;
        case (Cell::Type::TETRA_4) : return 4; break;
        case (Cell::Type::HEXA_8) : return 8; break;
        default : 
            cout << "Unsupported VTK cell type " << Cell::CellTypeName(cellType) << endl; 
            return -1; 
    }
}

int VTK::GetVTKCellType(const Cell::Type &cellType) {
    switch (cellType) {
        case (Cell::Type::BAR_2) : return 3; break;
        case (Cell::Type::TRI_3) : return 5; break;
        case (Cell::Type::QUAD_4) : return 9; break;
        case (Cell::Type::TETRA_4) : return 10; break;
        case (Cell::Type::HEXA_8) : return 12; break;
        default : 
            cout << "Unsupported VTK cell type " << Cell::CellTypeName(cellType) << endl; 
            return -1; 
    }
}
