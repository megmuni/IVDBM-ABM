/*
 * File: World.cpp 
 *
 * File Contents: Contains World class
 *
 * Author: alireza
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
using namespace std;

#include "World.h"

double World::clock = 0;
 
World::World() {
    this->nx = 0; 
    this->ny = 0; 
    this->nz = 0;
    this->seed = 27000; //srand(seed);
    srand(time(0));
}

World::World(const World& orig) {}
 
World::~World() {}

void World::setupGrid(int nx, int ny, int nz, REAL x_min, REAL x_max, REAL y_min, REAL y_max, REAL z_min, REAL z_max) {
    cout << "Setting Grid .." << endl;

    // Number of grid points (patches) in x,y,z dimensions
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;

    // Min and max coordinates (mm)
    this->x_min = x_min;
    this->x_max = x_max;
    this->y_min = y_min;
    this->y_max = y_max;
    this->z_min = z_min;
    this->z_max = z_max;

    //Length (mm)
    REAL L_x = x_max - x_min;
    REAL L_y = y_max - y_min;
    REAL L_z = z_max - z_min;

    //Length of each grid/patch (mm)     
    this->dx = (nx > 1)? L_x/(nx - 1) : 0.0;
    this->dy = (ny > 1)? L_y/(ny - 1) : 0.0;
    this->dz = (nz > 1)? L_z/(nz - 1) : 0.0;
    
    this->x.resize(nx*ny*nz);
    this->y.resize(nx*ny*nz);
    this->z.resize(nx*ny*nz);
    
    //ofstream outfile ("output/initGridData");

    // At every patch, calculate x, y, z coordinates (mm)
    for (int iz = 0; iz < nz; iz++)
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++) {
                int in = ix + iy*nx + iz*nx*ny;  //Patch row major index
                x[in] = dx*ix + x_min;
                y[in] = dy*iy + y_min;
                z[in] = dz*iz + z_min;
                //outfile << " x=" << x[in] << "  y=" << y[in] << "  z=" << z[in] << "  in=" << in << endl;
            }
            //outfile.close();
    cout << " Setting Grid completed." << endl;
}

void World::outputWorld_VTK_binary(const char* filename, double t) {
   cout << endl << setfill('-') << setw(80) << "-" << endl;
   cout << setw(15) << " Writing mesh in VTK binary format to file:" << filename << endl;
   int NPoints = nx*ny*nz;
   ofstream outfile(filename, ofstream::binary);

   int x1 = 0;
   int x2 = nx - 1;
   int y1 = 0;
   int y2 = ny - 1;
   int z1 = 0;
   int z2 = nz - 1;
   
   outfile << "<?xml version=\"1.0\"?>" << endl;
   outfile << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
   outfile << "<StructuredGrid WholeExtent=\"" << x1 << " " << x2 << " " << y1 << " " << y2 << " " << z1 << " " << z2 << "\">\n";
   outfile << "<Piece Extent=\"" << x1 << " " << x2 << " " << y1 << " " << y2 << " " << z1 << " " << z2 << "\">\n";
   outfile << "<PointData Scalars=\"default_fv\">\n";
   
   int offset = 0;
   outfile << "<DataArray type=\"Float64\" Name=\"default_fv\" format=\"appended\" encoding=\"raw\" offset=\"" << offset << "\"/>\n";
   offset += sizeof(int) + NPoints*sizeof(double);
   outfile << "</PointData>\n";
   outfile << "<CellData></CellData>\n";
   outfile << "<Points>\n";
   
   int offset_cord = offset;
   outfile << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" encoding=\"raw\" offset=\"" << offset_cord << "\"/>\n";
   outfile << "</Points>\n";
   outfile << "</Piece>\n";
   outfile << "</StructuredGrid>\n";
   outfile << "<AppendedData encoding=\"raw\">\n";
   outfile << "_";
   
   int size = NPoints*sizeof(double);
   outfile.write((char*)&size, sizeof(int));
   for (int in = 0; in < NPoints; in++)
       outfile.write(reinterpret_cast<char*>(&this->field_var[0][in]),sizeof(double));
   
   size = 3*nx*ny*nx*sizeof(double);
   cout << " NPoints=" << NPoints << endl;
   outfile.write((char*)&size, sizeof(int));
   for (int in = 0; in < NPoints; in++) {
    outfile.write((char*)&this->x[in], sizeof(double));
    outfile.write((char*)&this->y[in], sizeof(double));
    outfile.write((char*)&this->z[in], sizeof(double));
    //cout << " x= " << x[in] << "  y=" << y[in] << "  z=" << z[in] << endl;
   }
   
   outfile << "\n";
   outfile << "</AppendedData>\n";
   outfile << "</VTKFile>\n";
   outfile.close();

   cout << setw(15) << " Writing mesh in VTK binary format completed." << endl;
   cout << setfill('-') << setw(80) << "-" << endl;
}

void outputWorld_csv() {
    if (this->clock == 0) {
        remove(get_output_filename().c_str());

        ofstream output_file(get_output_filename(), ios::app);
        write_csv_header(output_file);
        output_file.close();

        write_auxiliary_header();
    }

    ofstream output_file(get_output_filename(), ios::app);

    // agent counting
    map<string, int> agent_counts;
    for (const auto& name : get_agent_type_names())
        agent_counts[name] = 0;

    count_agent_types(agent_counts);

    // derived class owns agent (e.g. cells) container, total comes via getter
    int totalAgents = get_total_agent_count();

    // env element (e.g. ecm) type counting
    map<string, float> env_counts;

    for (const auto& name : get_env_type_names())
        env_counts[name] = 0;

    count_env(env_counts);

    // print counts generically
    cout << " total agents: " << totalAgents << endl;
    for (const auto& [name, count] : agent_counts)
        cout << " " << name << " agents: " << count << endl;

    for (const auto& [name, amount] : ecm_counts)
        cout << " " << name << ": " << amount << endl;

    // shared columns — clock and day always written first
    output_file << this->clock << ","
        << this->clock / 48 << ",";

    // hook — remaining columns differ per world type
    write_data_row(output_file, agent_counts, ecm_counts);

    output_file.close();

    // hook — auxiliary outputs e.g. tgf_line, empty in base
    write_auxiliary_outputs();

    // hook — derived class updates its own prev agent count
    update_prev_agents();
}

// DEFAULT OUTPUT HOOK DEFINITIONS ///
void World::write_csv_header(ofstream& file) { }

void World::write_data_row(ofstream& file,
    map<string, int>& agent_counts,
    map<string, float>& env_counts) { }

void World::write_auxiliary_header() { }

void World::write_auxiliary_outputs() { }

vector<string> World::get_agent_type_names() { return {}; }

void World::count_agent_types(map<string, int>& counts) { }

vector<string> World::get_env_type_names() { return {}; }

void World::count_env(map<string, float>& counts) { }

int  World::get_total_agent_count() { return 0; }