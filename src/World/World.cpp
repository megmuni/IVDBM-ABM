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
using namespace std;

#include "World.h"
 
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