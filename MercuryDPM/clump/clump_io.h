// This module loads clump configuration produced by MClump tool

#ifndef CLUMP_IO_H
#define CLUMP_IO_H

#include <dirent.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include <cmath>
#include "Math/Matrix.h"
#include<CMakeDefinitions.h>

typedef std::vector<double> dvec;
typedef std::vector<dvec> ddvec;
typedef std::vector<std::string> svec;

struct clump_data     // structure for storing clump instances' parameters
{
public:
    
	std::string path;	// Path to MClump working directory
	svec clump_names;	// Array of names of clumps that will be used 
	dvec mass;		//  clump mass
	ddvec  pebbles_x;	//  Pebbles geometry (outer index goes over clumps, inner - over pebbles)
	ddvec  pebbles_y;
	ddvec  pebbles_z;
	ddvec  pebbles_r;
	
	ddvec  toi;		// Clump tensor of inertia (I11, I12, I13, I21..I33)
	ddvec  pd;		// Clump principal directions v1, v2, v3
	
};


void load_conf(clump_data &a)
{
    // Path to MCLump tool
    a.path = getMercurySourceDir() + "/Tools/MClump/clumps/";

    a.clump_names.push_back("clump01");
    a.clump_names.push_back("clump02");
    a.clump_names.push_back("clump03");
    a.clump_names.push_back("clumpCarlos");
}


void load_pebbles(clump_data &a)
{
    std::cout<<"Loading clump pebbles...";
    
    a.pebbles_x.resize(a.clump_names.size());
    a.pebbles_y.resize(a.clump_names.size());
    a.pebbles_z.resize(a.clump_names.size());
    a.pebbles_r.resize(a.clump_names.size());
    
    for (int i = 0; i < a.clump_names.size(); i++ ){
    	//std::cout<<a.path + a.clump_names[i]+"/clump/clump.txt"<<std::endl;
    	std::ifstream infile((a.path + a.clump_names[i]+"/clump/clump.txt").c_str(), std::ios::in | std::ios::binary);
    	char lin[256];
    	while (infile.getline(lin, 256, '\n')){
    		svec val;
    		std::string line(lin);
    		line+=",";
    		std::string buffer = "";
    		for(int j=0; j < line.size(); j++){
    			if (line[j] != ',') {buffer += line[j];}
        		else {val.push_back(buffer); buffer = "";}
        	}
        	a.pebbles_x[i].push_back(std::stof(val[0]));
        	a.pebbles_y[i].push_back(std::stof(val[1]));
        	a.pebbles_z[i].push_back(std::stof(val[2]));
        	a.pebbles_r[i].push_back(std::stof(val[3]));	
    	}
    	infile.close();	
    }
    std::cout<<"\t OK"<<std::endl;
}


void load_mass(clump_data &a)
{
    std::cout<<"Loading clump masses..";
    
    a.mass.resize(a.clump_names.size());
    
    for (int i = 0; i < a.clump_names.size(); i++ ){
    	//std::cout<<a.path + a.clump_names[i]+"/inertia/mass.txt"<<std::endl;
    	std::ifstream infile((a.path + a.clump_names[i]+"/inertia/mass.txt").c_str(), std::ios::in | std::ios::binary);
    	char lin[256];
    	infile.getline(lin, 256, '\n');
    	std::string line(lin);
    	a.mass[i] = std::stof(line);		
    	infile.close();	
    }
    std::cout<<"\t OK"<<std::endl;
}


void load_toi(clump_data &a)
{
    std::cout<<"Loading clump TOI..";
    a.toi.resize(a.clump_names.size());
    for (int i = 0; i < a.clump_names.size(); i++ ){
    	//std::cout<<a.path + a.clump_names[i]+"/inertia/toi.txt"<<std::endl;
    	std::ifstream infile((a.path + a.clump_names[i]+"/inertia/toi.txt").c_str(), std::ios::in | std::ios::binary);
    	char lin[256];
    	for (int row = 0; row<3; row++){
    		infile.getline(lin, 256, '\n');
    		svec val;
    		std::string line(lin);
    		line+=",";
    		std::string buffer = "";
    		for(int j=0; j < line.size(); j++){
    			if (line[j] != ',') {buffer += line[j];}
        		else { val.push_back(buffer); buffer = "";}
        	}
        	for (int k = 0; k < 3; k++) a.toi[i].push_back(std::stof(val[k]));
        }
    	infile.close();	
    }
    std::cout<<"\t OK"<<std::endl;
}


void load_pd(clump_data &a)
{
    std::cout<<"Loading clump PD..";
    a.pd.resize(a.clump_names.size());
    for (int i = 0; i < a.clump_names.size(); i++ ){
    	//std::cout<<a.path + a.clump_names[i]+"/inertia/pd.txt"<<std::endl;
    	std::ifstream infile((a.path + a.clump_names[i]+"/inertia/pd.txt").c_str(), std::ios::in | std::ios::binary);
    	char lin[256];
    	for (int row = 0; row<3; row++){
    		infile.getline(lin, 256, '\n');	
    		svec val;
    		std::string line(lin);
    		line+=",";
    		std::string buffer = "";
    		for(int j=0; j < line.size(); j++){
    			if (line[j] != ',') {buffer += line[j];}
        		else { val.push_back(buffer); buffer = "";}
        	}
        	for (int k = 0; k < 3; k++) a.pd[i].push_back(std::stof(val[k]));
        }
    	infile.close();	
    }
    std::cout<<"\t OK"<<std::endl;
}

void load_clumps(clump_data &data)
{
	// Umbrela function that loads all the necessary files to initiate clumps
	std::cout<<"LOAD CLUMP DATA"<<std::endl;
	load_conf(data);
	load_pebbles(data);
	load_mass(data);
	load_toi(data);
	load_pd(data);

    std::cout<<"LOADED CLUMPS"<<std::endl;
    for (int i = 0; i < data.pebbles_x.size(); i++){
        std::cout<<data.clump_names[i]<<" mass:"<<data.mass[i]<<std::endl;
        std::cout<<data.clump_names[i]<<" list of pebbles:"<<std::endl;
        for (int j = 0; j< data.pebbles_x[i].size(); j++){
            std::cout<<"Pebble "<<j<<": ("<<data.pebbles_x[i][j]<<","<<data.pebbles_y[i][j]<<","<<data.pebbles_z[i][j]<<"),"<<data.pebbles_r[i][j]<<std::endl;
        }

        std::cout<<data.clump_names[i]<<" TOI:"<<std::endl;
        std::cout<<data.toi[i][0]<<","<<data.toi[i][1]<<","<<data.toi[i][2]<<std::endl;
        std::cout<<data.toi[i][3]<<","<<data.toi[i][4]<<","<<data.toi[i][5]<<std::endl;
        std::cout<<data.toi[i][6]<<","<<data.toi[i][7]<<","<<data.toi[i][8]<<std::endl;

        std::cout<<data.clump_names[i]<<" Principal directions:"<<std::endl;
        std::cout<<data.pd[i][0]<<","<<data.pd[i][1]<<","<<data.pd[i][2]<<std::endl;
        std::cout<<data.pd[i][3]<<","<<data.pd[i][4]<<","<<data.pd[i][5]<<std::endl;
        std::cout<<data.pd[i][6]<<","<<data.pd[i][7]<<","<<data.pd[i][8]<<std::endl;


    }

}

Matrix3D transpose(Matrix3D M) { return Matrix3D(M.XX, M.YX, M.ZX, M.XY, M.YY, M.ZY, M.XZ, M.YZ, M.ZZ);}

clump_data rotate_clump(clump_data data, int clump_index, dvec new_pd)
{
    clump_data new_data = data;

    // Prepare rotation matrices
    Vec3D e10 = Vec3D(data.pd[clump_index][0], data.pd[clump_index][1], data.pd[clump_index][2]);
    Vec3D e20 = Vec3D(data.pd[clump_index][3], data.pd[clump_index][4], data.pd[clump_index][5]);
    Vec3D e30 = Vec3D(data.pd[clump_index][6], data.pd[clump_index][7], data.pd[clump_index][8]);

    Vec3D e1  = Vec3D(new_pd[0], new_pd[1], new_pd[2]);
    Vec3D e2  = Vec3D(new_pd[3], new_pd[4], new_pd[5]);
    Vec3D e3  = Vec3D(new_pd[6], new_pd[7], new_pd[8]);

    Matrix3D Q(Vec3D::dot(e10, e1), Vec3D::dot(e10, e2), Vec3D::dot(e10, e3),
               Vec3D::dot(e20, e1), Vec3D::dot(e20, e2), Vec3D::dot(e20, e3),
               Vec3D::dot(e30, e1), Vec3D::dot(e30, e2), Vec3D::dot(e30, e3));

    Matrix3D Qt = transpose(Q);

    // Set new pd's
    new_data.pd[clump_index] = new_pd;

    // Rotate tensor of inertia
    Matrix3D iI = Matrix3D(data.toi[clump_index][0], data.toi[clump_index][1], data.toi[clump_index][2],
                                data.toi[clump_index][3], data.toi[clump_index][4], data.toi[clump_index][5],
                                data.toi[clump_index][6], data.toi[clump_index][7], data.toi[clump_index][8]);

    Matrix3D nI = Q * iI * Qt;
    dvec d{ nI.XX, nI.XY, nI.XZ, nI.YX, nI.YY, nI.YZ, nI.ZX, nI.ZY, nI.ZZ };
    new_data.toi[clump_index] = d;

    return new_data; // All loaded clumps in new_data remain unchanged except the clump_index one that is rotated to new_pd
}

double random_double(double Max)
{return Max*((double) rand() / (RAND_MAX));}

dvec uniform_random_pds(){

    Vec3D n1, n2, n3, ref;

    // basis vector n1
    double r1 = random_double(2) - 1.0;
    double r2 = random_double(1);

    double theta = acos(r1); // Note that for isotropy of n1 theta is NOT uniform!
    double phi = 2 * M_PI * r2;

    n1 = Vec3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

    // basis vector n2
    r1 = random_double(1);
    r2 = random_double(1);

    theta = acos(r1); // Note that for isotropy of n1 theta is NOT uniform!
    phi = 2 * M_PI * r2;

    ref = Vec3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    n2 = Vec3D::cross(ref, n1); n2.normalise();
    n3 = Vec3D::cross(n1, n2); n3.normalise();
    return dvec{n1.X, n1.Y, n1.Z, n2.X, n2.Y, n2.Z, n3.X, n3.Y, n3.Z};

}


#endif // CLUMP_IO_H
