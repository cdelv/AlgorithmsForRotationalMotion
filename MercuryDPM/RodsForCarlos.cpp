//Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Example 6 - Single clump in the periodic box

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Math/Helpers.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Particles/MultiParticle.h"
#include "clump/clump_io.h"
#include "clump/mercury3Dclump.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

Mdouble f_min = -0.5; Mdouble f_max = 1.5; // Size of the box and the margin/clearance for clump seeds

int N_att = 30;   // Number of attempts to add particle

class multiParticleT1 : public Mercury3Dclump
{
public:
    explicit  multiParticleT1()
    {
        setGravity(Vec3D(0.0, 0.0, -9.81));
        setName("RodsForCarlos");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(f_max);
        setYMax(f_max);
        setZMax(f_max); // Unbounded domain
        setXMin(f_min);
        setYMin(f_min);
        setZMin(f_min);
        load_clumps(data);
        clump_mass = data.mass[clump_index]* MatDensity;
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
        auto species = speciesHandler.copyAndAddObject(HertzianViscoelasticMindlinSpecies());
        species->setDensity(MatDensity); // sets the species type-0 density
        species->setDissipation(0.0);
        species->setSlidingDissipation(0.0);
        Mdouble ElasticModulus = 1e9; // taking 100 times that of soft material
        Mdouble PoissonRatio = 0.30;
        Mdouble SlidingFrictionCoefficient = 0.0;
        Mdouble ShearModulus = ElasticModulus/ 2 /(1+PoissonRatio);
        Mdouble EffectiveShearModulus = ShearModulus/(2-PoissonRatio)/2;
        Mdouble EffectiveElasticModulus = ElasticModulus/(1-PoissonRatio*PoissonRatio)/2;
        species->setEffectiveElasticModulus(EffectiveElasticModulus);
        species->setEffectiveShearModulus(EffectiveShearModulus);
        species->setSlidingFrictionCoefficient(SlidingFrictionCoefficient);
        
        // Rayleigh timestep
        Mdouble MinRadius = *min_element(data.pebbles_r[clump_index].begin(), data.pebbles_r[clump_index].end()); //in m
        Mdouble timeStep = 0.1 * helpers::getRayleighTime(MinRadius,EffectiveShearModulus, PoissonRatio, MatDensity);
        std::cout<<"time step"<< timeStep <<std::endl;
        Mdouble DT = 3.795117781647068e-05;
        timeStep = 0.7*DT;
        setClumpDamping(0);
        setTimeStep(timeStep);
        std::cout<<"time step" << getTimeStep()<<std::endl;

        /* Double periodic + bottom wall + unlimited top
        auto per_x = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_x->set(Vec3D(1, 0, 0), getXMin(), getXMax());
        auto per_y = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
        per_y->set(Vec3D(0, 1, 0), getYMin(), getYMax());
        wallHandler.clear();
        InfiniteWall w0;
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        */

         //Rectangular box
        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);
        
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);

        // Generate a dense packing of clumps
        setClumpIndex(3);

        string fname = "Clump.csv";

        vector<vector<std::string>> content;
        vector<string> row;
        string line, word;

        fstream file (fname, ios::in);
        if(file.is_open())
        {
            while(getline(file, line))
            {
                row.clear();

                stringstream str(line);

                while(getline(str, word, ','))
                    row.push_back(word);
                content.push_back(row);
            }
        }

        for (int part = 0; part<N_att; part++) {
            MultiParticle p0;
            p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to MultiParticle 1
            p0.setMaster();
            Mdouble qw = std::atof(content[part+1][9].c_str());
            Mdouble q1 = std::atof(content[part+1][10].c_str());
            Mdouble q2 = std::atof(content[part+1][11].c_str());
            Mdouble q3 = std::atof(content[part+1][12].c_str());
            Quaternion q(qw,q1,q2,q3);
            Vec3D v(1.0,1.0,1.0);
            q.rotate(v);
            q.rotateBack(v);
            clump_data rdata = rotate_clump(data, clump_index, uniform_random_pds()); // Rotate clump arbitrarily


            p0.setRadius(rdata.pebbles_r[clump_index][0]);

            for (int j = 0; j < rdata.pebbles_r[clump_index].size(); j++) {
                p0.addSlave(Vec3D(rdata.pebbles_x[clump_index][j],
                                  rdata.pebbles_y[clump_index][j],
                                  rdata.pebbles_z[clump_index][j]),
                            rdata.pebbles_r[clump_index][j]);
            }
            p0.setPrincipalDirections(
                    Matrix3D(rdata.pd[clump_index][0], rdata.pd[clump_index][1], rdata.pd[clump_index][2],
                             rdata.pd[clump_index][3], rdata.pd[clump_index][4], rdata.pd[clump_index][5],
                             rdata.pd[clump_index][6], rdata.pd[clump_index][7], rdata.pd[clump_index][8]));
            p0.setInitInertia(
                    MatrixSymmetric3D(rdata.toi[clump_index][0] * MatDensity, rdata.toi[clump_index][1] * MatDensity, rdata.toi[clump_index][2] * MatDensity,
                                      rdata.toi[clump_index][4] * MatDensity, rdata.toi[clump_index][5] * MatDensity,
                                      rdata.toi[clump_index][8] * MatDensity));
            p0.setMassMultiparticle(rdata.mass[clump_index] * MatDensity);

            p0.setDamping(clump_damping);


            Vec3D Pos1 = (Vec3D(std::atof(content[part+1][0].c_str()), std::atof(content[part+1][1].c_str()), std::atof(content[part+1][2].c_str()))); //add them to vector3d


            p0.setPosition(Pos1);


            Vec3D vel = (Vec3D(std::atof(content[part+1][3].c_str()), std::atof(content[part+1][4].c_str()), std::atof(content[part+1][5].c_str()))); //add them to vector3d
            Vec3D angVel = (Vec3D(std::atof(content[part+1][6].c_str()), std::atof(content[part+1][7].c_str()), std::atof(content[part+1][8].c_str()))); //add them to vector3d

            p0.setAngularVelocity(angVel);
            p0.setVelocity(vel);


            if (checkClumpForInteractionPeriodic(p0)) {
                particleHandler.copyAndAddObject(p0);
            }
        }
    }
    void printTime() const
    {
        // std::cout << "t=" << getTime() << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
        // // Calculate Kinetic Energy of Clumps
        Mdouble KineticEnergyClump = 0.0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->IsMaster())
            {
                Mdouble KineticTranslational = .5 * p->getMass() * p->getVelocity().getLengthSquared();
                Mdouble KineticRotational = .5 * Vec3D::dot(p->getAngularVelocity(), p->getInertia() * p->getAngularVelocity());
                KineticEnergyClump += (KineticTranslational+KineticRotational);
                // std::cout << "Ene" << KineticRotational/KineticTranslational << std::endl;
            }
            // std::cout << "KineticEnergyClump" << KineticEnergyClump << std::endl;
        } 

        Mdouble ElasticEnergy = 0.0;
        for (BaseInteraction* c : interactionHandler)
        {
            ElasticEnergy += c->getElasticEnergy();
            // std::cout << "ElasticEnergy" << ElasticEnergy << std::endl;
        } 
        std::cout << "t=" << getTime() << " Ene " << KineticEnergyClump/getElasticEnergy() << std::endl;
    }     

private:
    int clump_index;
    clump_data data;
    Mdouble clump_mass;
    Mdouble clump_damping = 0;
    Mdouble MatDensity = 2500;
    bool writeVTK = true;
    unsigned int wallsBoundary = 0; //0->PBC , 1->wall

};

int main(int argc, char* argv[])
{

    // for(int i=0;i<content.size();i++)
    // {
    //     for(int j=0;j<content[i].size();j++)
    //     {
    //         cout<<content[i][j]<<" ";
    //     }
    //     cout<<"\n";
    // }

    // return 0;
    multiParticleT1 problem;

    // Quick demonstration
    problem.setSaveCount(200000);
    problem.setTimeMax(250.0);

    problem.removeOldFiles();
    problem.solve();
    return 0;
}
