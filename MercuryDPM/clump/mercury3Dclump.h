#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/MultiParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include <stdlib.h>
// This class contains necessary modifications of Mercury3D framework to enable rigid clumps
class Mercury3Dclump : public Mercury3D
{
public:
    explicit  Mercury3Dclump()
    {

    }
    // Redefine force computation for clumps
    void computeInternalForce(BaseParticle* P1, BaseParticle* P2) override
    {
        // Quit if:
        // 1) at least one particle is master
        // 2) pebbles of the same clump
        // 3) both particles are fixed
        // 4) both particles are "ghosts".
        if (P1->IsMaster() || P2->IsMaster()) return;
        if (P1->IsSlave() && P2->IsSlave() && P1->getMaster() == P2->getMaster()) return;
        if (P1->isFixed() && P2->isFixed()) return;
        if ((P1->getPeriodicFromParticle() != nullptr) && (P2->getPeriodicFromParticle() != nullptr)) return;

        // Assign pointers PI, PJ such that PI has the lower id than PJ
        BaseParticle * PI, * PJ;
        if (P1->getId() > P2->getId()) {PI = P2; PJ = P1;}
        else                           {PI = P1; PJ = P2;}

        // Check for the interaction between  PI and PJ
        BaseInteraction* i = PJ->getInteractionWith(PI, getNumberOfTimeSteps() + 1,
                                                    &interactionHandler);
        if (i!= nullptr) {
            i->computeForce(); // Compute and apply forces
            PI->addForce(i->getForce());
            PJ->addForce(-i->getForce());

            // Add extra torques to master particles arising from the forces acting on pebbles
            if (getRotation()) {
                PI->addTorque(i->getTorque() - Vec3D::cross(PI->getPosition() - i->getContactPoint(), i->getForce()));
                PJ->addTorque(-i->getTorque() + Vec3D::cross(PJ->getPosition() - i->getContactPoint(), i->getForce()));
            }
        }
    }

    void computeForcesDueToWalls(BaseParticle* pI, BaseWall* w) override
    {
        // No interaction between walls and master particles
        if (pI->IsMaster())
            return;

        //No need to compute interactions between periodic particle images and walls
        if (pI->getPeriodicFromParticle() != nullptr)
            return;

        //Checks if the particle is interacting with the current wall
        BaseInteraction* i = w->getInteractionWith(pI, getNumberOfTimeSteps() + 1,
                                                   &interactionHandler);
        if (i!=nullptr) {
            //...calculates the forces between the two objects...
            i->computeForce();

            //...and applies them to each of the two objects (wall and particle).
            pI->addForce(i->getForce());
            w->addForce(-i->getForce());

            //If the rotation flag is on, also applies the relevant torques
            //(getRotation() returns a boolean).
            if (getRotation()) // getRotation() returns a boolean.
            {
                pI->addTorque(i->getTorque() - Vec3D::cross(pI->getPosition() - i->getContactPoint(), i->getForce()));
                ///\todo TW: I think this torque has the wrong sign
                w->addTorque(-i->getTorque() + Vec3D::cross(w->getPosition() - i->getContactPoint(), i->getForce()));
            }
        }
    }

    void computeAllForces() override
    {
        //Resetting all forces on both particles and walls to zero
        #pragma omp parallel num_threads(getNumberOfOMPThreads())
        {
            #pragma omp for
            for (int k = 0; k < particleHandler.getSize(); ++k) {
                particleHandler.getObject(k)->resetForceTorque(getNumberOfOMPThreads());
            }
            #pragma omp for
            for (int k = 0; k < wallHandler.getSize(); k++) {
                wallHandler.getObject(k)->resetForceTorque(getNumberOfOMPThreads());
            }
        }
        logger(DEBUG,"All forces set to zero");


        // compute all internal and external forces; for omp simulations, this can be done in parallel
        #pragma omp parallel num_threads(getNumberOfOMPThreads())
        {
            //logger(INFO, "Number of omp threads = %", getNumberOfOMPThreads());
            ///Now loop over all particles contacts computing force contributions
            #pragma omp for schedule(dynamic)
            for (int k = 0; k < particleHandler.getSize(); ++k) {
                BaseParticle *p = particleHandler.getObject(k);
                //computing both internal forces (e.g. due to collisions)
                //and external forces (e.g. gravity)
                //(compute internal forces compares the current particle p
                //with all others in the handler!)
                computeInternalForces(p);
                // body forces
                computeExternalForces(p);
            }

            // wall-forces
            #pragma omp for schedule(dynamic)
            for (int k = 0; k < wallHandler.getSize(); k++) {
                BaseWall *w = wallHandler.getObject(k);
                computeWallForces(w);
            }
        }
        #pragma omp for schedule(dynamic)
        ///Loop on pebbles to add forces/moments to masters...
        for (BaseParticle* p : particleHandler)
        {
            if ((p->IsSlave())&&(p->getPeriodicFromParticle() == nullptr))
                {
                    p->getMaster()->addForce(p->getForce());
                    Vec3D lever = p->getPosition()-p->getMaster()->getPosition();

                    // Patch to fix lever - we do not create ghost masters, that's why this workaround is needed
                    // Should affect all boundaries derived from PeriodicBoundary and do not affect other types
                    for (auto p : boundaryHandler) {
                        auto pb = dynamic_cast<PeriodicBoundary*> (p);
                        if (pb){
                            Vec3D shift = pb->getShift();
                            if (lever.getLength() > (lever - shift).getLength()) lever -= shift;
                            if (lever.getLength() > (lever + shift).getLength()) lever += shift;
                        }
                    }
                    // End patch

                    Vec3D torque = p->getTorque();
                    torque += Vec3D::cross(lever, p->getForce());
                    p->getMaster()->addTorque(torque);
            }
        }
    }

    bool checkClumpForInteraction(BaseParticle& particle)
    {
        MultiParticle* mp = dynamic_cast<MultiParticle*>(&particle);
        if (mp->IsMaster()){
            for (int i = 0; i<mp->NSlave(); i++){
                for (BaseParticle *q: particleHandler) {
                    if (!q->IsMaster()) { // check every slave vs every non-master intersection
                        if (Vec3D::getDistanceSquared(q->getPosition(), mp->getSlavePositions()[i]) < q->getRadius() + mp->getSlaveRadii()[i]) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    bool checkClumpForInteractionPeriodic(BaseParticle& particle)
    {
        // Note that this implementation only check for interaction with particles
        bool NP = true;
        // Periodic case
        for (BaseBoundary* b : boundaryHandler)
        {
            PeriodicBoundary* pb = dynamic_cast<PeriodicBoundary*>(b);
            if (pb) NP = false;
            if (pb && particleHandler.getNumberOfObjects() > 0 )
            {
                MultiParticle* mp = dynamic_cast<MultiParticle*>(&particle);
                if (mp->IsMaster()){
                    for (int i = 0; i<mp->NSlave(); i++){
                        for (BaseParticle *q: particleHandler) {
                            if (!q->IsMaster()) { // check every slave vs every non-master intersection
                                if (Vec3D::getDistanceSquared(q->getPosition(), mp->getSlavePositions()[i]) < (q->getMaxInteractionRadius() + mp->getSlaveRadii()[i]) * (q->getMaxInteractionRadius()+ mp->getSlaveRadii()[i]) ) {
                                    return false;
                                }
                                BaseParticle *q1 = q->copy(); // check every slave vs non-master ghost intersection
                                pb->shiftPosition(q1);
                                Mdouble dist2 = Vec3D::getDistanceSquared(q1->getPosition(), mp->getSlavePositions()[i]);
                                Mdouble dist02 = (q1->getMaxInteractionRadius() + mp->getSlaveRadii()[i]) * (q1->getMaxInteractionRadius() + mp->getSlaveRadii()[i]);
                                delete q1;
                                if (dist2 < dist02) {
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }
        // Non-periodic case
        if ((NP)&&(particleHandler.getNumberOfObjects() > 0 ))
        {
            MultiParticle* mp = dynamic_cast<MultiParticle*>(&particle);
            if (mp->IsMaster()){
                for (int i = 0; i<mp->NSlave(); i++){
                    for (BaseParticle *q: particleHandler) {
                        if (!q->IsMaster()) { // check every slave vs every non-master intersection
                            if (Vec3D::getDistanceSquared(q->getPosition(), mp->getSlavePositions()[i]) < (q->getRadius() + mp->getSlaveRadii()[i]) * (q->getRadius() + mp->getSlaveRadii()[i])) {
                                return false;
                            }
                        }
                    }
                }
            }
        }
        return true;
    }
};


// mathsFunc::square(sp->getSumOfInteractionRadii(q))

