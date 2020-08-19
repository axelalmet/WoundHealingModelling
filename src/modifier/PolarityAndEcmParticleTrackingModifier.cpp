/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PolarityAndEcmParticleTrackingModifier.hpp"
#include <algorithm>
#include "ExtracellularMatrixCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "BloodCellProliferativeType.hpp"
#include "CollagenCellMutationState.hpp"
#include "FibrinCellMutationState.hpp"
#include "FibroblastStateDependentCollagenSrnModel.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PolarityAndEcmParticleTrackingModifier<DIM>::PolarityAndEcmParticleTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mVelocityReorientationStrength(DOUBLE_UNSET),
      mFibreReorientationStrength(DOUBLE_UNSET),
      mNeighbourhoodRadius(DOUBLE_UNSET),
      mMorphogenThreshold(DOUBLE_UNSET),
      mMeanActivationLifetime(DOUBLE_UNSET),
      mMeanLeukocyteDeathTime(DOUBLE_UNSET),
      mFibreDepositionProbability(DOUBLE_UNSET),
      mCollagenProductionRate(DOUBLE_UNSET),
      mCollagenDegradationRate(DOUBLE_UNSET),
      mEcmParticleInformation()
{
}

template<unsigned DIM>
PolarityAndEcmParticleTrackingModifier<DIM>::~PolarityAndEcmParticleTrackingModifier()
{
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetVelocityReorientationStrength()
{
    return mVelocityReorientationStrength;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetVelocityReorientationStrength(double velocityReorientationStrength)
{
    mVelocityReorientationStrength = velocityReorientationStrength;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetFibreReorientationStrength()
{
    return mFibreReorientationStrength;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetFibreReorientationStrength(double fibreReorientationStrength)
{
    mFibreReorientationStrength = fibreReorientationStrength;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetNeighbourhoodRadius()
{
    return mNeighbourhoodRadius;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetNeighbourhoodRadius(double neighbourhoodRadius)
{
    mNeighbourhoodRadius = neighbourhoodRadius;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetMorphogenThreshold()
{
    return mMorphogenThreshold;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetMorphogenThreshold(double morphogenThreshold)
{
    mMorphogenThreshold = morphogenThreshold;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetMeanActivationLifetime()
{
    return mMeanActivationLifetime;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetMeanLeukocyteDeathTime(double meanLeukocyteDeathTime)
{
    mMeanLeukocyteDeathTime = meanLeukocyteDeathTime;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetMeanLeukocyteDeathTime()
{
    return mMeanLeukocyteDeathTime;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetMeanActivationLifetime(double meanActivationLifetime)
{
    mMeanActivationLifetime = meanActivationLifetime;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetFibreDepositionProbability()
{
    return mFibreDepositionProbability;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetFibreDepositionProbability(double fibreDepositionProbability)
{
    mFibreDepositionProbability = fibreDepositionProbability;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetCollagenProductionRate()
{
    return mCollagenProductionRate;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetCollagenProductionRate(double collagenProductionRate)
{
    mCollagenProductionRate = collagenProductionRate;
}

template<unsigned DIM>
double PolarityAndEcmParticleTrackingModifier<DIM>::GetCollagenDegradationRate()
{
    return mCollagenDegradationRate;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetCollagenDegradationRate(double collagenDegradationRate)
{
    mCollagenDegradationRate = collagenDegradationRate;
}

template <unsigned DIM>
std::map<unsigned, c_vector<double, 3> > PolarityAndEcmParticleTrackingModifier<DIM>::GetEcmParticleInformation()
{
    return mEcmParticleInformation;
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetEcmParticleInformation(std::map<unsigned, c_vector<double, 3> > ecmParticleInformation)
{
    mEcmParticleInformation = ecmParticleInformation;
}

template<unsigned DIM>
c_vector<double, 2> PolarityAndEcmParticleTrackingModifier<DIM>::GetLocalFibreOrientationAndDensity(unsigned fibroblastIndex, AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{    
    c_vector<double, 2> fibre_angle_and_density = zero_vector<double>(2); // Initialise the vector storing the angle and density

    c_vector<double, DIM> fibroblast_location = rCellPopulation.GetNode(fibroblastIndex)->rGetLocation();

    double normalisation_factor = 0.0;

    // This really only works with node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulationWithParticles<DIM>* p_cell_population = static_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation);

    std::set<unsigned> neighbour_indices = p_cell_population->GetNodesWithinNeighbourhoodRadius(fibroblastIndex, mNeighbourhoodRadius);

    // As there may be non-zero fibrin and collagen densities,
    // we will just sum the amounts
    for (auto iter = neighbour_indices.begin();
    iter != neighbour_indices.end();
    ++iter)
    {

        Node<DIM>* p_node = rCellPopulation.GetNode(*iter);

        if ( (p_node->IsParticle())&&(mEcmParticleInformation.find(*iter) != mEcmParticleInformation.end()) ) // ECM nodes should be particles only (and should be in the map)
        {
            c_vector<double, DIM> neighbour_location = p_node->rGetLocation(); // Get the neighbour location

            c_vector<double, 3> ecm_node_information = mEcmParticleInformation[*iter]; // Get the density and orientation for the ECM particle

            double fibre_angle = ecm_node_information[1];
            double fibre_density = ecm_node_information[2];

            // Get the distance between the neighbour and the fibroblast
            c_vector<double, DIM> fibroblast_to_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(fibroblast_location, neighbour_location);
            double distance_to_neighbour = norm_2(fibroblast_to_neighbour);

            if (distance_to_neighbour < 1e-4) // If we're sitting on an ECM fibre, we take that location/density without interpolation
            {
                fibre_angle_and_density[0] = fibre_angle;
                fibre_angle_and_density[1] = fibre_density;
                break; // Stop the loop.
            }
            else
            {
                // Update the fibre angle
                fibre_angle_and_density[0] += pow(mNeighbourhoodRadius - distance_to_neighbour, 2.0)/pow(mNeighbourhoodRadius * distance_to_neighbour, 2.0) * fibre_angle;

                // Update the fibre density
                fibre_angle_and_density[1] += pow(mNeighbourhoodRadius - distance_to_neighbour, 2.0)/pow(mNeighbourhoodRadius * distance_to_neighbour, 2.0) * fibre_density;

                // Update the normalisation factor
                normalisation_factor += pow(mNeighbourhoodRadius - distance_to_neighbour, 2.0)/pow(mNeighbourhoodRadius * distance_to_neighbour, 2.0);
            }
        }
    }

    // Normalise the angle and density
    fibre_angle_and_density /= normalisation_factor;

    return fibre_angle_and_density;
}

template<unsigned DIM>
bool PolarityAndEcmParticleTrackingModifier<DIM>::DoesCellIntersectWithFibre(CellPtr pFibroblastCell, c_vector<double, DIM> fibroblastLocation,
                                    unsigned fibreIndex, c_vector<double, DIM> fibreLocation)
{
    double cell_direction = pFibroblastCell->GetCellData()->GetItem("direction"); // Get the cell direction

    c_vector<double, 3> fibre_information = mEcmParticleInformation[fibreIndex];
    double fibre_direction = fibre_information[1]; // Get the fibre directions

    // If the angles are the same, then they'll only intersect if they span from teh same point
    if (cell_direction == fibre_direction)
    {
        return ( (fibroblastLocation[0] == fibreLocation[0])&&(fibroblastLocation[1] == fibreLocation[1]) );
    }
    else // Otherwise, we can solve for the intersection, and if that falls within the neighbourhood radius, we say they've intersected
    {
        return ( (fabs(((fibroblastLocation[1] - fibreLocation[1])*cos(fibre_direction) 
                - (fibroblastLocation[0] - fibreLocation[0])*sin(fibre_direction))
                    /sin(fibre_direction - cell_direction)) <= mNeighbourhoodRadius)
                ||(fabs(((fibroblastLocation[1] - fibreLocation[1])*cos(fibre_direction)
                + (fibroblastLocation[0] - fibreLocation[0])*sin(fibre_direction))
                /sin(fibre_direction - cell_direction)) <= mNeighbourhoodRadius) );
    }
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::UpdateEcmParticleIndices(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    rCellPopulation.Update();

    // This really only works with node-based cell populations with particles
    assert(dynamic_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation));

    NodeBasedCellPopulationWithParticles<DIM>* p_cell_population = static_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation);

    std::set<unsigned> particle_set = p_cell_population->GetParticleIndices();

    // Convert the set into a vector
    std::vector<unsigned> particle_indices(particle_set.begin(), particle_set.end());

    // We'll initialise a new map and fill that with the indices. If there's a better way (pre-C++17),
    // I'd love to know.
    std::map<unsigned, c_vector<double, 3> > new_ecm_map; 

    for (auto map_iter = mEcmParticleInformation.begin();
        map_iter != mEcmParticleInformation.end(); map_iter++)
    {
        unsigned particle_iter = std::distance(mEcmParticleInformation.begin(), map_iter); // Convert the iter to unsigned for the vector
        unsigned particle_index = particle_indices[particle_iter];

        new_ecm_map[particle_index] = map_iter->second; // Assign the new key

    }

    mEcmParticleInformation = new_ecm_map; // Assign the new ECM map
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::UpdateEcmParticleInformation(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    rCellPopulation.Update();

    // This really only works with node-based cell populations with particles
    assert(dynamic_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation));

    NodeBasedCellPopulationWithParticles<DIM>* p_cell_population = static_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation);

    std::set<unsigned> particle_indices = p_cell_population->GetParticleIndices();

    double dt = SimulationTime::Instance()->GetTimeStep(); // Get the timestep to update the orientation and densities of the ECM fibre

    for (auto map_iter = mEcmParticleInformation.begin();
        map_iter != mEcmParticleInformation.end(); map_iter++)
    {
        // Get the index and particle information
        unsigned ecm_node_index = map_iter->first;
        c_vector<double, 3> ecm_information = map_iter->second;

        // Get the location of the particle
        c_vector<double, DIM> ecm_node_location = rCellPopulation.GetNode(ecm_node_index)->rGetLocation();

        double ecm_type = ecm_information[0]; // 6 if collagen state, 7 if fibrin state (these values are chosen to interface with previous visualisations)
        double fibre_angle = ecm_information[1]; // Orientation if ECM fibre
        double fibre_density = ecm_information[2]; // Current density

        // Intiialise the variables for the new 
        double ecm_type_new = ecm_type;
        double fibre_angle_new = fibre_angle;
        double fibre_density_new = fibre_density;

        // We first found whether the ECM node is surrounded by activated fibroblasts.
        // If there are activated fibroblast neighbours, we update the fibre angle and density
        unsigned num_activated_fibroblast_neighbours = 0;

        // Get the neighbouring indices
        std::set<unsigned> neighbour_indices = p_cell_population->GetNodesWithinNeighbourhoodRadius(ecm_node_index, mNeighbourhoodRadius);

        for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
        iter != neighbour_indices.end();
        ++iter)
        {
            // Only consider non-particles, as these will correspond to actual cells
            Node<DIM>* p_neighbour_node = rCellPopulation.GetNode(*iter);

            if (!p_neighbour_node->IsParticle())
            {
                CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(*iter);
                boost::shared_ptr<AbstractCellProliferativeType> p_neighbour_cell_type = p_cell->GetCellProliferativeType();

                // Only consider fibrin or collagen neighbours (obviously)
                if (p_neighbour_cell_type->IsType<FibroblastCellProliferativeType>())
                { 
                    unsigned fibroblast_index = rCellPopulation.GetLocationIndexUsingCell(p_cell);
                    c_vector<double, DIM> fibroblast_location = rCellPopulation.GetNode(fibroblast_index)->rGetLocation();

                    double neighbour_activation_status = p_cell->GetCellData()->GetItem("activated");

                    if ( (neighbour_activation_status == 1.0) // Fibroblast has been activated
                            &&(DoesCellIntersectWithFibre(p_cell, fibroblast_location, ecm_node_index, ecm_node_location)) ) // Fibroblast intersects with fibre
                    {
                        num_activated_fibroblast_neighbours += 1;

                        double fibroblast_direction = p_cell->GetCellData()->GetItem("direction");

                        // Adjust the fibroblast direction so that the change in the fibre angle is minimised.
                        // This may need to be changed in the future as we're making an implicit assumption
                        // that we're only considering the principle angle for the fibre orientation here.
                        if ( (fibroblast_direction > 0.5*M_PI)&&(fibroblast_direction < 1.5*M_PI) ) // Q2 or Q3
                        {
                            fibroblast_direction -= M_PI;
                        }
                        else if (fibroblast_direction >= 1.5*M_PI) // Q4
                        {
                            fibroblast_direction -= 2.0 * M_PI;
                        }

                        fibre_angle_new += exp(-5.0 * norm_2(fibroblast_location - ecm_node_location))* mFibreReorientationStrength * dt * sin(fibroblast_direction - fibre_angle); // Update the fibre orientation
                    }
                }
            }
        }

        // Update the and type ECM density if there are activated fibroblasts near it
        if (num_activated_fibroblast_neighbours > 0)
        {
            ecm_type_new = 5; // If it's a fibrin state, alter it to collagen
            fibre_density_new += dt * (mCollagenProductionRate - mCollagenDegradationRate * fibre_density); // Update collagen production 
        }

        ecm_information[0] = ecm_type_new;
        ecm_information[1] = fibre_angle_new;
        ecm_information[2] = fibre_density_new;

        // Update the ECM particle information
        mEcmParticleInformation[ecm_node_index] = ecm_information; 
    }

}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::UpdateEcmParticleInformationAfterCellDeath(std::vector<unsigned> removedNodeIndices)
{

    // Here, for every removed node index, the particle indices AFTER that index will shift
    // to the left by one, and we need to account for this specifically in the particle indices. The cell
    // population should already account for this, (I think).

    std::sort(removedNodeIndices.begin(), removedNodeIndices.end()); // First sort the indices in increasing order

    unsigned previous_size = mEcmParticleInformation.size();

    for (unsigned i = 0; i < removedNodeIndices.size(); i++)
    {
        if (i == removedNodeIndices.size() - 1)
        {
            for (auto map_iter = mEcmParticleInformation.upper_bound(removedNodeIndices[i]); 
                    map_iter != mEcmParticleInformation.end(); map_iter++)
            {
                unsigned ecm_node_index = map_iter->first;
                c_vector<double, 3> ecm_information = map_iter->second;

                // We update the map now to shift one to the left
                mEcmParticleInformation[ecm_node_index - 1] = ecm_information;

                // Erase the current element now
                mEcmParticleInformation.erase(map_iter);
            }
        }
        else
        {
            for (auto map_iter = mEcmParticleInformation.upper_bound(removedNodeIndices[i]); 
                    map_iter != mEcmParticleInformation.lower_bound(removedNodeIndices[i + 1]); map_iter++)
            {
                unsigned ecm_node_index = map_iter->first;
                c_vector<double, 3> ecm_information = map_iter->second;

                // We update the map now to shift one to the left
                mEcmParticleInformation[ecm_node_index - 1] = ecm_information;

                // Erase the current element now
                mEcmParticleInformation.erase(map_iter);
            }
        }
    }

    unsigned current_size = mEcmParticleInformation.size();

    assert(previous_size == current_size); // Need to make sure it's doing what it's supposed to.
}


template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    UpdateEcmParticleIndices(rCellPopulation);
    UpdateCellData(rCellPopulation);
    UpdateEcmParticleInformation(rCellPopulation);

    // rCellPopulation.Update();
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);

    double time = SimulationTime::Instance()->GetTime();

    *mpCellAndEcmDataFile << time << " ";

    // Write the data to file
    WriteCellAndEcmData(rCellPopulation);

    *mpCellAndEcmDataFile << "\n";
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
    
    // Create output file
    OutputFileHandler output_file_handler(outputDirectory + "/", false);
    mpCellAndEcmDataFile = output_file_handler.OpenOutputFile("cellandecmpolarities.dat");

    double time = SimulationTime::Instance()->GetTime();

    *mpCellAndEcmDataFile << time << " ";

    // Write the data to file
    WriteCellAndEcmData(rCellPopulation);

    *mpCellAndEcmDataFile << "\n";
}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Static-cast the population, so that we can add nodes.
    NodeBasedCellPopulationWithParticles<DIM>* p_cell_population = static_cast<NodeBasedCellPopulationWithParticles<DIM>*>(&rCellPopulation);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double dt = p_simulation_time->GetTimeStep();
    double current_time = p_simulation_time->GetTime();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        // Only consider non-ECM cells for polarity alignment
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get node index
        Node<DIM>* p_node = rCellPopulation.GetNode(node_index); // Get node
        c_vector<double, DIM> node_location = p_node->rGetLocation();

        double phi = cell_iter->GetCellData()->GetItem("direction"); // Get the current cell direction

        // Get the morphogen levels, activation status, etc
        double activation_time = cell_iter->GetCellData()->GetItem("activation_time"); // Activation lifetime

        if (p_cell_type->template IsType<FibroblastCellProliferativeType>())  // For fibroblasts, we align cells with the local fibre orientation
        { 
            // double activated_status = cell_iter->GetCellData()->GetItem("activated");
            double morphogen = cell_iter->GetCellData()->GetItem("morphogen");

            // If the morphogen is above a critical threshold, keep the status as activated
            if (morphogen > mMorphogenThreshold)
            {
                double activated_status = cell_iter->GetCellData()->GetItem("activated");
                double activation_time = cell_iter->GetCellData()->GetItem("activation_time");

                if (activated_status == 0.0) // If the fibroblast needs to be reactivated
                {
                    cell_iter->GetCellData()->SetItem("activated", 1.0); // Activate the fibroblast
                }
                if (activation_time > 0.0) // If the fibroblast needs to have its reactivation time reset
                {
                    cell_iter->GetCellData()->SetItem("activation_time", 0.0); // Reset the fibroblast lifetime
                }

            }
            else
            {
                // If we haven't set an activation_time 
                if (activation_time == 0.0)
                {
                    double exponential_random_number = p_gen->ExponentialRandomDeviate(mMeanActivationLifetime);
                    double time_of_deactivation = current_time + exponential_random_number;
                    cell_iter->GetCellData()->SetItem("activation_time", time_of_deactivation); // Reset the activation_time
                    cell_iter->GetCellData()->SetItem("activated", 0.0); // Reset the activation_time

                }
                // else
                // {
                //     // If the current time has passed the activation_time
                //     if (current_time > activation_time)
                //     {
                //         cell_iter->GetCellData()->SetItem("activated", 0.0); // Deactivate the fibroblast to stop it proliferating
                //         cell_iter->GetCellData()->SetItem("activation_time", 0.0); // Reset the activation_time
                //     }
                // }
            }
                
            // Remodel the fibroblast direction to the local fibre angle
            c_vector<double, 2> fibre_angle_and_density = GetLocalFibreOrientationAndDensity(node_index, rCellPopulation);
            double fibre_angle = fabs(fibre_angle_and_density[0]); // Get the local fibre angle (we're going to consider all )
            double fibre_density = fibre_angle_and_density[1]; // Get the fibre density too

            if (!std::isnan(fibre_angle)) // In case there's no ECM particles in the neighbourhood set.
            {
                // We need to calculate the spanning fibre direction that minimises the change in the direction
                double change_in_phi = fibre_angle - phi;

                if (fabs(M_PI - fibre_angle - phi) < fabs(change_in_phi)) // Q2 adjustment
                {
                    change_in_phi = M_PI - fibre_angle - phi; 
                }
                else if (fabs(M_PI + fibre_angle - phi) < fabs(change_in_phi)) // QIII
                {
                    change_in_phi = M_PI + fibre_angle - phi;
                }
                else if (fabs(2.0 * M_PI - fibre_angle - phi ) < fabs(change_in_phi))
                {
                    change_in_phi = 2.0 * M_PI - fibre_angle - phi;
                }

                // We can finally update the orientation now
                phi += mFibreReorientationStrength * fibre_density * dt * sin(change_in_phi);  
            }              

            c_vector<double, DIM> velocity = p_node->rGetAppliedForce(); // The applied force is proportional to the vector

            if ( (norm_2(velocity) > 0.0)&&(!std::isnan(norm_2(velocity))) )
            {
                velocity /= norm_2(velocity);

                // Determine the angle the velocity direction makes with the x-axis
                double velocity_angle = atan(velocity[1]/velocity[0]);

                // Adjust for the relevant quadrants
                if (velocity[0] < 0.0) // Quadrant II or III
                {
                    velocity_angle += M_PI;
                }
                else
                {
                    if (velocity[1] < 0.0) // Quadrant IV
                    {
                        velocity_angle += 2.0 * M_PI;
                    }
                }

                phi += mVelocityReorientationStrength * dt * sin(velocity_angle - phi); // Update the cell polarity
                cell_iter->GetCellData()->SetItem("direction", phi);
            }

        }
        // For blood cells, we'll kill them if they're in the neighbourhood
        else if ((!cell_iter->HasApoptosisBegun())&&(p_cell_type->template IsType<BloodCellProliferativeType>()) )
        { 
            if (activation_time == 0.0)
            {
                std::set<unsigned> neighbouring_indices = p_cell_population->GetNeighbouringNodeIndices(node_index);

                // Now check whether there are any fibroblast neighbours.
                for (std::set<unsigned>::iterator neighbour_iter=neighbouring_indices.begin();
                        neighbour_iter != neighbouring_indices.end();
                        ++neighbour_iter)
                {
                    if (!rCellPopulation.GetNode(*neighbour_iter)->IsParticle())
                    {
                        // Only consider fibroblast neighbours
                        if (rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType()->template IsType<FibroblastCellProliferativeType>() )
                        {
                            // We kill if the fibroblast is sufficiently activated during the wounding, i.e. it has been exposed to a sufficient 
                            // amount of growth factor.
                            double activation_status = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter)->GetCellData()->GetItem("activated");
                            
                            if ( (activation_status == 1.0) )
                            {
                                // Set the death time
                                double exponential_random_number = p_gen->ExponentialRandomDeviate(mMeanLeukocyteDeathTime);
                                double time_of_deactivation = current_time + exponential_random_number;
                                cell_iter->GetCellData()->SetItem("activation_time", time_of_deactivation); // Reset the activation_time
                                cell_iter->GetCellData()->SetItem("activated", 0.0); // Reset the activation_time
                                // cell_iter->StartApoptosis();
                                // cell_iter->Kill();
                            }
                        }
                    }
                }
            }
            else
            {
                if (current_time > activation_time)
                {
                    cell_iter->Kill();
                }
            }
        }
    }
}

template <unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::WriteCellAndEcmData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    // First write the data for the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        // Only consider non-ECM cells for polarity alignment
        boost::shared_ptr<AbstractCellProliferativeType> p_cell_type = cell_iter->GetCellProliferativeType();
        unsigned cell_colour = p_cell_type->GetColour(); // Get the colour ID to determine cell type

        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get node index
        c_vector<double, DIM> node_location = rCellPopulation.GetNode(node_index)->rGetLocation(); // Get node

        double cell_orientation = cell_iter->GetCellData()->GetItem("direction"); // Get the current cell direction

        // We write the following data:
        // node index, cell colour, x, y, orientation, 0 (cells have 0 density)
        *mpCellAndEcmDataFile << node_index << " " << cell_colour << " " 
                        << node_location[0] << " " << node_location[1] << " "
                        << cell_orientation << " " << 0.0 << " ";

    }
    
    for (auto map_iter = mEcmParticleInformation.begin();
    map_iter != mEcmParticleInformation.end(); map_iter++)
    {
        // Get the index and particle information
        unsigned ecm_node_index = map_iter->first;
        c_vector<double, 3> ecm_information = map_iter->second;

        unsigned ecm_cell_colour = (unsigned)ecm_information[0];
        double ecm_fibre_orientation = ecm_information[1];
        double ecm_fibre_density = ecm_information[2];

        // Get the node location too
        c_vector<double, DIM> ecm_node_location = rCellPopulation.GetNode(ecm_node_index)->rGetLocation();

        // We write the following data:
        // node index, cell colour, x, y, orientation, fibre density
        *mpCellAndEcmDataFile << ecm_node_index << " " << ecm_cell_colour << " " 
                        << ecm_node_location[0] << " " << ecm_node_location[1] << " "
                        << ecm_fibre_orientation << " " << ecm_fibre_density << " ";
    }

}

template<unsigned DIM>
void PolarityAndEcmParticleTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityAndEcmParticleTrackingModifier<1>;
template class PolarityAndEcmParticleTrackingModifier<2>;
template class PolarityAndEcmParticleTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityAndEcmParticleTrackingModifier)
