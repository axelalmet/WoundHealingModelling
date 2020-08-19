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

#include "PolarityTrackingModifier.hpp"
#include <algorithm>
#include "ExtracellularMatrixCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "BloodCellProliferativeType.hpp"
#include "CollagenCellMutationState.hpp"
#include "FibrinCellMutationState.hpp"
#include "FibroblastStateDependentCollagenSrnModel.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PolarityTrackingModifier<DIM>::PolarityTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mVelocityReorientationStrength(DOUBLE_UNSET),
      mFibreReorientationStrength(DOUBLE_UNSET),
      mNeighbourhoodRadius(DOUBLE_UNSET),
      mMorphogenThreshold(DOUBLE_UNSET),
      mMeanActivationLifetime(DOUBLE_UNSET),
      mFibreDepositionProbability(DOUBLE_UNSET),
      mCollagenProductionRate(DOUBLE_UNSET),
      mCollagenDegradationRate(DOUBLE_UNSET)
{
}

template<unsigned DIM>
PolarityTrackingModifier<DIM>::~PolarityTrackingModifier()
{
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetVelocityReorientationStrength()
{
    return mVelocityReorientationStrength;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetVelocityReorientationStrength(double velocityReorientationStrength)
{
    mVelocityReorientationStrength = velocityReorientationStrength;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetFibreReorientationStrength()
{
    return mFibreReorientationStrength;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetFibreReorientationStrength(double fibreReorientationStrength)
{
    mFibreReorientationStrength = fibreReorientationStrength;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetNeighbourhoodRadius()
{
    return mNeighbourhoodRadius;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetNeighbourhoodRadius(double neighbourhoodRadius)
{
    mNeighbourhoodRadius = neighbourhoodRadius;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetMorphogenThreshold()
{
    return mMorphogenThreshold;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetMorphogenThreshold(double morphogenThreshold)
{
    mMorphogenThreshold = morphogenThreshold;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetMeanActivationLifetime()
{
    return mMeanActivationLifetime;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetMeanActivationLifetime(double meanActivationLifetime)
{
    mMeanActivationLifetime = meanActivationLifetime;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetFibreDepositionProbability()
{
    return mFibreDepositionProbability;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetFibreDepositionProbability(double fibreDepositionProbability)
{
    mFibreDepositionProbability = fibreDepositionProbability;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetCollagenProductionRate()
{
    return mCollagenProductionRate;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetCollagenProductionRate(double collagenProductionRate)
{
    mCollagenProductionRate = collagenProductionRate;
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetCollagenDegradationRate()
{
    return mCollagenDegradationRate;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetCollagenDegradationRate(double collagenDegradationRate)
{
    mCollagenDegradationRate = collagenDegradationRate;
}

template<unsigned DIM>
c_vector<double, 2> PolarityTrackingModifier<DIM>::GetLocalFibreOrientationAndDensity(unsigned fibroblastIndex, AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{    
    c_vector<double, 2> fibre_angle_and_density = zero_vector<double>(2); // Initialise the vector storing the angle and density

    c_vector<double, DIM> fibroblast_location = rCellPopulation.GetNode(fibroblastIndex)->rGetLocation();

    double normalisation_factor = 0.0;

    // This really only works with node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

    std::set<unsigned> neighbour_indices = p_cell_population->GetNodesWithinNeighbourhoodRadius(fibroblastIndex, mNeighbourhoodRadius);

    // As there may be non-zero fibrin and collagen densities,
    // we will just sum the amounts
    for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
    iter != neighbour_indices.end();
    ++iter)
    {
        if (p_cell_population->IsCellAttachedToLocationIndex(*iter))
        {
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(*iter);
            boost::shared_ptr<AbstractCellProliferativeType> p_neighbour_cell_type = p_cell->GetCellProliferativeType();

            // Only consider fibrin or collagen neighbours (obviously)
            if (p_neighbour_cell_type->IsType<ExtracellularMatrixCellProliferativeType>())
            {
                boost::shared_ptr<CellData> cell_data = p_cell->GetCellData();

                double current_angle = cell_data->GetItem("direction"); // Get the fibre orientation
                double fibre_density = cell_data->GetItem("density"); // Fibre density is just the sum of the amounts
                c_vector<double, DIM> neighbour_location = rCellPopulation.GetNode(*iter)->rGetLocation(); // Get the neighbour location

                // Get the distance between the neighbour and the fibroblast
                c_vector<double, DIM> fibroblast_to_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(fibroblast_location, neighbour_location);
                double distance_to_neighbour = norm_2(fibroblast_to_neighbour);

                if (distance_to_neighbour < 1e-4) // If we're sitting on an ECM fibre, we take that location/density without interpolation
                {
                    fibre_angle_and_density[0] = current_angle;
                    fibre_angle_and_density[1] = fibre_density;
                    break; // Stop the loop.
                }
                else
                {
                    // Update the fibre angle
                    fibre_angle_and_density[0] += pow(mNeighbourhoodRadius - distance_to_neighbour, 2.0)/pow(mNeighbourhoodRadius * distance_to_neighbour, 2.0) * current_angle;

                    // Update the fibre density
                    fibre_angle_and_density[1] += pow(mNeighbourhoodRadius - distance_to_neighbour, 2.0)/pow(mNeighbourhoodRadius * distance_to_neighbour, 2.0) * fibre_density;

                    // Update the normalisation factor
                    normalisation_factor += pow(mNeighbourhoodRadius - distance_to_neighbour, 2.0)/pow(mNeighbourhoodRadius * distance_to_neighbour, 2.0);

                }
            }
        }
    }

    // Normalise the angle and density
    fibre_angle_and_density /= normalisation_factor;

    return fibre_angle_and_density;
}

template<unsigned DIM>
bool PolarityTrackingModifier<DIM>::DoesCellIntersectWithFibre(CellPtr pFibroblastCell, c_vector<double, DIM> fibroblastLocation,
                                    CellPtr pMatrixCell, c_vector<double, DIM> matrixCellLocation)
{
    double cell_direction = pFibroblastCell->GetCellData()->GetItem("direction"); // Get the cell direction
    double fibre_direction = pMatrixCell->GetCellData()->GetItem("direction"); // Get the fibre directions

    // If the angles are the same, then they'll only intersect if they span from teh same point
    if (cell_direction == fibre_direction)
    {
        return ( (fibroblastLocation[0] == matrixCellLocation[0])&&(fibroblastLocation[1] == matrixCellLocation[1]) );
    }
    else
    {
        return ( (fabs(((fibroblastLocation[1] - matrixCellLocation[1])*cos(fibre_direction) 
                - (fibroblastLocation[0] - matrixCellLocation[0])*sin(fibre_direction))
                    /sin(fibre_direction - cell_direction)) <= mNeighbourhoodRadius)
                ||(fabs(((fibroblastLocation[1] - matrixCellLocation[1])*cos(fibre_direction)
                + (fibroblastLocation[0] - matrixCellLocation[0])*sin(fibre_direction))
                /sin(fibre_direction - cell_direction)) <= mNeighbourhoodRadius) );
    }
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    // UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    boost::shared_ptr<AbstractCellProperty> p_collagen_state(CellPropertyRegistry::Instance()->Get<CollagenCellMutationState>()); // ECM collagen mutation state
    boost::shared_ptr<AbstractCellProperty> p_ecm_type(CellPropertyRegistry::Instance()->Get<ExtracellularMatrixCellProliferativeType>()); // ECM collagen mutation state

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only consider non-ECM cells for polarity alignment
        boost::shared_ptr<AbstractCellProliferativeType> p_cell_type = cell_iter->GetCellProliferativeType();

        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get node index
        Node<DIM>* p_node = rCellPopulation.GetNode(node_index); // Get node
        c_vector<double, DIM> node_location = p_node->rGetLocation();

        boost::shared_ptr<CellData> cell_data = cell_iter->GetCellData();
        double phi = cell_data->GetItem("direction"); // Get the current cell direction
        double dt = SimulationTime::Instance()->GetTimeStep(); // Get dt

        // Get the morphogen levels, activation status, etc
        double activation_time = cell_data->GetItem("activation time"); // Activation lifetime
        double current_time = SimulationTime::Instance()->GetTime(); // Current time

        if (!p_cell_type->IsType<ExtracellularMatrixCellProliferativeType>())       
        {
            if (p_cell_type->IsType<FibroblastCellProliferativeType>()) // For fibroblasts, we align cells with the local fibre orientation
            { 

                double activated_status = cell_data->GetItem("activated");
                double morphogen = cell_data->GetItem("morphogen");

                if (activated_status == 1.0) // If the fibroblast is activated
                {
                    if (morphogen < mMorphogenThreshold) // We should initiative the deactivation of the fibroblast
                    {
                        // Generate a random activation lifetime
                        double uniform_random_number = RandomNumberGenerator::Instance()->ranf();
                        double time_of_deactivation = current_time - mMeanActivationLifetime * log(uniform_random_number);

                        cell_iter->GetCellData()->SetItem("activation time", time_of_deactivation); 
                    }
                    else // The fibroblast is still activated and we should reset any deactivation time
                    {
                        if (activation_time > 0.0)
                        {
                            cell_iter->GetCellData()->SetItem("activation time", 0.0); // Reset the activation time
                        }
                    }

                    // // Determine if we should deposit an ECM fibre
                    // double probability_of_deposition = RandomNumberGenerator::Instance()->ranf();

                    // if (probability_of_deposition <= mFibreDepositionProbability)
                    // {
                    //     // Create the SRN model and cell cycle model
                    //     FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
                    //     p_srn_model->SetOdeParameters(mCollagenProductionRate, mCollagenDegradationRate);

                    //     NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Place-holder cell cycle model
                    //     p_cycle_model->SetDimension(2);

                    //     // Let's copy the cell data from the fibroblast cell to this cell
                    //     CellPropertyCollection p_new_fibre_collection;

                    //     // Create a new cell data object using the copy constructor and add this to the daughter cell
                    //     MAKE_PTR_ARGS(CellData, p_new_fibre_cell_data, (*cell_data));
                    //     p_new_fibre_collection.AddProperty(p_new_fibre_cell_data);

                    //     // Copy all cell Vec data (note we create a new object not just copying the pointer)
                    //     if (cell_iter->rGetCellPropertyCollection().template HasPropertyType<CellVecData>())
                    //     {
                    //         // Get the existing copy of the cell data and remove it from the daughter cell
                    //         boost::shared_ptr<CellVecData> cell_vec_data = cell_iter->GetCellVecData();

                    //         // Create a new cell data object using the copy constructor and add this to the daughter cell
                    //         MAKE_PTR_ARGS(CellVecData, p_new_fibre_cell_vec_data, (*cell_vec_data));
                    //         p_new_fibre_collection.AddProperty(p_new_fibre_cell_vec_data);
                    //     }

                    //     CellPtr p_new_fibre_cell(new Cell(p_collagen_state, p_cycle_model, p_srn_model, false, p_new_fibre_collection));
                    //     p_new_fibre_cell->InitialiseSrnModel();
                    //     p_new_fibre_cell->SetCellProliferativeType(p_ecm_type);
                    //     p_new_fibre_cell->GetCellData()->SetItem("scale", 1.0);
                    //     p_new_fibre_cell->GetCellData()->SetItem("activated", 1.0);
                    //     p_new_fibre_cell->GetCellData()->SetItem("activation time", 0.0);

                    //     // rCellPopulation.AddCell(p_new_fibre_cell, *cell_iter);
                    //     // unsigned node_index = rCellPopulation.GetNumAllCells();
                    //     // unsigned node_index = rCellPopulation.GetNumNodes();
                    //     // unsigned node_index = rCellPopulation.GetNumAliveCells();
                    //     rCellPopulation.AddCell(p_new_fibre_cell, *cell_iter); 
                    //     rCellPopulation.Update();
                    // }

                    // If there's a non-zero activation time, it means that the fibroblast has been primed to deactivate
                    if (activation_time > 0.0)
                    {
                        if (current_time > activation_time)
                        {
                            cell_iter->GetCellData()->SetItem("activated", 0.0); // Deactivate the fibroblast
                            cell_iter->GetCellData()->SetItem("activation time", 0.0); // Reset the activation time
                        }
                    }

                }
                else // Else, we check if we should activate the fibroblast
                {
                    if (morphogen >= mMorphogenThreshold)
                    {
                        cell_iter->GetCellData()->SetItem("activated", 1.0); // Activate the fibroblast
                        cell_iter->GetCellData()->SetItem("activation time", 0.0); // Reset the activation time (if we need to)
                    }
                }
                
                // Remodel the fibroblast direction to the local fibre angle
                c_vector<double, 2> fibre_angle_and_density = GetLocalFibreOrientationAndDensity(node_index, rCellPopulation);
                double fibre_angle = fabs(fibre_angle_and_density[0]); // Get the local fibre angle (we're going to consider all )
                double fibre_density = fibre_angle_and_density[1]; // Get the fibre density too

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
            else if (p_cell_type->IsType<BloodCellProliferativeType>() )
            {
                if (current_time > activation_time)
                {
                    cell_iter->Kill();
                    // rCellPopulation.Update();
                }
            }

            c_vector<double, DIM> velocity = p_node->rGetAppliedForce(); // The applied force is proportional to the vector

            if (norm_2(velocity) > 0.0)
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
        else if (p_cell_type->IsType<ExtracellularMatrixCellProliferativeType>() ) // ECM fibres are locally aligned by activated (myo)fibroblasts
        {
            double phi_new = phi; // We need to initialise a new fibre angle to account for multiple fibroblasts remodellign the local orientation

            unsigned num_activated_fibroblast_neighbours = 0; 

            // This really only works with node-based cell populations
            assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

            NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

            std::set<unsigned> neighbour_indices = p_cell_population->GetNodesWithinNeighbourhoodRadius(node_index, mNeighbourhoodRadius);

            // As there may be non-zero fibrin and collagen densities,
            // we will just sum the amounts
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
            iter != neighbour_indices.end();
            ++iter)
            {
                if (p_cell_population->IsCellAttachedToLocationIndex(*iter))
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
                                &&(DoesCellIntersectWithFibre(p_cell, fibroblast_location, *cell_iter, node_location)) ) // Fibroblast intersects with fibre
                        {
                            num_activated_fibroblast_neighbours += 1;

                            // If it's a fibrin mutation state, alter it to collagen.
                            if (cell_iter->GetMutationState()->template IsType<FibrinCellMutationState>() )
                            {
                                cell_iter->SetMutationState(p_collagen_state);
                            }

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

                            phi_new += exp(-5.0 * norm_2(fibroblast_location - node_location))* mFibreReorientationStrength * dt * sin(fibroblast_direction - phi); // Update the fibre orientation
                        }
                    }
                }
            }

            // If the ECM node has interacted with an activated fibroblast, it will begin to produce/degrade in ECM fibres
            // Otherwise, switch it off
            if (num_activated_fibroblast_neighbours > 0)
            {
                cell_iter->GetCellData()->SetItem("activated", 1.0);
            }  
            else
            {
                cell_iter->GetCellData()->SetItem("activated", 0.0);
            }
        
            cell_iter->GetCellData()->SetItem("direction", phi_new);

        }
    }

}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityTrackingModifier<1>;
template class PolarityTrackingModifier<2>;
template class PolarityTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityTrackingModifier)
