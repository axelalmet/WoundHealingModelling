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

#include "CollagenAlignmentTrackingModifier.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CollagenAlignmentTrackingModifier<DIM>::CollagenAlignmentTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mNeighbourhoodRadius(1.5),
      mReorientationStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
CollagenAlignmentTrackingModifier<DIM>::~CollagenAlignmentTrackingModifier()
{
}

template<unsigned DIM>
double CollagenAlignmentTrackingModifier<DIM>::GetNeighbourhoodRadius()
{
    return mNeighbourhoodRadius;
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::SetNeighbourhoodRadius(double neighbourhoodRadius)
{
    mNeighbourhoodRadius = neighbourhoodRadius;
}

template<unsigned DIM>
double CollagenAlignmentTrackingModifier<DIM>::GetReorientationStrength()
{
    return mReorientationStrength;
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::SetReorientationStrength(double reorientationStrength)
{
    mReorientationStrength = reorientationStrength;
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // This really only works with node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

    double neighbourhood_radius = GetNeighbourhoodRadius();

    // First check the current attachments to the basement membrane, which depends on the cell type.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only consider fibroblasts
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        if (p_cell_type->template IsType<FibroblastCellProliferativeType>())
        {

            // Get the node index
            unsigned node_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);
                
            //Sort the neighbours by collagen amounts.
            std::vector<std::pair<double, unsigned> > sorted_neighbours_by_collagen;

            // Get the set of neighbouring location indices within a neighbourhood radius
            std::set<unsigned> neighbour_indices = p_cell_population->GetNodesWithinNeighbourhoodRadius(node_index, neighbourhood_radius);

            // We only update the orientation if, well, we can.
            if (!neighbour_indices.empty())
            {
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    boost::shared_ptr<AbstractCellProperty> p_neighbour_cell_type = p_cell->GetCellProliferativeType();

                    // Only look at fibroblast neighbours.
                    if (p_neighbour_cell_type->template IsType<FibroblastCellProliferativeType>())
                    {
                        double neighbour_collagen = p_cell->GetCellData()->GetItem("collagen");

                        if (neighbour_collagen > 0.0) // We should only orient towards neighbours which actually are enriched with collagen
                        {
                            std::pair<double, unsigned> collagen_index = std::make_pair(neighbour_collagen, *iter);
                            sorted_neighbours_by_collagen.push_back(collagen_index); //Add the collagen amount and index
                        }

                    }
                }
            }

            //Sort the vector by the collagen amount in decreasing order
            std::sort(sorted_neighbours_by_collagen.rbegin(), sorted_neighbours_by_collagen.rend());

            // Once again, only do something if we have a sufficient number of collagen neighbours, two
            if (sorted_neighbours_by_collagen.size() > 1)
            {
                // Get the relevant neighbours, i.e. the neighbours with the highest and second highest
                // amount of collagen
                std::pair<double, unsigned> first_neighbour = sorted_neighbours_by_collagen[0];
                std::pair<double, unsigned> second_neighbour = sorted_neighbours_by_collagen[1];

                // Get the neighbour indices
                unsigned first_neighbour_index = first_neighbour.second;
                unsigned second_neighbour_index = second_neighbour.second;

                // Get their locations
                c_vector<double, DIM> first_neighbour_location = rCellPopulation.GetNode(first_neighbour_index)->rGetLocation();
                c_vector<double, DIM> second_neighbour_location = rCellPopulation.GetNode(second_neighbour_index)->rGetLocation();

                // Determine the collagen direction
                c_vector<double, DIM> collagen_direction = rCellPopulation.rGetMesh().GetVectorFromAtoB(second_neighbour_location, first_neighbour_location);
                
                // Compare this to the considered node location
                double orientation = cell_iter->GetCellData()->GetItem("orientation");

                c_vector<double, DIM> current_orientation;
                current_orientation[0] = cos(orientation);
                current_orientation[1] = sin(orientation);

                // If the inner product is negative, reverse the collagen direction, so as to minimise the change in orientation
                if (inner_prod(collagen_direction, current_orientation) < 0.0)
                {
                    collagen_direction *= -1.0;
                }

                collagen_direction /= norm_2(collagen_direction);

                // Get the current node location
                c_vector<double, DIM> current_location = rCellPopulation.GetNode(node_index)->rGetLocation();

                PRINT_2_VARIABLES(current_location[0], current_location[1]);
                PRINT_2_VARIABLES(collagen_direction[0], collagen_direction[1]);

                // Get the orientation of the constructed collagen fibre.
                double collagen_orientation = atan(collagen_direction[1]/collagen_direction[0]); //Get initial angle argument

                if (collagen_direction[0] < 0.0) //If the point is in the second quadrant or third quadrant
                {
                    collagen_orientation += M_PI;
                }
                else if ((collagen_direction[0]>=0.0)&&(collagen_direction[1]<0.0)) //Fourth quadrant
                {
                    collagen_orientation += M_PI;
                }

                // Finally, calculate the new orientation
                double reorientation_strength = GetReorientationStrength();
                double new_orientation = orientation + reorientation_strength * sin(collagen_orientation - orientation);
                
                // Update the cell data
                cell_iter->GetCellData()->SetItem("orientation", new_orientation);
            }
            else if (sorted_neighbours_by_collagen.size() == 1)
            {
                // Get the one collagen-expressing neighbour
                std::pair<double, unsigned> neighbour = sorted_neighbours_by_collagen[0];

                // Get the neighbour indices
                unsigned neighbour_index = neighbour.second;

                // Get their locations
                c_vector<double, DIM> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

                // Get the current node location
                c_vector<double, DIM> current_location = rCellPopulation.GetNode(node_index)->rGetLocation();

                // Determine the collagen direction
                c_vector<double, DIM> collagen_direction = rCellPopulation.rGetMesh().GetVectorFromAtoB(current_location, neighbour_location);

                // Compare this to the considered node location
                double orientation = cell_iter->GetCellData()->GetItem("orientation");

                c_vector<double, DIM> current_orientation;
                current_orientation[0] = cos(orientation);
                current_orientation[1] = sin(orientation);

                // If the inner product is negative, reverse the collagen direction, so as to minimise the change in orientation
                if (inner_prod(collagen_direction, current_orientation) < 0.0)
                {
                    collagen_direction *= -1.0;
                }

                collagen_direction /= norm_2(collagen_direction);

                PRINT_2_VARIABLES(current_location[0], current_location[1]);
                PRINT_2_VARIABLES(collagen_direction[0], collagen_direction[1]);

                // Get the orientation of the constructed collagen fibre.
                double collagen_orientation = atan(collagen_direction[1]/collagen_direction[0]); //Get initial angle argument

                if ((collagen_direction[0] < 0.0)&&(collagen_direction[1] >= 0.0) ) //If the point is in the second quadrant or third quadrant
                {
                    collagen_orientation += M_PI;
                }
                else if ((collagen_direction[0]>=0.0)&&(collagen_direction[1]<0.0)) //Fourth quadrant
                {
                    collagen_orientation += M_PI;
                }

                // Finally, calculate the new orientation
                double reorientation_strength = GetReorientationStrength();
                double new_orientation = orientation + reorientation_strength * sin(collagen_orientation - orientation);
                
                // Update the cell data
                cell_iter->GetCellData()->SetItem("orientation", new_orientation);
            }
        }

    }
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CollagenAlignmentTrackingModifier<1>;
template class CollagenAlignmentTrackingModifier<2>;
template class CollagenAlignmentTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CollagenAlignmentTrackingModifier)
