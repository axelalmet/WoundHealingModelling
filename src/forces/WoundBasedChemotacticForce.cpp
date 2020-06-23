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

#include "WoundBasedChemotacticForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlateletCellProliferativeType.hpp"
#include "CollagenCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
WoundBasedChemotacticForce<DIM>::WoundBasedChemotacticForce()
    : AbstractForce<DIM>(),
    mChemotacticStrength(DOUBLE_UNSET),
    mNeighbourhoodRadius(DOUBLE_UNSET)
{
}

template<unsigned DIM>
WoundBasedChemotacticForce<DIM>::~WoundBasedChemotacticForce()
{
}

// Method to get the chemotactic strength parameter, i.e. chi
template<unsigned DIM>
double WoundBasedChemotacticForce<DIM>::GetChemotacticStrength()
{
    return mChemotacticStrength;
}

// Method to set the chemotactic strength parameter, chi
template<unsigned DIM>
void WoundBasedChemotacticForce<DIM>::SetChemotacticStrength(double chemotacticStrength)
{
    mChemotacticStrength = chemotacticStrength;
}

// Method to get the neighbourhood radius of interaction
template<unsigned DIM>
double WoundBasedChemotacticForce<DIM>::GetNeighbourhoodRadius()
{
    return mNeighbourhoodRadius;
}

// Method to set the neighbourhood radius of interaction
template<unsigned DIM>
void WoundBasedChemotacticForce<DIM>::SetNeighbourhoodRadius(double neighbourhoodRadius)
{
    mNeighbourhoodRadius = neighbourhoodRadius;
}


template<unsigned DIM>
void WoundBasedChemotacticForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    // For now, this will only apply to node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_tissue = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    double neighbourhood_radius = GetNeighbourhoodRadius();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        // Should only look at non-platelet cells and non-collagen cells
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType(); // Get the cell type

        if ( (!p_cell_type->IsType<PlateletCellProliferativeType>())
            &&(!p_cell_type->IsType<CollagenCellProliferativeType>()) )
        {
            // Get the node index
            unsigned current_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Get the location of the node
            c_vector<double, DIM> current_location = rCellPopulation.GetNode(current_index)->rGetLocation();

            // Get the current morphogen concentration that we will compare for chemotaxis.
            double current_morphogen_concentration = cell_iter->GetCellData()->GetItem("morphogen");

            // Get the neighbouring node indices
            std::set<unsigned> neighbouring_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(current_index, neighbourhood_radius);

            // We determine which neighbours result in a direction of maximal change. There may be more than one, in which case
            // we will pick one at random.

            // First figure out what the maximal change in morphogen is. 
            double morphogen_max_grad = 0.0;

            // Iterate over these elements
            for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
                    elem_iter != neighbouring_indices.end();
                    ++elem_iter)
            {
                // Get the cell according to the index
                CellPtr neighbour_cell_iter = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

                // Get the neighbouring concentration
                double neighbour_morphogen_concentration = neighbour_cell_iter->GetCellData()->GetItem("morphogen");

                // Get the neighbour's location
                c_vector<double, DIM> neighbour_location = rCellPopulation.GetNode(*elem_iter)->rGetLocation();

                double grad = (neighbour_morphogen_concentration - current_morphogen_concentration)/norm_2(neighbour_location - current_location);

                if (grad > morphogen_max_grad)
                {
                    morphogen_max_grad = grad;
                }
            }

            // Iterate again, now determining which neighbours have a grad equal to max grad.
            std::vector<unsigned> maximal_gradient_indices; 

            // Iterate over these elements
            for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
                    elem_iter != neighbouring_indices.end();
                    ++elem_iter)
            {
                // Get the cell according to the index
                CellPtr neighbour_cell_iter = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

                // Get the neighbouring concentration
                double neighbour_morphogen_concentration = neighbour_cell_iter->GetCellData()->GetItem("morphogen");

                // Get the neighbour's location
                c_vector<double, DIM> neighbour_location = rCellPopulation.GetNode(*elem_iter)->rGetLocation();

                double grad = (neighbour_morphogen_concentration - current_morphogen_concentration)/norm_2(neighbour_location - current_location);

                if (grad == morphogen_max_grad)
                {
                    maximal_gradient_indices.push_back(*elem_iter);
            }
            }

            if (!maximal_gradient_indices.empty())
            {
                // Now choose an index at random
                unsigned chosen_index = floor(RandomNumberGenerator::Instance()->ranf() * maximal_gradient_indices.size());

                unsigned chosen_node_index = maximal_gradient_indices[chosen_index];
                
                c_vector<double, DIM> chosen_neighbour_location = rCellPopulation.GetNode(chosen_node_index)->rGetLocation();

                c_vector<double, DIM> gradient_direction = rCellPopulation.rGetMesh().GetVectorFromAtoB(current_location, chosen_neighbour_location);

                // Only add the chemotactic force when the change is positive
                if (morphogen_max_grad > 1e-4)
                {
                    // Get the chemotactic strength parameter 
                    double chemotactic_strength = GetChemotacticStrength();

                    // F += chi * gradC/|gradC|
                    c_vector<double, DIM> force = chemotactic_strength * morphogen_max_grad * gradient_direction / norm_2(gradient_direction);
                    rCellPopulation.GetNode(current_index)->AddAppliedForceContribution(force);

                    // Also need to update the fibroblast migration direction, which is stored as cell data.
                    double gradient_angle = atan(gradient_direction[1]/gradient_direction[0]);

                    // Correct for the quadrants
                    if ( (gradient_direction[0] < 0.0) ) // Second or third quadrant quadrant
                    {
                        gradient_angle += M_PI; 
                    }
                    else if ( (gradient_direction[0] > 0.0)&&(gradient_direction[1] < 0.0) ) // Third quadrant
                    {
                        gradient_angle += 2.0*M_PI;
                    }

                    cell_iter->GetCellData()->SetItem("direction", gradient_angle); // Update the fibroblast migration direction
                }
                // else Fc=0
            }
        }
    }
}

template<unsigned DIM>
void WoundBasedChemotacticForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<ChemotacticStrength>"<<  mChemotacticStrength << "</ChemotacticStrength> \n" ;
	*rParamsFile <<  "\t\t\t<NeighbourhoodRadius>" << mNeighbourhoodRadius << "</NeighbourhoodRadius> \n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class WoundBasedChemotacticForce<1>;
template class WoundBasedChemotacticForce<2>;
template class WoundBasedChemotacticForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WoundBasedChemotacticForce)
