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

#include "BasementMembraneAttachmentTrackingModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
BasementMembraneAttachmentTrackingModifier<DIM>::BasementMembraneAttachmentTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mNeighbourhoodRadius(1.5)
{
}

template<unsigned DIM>
BasementMembraneAttachmentTrackingModifier<DIM>::~BasementMembraneAttachmentTrackingModifier()
{
}

template<unsigned DIM>
double BasementMembraneAttachmentTrackingModifier<DIM>::GetNeighbourhoodRadius()
{
    return mNeighbourhoodRadius;
}

template<unsigned DIM>
void BasementMembraneAttachmentTrackingModifier<DIM>::SetNeighbourhoodRadius(double neighbourhoodRadius)
{
    mNeighbourhoodRadius = neighbourhoodRadius;
}

template<unsigned DIM>
void BasementMembraneAttachmentTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void BasementMembraneAttachmentTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void BasementMembraneAttachmentTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        // Only consider stem cells
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        if (!p_cell_type->template IsType<FibroblastCellProliferativeType>())
        {
            // Get the node index
            unsigned node_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);
                
            // Initialise number of fibroblast neighbours
            unsigned num_fibroblast_neighbours = 0;

            // Get the set of neighbouring location indices within a neighbourhood radius
            std::set<unsigned> neighbour_indices = p_cell_population->GetNodesWithinNeighbourhoodRadius(node_index, neighbourhood_radius);

            // Count the number of fibroblast neighbours
            if (!neighbour_indices.empty())
            {
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    boost::shared_ptr<AbstractCellProperty> p_neighbour_cell_type = p_cell->GetCellProliferativeType();

                    // If there are ANY fibroblast neighbours, we immediately know that it's attached to the BM and can stop the iterations.
                    if (p_neighbour_cell_type->template IsType<FibroblastCellProliferativeType>())
                    {
                        num_fibroblast_neighbours += 1;
                        cell_iter->GetCellData()->SetItem("attachment", 1.0);
                        break;
                    }
                }
            }

            // Else the cell has been detached from the basement membrane.
            if (num_fibroblast_neighbours == 0)
            {
                // If this cell has no fibroblast neighbours, store 0.0 for the cell data
                cell_iter->GetCellData()->SetItem("attachment", 0.0);
            }

        }
        else
        {
            cell_iter->GetCellData()->SetItem("attachment", -1.0);
        }
        

    }
}

template<unsigned DIM>
void BasementMembraneAttachmentTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class BasementMembraneAttachmentTrackingModifier<1>;
template class BasementMembraneAttachmentTrackingModifier<2>;
template class BasementMembraneAttachmentTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BasementMembraneAttachmentTrackingModifier)
