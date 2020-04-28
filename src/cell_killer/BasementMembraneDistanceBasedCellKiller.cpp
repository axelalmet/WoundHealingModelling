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
#include "BasementMembraneDistanceBasedCellKiller.hpp"
#include "NodeBasedCellPOpulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "PlateletCellProliferativeType.hpp"

BasementMembraneDistanceBasedCellKiller::BasementMembraneDistanceBasedCellKiller(AbstractCellPopulation<2>* pCellPopulation, double radius, double maxHeight)
    : AbstractCellKiller<2>(pCellPopulation),
      mRadius(radius),
      mMaxHeight(maxHeight)
{
}

double BasementMembraneDistanceBasedCellKiller::GetRadius() const
{
    return mRadius;
}

double BasementMembraneDistanceBasedCellKiller::GetMaxHeight() const
{
    return mMaxHeight;
}

void BasementMembraneDistanceBasedCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{

    // This method assumes that we're using a NodeBasedCellPopulation
    assert(dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation));
    NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

    for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {

        // Only consider differentiated cells (these are assumed to be epidermal cells)
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        if (p_cell_type->IsType<DifferentiatedCellProliferativeType>())
        {
            // Get the node index
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

            double y = this->mpCellPopulation->GetNode(node_index)->rGetLocation()[1];
            if (y > mMaxHeight)
            {
                cell_iter->Kill();
            }
            else
            {
                unsigned num_fibroblast_neighbours = 0;
                // Now get the neighbours within mRadius
                std::set<unsigned> neighbour_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(node_index, mRadius);
                
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*iter);
                    boost::shared_ptr<AbstractCellProperty> p_neighbour_cell_type = p_cell->GetCellProliferativeType();

                    // If there are ANY fibroblast neighbours, we immediately know that it's attached to the BM and can stop the iterations.
                    if (p_neighbour_cell_type->IsType<FibroblastCellProliferativeType>())
                    {
                        num_fibroblast_neighbours += 1;
                        break;
                    }
                }

                if (num_fibroblast_neighbours == 0)
                {
                    cell_iter->Kill();
                }
            }
        }
        else if (p_cell_type->IsType<PlateletCellProliferativeType>()) // Slough off platelet cells if they're too high up
        {
            // Get the node index
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

            double y = this->mpCellPopulation->GetNode(node_index)->rGetLocation()[1];

            if (y > mMaxHeight)
            {
                cell_iter->Kill();
            }
        }
    }
}

void BasementMembraneDistanceBasedCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<mRadius>" << mRadius << "</mRadius>\n";
    *rParamsFile << "\t\t\t<mMaxHeight>" << mMaxHeight << "</mMaxHeight>\n";

    // Call method on direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BasementMembraneDistanceBasedCellKiller)
