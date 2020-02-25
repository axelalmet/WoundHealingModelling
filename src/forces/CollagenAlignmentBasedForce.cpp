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

#include "CollagenAlignmentBasedForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FibroblastCellProliferativeType.hpp"

template<unsigned DIM>
CollagenAlignmentBasedForce<DIM>::CollagenAlignmentBasedForce()
    : AbstractForce<DIM>(),
    mForceMagnitude(DOUBLE_UNSET),
    mNeighbourhoodRadius(DOUBLE_UNSET)
{
}

template<unsigned DIM>
CollagenAlignmentBasedForce<DIM>::~CollagenAlignmentBasedForce()
{
}

// Method to get the chemotactic strength parameter, i.e. chi
template<unsigned DIM>
double CollagenAlignmentBasedForce<DIM>::GetForceMagnitude()
{
    return mForceMagnitude;
}

// Method to set the chemotactic strength parameter, chi
template<unsigned DIM>
void CollagenAlignmentBasedForce<DIM>::SetForceMagnitude(double ForceMagnitude)
{
    mForceMagnitude = ForceMagnitude;
}

template<unsigned DIM>
void CollagenAlignmentBasedForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    // For now, this will only apply to node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_tissue = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        // Only look at fibroblasts
        if (cell_iter->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>())
        {
            // Get the node index
            unsigned current_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Get the current orientation of the cell
            double orientation = cell_iter->GetCellData()->GetItem("orientation");

            // Determine the force direction using the cell orientation
            c_vector<double, DIM> force;

            // NOTE: This really only applies in 2D, despite the DIM template.
            force[0] = cos(orientation);
            force[1] = sin(orientation);

            double force_magnitude = GetForceMagnitude();
            
            // Multiply by the force magnitude
            force *= force_magnitude;

            // Add the force contribution
            rCellPopulation.GetNode(current_index)->AddAppliedForceContribution(force);
        }

    }
}

template<unsigned DIM>
void CollagenAlignmentBasedForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<ForceMagnitude>"<<  mForceMagnitude << "</ForceMagnitude> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CollagenAlignmentBasedForce<1>;
template class CollagenAlignmentBasedForce<2>;
template class CollagenAlignmentBasedForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CollagenAlignmentBasedForce)
