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

#include "PdeBasedChemotacticForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "BloodCellProliferativeType.hpp"
#include "ExtracellularMatrixCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
PdeBasedChemotacticForce<DIM>::PdeBasedChemotacticForce()
    : AbstractForce<DIM>(),
    mChemotacticStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
PdeBasedChemotacticForce<DIM>::~PdeBasedChemotacticForce()
{
}

// Method to get the chemotactic strength parameter, i.e. chi
template<unsigned DIM>
double PdeBasedChemotacticForce<DIM>::GetChemotacticStrength()
{
    return mChemotacticStrength;
}

// Method to set the chemotactic strength parameter, chi
template<unsigned DIM>
void PdeBasedChemotacticForce<DIM>::SetChemotacticStrength(double chemotacticStrength)
{
    mChemotacticStrength = chemotacticStrength;
}

template<unsigned DIM>
void PdeBasedChemotacticForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    // Get the current time
    double current_time = SimulationTime::Instance()->GetTime();

    if (current_time > 0.0) // I want to see if we need a time step to happen in order for this to work.
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {

            // Should only look at non-platelet cells and non-ECM cells
            boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType(); // Get the cell type

            if ( (!p_cell_type->IsType<BloodCellProliferativeType>())
                &&(!p_cell_type->IsType<ExtracellularMatrixCellProliferativeType>()) )
            {
                // Get the node index
                unsigned current_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

                Node<DIM>* p_node = rCellPopulation.GetNode(current_index);

                if (!p_node->IsParticle())
                {
                    
                    c_vector<double, DIM> force_direction = zero_vector<double>(DIM);

                    // The grad direction vectors depend on the dimension.
                    switch (DIM)
                    {
                        case 1:
                            force_direction[0] = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
                            break;
                        case 2:
                            force_direction[0] = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
                            force_direction[1] = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
                            break;
                        case 3:
                            force_direction[0] = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
                            force_direction[1] = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
                            force_direction[2] = cell_iter->GetCellData()->GetItem("morphogen_grad_z");                        
                            break;
                    }

                    c_vector<double, DIM> chemotactic_force = mChemotacticStrength * force_direction; // With the grad vector, the chemotactic force is easily calculated

                    p_node->AddAppliedForceContribution(chemotactic_force); // Add the force           
                    
                }
            }
        }
    }
}

template<unsigned DIM>
void PdeBasedChemotacticForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<ChemotacticStrength>"<<  mChemotacticStrength << "</ChemotacticStrength> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PdeBasedChemotacticForce<1>;
template class PdeBasedChemotacticForce<2>;
template class PdeBasedChemotacticForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PdeBasedChemotacticForce)
