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

#include "PlateletDerivedGrowthFactorCellwiseSourceParabolicPde.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlateletCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<DIM>::PlateletDerivedGrowthFactorCellwiseSourceParabolicPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double duDtCoefficient,
                                                            double diffusionCoefficient,
                                                            double productionCoefficient,
                                                            double degradationCoefficient)
    : mrCellPopulation(rCellPopulation),
      mDuDtCoefficient(duDtCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mProductionCoefficient(productionCoefficient),
      mDegradationCoefficient(degradationCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

// LCOV_EXCL_START
template<unsigned DIM>
double PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
double PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    double source_coefficient = 0.0;

    unsigned node_index = rNode.GetIndex();
    
    if (mrCellPopulation.IsPdeNodeAssociatedWithNonApoptoticCell(node_index))
    {

        // Get the cell proliferative type
        CellPtr p_cell = mrCellPopulation.GetCellUsingLocationIndex(node_index); // Get the cell
        boost::shared_ptr<AbstractCellProperty> p_cell_type = p_cell->GetCellProliferativeType(); // Get the mutation state

        if (p_cell_type->IsType<PlateletCellProliferativeType>())
        {
            // source_coefficient = mProductionCoefficient*u/(1.0 + mProductionCoefficient*u) - mDegradationCoefficient*u;
            source_coefficient = mProductionCoefficient*u - mDegradationCoefficient*u;
        }
        else
        {
            source_coefficient = -mDegradationCoefficient*u;
        }
    }

    // The source term is C*u
    return source_coefficient;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

// Explicit instantiation
template class PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<1>;
template class PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<2>;
template class PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlateletDerivedGrowthFactorCellwiseSourceParabolicPde)
