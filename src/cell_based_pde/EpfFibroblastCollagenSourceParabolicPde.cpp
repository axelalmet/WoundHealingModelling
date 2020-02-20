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

#include "EpfFibroblastCollagenSourceParabolicPde.hpp"
#include "EpfFibroblastCellMutationState.hpp"

template<unsigned DIM>
EpfFibroblastCollagenSourceParabolicPde<DIM>::EpfFibroblastCollagenSourceParabolicPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
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
const AbstractCellPopulation<DIM,DIM>& EpfFibroblastCollagenSourceParabolicPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double EpfFibroblastCollagenSourceParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

// LCOV_EXCL_START
template<unsigned DIM>
double EpfFibroblastCollagenSourceParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
double EpfFibroblastCollagenSourceParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    double source_coefficient = 0.0;

    if (mrCellPopulation.IsPdeNodeAssociatedWithNonApoptoticCell(rNode.GetIndex()))
    {
      // Get the cell mutation state.
      CellPtr p_cell = mrCellPopulation.GetCellUsingLocationIndex(rNode.GetIndex()); // Get the cell
      boost::shared_ptr<AbstractCellProperty> p_mutation_state = p_cell->GetMutationState(); // Get the mutation state
      
      if (p_mutation_state->IsType<EpfFibroblastCellMutationState>())
      {
        source_coefficient = mProductionCoefficient - mDegradationCoefficient;
      }
      else
      {
        source_coefficient = -mDegradationCoefficient;
      }
      
      
    }

    // The source term is C*u
    return source_coefficient*u;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> EpfFibroblastCollagenSourceParabolicPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

// Explicit instantiation
template class EpfFibroblastCollagenSourceParabolicPde<1>;
template class EpfFibroblastCollagenSourceParabolicPde<2>;
template class EpfFibroblastCollagenSourceParabolicPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EpfFibroblastCollagenSourceParabolicPde)
