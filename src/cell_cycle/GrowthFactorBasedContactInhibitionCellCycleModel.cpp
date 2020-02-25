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

#include "GrowthFactorBasedContactInhibitionCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "PlateletCellProliferativeType.hpp"
#include "Debug.hpp"

GrowthFactorBasedContactInhibitionCellCycleModel::GrowthFactorBasedContactInhibitionCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel(),
      mQuiescentVolumeFraction(DOUBLE_UNSET),
      mEquilibriumVolume(DOUBLE_UNSET),
      mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
      mCurrentQuiescentDuration(0.0),
      mGrowthFactorThreshold(0.5)
{
}

GrowthFactorBasedContactInhibitionCellCycleModel::GrowthFactorBasedContactInhibitionCellCycleModel(const GrowthFactorBasedContactInhibitionCellCycleModel& rModel)
    : AbstractSimplePhaseBasedCellCycleModel(rModel),
      mQuiescentVolumeFraction(rModel.mQuiescentVolumeFraction),
      mEquilibriumVolume(rModel.mEquilibriumVolume),
      mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
      mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration),
      mGrowthFactorThreshold(rModel.mGrowthFactorThreshold)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-CycleModel model.
     */
}

void GrowthFactorBasedContactInhibitionCellCycleModel::UpdateCellCyclePhase()
{
    if ((mQuiescentVolumeFraction == DOUBLE_UNSET) || (mEquilibriumVolume == DOUBLE_UNSET))
    {
        EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
    }

    // Get cell volume
    double cell_volume = mpCell->GetCellData()->GetItem("volume");

    // Get the growth factor level and threshold
    double growth_factor_level = mpCell->GetCellData()->GetItem("morphogen");

    // Get the threshold to determine proliferation
    double growth_factor_threshold = GetGrowthFactorThreshold();

    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell volume
        double dt = SimulationTime::Instance()->GetTimeStep();
        double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

        // Pause proliferation if the cell is either sufficiently stressed or is not exposed to a sufficiently
        // high morphogen concentration
        if ( (cell_volume < quiescent_volume)||(growth_factor_level < growth_factor_threshold) )
        {
            // Update the duration of the current period of contact inhibition.
            mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;
            mG1Duration += dt;

            PRINT_VARIABLE(1.0);

            /*
             * This method is usually called within a CellBasedSimulation, after the CellPopulation
             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
             * CellPropertyRegistry::Instance() here when adding the CellLabel, we would be creating
             * a new CellPropertyRegistry. In this case the CellLabel's cell count would be incorrect.
             * We must therefore access the CellLabel via the cell's CellPropertyCollection.
             */
            // boost::shared_ptr<AbstractCellProperty> p_label =
            //     mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
            // mpCell->AddCellProperty(p_label);
        }
        else
        {
            // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
        }
    }


    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    // Differentiated cells, platelet cells, or labelled cells (which represent fixed boundary cells)
    // don't proliferate
    if ((mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        ||(mpCell->GetCellProliferativeType()->IsType<PlateletCellProliferativeType>())
        ||(mpCell->HasCellProperty<CellLabel>()) )
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (time_since_birth < GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
}

void GrowthFactorBasedContactInhibitionCellCycleModel::SetG1Duration()
{
    // We need to override this inherited function to account for fibroblasts and platelets.

    assert(mpCell != nullptr);

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = GetTransitCellG1Duration();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>())
    {
        mG1Duration = GetTransitCellG1Duration();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<PlateletCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else if (mpCell->HasCellProperty<CellLabel>()) // Labelled boundary cells need to be accounted for.
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

AbstractCellCycleModel* GrowthFactorBasedContactInhibitionCellCycleModel::CreateCellCycleModel()
{
    return new GrowthFactorBasedContactInhibitionCellCycleModel(*this);
}

void GrowthFactorBasedContactInhibitionCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double GrowthFactorBasedContactInhibitionCellCycleModel::GetQuiescentVolumeFraction() const
{
    return mQuiescentVolumeFraction;
}

void GrowthFactorBasedContactInhibitionCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double GrowthFactorBasedContactInhibitionCellCycleModel::GetEquilibriumVolume() const
{
    return mEquilibriumVolume;
}

double GrowthFactorBasedContactInhibitionCellCycleModel::GetCurrentQuiescentDuration() const
{
    return mCurrentQuiescentDuration;
}

double GrowthFactorBasedContactInhibitionCellCycleModel::GetCurrentQuiescentOnsetTime() const
{
    return mCurrentQuiescentOnsetTime;
}

double GrowthFactorBasedContactInhibitionCellCycleModel::GetGrowthFactorThreshold() const
{
    return mGrowthFactorThreshold;
}

void GrowthFactorBasedContactInhibitionCellCycleModel::SetGrowthFactorThreshold(double growthFactorThreshold)
{
    mGrowthFactorThreshold = growthFactorThreshold;
}

void GrowthFactorBasedContactInhibitionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";
    *rParamsFile << "\t\t\t<GrowthFactorThreshold>" << mGrowthFactorThreshold << "</GrowthFactorThresholde>\n";


    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(GrowthFactorBasedContactInhibitionCellCycleModel)
