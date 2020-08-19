/*
 * LAST MODIFIED: 21/12/2014
 * Platelet cell killer created for wound healing model. Removes platelets that are in contact with 
 * fibroblasts to simulate degradation of the provisional matrix that forms during wound healing.
 *
 * Created on: Dec 21 2014
 * Last modified:
 * 			Author: Axel Almet
 */

#include "PlateletCellKiller.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellProperty.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "BloodCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"

PlateletCellKiller::PlateletCellKiller(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellKiller<2>(pCellPopulation),
    mCutOffRadius(1.5),
    mMeanDeathTime(1.0),
	mVolumeThreshold(0.5*0.25*M_PI)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "PlateletData/", false);
//	mPlateletOutputFile = output_file_handler.OpenOutputFile("results.Platelet");
}

//Method to get mGrowthFactorThreshold
double PlateletCellKiller::GetMeanDeathTime()
{
	return mMeanDeathTime;
}

//Method to set mGrowthFactorThreshold
void PlateletCellKiller::SetMeanDeathTime(double meanDeathTime)
{
	mMeanDeathTime = meanDeathTime;
}

//Method to get mVolumeThreshold
double PlateletCellKiller::GetVolumeThreshold()
{
	return mVolumeThreshold;
}

//Method to set mVolumeThreshold
void PlateletCellKiller::SetVolumeThreshold(double volumeThreshold)
{
	mVolumeThreshold = volumeThreshold;
}

//Method to get mCutOffRadius
double PlateletCellKiller::GetCutOffRadius()
{
	return mCutOffRadius;
}

//Method to set mCutOffRadius
void PlateletCellKiller::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

// Wrapper for GetNodesWithinNeighbourhoodRadius
std::set<unsigned> PlateletCellKiller::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	// Need access to the mesh but can't get to it because the cell killer only owns a
	// pointer to an AbstractCellPopulation
	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

	neighbouring_node_indices = p_tissue->GetNeighbouringNodeIndices(nodeIndex);

    return neighbouring_node_indices;
}

/** Method to determine whether or not cell should be removed
 */
bool PlateletCellKiller::ShouldCellBeRemoved(CellPtr pCell)
{
	bool should_cell_be_removed = false;	// Initialising

    double mean_death_time = GetMeanDeathTime(); // Get the growth factor threshold
	
	// double volume_threshold = GetVolumeThreshold();

	unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(pCell);
	// double current_volume = pCell->GetCellData()->GetItem("volume"); // Get the cell's volume.
	double death_time = pCell->GetCellData()->GetItem("activation_time"); // Get the activation_time
	double current_time =  SimulationTime::Instance()->GetTime(); // Get the current simulation time

	// Only consider platelet cells
	if (pCell->GetCellProliferativeType()->IsType<BloodCellProliferativeType>() )
	{
		std::set<unsigned> neighbours = GetNeighbouringNodeIndices(node_index);

		// Now check whether there are any fibroblast neighbours.
		for (std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
				neighbour_iter != neighbours.end();
				++neighbour_iter)
		{
			if (!this->mpCellPopulation->GetNode(*neighbour_iter)->IsParticle())
			{
				// Only consider fibroblast neighbours
				if (this->mpCellPopulation->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>() )
				{
					// We kill if the fibroblast is sufficiently activated during the wounding, i.e. it has been exposed to a sufficient 
					// amount of growth factor.
					double activation_status = this->mpCellPopulation->GetCellUsingLocationIndex(*neighbour_iter)->GetCellData()->GetItem("activated");
					
					// if ( (activation_status == 1.0)&&(current_volume < volume_threshold)&&(death_time == 0.0) )
					if ( (activation_status == 1.0)&&(death_time == 0.0) )
					{
						double exponential_random_number = RandomNumberGenerator::Instance()->ExponentialRandomDeviate(mean_death_time);
						double death_time = current_time + exponential_random_number;

						// Set the death time
						pCell->GetCellData()->SetItem("activation_time", death_time);
						break;
					}
				}
			}
		}

		// Also check if the platelet should be killed via natural cell death, so to speak
		if ((death_time > 0.0)&&(current_time > death_time))
		{
			should_cell_be_removed = true;
		}

	}

	return should_cell_be_removed;
}

/*
 * Cell Killer that kills platelet cells that are in contact with
 * activated fibroblasts as a result of the wound-induced growth factors.
 */
void PlateletCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
	for (typename AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
	cell_iter != this->mpCellPopulation->End();
	++cell_iter)
    {
		if (cell_iter->GetCellProliferativeType()->IsType<BloodCellProliferativeType>()) // Paranoia
		{
			if(ShouldCellBeRemoved(*cell_iter))
			{
				if (!cell_iter->HasApoptosisBegun())
				{
					// cell_iter->StartApoptosis();
					cell_iter->Kill();
				}
			}
		}
    }

}

void PlateletCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CutOffRadius>" << mCutOffRadius << "</CutOffRadius> \n";
    *rParamsFile << "\t\t\t<MeanDeathTime>" << mMeanDeathTime << "</MeanDeathTime> \n";
	*rParamsFile << "\t\t\t<VolumeThreshold>" << mVolumeThreshold << "</VolumeThreshold \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PlateletCellKiller)
