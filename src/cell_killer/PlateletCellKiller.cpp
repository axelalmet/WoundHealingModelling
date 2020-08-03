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
    mCellsRemovedByFibroblasts(0),
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

    // if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	// {
	// Need access to the mesh but can't get to it because the cell killer only owns a
	// pointer to an AbstractCellPopulation
	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

	// //Update cell population
	// p_tissue->Update();

	// double radius = GetCutOffRadius;


	neighbouring_node_indices = p_tissue->GetNeighbouringNodeIndices(nodeIndex);

	// neighbouring_node_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(nodeIndex, radius);
	// }

    return neighbouring_node_indices;
}

/** Method to determine whether or not cell should be removed
 */
bool PlateletCellKiller::ShouldCellBeRemoved(unsigned nodeIndex)
{
	bool should_cell_be_removed = false;	// Initialising

    double mean_death_time = GetMeanDeathTime(); // Get the growth factor threshold
	
	double volume_threshold = GetVolumeThreshold();

	CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(nodeIndex);
	double current_volume = p_cell->GetCellData()->GetItem("volume"); // Get the cell's volume.
	double death_time = p_cell->GetCellData()->GetItem("activation time"); // Get the activation time
	double current_time =  SimulationTime::Instance()->GetTime(); // Get the current simulation time

	// Only consider platelet cells
	if (p_cell->GetCellProliferativeType()->IsType<BloodCellProliferativeType>() )
	{
		std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

		// Now check whether there are any fibroblast neighbours.
		for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
				neighbour_iter != neighbours.end();
				++neighbour_iter)
		{
			// Only consider fibroblast neighbours
			if (this->mpCellPopulation->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>() )
			{
				// We kill if the fibroblast is sufficiently activated during the wounding, i.e. it has been exposed to a sufficient 
				// amount of growth factor.
				double activation_status = this->mpCellPopulation->GetCellUsingLocationIndex(*neighbour_iter)->GetCellData()->GetItem("activated");
				
				if ( (activation_status == 1.0)&&(current_volume < volume_threshold)&&(death_time == 0.0) )
				{
					double random_uniform_number = RandomNumberGenerator::Instance()->ranf();
					double death_time = current_time - mean_death_time * log(random_uniform_number);

					// Set the death time
					p_cell->GetCellData()->SetItem("activation time", death_time);

					break;
				}
			}
		}

		// Also check if the platelet should be killed via natural cell death, so to speak
		if ((death_time > 0.0)&&(current_time > death_time))
		{
			should_cell_be_removed = true;
		}

	}
	// }

	return should_cell_be_removed;
}

/** A method to return a vector that indicates which cells should be killed by the fibroblasts
 */
std::vector<c_vector<unsigned,2> > PlateletCellKiller::RemoveByFibroblasts()
{

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    // if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
    // {
	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

	c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
			cell_iter != this->mpCellPopulation->End();
			++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();

		// Initialise
		individual_node_information[0] = node_index;
		individual_node_information[1] = 0;

		// Examine each epithelial node to see if it should be removed by Platelet and then if it
		// should be removed by compression-driven apoptosis
		if (cell_iter->GetCellProliferativeType()->IsType<BloodCellProliferativeType>())
		{
			// Determining whether to remove this cell by fibroblasts

			if (this->mpCellPopulation->IsCellAttachedToLocationIndex(node_index))
			{
				if(this->ShouldCellBeRemoved(node_index))
				{
					individual_node_information[1] = 1;
				}
			}
		}

		cells_to_remove.push_back(individual_node_information);
	}
    // }

	return cells_to_remove;
}


/*
 * Cell Killer that kills platelet cells that are in contact with
 * activated fibroblasts as a result of the wound-induced growth factors.
 */
void PlateletCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
    // if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	// {
	// NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

	// Get the information at this timestep for each node index that says whether to remove by fibroblasts
	std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByFibroblasts();

	// Keep a record of how many cells have been removed at this timestep
	this->SetNumberCellsRemoved(cells_to_remove);
	this->SetLocationsOfCellsRemovedByFibroblasts(cells_to_remove);

	// Need to avoid trying to kill any cells twice (i.e. both by fibroblasts)
	// Loop over these vectors individually and kill any cells that they tell you to

	for (unsigned i=0; i<cells_to_remove.size(); i++)
	{
		if (cells_to_remove[i][1] == 1)
		{
			// Get cell associated to this node
			CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(cells_to_remove[i][0]);

			// if (!p_cell->HasApoptosisBegun())
			// {
				// p_cell->StartApoptosis();
			p_cell->Kill();
			// }
		}
	}

}

void PlateletCellKiller::SetNumberCellsRemoved(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	unsigned num_removed_by_fibroblasts = 0;

    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_fibroblasts+=1;
    	}
    }

    mCellsRemovedByFibroblasts += num_removed_by_fibroblasts;
}

unsigned PlateletCellKiller::GetNumberCellsRemoved()
{
	return mCellsRemovedByFibroblasts;
}

void PlateletCellKiller::SetLocationsOfCellsRemovedByFibroblasts(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
 
    // if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	// {
	// NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);
	double x_location, y_location;
	c_vector<double, 3> time_and_location;

	// Need to use the node indices to store the locations of where cells are removed
	for (unsigned i=0; i<cellsRemoved.size(); i++)
	{
		if (cellsRemoved[i][1] == 1)		// This cell has been removed by Platelet
		{
			time_and_location[0] = SimulationTime::Instance()->GetTime();

			unsigned node_index = cellsRemoved[i][0];

			CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(node_index);
			x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
			y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];

			time_and_location[1] = x_location;
			time_and_location[2] = y_location;

			mLocationsOfPlateletCells.push_back(time_and_location);
		}
	}
	// }
}

std::vector<c_vector<double,3> > PlateletCellKiller::GetLocationsOfCellsRemovedByFibroblasts()
{
	return mLocationsOfPlateletCells;
}

void PlateletCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByFibroblasts>" << mCellsRemovedByFibroblasts << "</CellsRemovedByFibroblasts> \n";
    *rParamsFile << "\t\t\t<CutOffRadius>" << mCutOffRadius << "</CutOffRadius> \n";
    *rParamsFile << "\t\t\t<MeanDeathTime>" << mMeanDeathTime << "</MeanDeathTime> \n";
	*rParamsFile << "\t\t\t<VolumeThreshold>" << mVolumeThreshold << "</VolumeThreshold \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PlateletCellKiller)
