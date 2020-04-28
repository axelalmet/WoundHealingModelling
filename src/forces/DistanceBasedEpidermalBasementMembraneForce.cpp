#include "DistanceBasedEpidermalBasementMembraneForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "AbstractCellProperty.hpp"
#include "Debug.hpp"
#include "Exception.hpp"

/*
 * This is a force class to simulate the basement membrane,
 * based off that originally proposed by Dunn et al. (2012).
 */


/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
DistanceBasedEpidermalBasementMembraneForce::DistanceBasedEpidermalBasementMembraneForce()
:  AbstractForce<2>(),
   mBasementMembraneParameter(DOUBLE_UNSET),
   mCutOffRadius(1.5)  
{
	// Sets up output file
	//	OutputFileHandler output_file_handler("CurvatureData/", false);
	//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

DistanceBasedEpidermalBasementMembraneForce::~DistanceBasedEpidermalBasementMembraneForce()
{
	//    mMeinekeOutputFile->close();
}

void DistanceBasedEpidermalBasementMembraneForce::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double DistanceBasedEpidermalBasementMembraneForce::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}

double DistanceBasedEpidermalBasementMembraneForce::GetCutOffRadius()
{
	return mCutOffRadius;
}

void DistanceBasedEpidermalBasementMembraneForce::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

/*
 * Method to find all the Epidermal cells that make up the monolayer
 */
std::vector<unsigned> DistanceBasedEpidermalBasementMembraneForce::GetEpidermalIndices(AbstractCellPopulation<2>& rCellPopulation)
{
	//Create the vector of Epidermal cell indices
	std::vector<unsigned> epidermal_indices;

	std::vector<std::pair<double, unsigned> > epidermal_indices_and_x_coordinates; //Define vector of pairs, so that we may sort by the x-coordinate

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
				cell_iter != p_tissue->End();
				++cell_iter)
		{
			//Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			if(p_type->IsType<StemCellProliferativeType>()) //If we have an Epidermal cell
			{
				unsigned node_index = p_tissue->GetLocationIndexUsingCell(*cell_iter);
				double x_coordinate = rCellPopulation.GetNode(node_index)->rGetLocation()[0];

				// Create the pair of index and x-coordinate
				std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(x_coordinate, node_index);
				epidermal_indices_and_x_coordinates.push_back(x_coordinate_and_index);
			}
		}
	}

	//Sort indices by the x-coordinate
	std::sort(epidermal_indices_and_x_coordinates.begin(), epidermal_indices_and_x_coordinates.end());

	// Push back the indices into the vector
	for (unsigned i = 0; i < epidermal_indices_and_x_coordinates.size(); i++)
	{
		unsigned epidermal_index = epidermal_indices_and_x_coordinates[i].second; //Get node index

		epidermal_indices.push_back(epidermal_index);
	}

	return epidermal_indices;
}

/*
 * Method to get the neighbouring nodes that are epidermal
 */
std::vector<unsigned> DistanceBasedEpidermalBasementMembraneForce::GetNeighbouringFibroblastIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned epidermalIndex)
{
	// Create a set of neighbouring node indices
	std::vector<unsigned> neighbouring_fibroblast_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		//Get cut off radius for defining neighbourhood
		double radius = GetCutOffRadius();

		// Find the indices of the elements owned by this node
		// std::set<unsigned> neighbouring_indices = p_tissue->GetNeighbouringNodeIndices(epidermalIndex);
		std::set<unsigned> neighbouring_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(epidermalIndex, radius);

		// Iterate over these elements
		for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
				elem_iter != neighbouring_indices.end();
				++elem_iter)
		{
			//Get the cell according to the index
			CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

			//Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			//if the cell is not differentiated and thus an Epidermal cell, we add it to the vector
			if(p_type->IsType<FibroblastCellProliferativeType>())
			{
				neighbouring_fibroblast_indices.push_back(*elem_iter);
			}

		}
	}

	return neighbouring_fibroblast_indices;
}


// c_vector<double, 2> DistanceBasedEpidermalBasementMembraneForce::CalculateForceDirection(AbstractCellPopulation<2>& rCellPopulation, 
//                                                                     std::vector<unsigned> fibroblastIndices,
//                                                                     unsigned epidermalIndex)
// {
// 	// Get the location of the epidermal index
// 	c_vector<double, 2> epidermal_location = rCellPopulation.GetNode(epidermalIndex)->rGetLocation();

// 	// Get the number of fibroblast indices
// 	unsigned num_fibroblast_indices = fibroblastIndices.size(); 

// 	// Initialise the force direction vector
// 	c_vector<double, 2> force_direction;

// 	for (unsigned i = 0; i < num_fibroblast_indices; i++)
// 	{
// 		unsigned fibroblast_index = fibroblastIndices[i];
// 		c_vector<double, 2> fibroblast_location = rCellPopulation.GetNode(fibroblast_index)->rGetLocation();

// 		// Get the vector between the fibroblast and the epidermal cell
// 		c_vector<double, 2> epidermal_to_fibroblast = rCellPopulation.rGetMesh().GetVectorFromAtoB(epidermal_location, fibroblast_location);

// 		force_direction += epidermal_to_fibroblast / num_fibroblast_indices; 
// 	}

// 	force_direction /= norm_2(force_direction);

// 	// If there are more than 2 fibroblast neighbours, it's likely the epidermal cell is too far into the dermis
// 	if (num_fibroblast_indices > 2)
// 	{
// 		force_direction *= -1.0;
// 	}

// 	return force_direction;
// }

//Method to calculate the force due to the basement membrane on an Epidermal cell
c_vector<double, 2> DistanceBasedEpidermalBasementMembraneForce::CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, 
																						std::vector<unsigned> epidermalIndices,
																						unsigned nodeIndex)
{
	// Get the basement membrane stiffness
	double basement_membrane_parameter = GetBasementMembraneParameter(); 

	//Initialise the force vector
	c_vector<double, 2> force_due_to_basement_membrane;
	force_due_to_basement_membrane[0] = 0.0;
	force_due_to_basement_membrane[1] = 0.0;

	// Get the considered node's location
	c_vector<double, 2> epidermal_location = rCellPopulation.GetNode(nodeIndex)->rGetLocation();

	// Get the epidermal node radius
	double epidermal_radius = rCellPopulation.GetNode(nodeIndex)->GetRadius();

	// Get the fibroblast-to-centre vector
	std::vector<unsigned> closest_fibroblast_indices = GetNeighbouringFibroblastIndices(rCellPopulation, nodeIndex);

	 // If there's no fibroblast neighbours, the cell has detached from the membrane.
	for (unsigned i = 0; i < closest_fibroblast_indices.size(); i++)
	{
		unsigned fibroblast_index = closest_fibroblast_indices[i];
		c_vector<double, 2> fibroblast_location = rCellPopulation.GetNode(fibroblast_index)->rGetLocation();
		// Get the fibroblast node radius too
		double fibroblast_radius = rCellPopulation.GetNode(fibroblast_index)->GetRadius();

		// Get the vector from the epidermal cell to the fibroblast
		c_vector<double, 2> epidermal_to_fibroblast = rCellPopulation.rGetMesh().GetVectorFromAtoB(epidermal_location, fibroblast_location);
		double distance_to_bm = norm_2(epidermal_to_fibroblast);

		// Only add a non-zero force if the overlap is positive, i.e. it's not too close.
		if (distance_to_bm > epidermal_radius + fibroblast_radius)
		{
			c_vector<double, 2> force_direction = epidermal_to_fibroblast / distance_to_bm;

			force_due_to_basement_membrane += basement_membrane_parameter/(distance_to_bm + 0.01) * force_direction;
		}

		
	}

	if (closest_fibroblast_indices.size() > 2) // If there are too many neighbours, it's too far into the dermis
	{
		force_due_to_basement_membrane *= -1.0;
	}

	// c_vector<double, 2> force_due_to_basement_membrane = zero_vector<double>(2);
	return force_due_to_basement_membrane;
}

//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void DistanceBasedEpidermalBasementMembraneForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{

	//Get the epidermal indices
	std::vector<unsigned> epidermal_indices = GetEpidermalIndices(rCellPopulation);

	// Iterate over each node and apply the force
	for (unsigned i = 0; i < epidermal_indices.size(); i++)
	{
		unsigned epidermal_index = epidermal_indices[i];

		// Calculate the basement membrane force
		c_vector<double, 2> force_on_node = CalculateForceDueToBasementMembrane(rCellPopulation, epidermal_indices, epidermal_index);

		// // Apply the force
		rCellPopulation.GetNode(epidermal_index)->AddAppliedForceContribution(force_on_node);
	}
}

void DistanceBasedEpidermalBasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DistanceBasedEpidermalBasementMembraneForce)
