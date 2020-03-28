	#include "EpidermalBasementMembraneForce.hpp"
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
EpidermalBasementMembraneForce::EpidermalBasementMembraneForce()
:  AbstractForce<2>(),
   mBasementMembraneParameter(DOUBLE_UNSET),
   mTargetCurvature(DOUBLE_UNSET),
   mCutOffRadius(1.5)  
{
	// Sets up output file
	//	OutputFileHandler output_file_handler("CurvatureData/", false);
	//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

EpidermalBasementMembraneForce::~EpidermalBasementMembraneForce()
{
	//    mMeinekeOutputFile->close();
}

void EpidermalBasementMembraneForce::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double EpidermalBasementMembraneForce::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}


void EpidermalBasementMembraneForce::SetTargetCurvature(double targetCurvature)
{
	mTargetCurvature = targetCurvature;
}


double EpidermalBasementMembraneForce::GetTargetCurvature()
{
	return mTargetCurvature;
}

double EpidermalBasementMembraneForce::GetCutOffRadius()
{
	return mCutOffRadius;
}

void EpidermalBasementMembraneForce::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

c_vector<double,2> EpidermalBasementMembraneForce::GetEpidermisHeightExtremes(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> epidermalWidthExtremes)
{
	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

	// Create a vector to store the y-coordinates of the lowest point of the Epidermis base and the highest point of the
	// Epidermis orifice
	c_vector<double,2> height_extremes;

	double max_height = 0.0;
	double min_height = 0.0; // This will make sense in a bit

	// We will need the width extremes, actually
	double max_width = epidermalWidthExtremes[0];

	// To account for injury, we will consider the dermal fibroblast positions
	// instead. The maximum height is found by 
	for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{

		// Need these to not be labelled cells
		if ( (cell_iter->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>()) )
		{
			Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

			double x = p_node->rGetLocation()[0];
			double y = p_node->rGetLocation()[1];

			if (y > max_height)
			{
				max_height = y;
			}

			// We determine the minimum height by looking for the dermal cell
			// by considering the cells around x = L/2 and taking the maximal
			// dermal height of cells at x = 1/2. The ball of 0.25 is to just
			// give us a bit of a buffer.
			if ( (x > 0.5*max_width - 0.25)&&(x < 0.5*max_width + 0.25) )
			{
				if (y > min_height)
				{
					min_height = y;
				}
			}
		}
	}

	height_extremes[0] = max_height;
	height_extremes[1] = min_height;

	return height_extremes;
}

c_vector<double,2> EpidermalBasementMembraneForce::GetEpidermisWidthExtremes(AbstractCellPopulation<2>& rCellPopulation)
{
	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

	// Create a vector to store the y-coordinates of the lowest point of the Epidermis base and the highest point of the
	// Epidermis orifice
	c_vector<double,2> width_extremes;

	double max_width = 0.0;
	double min_width = DBL_MAX;

	double current_width_coordinate;

	// We iterate over all cells in the tissue, and deal only with those that are Epidermal cells
	for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{

		// Need these to not be labelled cells
		if ( (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) )
		{
			Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

			current_width_coordinate = p_node->rGetLocation()[0];

			if (current_width_coordinate > max_width)
			{
				max_width = current_width_coordinate;
			}
			else if (current_width_coordinate < min_width)
			{
				min_width = current_width_coordinate;
			}
		}
	}

	width_extremes[0] = max_width;
	width_extremes[1] = min_width;

	return width_extremes;
}

/*
 * Method to calculate the tangent line at a point, based on approximating the epithelium positions
 * by a cosine curve.
 */
c_vector<double, 2> EpidermalBasementMembraneForce::GetCosineBasedTangentVector(AbstractCellPopulation<2>& rCellPopulation, 
																				c_vector<double, 2> epidermalHeightExtremes,
                                                            					c_vector<double, 2> epidermalWidthExtremes,
																				c_vector<double, 2> point)
{
	// Get the minimal and maximal x-coordinates of epithelium
	double max_width = epidermalWidthExtremes[0];
	double min_width = epidermalWidthExtremes[1];

	// Get the minimal and maximal y-coordinates of epithelium
	double max_height = epidermalHeightExtremes[0];
	double min_height = epidermalHeightExtremes[1];

	// Define width of epidermis
	double epidermis_width = max_width - min_width;

	// Define epidermis height
	double epidermis_height = max_height - min_height;

	// Tangent defined by derivative of approximated cosine function
	c_vector<double, 2> tangent_vector;
	double x = point[0];

	// y-component is derivative of cos function approximation of epithelium, with appropriate scalings.
	tangent_vector[0] = 1.0;
	tangent_vector[1] = -(M_PI/epidermis_width)*epidermis_height*sin(2.0*M_PI/epidermis_width*(x - min_width));

	tangent_vector /= norm_2(tangent_vector); // Convert to unit vector

	// Return a double-length vector---this works best for tracking neighbours
	tangent_vector *= 2.0;

	return tangent_vector;
}

/*
 * Return vector of Epidermal indices that are close to the considered Epidermal node,
 * but based on an approximated cosine approximation
 */
std::vector<unsigned> EpidermalBasementMembraneForce::GetClosestNeighboursBasedOnCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, 
																									std::vector<unsigned> epidermalIndices,
																									c_vector<double, 2> epidermalHeightExtremes,
                                                            										c_vector<double, 2> epidermalWidthExtremes,
																									unsigned epidermalIndex,
																									double leftOrRight)
{
	// Initialise vector
	std::vector<unsigned> closest_neighbours;

	// Get the location of the Epidermal node
	c_vector<double, 2> epidermal_location = rCellPopulation.GetNode(epidermalIndex)->rGetLocation();

	// Get the tangent vector based on the cosine approximation
	c_vector<double, 2> tangent_vector = GetCosineBasedTangentVector(rCellPopulation, epidermalHeightExtremes, epidermalWidthExtremes, epidermal_location);

	// We multiply the tangent vector by leftOrRight to account for whether we're looking at 
	// 'left' or 'right' neighbours
	tangent_vector *= leftOrRight;

	// Get the cut-off radius to determine the nearest neighbours
	double neighbourhood_radius = GetCutOffRadius();

	// Sweep through the indices
	for (unsigned i = 0; i < epidermalIndices.size(); i++)
	{
		unsigned neighbour_index = epidermalIndices[i];

		if (neighbour_index != epidermalIndex)
		{
			c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

			if ((norm_2(neighbour_location - epidermal_location) < neighbourhood_radius + 0.5) // Is the neighbour reasonably close
					&&(inner_prod(neighbour_location - epidermal_location, tangent_vector) > 0.0 ) ) // Is the neighbour pointing in the same direction as the tangent vector
			{
				closest_neighbours.push_back(neighbour_index);
			}

		}
	}

	return closest_neighbours;
}

/*
 * Return the nearest neighbour based on vector projections from the tangent vector
 * at a point along the cosine approximation of the epithelium.
 */
unsigned EpidermalBasementMembraneForce::GetNearestNeighbourAlongCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, 
																					std::vector<unsigned> epidermalIndices,
																					c_vector<double, 2> epidermalHeightExtremes,
																					c_vector<double, 2> epidermalWidthExtremes,
																					unsigned epidermalIndex,
																					double leftOrRight)
{
	double min_scalar_projection = DBL_MAX;
	double min_projection_index = 0;

	// Get the closest neighbours, based on the cosine approximation
	std::vector<unsigned> closest_neighbours = GetClosestNeighboursBasedOnCosineApproximation(rCellPopulation, 
																								epidermalIndices,
																								epidermalHeightExtremes,
																								epidermalWidthExtremes,
																								epidermalIndex,
																								leftOrRight);

	// Get the location of the considered Epidermal node
	c_vector<double, 2> epidermal_location = rCellPopulation.GetNode(epidermalIndex)->rGetLocation();

	// Calculate the tangent vector at the Epidermal node
	c_vector<double, 2> tangent_vector = GetCosineBasedTangentVector(rCellPopulation, epidermalHeightExtremes, epidermalWidthExtremes, epidermal_location);

	for (unsigned i = 0; i < closest_neighbours.size(); i++)
	{
		// Get the location of the neighbour
		unsigned neighbour_index = closest_neighbours[i];
		c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

		// Get the scalar projection of the relative position vector onto the tangent vector
		double scalar_projection = inner_prod(neighbour_location - epidermal_location, tangent_vector)/norm_2(tangent_vector);

		if (scalar_projection < min_scalar_projection)
		{
			min_scalar_projection = scalar_projection;
			min_projection_index = neighbour_index;
		}
	}

	return min_projection_index;
}

/*
 * Return the nearest neighbour based on vector projections from the tangent vector
 * at a point along the cosine approximation of the epithelium.
 */
unsigned EpidermalBasementMembraneForce::GetNearestFibroblastNeighbour(AbstractCellPopulation<2>& rCellPopulation, unsigned epidermalIndex)
{
	double min_fibroblast_distance = DBL_MAX;
	unsigned closest_fibroblast_index = 0;

	// Get the location of the considered Epidermal node
	c_vector<double, 2> epidermal_location = rCellPopulation.GetNode(epidermalIndex)->rGetLocation();

	// This method really only applies to NodeBasedCellPopulations
	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		//Get cut off radius for defining neighbourhood
		double radius = GetCutOffRadius();

		// Find the indices of the elements owned by this node
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

			// Only consider fibroblasts
			if(p_type->IsType<FibroblastCellProliferativeType>())
			{
				
				c_vector<double, 2> fibroblast_location = rCellPopulation.GetNode(*elem_iter)->rGetLocation();

				// By the end of this loop, we should end up with the closest fibroblast neighbour
				if (norm_2(fibroblast_location - epidermal_location) < min_fibroblast_distance)
				{
					min_fibroblast_distance = norm_2(fibroblast_location - epidermal_location);
					closest_fibroblast_index = rCellPopulation.GetNode(*elem_iter)->GetIndex(); // Roundabout way of storing the index, as we can't convert the iterator to an unsigned
				}
			}

		}
	}

	return closest_fibroblast_index;
}


/*
 * Function to find the curvature along three points, using the method previously described by SJD
 * Method has been adjusted to account for periodic meshes etc, i.e. heavy use of GetVectorFromAtoB
 */
double EpidermalBasementMembraneForce::FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
		c_vector<double, 2> leftPoint,
		c_vector<double, 2> centrePoint,
		c_vector<double, 2> rightPoint)
{

	//Get the relevant vectors (all possible differences)
	c_vector<double, 2> left_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftPoint, centrePoint);
	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centrePoint, rightPoint);
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftPoint, rightPoint);

	// Firstly find the parametric intervals
	double left_s = pow(pow(left_to_centre[0],2.0) + pow(left_to_centre[1],2.0), 0.5);
	double right_s = pow(pow(centre_to_right[0],2.0) + pow(centre_to_right[1],2.0), 0.5);
	//	(left_s, right_s);

	double sum_intervals = left_s + right_s;

	//Calculate finite difference of first derivatives
	double x_prime = (left_to_right[0])/sum_intervals;
	double y_prime = (left_to_right[1])/sum_intervals;

	//Calculate finite difference of second derivatives
	double x_double_prime = 2*(left_s*centre_to_right[0] - right_s*left_to_centre[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*centre_to_right[1] - right_s*left_to_centre[1])/(left_s*right_s*sum_intervals);

	//Calculate curvature using formula
	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2.0) + pow(y_prime,2.0)), 1.5);

	return curvature;
}

/*
 * Method to find all the Epidermal cells that make up the monolayer
 */
std::vector<unsigned> EpidermalBasementMembraneForce::GetEpidermalIndices(AbstractCellPopulation<2>& rCellPopulation)
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
std::vector<unsigned> EpidermalBasementMembraneForce::GetNeighbouringEpidermalIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::vector<unsigned> neighbouring_epidermal_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		//Get cut off radius for defining neighbourhood
		double radius = GetCutOffRadius();

		// Find the indices of the elements owned by this node
		std::set<unsigned> neighbouring_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(nodeIndex, radius);

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
			if(p_type->IsType<StemCellProliferativeType>())
			{
				neighbouring_epidermal_indices.push_back(*elem_iter);
			}

		}
	}

	return neighbouring_epidermal_indices;
}

//Method to calculate the force due to the basement membrane on an Epidermal cell
c_vector<double, 2> EpidermalBasementMembraneForce::CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, 
																						std::vector<unsigned> epidermalIndices,
																						c_vector<double, 2> epidermalHeightExtremes,
                                                            							c_vector<double, 2> epidermalWidthExtremes, 
																						unsigned nodeIndex)
{

	// Get the considered node's location
	c_vector<double, 2> centre_point = rCellPopulation.GetNode(nodeIndex)->rGetLocation();

	//Get the left and neighbours of the node index, based on the cosine approximation of the dermal interface.
	std::vector<unsigned> closest_left_neighbours = GetClosestNeighboursBasedOnCosineApproximation(rCellPopulation, 
																									epidermalIndices,
																									epidermalHeightExtremes,
																									epidermalWidthExtremes,
																									nodeIndex,
																									-1.0);
	std::vector<unsigned> closest_right_neighbours = GetClosestNeighboursBasedOnCosineApproximation(rCellPopulation,
																									epidermalIndices,
																									epidermalHeightExtremes,
																									epidermalWidthExtremes,
																									nodeIndex,
																									1.0);

	// Initialise the left and right node indices and left and right points
	c_vector<double, 2> left_point, right_point; 
	// If there is a closest left and right neighbour, this is straightforward
	if ( (!closest_left_neighbours.empty())&&(!closest_right_neighbours.empty()) )
	{
		// Get the indices
		unsigned left_node_index = GetNearestNeighbourAlongCosineApproximation(rCellPopulation, 
																				epidermalIndices,
																				epidermalHeightExtremes,
																				epidermalWidthExtremes,
																				nodeIndex,
																				-1.0);
		unsigned right_node_index = GetNearestNeighbourAlongCosineApproximation(rCellPopulation,
																				epidermalIndices,
																				epidermalHeightExtremes,
																				epidermalWidthExtremes,
																				nodeIndex,
																				1.0);

		// Get the points
		left_point = rCellPopulation.GetNode(left_node_index)->rGetLocation();
		right_point = rCellPopulation.GetNode(right_node_index)->rGetLocation();
	}
	else if ( (!closest_left_neighbours.empty())&&(closest_right_neighbours.empty()) )// No right neighbour, but there is a left neighbour
	{
		// Get the left index and point
		unsigned left_node_index = GetNearestNeighbourAlongCosineApproximation(rCellPopulation, 
																				epidermalIndices,
																				epidermalHeightExtremes,
																				epidermalWidthExtremes,
																				nodeIndex,
																				-1.0);
		left_point = rCellPopulation.GetNode(left_node_index)->rGetLocation();

		// Obtain the right neighbour by rotating the vector from the centre node to the left node.
		// We define the angle via the inner product between the centre-to-left vector and the vector from
		// the closest fibroblast neighbour to the centre node.

		// Get the centre-to-left vector
		c_vector<double, 2> centre_to_left = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, left_point);
		double centre_to_left_norm = norm_2(centre_to_left);
		centre_to_left /= centre_to_left_norm; // Normalise the vector length
		
		// Get the fibroblast-to-centre vector
		unsigned closest_fibroblast_index = GetNearestFibroblastNeighbour(rCellPopulation, nodeIndex);
		c_vector<double, 2> closest_fibroblast_point = rCellPopulation.GetNode(closest_fibroblast_index)->rGetLocation();

		c_vector<double, 2> fibroblast_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(closest_fibroblast_point, centre_point);
		fibroblast_to_centre /= norm_2(fibroblast_to_centre); // Normalise the vector length

		// Reflect the fibroblast-to-centre vector so they're pointing in the same quadrant
		if (inner_prod(fibroblast_to_centre, centre_to_left) < 0.0)
		{
			fibroblast_to_centre *= -1.0;
		}

		// We can now define the right point via a clockwise rotation
		double theta = acos(inner_prod(fibroblast_to_centre, centre_to_left));

		// Rotate the centre point
		right_point[0] = centre_to_left_norm*(fibroblast_to_centre[0]*cos(theta) + fibroblast_to_centre[1]*sin(theta));
		right_point[1] = centre_to_left_norm*(fibroblast_to_centre[1]*cos(theta) - fibroblast_to_centre[0]*sin(theta));

		// // Define the right neighbour by reflecting the vector from the considered node to the right neighbour
		// right_point = centre_point + rCellPopulation.rGetMesh().GetVectorFromAtoB(left_point, centre_point);

	}
	else if ( (closest_left_neighbours.empty())&&(!closest_right_neighbours.empty()) )// No left neighbour, but there is a right neighbour
	{
		// Get the right index and point
		unsigned right_node_index = GetNearestNeighbourAlongCosineApproximation(rCellPopulation,
																				epidermalIndices,
																				epidermalHeightExtremes,
																				epidermalWidthExtremes,
																				nodeIndex,
																				1.0);
		right_point = rCellPopulation.GetNode(right_node_index)->rGetLocation();

		// Obtain the left neighbour by rotating the vector from the centre node to the left node.
		// We define the angle via the inner product between the centre-to-right vector and the vector from
		// the closest fibroblast neighbour to the centre node.

		// Get the centre-to-right vector
		c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, right_point);
		double centre_to_right_norm = norm_2(centre_to_right);
		centre_to_right /= centre_to_right_norm; // Normalise the vector length

		// Get the fibroblast-to-centre vector
		unsigned closest_fibroblast_index = GetNearestFibroblastNeighbour(rCellPopulation, nodeIndex);
		c_vector<double, 2> closest_fibroblast_point = rCellPopulation.GetNode(closest_fibroblast_index)->rGetLocation();

		c_vector<double, 2> fibroblast_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(closest_fibroblast_point, centre_point);
		fibroblast_to_centre /= norm_2(fibroblast_to_centre); // Normalise the vector length

		// Reflect the fibroblast-to-centre vector so they're pointing in the same quadrant
		if (inner_prod(fibroblast_to_centre, centre_to_right) < 0.0)
		{
			fibroblast_to_centre *= -1.0;
		}

		// We can now define the right point via a clockwise rotation
		double theta = acos(inner_prod(fibroblast_to_centre, centre_to_right));

		// Rotate the centre point
		left_point[0] = centre_to_right_norm*(fibroblast_to_centre[0]*cos(theta) - fibroblast_to_centre[1]*sin(theta));
		left_point[1] = centre_to_right_norm*(fibroblast_to_centre[1]*cos(theta) + fibroblast_to_centre[0]*sin(theta));

		// // Define the right neighbour by reflecting the vector from the considered node to the right neighbour
		// left_point = centre_point + rCellPopulation.rGetMesh().GetVectorFromAtoB(right_point, centre_point);

	}
	else // Hopefully this never happens. We COULD deal with this, but it really shouldn't happen.
	{
		EXCEPTION("Epidermal cell is isolated and has no neighbours.");
	}

	double curvature = FindParametricCurvature(rCellPopulation, left_point, centre_point, right_point);

	//Get the unit vectors from the centre points to its left and right neighbours
	// c_vector<double, 2> centre_to_left = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, left_point);
	// centre_to_left /= norm_2(centre_to_left); //Normalise vector

	// c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, right_point);
	// centre_to_right /= norm_2(centre_to_right); //Normalise vector

	/* Define direction as the perpendicular bisector from the left to right point
		*/
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(left_point, right_point);


	assert(norm_2(left_to_right) != 0.0);

	//Define force direction
	c_vector<double, 2> force_direction;
	force_direction(0) = left_to_right[1];
	force_direction(1) = -left_to_right[0];

	force_direction /= norm_2(left_to_right);

	double target_curvature = GetTargetCurvature(); // Get the target curvature

	/* We now ensure the vector is pointing in the appropriate direction
	* (it will always point "down" and "outwards" initially).
	* */
	if (target_curvature > 0.0) //If the force has overshot the target curvature, we need to reverse the force direction
	{
		if ( (curvature > 0.0)&&(curvature - target_curvature > 0.0) ) //If points look like V and the 'v' is too pointy, we send it away from the CoM
		{
			force_direction *= -1.0;
		}

	}
	else if (target_curvature < 0.0) //Similar situation but with "/\"
	{
		if ( (curvature < 0.0)&&(curvature - target_curvature < 0.0) )
		{
			force_direction *= -1.0;
		}
	}
	else //Reverse the force direction if we get a "V"
	{
		if (curvature > 0.0)
		{
			force_direction *= -1.0;
		}
	}

	//Initialise the force vector
	c_vector<double, 2> force_due_to_basement_membrane;

	// Get the basement membrane stiffness and target curvature
	double basement_membrane_parameter = GetBasementMembraneParameter(); //Get the basement membrane stiffness

	force_due_to_basement_membrane = basement_membrane_parameter*( fabs(curvature - target_curvature) )*force_direction;

	// c_vector<double, 2> force_due_to_basement_membrane = zero_vector<double>(2);
	return force_due_to_basement_membrane;
}

//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void EpidermalBasementMembraneForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{

	//Get the epidermal indices
	std::vector<unsigned> epidermal_indices = GetEpidermalIndices(rCellPopulation);

	// Get the epidermal width and height extremes
	c_vector<double, 2> epidermis_widths = GetEpidermisWidthExtremes(rCellPopulation);
	c_vector<double, 2> epidermis_heights = GetEpidermisHeightExtremes(rCellPopulation, epidermis_widths);

	// Iterate over each node and apply the force
	for (unsigned i = 0; i < epidermal_indices.size(); i++)
	{
		unsigned epidermal_index = epidermal_indices[i];

		// Calculate the basement membrane force
		c_vector<double, 2> force_on_node = CalculateForceDueToBasementMembrane(rCellPopulation, epidermal_indices, epidermis_heights, epidermis_widths, epidermal_index);

		// // Apply the force
		rCellPopulation.GetNode(epidermal_index)->AddAppliedForceContribution(force_on_node);
	}
}

void EpidermalBasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<TargetCurvature>" << mTargetCurvature << "</TargetCurvature> \n";

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EpidermalBasementMembraneForce)
