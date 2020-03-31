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

#include "RandomSymmetricDivisionBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RandomSymmetricDivisionBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::RandomSymmetricDivisionBasedDivisionRule(double& rSymmetricDivisionProbability)
{
    mSymmetricDivisionProbability = rSymmetricDivisionProbability;
}

/*
* Get probability of symmetric division
*/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const double& RandomSymmetricDivisionBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::rGetSymmetricDivisionProbability() const
{
    return mSymmetricDivisionProbability;
}

/*
* Return closest fibroblast index to the considered epidermal index
* 
* @param rCellPopulation the cell population
* @param epidermalIndex the considered epidermal node index
*/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned RandomSymmetricDivisionBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::GetNearestFibroblastNeighbour(AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, unsigned epidermalIndex)
{
    double min_fibroblast_distance = DBL_MAX;
	unsigned closest_fibroblast_index = 0;

	// Get the location of the considered Epidermal node
	c_vector<double, SPACE_DIM> epidermal_location = rCellPopulation.GetNode(epidermalIndex)->rGetLocation();

	// This method really only applies to NodeBasedCellPopulations
	if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
	{
		// NodeBasedCellPopulation<SPACE_DIM>* p_tissue = static_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);

		// Find the indices of the elements owned by this node
		std::set<unsigned> neighbouring_indices = rCellPopulation.GetNeighbouringNodeIndices(epidermalIndex);

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > RandomSymmetricDivisionBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Get separation parameter
    double separation = rCellPopulation.GetMeinekeDivisionSeparation();
    
    // Initialise division vector
    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions;

    //Get the cell type
    boost::shared_ptr<AbstractCellProperty> p_type = pParentCell->GetCellProliferativeType();

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        if (p_type->IsType<StemCellProliferativeType>()) // Epidermal stem cells
        {
            // Get the cell location index
            unsigned parent_cell_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);

            // Determine where the nearest fibroblast neighbour is
            unsigned nearest_fibroblast_index = GetNearestFibroblastNeighbour(rCellPopulation, parent_cell_index);

            // Get the locations of the cells
            c_vector<double, SPACE_DIM> stem_position = rCellPopulation.GetLocationOfCellCentre(pParentCell); // Parent cell

            c_vector<double, SPACE_DIM> fibroblast_position = rCellPopulation.GetNode(nearest_fibroblast_index)->rGetLocation(); // Fibroblast

            // Get the vector between the two
            c_vector<double, SPACE_DIM> fibroblast_to_stem = stem_position - fibroblast_position;
            
            fibroblast_to_stem /= norm_2(fibroblast_to_stem); // Normalise the vector

            /* If the probability is less than mSymmetricDivisionProbability, the division direction is
            * roughly parallel to the basement membrane. If not, the division direction is roughly
            * perpendicular to the basement membrane.
            */

            double division_probability = RandomNumberGenerator::Instance()->ranf();

            double symmetric_division_probability = mSymmetricDivisionProbability;

            // Initialise parent and daughter positions
            c_vector<double, SPACE_DIM> parent_position, daughter_position;


            if (division_probability < symmetric_division_probability)
            {
                // Rotate the division direction 90 degrees
                c_vector<double, SPACE_DIM> division_vector;

                division_vector(0) = -fibroblast_to_stem(1);
                division_vector(1) = fibroblast_to_stem(0);

                parent_position = stem_position - 0.5*division_vector;
                daughter_position = stem_position + 0.5*division_vector;

            }
            else // Division direction is perpendicular to basement membrane, i.e. from fibroblast to stem cell
            {
                parent_position = stem_position;

                daughter_position = stem_position + separation*fibroblast_to_stem;
            }

            // Division division pair
            positions.first = parent_position;
            positions.second = daughter_position;
        }
        else // Essentially for anything not a stem cell 
        {

            // Make a random direction vector of the required length
            c_vector<double, SPACE_DIM> random_vector;

            /*
            * Pick a random direction and move the parent cell backwards by 0.5*separation
            * in that direction and return the position of the daughter cell 0.5*separation
            * forwards in that direction.
            */
            switch (SPACE_DIM)
            {
                case 1:
                {
                    double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                    random_vector(0) = 0.5*separation*random_direction;
                    break;
                }
                case 2:
                {
                    double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                    random_vector(0) = 0.5*separation*cos(random_angle);
                    random_vector(1) = 0.5*separation*sin(random_angle);
                    break;
                }
                case 3:
                {
                    /*
                    * Note that to pick a random point on the surface of a sphere, it is incorrect
                    * to select spherical coordinates from uniform distributions on [0, 2*pi) and
                    * [0, pi) respectively, since points picked in this way will be 'bunched' near
                    * the poles. See #2230.
                    */
                    double u = RandomNumberGenerator::Instance()->ranf();
                    double v = RandomNumberGenerator::Instance()->ranf();

                    double random_azimuth_angle = 2*M_PI*u;
                    double random_zenith_angle = std::acos(2*v - 1);

                    random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
                    random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
                    random_vector(2) = 0.5*separation*cos(random_zenith_angle);
                    break;
                }
                default:
                    // This can't happen
                    NEVER_REACHED;
            }

            c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - random_vector;
            c_vector<double, SPACE_DIM> daughter_position = parent_position + random_vector;

            positions.first = parent_position;
            positions.second = daughter_position;

        }
    }

    return positions;
}

// Explicit instantiation
template class RandomSymmetricDivisionBasedDivisionRule<1,1>;
template class RandomSymmetricDivisionBasedDivisionRule<1,2>;
template class RandomSymmetricDivisionBasedDivisionRule<2,2>;
template class RandomSymmetricDivisionBasedDivisionRule<1,3>;
template class RandomSymmetricDivisionBasedDivisionRule<2,3>;
template class RandomSymmetricDivisionBasedDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RandomSymmetricDivisionBasedDivisionRule)
