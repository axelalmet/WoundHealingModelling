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

#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp"
#include "ExtracellularMatrixCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::GeneralisedLinearSpringForceWithVariableInteractionDistance()
   : GeneralisedLinearSpringForce<ELEMENT_DIM,SPACE_DIM>(),
     mSpringStiffness(15.0),        // denoted by mu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     mMeinekeDivisionRestingSpringLength(0.5),
     mMeinekeSpringGrowthDuration(1.0),
     mApplyCellToEcmForce(false),
     mApplyForceOnMarkedSprings(false)
{
    if (SPACE_DIM == 1)
    {
        mSpringStiffness = 30.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::GetSpringConstantMultiplicationFactor(unsigned nodeAColour,
                                                                                     unsigned nodeBColour)
{
    /* 
     * As a first step, we just take the average of the two stiffness factors,
     * because it will reduce to the correct stiffness when the cell types are the same. 
     */
    double stiffness_factor_A = mStiffnessMultiplicationFactors[nodeAColour]; 
    double stiffness_factor_B = mStiffnessMultiplicationFactors[nodeBColour];

    double stiffness_mulltiplier = 0.5 * (stiffness_factor_A + stiffness_factor_B);

    return stiffness_mulltiplier; 

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::~GeneralisedLinearSpringForceWithVariableInteractionDistance()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    return zero_vector<double>(SPACE_DIM); // Just return zero here.
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::CalculateRestLength(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    c_vector<double, SPACE_DIM> unitDifference,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // Initialise the rest length
    double rest_length;

    // Get the cells (need them for cell type and cell data)
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    // Get the cell types
    boost::shared_ptr<AbstractCellProperty> p_cell_type_A = p_cell_A->GetCellProliferativeType();
    boost::shared_ptr<AbstractCellProperty> p_cell_type_B = p_cell_B->GetCellProliferativeType();

    // Get the node radii as well
    double node_A_radius = rCellPopulation.GetNode(nodeAGlobalIndex)->GetRadius();
    double node_B_radius = rCellPopulation.GetNode(nodeBGlobalIndex)->GetRadius();

    // Define the length scales for interaction
    double scale_A = p_cell_A->GetCellData()->GetItem("scale");
    double scale_B = p_cell_B->GetCellData()->GetItem("scale");

    // Define the values of sigma needed to calculate the rest length.
    // The sqrt(2) factors are to ensure that when sigmas are all the same,
    // it reduces to the case of two circles, i.e. l = 2R.
    double sigma_A_pa = sqrt(2.0) * node_A_radius;
    double sigma_A_pe = scale_A * sigma_A_pa;
    double sigma_B_pa = sqrt(2.0) * node_B_radius;
    double sigma_B_pe = scale_B * sigma_B_pa;

    // Get the migration directions
    double phi_A = p_cell_A->GetCellData()->GetItem("direction");
    double phi_B = p_cell_B->GetCellData()->GetItem("direction");
    
    // Define the orientation vectors from these angles phi
    c_vector<double, 2> direction_A, direction_B;
    direction_A[0] = cos(phi_A);
    direction_A[1] = sin(phi_A);

    direction_B[0] = cos(phi_B);
    direction_B[1] = sin(phi_B);

    // Calculate the numerator and denominator separately
    double sigma_numerator = (sigma_A_pa*sigma_A_pa + sigma_B_pe*sigma_B_pe) * (sigma_A_pe*sigma_A_pe + sigma_B_pa*sigma_B_pa)
                                - (sigma_A_pa*sigma_A_pa - sigma_A_pe*sigma_A_pe) * (sigma_B_pa*sigma_B_pa - sigma_B_pe*sigma_B_pe) * inner_prod(direction_A, direction_B) * inner_prod(direction_A, direction_B);
    double sigma_denominator = sigma_A_pa*sigma_A_pa + sigma_B_pa*sigma_B_pa
                                - (sigma_A_pa*sigma_A_pa - sigma_A_pe*sigma_A_pe)*inner_prod(unitDifference, direction_A)*inner_prod(unitDifference, direction_A)
                                - (sigma_B_pa*sigma_B_pa - sigma_B_pe*sigma_B_pe)*inner_prod(unitDifference, direction_B)*inner_prod(unitDifference, direction_B);

    // Finally calculate the rest length
    rest_length = pow(sigma_numerator/sigma_denominator, 0.5);

    return rest_length;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::DoesCellIntersectWithFibre(CellPtr pCell, c_vector<double, SPACE_DIM> cellLocation,
                                    CellPtr pMatrixCell, c_vector<double, SPACE_DIM> matrixCellLocation)
{
    double cell_direction = pCell->GetCellData()->GetItem("direction"); // Get the cell direction
    double fibre_direction = pMatrixCell->GetCellData()->GetItem("direction"); // Get the fibre directions

    // If the angles are the same, then they'll only intersect if they span from the same point
    if (cell_direction == fibre_direction)
    {
        return ( (cellLocation[0] == matrixCellLocation[0])&&(cellLocation[1] == matrixCellLocation[1]) );
    }
    else
    {
        return (fabs(((cellLocation[1] - matrixCellLocation[1])*cos(fibre_direction) 
                - (cellLocation[0] - matrixCellLocation[0])*sin(fibre_direction))
                    /sin(fibre_direction - cell_direction)) <= this->GetCutOffLength());
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::GetSpringStiffness()
{
    return mSpringStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::GetMeinekeDivisionRestingSpringLength()
{
    return mMeinekeDivisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::GetMeinekeSpringGrowthDuration()
{
    return mMeinekeSpringGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::SetSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mSpringStiffness = springStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mMeinekeDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mMeinekeSpringGrowthDuration = springGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::SetStiffnessMultiplicationFactors(std::map<unsigned, double> stiffnessMultiplicationFactors)
{
    mStiffnessMultiplicationFactors = stiffnessMultiplicationFactors;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::ApplyCellToEcmForce(bool applyCellToEcmForce)
{
    mApplyCellToEcmForce = applyCellToEcmForce;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::ApplyForceOnMarkedSprings(bool applyForceOnMarkedSprings)
{
    mApplyForceOnMarkedSprings = applyForceOnMarkedSprings;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableInteractionDistance<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{

    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

    std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& r_node_pairs = p_static_cast_cell_population->rGetNodePairs();

    // We iterate over the node pairs because it's easier to control which cell pairs undergo both adhesion or repulsion 
    // or only one of adhesion/repulsion
    for (typename std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >::iterator iter = r_node_pairs.begin();
        iter != r_node_pairs.end();
        iter++)
    {
        // Get the node pair and their indices
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > pair = *iter;

        unsigned node_a_index = pair.first->GetIndex();
        unsigned node_b_index = pair.second->GetIndex();
        
        // Get the cell types
        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(node_a_index);
        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(node_b_index);

        // Get the cell types
        boost::shared_ptr<AbstractCellProliferativeType> p_cell_type_A = p_cell_A->GetCellProliferativeType();
        boost::shared_ptr<AbstractCellProliferativeType> p_cell_type_B = p_cell_B->GetCellProliferativeType();

        // Get the colours associated with cell types
        unsigned node_a_colour = p_cell_type_A->GetColour();
        unsigned node_b_colour = p_cell_type_B->GetColour();

        // Get the node locations
        const c_vector<double, SPACE_DIM>& r_node_a_location = pair.first->rGetLocation();
        const c_vector<double, SPACE_DIM>& r_node_b_location = pair.second->rGetLocation();

        // Get the unit vector parallel to the line joining the two nodes
        c_vector<double, SPACE_DIM> unit_difference;

        c_vector<double, SPACE_DIM> force; // Initialise the force vector

        /*
        * We use the mesh method GetVectorFromAtoB() to compute the direction of the
        * unit vector along the line joining the two nodes, rather than simply subtract
        * their positions, because this method can be overloaded (e.g. to enforce a
        * periodic boundary in Cylindrical2dMesh).
        */
        unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

        // Calculate the distance between the two nodes
        double distance_between_nodes = norm_2(unit_difference);

        // Really only do anything if we're at non-zero and non-NAN distance
        if ( (distance_between_nodes > 0)&&(!std::isnan(distance_between_nodes)) )
        {
            unit_difference /= distance_between_nodes; // Normalilse the unit difference

            /*
            * Calculate the rest length of the spring connecting the two nodes using the new method
            * that accounts for the differences in shapes between epidermal and fibroblats
            */
            double rest_length_final = CalculateRestLength(node_a_index, node_b_index, unit_difference, rCellPopulation);

            /*
            * If mUseCutOffLength has been set, then there is zero force between
            * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
            */
            if (this->mUseCutOffLength)
            {
                if (distance_between_nodes <= (rest_length_final) * (this->GetCutOffLength()) )
                {
                    double rest_length = rest_length_final;

                    double ageA = p_cell_A->GetAge();
                    double ageB = p_cell_B->GetAge();

                    assert(!std::isnan(ageA));
                    assert(!std::isnan(ageB));

                    /*
                    * If the cells are both newly divided, then the rest length of the spring
                    * connecting them grows linearly with time, until 1 hour after division.
                    */
                    if (ageA < mMeinekeSpringGrowthDuration && ageB < mMeinekeSpringGrowthDuration)
                    {

                        // if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
                        // {
                        // Spring rest length increases from a small value to the normal rest length over 1 hour
                        double lambda = mMeinekeDivisionRestingSpringLength;
                        rest_length = lambda + (rest_length_final - lambda) * ageA/mMeinekeSpringGrowthDuration;
                        // }
                        if (ageA + SimulationTime::Instance()->GetTimeStep() >= mMeinekeSpringGrowthDuration)
                        {
                            // This spring is about to go out of scope
                            // p_static_cast_cell_population->UnmarkSpring(cell_pair);
                        }
                    }

                    /*
                    * For apoptosis, progressively reduce the radius of the cell
                    */
                    double a_rest_length = rest_length*0.5;
                    double b_rest_length = a_rest_length;

                    /*
                    * If either of the cells has begun apoptosis, then the length of the spring
                    * connecting them decreases linearly with time.
                    */
                    if (p_cell_A->HasApoptosisBegun())
                    {
                        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
                        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
                    }
                    if (p_cell_B->HasApoptosisBegun())
                    {
                        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
                        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
                    }

                    rest_length = a_rest_length + b_rest_length;

                            //assert(rest_length <= 1.0+1e-12); ///\todo #1884 Magic number: would "<= 1.0" do?

                    // Although in this class the 'spring constant' is a constant parameter, in
                    // subclasses it can depend on properties of each of the cells
                    double overlap = distance_between_nodes - rest_length;
                    bool is_closer_than_rest_length = (overlap <= 0);
                    double multiplication_factor = GetSpringConstantMultiplicationFactor(node_a_colour, node_b_colour);
                    double spring_stiffness = mSpringStiffness;

                    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
                    {
                        force = multiplication_factor * spring_stiffness * unit_difference * overlap;
                    }
                    else
                    {
                        // A reasonably stable simple force law
                        if (is_closer_than_rest_length) //overlap is negative
                        {
                            //log(x+1) is undefined for x<=-1
                            assert(overlap > -rest_length_final);
                            force = multiplication_factor*spring_stiffness * unit_difference * rest_length_final* log(1.0 + overlap/rest_length_final);
                        }
                        else
                        {
                            double alpha = 5.0;
                            force = multiplication_factor*spring_stiffness * unit_difference * overlap * exp(-alpha * overlap/rest_length_final);
                        }
                    }

                    // We now check for the cell types, because there's a bunch of conditions that we want to implement
                    if ( (p_cell_type_A->template IsType<FibroblastCellProliferativeType>()) 
                        &&(p_cell_type_B->template IsType<FibroblastCellProliferativeType>()) )
                    {
                        if (is_closer_than_rest_length) // Only apply a repulsive force
                        {
                            // Add the force contribution to each node
                            pair.first->AddAppliedForceContribution(force);
                            pair.second->AddAppliedForceContribution(-1.0*force);
                        }
                    }
                    else if ( (p_cell_type_A->template IsType<ExtracellularMatrixCellProliferativeType>()) // Fibrin fibres should only be under tension, if there's a force.
                        &&(p_cell_type_B->template IsType<ExtracellularMatrixCellProliferativeType>()) )
                    {
                        if (!is_closer_than_rest_length) // Only apply if the spring is under tension
                        {                   // If we only apply on marked springs, we should check
                            if (mApplyForceOnMarkedSprings)
                            {
                                std::pair<CellPtr, CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

                                if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
                                {
                                    pair.first->AddAppliedForceContribution(force);
                                    pair.second->AddAppliedForceContribution(-1.0*force);
                                }
                            }
                            else
                            {
                                pair.first->AddAppliedForceContribution(force);
                                pair.second->AddAppliedForceContribution(-1.0*force);
                            }
                        }
                    }
                    else if ( (p_cell_type_A->template IsType<ExtracellularMatrixCellProliferativeType>()) // Cell-ecm (fibrin) forces need to be considered
                        &&(!p_cell_type_B->template IsType<ExtracellularMatrixCellProliferativeType>()) )
                    {
                        // Only consider cells under tension and intersecting fibres/cells
                        if ((!is_closer_than_rest_length)
                                &&(DoesCellIntersectWithFibre(p_cell_B, r_node_b_location, p_cell_A, r_node_a_location)) )
                        {
                            if (mApplyCellToEcmForce)
                            {
                                    pair.first->AddAppliedForceContribution(force);
                                    pair.second->AddAppliedForceContribution(-1.0*force);
                            }
                            else // Else only apply the force on the cell
                            {
                                // Cell B is the non-fibrin cell
                                pair.second->AddAppliedForceContribution(-1.0*force);
                            }
                        }
                    }
                    else if ( (!p_cell_type_A->template IsType<ExtracellularMatrixCellProliferativeType>()) 
                        &&(p_cell_type_B->template IsType<ExtracellularMatrixCellProliferativeType>()) )
                    {
                        // Only consider cells under tension and intersecting fibres/cells
                        if ((!is_closer_than_rest_length)
                            &&(DoesCellIntersectWithFibre(p_cell_A, r_node_a_location, p_cell_B, r_node_b_location)) )
                        {
                            if (mApplyCellToEcmForce)
                            {
                                    pair.first->AddAppliedForceContribution(force);
                                    pair.second->AddAppliedForceContribution(-1.0*force);
                            }
                            else // Else only apply the force on the cell
                            {
                                // Cell A is the non-fibrin cell
                                pair.first->AddAppliedForceContribution(force); 
                            }
                        }
                    }
                    else
                    {
                        // Add the force contribution to each node
                        pair.first->AddAppliedForceContribution(force);
                        pair.second->AddAppliedForceContribution(-1.0*force);
                    }

                }
            }
        }
    }
}

// Explicit instantiation
template class GeneralisedLinearSpringForceWithVariableInteractionDistance<1,1>;
template class GeneralisedLinearSpringForceWithVariableInteractionDistance<1,2>;
template class GeneralisedLinearSpringForceWithVariableInteractionDistance<2,2>;
template class GeneralisedLinearSpringForceWithVariableInteractionDistance<1,3>;
template class GeneralisedLinearSpringForceWithVariableInteractionDistance<2,3>;
template class GeneralisedLinearSpringForceWithVariableInteractionDistance<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GeneralisedLinearSpringForceWithVariableInteractionDistance)
