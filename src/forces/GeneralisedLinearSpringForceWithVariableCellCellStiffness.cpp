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

#include "GeneralisedLinearSpringForceWithVariableCellCellStiffness.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "BloodCellProliferativeType.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GeneralisedLinearSpringForceWithVariableCellCellStiffness()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mMeinekeSpringStiffness(15.0),        // denoted by mu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     mMeinekeDivisionRestingSpringLength(0.5),
     mMeinekeSpringGrowthDuration(1.0),
     mStemStemMultiplicationFactor(1.0),
     mStemDifferentiatedMultiplicationFactor(1.0),
     mStemFibroblastMultiplicationFactor(1.0),
     mDifferentiatedDifferentiatedMultiplicationFactor(1.0),
     mDifferentiatedFibroblastMultiplicationFactor(1.0),
     mFibroblastFibroblastMultiplicationFactor(1.0),
     mStemPlateletMultiplicationFactor(1.0),
     mDifferentiatedPlateletMultiplicationFactor(1.0),
     mFibroblastPlateletMultiplicationFactor(1.0),
     mPlateletPlateletMultiplicationFactor(1.0)

{
    if (SPACE_DIM == 1)
    {
        mMeinekeSpringStiffness = 30.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetStemStemMultiplicationFactor(double stemStemMultiplicationFactor)
{
    mStemStemMultiplicationFactor = stemStemMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetStemStemMultiplicationFactor()
{
    return mStemStemMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetStemDifferentiatedMultiplicationFactor(double stemDifferentiatedMultiplicationFactor)
{
    mStemDifferentiatedMultiplicationFactor = stemDifferentiatedMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetStemDifferentiatedMultiplicationFactor()
{
    return mStemDifferentiatedMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetStemFibroblastMultiplicationFactor(double stemFibroblastMultiplicationFactor)
{
    mStemFibroblastMultiplicationFactor = stemFibroblastMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetStemFibroblastMultiplicationFactor()
{
    return mStemFibroblastMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetDifferentiatedDifferentiatedMultiplicationFactor(double differentiatedDifferentiatedMultiplicationFactor)
{
    mDifferentiatedDifferentiatedMultiplicationFactor = differentiatedDifferentiatedMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetDifferentiatedDifferentiatedMultiplicationFactor()
{
    return mDifferentiatedDifferentiatedMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetDifferentiatedFibroblastMultiplicationFactor(double differentiatedFibroblastMultiplicationFactor)
{
    mDifferentiatedFibroblastMultiplicationFactor = differentiatedFibroblastMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetDifferentiatedFibroblastMultiplicationFactor()
{
    return mDifferentiatedFibroblastMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetFibroblastFibroblastMultiplicationFactor(double fibroblastFibroblastMultiplicationFactor)
{
    mFibroblastFibroblastMultiplicationFactor = fibroblastFibroblastMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetFibroblastFibroblastMultiplicationFactor()
{
    return mFibroblastFibroblastMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetStemPlateletMultiplicationFactor(double stemPlateletMultiplicationFactor)
{
    mStemPlateletMultiplicationFactor = stemPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetStemPlateletMultiplicationFactor()
{
    return mStemPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetDifferentiatedPlateletMultiplicationFactor(double differentiatedPlateletMultiplicationFactor)
{
    mDifferentiatedPlateletMultiplicationFactor = differentiatedPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetDifferentiatedPlateletMultiplicationFactor()
{
    return mDifferentiatedPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetFibroblastPlateletMultiplicationFactor(double fibroblastPlateletMultiplicationFactor)
{
    mFibroblastPlateletMultiplicationFactor = fibroblastPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetFibroblastPlateletMultiplicationFactor()
{
    return mFibroblastPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetPlateletPlateletMultiplicationFactor(double plateletPlateletMultiplicationFactor)
{
    mPlateletPlateletMultiplicationFactor = plateletPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetPlateletPlateletMultiplicationFactor()
{
    return mPlateletPlateletMultiplicationFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                                                     unsigned nodeBGlobalIndex,
                                                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                                                     bool isCloserThanRestLength)
{

    double multiplication_factor;

    // Get the cell proliferative types
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    if ( (p_cell_A->GetCellProliferativeType()->IsType<StemCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) ) // Stem-stem connection
    {
        multiplication_factor = mStemStemMultiplicationFactor;
    }
    else if ( ( (p_cell_A->GetCellProliferativeType()->IsType<StemCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) )
            ||  ( (p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) ) )// Stem-differentiated connection
    {
        multiplication_factor = mStemDifferentiatedMultiplicationFactor;
    }
    else if ( ( (p_cell_A->GetCellProliferativeType()->IsType<StemCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>()) )
            ||  ( (p_cell_A->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) ) )// Stem-fibroblast connection
    {
        multiplication_factor = mStemFibroblastMultiplicationFactor;
    }
    else if ( (p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) ) // Differentiated-differentiated
    {
        multiplication_factor = mDifferentiatedDifferentiatedMultiplicationFactor;
    }
    else if ( ( (p_cell_A->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) )
            ||  ( (p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>()) ) )// Fibroblast-differentiated connection
    {
        multiplication_factor =  mDifferentiatedFibroblastMultiplicationFactor;
    }
    else if ( (p_cell_A->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>())&&(p_cell_A->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>()) ) // Fibroblast-fibroblast
    {
        multiplication_factor =  mFibroblastFibroblastMultiplicationFactor;
    }
    else if ( ( (p_cell_A->GetCellProliferativeType()->IsType<StemCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<BloodCellProliferativeType>()) )
            ||  ( (p_cell_A->GetCellProliferativeType()->IsType<BloodCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) ) )// Stem-platelet connection
    {
        multiplication_factor = mStemPlateletMultiplicationFactor;
    }
    else if ( ( (p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<BloodCellProliferativeType>()) )
            ||  ( (p_cell_A->GetCellProliferativeType()->IsType<BloodCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) ) )// Differentiated-platelet connection
    {
        multiplication_factor = mDifferentiatedPlateletMultiplicationFactor;
    }
    else if ( ( (p_cell_A->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<BloodCellProliferativeType>()) )
            ||  ( (p_cell_A->GetCellProliferativeType()->IsType<BloodCellProliferativeType>())&&(p_cell_B->GetCellProliferativeType()->IsType<FibroblastCellProliferativeType>()) ) )// Fibroblast-platelet connection
    {
        multiplication_factor =  mFibroblastPlateletMultiplicationFactor;
    }
    else if ( (p_cell_A->GetCellProliferativeType()->IsType<BloodCellProliferativeType>())&&(p_cell_A->GetCellProliferativeType()->IsType<BloodCellProliferativeType>()) ) // Platelet-platelet
    {
        multiplication_factor =  mPlateletPlateletMultiplicationFactor;
    }
    else
    {
        multiplication_factor = 1.0;
    } 
    
    return multiplication_factor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::~GeneralisedLinearSpringForceWithVariableCellCellStiffness()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    const c_vector<double, SPACE_DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_b_location = p_node_b->rGetLocation();

    // Get the node radii for a NodeBasedCellPopulation
    double node_a_radius = 0.0;
    double node_b_radius = 0.0;

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        node_a_radius = p_node_a->GetRadius();
        node_b_radius = p_node_b->GetRadius();
    }

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, SPACE_DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutOffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(SPACE_DIM); // c_vector<double,SPACE_DIM>() is not guaranteed to be fresh memory
        }
    }

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     * value of 1.0.
     */
    double rest_length_final = 1.0;

    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        rest_length_final = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);
    }
    else if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        rest_length_final = node_a_radius+node_b_radius;
    }

    double rest_length = rest_length_final;

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

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
        AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        // if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        // {
        // Spring rest length increases from a small value to the normal rest length over 1 hour
        double lambda = mMeinekeDivisionRestingSpringLength;
        rest_length = lambda + (rest_length_final - lambda) * ageA/mMeinekeSpringGrowthDuration;
        // }
        if (ageA + SimulationTime::Instance()->GetTimeStep() >= mMeinekeSpringGrowthDuration)
        {
            // This spring is about to go out of scope
            p_static_cast_cell_population->UnmarkSpring(cell_pair);
        }
    }

    /*
     * For apoptosis, progressively reduce the radius of the cell
     */
    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        a_rest_length = (node_a_radius/(node_a_radius+node_b_radius))*rest_length;
        b_rest_length = (node_b_radius/(node_a_radius+node_b_radius))*rest_length;
    }

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
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = mMeinekeSpringStiffness;

    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        return multiplication_factor * spring_stiffness * unit_difference * overlap;
    }
    else
    {
        // A reasonably stable simple force law
        if (is_closer_than_rest_length) //overlap is negative
        {
            //log(x+1) is undefined for x<=-1
            assert(overlap > -rest_length_final);
            c_vector<double, SPACE_DIM> temp = multiplication_factor*spring_stiffness * unit_difference * rest_length_final* log(1.0 + overlap/rest_length_final);
            return temp;
        }
        else
        {
            double alpha = 5.0;
            c_vector<double, SPACE_DIM> temp = multiplication_factor*spring_stiffness * unit_difference * overlap * exp(-alpha * overlap/rest_length_final);
            return temp;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetMeinekeSpringStiffness()
{
    return mMeinekeSpringStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetMeinekeDivisionRestingSpringLength()
{
    return mMeinekeDivisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::GetMeinekeSpringGrowthDuration()
{
    return mMeinekeSpringGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mMeinekeSpringStiffness = springStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mMeinekeDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mMeinekeSpringGrowthDuration = springGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GeneralisedLinearSpringForceWithVariableCellCellStiffness<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MeinekeSpringStiffness>" << mMeinekeSpringStiffness << "</MeinekeSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MeinekeDivisionRestingSpringLength>" << mMeinekeDivisionRestingSpringLength << "</MeinekeDivisionRestingSpringLength>\n";
    *rParamsFile << "\t\t\t<MeinekeSpringGrowthDuration>" << mMeinekeSpringGrowthDuration << "</MeinekeSpringGrowthDuration>\n";
    *rParamsFile << "\t\t\t<StemStemMultiplicationFactor>" << mStemStemMultiplicationFactor << "</StemStemMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<StemDifferentiatedMultiplicationFactor>" << mStemDifferentiatedMultiplicationFactor << "</StemDifferentiatedMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<StemFibroblastMultiplicationFactor>" << mStemFibroblastMultiplicationFactor << "</StemFibroblastMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<DifferentiatedFibroblastMultiplicationFactor>" << mDifferentiatedFibroblastMultiplicationFactor << "</DifferentiatedFibroblastMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<DifferentiatedDifferentiatedMultiplicationFactor>" << mDifferentiatedDifferentiatedMultiplicationFactor << "</DifferentiatedDifferentiatedMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<FibroblastFibroblastMultiplicationFactor>" << mFibroblastFibroblastMultiplicationFactor << "</FibroblastFibroblastMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<StemPlateletMultiplicationFactor>" << mStemPlateletMultiplicationFactor << "</StemPlateletMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<DifferentiatedPlateletMultiplicationFactor>" << mDifferentiatedPlateletMultiplicationFactor << "</DifferentiatedPlateletMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<FibroblastPlateletMultiplicationFactor>" << mFibroblastPlateletMultiplicationFactor << "</FibroblastPlateletMultiplicationFactor>\n";
    *rParamsFile << "\t\t\t<PlateletPlateletMultiplicationFactor>" << mPlateletPlateletMultiplicationFactor << "</PlateletPlateletMultiplicationFactor>\n";


    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class GeneralisedLinearSpringForceWithVariableCellCellStiffness<1,1>;
template class GeneralisedLinearSpringForceWithVariableCellCellStiffness<1,2>;
template class GeneralisedLinearSpringForceWithVariableCellCellStiffness<2,2>;
template class GeneralisedLinearSpringForceWithVariableCellCellStiffness<1,3>;
template class GeneralisedLinearSpringForceWithVariableCellCellStiffness<2,3>;
template class GeneralisedLinearSpringForceWithVariableCellCellStiffness<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GeneralisedLinearSpringForceWithVariableCellCellStiffness)
