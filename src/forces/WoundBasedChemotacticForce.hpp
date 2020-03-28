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

#ifndef WOUNDBASEDCHEMOTACTICFORCE_HPP_
#define WOUNDBASEDCHEMOTACTICFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/**
 * A WoundBasedChemotactic force class.
 */
template<unsigned DIM>
class WoundBasedChemotacticForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    /*
     * Chemotactic Force strenghth parameter
     */
    double mChemotacticStrength;

    /*
     * Radius of interaction to determine which neighbours we compare the concentrations to.
     */ 
    double mNeighbourhoodRadius;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mChemotacticStrength;
        archive & mNeighbourhoodRadius;
    }

public:

    /**
     * Constructor.
     */
    WoundBasedChemotacticForce();

    /**
     * Destructor.
     */
    ~WoundBasedChemotacticForce();

    /*
     * Get chemotactic force strength
     * 
     * @return mChemotacticStrength
     */
    double GetChemotacticStrength();

    /*
     * Set the chemotactic force strength
     * 
     * @param chemotacticStrength
     */
    void SetChemotacticStrength(double chemotacticStrength);


    /*
     * Get neighbourhood radius
     * 
     * @return mNeighbourhoodRadius
     */
    double GetNeighbourhoodRadius();

    /*
     * Set the neighbourhood radius
     * 
     * @param neighbourhoodRadius
     */
    void SetNeighbourhoodRadius(double neighbourhoodRadius);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     * Fc = chi(C,|gradC|) gradC/|gradC|  (if |gradC|>0, else Fc = 0)
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WoundBasedChemotacticForce)

#endif /*WOUNDBASEDChemotacticForce_HPP_*/
