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

#ifndef BASEMENTMEMBRANEDISTANCEBASEDCELLKILLER_HPP_
#define BASEMENTMEMBRANEDISTANCEBASEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 *  Distance-based cell killer, where differentiated cells that
 *  are a specified distance away from the basement membrane are 
 *  sloughed.
 */
class BasementMembraneDistanceBasedCellKiller : public AbstractCellKiller<2>
{
private:


    /** Radius of death. */
    double mRadius;

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
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population.
     * @param radius the radius of death.
     */
    BasementMembraneDistanceBasedCellKiller(AbstractCellPopulation<2>* pCellPopulation,
                              double radius);

    /**
     * @return mRadius.
     */
    double GetRadius() const;

    /**
     * Loop over cells and kills cells outside boundary.
     */
    virtual void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BasementMembraneDistanceBasedCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BasementMembraneDistanceBasedCellKiller.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const BasementMembraneDistanceBasedCellKiller * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
    ar & p_cell_population;
    double radius = t->GetRadius();
    ar & radius;
}

/**
 * De-serialize constructor parameters and initialise a BasementMembraneDistanceBasedCellKiller.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, BasementMembraneDistanceBasedCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar & p_cell_population;
    double radius;
    ar & radius;

    // Invoke inplace constructor to initialise instance
    ::new(t)BasementMembraneDistanceBasedCellKiller(p_cell_population, radius);
}
}
} // namespace ...


#endif /*BASEMENTMEMBRANEDISTANCEBASEDCELLKILLER_HPP_*/

