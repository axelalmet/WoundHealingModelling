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

#ifndef POLARITYTRACKINGMODIFIER_HPP_
#define POLARITYTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class in which the polarity of non-ECM-type cells
 * are updated. Polarity is assumed to align to cell velocity,
 * which is influenced by various mechanical forces.
 */
template<unsigned DIM>
class PolarityTrackingModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mVelocityReorientationStrength;
        archive & mFibreReorientationStrength;
        archive & mNeighbourhoodRadius;
        archive & mMorphogenThreshold;
        archive & mMeanActivationLifetime;
        archive & mFibreDepositionProbability;
        archive & mCollagenProductionRate;
        archive & mCollagenDegradationRate;
    }
protected:

    /*
     *  The remodelling rate of the polarity to the cell velocity.
     */
    double mVelocityReorientationStrength;

    /*
     * The remodelling rate of the polarity to the 'average' fibre orientation at that position
     */
    double mFibreReorientationStrength;

    /*
     * Neighbourhood radius
     */
    double mNeighbourhoodRadius;

    /*
     * Morphogen threshold for determining whether or not a fibroblast rearranges the local
     * ECM fibres
     */
    double mMorphogenThreshold;

    /*
     * Mean activation lifetime for a fibroblast once it is no longer exposed to a sufficient
     * concentration of morphogen
     */
    double mMeanActivationLifetime;

    /*
     * Probability of an activated fibroblast depositing a collagen fibre in front of it.
     */
    double mFibreDepositionProbability; 

    /* 
     * Rate of collagen production (needed for newly-deposited collagen fibres)
     */
    double mCollagenProductionRate;

    /* 
     * Rate of collagen degradation (needed for newly-deposited collagen fibres)
     */
    double mCollagenDegradationRate;

public:

    /**
     * Default constructor.
     */
    PolarityTrackingModifier();

    /**
     * Destructor.
     */
    virtual ~PolarityTrackingModifier();

    /*
     * Get the reorientation strength parameters
     */
    double GetVelocityReorientationStrength();


    /*
     * Set the reorientation strength parameters
     * 
     * @param reorientationStrength the new strength of reorientation
     */
    void SetVelocityReorientationStrength(double velocityReorientationStrength);

    /*
     * Get the reorientation strength parameters
     */
    double GetFibreReorientationStrength();

    /*
     * Set the reorientation strength parameters
     * 
     * @param reorientationStrength the new strength of reorientation
     */
    void SetFibreReorientationStrength(double fibreReorientationStrength);

    /* 
     * Get the neighbourhood radius used to calculate the local fibre orientation
     */
    double GetNeighbourhoodRadius();

    /*
     * Set the neighbourhood radius
     * 
     * @param neighbourhoodRadius the new neighbourhood radius
     */
    void SetNeighbourhoodRadius(double neighbourhoodRadius);

    /* 
     * Get the morphogen threshold used to determine whether or not a fibroblast will
     * sculpt the local ECM structure
     */
    double GetMorphogenThreshold();

    /*
     * Set the morphogen threshold
     * 
     * @param morphogenThreshold the new morphogen threshold
     */
    void SetMorphogenThreshold(double morphogenThreshold);

    /*
     * Get the mean activation lifetime parameter for fibroblasts
     * that are no longer 
     */
    double GetMeanActivationLifetime();

    /*
     * Set the mean activation lifetime for fibroblasts
     * 
     * @param meanActivationLifetime the new mean activationlifetime
     */
    void SetMeanActivationLifetime(double meanActivationLifetime);

    /*
     * Get the probability of depositing a collagen fibre (activated fibroblasts only)
     */
    double GetFibreDepositionProbability();

    /*
     * Set the probability of ECM fibre deposition
     * 
     * @param fibreDepositionProbability
     */
    void SetFibreDepositionProbability(double fibreDepositionProbability);

    /*
     * Get the collagen production rate (for newly-created fibres)
     */
    double GetCollagenProductionRate();

    /*
     * Set the rate of collagen production
     * 
     * @param collagenProducitonRate
     */
    void SetCollagenProductionRate(double collagenProductionRate);

    /*
     * Get the collagen degradation rate (for newly-created fibres)
     */
    double GetCollagenDegradationRate();

    /*
     * Set the rate of collagen degradation
     * 
     * @param collagenDegradationRate
     */
    void SetCollagenDegradationRate(double collagenDegradationRate);

      
     /* Function to calculate the local fibre orientation and density
     * based on an inverse distance weighting sampling method,
     * so that the fibre orientation and density is more similar
     * to the ECM nodes that are closer.
     *
     * @param fibroblastIndex the index of the fibroblast cell
     * @param rCellPopulation the cell population
     */
    c_vector<double, 2> GetLocalFibreOrientationAndDensity(unsigned fibroblastIndex, AbstractCellPopulation<DIM,DIM>& rCellPopulation);

      /* 
     * This method determines whether or not a spanning ECM (fibrin or collagen)
     * will intersect with the cell within a reasonable range (determined by cut-off length),
     * which will lthen determine whether or not we apply the cell-matrix force.
     * 
     * @param CellPtr pCell, the non-ECM cell
     * @param cellLocation, position of the cell
     * @param CellPtr pMatrixCell, the ECM cell
     * @param matrixCellLocation, the ECM node location
     * 
     * @return Whether or not the fibre intersects with the cell
     */
    bool DoesCellIntersectWithFibre(CellPtr pFibroblastCell, c_vector<double, DIM> fibroblastLocation,
                                    CellPtr pMatrixCell, c_vector<double, DIM> matrixCellLocation);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the mean level of Delta in each cell's neighbours and store these in the CellData.
     *
     * Note: If using a CaBasedCellPopulation, we assume a Moore neighbourhood and unit carrying capacity.
     * If a cell has no neighbours (such as an isolated cell in a CaBasedCellPopulation), we store the
     * value -1 in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityTrackingModifier)

#endif /*PolarityTrackingMODIFIER_HPP_*/
