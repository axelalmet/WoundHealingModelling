#ifndef DISTANCEBASEDEPIDERMALBASEMENTMEMBRANEFORCE_HPP_
#define DISTANCEBASEDEPIDERMALBASEMENTMEMBRANEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "RandomNumberGenerator.hpp"

#include <cmath>
#include <list>
#include <fstream>

/**
 * MODIFIED FOR PERSONAL USE BY AXEL ALMET
 * A force class that defines the force due to the basement membrane, where
 * the strength of attachment is based on the distance to the basement membrane,
 * as based off the Du et al. (2018) and Wang et al. (2019) papers..
 */

class DistanceBasedEpidermalBasementMembraneForce : public AbstractForce<2>
{
    friend class TestCrossSectionModelInteractionForce;

private :

    /** Parameter that multiplies the curvature to give the basement membrane force */
    double mBasementMembraneParameter;

    /** The cut off radius for defining neighbouring nodes */
    double mCutOffRadius;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mBasementMembraneParameter;
        archive & mCutOffRadius;
    }

public :

    /**
     * Constructor.
     */
	DistanceBasedEpidermalBasementMembraneForce();

    /**
     * Destructor.
     */
    ~DistanceBasedEpidermalBasementMembraneForce();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Get method for cut off radius */
    double GetCutOffRadius();

    /* Set method for cut off radius */
    void SetCutOffRadius(double cutOffRadius);

    /*
     * Method to get the indices of the monolayer
     */
    std::vector<unsigned> GetEpidermalIndices(AbstractCellPopulation<2>& rCellPopulation);


    /*
     * Return vector of Collagen indices that are close to the considered Epidermal node,
     * based on a specified neighbourhood radius
     * 
     * @param rCellPopulation the cell population
     * @param epidermalIndex the considered epidermal node index
     */
    std::vector<unsigned> GetNeighbouringCollagenIndices(AbstractCellPopulation<2>& rCellPopulation, 
                                                                        unsigned epidermalIndex);

    /*
     * Return force direction with respect ot the basement membrane, which accounts for whether
     * or not the direction from the epidermal index to the average Collagen position should be 
     * reflected.
     * 
     * @param rCellPopulation the cell population
     * @param CollagenIndices the neighbouring Collagen indices to the epidermal cell
     * @param epidermalIndex the epidermal cell
     */
    // c_vector<double, 2> CalculateForceDirection(AbstractCellPopulation<2>& rCellPopulation, 
    //                                                                 std::vector<unsigned> CollagenIndices,
    //                                                                 unsigned epidermalIndex);

    /*
     * Method to calculate the force due to basement membrane on an Epidermal cell
     */
    c_vector<double, 2> CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, 
                                                            std::vector<unsigned> epidermalIndices,
                                                            unsigned nodeIndex);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DistanceBasedEpidermalBasementMembraneForce)

#endif /*DISTANCEBASEDEPIDERMALBASEMENTMEMBRANEFORCE_HPP_*/
