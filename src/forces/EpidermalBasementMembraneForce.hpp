#ifndef EPIDERMALBASEMENTMEMBRANEFORCE_HPP_
#define EPIDERMALBASEMENTMEMBRANEFORCE_HPP_

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
 * A force class that defines the force due to the basement membrane.
 */

class EpidermalBasementMembraneForce : public AbstractForce<2>
{
    friend class TestCrossSectionModelInteractionForce;

private :

    /** Parameter that multiplies the curvature to give the basement membrane force */
    double mBasementMembraneParameter;

    /** Target curvature for the ring of cells (NodeBased) */
    double mTargetCurvature;

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
        archive & mTargetCurvature;
        archive & mCutOffRadius;
    }

public :

    /**
     * Constructor.
     */
	EpidermalBasementMembraneForce();

    /**
     * Destructor.
     */
    ~EpidermalBasementMembraneForce();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Value of curvature at Epidermis base */
    void SetTargetCurvature(double targetCurvature = 0.0);

    /*
     * Get method for Target Curvature Parameters
     */
    double GetTargetCurvature();

    /* Returns Epidermis height extremes to apply non-zero curvature to base */
    c_vector<double, 2> GetEpidermisHeightExtremes(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> epidermalWidthExtremes);

    /* Returns Epidermis width extremes */
    c_vector<double, 2> GetEpidermisWidthExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /* Return tangent vector at point based on cosine parametrisation */
    c_vector<double, 2> GetCosineBasedTangentVector(AbstractCellPopulation<2>& rCellPopulation,
                                                    c_vector<double, 2> epidermalHeightExtremes,
                                                    c_vector<double, 2> epidermalWidthExtremes,
                                                    c_vector<double, 2> point);

    /*
     * Return vector of Epidermal indices that are close to the considered Epidermal node,
     * but based on an approximated cosine approximation
     * 
     * @param rCellPopulation the cell population
     * @param epidermalIndex the considered epidermal node index
     * @param leftOrRight, should be 1.0 or -1.0, determines whether we consider 'left' or 'right' neighbours
     */
    std::vector<unsigned> GetClosestNeighboursBasedOnCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, 
                                                                        std::vector<unsigned> epidermalIndices,
                                                                        c_vector<double, 2> epidermalHeightExtremes,
                                                                        c_vector<double, 2> epidermalWidthExtremes,
                                                                        unsigned epidermalIndex,
                                                                        double leftOrRight);

    /*
     * Return closest epidermal node indice to the considered Epidermal node,
     * 
     * @param rCellPopulation the cell population
     * @param epidermalIndex the considered epidermal node index
     * @param leftOrRight, should be 1.0 or -1.0, determines whether we consider 'left' or 'right' neighbours
     */
    unsigned GetNearestNeighbourAlongCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, 
                                                        std::vector<unsigned> epidermalIndices,
                                                        c_vector<double, 2> epidermalHeightExtremes,
                                                        c_vector<double, 2> epidermalWidthExtremes,
                                                        unsigned epidermalIndex,
                                                        double leftOrRight);

    /*
     * Return closest Collagen index to the considered epidermal index
     * 
     * @param rCellPopulation the cell population
     * @param epidermalIndex the considered epidermal node index
     */
    c_vector<unsigned,2> GetTwoNearestCollagenNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned epidermalIndex);

    /* Get method for cut off radius */
    double GetCutOffRadius();

    /* Set method for cut off radius */
    void SetCutOffRadius(double cutOffRadius);

    /* Method to calculate parametre */ 
    double FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
    								c_vector<double, 2> leftPoint,
									c_vector<double, 2> centrePoint,
									c_vector<double, 2> rightPoint);

    /*
     * Method to get the indices of the monolayer
     */
    std::vector<unsigned> GetEpidermalIndices(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to return the nodes connected to a particular node within a defined cut-off
     * radius
     */
    std::vector<unsigned> GetNeighbouringEpidermalIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /*
     * Method to calculate the force due to basement membrane on an Epidermal cell
     */
    c_vector<double, 2> CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, 
                                                            std::vector<unsigned> epidermalIndices,
                                                            c_vector<double, 2> epidermalHeightExtremes,
                                                            c_vector<double, 2> epidermalWidthExtremes,
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
CHASTE_CLASS_EXPORT(EpidermalBasementMembraneForce)

#endif /*EPIDERMALBASEMENTMEMBRANEFORCE_HPP_*/
