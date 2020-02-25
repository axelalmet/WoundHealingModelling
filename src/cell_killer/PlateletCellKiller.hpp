#ifndef PLATELETCELLKILLER_HPP_
#define PLATELETCELLKILLER_HPP_
#include "CheckpointArchiveTypes.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellKiller.hpp"

/*
 * Cell killer that removes platelet cells due to contact with
 * activated fibroblasts during wound healing.
 */

class PlateletCellKiller : public AbstractCellKiller<2>
{
private:

	// Number of cells removed by fibroboasts
	unsigned mCellsRemovedByFibroblasts;

    std::vector<c_vector<double,3> > mLocationsOfPlateletCells;

    //Cut off radius for NodeBasedCellPopulations
    double mCutOffRadius;

    // Growth factor threshold in order to induce anoikis
    double mGrowthFactorThreshold;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by Platelet
    out_stream mCellKillerOutputFile;

    std::string mOutputDirectory;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        archive & mCellsRemovedByFibroblasts;
        archive & mCutOffRadius;
        archive & mGrowthFactorThreshold;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     */
	PlateletCellKiller(AbstractCellPopulation<2>* pCellPopulation);

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    /*
     * @return mGrowthFactorThreshold for cell killing
     */
    double GetGrowthFactorThreshold();

    /*
     * Method to set the growth factor threshold that determines platelet cell degradation
     */
    void SetGrowthFactorThreshold(double growthFactorThreshold);

    /*
     * @return mCutOffRadius
     */
    double GetCutOffRadius();

    /*
     * Method to defin mCutOffRadius by
     * cutOffRadius
     */
    void SetCutOffRadius(double cutOffRadius);

    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    bool ShouldCellBeRemoved(unsigned nodeIndex);

    std::vector<c_vector<unsigned,2> > RemoveByFibroblasts();

    /**
     *  Loops over and kills cells by Platelet or at the orifice if instructed.
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /* After each event of cell killing in CheckAndLabelCellsForApoptosisOrDeath(), the information of whether to kill each cell
     * or not is passed to this method which then increments the member variables corresponding to the total number of cells
     * killed by Platelet or apoptosis through compression
     */
    void SetNumberCellsRemoved(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the total number of cells removed by Platelet ([0]) and by compression ([1])
     *
     */
    unsigned GetNumberCellsRemoved();

    /* Storing the x-locations of those epithelial cells that get removed by Platelet
     *
     */
    void SetLocationsOfCellsRemovedByFibroblasts(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the x-coordinates of those cells removed by Platelet
     *
     */
    std::vector<c_vector<double,3> > GetLocationsOfCellsRemovedByFibroblasts();

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
CHASTE_CLASS_EXPORT(PlateletCellKiller)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const PlateletCellKiller * t, const unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, PlateletCellKiller * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)PlateletCellKiller(p_cell_population);
        }
    }
}

#endif /* PlateletCELLKILLER_HPP_ */
