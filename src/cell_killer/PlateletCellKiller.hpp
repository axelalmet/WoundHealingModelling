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

    //Cut off radius for NodeBasedCellPopulations
    double mCutOffRadius;

    // Mean death time for platelets when subject to apoptosis
    double mMeanDeathTime;

    // Threshold volume under which to induce apoptosis
    double mVolumeThreshold;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by Platelet
    out_stream mCellKillerOutputFile;

    std::string mOutputDirectory;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        archive & mCutOffRadius;
        archive & mMeanDeathTime;
        archive & mVolumeThreshold;
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
     * @return mDeathTime
     */
    double GetMeanDeathTime();

    /*
     * Method to set the mean death time for platelets when subject to apoptosis
     */
    void SetMeanDeathTime(double meanDeathTime);

        /*
     * @return mVolumeThreshold for cell killing
     */
    double GetVolumeThreshold();

    /*
     * Method to set the volume threshold that determines platelet cell degradation
     */
    void SetVolumeThreshold(double volumeThreshold);

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

    bool ShouldCellBeRemoved(CellPtr pCell);

    /**
     *  Loops over and kills cells by Platelet or at the orifice if instructed.
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

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
