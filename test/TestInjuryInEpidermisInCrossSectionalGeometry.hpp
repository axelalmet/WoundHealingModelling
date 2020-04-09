#ifndef TESTINJURYINEPIDERMISINCROSSSECTIONALGEOMETRY_HPP_
#define TESTINJURYINEPIDERMISINCROSSSECTIONALGEOMETRY_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "GeneralisedLinearSpringForce.hpp" // Standard spring force that implements logarithmic repulsion and exponential attraction for OS models
#include "EpidermalBasementMembraneForce.hpp" // Force to anchor basal stem cells to dermis
#include "WoundBasedChemotacticForce.hpp" // Individual-based chemotactic force to induce migration.
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "PlaneBoundaryCondition.hpp" // Plane-based boundary condition
#include "HoneycombMeshGenerator.hpp" // Generates mesh
#include "NoCellCycleModel.hpp" // Useful for running tests where cell proliferation isn't needed.
#include "BasementMembraneBasedContactInhibitionCellCycleModel.hpp" // Cell cycle for epidermal cell, where proliferative capacity is determined by attachment to the basement membrane
#include "GrowthFactorBasedContactInhibitionCellCycleModel.hpp" // Cell cycle for fibroblasts that is dependent on exposure to wound-derived growth factors
#include "FibroblastStateDependentCollagenSrnModel.hpp"
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "RandomSymmetricDivisionBasedDivisionRule.hpp" // Symmetric and asymmetric division rule
#include "Cylindrical2dNodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "CellDataItemWriter.hpp" // Allows us to track different cell data items
#include "ParabolicGrowingDomainWithCellDeathPdeModifier.hpp" // Modifier to track PDE solutions
#include "PlateletDerivedGrowthFactorCellwiseSourceParabolicPde.hpp" // Cellwise-source-based PDE to simulate PDGF due to wound healing
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CellLabel.hpp" // What we use to mark cells along the bottom boundary
#include "StemCellProliferativeType.hpp" // Epidermal basal stem cell type
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "DifferentiatedCellProliferativeType.hpp" // Differentiated cell type
#include "PlateletCellProliferativeType.hpp" // Blood platelet cell type
#include "EnfFibroblastCellMutationState.hpp" // ENF fibroblast mutation state
#include "EpfFibroblastCellMutationState.hpp" // EPF fibroblast mutation state
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "PlateletCellMutationState.hpp" //Platelet cell mutation state
#include "BasementMembraneAttachmentTrackingModifier.hpp" // Modifier to track stem cell attachment to the basement membrane
#include "CollagenAlignmentTrackingModifier.hpp" // Modifier to align fibroblasts with local collagen fibre orientation
#include "VolumeTrackingModifier.hpp" // Modifier to track cell volume
#include "PlateletCellKiller.hpp" // Cell killer to remove platelets upon wound healing
#include "BasementMembraneDistanceBasedCellKiller.hpp" // Cell killer based on distance to the basement membrane
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_STEADY_STATE_DIRECTORY = "WoundHealingModel/CrossSection/Epidermis/SteadyState";
static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel/CrossSection";
static const double M_DT = 0.005;
static const double M_STEADY_STATE_TIME = 100.0;
static const double M_END_TIME = M_STEADY_STATE_TIME + 100.0;
static const double M_SAMPLING_TIMESTEP = 10.0/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into modelling the epidermis in wound healing. One of the main obstacles is modelling
* the effect of hte basement membrane in maintaining a basal layer of epidermal cells.
*/
class TestCrossSectionalEpidermis : public AbstractCellBasedTestSuite
{
public:
    void TestInjuryInEpidermis()
    {
        // Load the steady state directory
		OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(M_STEADY_STATE_DIRECTORY, 100.0);

        //Get the tissue geometry so we may correctly wound the model.
        double max_width = 0.0;
        double max_height = 0.0;

        for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
                cell_iter != p_simulator->rGetCellPopulation().End();
                ++cell_iter)
        {
            //Get location of cell
            double x = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[0];
            double y = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

            if (x > max_width)
            {
                max_width = x;
            }
            
            if (y > max_height)
            {
                max_height = y;
            }
        }

        // Now we can define the wound geometry
        double wound_centre = 0.5*max_width;
        double wound_width = 0.25*max_width;
        double wound_base_height = 0.375*max_height;


        // Kill the cells within the wound area
        for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
                cell_iter != p_simulator->rGetCellPopulation().End();
                ++cell_iter)
        {
            //Get location of cell
            double x = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[0];
            double y = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

            //If the cell is within the 'wound area', we kill it.
            if ( (x > (wound_centre - 0.5*wound_width))&&(x < (wound_centre + 0.5*wound_width))&&(y > wound_base_height) )
            {
                cell_iter->Kill();
            }
        
        }

        //Set output directory
        std::stringstream out;
        out << "/Epidermis/Injury/";
        std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
        p_simulator->SetOutputDirectory(output_directory);

        p_simulator->SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        p_simulator->SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        p_simulator->Solve(); // Run the simulation.

        delete p_simulator;

    }
};

#endif /* TESTINJURYINEPIDERMISINCROSSSECTIONALGEOMETRY */