#ifndef TESTCLOTFORMATIONINCROSSSECTIONALGEOMETRY_HPP_
#define TESTCLOTFORMATIONINCROSSSECTIONALGEOMETRY_HPP_

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
#include "CellMigrationDirectionWriter.hpp" // Cell writer for migration direction
#include "CellCollagenFibreOrientationWriter.hpp" // Cell writer for collagen fibre orientations
#include "PlateletCellKiller.hpp" // Cell killer to remove platelets upon wound healing
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel/CrossSection";
static const double M_DT = 0.005;
static const double M_END_TIME = 5.0;
// static const double M_SAMPLING_TIMESTEP = M_END_TIME / M_DT;
static const double M_SAMPLING_TIMESTEP = 0.1/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into the clot formation process that induces the majority of wound healing. Here,
* the main challenge is modelling the spread of growth factors induced by the clot
* and the subsequent degradation of the clot by activated fibroblasts.
*/
class TestCrossSectionalClotFormation : public AbstractCellBasedTestSuite
{
public:
    void TestClotFormation()
    {

        //Set the number of cells across and down for the array
        unsigned cells_across = 20;
        unsigned cells_up = 8;

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Mechanical parameters
        double spring_stiffness = 15.0; // Spring stiffness
        // double bm_stiffness = 6.0; // Basement membrane attachment strength
        // double target_curvature = 0.0; // Target curvature

        // Set the probability of being an EPF fibroblast.
        double epf_fibroblast_probability = 0.75;

        HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

        // Construct a periodic mesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(1.0*cells_across);
		p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0); //Construct mesh

        // Create a non-periodic mesh.
		// NodesOnlyMesh<2> mesh;
        // mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 2.0); //Construct mesh

        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_fibroblast_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>());        
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_enf_state(CellPropertyRegistry::Instance()->Get<EnfFibroblastCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_epf_state(CellPropertyRegistry::Instance()->Get<EpfFibroblastCellMutationState>());

        //Create tissue of cells. Initially we set them all to be differentiated
        std::vector<CellPtr> cells; //Create vector of cells

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++) // Iterator for periodic mesh
        // for (unsigned i = 0; i < mesh.GetNumNodes(); i++) // Iterator for non-periodic mesh
        {
            // Set contact inhibition based cell cycle
            GrowthFactorBasedContactInhibitionCellCycleModel* p_cycle_model = new GrowthFactorBasedContactInhibitionCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetGrowthFactorThreshold(1.0);
            p_cycle_model->SetDimension(2);

            // Randomly fill the fibroblast population with EPF and ENF fibroblasts, according to proportions
            // from the Rinkevich et al. (2018) paper.
            double fibroblast_state = RandomNumberGenerator::Instance()->ranf();

            // Set a random fibroblast direction
            double fibroblast_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

            if (fibroblast_state < epf_fibroblast_probability) // Roughly in line with the Rinkevich et al. (2018) paper.
            {
                CellPtr p_cell(new Cell(p_epf_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated

                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = - RandomNumberGenerator::Instance()->ranf() * 1.0;
                p_cell->SetBirthTime(birth_time);

                // For completeness in the stupid contact inhibition model.
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);

                // Initialise morphogen for later.
                p_cell->GetCellData()->SetItem("morphogen", 0.0);

                // Set EPF state
                p_cell->GetCellData()->SetItem("epf", 1.0);

                // Set collagen amount
                p_cell->GetCellData()->SetItem("collagen", 0.1);

                // Set fibroblast migration direction
                p_cell->GetCellData()->SetItem("direction", fibroblast_direction);

                cells.push_back(p_cell);
            }   
            else
            {
                CellPtr p_cell(new Cell(p_enf_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated

                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = - RandomNumberGenerator::Instance()->ranf() * 1.0;
                p_cell->SetBirthTime(birth_time);

                // For completeness in the stupid contact inhibition model.
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);

                //Initialise morphogen for later.
                p_cell->GetCellData()->SetItem("morphogen", 0.0);

                // Set EPF state
                p_cell->GetCellData()->SetItem("epf", 0.0);

                // Set collagen amount
                p_cell->GetCellData()->SetItem("collagen", 1e-2);

                // Set fibroblasti migration direction
                p_cell->GetCellData()->SetItem("direction", fibroblast_direction);

                cells.push_back(p_cell);
            }
        }

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells); // Used for periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        // Add cell writers
        cell_population.AddCellWriter<CellMigrationDirectionWriter>();
        // cell_population.AddCellWriter<CellCollagenFibreOrientationWriter>();

        //Get the maximum width so we know where to apply the right BC.
        double min_width = 0.0;
        double max_width = 0.0;
        double min_height = 0.0;
        double max_height = 0.0;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            if (x > max_width)
            {
                max_width = x;
            }
            else if (x < min_width)
            {
                min_width = x;
            }
            
            if (y > max_height)
            {
                max_height = y;
            }
            else if (y < min_height)
            {
                min_height = y;
            }
        }

        // Let's set the cell populations for the boundary cells and the epithelium

        boost::shared_ptr<AbstractCellProperty> p_cell_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        cell_iter != cell_population.End(); ++cell_iter)
        {           
            // double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            
            if (y > max_height - 1.5)
            {
                cell_iter->SetCellProliferativeType(p_stem_type);
                cell_iter->SetMutationState(p_wildtype_state);                       
            }
            if (y == min_height)
            {
                cell_iter->AddCellProperty(p_cell_label);
            }

        }


        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory
        std::stringstream out;
        out << "/ClotFormation/";
        std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
        simulator.SetOutputDirectory(output_directory);

        simulator.SetDt(M_DT);
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(radius_of_interaction);
        simulator.AddForce(p_spring_force);

        // Add the chemotactic force
        MAKE_PTR(WoundBasedChemotacticForce<2>, p_chemotactic_force);
        p_chemotactic_force->SetNeighbourhoodRadius(radius_of_interaction);
        p_chemotactic_force->SetChemotacticStrength(5.0*M_DT);
        simulator.AddForce(p_chemotactic_force);

        // Define a fixed-regions boundary condition so that cells can't move past y = 0
        c_vector<double, 2> point, normal;

        // Bottom boundary
        point(0) = 0.0;
        point(1) = 0.25;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc_bottom, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_bottom);

        // Create a modifier to track which cells are attached to the basement membrane.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // Define the reaction-diffusion PDE, using the value's from YangYang's paper.
        MAKE_PTR_ARGS(PlateletDerivedGrowthFactorCellwiseSourceParabolicPde<2>, p_pde, (simulator.rGetCellPopulation(), 1.0, 0.36, 1.0, 0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // Create a PDE Modifier object using this pde and bcs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainWithCellDeathPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true));
        p_pde_modifier->SetDependentVariableName("morphogen");
        simulator.AddSimulationModifier(p_pde_modifier);

        // Wound the model. 
        double wound_centre = 0.5*max_width;
        double wound_width = 0.25*max_width;
        double wound_base_height = 0.4*max_height;

        boost::shared_ptr<AbstractCellProperty> p_platelet_type(CellPropertyRegistry::Instance()->Get<PlateletCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_platelet_state(CellPropertyRegistry::Instance()->Get<PlateletCellMutationState>());

        //Obtain the proliferative cells
        for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                cell_iter != simulator.rGetCellPopulation().End();
                ++cell_iter)
        {
            //Get location of cell
            double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[0];
            double y = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

            //If the cell is within the 'wound area', we kill it.
            if ( (x > (wound_centre - 0.5*wound_width))&&(x < (wound_centre + 0.5*wound_width))&&(y > wound_base_height) )
            // if ( ( pow(x - wound_centre, 2.0) + pow(y - min_height, 2.0) < pow(wound_width, 2.0) ))
            {
                cell_iter->SetMutationState(p_platelet_state);
                cell_iter->SetCellProliferativeType(p_platelet_type);
                cell_iter->GetCellData()->SetItem("morphogen", 3.3);
            }
            else
            {
                cell_iter->GetCellData()->SetItem("morphogen", 0.0);
            }
        
        }

        // // Add the platelet cell killer
        MAKE_PTR_ARGS(PlateletCellKiller, p_platelet_cell_killer, (&cell_population));
        p_platelet_cell_killer->SetCutOffRadius(radius_of_interaction);
        p_platelet_cell_killer->SetGrowthFactorThreshold(0.1);
        p_platelet_cell_killer->SetVolumeThreshold(0.95*0.25*M_PI); // Platelet cells can be compressed to half their size before dying
        simulator.AddCellKiller(p_platelet_cell_killer);

        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTCLOTFORMATIONINCROSSSECTIONALGEOMETRY */