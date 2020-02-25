#ifndef TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY_HPP_
#define TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY_HPP_

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
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "Cylindrical2dNodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "CellDataItemWriter.hpp" // Allows us to track different cell data items
#include "ParabolicGrowingDomainPdeModifier.hpp" // Modifier to track PDE solutions
#include "EpfFibroblastCollagenSourceParabolicPde.hpp" // Cellwise-source-based PDE to simulate collagen production due to EPF fibroblasts
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
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel";
static const double M_DT = 0.005;
static const double M_END_TIME = 1.0;
// static const double M_SAMPLING_TIMESTEP = M_END_TIME / M_DT;
static const double M_SAMPLING_TIMESTEP = 1.0/M_DT;

class TestCrossSectionalWoundHealing : public AbstractCellBasedTestSuite
{
public:
    void TestWounding()
    {

        //Set the number of cells across and down for the array
        unsigned cells_across = 20;
        unsigned cells_up = 10;

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Mechanical parameters
        double spring_stiffness = 50.0; // Spring stiffness
        // double bm_stiffness = 6.0; // Basement membrane attachment strength
        // double target_curvature = 0.0; // Target curvature

        // Set the probability of being an EPF fibroblast.
        double epf_fibroblast_probability = 0.5;

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
            p_cycle_model->SetQuiescentVolumeFraction(0.8);
            p_cycle_model->SetDimension(2);

            // Randomly fill the fibroblast population with EPF and ENF fibroblasts, according to proportions
            // from the Rinkevich et al. (2018) paper.
            double fibroblast_state = RandomNumberGenerator::Instance()->ranf();

            // Randomly initiate a collagen orientation
            double collagen_orientation = 2.0* M_PI * RandomNumberGenerator::Instance()->ranf();

            if (fibroblast_state < epf_fibroblast_probability) // Roughly in line with the Rinkevich et al. (2018) paper.
            {
                CellPtr p_cell(new Cell(p_epf_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated

                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = - RandomNumberGenerator::Instance()->ranf() * 0.0;
                p_cell->SetBirthTime(birth_time);

                // For completeness in the stupid contact inhibition model.
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);

                // Initialise morphogen for later.
                p_cell->GetCellData()->SetItem("morphogen", 1e-4);

                // Initialise cell data to describe BM attachment.
                p_cell->GetCellData()->SetItem("attachment", -1.0);

                // Initialise collagen orientation
                p_cell->GetCellData()->SetItem("orientation", collagen_orientation);

                cells.push_back(p_cell);
            }
            else
            {
                CellPtr p_cell(new Cell(p_enf_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated

                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = - RandomNumberGenerator::Instance()->ranf() * 0.0;
                p_cell->SetBirthTime(birth_time);

                // For completeness in the stupid contact inhibition model.
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);

                //Initialise morphogen for later.
                p_cell->GetCellData()->SetItem("morphogen", 1e-4);

                // Initialise cell data to describe BM attachment.
                p_cell->GetCellData()->SetItem("attachment", -1.0);

                // Initialise collagen orientation
                p_cell->GetCellData()->SetItem("orientation", collagen_orientation);

                cells.push_back(p_cell);
            }
        }

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells); // Used for periodic
        // NodeBasedCellPopulation<2> cell_population(mesh, cells); // Used for non-periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        // // Add a cell item writer to track the morphogen that we will encode to evolve due to wounding.
        // boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        // cell_population.AddCellWriter(p_cell_data_item_writer);

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
        {            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            
            // // Turn the 'upper' part of the tissue epidermis
            // if (y > (max_height - 1.25))
            // // if (y == max_height)
            // {
            //     //Change the cell cycle to one that is based on the basement membrane attachment
            //     BasementMembraneBasedContactInhibitionCellCycleModel* p_cycle_model = new BasementMembraneBasedContactInhibitionCellCycleModel(); //Contact-inhibition-based cycle model yet.
            //     p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
            //     p_cycle_model->SetQuiescentVolumeFraction(0.8);
            //     p_cycle_model->SetDimension(2);
            
            //     cell_iter->SetCellCycleModel(p_cycle_model);
            //     cell_iter->SetCellProliferativeType(p_stem_type);
            //     cell_iter->SetMutationState(p_wildtype_state);
            // }
            if (y == min_height)
            {
                cell_iter->AddCellProperty(p_cell_label);
            }

        }

        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory
        std::stringstream out;
        out << "/CrossSection/";
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

        // Add basement membrane force
        // MAKE_PTR(EpidermalBasementMembraneForce, p_bm_force);
        // p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        // p_bm_force->SetTargetCurvature(target_curvature);
        // p_bm_force->ApplyPeriodicForce(false);
        // simulator.AddForce(p_bm_force);

        // Define a fixed-regions boundary condition so that cells can't move past y = 0
        c_vector<double, 2> point, normal;

        // Bottom boundary
        point(0) = 0.0;
        point(1) = 0.25;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc_bottom, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_bottom);

        // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_bottom, (&cell_population, point, normal));
        // simulator.AddCellPopulationBoundaryCondition(p_bc_bottom);

        // Create a modifier to track which cells are attached to the basement membrane.
        // MAKE_PTR(BasementMembraneAttachmentTrackingModifier<2>, p_bm_attachment_tracking_modifier);
        // p_bm_attachment_tracking_modifier->SetNeighbourhoodRadius(radius_of_interaction);
		// simulator.AddSimulationModifier(p_bm_attachment_tracking_modifier);

        // Create a modifier to track which cells are attached to the basement membrane.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // // Define the reaction-diffusion PDE, using the value's from YangYang's paper.
        // MAKE_PTR_ARGS(EpfFibroblastCollagenSourceParabolicPde<2>, p_pde, (simulator.rGetCellPopulation(), 1.0, 0.3537, 1.0));
        // MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // // Create a PDE Modifier object using this pde and bcs object
        // MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true));
        // p_pde_modifier->SetDependentVariableName("morphogen");
        // simulator.AddSimulationModifier(p_pde_modifier);

        // Wound the model. 
        double wound_centre = 0.5*max_width;
        double wound_width = 0.5*max_width;
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
            {
                cell_iter->SetMutationState(p_platelet_state);
                cell_iter->SetCellProliferativeType(p_platelet_type);
                cell_iter->GetCellData()->SetItem("morphogen", 1.0);
            }
            else
            {
                cell_iter->GetCellData()->SetItem("morphogen", 0.0);
            }
        
        }

        // Add the platelet cell killer
        MAKE_PTR_ARGS(PlateletCellKiller, p_platelet_cell_killer, (&cell_population));
        p_platelet_cell_killer->SetCutOffRadius(radius_of_interaction);
        p_platelet_cell_killer->SetGrowthFactorThreshold(0.5);
        simulator.AddCellKiller(p_platelet_cell_killer);

        // Add a chemotactic force
        MAKE_PTR(WoundBasedChemotacticForce<2>, p_chemotactic_force);
        p_chemotactic_force->SetChemotacticStrength(2.0);
        p_chemotactic_force->SetNeighbourhoodRadius(radius_of_interaction);
        simulator.AddForce(p_chemotactic_force);

        simulator.Solve(); // Run the simulation.

        // Back-up code when we need to start saving these steady states for the wounding experiments
        // Save simulation in steady state
		// CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Wound the model. 
        // double wound_centre = 0.5*max_width;
        // double wound_width = 0.5*max_width;
        // double wound_base_height = 0.4*max_height;

        // //Obtain the proliferative cells
        // for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
        //         cell_iter != simulator.rGetCellPopulation().End();
        //         ++cell_iter)
        // {
        //     //Get location of cell
        //     double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[0];
        //     double y = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

        //     //If the cell is within the 'wound area', we kill it.
        //     if ( (x > (wound_centre - 0.5*wound_width))&&(x < (wound_centre + 0.5*wound_width))&&(y > wound_base_height) )
        //     {
        //         cell_iter->SetCellProliferativeType(p_diff_type);
        //         cell_iter->GetCellData()->SetItem("morphogen", 1.0);
        //     }

        // }

        // // Define the reaction-diffusion PDE, using the value's from YangYang's paper.
        // MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (simulator.rGetCellPopulation(), 1.0, 0.1*0.3537, -0.1));
        // MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // // Create a PDE Modifier object using this pde and bcs object
        // MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true));
        // p_pde_modifier->SetDependentVariableName("morphogen");
        // simulator.AddSimulationModifier(p_pde_modifier);

        // simulator.SetSamplingTimestepMultiple(0.25*M_SAMPLING_TIMESTEP);
        // simulator.SetEndTime(2.0*M_END_TIME);

        // simulator.Solve(); // Run the simulation again

        // // Tidying up
        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /* TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY */