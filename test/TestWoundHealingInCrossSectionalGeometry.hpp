#ifndef TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY_HPP_
#define TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "GeneralisedLinearSpringForce.hpp" // Standard spring force that implements logarithmic repulsion and exponential attraction for OS models
#include "GeneralisedLinearSpringForceWithVariableCellCellStiffness.hpp" // Version of generalised linear spring force where the spring stiffness depends on cell-cell types
#include "DistanceBasedEpidermalBasementMembraneForce.hpp" // Force to anchor basal stem cells to dermis
#include "WoundBasedChemotacticForce.hpp" // Individual-based chemotactic force to induce migration.
#include "FibreAlignmentBasedMigrationForce.hpp" // Migration force that aligns each fibroblast with surrounding collagen fibres
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "BasementMembraneBasedContactInhibitionCellCycleModel.hpp" // Cell cycle for epidermal cell, where proliferative capacity is determined by attachment to the basement membrane
#include "GrowthFactorBasedContactInhibitionCellCycleModel.hpp" // Cell cycle for fibroblasts that is dependent on exposure to wound-derived growth factors
#include "FibroblastStateDependentCollagenSrnModel.hpp"
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "Cylindrical2dNodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "ParabolicGrowingDomainWithCellDeathPdeModifier.hpp" // Modifier to track PDE solutions
#include "ModifiedParabolicBoxDomainPdeModifier.hpp" // Modifier to track PDE solutions in box domain
#include "PlateletDerivedGrowthFactorAveragedSourceParabolicPde.hpp" // Averaged-source-based PDE to simulate PDGF due to wound healing
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
#include "BasementMembraneDistanceBasedCellKiller.hpp" // Random cell killer
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel/CrossSection";
static const double M_DT = 0.005;
static const double M_END_TIME = 168.0;
static const double M_SAMPLING_TIMESTEP = 12.0/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into modelling scar formation in wound healing. That is, how does the presence
* of collagen produced by EPF fibroblasts during wound healing affect the orientation
* of collagen and fibroblast behaviour?
*/
class TestCrossSectionalWoundHealing : public AbstractCellBasedTestSuite
{
public:
    void TestWoundHealing()
    {

        //Set the number of cells across and down for the array
        unsigned cells_across = 25;
        unsigned cells_up = 13;

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Mechanical parameters
        double spring_stiffness = 30.0; // Spring stiffness
        double bm_stiffness = 0.1; // Basement membrane attachment strength
        // double target_curvature = 0.0; // Target curvature

        // Reseed the random number generator
		RandomNumberGenerator::Instance()->Reseed(100);

        // Set the probability of being an EPF fibroblast.
        double epf_fibroblast_probability = 0.25;

        // Morphogen threshold for fibroblast proliferation and collagen activation
        double morphogen_threshold = 2.5;

        HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

        // Construct a periodic mesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(1.0*cells_across);
		p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 25.0); //Construct mesh

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
        {
            // Set contact inhibition based cell cycle
            GrowthFactorBasedContactInhibitionCellCycleModel* p_cycle_model = new GrowthFactorBasedContactInhibitionCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetGrowthFactorThreshold(morphogen_threshold);
            p_cycle_model->SetStemCellG1Duration(14.0);
            p_cycle_model->SetDimension(2);

            // Set collagen-based SRN model
            FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
            p_srn_model->SetMorphogenThreshold(morphogen_threshold);

            // Randomly fill the fibroblast population with EPF and ENF fibroblasts, according to proportions
            // from the Rinkevich et al. (2018) paper (0.75 EPF, 0.25 ENF)
            double fibroblast_state = RandomNumberGenerator::Instance()->ranf();

            // Randomly initiate a collagen orientation
            double collagen_orientation = -0.5*M_PI + (M_PI * RandomNumberGenerator::Instance()->ranf());

            // Random initiate a fibroblast direction
            double fibroblast_direction = 2.0*M_PI * RandomNumberGenerator::Instance()->ranf(); 

            // Randomly initiate collagen amount
            double collagen_amount = RandomNumberGenerator::Instance()->ranf();

            if (fibroblast_state < epf_fibroblast_probability) 
            {
                CellPtr p_cell(new Cell(p_epf_state, p_cycle_model, p_srn_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated
                // p_cell->InitialiseCellCycleModel();
                // p_cell->InitialiseSrnModel();

                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = - 10.0* RandomNumberGenerator::Instance()->ranf() * 1.0;
                p_cell->SetBirthTime(birth_time);

                // Initialise all the scalar cell data

                // For completeness in the contact inhibition model.
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);
                // Initialise morphogen for later.
                p_cell->GetCellData()->SetItem("morphogen", 0.1);

                // Initialise cell data to describe BM attachment.
                p_cell->GetCellData()->SetItem("attachment", -1.0);

                // Set EPF state
                p_cell->GetCellData()->SetItem("epf", 1.0);

                // Initialise fibroblast direction
                p_cell->GetCellData()->SetItem("direction", fibroblast_direction);

                // Initialise collagen orientation
                p_cell->GetCellData()->SetItem("orientation", collagen_orientation);

                // Set collagen amount
                p_cell->GetCellData()->SetItem("collagen", collagen_amount);

                // Set collagen amount
                p_cell->GetCellData()->SetItem("collagen", collagen_amount);

                cells.push_back(p_cell);
            }
            else
            {
                CellPtr p_cell(new Cell(p_enf_state, p_cycle_model, p_srn_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated
                // p_cell->InitialiseCellCycleModel();
                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = - 10.0 * RandomNumberGenerator::Instance()->ranf() * 1.0;
                p_cell->SetBirthTime(birth_time);

                // For completeness in the stupid contact inhibition model.
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);

                //Initialise morphogen for later.
                p_cell->GetCellData()->SetItem("morphogen", 0.1);

                // Initialise cell data to describe BM attachment.
                p_cell->GetCellData()->SetItem("attachment", -1.0);

                // Initialise fibroblast direction
                p_cell->GetCellData()->SetItem("direction", fibroblast_direction);

                // Initialise collagen orientation
                p_cell->GetCellData()->SetItem("orientation", collagen_orientation);

                // Set EPF state
                p_cell->GetCellData()->SetItem("epf", 0.0);

                // Set collagen amount
                p_cell->GetCellData()->SetItem("collagen", collagen_amount);

                cells.push_back(p_cell);
            }
        }

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells); // Used for non-periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        // Add cell writers
        cell_population.AddCellWriter<CellMigrationDirectionWriter>();
        cell_population.AddCellWriter<CellCollagenFibreOrientationWriter>();

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

        double max_fibroblast_height = 0.0;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        cell_iter != cell_population.End(); ++cell_iter)
        {           
            // double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            
            if (y == min_height)
            {
                cell_iter->AddCellProperty(p_cell_label);
            }
            if (y > max_height - 1.5*sqrt(3.0) - 0.1)
            {
                BasementMembraneBasedContactInhibitionCellCycleModel* p_cycle_model = new BasementMembraneBasedContactInhibitionCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
                p_cycle_model->SetQuiescentVolumeFraction(0.9);
                p_cycle_model->SetStemCellG1Duration(14.0);
                p_cycle_model->SetDimension(2);

                cell_iter->SetCellCycleModel(p_cycle_model);
                cell_iter->SetCellProliferativeType(p_stem_type);
                cell_iter->SetMutationState(p_wildtype_state);

                // Initialise cell data to describe BM attachment.
                cell_iter->GetCellData()->SetItem("attachment", 1.0);

                // Should turn of the EPF status
                cell_iter->GetCellData()->SetItem("epf", 0.0);

                // Should turn off collagen as well
                cell_iter->GetCellData()->SetItem("collagen", 0.0);

            }
            else
            {
                if (y > max_fibroblast_height)
                {
                    max_fibroblast_height = y;
                }
            }

        }

        // Wound the model. 
        double wound_centre = 0.5*max_width;
        double wound_width = 0.5*max_width;
        double wound_base_height = 0.3*max_height;

        boost::shared_ptr<AbstractCellProperty> p_platelet_type(CellPropertyRegistry::Instance()->Get<PlateletCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_platelet_state(CellPropertyRegistry::Instance()->Get<PlateletCellMutationState>());

        //Obtain the proliferative cells
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            //Get location of cell
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            //If the cell is within the 'wound area', we kill it.
            if ( (x > (wound_centre - 0.5*wound_width))&&(x < (wound_centre + 0.5*wound_width))&&(y > wound_base_height) )
            // if ( ( pow(x - wound_centre, 2.0) + pow(y - min_height, 2.0) < pow(wound_width, 2.0) ))
            {
                cell_iter->SetMutationState(p_platelet_state);
                cell_iter->SetCellProliferativeType(p_platelet_type);
                cell_iter->GetCellData()->SetItem("collagen", 0.0);
                cell_iter->GetCellData()->SetItem("epf", 0.0);
                
            }
        
        }

        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory
        std::stringstream out;
        out << "/FullModel/EPF_" << epf_fibroblast_probability << "/";
        std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
        simulator.SetOutputDirectory(output_directory);

        simulator.SetDt(M_DT);
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(radius_of_interaction);
        // p_spring_force->SetFibroblastFibroblastMultiplicationFactor(1.0);
        // p_spring_force->SetStemFibroblastMultiplicationFactor(2.0);
        // p_spring_force->SetStemStemMultiplicationFactor(2.0);
        // p_spring_force->SetStemDifferentiatedMultiplicationFactor(4.0);
        // p_spring_force->SetDifferentiatedDifferentiatedMultiplicationFactor(4.0);
        // p_spring_force->SetDifferentiatedFibroblastMultiplicationFactor(4.0);
        // p_spring_force->SetPlateletPlateletMultiplicationFactor(1.0);
        // p_spring_force->SetStemPlateletMultiplicationFactor(2.0);
        // p_spring_force->SetDifferentiatedPlateletMultiplicationFactor(4.0);
        // p_spring_force->SetFibroblastPlateletMultiplicationFactor(1.0);
        simulator.AddForce(p_spring_force);


        // Add basement membrane force
        MAKE_PTR(DistanceBasedEpidermalBasementMembraneForce, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        simulator.AddForce(p_bm_force);

        // Add the chemotactic force
        MAKE_PTR(WoundBasedChemotacticForce<2>, p_chemotactic_force);
        p_chemotactic_force->SetNeighbourhoodRadius(radius_of_interaction);
        p_chemotactic_force->SetChemotacticStrength(2.5*M_DT);
        simulator.AddForce(p_chemotactic_force);

        // // Add fibre-alignment-based migration force
        // MAKE_PTR(FibreAlignmentBasedMigrationForce<2>, p_migration_force);
        // p_migration_force->SetMigrationForceStrength(0.0);
        // p_migration_force->SetReorientationStrength(2.5*M_DT);
        // simulator.AddForce(p_migration_force);

        // Define a fixed-regions boundary condition so that cells can't move past y = 0
        c_vector<double, 2> point, normal;

        // Bottom boundary
        point(0) = 0.0;
        point(1) = 0.25;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc_bottom, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_bottom);

        // Define the reaction-diffusion PDE, using the value's from YangYang's paper.
        MAKE_PTR_ARGS(PlateletDerivedGrowthFactorAveragedSourceParabolicPde<2>, p_pde, (cell_population, 1.0, 0.36, 1.0, 0.1));
        // MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, 1.0, 0.36, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // Define the box domain for the PDE
        ChastePoint<2> lower(-1.0, -1.0);
        ChastePoint<2> upper(1.0*(cells_across + 1), 0.5*sqrt(3)*(cells_up + 2));
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_box_domain, (lower, upper));

        // Create a PDE Modifier object using this pde and bcs object
        // MAKE_PTR_ARGS(ParabolicGrowingDomainWithCellDeathPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true));
        MAKE_PTR_ARGS(ModifiedParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true, p_box_domain));
        p_pde_modifier->SetDependentVariableName("morphogen");
        simulator.AddSimulationModifier(p_pde_modifier);

        // // Create a modifier to track which cells are attached to the basement membrane.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // Create a modifier to realign cell orientations with collagen.
        MAKE_PTR(CollagenAlignmentTrackingModifier<2>, p_collagen_alignment_modifier);
        p_collagen_alignment_modifier->SetNeighbourhoodRadius(radius_of_interaction);
        p_collagen_alignment_modifier->SetReorientationStrength(1.0*M_DT);
		simulator.AddSimulationModifier(p_collagen_alignment_modifier);

        // Create a modifier to track which cells are attached to the basement membrane.
        MAKE_PTR(BasementMembraneAttachmentTrackingModifier<2>, p_bm_attachment_tracking_modifier);
        p_bm_attachment_tracking_modifier->SetNeighbourhoodRadius(radius_of_interaction);
        simulator.AddSimulationModifier(p_bm_attachment_tracking_modifier);

        // Add a cell killer to remove differentiated cells that are too far from the basement membrane.
        MAKE_PTR_ARGS(BasementMembraneDistanceBasedCellKiller, p_cell_killer, (&cell_population, max_height - max_fibroblast_height + 0.1, max_height));
        simulator.AddCellKiller(p_cell_killer);

        // // Add the platelet cell killer
        MAKE_PTR_ARGS(PlateletCellKiller, p_platelet_cell_killer, (&cell_population));
        p_platelet_cell_killer->SetCutOffRadius(radius_of_interaction);
        p_platelet_cell_killer->SetGrowthFactorThreshold(morphogen_threshold);
        p_platelet_cell_killer->SetVolumeThreshold(0.95*0.25*M_PI); // Platelet cells can be compressed to half their size before dying
        simulator.AddCellKiller(p_platelet_cell_killer);

        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY */