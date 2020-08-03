#ifndef TESTEPIDERMISONCOLLAGENLATTICEINCROSSSECTIONALGEOMETRY_HPP_
#define TESTEPIDERMISONCOLLAGENLATTICEINCROSSSECTIONALGEOMETRY_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp" // Standard spring force that implements logarithmic repulsion and exponential attraction for OS models
#include "EpidermalBasementMembraneForce.hpp" // Force to anchor basal stem cells to dermis (based off Du et al. (2018) model)
#include "DistanceBasedEpidermalBasementMembraneForce.hpp" // Force to anchor basal stem cells to dermis (based off Dunn et al. (2012) model)
#include "WoundBasedChemotacticForce.hpp" // Individual-based chemotactic force to induce migration.
#include "FibreAlignmentBasedMigrationForce.hpp" // Fibre-alignment-based migration force
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "EcmBasedPlaneBoundaryCondition.hpp" // Plane boundary condition for ecm nodes only
#include "EcmBasedPlaneBoundaryCondition.hpp" // Plane boundary condition for ecm nodes only
#include "NoCellCycleModel.hpp" // Useful for running tests where cell proliferation isn't needed.
#include "BasementMembraneBasedContactInhibitionCellCycleModel.hpp" // Cell cycle for epidermal cell, where proliferative capacity is determined by attachment to the basement membrane
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "Cylindrical2dNodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "CellMigrationDirectionWriter.hpp" // Allows us to track migration direction
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CellLabel.hpp" // What we use to mark cells along the bottom boundary
#include "StemCellProliferativeType.hpp" // Epidermal basal stem cell type
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "CollagenCellProliferativeType.hpp" // Dermal cell type
#include "DifferentiatedCellProliferativeType.hpp" // Differentiated cell type
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "CellMigrationDirectionWriter.hpp" // Cell writer for migration direction
#include "BasementMembraneAttachmentTrackingModifier.hpp" // Modifier to track stem cell attachment to the basement membrane
#include "VolumeTrackingModifier.hpp" // Modifier to track cell volume
#include "CollagenFibreTrackingModifier.hpp" // Modifier to track collagen fibres, defined by marked springs
#include "PolarityTrackingModifier.hpp" // Modifier to update cell polarities
#include "SloughingCellKiller.hpp" // Cell killer that sloughs of cells past a certain height
#include "MathsCustomFunctions.hpp" // So that I can use SmallPow
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel/CrossSection/EpidermisOnCollagenLattice";
static const double M_DT = 0.005;
static const double M_END_TIME = 5.0;
static const double M_SS_SAMPLING_TIMESTEP = 0.1*M_END_TIME/M_DT;
// static const double M_INJURY_SAMPLING_TIMESTEP = 0.05*M_END_TIME/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into modelling the epidermis in wound healing. One of the main obstacles is modelling
* the effect of hte basement membrane in maintaining a basal layer of epidermal cells.
*/
class TestCrossSectionalEpidermisOnCollagenLattice : public AbstractCellBasedTestSuite
{
private:

    // Generator function to assemble the nodes and cells.
    void GenerateNodesAndCells(MutableMesh<2, 2>& rMesh, std::vector<CellPtr>& rCells,
                                std::vector<unsigned>& rCollagenIndices, 
                                unsigned cellsAcross, unsigned cellsUp, 
                                double collagenRadius,
                                double reticularDermisHeight, double papillaryDermisHeight,
                                unsigned reticularFibroblastsAcross, unsigned reticularFibroblastsUp, 
                                unsigned papillaryFibroblastsAcross, unsigned papillaryFibroblastsUp, 
                                double fibroblastScale)
    {
        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>()); // Epidermal stem cell type
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>()); // Epidermal differentiated cell type
        boost::shared_ptr<AbstractCellProperty> p_fibroblast_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>()); // Dermal fibroblast cell
        boost::shared_ptr<AbstractCellProperty> p_collagen_type(CellPropertyRegistry::Instance()->Get<CollagenCellProliferativeType>()); // ECM collagen cell
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>()); // Generic WT mutation state
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>()); // Cell label to mark bottom row of cells

        double mesh_width = (double)cellsAcross;

        // Generate the epidermis first, taking care to make sure the bottom row of cells are no longer considered boundary notes
        for (unsigned i = 0; i < rMesh.GetNumNodes(); i++)
        {
            Node<2>* p_node = rMesh.GetNode(i); // Get the node
            // Get the node location
            double x = p_node->rGetLocation()[0];
            double y = p_node->rGetLocation()[1];
    
            // Adjust the interior nodes of the bottom row to not be boundary nodes
            if ( (x > 0.0)&&(x < mesh_width)&&(y == papillaryDermisHeight) )
            {
                p_node->SetAsBoundaryNode(false);
            }

            // Create the cell
            // NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
            BasementMembraneBasedContactInhibitionCellCycleModel* p_cycle_model = new BasementMembraneBasedContactInhibitionCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetDimension(2);

            // Initialise CellPtr
            CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));

            // Set a random birth time for each cell so that you don't get synchronised division.
            double birth_time = -12.0 * RandomNumberGenerator::Instance()->ranf();
            p_cell->SetBirthTime(birth_time);

            double migration_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf(); // Set a random migration direction

            if (y >= papillaryDermisHeight + 0.5)
            {
                p_cell->SetCellProliferativeType(p_diff_type); // Make cellstem type
                p_cell->GetCellData()->SetItem("attachment", 0.0);  // Initialise cell data to describe BM attachment.
            }
            else
            {
                p_cell->SetCellProliferativeType(p_stem_type); // Make cellstem type
                p_cell->GetCellData()->SetItem("attachment", 1.0);  // Initialise cell data to describe BM attachment.
            }

            p_cell->GetCellData()->SetItem("scale", 1.0); // Set the scale parameter needed to calculate the forces
            p_cell->GetCellData()->SetItem("volume", 0.25*M_PI); // Set the current volume 
            p_cell->GetCellData()->SetItem("direction", migration_direction);

            rCells.push_back(p_cell);
        }
    
        unsigned node_index = rMesh.GetNumNodes(); // Initialise node index counter, which is just the number of epidermal cells

        double horizontal_spacing = 2.0 * collagenRadius;
        double vertical_spacing = collagenRadius * sqrt(3.0);

        unsigned collagen_cells_across = (unsigned)(mesh_width/horizontal_spacing);
        unsigned collagen_cells_up = (unsigned)(papillaryDermisHeight/vertical_spacing) + 1;

        bool is_boundary = false;

        for (unsigned i = 0; i < collagen_cells_up; i++)
        {
            for (unsigned j = 0; j < collagen_cells_across; j++)
            {
                // This essentially is borrowed from HoneycombMeshGenerator, so that we can create a triangular lattice
                double x = horizontal_spacing * ((double)j + 0.25 * (1.0 + SmallPow(-1.0, i + 1)));
                double y = vertical_spacing * (double)i;

                // Mark the appropriate nodes as boundary nodes
                if ( (i == 0)||(i == collagen_cells_up - 1)||(j == 0)||(j == collagen_cells_across - 1))
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false;
                }

                rCollagenIndices.push_back(node_index);

                // Create the node
                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y));
                p_node->SetRadius(collagenRadius);

                rMesh.AddNode(p_node); // Add the node
                node_index += 1;
                
                // Now create the cell
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Place-holder cell cycle model
                p_cycle_model->SetDimension(2);

                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                p_cell->GetCellData()->SetItem("volume", M_PI * collagenRadius * collagenRadius);
                p_cell->GetCellData()->SetItem("attachment", -1.0);
                p_cell->GetCellData()->SetItem("direction", 0.0);
                p_cell->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

                rCells.push_back(p_cell); // Add the cell

            }
        }

        rMesh.ReMesh(); // Re-mesh, as we've added several nodes
       
        // We now populate the dermis with fibroblasts.
        // Fibroblasts are more scarce in the reticular dermis.
        double reticular_horizontal_spacing = mesh_width/( (double)reticularFibroblastsAcross); // Rough horizontal spacing
        double reticular_vertical_spacing = reticularDermisHeight/( (double)reticularFibroblastsUp); // Rough vertical spacing

        // We fill in the dermis by breaking up each layer into 'cells', and randomly place a 
        // fibroblast in each of the cells, so that fibroblasts are distributed more evenly
        for (unsigned i = 0; i < reticularFibroblastsUp; i++)
        {
            for (unsigned j = 0; j < reticularFibroblastsAcross; j++)
            {
                double x = (double)j * reticular_horizontal_spacing + reticular_horizontal_spacing * RandomNumberGenerator::Instance()->ranf();
                double y = 0.5 + (double)i * reticular_vertical_spacing + reticular_vertical_spacing * RandomNumberGenerator::Instance()->ranf();

                if ( (x == 0.0)||(x == mesh_width) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false;
                }  

                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y));  // Define the node

                rMesh.AddNode(p_node); // Add the node
                node_index += 1;

                // Set contact inhibition based cell cycle
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                // p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
                // p_cycle_model->SetQuiescentVolumeFraction(0.8);
                // p_cycle_model->SetGrowthFactorThreshold(0.25);
                p_cycle_model->SetDimension(2);

                // Initialise CellPtr
                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type);
                
                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = -12.0 * RandomNumberGenerator::Instance()->ranf();
                p_cell->SetBirthTime(birth_time);

                // Generate hte migration direction
                double migration_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

                // Add all the relevant cell data
                p_cell->GetCellData()->SetItem("scale", fibroblastScale); // Set the scale parameter needed to calculate the forces
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI); // Mature cell volume
                p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                p_cell->GetCellData()->SetItem("direction", migration_direction); // Migration direciton

                rCells.push_back(p_cell); // Add the cell
            }
        }

        rMesh.ReMesh(); // Re-mesh, as we've added several nodes

        // Papillary dermis contains more fibroblasts
        double papillary_horizontal_spacing = mesh_width/( (double)papillaryFibroblastsAcross); 
        double papillary_vertical_spacing = (papillaryDermisHeight - reticularDermisHeight - 0.5)/( (double)papillaryFibroblastsUp);

        for (unsigned i = 0; i < papillaryFibroblastsUp; i++)
        {
            for (unsigned j = 0; j < papillaryFibroblastsAcross; j++)
            {
                double x = (double)j * papillary_horizontal_spacing + papillary_horizontal_spacing * RandomNumberGenerator::Instance()->ranf();
                double y = reticularDermisHeight + (double)i * papillary_vertical_spacing + papillary_vertical_spacing * RandomNumberGenerator::Instance()->ranf();

                if ( (x == 0.0)||(x == mesh_width) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false;
                } 

                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y));  // Define the node
                
                rMesh.AddNode(p_node); // Add the node
                node_index += 1;

                // Set contact inhibition based cell cycle
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                // p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
                // p_cycle_model->SetQuiescentVolumeFraction(0.8);
                // p_cycle_model->SetGrowthFactorThreshold(0.25);
                p_cycle_model->SetDimension(2);

                // Initialise CellPtr
                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type);
                
                // Set a random birth time for each cell so that you don't get synchronised division.
                double birth_time = -12.0 * RandomNumberGenerator::Instance()->ranf();
                p_cell->SetBirthTime(birth_time);

                // Generate hte migration direction
                double migration_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

                // Add all the relevant cell data
                p_cell->GetCellData()->SetItem("scale", fibroblastScale); // Set the scale parameter needed to calculate the forces
                p_cell->GetCellData()->SetItem("volume", 0.25*M_PI); // Mature cell volume
                p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                p_cell->GetCellData()->SetItem("direction", migration_direction); // Migration direciton

                rCells.push_back(p_cell); // Add the cell
            }
        }

        rMesh.ReMesh(); // Re-mesh, as we've added several nodes

    }
    
    // Function to mark collagen fibres
    void MarkCollagenFibres(NodeBasedCellPopulation<2>& rCellPopulation, std::vector<unsigned>& rCollagenIndices, double collagenRadius, double crosslinkProbability)
    {
        std::random_shuffle(rCollagenIndices.begin(), rCollagenIndices.end()); // Shuffle the collagen indices

        for (unsigned i = 0; i < rCollagenIndices.size(); i++)
        {
            unsigned node_index = rCollagenIndices[i]; // Get the node index

            // Need to set the node radius
            Node<2>* p_node = rCellPopulation.GetNode(node_index);
            p_node->SetRadius(collagenRadius);

            c_vector<double, 2> node_location = p_node->rGetLocation(); // Will need the current node's location to check for periodicity

            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

            // Get the neighbours
            std::set<unsigned> neighbouring_nodes_set = rCellPopulation.GetNodesWithinNeighbourhoodRadius(node_index, 2.1 * collagenRadius);
            std::vector<unsigned> neighbouring_nodes;

            unsigned fibres_to_mark = 3;
            double probability = RandomNumberGenerator::Instance()->ranf();

            // If p < crossLinkProbability, we choose 3 neighbours. Otherwise, we choose 4.
            if (probability > crosslinkProbability)
            {
                fibres_to_mark += 1;
            }

            unsigned marked_fibres = 0; // Count we initialise for how many fibres we will mark;

            // First check if any fibres have been marked already
            for (std::set<unsigned>::iterator elem_iter = neighbouring_nodes_set.begin();
            elem_iter != neighbouring_nodes_set.end();
            ++elem_iter)
            {                
                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

                if (p_neighbour_cell->GetCellProliferativeType()->IsType<CollagenCellProliferativeType>() )
                {
                    // Chceck if the considered cell pair has been marked or not.
                    std::pair<CellPtr, CellPtr> cell_pair = rCellPopulation.CreateCellPair(p_cell, p_neighbour_cell);

                    if (rCellPopulation.IsMarkedSpring(cell_pair)) // If unmarked, mark it
                    {
                        marked_fibres += 1;
                    }
                    else
                    {
                        // We also need to check for periodicity. We do this by considering the absollute difference,
                        // rather than using GetVectorFromAToB, as GetNodesWithinHenighbourhoodRadius does.

                        c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(*elem_iter)->rGetLocation(); 

                        if (norm_2(neighbour_location - node_location) <= 2.1 * collagenRadius)
                        {
                            neighbouring_nodes.push_back(*elem_iter);
                        }
                    }

                    if (marked_fibres == fibres_to_mark)
                    {
                        break;
                }
                }
            }

            if ( (marked_fibres < fibres_to_mark)&&(neighbouring_nodes.size() > 0) ) // If there are still fibres to mark
            {
                // Randomly shuffle the indices
                std::random_shuffle(neighbouring_nodes.begin(), neighbouring_nodes.end());

                for (unsigned j = 0; j < neighbouring_nodes.size(); j++)
                {
                    unsigned neighbour_index = neighbouring_nodes[j];
                    
                    CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);
                    
                    // Chceck if the considered cell pair has been marked or not.
                    std::pair<CellPtr, CellPtr> cell_pair = rCellPopulation.CreateCellPair(p_cell, p_neighbour_cell);

                    if (!rCellPopulation.IsMarkedSpring(cell_pair)) // If unmarked, mark it
                    {
                        rCellPopulation.MarkSpring(cell_pair);
                    }

                    marked_fibres += 1; // Update the counter

                    if (marked_fibres == fibres_to_mark)
                    {
                        break;
                    }
                }
            }
        }
    }

public:
    void TestEpidermisOnCollagenLattice()
    {
        //Set the number of cells across and down for the array
        unsigned cells_across = 12;
        unsigned epidermal_cell_layers = 3;
        
        // Set the number of reticular and papillary dermis layers
        unsigned reticular_dermis_layers = 5;
        unsigned papillary_dermis_layers = 3;

        // Need to define the total cell layers for later
        unsigned cells_up = epidermal_cell_layers + reticular_dermis_layers + papillary_dermis_layers;

        // From this we can define the relevant heights
        // double max_height = cells_up * 0.5 * sqrt(3.0);
        double reticular_dermis_height = (double)reticular_dermis_layers * 0.5 * sqrt(3.0); // Five rows will be RD
        double papillary_dermis_height = reticular_dermis_height + (double)papillary_dermis_layers * 0.5 * sqrt(3.0); // Next three rows will be PD

        // Some numbers so that we can add additional collagen fibres and fibroblasts
        unsigned papillary_fibroblasts_across = 4;
        unsigned papillary_fibroblasts_up = 2;
        unsigned reticular_fibroblasts_across  = 3;
        unsigned reticular_fibroblasts_up = 3;

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Mechanical parameters
        double spring_stiffness = 30.0; // Spring stiffness
        // double bm_stiffness = 5.0; // Basement membrane attachment strength
        // double target_curvature = 0.0; // Target curvature
        double reorientation_strength = 2.5;
        bool apply_cell_to_ecm_force = true;
        bool apply_force_on_marked_springs = true;

        // Set the stiffness multiplication factors
        std::map<unsigned, double> stiffness_multiplication_factors;

        stiffness_multiplication_factors[0] = 1.0; // Stem cells 
        stiffness_multiplication_factors[1] = 1.0; // TA cells 
        stiffness_multiplication_factors[2] = 1.0; // Differentiated cells 
        stiffness_multiplication_factors[3] = 1.0; // Fibroblasts cells 
        stiffness_multiplication_factors[5] = 0.01; // Collagen cells 


        // Shape scale of fibroblasts and collagen fibres
        double fibroblast_scale = 0.5;
        double collagen_radius = 0.5;
        double crosslink_probability = 3.2;

        std::vector<CellPtr> cells; // Vector of cells
        std::vector<unsigned> collagen_indices; // Indices of collagen cells that we will track
        HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh(); //Generate mesh
        c_vector<double, 2> translate = zero_vector<double>(2);
        translate[1] = 0.5;
        p_generating_mesh->Translate(translate);

        // Remove the nodes that are below the papillary dermis height
        for (MutableMesh<2, 2>::NodeIterator node_iter = p_generating_mesh->GetNodeIteratorBegin(); 
            node_iter != p_generating_mesh->GetNodeIteratorEnd(); 
            ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            double y = node_iter->rGetLocation()[1];

            if (y < papillary_dermis_height)
            {
                p_generating_mesh->DeleteNodePriorToReMesh(node_index);
            }
        }

        p_generating_mesh->ReMesh(); // Remesh, as we've removed several cells.

        GenerateNodesAndCells(*p_generating_mesh, cells, collagen_indices,
                                cells_across, cells_up, collagen_radius,
                                reticular_dermis_height, papillary_dermis_height,
                                reticular_fibroblasts_across, reticular_fibroblasts_up, 
                                papillary_fibroblasts_across, papillary_fibroblasts_up, 
                                fibroblast_scale); // Generate the initial tissue

        Cylindrical2dNodesOnlyMesh* p_cylindrical_mesh = new Cylindrical2dNodesOnlyMesh(1.0*cells_across);
		p_cylindrical_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0); //Construct mesh;

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_cylindrical_mesh, cells); // Used for periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        cell_population.Update();

        // Mark the collagen fibres now
        MarkCollagenFibres(cell_population, collagen_indices, collagen_radius, crosslink_probability);

        // Add cell writers
        cell_population.AddCellWriter<CellMigrationDirectionWriter>(); // Visualise cell migration directions
 
        OffLatticeSimulation<2> simulator(cell_population);

        // Set output directory
        std::string output_directory = M_OUTPUT_DIRECTORY + "/SteadyState";
        simulator.SetOutputDirectory(M_OUTPUT_DIRECTORY);

        simulator.SetDt(M_DT);
        simulator.SetSamplingTimestepMultiple(M_SS_SAMPLING_TIMESTEP); // Set the sampling frequency
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        MAKE_PTR(GeneralisedLinearSpringForceWithVariableInteractionDistance<2>, p_spring_force);
        p_spring_force->SetSpringStiffness(spring_stiffness);
        p_spring_force->SetStiffnessMultiplicationFactors(stiffness_multiplication_factors);
        p_spring_force->ApplyCellToEcmForce(apply_cell_to_ecm_force);
        p_spring_force->ApplyForceOnMarkedSprings(apply_force_on_marked_springs);
        p_spring_force->SetCutOffLength(radius_of_interaction);
        simulator.AddForce(p_spring_force);

        // // // Add basement membrane force
        // // MAKE_PTR(DistanceBasedEpidermalBasementMembraneForce, p_bm_force);
        // MAKE_PTR(EpidermalBasementMembraneForce, p_bm_force);
        // p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        // p_bm_force->SetTargetCurvature(target_curvature);
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

        // Left boundary
        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(EcmBasedPlaneBoundaryCondition<2>, p_bc_left, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_left);

        // Bottom boundary
        point(0) = (double)cells_across;
        point(1) = 0.0;
        normal(0) = 1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(EcmBasedPlaneBoundaryCondition<2>, p_bc_right, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_right);

        // // Create a modifier to track which cells are attached to the basement membrane.
        // MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		// simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // // // // // // Create a modifier to track which cells are attached to the basement membrane.
        // MAKE_PTR(BasementMembraneAttachmentTrackingModifier<2>, p_bm_attachment_tracking_modifier);
        // p_bm_attachment_tracking_modifier->SetNeighbourhoodRadius(radius_of_interaction);
		// simulator.AddSimulationModifier(p_bm_attachment_tracking_modifier);

        // double max_height = 0.5 * sqrt(3.0) * (double)cells_up;

        // // Add a cell killer to remove differentiated cells that are too far from the basement membrane.
        // MAKE_PTR_ARGS(SloughingCellKiller<2>, p_sloughing_killer, (&cell_population, max_height + 0.5));
        // simulator.AddCellKiller(p_sloughing_killer);

        // Modifier to track the collagen fibres
        MAKE_PTR(CollagenFibreTrackingModifier<2>, p_collagen_fibre_tracking_modifier);
        simulator.AddSimulationModifier(p_collagen_fibre_tracking_modifier);

        // Modifier to align the cell polarities with their velocities
        MAKE_PTR(PolarityTrackingModifier<2>, p_polarity_tracking_modifier);
        p_polarity_tracking_modifier->SetReorientationStrength(reorientation_strength);
        simulator.AddSimulationModifier(p_polarity_tracking_modifier);

        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTEPIDERMISINCROSSSECTIONALGEOMETRY */