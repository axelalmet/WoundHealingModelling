#ifndef TESTSCARFORMATIONWITHOUTEPIDERMIS_HPP_
#define TESTSCARFORMATIONWITHOUTEPIDERMIS_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp" // Standard spring force that implements logarithmic repulsion and exponential attraction for OS models
#include "WoundBasedChemotacticForce.hpp" // Individual-based chemotactic force to induce migration.
#include "FibreAlignmentBasedMigrationForce.hpp" // Fibre-alignment-based migration force
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "EcmBasedPlaneBoundaryCondition.hpp" // Plane boundary condition that only applies to ECM cells
#include "NoCellCycleModel.hpp" // Useful for running tests where cell proliferation isn't needed.
#include "FibroblastStateDependentCollagenSrnModel.hpp" // SRN model that will model collagen production/degradation (when activated)
#include "GrowthFactorBasedContactInhibitionCellCycleModel.hpp" 
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "Cylindrical2dNodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "CellMigrationDirectionWriter.hpp" // Allows us to track migration direction
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CellLabel.hpp" // What we use to mark cells along the bottom boundary
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "ExtracellularMatrixCellProliferativeType.hpp" // Collagen (ECM) cell type
#include "FibrinCellMutationState.hpp" // Fibrin matrix mutation state
#include "CollagenCellMutationState.hpp" // Collagen mutation state
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "PlateletCellMutationState.hpp" // Platelet mutation state
#include "BloodCellProliferativeType.hpp" // Platelet cell type
#include "CellDataItemWriter.hpp" // Generic cell data item writer, in case we need it
#include "CellMigrationDirectionWriter.hpp" // Cell writer for migration direction
#include "CellMutationStatesWriter.hpp" // Cell writer for migration direction
#include "BasementMembraneAttachmentTrackingModifier.hpp" // Modifier to track stem cell attachment to the basement membrane
#include "ModifiedVolumeTrackingModifier.hpp" // Modifier to track cell volume
#include "CollagenFibreTrackingModifier.hpp" // Modifier to track collagen fibres, defined by marked springs
#include "PolarityTrackingModifier.hpp" // Modifier to update cell polarities
#include "SloughingCellKiller.hpp" // Cell killer that sloughs of cells past a certain height
#include "PlateletCellKiller.hpp" // Cell killer that removes platelets once fibroblasts reach them
#include "ModifiedParabolicBoxDomainPdeModifier.hpp" // Modifier to track PDE solutions in box domain
#include "PlateletDerivedGrowthFactorAveragedSourceParabolicPde.hpp" // Averaged-source-based PDE to simulate PDGF due to wound healing
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel/CrossSection/ScarFormationWithoutEpidermis";
static const double M_DT = 0.005;
static const double M_END_TIME = 10.0;
static const double M_SAMPLING_TIMESTEP = 0.1*M_END_TIME/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into modelling the epidermis in wound healing. One of the main obstacles is modelling
* the effect of hte basement membrane in maintaining a basal layer of epidermal cells.
*/
class TestScarFormationWithoutEpidermis : public AbstractCellBasedTestSuite
{
private:

    // Generator function to assemble the nodes and cells.
    void GenerateNodesAndCells(std::vector<Node<2>* >& rNodes, std::vector<CellPtr>& rCells,
                                std::vector<unsigned>& rCollagenIndices, double maxHeight,
                                unsigned cellsAcross, double collagenRadius,
                                double collagenProductionRate, double collagenDegradationRate,
                                double reticularDermisHeight, double papillaryDermisHeight,
                                unsigned reticularFibroblastsAcross, unsigned reticularFibroblastsUp, 
                                unsigned papillaryFibroblastsAcross, unsigned papillaryFibroblastsUp, 
                                double fibroblastScale, double jiggleRadius,
                                double morphogenThreshold, double quiescentVolumeFaction,
                                double woundBaseCentre, double woundBaseWidth, double woundBaseHeight,
                                unsigned numBloodCells, double meanBloodCellDeathTime)
    {
        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_fibroblast_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>()); // Dermal fibroblast cell
        boost::shared_ptr<AbstractCellProperty> p_ecm_type(CellPropertyRegistry::Instance()->Get<ExtracellularMatrixCellProliferativeType>()); // ECM collagen cell
        boost::shared_ptr<AbstractCellProperty> p_blood_type(CellPropertyRegistry::Instance()->Get<BloodCellProliferativeType>()); // Leukocyte type
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>()); // Generic WT mutation state
        boost::shared_ptr<AbstractCellProperty> p_collagen_state(CellPropertyRegistry::Instance()->Get<CollagenCellMutationState>()); // ECM collagen mutation state
        boost::shared_ptr<AbstractCellProperty> p_fibrin_state(CellPropertyRegistry::Instance()->Get<FibrinCellMutationState>()); // ECM fibrin mutation state
        boost::shared_ptr<AbstractCellProperty> p_platelet_state(CellPropertyRegistry::Instance()->Get<PlateletCellMutationState>()); // Platelet mutation state
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>()); // Cell label to mark bottom row of cells

        double mesh_width = (double)cellsAcross;

        unsigned node_index = 0; // Initialise node index counter, which is just the number of epidermal cells

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

                // We will perturb the collagen node positions as well
                if (i == 0) // Bottom row should only be perturbed horizontally
                {
                    x += -1.0 * jiggleRadius + 2.0 * jiggleRadius * RandomNumberGenerator::Instance()->ranf();
                }
                else if ( (j == 0)||(j == collagen_cells_across - 1) ) // The Left and right boundaries should only be perturbed vertically
                {
                    y += -1.0 * jiggleRadius + 2.0 * jiggleRadius * RandomNumberGenerator::Instance()->ranf();
                }
                else // Otherwise perturb them with a random angle and radius
                {
                    double perturbation_angle = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();
                    double perturbation_radius = jiggleRadius * RandomNumberGenerator::Instance()->ranf();

                    x += perturbation_radius * cos(perturbation_angle);
                    y += perturbation_radius * sin(perturbation_angle);
                }

                rCollagenIndices.push_back(node_index);

                // Create the node
                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y));
                p_node->SetRadius(collagenRadius);

                rNodes.push_back(p_node); // Add the node
                node_index += 1;
                
                // Now create the cell
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Place-holder cell cycle model
                p_cycle_model->SetDimension(2);

                // Set the SRN model
                FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
                p_srn_model->SetOdeParameters(collagenProductionRate, collagenDegradationRate);

                CellPtr p_cell(new Cell(p_collagen_state, p_cycle_model, p_srn_model));

                // If we're in the wounded area, we will set the mutation state to be a fibrin
                p_cell->SetCellProliferativeType(p_ecm_type); // Set cell to be collagen cells

                p_cell->GetCellData()->SetItem("volume", M_PI * collagenRadius * collagenRadius);
                p_cell->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)
                p_cell->GetCellData()->SetItem("morphogen", 0.1); // Morphogen amount for PDE
                p_cell->GetCellData()->SetItem("activated", 0.0); // Activation status of cell (for ECM cells mainly)
                p_cell->GetCellData()->SetItem("activation time", 0.0); // Activation time of cell (for fibroblast cells mainly)

                // Set the fibre orientation and density
                double theta;

                double density = 0.5 * RandomNumberGenerator::Instance()->ranf(); // Randomly initialise density

                // If we're in the wound area  
                if ( (x >= woundBaseCentre - woundBaseWidth)&&(x <= woundBaseCentre + woundBaseWidth)&&(y >= woundBaseHeight) )
                {
                    p_cell->SetMutationState(p_fibrin_state); // First change the mutation state to a fibrin cell

                    theta = 0.25 * M_PI * RandomNumberGenerator::Instance()->ranf();
                }
                else
                {
                    if (y > reticularDermisHeight)
                    {
                        theta = pow(-1.0, (double)j) * (0.2 + 0.1 * RandomNumberGenerator::Instance()->ranf()) * M_PI;

                        if ( (i == collagen_cells_up - 1)&&(theta > 0.0) )
                        {
                            theta = -0.25 * M_PI * RandomNumberGenerator::Instance()->ranf();
                        }
                    }
                    else
                    {
                        theta = -0.2 + 0.4 * RandomNumberGenerator::Instance()->ranf();
                        density += 0.5;
                    }
                }

                p_cell->GetCellData()->SetItem("direction", theta);
                p_cell->GetCellData()->SetItem("density", density); // Collagen amount

                rCells.push_back(p_cell); // Add the cell

            }
        }

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

                rNodes.push_back(p_node); // Add the node
                node_index += 1;

                if ( (x >= woundBaseCentre - woundBaseWidth)&&(x <= woundBaseCentre + woundBaseWidth)&&(y >= woundBaseHeight) )
                {
                    // Set contact inhibition based cell cycle
                    NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Growth-factor dependent cell cycle model
                    // p_cycle_model->SetStemCellG1Duration(7.0); // I want to speed division up a bit for this simulation
                    p_cycle_model->SetDimension(2);

                    // Set the SRN model
                    FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
                    p_srn_model->SetOdeParameters(0.0, 0.0);

                    CellPtr p_cell(new Cell(p_fibrin_state, p_cycle_model, p_srn_model)); // Initialise CellPtr
                    p_cell->SetCellProliferativeType(p_ecm_type);
                    
                    // Set a random birth time for each cell so that you don't get synchronised division.
                    double birth_time = -24.0 * RandomNumberGenerator::Instance()->ranf();
                    p_cell->SetBirthTime(birth_time);

                    // Generate hte migration direction
                    double migration_direction = 0.25 * M_PI * RandomNumberGenerator::Instance()->ranf();

                    // Add all the relevant cell data
                    p_cell->GetCellData()->SetItem("scale", 1.0); // Set the scale parameter needed to calculate the forces
                    p_cell->GetCellData()->SetItem("volume", 0.25 * M_PI); // Mature cell volume
                    p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                    p_cell->GetCellData()->SetItem("direction", migration_direction); // Migration direciton
                    p_cell->GetCellData()->SetItem("density", 0.0); // Collagen amount
                    p_cell->GetCellData()->SetItem("morphogen", 0.1); // Set the morphogen variable for the PDE
                    p_cell->GetCellData()->SetItem("activated", 0.0); // Activation status of cell (for ECM cells mainly)
                    p_cell->GetCellData()->SetItem("activation time", 0.0); // Activation time of cell (for fibroblast cells mainly)

                    rCells.push_back(p_cell); // Add the cell
                }
                else
                {
                    // Set contact inhibition based cell cycle
                    GrowthFactorBasedContactInhibitionCellCycleModel* p_cycle_model = new GrowthFactorBasedContactInhibitionCellCycleModel(); // Growth-factor dependent cell cycle model
                    // p_cycle_model->SetStemCellG1Duration(7.0); // I want to speed division up a bit for this simulation
                    p_cycle_model->SetEquilibriumVolume(0.5 * fibroblastScale * M_PI);
                    p_cycle_model->SetQuiescentVolumeFraction(quiescentVolumeFaction);
                    p_cycle_model->SetGrowthFactorThreshold(morphogenThreshold);
                    p_cycle_model->SetDimension(2);

                    // Set the SRN model
                    FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
                    p_srn_model->SetOdeParameters(collagenProductionRate, collagenDegradationRate);

                    CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model, p_srn_model)); // Initialise CellPtr
                    p_cell->SetCellProliferativeType(p_fibroblast_type);
                    
                    // Set a random birth time for each cell so that you don't get synchronised division.
                    double birth_time = -24.0 * RandomNumberGenerator::Instance()->ranf();
                    p_cell->SetBirthTime(birth_time);

                    // Generate hte migration direction
                    double migration_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

                    // Add all the relevant cell data
                    p_cell->GetCellData()->SetItem("scale", fibroblastScale); // Set the scale parameter needed to calculate the forces
                    p_cell->GetCellData()->SetItem("volume", 0.5 * 0.5 * fibroblastScale * M_PI); // Mature cell volume
                    p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                    p_cell->GetCellData()->SetItem("direction", migration_direction); // Migration direciton
                    p_cell->GetCellData()->SetItem("density", 0.0); // Collagen amount
                    p_cell->GetCellData()->SetItem("morphogen", 0.1); // Set the morphogen variable for the PDE
                    p_cell->GetCellData()->SetItem("activated", 0.0); // Activation status of cell (for ECM cells mainly)
                    p_cell->GetCellData()->SetItem("activation time", 0.0); // Activation time of cell (for fibroblast cells mainly)

                    rCells.push_back(p_cell); // Add the cell
                }
            }
        }

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
                
                rNodes.push_back(p_node); // Add the node
                node_index += 1;

                if ( (x >= woundBaseCentre - woundBaseWidth)&&(x <= woundBaseCentre + woundBaseWidth)&&(y >= woundBaseHeight) )
                {
                    // Set contact inhibition based cell cycle
                    NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Growth-factor dependent cell cycle model
                    // p_cycle_model->SetStemCellG1Duration(7.0); // I want to speed division up a bit for this simulation
                    p_cycle_model->SetDimension(2);

                    // Set the SRN model
                    FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
                    p_srn_model->SetOdeParameters(0.0, 0.0);

                    CellPtr p_cell(new Cell(p_fibrin_state, p_cycle_model, p_srn_model)); // Initialise CellPtr
                    p_cell->SetCellProliferativeType(p_ecm_type);
                    
                    // Set a random birth time for each cell so that you don't get synchronised division.
                    double birth_time = -24.0 * RandomNumberGenerator::Instance()->ranf();
                    p_cell->SetBirthTime(birth_time);

                    // Generate hte migration direction
                    double migration_direction = 0.25 * M_PI * RandomNumberGenerator::Instance()->ranf();

                    // Add all the relevant cell data
                    p_cell->GetCellData()->SetItem("scale", 1.0); // Set the scale parameter needed to calculate the forces
                    p_cell->GetCellData()->SetItem("volume", 0.25 * M_PI); // Mature cell volume
                    p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                    p_cell->GetCellData()->SetItem("direction", migration_direction); // Migration direciton
                    p_cell->GetCellData()->SetItem("density", 0.0); // Collagen amount
                    p_cell->GetCellData()->SetItem("morphogen", 0.1); // Set the morphogen variable for the PDE
                    p_cell->GetCellData()->SetItem("activated", 0.0); // Activation status of cell (for ECM cells mainly)
                    p_cell->GetCellData()->SetItem("activation time", 0.0); // Activation time of cell (for fibroblast cells mainly)

                    rCells.push_back(p_cell); // Add the cell
                }
                else
                {
                    // Set contact inhibition based cell cycle
                    GrowthFactorBasedContactInhibitionCellCycleModel* p_cycle_model = new GrowthFactorBasedContactInhibitionCellCycleModel(); // Growth-factor dependent cell cycle model
                    // p_cycle_model->SetStemCellG1Duration(7.0); // I want to speed division up a bit for this simulation
                    p_cycle_model->SetEquilibriumVolume(0.5 * fibroblastScale * M_PI);
                    p_cycle_model->SetQuiescentVolumeFraction(quiescentVolumeFaction);
                    p_cycle_model->SetGrowthFactorThreshold(morphogenThreshold);
                    p_cycle_model->SetDimension(2);

                    // Set the SRN model
                    FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
                    p_srn_model->SetOdeParameters(collagenProductionRate, collagenDegradationRate);

                    CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model, p_srn_model)); // Initialise CellPtr
                    p_cell->SetCellProliferativeType(p_fibroblast_type);
                    
                    // Set a random birth time for each cell so that you don't get synchronised division.
                    double birth_time = -24.0 * RandomNumberGenerator::Instance()->ranf();
                    p_cell->SetBirthTime(birth_time);

                    // Generate hte migration direction
                    double migration_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

                    // Add all the relevant cell data
                    p_cell->GetCellData()->SetItem("scale", fibroblastScale); // Set the scale parameter needed to calculate the forces
                    p_cell->GetCellData()->SetItem("volume", 0.5 * 0.5 * fibroblastScale * M_PI); // Mature cell volume
                    p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                    p_cell->GetCellData()->SetItem("direction", migration_direction); // Migration direciton
                    p_cell->GetCellData()->SetItem("density", 0.0); // Collagen amount
                    p_cell->GetCellData()->SetItem("morphogen", 0.1); // Set the morphogen variable for the PDE
                    p_cell->GetCellData()->SetItem("activated", 0.0); // Activation status of cell (for ECM cells mainly)
                    p_cell->GetCellData()->SetItem("activation time", 0.0); // Activation time of cell (for fibroblast cells mainly)

                    rCells.push_back(p_cell); // Add the cell
                }
            }
        }

        for (unsigned i = 0; i < numBloodCells; i++)
        {
            // Randomly place the node in the wound area.
            double x = woundBaseCentre - woundBaseWidth + 2.0 * woundBaseWidth * RandomNumberGenerator::Instance()->ranf();
            double y = woundBaseHeight + (maxHeight - woundBaseHeight) * RandomNumberGenerator::Instance()->ranf();

            Node<2>* p_node(new Node<2>(node_index, false, x, y));  // Define the node

            rNodes.push_back(p_node); // Add the node
            node_index += 1;

            // Create the cell
            FibroblastStateDependentCollagenSrnModel* p_srn_model = new FibroblastStateDependentCollagenSrnModel(); //Fibroblast-state-dependent collagen SRN model
            p_srn_model->SetOdeParameters(0.0, 0.0);

            // Also need to change the cell cycle (if we converted a fibroblast)
            NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Place-holder cell cycle model
            p_cycle_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_platelet_state, p_cycle_model, p_srn_model)); // Initialise CellPtr
            p_cell->SetCellProliferativeType(p_blood_type);

            // Add all the relevant cell data
            p_cell->GetCellData()->SetItem("scale", 1.0); // Set the scale parameter needed to calculate the forces
            p_cell->GetCellData()->SetItem("volume", 0.25 * M_PI); // Mature cell volume
            p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
            p_cell->GetCellData()->SetItem("direction", 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf()); // Migration direciton
            p_cell->GetCellData()->SetItem("density", 0.0); // Collagen amount
            p_cell->GetCellData()->SetItem("morphogen", 0.1); // Set the morphogen variable for the PDE
            p_cell->GetCellData()->SetItem("activated", 0.0); // Activation status of cell (for ECM cells mainly)

            // Set the lifetime to follow an exponential distributions
            double activation_time = - meanBloodCellDeathTime * log(RandomNumberGenerator::Instance()->ranf());
            p_cell->GetCellData()->SetItem("activation time", activation_time); // Activation time of cell (for fibroblast cells mainly)
            
            rCells.push_back(p_cell);
        }
    }

public:
    void TestScarFormation()
    {
        //Set the number of cells across and down for the array
        unsigned cells_across = 20;
        
        // Set the number of reticular and papillary dermis layers
        unsigned reticular_dermis_layers = 5;
        unsigned papillary_dermis_layers = 3;

        // Need to define the total cell layers for later
        unsigned cells_up = reticular_dermis_layers + papillary_dermis_layers;
        double max_height = cells_up * 0.5 * sqrt(3.0);

        // Define the wound base height and width
        double wound_base = 0.3 * max_height;
        double wound_centre = 0.5 * (double)cells_across;
        double wound_width = 0.4 * wound_centre;

        // From this we can define the relevant heights
        double reticular_dermis_height = (double)reticular_dermis_layers * 0.5 * sqrt(3.0); // Five rows will be RD
        double papillary_dermis_height = reticular_dermis_height + (double)papillary_dermis_layers * 0.5 * sqrt(3.0); // Next three rows will be PD

        // Some numbers so that we can add additional collagen fibres and fibroblasts
        unsigned papillary_fibroblasts_across = 5;
        unsigned papillary_fibroblasts_up = 2;
        unsigned reticular_fibroblasts_across  = 4;
        unsigned reticular_fibroblasts_up = 3;

        unsigned num_leukocytes = 25; // Add some platelets to act as the source for PDEs

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Mechanical parameters
        double spring_stiffness = 30.0; // Spring stiffness
        double reorientation_strength = 2.5;
        bool apply_cell_to_ecm_force = false;
        bool apply_force_on_marked_springs = false;

        // Some extra parameters
        double chemotactic_force_strength = 2.0; // Strength of the chemotactic force
        double morphogen_threshold = 0.5; // Threshold of the morphogen to trigger proliferation
        double quiescent_fraction = 0.8; // Volume fraction to trigger contact inhibition

        // Collagen remodelling rates for fibres
        double collagen_production_rate = 1.0;
        double collagen_degradation_rate = 0.5;

        double mean_activation_lifetime = 10.0; // Mean fibroblast lifetime when no longer exposed to morphogen
        double mean_platelet_death_time = 5.0; // Mean death time for platelets

        // Note that the probability of fibre deposition essentially follows a Binomial distribution (you deposit or you don't), 
        // so we can backtrack from the expected number of fibres we want to deposit over [0, T] to work out the probability.
        double fibre_deposition_probability = M_DT/M_END_TIME * 5.0; // Probability of depositing ECM fibre. 

        // Set the stiffness multiplication factors
        std::map<unsigned, double> stiffness_multiplication_factors;

        stiffness_multiplication_factors[0] = 1.0; // Stem cells 
        stiffness_multiplication_factors[1] = 1.0; // TA cells 
        stiffness_multiplication_factors[2] = 1.0; // Differentiated cells 
        stiffness_multiplication_factors[3] = 1.0; // Fibroblasts cells 
        stiffness_multiplication_factors[4] = 0.1; // Leukocytes (platelets, for now)
        stiffness_multiplication_factors[5] = 0.0; // Collagen cells 
        stiffness_multiplication_factors[6] = 0.0; // Fibrin cells 

        // Shape scale of fibroblasts and collagen fibres
        double fibroblast_scale = 0.5;
        double collagen_radius = 0.5;
        double jiggle_radius = 0.1; // Radius to perturb the node positions for the fibres

        std::vector<Node<2>* > nodes; // Vector of nodes
        std::vector<CellPtr> cells; // Vector of cells
        std::vector<unsigned> collagen_indices; // Indices of collagen cells that we will track

        GenerateNodesAndCells(nodes, cells, collagen_indices, 
                                max_height, cells_across, collagen_radius,
                                collagen_production_rate, collagen_degradation_rate, 
                                reticular_dermis_height, papillary_dermis_height,
                                reticular_fibroblasts_across, reticular_fibroblasts_up, 
                                papillary_fibroblasts_across, papillary_fibroblasts_up, 
                                fibroblast_scale, jiggle_radius,
                                morphogen_threshold, quiescent_fraction,
                                wound_centre, wound_width, wound_base, 
                                num_leukocytes, mean_platelet_death_time); // Generate the initial tissue

        Cylindrical2dNodesOnlyMesh* p_cylindrical_mesh = new Cylindrical2dNodesOnlyMesh((double)cells_across);
		p_cylindrical_mesh->ConstructNodesWithoutMesh(nodes, 2.0); //Construct mesh;

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_cylindrical_mesh, cells); // Used for periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        // Add cell writers
        cell_population.AddCellWriter<CellMigrationDirectionWriter>(); // Visualise cell migration directions
        // cell_population.AddCellWriter<CellMutationStatesWriter>(); // Count the mutation states

        // // Add writer for morphogen so we can post-process the PDE field
        // boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        // cell_population.AddCellWriter(p_cell_data_item_writer);

        OffLatticeSimulation<2> simulator(cell_population);

        // Set output directory
        simulator.SetOutputDirectory(M_OUTPUT_DIRECTORY); // Output directory
        simulator.SetDt(M_DT); // Set dt
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); // Set the sampling frequency
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        MAKE_PTR(GeneralisedLinearSpringForceWithVariableInteractionDistance<2>, p_spring_force);
        p_spring_force->SetSpringStiffness(spring_stiffness);
        p_spring_force->SetStiffnessMultiplicationFactors(stiffness_multiplication_factors);
        p_spring_force->ApplyCellToEcmForce(apply_cell_to_ecm_force);
        p_spring_force->ApplyForceOnMarkedSprings(apply_force_on_marked_springs);
        p_spring_force->SetCutOffLength(radius_of_interaction);
        simulator.AddForce(p_spring_force);

        // Add the chemotactic force
        MAKE_PTR(WoundBasedChemotacticForce<2>, p_chemotactic_force);
        p_chemotactic_force->SetNeighbourhoodRadius(radius_of_interaction);
        p_chemotactic_force->SetChemotacticStrength(chemotactic_force_strength);
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

        // MAKE_PTR_ARGS(PlateletCellKiller, p_platelet_cell_killer, (&cell_population));
        // p_platelet_cell_killer->SetCutOffRadius(radius_of_interaction);
        // p_platelet_cell_killer->SetMeanDeathTime(mean_platelet_death_time);
        // p_platelet_cell_killer->SetVolumeThreshold(0.9*0.25*M_PI); // Platelet cells can be compressed to half their size before dying
        // simulator.AddCellKiller(p_platelet_cell_killer);

        // Modified to track cell volumes to apply contact inhibition
        MAKE_PTR(ModifiedVolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // Modifier to align the cell polarities with their velocities
        MAKE_PTR(PolarityTrackingModifier<2>, p_polarity_tracking_modifier);
        p_polarity_tracking_modifier->SetFibreReorientationStrength(reorientation_strength);
        p_polarity_tracking_modifier->SetVelocityReorientationStrength(reorientation_strength);
        p_polarity_tracking_modifier->SetNeighbourhoodRadius(radius_of_interaction);
        p_polarity_tracking_modifier->SetMorphogenThreshold(morphogen_threshold);
        p_polarity_tracking_modifier->SetMeanActivationLifetime(mean_activation_lifetime);
        p_polarity_tracking_modifier->SetFibreDepositionProbability(fibre_deposition_probability);
        p_polarity_tracking_modifier->SetCollagenProductionRate(collagen_production_rate);
        p_polarity_tracking_modifier->SetCollagenProductionRate(collagen_degradation_rate);
        simulator.AddSimulationModifier(p_polarity_tracking_modifier);

        // Define the reaction-diffusion PDE, using the value's from YangYang's paper.
        MAKE_PTR_ARGS(PlateletDerivedGrowthFactorAveragedSourceParabolicPde<2>, p_pde, (cell_population, 1.0, 0.36, 1.0, 0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_pde_bc, (0.0));

        // Define the box domain for the PDE
        ChastePoint<2> lower(-1.0, -1.0);
        ChastePoint<2> upper((double)(2 + cells_across), std::ceil(2.0 + max_height));
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_box_domain, (lower, upper));

        // Create a PDE Modifier object using this pde and bcs object
        MAKE_PTR_ARGS(ModifiedParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_pde_bc, true, p_box_domain));
        p_pde_modifier->SetDependentVariableName("morphogen");
        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve(); // Run the simulation.
    }
};

#endif /* TESTSCARFORMATIONWITHOUTEPIDERMIS_HPP_ */