#ifndef TESTEPIDERMISINCROSSSECTIONALGEOMETRY_HPP_
#define TESTEPIDERMISINCROSSSECTIONALGEOMETRY_HPP_

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
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "WoundHealingModel/CrossSection/Epidermis";
static const double M_DT = 0.005;
static const double M_END_TIME = 1.0;
static const double M_SS_SAMPLING_TIMESTEP = 0.1*M_END_TIME/M_DT;
// static const double M_INJURY_SAMPLING_TIMESTEP = 0.05*M_END_TIME/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into modelling the epidermis in wound healing. One of the main obstacles is modelling
* the effect of hte basement membrane in maintaining a basal layer of epidermal cells.
*/
class TestCrossSectionalEpidermis : public AbstractCellBasedTestSuite
{
private:

    // Generator function to assemble the nodes and cells.
    void GenerateNodesAndCells(MutableMesh<2, 2>& rMesh, std::vector<CellPtr>& rCells,
                                std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> >& rMarkedPairs,
                                std::vector<c_vector<double, 3> >& rPairsAndLowerEndpoint,
                                unsigned cellsAcross, unsigned cellsUp, 
                                double reticularDermisHeight, double papillaryDermisHeight,
                                double collagenRadius, double collagenFibreLength,
                                unsigned numReticularFibres, unsigned numPapillaryFibres,
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
        bool is_boundary = false;
        unsigned start_index, end_index;

        // Generate the reticular dermis fibres, which are largely horizontal
        for (unsigned i = 0; i < numReticularFibres; i++)
        { 
            // Generate the start and endpoints of the fibre.
            double xStart = mesh_width * RandomNumberGenerator::Instance()->ranf();
            double yStart = reticularDermisHeight * RandomNumberGenerator::Instance()->ranf();

            // We want reticular fibres to be roughly horizontal
            double theta = -0.05 + 0.1 * M_PI * RandomNumberGenerator::Instance()->ranf();

            double xEnd = xStart + collagenFibreLength * cos(theta);
            double yEnd = yStart + collagenFibreLength * sin(theta);

            // Adjust the endpoints if they fall outside of the domain
            if ( (xEnd > mesh_width) ) // Falls outside the vertical boundaries
            {
                double phi = theta;
                theta = M_PI - phi;

                xEnd = xStart + collagenFibreLength * cos(theta);
            }
            if (yEnd < 0.0) // Falls below the mesh
            {
                theta *= -1.0; // Reflect theta from Q4 (likely) to Q1
                yEnd = yStart + collagenFibreLength * sin(theta);
            }

            if ( (xStart == 0.0)||(xStart == mesh_width)||(yStart == 0.0) )
            {
                is_boundary = true;
            }
            else
            {
                is_boundary = false;
            }
            
            start_index = node_index;

            // Add the start point to the vector of nodes and cells
            Node<2>* p_node_start(new Node<2>(start_index, is_boundary, xStart, yStart)); // Define the node
            p_node_start->SetRadius(collagenRadius);

            rMesh.AddNode(p_node_start); // Add the node
            node_index += 1;

            // Set placeholder cell cycle model (i.e. no proliferation)
            NoCellCycleModel* p_cycle_model_start = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model_start->SetDimension(2);

            CellPtr p_cell_start(new Cell(p_wildtype_state, p_cycle_model_start));
            p_cell_start->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
            p_cell_start->GetCellData()->SetItem("attachment", -1.0); // BM attachment
            p_cell_start->GetCellData()->SetItem("volume", M_PI * pow(collagenRadius, 2.0)); // Mature cell volume
            p_cell_start->GetCellData()->SetItem("direction", theta);
            p_cell_start->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

            if (yStart < 0.25) // If endpoint falls into 'bottom row' of tissue, we fix it via a cell label
            {
                p_cell_start->AddCellProperty(p_label);
            }

            rCells.push_back(p_cell_start); // Add the cell

            if ( (xEnd == 0.0)||(xEnd == mesh_width)||(yEnd == 0.0) )
            {
                is_boundary = true;
            }
            else
            {
                is_boundary = false;
            }

            // Add the endpoint to the vector of nodes and cells
            end_index = node_index;

            Node<2>* p_node_end(new Node<2>(end_index, is_boundary, xEnd, yEnd)); // Define the node
            p_node_end->SetRadius(collagenRadius);

            rMesh.AddNode(p_node_end); // Add the node
            node_index += 1;

            // Set placeholder cell cycle model (i.e. no proliferation)
            NoCellCycleModel* p_cycle_model_end = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model_end->SetDimension(2);

            CellPtr p_cell_end(new Cell(p_wildtype_state, p_cycle_model_end));
            p_cell_end->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
            p_cell_end->GetCellData()->SetItem("attachment", -1.0); // BM attachment
            p_cell_end->GetCellData()->SetItem("volume", M_PI * pow(collagenRadius, 2.0)); // Mature cell volume
            p_cell_end->GetCellData()->SetItem("direction", theta);
            p_cell_end->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

            if (yEnd < 0.25) // If endpoint falls into 'bottom row' of tissue, we fix it via a cell label
            {
                p_cell_end->AddCellProperty(p_label);
            }

            rCells.push_back(p_cell_end); // Add the cell

            // We will now store this index pair for future use when we mark cross-links
            std::pair<unsigned, unsigned> node_pair = std::make_pair(start_index, end_index); 
            rMarkedPairs[node_pair] = {start_index, end_index}; 

            c_vector<double, 3> fibre_and_lower_y_coord;
            fibre_and_lower_y_coord[0] = (double)start_index;
            fibre_and_lower_y_coord[1] = (double)end_index;
            fibre_and_lower_y_coord[2] = std::min(yStart, yEnd);

            rPairsAndLowerEndpoint.push_back(fibre_and_lower_y_coord);
        }
        
        rMesh.ReMesh(); // Re-mesh, as we've added several nodes
                
        // Generate papillar dermis fibres, which are more criss-crossed
        for (unsigned i = 0; i < numPapillaryFibres; i++)
        {
            // Generate the start and endpoints of the fibre.
            double xStart = mesh_width * RandomNumberGenerator::Instance()->ranf();
            double yStart = reticularDermisHeight 
                            + (papillaryDermisHeight - reticularDermisHeight) * RandomNumberGenerator::Instance()->ranf();

            // We want papillary fibres to be criss-crossed.
            double theta = (0.225 + 0.05 * RandomNumberGenerator::Instance()->ranf())* M_PI + 
                            0.5 * M_PI + std::round(RandomNumberGenerator::Instance()->ranf());

            double xEnd = xStart + collagenFibreLength * cos(theta);
            double yEnd = yStart + collagenFibreLength * sin(theta);

            // Adjust the endpoints if they fall outside of the domain
            if ( (xEnd < 0.0)||(xEnd > mesh_width) ) // Falls outside the vertical boundaries
            {
                theta = M_PI - theta;
                xEnd = xStart + collagenFibreLength * cos(theta);
            }
            if (yEnd > papillaryDermisHeight) // Falls above the epidermis
            {
                if (theta < 0.5*M_PI) // Q1
                {   
                    double phi = theta;
                    theta = 2.0*M_PI - phi;
                }
                else
                {
                    double phi = M_PI - theta;
                    theta = M_PI + phi;
                }

                yEnd = yStart + collagenFibreLength * sin(theta);
            }

            if ( (xStart == 0.0)||(xStart == mesh_width)||(yStart == 0.0) )
            {
                is_boundary = true;
            }
            else
            {
                is_boundary = false;
            }
            
            start_index = node_index;

            // Add the start point to the vector of nodes and cells
            Node<2>* p_node_start(new Node<2>(start_index, is_boundary, xStart, yStart)); // Define the node
            p_node_start->SetRadius(collagenRadius);

            rMesh.AddNode(p_node_start); // Add the node
            node_index += 1;

            // Set placeholder cell cycle model (i.e. no proliferation)
            NoCellCycleModel* p_cycle_model_start = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model_start->SetDimension(2);

            CellPtr p_cell_start(new Cell(p_wildtype_state, p_cycle_model_start));
            p_cell_start->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
            p_cell_start->GetCellData()->SetItem("volume", M_PI * pow(collagenRadius, 2.0)); // Mature cell volume
            p_cell_start->GetCellData()->SetItem("attachment", -1.0); // BM attachment
            p_cell_start->GetCellData()->SetItem("direction", theta);
            p_cell_start->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

            rCells.push_back(p_cell_start); // Add the cell

            // Add the endpoint to the vector of nodes and cells
            end_index = node_index;

            Node<2>* p_node_end(new Node<2>(end_index, is_boundary, xEnd, yEnd)); // Define the node
            p_node_end->SetRadius(collagenRadius);

            rMesh.AddNode(p_node_end); // Add the node
            node_index += 1;

            // Set placeholder cell cycle model (i.e. no proliferation)
            NoCellCycleModel* p_cycle_model_end = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
            p_cycle_model_end->SetDimension(2);

            CellPtr p_cell_end(new Cell(p_wildtype_state, p_cycle_model_end));
            p_cell_end->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
            p_cell_end->GetCellData()->SetItem("volume", M_PI * pow(collagenRadius, 2.0)); // Mature cell volume
            p_cell_end->GetCellData()->SetItem("attachment", -1.0); // BM attachment
            p_cell_end->GetCellData()->SetItem("direction", theta);
            p_cell_end->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

            rCells.push_back(p_cell_end); // Add the cell

            // We will now store this index pair for future use when we mark cross-links
            std::pair<unsigned, unsigned> node_pair = std::make_pair(start_index, end_index); 
            rMarkedPairs[node_pair] = {start_index, end_index}; 

            c_vector<double, 3> fibre_and_lower_y_coord;
            fibre_and_lower_y_coord[0] = (double)start_index;
            fibre_and_lower_y_coord[1] = (double)end_index;
            fibre_and_lower_y_coord[2] = std::min(yStart, yEnd);

            rPairsAndLowerEndpoint.push_back(fibre_and_lower_y_coord);
        }

        rMesh.ReMesh(); // Re-mesh, as we've added several nodes

        // Sort the pairs by the lower y-coordinate
        std::sort(rPairsAndLowerEndpoint.begin(), rPairsAndLowerEndpoint.end(), SortByCoordinate);

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

    // Function to determine all of the fibre intersections
    void DetermineFibreIntersections(std::vector<c_vector<unsigned, 5> >& rFibreIntersections,
                                        std::vector<c_vector<double, 3> >& rPairsAndLowerEndpoint,
                                        std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> >& rMarkedPairs,
                                        Cylindrical2dNodesOnlyMesh& rMesh,
                                        std::vector<CellPtr>& rCells,
                                        double collagenRadius)
    {

        // Initialise the node index counter to add crosslinks
        unsigned crosslink_node_index = rMesh.GetNumNodes();

        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_collagen_type(CellPropertyRegistry::Instance()->Get<CollagenCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        // Iterate through, based on a rough plane sweeping algorithm
        for (unsigned i = 0; i < rPairsAndLowerEndpoint.size(); i++)
        {
            // Determine the upper y-coordinate, which determines how many segments we may need to consider.
            c_vector<double, 3> node_pair_and_lower_endpoint = rPairsAndLowerEndpoint[i];
            unsigned index_A = (unsigned)node_pair_and_lower_endpoint[0]; // First index of pair
            unsigned index_B = (unsigned)node_pair_and_lower_endpoint[1]; // Second index of pair
            double lower_y_coord = node_pair_and_lower_endpoint[2];

            // Get the locations
            c_vector<double, 2> node_location_A = rMesh.GetNode(index_A)->rGetLocation();
            c_vector<double, 2> node_location_B = rMesh.GetNode(index_B)->rGetLocation();

            double upper_y_coord = std::max(node_location_A[1], node_location_B[1]); // The y-coordinate of the upper end point should tell us when to stop considering neighbours

            std::vector<c_vector<double, 3> > considered_pairs_and_x_coords; // We will sort by the x-coordinate of the intersection with the lower y-coordinate

            for (unsigned j = 0; j < rPairsAndLowerEndpoint.size(); j++) // We arguably shouldn't have to consider anything beyond the current segment
            {
                c_vector<double, 3> current_pair_and_lower_endpoint = rPairsAndLowerEndpoint[j];

                unsigned current_index_A = (unsigned)current_pair_and_lower_endpoint[0]; // First index of pair
                unsigned current_index_B = (unsigned)current_pair_and_lower_endpoint[1]; // Second index of pair
                double current_lower_endpoint = current_pair_and_lower_endpoint[2];

                if (current_lower_endpoint > upper_y_coord)
                {
                    break;
                }
                else
                {
                    // Get the locations
                    c_vector<double, 2> current_node_location_A = rMesh.GetNode(current_index_A)->rGetLocation();
                    c_vector<double, 2> current_node_location_B = rMesh.GetNode(current_index_B)->rGetLocation();

                    double min_y = std::min(current_node_location_A[1], current_node_location_B[1]);
                    double max_y = std::max(current_node_location_A[1], current_node_location_B[1]);

                    if ( (lower_y_coord >= min_y)&&(lower_y_coord <= max_y) ) // If the current lower endpoint falls in the interior of a segment, we consider it.
                    {
                        c_vector<double, 3> pair_and_x_coord; // Define new pair and x-coord
                        pair_and_x_coord[0] = current_index_A;
                        pair_and_x_coord[1] = current_index_B;

                        // We will sort by the x-coordinate of the intersection with the current lower endpoint
                        double x_intersection = current_node_location_A[0] + (lower_y_coord - current_node_location_A[1])*(current_node_location_B[0] - current_node_location_A[0])/(current_node_location_B[1] - current_node_location_A[1]);
                        pair_and_x_coord[2] = x_intersection;

                        considered_pairs_and_x_coords.push_back(pair_and_x_coord);
                    }
                }
            }

            // If non-empty, we can add the current index as well and then determine any intersections
            // on the sorted vector
            if (considered_pairs_and_x_coords.size() > 1) 
            {

                // Sort the vector
                std::sort(considered_pairs_and_x_coords.begin(), considered_pairs_and_x_coords.end(), SortByCoordinate);

                // Now iterate through the vectors and consider only consecutive fibres
                for (unsigned k = 0; k < (considered_pairs_and_x_coords.size() - 1); k++)
                {
                    // Get the start and endpoints of the first fibre
                    c_vector<double, 3> fibre_A = considered_pairs_and_x_coords[k];

                    unsigned start_index_A = (unsigned)fibre_A[0];
                    unsigned end_index_A = (unsigned)fibre_A[1];

                    c_vector<double, 2> fibre_start_A = rMesh.GetNode(start_index_A)->rGetLocation();
                    c_vector<double, 2> fibre_end_A = rMesh.GetNode(end_index_A)->rGetLocation();

                    // Get the start and endpoints of the second fibre
                    c_vector<double, 3> fibre_B = considered_pairs_and_x_coords[k + 1];

                    unsigned start_index_B = (unsigned)fibre_B[0];
                    unsigned end_index_B = (unsigned)fibre_B[1];

                    c_vector<double, 2> fibre_start_B = rMesh.GetNode(start_index_B)->rGetLocation();
                    c_vector<double, 2> fibre_end_B = rMesh.GetNode(end_index_B)->rGetLocation();
                    
                    // Check if they intersect
                    if (DoFibresIntersect(fibre_start_A, fibre_end_A, fibre_start_B, fibre_end_B))
                    {
                        // We've found an intersection!
                        c_vector<unsigned, 5> fibre_intersection;

                        // We always order fibre intersections by the start index, so that it makes it easier to determine
                        // whether or not we've already encountered this intersection
                        if (start_index_A < start_index_B)
                        {
                            fibre_intersection[0] = start_index_A;
                            fibre_intersection[1] = end_index_A;
                            fibre_intersection[2] = start_index_B;
                            fibre_intersection[3] = end_index_B;
                        }
                        else
                        {
                            fibre_intersection[0] = start_index_B;
                            fibre_intersection[1] = end_index_B;
                            fibre_intersection[2] = start_index_A;
                            fibre_intersection[3] = end_index_A;
                        }

                        fibre_intersection[4] = crosslink_node_index;

                        // We only add the intersection if it's new.
                        unsigned count = 0;
                        for (unsigned l = 0; l < rFibreIntersections.size(); l++)
                        {
                            c_vector<unsigned, 5> current_intersection = rFibreIntersections[l];

                            if ( AreFibreIntersectionsEqual(current_intersection, fibre_intersection) )
                            {
                                break;
                            }
                            else
                            {
                                count += 1;
                            }
                        }
                        if (count == rFibreIntersections.size()) // If we're at the end, then this intersection hasn't been accounted for.
                        {
                            rFibreIntersections.push_back(fibre_intersection);

                            // Add the node and cell
                            unsigned start_A = fibre_intersection[0];
                            unsigned end_A = fibre_intersection[1];
                            unsigned start_B = fibre_intersection[2];
                            unsigned end_B = fibre_intersection[3];

                            c_vector<double, 2> start_location_A = rMesh.GetNode(start_A)->rGetLocation();
                            c_vector<double, 2> end_location_A = rMesh.GetNode(end_A)->rGetLocation();
                            c_vector<double, 2> start_location_B = rMesh.GetNode(start_B)->rGetLocation();
                            c_vector<double, 2> end_location_B = rMesh.GetNode(end_B)->rGetLocation(); 

                            c_vector<double, 2> intersection_location = GetFibreIntersection(start_location_A, end_location_A, start_location_B, end_location_B);

                            // Add a new node to the mesh
                            Node<2>* p_node(new Node<2>(crosslink_node_index, false, intersection_location[0], intersection_location[1]));
                            unsigned new_crosslink_index = rMesh.AddNode(p_node);
                            
                            // Add a new cell
                            NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                            p_cycle_model->SetDimension(2);

                            double theta = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

                            CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                            p_cell->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                            p_cell->GetCellData()->SetItem("volume", M_PI * pow(collagenRadius, 2.0)); // Mature cell volume
                            p_cell->GetCellData()->SetItem("attachment", -1.0); // BM attachment
                            p_cell->GetCellData()->SetItem("direction", theta);
                            p_cell->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

                            rCells.push_back(p_cell);

                            // Also update the fibre and cross-link maps
                            std::pair<unsigned, unsigned> pair_A = std::make_pair(start_A, end_A);
                            std::pair<unsigned, unsigned> pair_B = std::make_pair(start_B, end_B);

                            std::vector<unsigned> marked_indices_A = rMarkedPairs[pair_A];
                            marked_indices_A.push_back(new_crosslink_index); // Add the crosslink node
                            rMarkedPairs[pair_A] = marked_indices_A; // Update the node pair map

                            std::vector<unsigned> marked_indices_B = rMarkedPairs[pair_B];
                            marked_indices_B.push_back(new_crosslink_index); // Add the crosslink node
                            rMarkedPairs[pair_B] = marked_indices_B; // Update the node pair map

                            crosslink_node_index += 1; // Update the crosslink node index
                        }
                    }
                    
                }
            }

        }
    }

    // Function to compare the coordinates of the fibres. The assumption is that we're comparing
    // x-coordinates to x-coordinates or y-coordinates to y-coordinates
    static bool SortByCoordinate(c_vector<double, 3>& fibreA, c_vector<double, 3>& fibreB)
    {
        bool cond = false;

        if (fibreA[2] != fibreB[2])
        {
            cond = (fibreA[2] < fibreB[2]);
        }
        else
        {
            cond = (fibreA[0] < fibreB[0]);
        }

        return cond;
    }

    static bool AreFibreIntersectionsEqual(c_vector<unsigned, 5>& fibreIntersectionA, c_vector<unsigned, 5>& fibreIntersectionB)
    {
        return ( (fibreIntersectionA[0]==fibreIntersectionB[0])&&(fibreIntersectionA[1]==fibreIntersectionB[1])
                    &&(fibreIntersectionA[2]==fibreIntersectionB[2])&&(fibreIntersectionA[3]==fibreIntersectionB[3]) );
    }

    // Function to determine if two fibres intersect, based on a parametrisation from the
    // respective start and endpoints. If both of the resulting parametrisations lie between
    // 0 and 1, they intersect.
    bool DoFibresIntersect(c_vector<double, 2> fibreStartA, c_vector<double, 2> fibreEndA, 
                            c_vector<double, 2> fibreStartB, c_vector<double, 2> fibreEndB)
    {
        // Define the respective x and y coordinates
        double xStartA = fibreStartA[0];
        double xEndA = fibreEndA[0];
        double yStartA = fibreStartA[1];
        double yEndA = fibreEndA[1];

        double xStartB = fibreStartB[0];
        double xEndB = fibreEndB[0];
        double yStartB = fibreStartB[1];
        double yEndB = fibreEndB[1];

        double t_A = ( xEndB*(yStartA - yStartB) + xStartA*(yStartB - yEndB) + xStartB*(yEndB - yStartA) )
                        /( -(xEndB - xStartB)*(yEndA - yStartA) + (xEndA - xStartA)*(yEndB - yStartB) );
        double t_B = ( xStartB*(yStartA - yEndA) + xStartA*(yEndA - yStartB) + xEndA*(yStartB - yStartA) )
                        /( (xEndB - xStartB)*(yEndA - yStartA) - (xEndA - xStartA)*(yEndB - yStartB) );

        return ( (t_A > 0.0)&&(t_A < 1.0)&&(t_B > 0.0)&&(t_B < 1.0) );
    }

    // Function that returns the fibre intersection
    c_vector<double, 2> GetFibreIntersection(c_vector<double, 2> fibreStartA, c_vector<double, 2> fibreEndA, 
                            c_vector<double, 2> fibreStartB, c_vector<double, 2> fibreEndB)
    {
        // Define the respective x and y coordinates
        double xStartA = fibreStartA[0];
        double xEndA = fibreEndA[0];
        double yStartA = fibreStartA[1];
        double yEndA = fibreEndA[1];

        double xStartB = fibreStartB[0];
        double xEndB = fibreEndB[0];
        double yStartB = fibreStartB[1];
        double yEndB = fibreEndB[1];

        double t_A = ( xEndB*(yStartA - yStartB) + xStartA*(yStartB - yEndB) + xStartB*(yEndB - yStartA) )
                        /( -(xEndB - xStartB)*(yEndA - yStartA) + (xEndA - xStartA)*(yEndB - yStartB) );

        c_vector<double, 2> intersection;
        intersection[0] = xStartA + t_A * (xEndA - xStartA);
        intersection[1] = yStartA + t_A * (yEndA - yStartA);

        return intersection;
    }

public:
    void TestEpidermis()
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

        unsigned num_reticular_fibres = 24 * reticular_fibroblasts_across * reticular_fibroblasts_up;
        unsigned num_papillary_fibres = 20 * papillary_fibroblasts_across * papillary_fibroblasts_up;

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Mechanical parameters
        double spring_stiffness = 15.0; // Spring stiffness
        double bm_stiffness = 10.0; // Basement membrane attachment strength
        double target_curvature = 0.0; // Target curvature
        double reorientation_strength = 2.5;
        bool apply_cell_to_ecm_force = false;
        bool apply_force_on_marked_springs = false;

        // Set the stiffness multiplication factors
        std::map<unsigned, double> stiffness_multiplication_factors;

        stiffness_multiplication_factors[0] = 1.0; // Stem cells 
        stiffness_multiplication_factors[1] = 1.0; // TA cells 
        stiffness_multiplication_factors[2] = 1.0; // Differentiated cells 
        stiffness_multiplication_factors[3] = 1.0; // Fibroblasts cells 
        stiffness_multiplication_factors[5] = 0.0; // Collagen cells 

        // Shape scale of fibroblasts and collagen fibres
        double fibroblast_scale = 0.5;
        double collagen_radius = 0.5;
        double collagen_fibre_length = 1.0;

        std::vector<CellPtr> cells; // Vector of cells
        std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> > fibres_to_mark; // How we will track cross-links
        std::vector<c_vector<double, 3> > fibres_and_lower_y_coord; // Sort the fibres by the y-coordinate of the lower endpoint, so that we can determine intersections

        HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

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

        GenerateNodesAndCells(*p_generating_mesh, cells, fibres_to_mark, fibres_and_lower_y_coord,
                                cells_across, cells_up, 
                                reticular_dermis_height, papillary_dermis_height,
                                collagen_radius, collagen_fibre_length,
                                num_reticular_fibres, num_papillary_fibres,
                                reticular_fibroblasts_across, reticular_fibroblasts_up, 
                                papillary_fibroblasts_across, papillary_fibroblasts_across, 
                                fibroblast_scale); // Generate the initial tissue

        Cylindrical2dNodesOnlyMesh* p_cylindrical_mesh = new Cylindrical2dNodesOnlyMesh(1.0*cells_across);
		p_cylindrical_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0); //Construct mesh

        // Let's now determine all the fibre intersections (and their crosslinks)
        std::vector<c_vector<unsigned, 5> > intersecting_fibres;

        DetermineFibreIntersections(intersecting_fibres, fibres_and_lower_y_coord, 
                                        fibres_to_mark, *p_cylindrical_mesh, cells, collagen_radius);

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_cylindrical_mesh, cells); // Used for periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        // Mark the ECM fibres
        // Mark all the relevant cell pairs
        for (std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> >::iterator map_iter = fibres_to_mark.begin();
        map_iter != fibres_to_mark.end();
        map_iter++)
        {
            std::vector<unsigned> nodes_along_fibre = map_iter->second; // These are all the possible nodes along the fibre

            // We will sort the node indices now by x-coordinate and then mark consecutive pairs
            std::vector<std::pair<double, unsigned> > x_coordinates_and_nodes;

            // Determine the x-coordinate of each node
            for (unsigned j = 0; j < nodes_along_fibre.size(); j++)
            {
                unsigned node_index = nodes_along_fibre[j];
                double x = cell_population.rGetMesh().GetNode(node_index)->rGetLocation()[0];

                std::pair<double, unsigned> x_and_index = std::make_pair(x, node_index);

                x_coordinates_and_nodes.push_back(x_and_index);
            }
            
            std::sort(x_coordinates_and_nodes.begin(), x_coordinates_and_nodes.end()); // Sort the vector by x-coordinate

            // Let's average out the new radii, otherwise putting in these new crosslinks will push everything out.
            double new_node_radius = collagen_radius;

            for (unsigned j = 0; j < (x_coordinates_and_nodes.size() - 1); j++)
            {
                unsigned start_index = x_coordinates_and_nodes[j].second;
                unsigned end_index = x_coordinates_and_nodes[j + 1].second;

                // Update the radii
                cell_population.rGetMesh().GetNode(start_index)->SetRadius(new_node_radius);
                cell_population.rGetMesh().GetNode(end_index)->SetRadius(new_node_radius);

                // Create the cell pair
                CellPtr p_cell_start = cell_population.GetCellUsingLocationIndex(start_index);
                CellPtr p_cell_end = cell_population.GetCellUsingLocationIndex(end_index);

                std::pair<CellPtr, CellPtr> cell_pair = cell_population.CreateCellPair(p_cell_start, p_cell_end);

                cell_population.MarkSpring(cell_pair);
            }

        }

        // Add cell writers
        cell_population.AddCellWriter<CellMigrationDirectionWriter>(); // Visualise cell migration directions
 
        OffLatticeSimulation<2> simulator(cell_population);

        // Set output directory
        std::string output_directory = M_OUTPUT_DIRECTORY + "/SteadyState";
        simulator.SetOutputDirectory(output_directory);

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

        // // Add basement membrane force
        // MAKE_PTR(DistanceBasedEpidermalBasementMembraneForce, p_bm_force);
        MAKE_PTR(EpidermalBasementMembraneForce, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        p_bm_force->SetTargetCurvature(target_curvature);
        simulator.AddForce(p_bm_force);

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

        // // // // // Create a modifier to track which cells are attached to the basement membrane.
        MAKE_PTR(BasementMembraneAttachmentTrackingModifier<2>, p_bm_attachment_tracking_modifier);
        p_bm_attachment_tracking_modifier->SetNeighbourhoodRadius(radius_of_interaction);
		simulator.AddSimulationModifier(p_bm_attachment_tracking_modifier);

        double max_height = 0.5 * sqrt(3.0) * (double)cells_up;

        // Add a cell killer to remove differentiated cells that are too far from the basement membrane.
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_sloughing_killer, (&cell_population, max_height + 1.0));
        simulator.AddCellKiller(p_sloughing_killer);

        // Modifier to track the collagen fibres
        MAKE_PTR(CollagenFibreTrackingModifier<2>, p_collagen_fibre_tracking_modifier);
        simulator.AddSimulationModifier(p_collagen_fibre_tracking_modifier);

        // Modifier to align the cell polarities with their velocities
        MAKE_PTR(PolarityTrackingModifier<2>, p_polarity_tracking_modifier);
        p_polarity_tracking_modifier->SetReorientationStrength(reorientation_strength);
        simulator.AddSimulationModifier(p_polarity_tracking_modifier);

        simulator.Solve(); // Run the simulation.

        // // Save simulation in steady state
		// CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

    }

    // void TestInjury()
    // {

    //     // Load the saved simulation
    //     OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(M_OUTPUT_DIRECTORY + "/SteadyState", M_END_TIME);

    //     double max_width = 0.0;
    //     double max_height = 0.0;
    //     // Kill the cells within the wound area
    //     for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
    //             cell_iter != p_simulator->rGetCellPopulation().End();
    //             ++cell_iter)
    //     {
    //         //Get location of cell
    //         double x = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[0];
    //         double y = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

    //         //If the cell is within the 'wound area', we kill it.
    //         if (x > max_width)
    //         {
    //             max_width = x;
    //         }
    //         if (y > max_width)
    //         {
    //             max_height = y;
    //         }
        
    //     }

    //     // Now we can define the wound geometry
    //     double wound_centre = 0.5*max_width;
    //     double wound_width = 0.25*max_width;
    //     double wound_base_height = 0.4*max_height;

    //     // Kill the cells within the wound area
    //     for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
    //             cell_iter != p_simulator->rGetCellPopulation().End();
    //             ++cell_iter)
    //     {
    //         //Get location of cell
    //         double x = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[0];
    //         double y = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

    //         //If the cell is within the 'wound area', we kill it.
    //         if ( (x > (wound_centre - 0.5*wound_width))&&(x < (wound_centre + 0.5*wound_width))&&(y > wound_base_height) )
    //         {
    //             cell_iter->Kill();
    //         }
    //         if (x < wound_centre)
    //         {
    //             double direction = 0.0;
    //             cell_iter->GetCellData()->SetItem("direction", direction);
    //             cell_iter->GetCellData()->SetItem("collagen", 0.0);
    //             cell_iter->GetCellData()->SetItem("orientation", 0.0);

    //         }
    //         else
    //         {
    //             double direction = M_PI;
    //             cell_iter->GetCellData()->SetItem("direction", direction);
    //             cell_iter->GetCellData()->SetItem("collagen", 0.0);
    //             cell_iter->GetCellData()->SetItem("orientation", 0.0);
    //         }
    //     }

    //     // Add cell writers
    //     p_simulator->rGetCellPopulation().AddCellWriter<CellMigrationDirectionWriter>();

    //     // // Add fibre-alignment-based migration force
    //     // MAKE_PTR(FibreAlignmentBasedMigrationForce<2>, p_migration_force);
    //     // p_migration_force->SetMigrationForceStrength(0.0*M_DT);
    //     // p_migration_force->SetReorientationStrength(0.0*M_DT);
    //     // p_simulator->AddForce(p_migration_force);

    //     //Set output directory
    //     std::string output_directory = M_OUTPUT_DIRECTORY + "/Epidermis/Injury";
    //     p_simulator->SetOutputDirectory(output_directory);

    //     p_simulator->SetSamplingTimestepMultiple(M_INJURY_SAMPLING_TIMESTEP); //Sample the simulation at every hour
    //     p_simulator->SetEndTime(2*M_END_TIME); //Hopefully this is long enough for a steady state

    //     p_simulator->Solve(); // Run the simulation.

    //     // Save simulation in after injury too
	// 	CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator);

    // }
};

#endif /* TESTEPIDERMISINCROSSSECTIONALGEOMETRY */