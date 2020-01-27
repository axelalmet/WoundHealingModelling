#ifndef TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY_HPP_
#define TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CellLabel.hpp"
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "MeasuringOSStresses";
static const double M_DT = 0.05;
static const double M_END_TIME = 1.0;
//static const double M_SECOND_END_TIME = 1.0;
static const double M_SAMPLING_TIMESTEP = M_END_TIME / M_DT;

class TestCrossSectionalWoundHealing : public AbstractCellBasedTestSuite
{
public:
    void TestWounding()
    {

        //Set all the spring stiffness variables
        double epithelial_epithelial_stiffness = 45.0;
        double epithelial_stromal_stiffness = 45.0;
        double stromal_stromal_stiffness = 45.0;

        //Set the number of cells across and down for the array
        unsigned cells_across = 20;
        unsigned cells_up = 6;

        //Set the basement membrane force parameters

        double epithelial_epithelial_resting_spring_length = 1.0;

        double radius_of_interaction = 1.5;
        double division_separation = 0.1;

        HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

        // Construct the mesh that we'll actually use.
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, radius_of_interaction); //Construct mesh

        //Get the real indices
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        //Create the vector of cells
        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

        //Create tissue of cells. Initially we set them all to be differentiated
        std::vector<CellPtr> cells; //Create vector of cells

        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            //Set stochastic duration based cell cycle
            NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
            p_cycle_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
            p_cell->InitialiseCellCycleModel(); // For paranoia really.

            p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated

            cells.push_back(p_cell);
        }

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        //Set the division separation
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        //Get the maximum width so we know where to apply the right BC.
        double max_width = 0.0;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];

            if (x > max_width)
            {
                max_width = x;
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

        //Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        // MAKE_PTR(LinearSpringForceWithVariableRestLength<2>, p_spring_force);
        // //			p_spring_force->SetCutOffLength(epithelial_epithelial_resting_spring_length);
        // p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness); //Default is 15
        // p_spring_force->SetEpithelialStromalSpringStiffness(epithelial_stromal_stiffness); //Default is 15
        // p_spring_force->SetStromalStromalSpringStiffness(stromal_stromal_stiffness); //Default is 15
        // p_spring_force->SetEpithelialEpithelialRestingSpringLength(epithelial_epithelial_resting_spring_length);
        // p_spring_force->SetCutOffLength(radius_of_interaction);
        // p_spring_force->SetMeinekeDivisionRestingSpringLength(division_separation);
        // simulator.AddForce(p_spring_force);

        // c_vector<double, 2> point, normal;

        // point(0) = 0.0;
        // point(1) = 0.25;
        // normal(0) = 0.0;
        // normal(1) = -1.0;
        // MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        // simulator.AddCellPopulationBoundaryCondition(p_bc1);

        simulator.Solve();
    }
};

#endif /* TESTWOUNDHEALINGINCROSSSECTIONALGEOMETRY