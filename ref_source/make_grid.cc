

#include "Phase_Field.h"

void PhaseField::make_grid()
{

    GridIn<2> gridin;
    gridin.attach_triangulation(triangulation);
    const std::string mesh_file = "mesh/geometry.msh";
    std::ifstream f(mesh_file);
    AssertThrow(f.is_open(),
                ExcMessage("Could not open mesh file: " + mesh_file));
    gridin.read_msh(f);
    //triangulation.refine_global(1);

    pcout << "No. of levels in triangulation: " << triangulation.n_global_levels () << std::endl;

    dof_handler_damage.distribute_dofs (fe_damage);
    dof_handler_elastic.distribute_dofs (fe_elastic);

    pcout << "   Number of locally owned cells on the process:       "
   << triangulation.n_locally_owned_active_cells () << std::endl;

    pcout << "Number of global cells:" << triangulation.n_global_active_cells ()
    << std::endl;

    pcout << "  Total Number of globally active cells:       "
    << triangulation.n_global_active_cells () << std::endl
    << "   Number of degrees of freedom for elasticity: "
    << dof_handler_elastic.n_dofs () << std::endl
    << " Number of degrees of freedom for damage: "
    << dof_handler_damage.n_dofs () << std::endl;

    //Initialising damage vectors
    locally_owned_dofs_damage = dof_handler_damage.locally_owned_dofs ();
    locally_relevant_dofs_damage = DoFTools::extract_locally_relevant_dofs (
        dof_handler_damage);

    completely_distributed_solution_damage_old.reinit (
        locally_owned_dofs_damage, mpi_communicator);
    completely_distributed_solution_damage.reinit (
        locally_owned_dofs_damage, mpi_communicator);
    locally_relevant_solution_damage.reinit (locally_owned_dofs_damage,
        locally_relevant_dofs_damage, mpi_communicator);

    locally_owned_dofs_elastic = dof_handler_elastic.locally_owned_dofs ();
    locally_relevant_dofs_elastic = DoFTools::extract_locally_relevant_dofs (
        dof_handler_elastic);

    completely_distributed_solution_elastic_old.reinit (
        locally_owned_dofs_elastic, mpi_communicator);


}
