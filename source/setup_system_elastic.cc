#include "Phase_Field.h"

void PhaseField::setup_system_elastic(double boundary_value)
{

    TimerOutput::Scope ts (computing_timer, "setup_system_elastic");

    LA::MPI::Vector saved_solution;
    LA::MPI::Vector saved_old_solution;
    const bool keep_current_solution =
        completely_distributed_solution_elastic.size() == dof_handler_elastic.n_dofs();
    const bool keep_old_solution =
        completely_distributed_solution_elastic_old.size() == dof_handler_elastic.n_dofs();

    if (keep_current_solution)
    {
        saved_solution.reinit(locally_owned_dofs_elastic, mpi_communicator);
        saved_solution = completely_distributed_solution_elastic;
    }
    if (keep_old_solution)
    {
        saved_old_solution.reinit(locally_owned_dofs_elastic, mpi_communicator);
        saved_old_solution = completely_distributed_solution_elastic_old;
    }

    locally_owned_dofs_elastic = dof_handler_elastic.locally_owned_dofs ();
    locally_relevant_dofs_elastic = DoFTools::extract_locally_relevant_dofs (
        dof_handler_elastic);

    locally_relevant_solution_elastic.reinit (locally_owned_dofs_elastic,
        locally_relevant_dofs_elastic, mpi_communicator);

    system_rhs_elastic.reinit (locally_owned_dofs_elastic, mpi_communicator);

    completely_distributed_solution_elastic.reinit (locally_owned_dofs_elastic,
        mpi_communicator);

    setup_constraints_elastic (boundary_value);

    DynamicSparsityPattern dsp (locally_relevant_dofs_elastic);

    DoFTools::make_sparsity_pattern (dof_handler_elastic, dsp,
        constraints_elastic, false);

    SparsityTools::distribute_sparsity_pattern (dsp,
        dof_handler_elastic.locally_owned_dofs (), mpi_communicator,
        locally_relevant_dofs_elastic);

    system_matrix_elastic.reinit (locally_owned_dofs_elastic,
        locally_owned_dofs_elastic, dsp, mpi_communicator);

    if (keep_current_solution)
    {
        completely_distributed_solution_elastic = saved_solution;
        locally_relevant_solution_elastic = saved_solution;
        locally_relevant_solution_elastic.update_ghost_values();
    }
    if (keep_old_solution)
        completely_distributed_solution_elastic_old = saved_old_solution;

}
