#include "Phase_Field.h"

void PhaseField::setup_system_damage()
{
    TimerOutput::Scope ts (computing_timer, "setup_system_damage");
    LA::MPI::Vector saved_solution;
    LA::MPI::Vector saved_old_solution;
    const bool keep_current_solution =
        completely_distributed_solution_damage.size() == dof_handler_damage.n_dofs();
    const bool keep_old_solution =
        completely_distributed_solution_damage_old.size() == dof_handler_damage.n_dofs();

    if (keep_current_solution)
    {
        saved_solution.reinit(locally_owned_dofs_damage, mpi_communicator);
        saved_solution = completely_distributed_solution_damage;
    }
    if (keep_old_solution)
    {
        saved_old_solution.reinit(locally_owned_dofs_damage, mpi_communicator);
        saved_old_solution = completely_distributed_solution_damage_old;
    }

    dof_handler_damage.distribute_dofs(fe_damage);
    locally_owned_dofs_damage = dof_handler_damage.locally_owned_dofs ();
    locally_relevant_dofs_damage = DoFTools::extract_locally_relevant_dofs (
        dof_handler_damage);

    locally_relevant_solution_damage.reinit (locally_owned_dofs_damage,
        locally_relevant_dofs_damage, mpi_communicator);

    system_rhs_damage.reinit (locally_owned_dofs_damage, mpi_communicator);

    completely_distributed_solution_damage.reinit (locally_owned_dofs_damage,
        mpi_communicator);

    DynamicSparsityPattern dsp (locally_relevant_dofs_damage);

    DoFTools::make_sparsity_pattern (dof_handler_damage, dsp,
        constraints_damage, false);
    SparsityTools::distribute_sparsity_pattern (dsp,
        dof_handler_damage.locally_owned_dofs (), mpi_communicator,
        locally_relevant_dofs_damage);

    system_matrix_damage.reinit (locally_owned_dofs_damage,
        locally_owned_dofs_damage, dsp, mpi_communicator);

    if (keep_current_solution)
    {
        completely_distributed_solution_damage = saved_solution;
        locally_relevant_solution_damage = saved_solution;
        locally_relevant_solution_damage.update_ghost_values();
    }
    if (keep_old_solution)
        completely_distributed_solution_damage_old = saved_old_solution;
}
