#include "Phase_Field.h"

void PhaseField::setup_system_damage()
{
    TimerOutput::Scope ts (computing_timer, "setup_system_damage");
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
}
