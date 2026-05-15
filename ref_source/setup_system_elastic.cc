#include "Phase_Field.h"

void PhaseField::setup_system_elastic(double boundary_value)
{

    TimerOutput::Scope ts (computing_timer, "setup_system_elastic");

    //dof_handler_elastic.distribute_dofs(fe_elastic);
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

}
