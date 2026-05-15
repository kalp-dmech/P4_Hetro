

#include "Phase_Field.h"


void PhaseField::solve_linear_system_damage()
{
    TimerOutput::Scope ts (computing_timer, "solve_linear_system_damage");

    SolverControl solver_control (10000,
        1e-12* system_rhs_damage.l2_norm());

    SolverCG<LA::MPI::Vector> solver (solver_control);

    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#else
    // Trilinos defaults are good
#endif
    LA::MPI::PreconditionAMG preconditioner;
    preconditioner.initialize (system_matrix_damage, data);

    solver.solve (system_matrix_damage, completely_distributed_solution_damage,
        system_rhs_damage, preconditioner);

    constraints_damage.distribute (completely_distributed_solution_damage);
    locally_relevant_solution_damage = completely_distributed_solution_damage;
    //locally_relevant_solution_damage.update_ghost_values();
}
