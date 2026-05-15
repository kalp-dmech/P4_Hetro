

#include "Phase_Field.h"


void PhaseField::solve_linear_system_elastic()
{
    // std::cout << "Solving for u is started" << std::endl;

    TimerOutput::Scope ts (computing_timer, "solve_linear_system_elastic");

    SolverControl solver_control (10000,
        1e-12* system_rhs_elastic.l2_norm());

    SolverCG<LA::MPI::Vector> solver (solver_control);

    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#else
    // Trilinos defaults are good
#endif
    LA::MPI::PreconditionAMG preconditioner;
    preconditioner.initialize (system_matrix_elastic, data);

    solver.solve (system_matrix_elastic,
        completely_distributed_solution_elastic, system_rhs_elastic,
        preconditioner);

    last_elastic_cg_iterations = solver_control.last_step();
    last_elastic_residual = solver_control.last_value();

    constraints_elastic.distribute (completely_distributed_solution_elastic);

    locally_relevant_solution_elastic = completely_distributed_solution_elastic;
    //locally_relevant_solution_elastic.update_ghost_values();
    
 /*  SparseDirectUMFPACK A_direct;
    A_direct.initialize(u_system_matrix);
    A_direct.vmult(u_solution, u_system_rhs);

    std::cout << "Solving for u is finished" << std::endl;*/
    
}
