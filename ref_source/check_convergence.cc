

#include "Phase_Field.h"


bool
  PhaseField::check_convergence ()
{
    double tol = 1e-3;
    LA::MPI::Vector solution_damage_difference (locally_owned_dofs_damage,
        mpi_communicator);
    LA::MPI::Vector solution_elastic_difference (locally_owned_dofs_elastic,
        mpi_communicator);
    LA::MPI::Vector solution_damage_difference_ghost (locally_owned_dofs_damage,
        locally_relevant_dofs_damage, mpi_communicator);
    LA::MPI::Vector solution_elastic_difference_ghost (
        locally_owned_dofs_elastic, locally_relevant_dofs_elastic,
        mpi_communicator);

    solution_damage_difference = locally_relevant_solution_damage;

    solution_damage_difference -= completely_distributed_solution_damage_old;

    solution_elastic_difference = locally_relevant_solution_elastic;

    solution_elastic_difference -= completely_distributed_solution_elastic_old;

    solution_damage_difference_ghost = solution_damage_difference;

    solution_elastic_difference_ghost = solution_elastic_difference;

    double error_elastic_solution_numerator, error_elastic_solution_denominator,
        error_damage_solution_numerator, error_damage_solution_denominator;

    error_damage_solution_numerator = solution_damage_difference.l2_norm ();
    error_elastic_solution_numerator = solution_elastic_difference.l2_norm ();
    error_damage_solution_denominator =
        completely_distributed_solution_damage.l2_norm ();
    error_elastic_solution_denominator =
        completely_distributed_solution_elastic.l2_norm ();

    double error_elastic_solution, error_damage_solution;
    error_damage_solution = error_damage_solution_numerator
        / error_damage_solution_denominator;
    error_elastic_solution = error_elastic_solution_numerator
        / error_elastic_solution_denominator;

    if ((error_elastic_solution < tol) && (error_damage_solution < tol))
        return true;
    else
        return false;
}
