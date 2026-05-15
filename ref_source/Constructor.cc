/**
* This class is the constructor for the Step6 class
*/

#include "Phase_Field.h"

namespace
{
    // Hardcoded structural and control constants moved from variable_constructor.cc
    const double lambda_val = 121.1538; // kN/mm^2
    const double mu_val = 80.7692;
    const double k_val = 1e-6;

    const unsigned int elastic_fe_degree_val = 1;
    const unsigned int damage_fe_degree_val = 1;
    const unsigned int top_boundary_id_val = 2;

    const unsigned int output_interval_val = 30;
    const unsigned int amr_frequency_val = 1;
    const unsigned int amr_stop_step_val = PhaseFieldVariables::n_steps;
    const unsigned int amr_max_level_val = 5;
    const double inc_large_val = 1.0e-3;
    const double inc_small_val = 1.0e-5;
    const double amr_refine_length_val = 0.25;
}

PhaseField::PhaseField()   // Constructor for initializing the objects
        : PhaseField(PhaseFieldVariables::pore_counts.front(),
                     1,
                     "output",
                     "force_displacement.csv")
    {}

PhaseField::PhaseField(const unsigned int pore_count_value,
                       const unsigned int current_trial_number,
                       const std::string &output_directory_value,
                       const std::string &force_displacement_file_value)

        : lambda (lambda_val),
          mu (mu_val),
          Gc (PhaseFieldVariables::Gc),
          lo (PhaseFieldVariables::lo),
          k(k_val),
          mpi_communicator (MPI_COMM_WORLD),
          pcout (std::cout,
      (Utilities::MPI::this_mpi_process (mpi_communicator) == 0)),
          computing_timer (mpi_communicator, pcout, TimerOutput::never,
      TimerOutput::wall_times),

         triangulation (mpi_communicator),

         fe_elastic(FE_Q<2>(elastic_fe_degree_val),2),
         dof_handler_elastic(triangulation),
         quadrature_formula_elastic(fe_elastic.degree + 1),

         fe_damage(FE_Q<2>(damage_fe_degree_val),1),
         dof_handler_damage(triangulation),
         quadrature_formula_damage (fe_damage.degree + 1),

         top_boundary_id(top_boundary_id_val),
         pore_count(pore_count_value),
         pore_trial_number(current_trial_number),
         total_trial_count(PhaseFieldVariables::pore_trial_number),
         case_number(PhaseFieldVariables::case_number),
         random_seed_base(PhaseFieldVariables::random_seed_base),
         min_pore_size(PhaseFieldVariables::min_pore_size),
         max_pore_size(PhaseFieldVariables::max_pore_size),
         crack_spread_radius(PhaseFieldVariables::crack_spread_radius),
         boundary_margin(PhaseFieldVariables::boundary_margin),
         output_directory(output_directory_value),
         force_displacement_file(force_displacement_file_value),
         n_steps(PhaseFieldVariables::n_steps),
         max_staggered_iterations_1(PhaseFieldVariables::max_staggered_iterations_1),
         max_staggered_iterations_2(PhaseFieldVariables::max_staggered_iterations_2),
         output_interval(output_interval_val),
         amr_frequency(amr_frequency_val),
         amr_stop_step(amr_stop_step_val),
         amr_max_level(amr_max_level_val),
         inc_large(inc_large_val),
         inc_small(inc_small_val),
         tol_1(PhaseFieldVariables::tol_1),
         tol_2(PhaseFieldVariables::tol_2),
         amr_refine_length(amr_refine_length_val)

    {}
