/**
* This class is the constructor for the Step6 class
*/

#include "Phase_Field.h"

namespace
{
    const double k_val = 1e-6;

    const unsigned int elastic_fe_degree_val = 1;
    const unsigned int damage_fe_degree_val = 1;

    const unsigned int output_interval_val = 10;
    const unsigned int amr_frequency_val = 1;
    const unsigned int amr_stop_step_val = PhaseFieldVariables::n_steps;
    const unsigned int amr_max_level_val = 5;
    const double amr_refine_length_val = 0.25;
}

PhaseField::PhaseField()   // Constructor for initializing the objects
        : PhaseField("output",
                     "force_displacement.txt",
                     false)
    {}

PhaseField::PhaseField(const std::string &output_directory_value,
                       const std::string &force_displacement_file_value,
                       const bool enable_amr_value)

        : lambda (PhaseFieldVariables::material_1_lambda),
          mu (PhaseFieldVariables::material_1_mu),
          Gc (PhaseFieldVariables::material_1_Gc),
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

         bottom_boundary_id(PhaseFieldVariables::bottom_boundary_id),
         top_boundary_id(PhaseFieldVariables::top_boundary_id),
         crack_boundary_id(PhaseFieldVariables::crack_boundary_id),
         applied_displacement(0.0),
         force_data(0.0),
         last_elastic_cg_iterations(0),
         last_damage_cg_iterations(0),
         last_elastic_residual(0.0),
         last_damage_residual(0.0),
         output_directory(output_directory_value),
         force_displacement_file(force_displacement_file_value),
         enable_amr(enable_amr_value),
         n_steps(PhaseFieldVariables::n_steps),
         max_staggered_iterations_1(PhaseFieldVariables::max_staggered_iterations_1),
         max_staggered_iterations_2(PhaseFieldVariables::max_staggered_iterations_2),
         output_interval(output_interval_val),
         amr_frequency(amr_frequency_val),
         amr_stop_step(amr_stop_step_val),
         amr_max_level(amr_max_level_val),
         inc_large(PhaseFieldVariables::inc_large),
         inc_small(PhaseFieldVariables::inc_small),
         tol_1(PhaseFieldVariables::tol_1),
         tol_2(PhaseFieldVariables::tol_2),
         amr_refine_length(amr_refine_length_val)

    {}
