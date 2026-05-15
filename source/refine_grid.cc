

#include "Phase_Field.h"


void PhaseField::refine_grid(double boundary_value)
{
    FEValues<2> fe_values_damage (fe_damage, quadrature_formula_damage,
         update_values | update_gradients | update_JxW_values
         | update_quadrature_points);


    parallel::distributed::ContinuousQuadratureDataTransfer<2, MyQData> data_transfer (
       fe_damage, quadrature_formula_damage, quadrature_formula_damage);


    Vector<float> estimated_error_per_cell (
        triangulation.n_locally_owned_active_cells ());

    KellyErrorEstimator<2>::estimate (dof_handler_damage,
        QGauss<1> (fe_damage.degree + 1),
          { }, locally_relevant_solution_damage, estimated_error_per_cell);

    // Initialize SolutionTransfer object
    parallel::distributed::SolutionTransfer<2, LA::MPI::Vector> soltransDamage (dof_handler_damage);

    // Initialize SolutionTransfer object
    parallel::distributed::SolutionTransfer<2, LA::MPI::Vector> soltransElastic (dof_handler_elastic);

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction (
        triangulation, estimated_error_per_cell, 0.01, // top 1% cells marked for refinement
        0.0); // bottom 0 % cells marked for coarsening

    if (triangulation.n_global_levels () >= 4)
    {
        for (const auto &cell : triangulation.active_cell_iterators_on_level (3))
            if (cell->is_locally_owned ())
                cell->clear_refine_flag ();
    }

    // prepare the triangulation,
    triangulation.prepare_coarsening_and_refinement ();

    // prepare CellDataStorage object for refinement
    data_transfer.prepare_for_coarsening_and_refinement (triangulation,
        quadrature_point_history_field);

    // prepare the SolutionTransfer object for coarsening and refinement
    // and give the solution vector that we intend to interpolate later,
    soltransDamage.prepare_for_coarsening_and_refinement (
        locally_relevant_solution_damage);
    soltransElastic.prepare_for_coarsening_and_refinement (
        locally_relevant_solution_elastic);

    triangulation.execute_coarsening_and_refinement ();

    // redistribute dofs,
    dof_handler_damage.distribute_dofs (fe_damage);
    dof_handler_elastic.distribute_dofs (fe_elastic);

    // Recreate locally_owned_dofs and locally_relevant_dofs index sets
    locally_owned_dofs_damage = dof_handler_damage.locally_owned_dofs ();
    locally_relevant_dofs_damage = DoFTools::extract_locally_relevant_dofs (
        dof_handler_damage);

    completely_distributed_solution_damage_old.reinit (
        locally_owned_dofs_damage, mpi_communicator);
    soltransDamage.interpolate (completely_distributed_solution_damage_old);

    // Apply constraints on the interpolated solution to make sure it conforms with the new mesh
    setup_constraints_damage();
    setup_system_damage();

    constraints_damage.distribute (completely_distributed_solution_damage_old);

    // Copy completely_distributed_solution_damage_old to locally_relevant_solution_damage
    locally_relevant_solution_damage.reinit (locally_owned_dofs_damage,
        locally_relevant_dofs_damage, mpi_communicator);
    locally_relevant_solution_damage =
        completely_distributed_solution_damage_old;

    // Interpolating elastic solution similarly
    locally_owned_dofs_elastic = dof_handler_elastic.locally_owned_dofs ();
    locally_relevant_dofs_elastic = DoFTools::extract_locally_relevant_dofs (
        dof_handler_elastic);

    completely_distributed_solution_elastic_old.reinit (
        locally_owned_dofs_elastic, mpi_communicator);
    soltransElastic.interpolate (completely_distributed_solution_elastic_old);

    // Apply constraints on the interpolated solution to make sure it conforms with the new mesh

    setup_constraints_elastic (boundary_value);
    constraints_elastic.distribute (
        completely_distributed_solution_elastic_old);
    setup_system_elastic(boundary_value);

    // Copy completely_distributed_solution_damage_old to locally_relevant_solution_damage
    locally_relevant_solution_elastic.reinit (locally_owned_dofs_elastic,
        locally_relevant_dofs_elastic, mpi_communicator);
    locally_relevant_solution_elastic =
        completely_distributed_solution_elastic_old;

    for (const auto &cell : triangulation.active_cell_iterators ())
      {
        if (cell->is_locally_owned ())
          quadrature_point_history_field.initialize (cell, quadrature_formula_damage.size());
      }
    data_transfer.interpolate ();
}
