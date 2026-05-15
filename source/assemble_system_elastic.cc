

#include "Phase_Field.h"
#include <algorithm>
   

void PhaseField::assemble_system_elastic()
{
    TimerOutput::Scope ts (computing_timer, "assembly_elastic");

    system_matrix_elastic = 0;
    system_rhs_elastic = 0;

    FEValues<2> fe_values_elastic (fe_elastic, quadrature_formula_elastic,
       update_values | update_gradients | update_quadrature_points
       | update_JxW_values);
    FEValues<2> fe_values_damage (fe_damage, quadrature_formula_damage,
        update_values | update_gradients | update_quadrature_points | update_JxW_values);


    const unsigned int dofs_per_cell = fe_elastic.n_dofs_per_cell();

    const unsigned int n_q_points = quadrature_formula_damage.size ();
    AssertDimension(quadrature_formula_elastic.size(), n_q_points);

    FullMatrix<double> cell_matrix_elastic(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs_elastic(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double>       damage_values (n_q_points); // phi values at gauss points

    const FEValuesExtractors::Scalar u_x(0);
    const FEValuesExtractors::Scalar u_y(1);

    for (const auto &cell : dof_handler_elastic.active_cell_iterators())
        if (cell->is_locally_owned ())
        {

        cell_matrix_elastic = 0;
        cell_rhs_elastic    = 0;

        // Find corresponding damage cell (same triangulation, different DoFHandler)
            const DoFHandler<2>::active_cell_iterator damage_cell =
                  Triangulation<2>::active_cell_iterator (cell)->as_dof_handler_iterator (
                      dof_handler_damage);


        fe_values_damage.reinit(damage_cell);
        fe_values_elastic.reinit(cell);
        const MaterialProperties material = material_properties(cell->material_id());

        fe_values_damage.get_function_values(locally_relevant_solution_damage, damage_values);

        for (const unsigned int q_index : fe_values_elastic.quadrature_point_indices())
        {
            const double gp_phi = damage_values[q_index];
            const double g = (1.0 - gp_phi) * (1.0 - gp_phi) + k;

            for (const unsigned int i : fe_values_elastic.dof_indices())
            {
                double N_wx_i_x = fe_values_elastic[u_x].gradient(i, q_index)[0];
                double N_wx_i_y = fe_values_elastic[u_x].gradient(i, q_index)[1];
                double N_wy_i_x = fe_values_elastic[u_y].gradient(i, q_index)[0];
                double N_wy_i_y = fe_values_elastic[u_y].gradient(i, q_index)[1];

               for (const unsigned int j : fe_values_elastic.dof_indices())
                {
                    double N_ux_j_x = fe_values_elastic[u_x].gradient(j, q_index)[0];
                    double N_ux_j_y = fe_values_elastic[u_x].gradient(j, q_index)[1];
                    double N_uy_j_x = fe_values_elastic[u_y].gradient(j, q_index)[0];
                    double N_uy_j_y = fe_values_elastic[u_y].gradient(j, q_index)[1];

                    // Isotropic elasticity with correct (lambda + 2*mu) diagonal
                    double Aij = (material.lambda + 2.0 * material.mu) * N_wx_i_x * N_ux_j_x
                               + material.mu * N_wx_i_y * N_ux_j_y;
                    double Bij = material.lambda * N_wx_i_x * N_uy_j_y
                               + material.mu * N_wx_i_y * N_uy_j_x;
                    double Cij = material.lambda * N_wy_i_y * N_ux_j_x
                               + material.mu * N_wy_i_x * N_ux_j_y;
                    double Dij = (material.lambda + 2.0 * material.mu) * N_wy_i_y * N_uy_j_y
                               + material.mu * N_wy_i_x * N_uy_j_x;

                    // Uniform degradation on full stiffness (Miehe 2010)
                    cell_matrix_elastic(i, j) +=
                        g * (Aij + Bij + Cij + Dij) *
                        fe_values_elastic.JxW(q_index);
                }
            }
        }
        cell->get_dof_indices(local_dof_indices);

            constraints_elastic.distribute_local_to_global (cell_matrix_elastic,
                  cell_rhs_elastic, local_dof_indices, system_matrix_elastic,
                  system_rhs_elastic);
    }

    system_matrix_elastic.compress (VectorOperation::add);
    system_rhs_elastic.compress (VectorOperation::add);


}
