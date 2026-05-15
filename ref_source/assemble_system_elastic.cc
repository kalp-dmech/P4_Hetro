

#include "Phase_Field.h"
   

void PhaseField::assemble_system_elastic()
{
    TimerOutput::Scope ts (computing_timer, "assembly_elastic");

    //  Reinitialize the system matrix and RHS for each boundary condition
    //system_matrix_elastic = 0;
    //system_rhs_elastic = 0;

    //QGauss<2> quadrature_formula(fe_elastic.degree + 1);
    FEValues<2> fe_values_elastic (fe_elastic, quadrature_formula_elastic,
       update_values | update_gradients | update_quadrature_points
       | update_JxW_values);
    FEValues<2> fe_values_damage (fe_damage, quadrature_formula_damage,
        update_values | update_gradients | update_quadrature_points | update_JxW_values);


    const unsigned int dofs_per_cell = fe_elastic.n_dofs_per_cell();

    const unsigned int n_q_points = quadrature_formula_damage.size ();

    FullMatrix<double> cell_matrix_elastic(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs_elastic(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double>       damage_values (n_q_points); // phi values at gauss points

     FEValuesExtractors::Scalar u_x(0);
     FEValuesExtractors::Scalar u_y(1);

     std::vector<Tensor<1, 2>> u_grad(n_q_points);
     std::vector<Tensor<1, 2>> v_grad(n_q_points);

     SymmetricTensor<2, 2> strains;

    // typename DoFHandler<2>::active_cell_iterator cell2 = phi_dof_handler.begin_active();

    for (const auto &cell : dof_handler_elastic.active_cell_iterators())
        if (cell->is_locally_owned ())
        {

        cell_matrix_elastic = 0;
        cell_rhs_elastic    = 0;

        // Find corresponding elastic cell (same triangulation, different DoFHandler)
            const DoFHandler<2>::active_cell_iterator damage_cell =
                  Triangulation<2>::active_cell_iterator (cell)->as_dof_handler_iterator (
                      dof_handler_damage);


        fe_values_damage.reinit(damage_cell);
        fe_values_elastic.reinit(cell);

        fe_values_damage.get_function_values(locally_relevant_solution_damage, damage_values);

        fe_values_elastic[u_x].get_function_gradients(locally_relevant_solution_elastic, u_grad);
        fe_values_elastic[u_y].get_function_gradients(locally_relevant_solution_elastic, v_grad);

        for (const unsigned int q_index : fe_values_elastic.quadrature_point_indices())
        {
           auto gp_phi = damage_values[q_index];

            //Imposing the phi constraint
	/*
            double exx = u_grad[q_index][0];
            double eyy = v_grad[q_index][1];
            double exy = 0.5 * (u_grad[q_index][1] + v_grad[q_index][0]);

            strains[0][0] = exx;
            strains[0][1] = exy;
            strains[1][0] = exy;
            strains[1][1] = eyy;

            double tr_strain = trace(strains);

            double Mac_tr_strain_plus  = (tr_strain > 0.0) ? tr_strain : 0.0;
            double Mac_tr_strain_minus = (tr_strain < 0.0) ? tr_strain : 0.0;

            std::array<double,2> Principal_Strains = eigenvalues(strains);

            double Mac_principal_strain_1 = (Principal_Strains[0] > 0.0) ? Principal_Strains[0] : 0.0;
            double Mac_principal_strain_2 = (Principal_Strains[1] > 0.0) ? Principal_Strains[1] : 0.0;

            double Min_principal_strain_1 = (Principal_Strains[0] < 0.0) ? Principal_Strains[0] : 0.0;
            double Min_principal_strain_2 = (Principal_Strains[1] < 0.0) ? Principal_Strains[1] : 0.0;

            double tr_sqr_Mac_principal = std::pow(Mac_principal_strain_1, 2.0)
                                        + std::pow(Mac_principal_strain_2, 2.0);

            double tr_sqr_Min_principal = std::pow(Min_principal_strain_1, 2.0)
                                        + std::pow(Min_principal_strain_2, 2.0);

            double strain_energy_density_plus  = 0.5 * lambda * std::pow(Mac_tr_strain_plus, 2.0)
                                         + mu * tr_sqr_Mac_principal;

            double strain_energy_density_minus = 0.5 * lambda * std::pow(Mac_tr_strain_minus, 2.0)
                                         + mu * tr_sqr_Min_principal;

            if (strain_energy_density_plus < strain_energy_density_minus)
            {

                gp_phi = 0;
            }
		*/

            for (const unsigned int i : fe_values_elastic.dof_indices())
            {

                
                // If it doesn't work out try the below code for getting gradients and also try using auto instead of double
                
                double N_wx_i_x = fe_values_elastic[u_x].gradient(i, q_index)[0];
                double N_wx_i_y = fe_values_elastic[u_x].gradient(i, q_index)[1];
                double N_wy_i_x = fe_values_elastic[u_y].gradient(i, q_index)[0];
                double N_wy_i_y = fe_values_elastic[u_y].gradient(i, q_index)[1];

                cell_rhs_elastic(i) += 0;


               for (const unsigned int j : fe_values_elastic.dof_indices())
                {

                    
                // If it doesn't work out try the below code for getting gradient and also use auto instead of double.
                    
                    double N_ux_j_x = fe_values_elastic[u_x].gradient(j, q_index)[0];
                    double N_ux_j_y = fe_values_elastic[u_x].gradient(j, q_index)[1];
                    double N_uy_j_x = fe_values_elastic[u_y].gradient(j, q_index)[0];
                    double N_uy_j_y = fe_values_elastic[u_y].gradient(j, q_index)[1];
                    

                   double  Aij = (lambda * N_wx_i_x * N_ux_j_x) + (mu * (N_wx_i_x*N_ux_j_x + N_wx_i_x*N_ux_j_x + N_wx_i_y*N_ux_j_y));
                   double  Bij = (lambda * N_wx_i_x * N_uy_j_y) + (mu * N_wx_i_y * N_uy_j_x);
                   double  Cij = (lambda * N_wy_i_y * N_ux_j_x) + ( mu * N_wy_i_x * N_ux_j_y );
                   double  Dij = (lambda * N_wy_i_y  * N_uy_j_y) + mu * ( N_wy_i_x*N_uy_j_x + N_wy_i_y*N_uy_j_y + N_wy_i_y*N_uy_j_y );

                   cell_matrix_elastic(i, j) += ((((1-gp_phi)*(1-gp_phi))+k)* ( Aij + Bij + Cij + Dij )) * fe_values_elastic.JxW(q_index);
                   //cell_matrix(i, j) += ( (1-(0))* ( Aij + Bij + Cij + Dij )) * fe_values.JxW(q_index);

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

