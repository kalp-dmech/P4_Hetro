

#include "Phase_Field.h"


void PhaseField::assemble_system_damage()
{
    TimerOutput::Scope ts (computing_timer, "assembly_damage");
    
    //system_matrix_damage = 0;
    //system_rhs_damage = 0;
        

    FEValues<2> fe_values_damage (fe_damage, quadrature_formula_damage,
         update_values | update_gradients | update_JxW_values
         | update_quadrature_points);
    FEValues<2> fe_values_elastic (fe_elastic, quadrature_formula_elastic,
        update_values | update_gradients | update_JxW_values
        | update_quadrature_points);

    const unsigned int dofs_per_cell = fe_damage.n_dofs_per_cell ();
    const unsigned int n_q_points = quadrature_formula_damage.size ();

    FullMatrix<double> cell_matrix_damage (dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs_damage (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    // Storing strain tensor for all Gauss points of a cell in vector strain_values.
    std::vector<SymmetricTensor<2, 2>> strain_values (n_q_points);

    for (const auto &cell : dof_handler_damage.active_cell_iterators())
      if (cell->is_locally_owned ())
       {

          std::vector<std::shared_ptr<MyQData>> qpdH =
                quadrature_point_history_field.get_data (cell);

          const DoFHandler<2>::active_cell_iterator elastic_cell =
              Triangulation<2>::active_cell_iterator (cell)->as_dof_handler_iterator (
                  dof_handler_elastic);

          fe_values_damage.reinit (cell);
          fe_values_elastic.reinit (elastic_cell);

          const FEValuesExtractors::Vector displacements (
              /* first vector component = */0);
          fe_values_elastic[displacements].get_function_symmetric_gradients (
              locally_relevant_solution_elastic, strain_values);
        // Get strain values at quadrature points from u_new_solution and store in strain_values and pass these
        // as input to
        // the H_plus function
          cell_matrix_damage = 0;
          cell_rhs_damage = 0;

        for (const unsigned int q_index : fe_values_damage.quadrature_point_indices())
        {

            // Compute positive strain energy
            double H_call = H_plus(strain_values[q_index]);

            // Update history (irreversibility)
            const double H = std::max(H_call, qpdH[q_index]->value_H);
            //pcout << "The H value at qp " << q_index << " is " << H << std::endl;
            qpdH[q_index]->value_H_new = H;

            for (const unsigned int i : fe_values_damage.dof_indices())
            {
                auto N_v_i = fe_values_damage.shape_value(i, q_index);
                auto grad_N_v_i = fe_values_damage.shape_grad(i, q_index);

                cell_rhs_damage(i) += (2*lo/Gc*H*N_v_i)*fe_values_damage.JxW(q_index);

               for (const unsigned int j : fe_values_damage.dof_indices())
                {
                   auto N_phi_j = fe_values_damage.shape_value(j, q_index);
                   auto grad_N_phi_j = fe_values_damage.shape_grad(j, q_index);
                    

                   cell_matrix_damage(i, j) += ( (lo*lo*grad_N_v_i*grad_N_phi_j) + (1 + (2*H*lo)/Gc) * (N_v_i*N_phi_j) ) * fe_values_damage.JxW(q_index);
                }
            }
        }


          cell->get_dof_indices (local_dof_indices);
          constraints_damage.distribute_local_to_global (cell_matrix_damage,
              cell_rhs_damage, local_dof_indices, system_matrix_damage,
              system_rhs_damage);
    }

    system_matrix_damage.compress (VectorOperation::add);
    system_rhs_damage.compress (VectorOperation::add);


}

