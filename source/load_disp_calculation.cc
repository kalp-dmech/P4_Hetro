//
// Created by narasimhan-swaminathan on 4/1/25.
// Modified for parallel computation
//

#include "Phase_Field.h"

void PhaseField::load_disp_calculation(double boundary_value)
{
    QGauss<1> face_quadrature_formula(fe_elastic.degree + 1);
    FEFaceValues<2> fe_face_values(fe_elastic,
                                   face_quadrature_formula,
                                   update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values | update_normal_vectors);

    FEFaceValues<2> fe_face_values_damage(fe_damage,
                                   face_quadrature_formula,
                                   update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values | update_normal_vectors);

    const unsigned int n_q_f_points = face_quadrature_formula.size();
    const FEValuesExtractors::Scalar u_x(0);
    const FEValuesExtractors::Scalar u_y(1);

    std::vector<double> u_val(n_q_f_points);
    std::vector<double> v_val(n_q_f_points);
    std::vector<double> phi_val(n_q_f_points);

    std::vector<Tensor<1,2>> u_grad(n_q_f_points);
    std::vector<Tensor<1,2>> v_grad(n_q_f_points);

    double force = 0.0;

    // Loop over locally owned cells only
    for (const auto &cell : dof_handler_elastic.active_cell_iterators())
    {
        if (cell->is_locally_owned ())
        {
            const DoFHandler<2>::active_cell_iterator damage_cell =
                      Triangulation<2>::active_cell_iterator (cell)->as_dof_handler_iterator (
                          dof_handler_damage);
            const MaterialProperties material = material_properties(cell->material_id());

        for (unsigned int face_number = 0; face_number < GeometryInfo<2>::faces_per_cell; ++face_number)
        {
            if (cell->face(face_number)->at_boundary() &&
                cell->face(face_number)->boundary_id() == top_boundary_id)
            {
                fe_face_values.reinit(cell, face_number);
                fe_face_values_damage.reinit(damage_cell, face_number);

                // Get function values and gradients using correct variable names
                fe_face_values[u_x].get_function_values(locally_relevant_solution_elastic, u_val);
                fe_face_values[u_y].get_function_values(locally_relevant_solution_elastic, v_val);
                fe_face_values_damage.get_function_values(locally_relevant_solution_damage, phi_val);

                fe_face_values[u_x].get_function_gradients(locally_relevant_solution_elastic, u_grad);
                fe_face_values[u_y].get_function_gradients(locally_relevant_solution_elastic, v_grad);

                for (unsigned int q = 0; q < n_q_f_points; ++q)
                {
                    double ux = u_grad[q][0];
                    double uy = u_grad[q][1];
                    double vx = v_grad[q][0];
                    double vy = v_grad[q][1];
                    double phi = phi_val[q];

                    double exx = ux;
                    double eyy = vy;
                    double exy = 0.5 * (uy + vx);

                    // Uniform degradation on full stress (Miehe 2010, ref_source)
                    double degradation = (1.0 - phi) * (1.0 - phi);

                    double sxy = degradation * 2.0 * material.mu * exy;
                    double syy = degradation * (material.lambda * (exx + eyy)
                               + 2.0 * material.mu * eyy);

                    force += (sxy * fe_face_values.normal_vector(q)[0] +
                              syy * fe_face_values.normal_vector(q)[1]) * fe_face_values.JxW(q);
                }
            }
        }
    }


    }


    // Sum force contributions from all processors
    force = Utilities::MPI::sum(force, mpi_communicator);

    // Store results
    applied_displacement =boundary_value;
    force_data = force;
}
