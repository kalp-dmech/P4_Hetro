#include "Phase_Field.h"
#include <algorithm>

void PhaseField::compute_energies(double &elastic_energy,
                                  double &fracture_energy) const
{
    FEValues<2> fe_values_elastic(fe_elastic,
                                  quadrature_formula_elastic,
                                  update_gradients | update_JxW_values);
    FEValues<2> fe_values_damage(fe_damage,
                                 quadrature_formula_damage,
                                 update_values | update_gradients |
                                     update_JxW_values);

    AssertDimension(quadrature_formula_elastic.size(),
                    quadrature_formula_damage.size());

    const FEValuesExtractors::Vector displacements(0);
    const unsigned int n_q_points = quadrature_formula_damage.size();
    std::vector<SymmetricTensor<2, 2>> strain_values(n_q_points);
    std::vector<double> damage_values(n_q_points);
    std::vector<Tensor<1, 2>> damage_gradients(n_q_points);

    double local_elastic_energy = 0.0;
    double local_fracture_energy = 0.0;

    for (const auto &cell : dof_handler_elastic.active_cell_iterators())
        if (cell->is_locally_owned())
        {
            const DoFHandler<2>::active_cell_iterator damage_cell =
                Triangulation<2>::active_cell_iterator(cell)
                    ->as_dof_handler_iterator(dof_handler_damage);

            fe_values_elastic.reinit(cell);
            fe_values_damage.reinit(damage_cell);

            fe_values_elastic[displacements].get_function_symmetric_gradients(
                locally_relevant_solution_elastic, strain_values);
            fe_values_damage.get_function_values(locally_relevant_solution_damage,
                                                 damage_values);
            fe_values_damage.get_function_gradients(locally_relevant_solution_damage,
                                                    damage_gradients);
            const MaterialProperties material =
                material_properties(cell->material_id());

            for (const unsigned int q : fe_values_damage.quadrature_point_indices())
            {
                const double damage =
                    std::min(1.0, std::max(0.0, damage_values[q]));
                const double degradation =
                    (1.0 - damage) * (1.0 - damage) + k;
                const double tr_strain = trace(strain_values[q]);
                const double Mac_tr_minus = std::min(tr_strain, 0.0);

                // Spectral (Miehe 2010) compressive energy:
                // psi- = lambda/2 <tr(eps)>-^2 + mu * sum_i <eps_i>-^2
                const std::array<double, 2> eigs = eigenvalues(strain_values[q]);
                const double Min_e1 = std::min(eigs[0], 0.0);
                const double Min_e2 = std::min(eigs[1], 0.0);
                const double compressive_energy =
                    0.5 * material.lambda * Mac_tr_minus * Mac_tr_minus
                    + material.mu * (Min_e1 * Min_e1 + Min_e2 * Min_e2);

                local_elastic_energy +=
                    (degradation *
                         tensile_energy_density(strain_values[q], material) +
                     compressive_energy) *
                    fe_values_elastic.JxW(q);

                local_fracture_energy +=
                    material.Gc *
                    (0.5 / lo * damage * damage +
                     0.5 * lo *
                         (damage_gradients[q] * damage_gradients[q])) *
                    fe_values_damage.JxW(q);
            }
        }

    elastic_energy = Utilities::MPI::sum(local_elastic_energy, mpi_communicator);
    fracture_energy = Utilities::MPI::sum(local_fracture_energy, mpi_communicator);
}
