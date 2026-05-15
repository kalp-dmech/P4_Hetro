

#include "strain_post_processor.h"
#include "vector_post_processor.h"

#include "Phase_Field.h"

#include <algorithm>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>


void PhaseField::output_results(const unsigned int step) const
{
    if (step != 0 && step != 1 && step != n_steps &&
        step % output_interval != 0)
        return;

    const int mkdir_status = mkdir(output_directory.c_str(), 0775);
    AssertThrow(mkdir_status == 0 || errno == EEXIST,
                ExcMessage("Could not create output directory in current working directory."));

    if (step == 0)
    {
        DataOut<2> initial_phi_out;
        initial_phi_out.add_data_vector(dof_handler_damage,
                                        locally_relevant_solution_damage,
                                        "PhaseField");

        LA::MPI::Vector initial_crack_indicator(locally_owned_dofs_damage,
                                                mpi_communicator);
        for (const auto index : locally_owned_dofs_damage)
            initial_crack_indicator[index] =
                (locally_relevant_solution_damage[index] >= 0.99 ? 1.0 : 0.0);
        initial_crack_indicator.compress(VectorOperation::insert);

        initial_phi_out.add_data_vector(dof_handler_damage,
                                        initial_crack_indicator,
                                        "crack_indicator");

        Vector<double> subdomain(triangulation.n_active_cells());
        Vector<double> cell_level(triangulation.n_active_cells());
        Vector<double> material_id(triangulation.n_active_cells());
        for (const auto &cell : triangulation.active_cell_iterators())
        {
            cell_level(cell->active_cell_index()) = cell->level();
            material_id(cell->active_cell_index()) = cell->material_id();
        }
        for (unsigned int i = 0; i < subdomain.size(); ++i)
            subdomain(i) = triangulation.locally_owned_subdomain();

        initial_phi_out.add_data_vector(cell_level, "cell_level");
        initial_phi_out.add_data_vector(material_id, "material_id");
        initial_phi_out.add_data_vector(subdomain, "subdomain");
        initial_phi_out.build_patches();
        initial_phi_out.write_vtu_with_pvtu_record(
            output_directory + "/", "initial_phi", step, mpi_communicator, 2, 8);
        return;
    }

    StrainPostProcessor strain_tensor;
    VectorPostProcessor solution1;

    // Compute cell-averaged degraded stress and strain energy density
    // using correct per-material properties and damage field
    const unsigned int n_active = triangulation.n_active_cells();
    Vector<double> cell_sigma_xx(n_active);
    Vector<double> cell_sigma_yy(n_active);
    Vector<double> cell_sigma_xy(n_active);
    Vector<double> cell_sed(n_active);

    {
        FEValues<2> fe_vals_e(fe_elastic, quadrature_formula_elastic,
                              update_gradients | update_JxW_values);
        FEValues<2> fe_vals_d(fe_damage, quadrature_formula_damage,
                              update_values);

        const FEValuesExtractors::Vector displacements(0);
        const unsigned int n_q = quadrature_formula_elastic.size();
        std::vector<SymmetricTensor<2, 2>> strain_vals(n_q);
        std::vector<double> damage_vals(n_q);

        for (const auto &cell : dof_handler_elastic.active_cell_iterators())
            if (cell->is_locally_owned())
            {
                const auto damage_cell =
                    Triangulation<2>::active_cell_iterator(cell)
                        ->as_dof_handler_iterator(dof_handler_damage);

                fe_vals_e.reinit(cell);
                fe_vals_d.reinit(damage_cell);

                const MaterialProperties mat =
                    material_properties(cell->material_id());

                fe_vals_e[displacements].get_function_symmetric_gradients(
                    locally_relevant_solution_elastic, strain_vals);
                fe_vals_d.get_function_values(
                    locally_relevant_solution_damage, damage_vals);

                SymmetricTensor<2, 2> avg_stress;
                double avg_sed = 0.0;
                double total_weight = 0.0;

                for (unsigned int q = 0; q < n_q; ++q)
                {
                    const double w = fe_vals_e.JxW(q);
                    const SymmetricTensor<2, 2> sig =
                        degraded_stress(strain_vals[q], damage_vals[q], mat);

                    avg_stress += w * sig;
                    avg_sed    += w * (sig * strain_vals[q]);
                    total_weight += w;
                }

                avg_stress /= total_weight;
                avg_sed    /= total_weight;

                const unsigned int idx = cell->active_cell_index();
                cell_sigma_xx(idx) = avg_stress[0][0];
                cell_sigma_yy(idx) = avg_stress[1][1];
                cell_sigma_xy(idx) = avg_stress[0][1];
                cell_sed(idx)      = 0.5 * avg_sed;
            }
    }

    DataOut<2> data_out_phasefield;

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
            interpretation (2,DataComponentInterpretation::component_is_part_of_vector);

    std::vector<std::string> u_solution_names;
    u_solution_names.emplace_back("x_displacement");
    u_solution_names.emplace_back("y_displacement");
    data_out_phasefield.add_data_vector(dof_handler_elastic, locally_relevant_solution_elastic, u_solution_names, interpretation);
    data_out_phasefield.add_data_vector(dof_handler_damage,locally_relevant_solution_damage, "PhaseField");

    LA::MPI::Vector crack_indicator(locally_owned_dofs_damage,
                                    mpi_communicator);
    for (const auto index : locally_owned_dofs_damage)
        crack_indicator[index] =
            (locally_relevant_solution_damage[index] >= 0.99 ? 1.0 : 0.0);
    crack_indicator.compress(VectorOperation::insert);

    data_out_phasefield.add_data_vector(dof_handler_damage,
                                        crack_indicator,
                                        "crack_indicator");

    data_out_phasefield.add_data_vector(dof_handler_elastic,locally_relevant_solution_elastic,strain_tensor);
    data_out_phasefield.add_data_vector(dof_handler_elastic,locally_relevant_solution_elastic,solution1);

    data_out_phasefield.add_data_vector(cell_sigma_xx, "sigma_xx");
    data_out_phasefield.add_data_vector(cell_sigma_yy, "sigma_yy");
    data_out_phasefield.add_data_vector(cell_sigma_xy, "sigma_xy");
    data_out_phasefield.add_data_vector(cell_sed, "strain_energy_density");

    Vector<double> subdomain (n_active);
    Vector<double> cell_level (n_active);
    Vector<double> material_id (n_active);
    for (const auto &cell : triangulation.active_cell_iterators())
    {
        cell_level(cell->active_cell_index()) = cell->level();
        material_id(cell->active_cell_index()) = cell->material_id();
    }
    for (unsigned int i = 0; i < subdomain.size (); ++i)
        subdomain (i) = triangulation.locally_owned_subdomain ();
    data_out_phasefield.add_data_vector (cell_level, "cell_level");
    data_out_phasefield.add_data_vector (material_id, "material_id");
    data_out_phasefield.add_data_vector (subdomain, "subdomain");
    data_out_phasefield.build_patches ();


    data_out_phasefield.write_vtu_with_pvtu_record(
      output_directory + "/", "solution",  step, mpi_communicator,2, 8);
   // std::cout << "output  for u is finished" << std::endl;

   /* DataOut<2> data_out_phasefield1;

    data_out_phasefield1.attach_dof_handler(dof_handler_damage);

    data_out_phasefield1.add_data_vector(dof_handler_damage,locally_relevant_solution_damage, "PhaseField");
    data_out_phasefield1.build_patches ();


    data_out_phasefield1.write_vtu_with_pvtu_record(
      "output/", "solution_phi",  step, mpi_communicator,2, 8);
    // std::cout << "output  for u is finished" << std::endl;
*/


}
