

#include "stress_post_processor.h"
#include "strain_post_processor.h"
#include "sed_data_post_processor.h"
#include "vector_post_processor.h"

#include "Phase_Field.h"

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

        Vector<double> subdomain(triangulation.n_active_cells());
        Vector<double> cell_level(triangulation.n_active_cells());
        for (const auto &cell : triangulation.active_cell_iterators())
            cell_level(cell->active_cell_index()) = cell->level();
        for (unsigned int i = 0; i < subdomain.size(); ++i)
            subdomain(i) = triangulation.locally_owned_subdomain();

        initial_phi_out.add_data_vector(cell_level, "cell_level");
        initial_phi_out.add_data_vector(subdomain, "subdomain");
        initial_phi_out.build_patches();
        initial_phi_out.write_vtu_with_pvtu_record(
            output_directory + "/", "initial_phi", step, mpi_communicator, 2, 8);
        return;
    }

    StressPostProcessor stress_tensor;
    StrainPostProcessor strain_tensor;
    SEDDataPostProcessor sed_data_post_processor;
    VectorPostProcessor solution1;


    DataOut<2> data_out_phasefield;

    //data_out_phasefield.attach_dof_handler(dof_handler_elastic);

    //data_out.add_data_vector(cell_id_averaged,"Averaged_strain_energy_density");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
            interpretation (2,DataComponentInterpretation::component_is_part_of_vector);

    std::vector<std::string> u_solution_names;
    u_solution_names.emplace_back("x_displacement");
    u_solution_names.emplace_back("y_displacement");
    data_out_phasefield.add_data_vector(dof_handler_elastic, locally_relevant_solution_elastic, u_solution_names, interpretation); //For outputting the u_solution
    data_out_phasefield.add_data_vector(dof_handler_damage,locally_relevant_solution_damage, "PhaseField");
    data_out_phasefield.add_data_vector(dof_handler_elastic,locally_relevant_solution_elastic,stress_tensor);
    data_out_phasefield.add_data_vector(dof_handler_elastic,locally_relevant_solution_elastic,strain_tensor);
    data_out_phasefield.add_data_vector(dof_handler_elastic,locally_relevant_solution_elastic,sed_data_post_processor);
    data_out_phasefield.add_data_vector(dof_handler_elastic,locally_relevant_solution_elastic,solution1);
    Vector<double> subdomain (triangulation.n_active_cells ());
    Vector<double> cell_level (triangulation.n_active_cells ());
    for (const auto &cell : triangulation.active_cell_iterators())
        cell_level(cell->active_cell_index()) = cell->level();
    for (unsigned int i = 0; i < subdomain.size (); ++i)
        subdomain (i) = triangulation.locally_owned_subdomain ();
    data_out_phasefield.add_data_vector (cell_level, "cell_level");
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
