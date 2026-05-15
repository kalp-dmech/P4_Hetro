
#include "Phase_Field.h"

#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

namespace
{
    void create_directory_recursive(const std::string &directory)
    {
        if (directory.empty())
            return;

        std::string current;
        for (const char character : directory)
        {
            current += character;
            if (character == '/')
            {
                if (current.size() > 1)
                {
                    const int mkdir_status = mkdir(current.c_str(), 0775);
                    AssertThrow(mkdir_status == 0 || errno == EEXIST,
                                ExcMessage("Could not create output directory."));
                }
            }
        }

        const int mkdir_status = mkdir(directory.c_str(), 0775);
        AssertThrow(mkdir_status == 0 || errno == EEXIST,
                    ExcMessage("Could not create output directory."));
    }

    void create_output_directory(MPI_Comm mpi_communicator,
                                 const std::string &directory)
    {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            create_directory_recursive(directory);
        MPI_Barrier(mpi_communicator);
    }

    std::string shell_quote(const std::string &value)
    {
        std::string quoted = "'";
        for (const char character : value)
        {
            if (character == '\'')
                quoted += "'\\''";
            else
                quoted += character;
        }
        quoted += "'";
        return quoted;
    }

    std::string format_force_value(const bool available,
                                   const double force_value)
    {
        if (!available)
            return "--";

        std::ostringstream out;
        out << std::fixed << std::setprecision(6) << force_value;
        return out.str();
    }
}


void PhaseField::run()
{
    Timer timer;
    timer.start();

    std::ostringstream material_1_line;
    material_1_line << std::fixed << std::setprecision(4)
                    << "Material 1 (steel): lambda="
                    << PhaseFieldVariables::material_1_lambda
                    << ", mu=" << PhaseFieldVariables::material_1_mu
                    << ", Gc=2.71e-3";

    std::ostringstream material_2_line;
    material_2_line << std::fixed << std::setprecision(4)
                    << "Material 2 (aluminum): lambda="
                    << PhaseFieldVariables::material_2_lambda
                    << ", mu=" << PhaseFieldVariables::material_2_mu
                    << ", Gc=8.0e-4";

    pcout << "\n============================================" << std::endl;
    pcout << "Phase-field fracture run (heterogeneous 2-material)" << std::endl;
    pcout << material_1_line.str() << std::endl;
    pcout << material_2_line.str() << std::endl;
    pcout << "Output directory: " << output_directory << std::endl;
    pcout << "Load steps: " << n_steps
          << " | strict tol = " << tol_1
          << " | relaxed tol = " << tol_2 << std::endl;
    pcout << "============================================" << std::endl;

    create_output_directory(mpi_communicator, output_directory);
    make_grid();
    setup_qph();
    setup_constraints_damage();

    completely_distributed_solution_damage = 0.0;
    completely_distributed_solution_damage_old = 0.0;
    locally_relevant_solution_damage = 0.0;

    output_results(0);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        std::cout << std::left
                  << std::setw(8) << "step" << "| "
                  << std::setw(12) << "disp" << "| "
                  << std::setw(22) << "staggered iterations" << "| "
                  << std::setw(12) << "AMR" << "| "
                  << std::setw(12) << "force" << "| "
                  << std::setw(12) << "cells" << std::endl;
        std::cout << std::string(88, '-') << std::endl;
    }

    double increment = 0;
    double increment1 = inc_small;
    double increment3 = inc_large;

    std::ofstream file;
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        file.open(force_displacement_file);
        AssertThrow(file, ExcMessage("Could not open force-displacement file."));
        file << "displacement,force" << std::endl;
        file << "0.0,0.0" << std::endl;
    }

    for (unsigned int i = 1; i <= n_steps; i++)
    {
        if (i <= 5)
            increment += increment3;
        else
            increment += increment1;

        double boundary_value = increment;

        setup_system_elastic(boundary_value);
        setup_constraints_damage();
        setup_system_damage();

        bool stoppingCriterion = false;
        unsigned int iteration = 0;
        unsigned int total_amr_passes = 0;
        const unsigned int cells_start_of_step = triangulation.n_global_active_cells();
        unsigned int displayed_cells = cells_start_of_step;
        std::string amr_report = "same";
        bool force_available = false;

        const auto print_progress_row =
            [&](const bool finalize_line)
            {
                if (Utilities::MPI::this_mpi_process(mpi_communicator) != 0)
                    return;
                std::ostringstream row;
                row << '\r'
                    << std::left
                    << std::setw(8) << i << "| "
                    << std::setw(12) << std::fixed << std::setprecision(6)
                    << boundary_value << "| "
                    << std::setw(22) << iteration << "| "
                    << std::setw(12) << amr_report << "| "
                    << std::setw(12) << format_force_value(force_available, force_data) << "| "
                    << std::setw(12) << displayed_cells;
                std::cout << row.str();
                if (finalize_line)
                    std::cout << std::endl;
                std::cout.flush();
            };

        print_progress_row(false);

        const unsigned int max_staggered_iterations =
            (i <= 5 ? max_staggered_iterations_1 : max_staggered_iterations_2);

        // AMR Refinement: Perform mesh refinement at the start of each load step
        // to ensure the grid is prepared for potential crack growth.
        if (enable_amr && 
            i % amr_frequency == 0 && 
            i <= amr_stop_step)
        {
            refine_grid(boundary_value);
            total_amr_passes++;
            displayed_cells = triangulation.n_global_active_cells();
            const int cell_change = static_cast<int>(displayed_cells) - static_cast<int>(cells_start_of_step);
            amr_report = (cell_change >= 0 ? "+" : "") + std::to_string(cell_change) + "/" + std::to_string(total_amr_passes);
            print_progress_row(false);
        }

        while (!stoppingCriterion && iteration < max_staggered_iterations)
        {
            assemble_system_elastic();
            solve_linear_system_elastic();
            locally_relevant_solution_elastic.update_ghost_values();

            assemble_system_damage();
            solve_linear_system_damage();
            locally_relevant_solution_damage.update_ghost_values();

            if (iteration > 0)
            {
                // Use relaxed tolerance (tol_2) if iterations exceed a threshold
                const double current_tol = (iteration > 50) ? tol_2 : tol_1;
                stoppingCriterion = check_convergence(current_tol);
            }

            completely_distributed_solution_elastic_old =
                locally_relevant_solution_elastic;
            completely_distributed_solution_damage_old =
                locally_relevant_solution_damage;

            if (!stoppingCriterion)
            {
                ++iteration;
            }
            else
                ++iteration;

            print_progress_row(false);
        }

        AssertThrow(stoppingCriterion,
                    ExcMessage("Staggered phase-field solve did not converge "
                               "within the configured iteration limit."));

        update_history_field();

        if (total_amr_passes == 0)
            amr_report = "same";
        displayed_cells = triangulation.n_global_active_cells();

        output_results(i);
        load_disp_calculation(boundary_value);
        force_available = true;
        print_progress_row(true);

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            file << applied_displacement << "," << force_data << std::endl;

        computing_timer.reset();

    } // for loop ends

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        file.close();
        const std::string comparison_plot =
            output_directory + "/force_displacement_comparison_hetero.png";
        const std::string plot_command =
            "python3 plot_force_displacement_comparison.py --generated " +
            shell_quote(force_displacement_file) + " --output " +
            shell_quote(comparison_plot);
        const int plot_status = std::system(plot_command.c_str());
        if (plot_status != 0)
            pcout << "Warning: force-displacement comparison plot failed."
                  << std::endl;
    }

    timer.stop();
    pcout << "\nRun finished in " << timer.wall_time() << " seconds." << std::endl;

} // run function
