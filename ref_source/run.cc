
#include "Phase_Field.h"

#include <cerrno>
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

    std::string comparison_plot_file_name(const unsigned int pore_count,
                                          const unsigned int trial_number)
    {
        return "force_displacement_comparison_Pore_" +
               std::to_string(pore_count) + "_trail_" +
               std::to_string(trial_number) + ".svg";
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
    timer.start ();
    const unsigned int total_pore_cases =
        PhaseFieldVariables::pore_counts.size();
    unsigned int pore_case_number = 1;
    for (unsigned int i = 0; i < total_pore_cases; ++i)
        if (PhaseFieldVariables::pore_counts[i] == pore_count)
        {
            pore_case_number = i + 1;
            break;
        }

    // ============================================================
    // 1. INITIAL SETUP & HEADER
    // ============================================================
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        std::cout << "\n============================================" << std::endl;
        std::cout << "Phase-field fracture run" << std::endl;
        std::cout << "Case: pore set " << pore_case_number
                  << "/" << total_pore_cases
                  << " | pores = " << pore_count
                  << " | trial = " << pore_trial_number
                  << "/" << total_trial_count << std::endl;
        std::cout << "Output directory: " << output_directory << std::endl;
        std::cout << "Load steps: " << n_steps
                  << " | strict tol = " << tol_1
                  << " | relaxed tol = " << tol_2 << std::endl;
        std::cout << "============================================" << std::endl;
    }

    create_output_directory(mpi_communicator, output_directory);
    make_grid();
    setup_qph();
    setup_constraints_damage();
    
    // Sustaining Gaussian Porosity Implementation
    initialize_random_porosity();
    pcout << "Initial damage field max(phi) = "
          << locally_relevant_solution_damage.linfty_norm() << std::endl;
    output_results(0);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        std::cout << std::left
                  << std::setw(8) << "step" << "| "
                  << std::setw(22) << "staggered iterations" << "| "
                  << std::setw(12) << "AMR" << "| "
                  << std::setw(12) << "force" << "| "
                  << std::setw(12) << "cells" << std::endl;
        std::cout << std::string(74, '-') << std::endl;
    }

    double increment = 0;
    double increment1 = inc_small;   // 1.0e-5 in variable_constructor.cc
    double increment3 = inc_large;   // 1.0e-3 in variable_constructor.cc

    std::ofstream file;
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        file.open(force_displacement_file);
        AssertThrow(file, ExcMessage("Could not open force-displacement file."));
        file << "displacement,force" << std::endl;
    }

    for (unsigned int i = 1; i <= n_steps; i++)
    {
       if (i <= 5)
        {
            increment = increment + increment3;
        }
        else
        {
            increment = increment + increment1;
        }

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

        while (stoppingCriterion == false)
        {
            assemble_system_elastic();
            solve_linear_system_elastic();
            locally_relevant_solution_elastic.update_ghost_values();

            assemble_system_damage();
            solve_linear_system_damage();
            locally_relevant_solution_damage.update_ghost_values();

            if (iteration > 0)
                stoppingCriterion = check_convergence ();

            completely_distributed_solution_elastic_old =
                locally_relevant_solution_elastic;
            completely_distributed_solution_damage_old =
                locally_relevant_solution_damage;

            if (stoppingCriterion == false)
            {
                ++iteration;
                refine_grid(boundary_value);
                total_amr_passes++;
                displayed_cells = triangulation.n_global_active_cells();
                const int cell_change =
                    static_cast<int>(displayed_cells) -
                    static_cast<int>(cells_start_of_step);
                amr_report = (cell_change >= 0 ? "+" : "") +
                             std::to_string(cell_change) + "/" +
                             std::to_string(total_amr_passes);
                print_progress_row(false);
            }

            if (stoppingCriterion == true)
                ++iteration;

            print_progress_row(false);
        } // while loop for converged solution ends

        update_history_field();

        const unsigned int cells_end_of_step = triangulation.n_global_active_cells();
        displayed_cells = cells_end_of_step;
        if (total_amr_passes == 0)
            amr_report = "same";

        output_results(i);
        load_disp_calculation(boundary_value);
        force_available = true;
        print_progress_row(true);
        
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            file << applied_displacement << "," << force_data << std::endl;

        computing_timer.reset ();

    } // for loop ends

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        file.close();
        const std::string comparison_plot =
            output_directory + "/" +
            comparison_plot_file_name(pore_count, pore_trial_number);
        const std::string plot_command =
            "cd /home/me24m086/PFM_final/P3 && "
            "python3 plot_force_displacement_comparison.py --generated " +
            shell_quote(force_displacement_file) + " --output " +
            shell_quote(comparison_plot);
        const int plot_status = std::system(plot_command.c_str());
        if (plot_status != 0)
            pcout << "Warning: force-displacement comparison plot failed."
                  << std::endl;
    }

    timer.stop ();
    pcout << "\nRun finished in " << timer.wall_time() << " seconds."
          << std::endl;

} // run function
