#include "allheaders.h"
#include "Phase_Field.h"

#include <cstdlib>
#include <string>

namespace
{
    bool smoke_test_enabled()
    {
        const char *value = std::getenv("PF_SMOKE_TEST");
        return value != nullptr && std::string(value) != "0";
    }
}

int
main (int argc,
      char *argv[])
{
    try
    {
        using namespace dealii;
        Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

        const bool enable_amr = true;
        const std::string output_directory =
            smoke_test_enabled() ? std::string("output_smoke") : std::string("output");

        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
            std::cout << "\nRuntime configuration" << std::endl;
            std::cout << "  MPI cores: "
                      << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl;
            std::cout << "  AMR: "
                      << (enable_amr ? "enabled" : "disabled")
                      << std::endl;
            std::cout << "  Output directory: " << output_directory << std::endl;
        }

        const std::string force_displacement_file =
            output_directory + "/force_displacement.csv";

        PhaseField phasefield(output_directory,
                              force_displacement_file,
                              enable_amr);
        phasefield.run();

        MPI_Barrier(MPI_COMM_WORLD);
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------" << std::endl;
        std::cerr << "Exception on processing: " << std::endl << exc.what ()
        << std::endl << "Aborting!" << std::endl
        << "----------------------------------------------------" << std::endl;

        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------" << std::endl;
        std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
        << "----------------------------------------------------" << std::endl;
        return 1;
    }

    return 0;

}
