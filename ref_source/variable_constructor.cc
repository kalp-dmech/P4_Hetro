#include "Phase_Field.h"

namespace PhaseFieldVariables
{
    // These variables are declared as parameters and can be modified via ParameterHandler
    const std::vector<unsigned int> pore_counts = {0,20,60,90,95};
    const unsigned int pore_trial_number = 1;
    const unsigned int case_number = 0;
    const unsigned int random_seed_base = 42;
    const double min_pore_size = 0.01;
    const double max_pore_size = 0.035;
    const double crack_spread_radius = 0.45;
    const double boundary_margin = 0.001;
    const double lo = 0.0075;
    const double Gc = 2.7e-3;

    const unsigned int n_steps = 200;
    const unsigned int max_staggered_iterations_1 = 5000;
    const unsigned int max_staggered_iterations_2 = 250;
    const double tol_1 = 1.0e-3;
    const double tol_2 = 1.0e-2;
}