#include "Phase_Field.h"

namespace PhaseFieldVariables
{
    // Core simulation constants.
    const double lo = 0.0075;
    const double Gc = 2.71e-3;

    const unsigned int material_1_id = 1;
    const unsigned int material_2_id = 2;
    
    // steel
    const double material_1_lambda = 121.1538;
    const double material_1_mu     = 80.7692;
    const double material_1_Gc     = 2.71e-3;

    // aluminum
    const double material_2_lambda = 53.8462;
    const double material_2_mu     = 26.3158;
    const double material_2_Gc     = 8.0e-4;

    const unsigned int bottom_boundary_id = 1;
    const unsigned int top_boundary_id = 3;
    const unsigned int crack_boundary_id = 5;

    const unsigned int n_steps = 100;
    const unsigned int max_staggered_iterations_1 = 1000 ;
    const unsigned int max_staggered_iterations_2 = 250;
    const double tol_1 = 1.0e-3;
    const double tol_2 = 1.0e-1;

    // Displacement increment sizes:
    // inc_large is used for the first 5 load steps (coarse ramp-up),
    // inc_small is used for all subsequent steps (fine crack-propagation phase).
    const double inc_large = 1.0e-3; // mm, large increment for initial steps (steps 1-5)
    const double inc_small = 1.0e-5; // mm, small increment for crack propagation (steps 6+)
}
