#include "Phase_Field.h"

namespace PhaseFieldVariables
{
    void declare_parameters(ParameterHandler &prm)
    {
        prm.declare_entry("pore_count", std::to_string(pore_counts.front()),
                          Patterns::Integer(0));
        prm.declare_entry("pore_trial_number",
                          std::to_string(pore_trial_number),
                          Patterns::Integer(0));
        prm.declare_entry("case_number", std::to_string(case_number),
                          Patterns::Integer(0));
        prm.declare_entry("random_seed_base", std::to_string(random_seed_base),
                          Patterns::Integer(0));
        prm.declare_entry("min_pore_size", std::to_string(min_pore_size),
                          Patterns::Double(0.0));
        prm.declare_entry("max_pore_size", std::to_string(max_pore_size),
                          Patterns::Double(0.0));
        prm.declare_entry("crack_spread_radius",
                          std::to_string(crack_spread_radius),
                          Patterns::Double(0.0));
        prm.declare_entry("boundary_margin", std::to_string(boundary_margin),
                          Patterns::Double(0.0));
        prm.declare_entry("lo", std::to_string(lo), Patterns::Double(0.0));
        prm.declare_entry("Gc", std::to_string(Gc), Patterns::Double(0.0));

        prm.declare_entry("n_steps", std::to_string(n_steps),
                          Patterns::Integer(1));
        prm.declare_entry("tol_1", std::to_string(tol_1),
                          Patterns::Double(0.0));
        prm.declare_entry("tol_2", std::to_string(tol_2),
                          Patterns::Double(0.0));
        prm.declare_entry("max_staggered_iterations_1",
                          std::to_string(max_staggered_iterations_1),
                          Patterns::Integer(1));
        prm.declare_entry("max_staggered_iterations_2",
                          std::to_string(max_staggered_iterations_2),
                          Patterns::Integer(1));
    }
}
