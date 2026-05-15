#include "Phase_Field.h"

void PhaseField::setup_constraints_damage()

{
    TimerOutput::Scope ts (computing_timer, "setup_bv_damage");

    constraints_damage.clear ();
    //constraints_damage.reinit (locally_owned_dofs_damage,locally_relevant_dofs_damage);
    constraints_damage.reinit (locally_relevant_dofs_damage);
    DoFTools::make_hanging_node_constraints (dof_handler_damage,
        constraints_damage);

    if (crack_boundary_id != numbers::invalid_boundary_id)
        VectorTools::interpolate_boundary_values(dof_handler_damage,
                                                  crack_boundary_id,
                                                  Functions::ConstantFunction<2>(1,1),
                                                    constraints_damage
                                                  ); //A mask that indicates which components of the solution vector should have the boundary condition applied

    constraints_damage.close ();

}
