#include "Phase_Field.h"

void PhaseField::setup_constraints_elastic(double boundary_value)
{
    constraints_elastic.clear ();
    //constraints_elastic.reinit (locally_owned_dofs_elastic,locally_relevant_dofs_elastic);
    constraints_elastic.reinit (locally_relevant_dofs_elastic);
    DoFTools::make_hanging_node_constraints (dof_handler_elastic,
       constraints_elastic);

    const FEValuesExtractors::Scalar u_x (0);
    const FEValuesExtractors::Scalar u_y (1);

    const ComponentMask u_x_mask = fe_elastic.component_mask (u_x);
    const ComponentMask u_y_mask = fe_elastic.component_mask (u_y);

    VectorTools::interpolate_boundary_values(dof_handler_elastic,
                                             bottom_boundary_id,
                                             Functions::ZeroFunction<2>(2),
                                             constraints_elastic,
                                             u_x_mask); //A mask that indicates which components of the solution vector should have the boundary condition applied

    VectorTools::interpolate_boundary_values(dof_handler_elastic,
                                             bottom_boundary_id,
                                             Functions::ZeroFunction<2>(2),
                                             constraints_elastic,
                                             u_y_mask); //A mask that indicates which components of the solution vector should have the boundary condition applied


    VectorTools::interpolate_boundary_values(dof_handler_elastic,
                                             top_boundary_id,
                                             Functions::ConstantFunction<2>(boundary_value,2), // <2> -- 2d problem (1.0e-5 --> boundayr value, 2 --> no of dof)
                                             constraints_elastic,
                                             u_y_mask);

    constraints_elastic.close();
}
