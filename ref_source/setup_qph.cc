//
// Created by vijay on 23/2/25.
//

#include "Phase_Field.h"

void PhaseField::setup_qph()
{


    for (const auto &cell : triangulation.active_cell_iterators())
        {
             if (cell -> is_locally_owned())
            {
            quadrature_point_history_field.initialize(cell, quadrature_formula_damage.size());
            }
       }

    FEValues<2> fe_values_damage (fe_damage, quadrature_formula_damage,
        update_values | update_gradients | update_JxW_values
        | update_quadrature_points);

    for (const auto &cell : triangulation.active_cell_iterators ())
        if (cell->is_locally_owned ())
        {
            const std::vector<std::shared_ptr<MyQData>> lqph =
                quadrature_point_history_field.get_data (cell);
            for (const unsigned int q_index : fe_values_damage.quadrature_point_indices ())
            {
                lqph[q_index]->value_H = 0.0;
                lqph[q_index]->value_H_new = 0.0;
            }
        }
}