 //
// Created by vijay on 23/2/25.
//

#include "Phase_Field.h"

 void
  PhaseField::update_history_field ()
{
	FEValues<2> fe_values_damage (fe_damage, quadrature_formula_damage,
		update_values | update_gradients | update_JxW_values
		| update_quadrature_points);

	for (const auto &cell : dof_handler_damage.active_cell_iterators ())
		if (cell->is_locally_owned ())
		{
			const std::vector<std::shared_ptr<MyQData>> lqph =
				quadrature_point_history_field.get_data (cell);
			for (unsigned int q_index = 0;
				q_index < quadrature_formula_damage.size (); ++q_index)
				lqph[q_index]->value_H = lqph[q_index]->value_H_new;
		}
}
