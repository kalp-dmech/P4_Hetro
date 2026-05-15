//
// Created by vijay on 6/3/25.
//

#include "strain_post_processor.h"

StrainPostProcessor::StrainPostProcessor()
    :
    DataPostprocessorTensor<2>("strain_tensor",update_gradients|update_values)
{ }

void StrainPostProcessor::evaluate_vector_field(const DataPostprocessorInputs::Vector<2> &inputs,
                                                            std::vector<Vector<double> > &computed_quantities) const
{
    Assert(computed_quantities.size() == inputs.solution_gradients.size(),
        ExcDimensionMismatch(computed_quantities.size(), inputs.solution_gradients.size()));


    for (unsigned int p = 0; p < inputs.solution_gradients.size(); ++p)
    {

        AssertDimension(computed_quantities[p].size(), (Tensor<2,2>::n_independent_components));


        double dux_dx = inputs.solution_gradients[p][0][0]; //  dux_dx = dux/dx --> derivative of ux with respect to x
        double dux_dy = inputs.solution_gradients[p][0][1];
        double duy_dx = inputs.solution_gradients[p][1][0];
        double duy_dy = inputs.solution_gradients[p][1][1];

        double exx = dux_dx;
        double eyy = duy_dy;
        double exy = 0.5*(dux_dy+duy_dx);
        double eyx = exy;

        
        computed_quantities[p][Tensor<2,2>::component_to_unrolled_index(TableIndices<2>(0,0))]=exx;
        computed_quantities[p][Tensor<2,2>::component_to_unrolled_index(TableIndices<2>(0,1))]=exy;
        computed_quantities[p][Tensor<2,2>::component_to_unrolled_index(TableIndices<2>(1,0))]=eyx;
        computed_quantities[p][Tensor<2,2>::component_to_unrolled_index(TableIndices<2>(1,1))]=eyy;


    }

}
