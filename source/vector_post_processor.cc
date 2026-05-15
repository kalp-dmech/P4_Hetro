//
// Created by vijay on 6/3/25.
//

#include "vector_post_processor.h"

VectorPostProcessor::VectorPostProcessor()
    :
    DataPostprocessorVector<2>("solution1",update_values )
{}

void VectorPostProcessor::evaluate_vector_field(const DataPostprocessorInputs::Vector<2> &inputs,
                                                std::vector<Vector<double> > &computed_quantities) const
{

    Assert(computed_quantities.size() == inputs.solution_values.size(),
            ExcDimensionMismatch (computed_quantities.size(), inputs.solution_values.size()));




    for (unsigned int p=0;p<inputs.solution_values.size();++p)
    {
        AssertDimension (computed_quantities[p].size(),2);

        double ux_val=inputs.solution_values[p][0];
        double uy_val=inputs.solution_values[p][1];

        computed_quantities[p][0]=ux_val;
        computed_quantities[p][1]=uy_val;
    }

}
