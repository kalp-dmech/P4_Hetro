//
// Created by vijay on 25/2/25.
//

#ifndef VECTOR_POST_PROCESSOR_H
#define VECTOR_POST_PROCESSOR_H

#include "allheaders.h"
#include "Phase_Field.h"

class VectorPostProcessor : public DataPostprocessorVector<2>
{
  public:
    VectorPostProcessor();
    virtual
    void evaluate_vector_field (const DataPostprocessorInputs::Vector<2> &inputs,std::vector<Vector<double>> &computed_quantities ) const;
};

#endif //VECTOR_POST_PROCESSOR_H
