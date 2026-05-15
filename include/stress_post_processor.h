//
// Created by vijay on 25/2/25.
//

#ifndef STRESS_POST_PROCESS_H
#define STRESS_POST_PROCESS_H

#include "allheaders.h"
#include "Phase_Field.h"

class StressPostProcessor: public DataPostprocessorTensor<2>
{
public:
    StressPostProcessor ();

    PhaseField constants;

    //double E = constants.E;
    //double nu = constants.nu;
    double lambda = constants.lambda;
    double mu = constants.mu;



    virtual
    void evaluate_vector_field (const DataPostprocessorInputs::Vector<2> &inputs, std::vector<Vector<double>> &computed_quantities ) const;
};



#endif //STRESS_POST_PROCESS_H
