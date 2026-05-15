//
// Created by vijay on 25/2/25.
//

#ifndef STRAIN_POST_PROCESS_H
#define STRAIN_POST_PROCESS_H

#include "allheaders.h"
#include "Phase_Field.h"

class StrainPostProcessor: public DataPostprocessorTensor<2>
{
public:
    StrainPostProcessor ();

    PhaseField constants;

    //double E = constants.E;
    //double nu = constants.nu;
    double lambda = constants.lambda;
    double mu = constants.mu;



    virtual
    void evaluate_vector_field (const DataPostprocessorInputs::Vector<2> &inputs, std::vector<Vector<double>> &computed_quantities ) const;
};



#endif //STRAIN_POST_PROCESS_H
