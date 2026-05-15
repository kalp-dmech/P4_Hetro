//
// Created by vijay on 25/2/25.
//

#ifndef SED_DATA_POST_PROCESSOR_H
#define SED_DATA_POST_PROCESSOR_H

#include "allheaders.h"
#include "Phase_Field.h"

class SEDDataPostProcessor: public DataPostprocessorScalar<2> // Inheriting the Datapostprocessorscalar class into user defined class
{
  public:

    SEDDataPostProcessor(); // Constructor of the class

    PhaseField constants; // To access the constants of the Elasticity class for doing strain energy density calculation.

    //double E = constants.E;
    //double nu = constants.nu;
    double mu = constants.mu;
    double lambda = constants.lambda;


    virtual         // Member function of the class.
    void evaluate_vector_field (const DataPostprocessorInputs::Vector<2> &inputs,
                                        std::vector<Vector<double>> &computed_quantities) const;
};



#endif //SED_DATA_POST_PROCESSOR_H
