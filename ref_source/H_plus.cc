#include "Phase_Field.h"


double
PhaseField::H_plus (const SymmetricTensor<2, 2> &strain)
{
    double Mac_tr_strain, Mac_first_principal_strain,
        Mac_second_principal_strain,
        tr_sqr_Mac_Principal_strain;

    const double tr_strain = trace (strain);

    Mac_tr_strain = tr_strain >0 ? tr_strain : 0;
    const std::array<double, 2> Principal_strains = eigenvalues (strain);

    Mac_first_principal_strain = (Principal_strains[0] > 0) ? Principal_strains[0] : 0;
    Mac_second_principal_strain = (Principal_strains[1] > 0) ? Principal_strains[1] : 0;


    tr_sqr_Mac_Principal_strain = pow (Mac_first_principal_strain, 2)
        + pow (Mac_second_principal_strain, 2);

    double H_plus_val;
    H_plus_val = 0.5 * lambda * pow (Mac_tr_strain, 2)
        + mu  * tr_sqr_Mac_Principal_strain;
    return H_plus_val;
}