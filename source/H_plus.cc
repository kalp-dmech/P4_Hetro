#include "Phase_Field.h"
#include <algorithm>

namespace
{
    constexpr unsigned int dim = 2;

    double macaulay_plus(const double value)
    {
        return std::max(value, 0.0);
    }

    double macaulay_minus(const double value)
    {
        return std::min(value, 0.0);
    }
}

PhaseField::MaterialProperties
PhaseField::material_properties(const types::material_id material_id) const
{
    if (material_id == PhaseFieldVariables::material_2_id)
        return {PhaseFieldVariables::material_2_lambda,
                PhaseFieldVariables::material_2_mu,
                PhaseFieldVariables::material_2_Gc};

    return {PhaseFieldVariables::material_1_lambda,
            PhaseFieldVariables::material_1_mu,
            PhaseFieldVariables::material_1_Gc};
}

double PhaseField::tensile_energy_density(
    const SymmetricTensor<2, 2> &strain) const
{
    return tensile_energy_density(
        strain, material_properties(PhaseFieldVariables::material_1_id));
}

double PhaseField::tensile_energy_density(
    const SymmetricTensor<2, 2> &strain,
    const MaterialProperties &material) const
{
    const double tr_strain = trace(strain);
    const double Mac_tr_strain = macaulay_plus(tr_strain);

    // Spectral decomposition (Miehe 2010):
    // psi+ = lambda/2 <tr(eps)>+^2 + mu * sum_i <eps_i>+^2
    const std::array<double, dim> eigs = eigenvalues(strain);
    const double Mac_e1 = macaulay_plus(eigs[0]);
    const double Mac_e2 = macaulay_plus(eigs[1]);

    return 0.5 * material.lambda * Mac_tr_strain * Mac_tr_strain
           + material.mu * (Mac_e1 * Mac_e1 + Mac_e2 * Mac_e2);
}

SymmetricTensor<2, 2>
PhaseField::tensile_stress(const SymmetricTensor<2, 2> &strain) const
{
    return tensile_stress(
        strain, material_properties(PhaseFieldVariables::material_1_id));
}

SymmetricTensor<2, 2>
PhaseField::tensile_stress(const SymmetricTensor<2, 2> &strain,
                           const MaterialProperties &material) const
{
    const SymmetricTensor<2, 2> identity = unit_symmetric_tensor<2>();
    const double tr_strain = trace(strain);

    // Spectral decomposition (Miehe 2010):
    // sigma+ = lambda <tr(eps)>+ I + 2 mu eps+
    // where eps+ = sum_i <eps_i>+ n_i (x) n_i
    const auto eig = eigenvectors(strain);
    SymmetricTensor<2, 2> strain_plus;
    for (unsigned int i = 0; i < dim; ++i)
    {
        const double ev = macaulay_plus(eig[i].first);
        strain_plus += ev * symmetrize(outer_product(eig[i].second,
                                                     eig[i].second));
    }

    return material.lambda * macaulay_plus(tr_strain) * identity +
           2.0 * material.mu * strain_plus;
}

SymmetricTensor<2, 2>
PhaseField::compressive_stress(const SymmetricTensor<2, 2> &strain) const
{
    return compressive_stress(
        strain, material_properties(PhaseFieldVariables::material_1_id));
}

SymmetricTensor<2, 2>
PhaseField::compressive_stress(const SymmetricTensor<2, 2> &strain,
                               const MaterialProperties &material) const
{
    const SymmetricTensor<2, 2> identity = unit_symmetric_tensor<2>();
    const double tr_strain = trace(strain);

    // Spectral decomposition (Miehe 2010):
    // sigma- = lambda <tr(eps)>- I + 2 mu eps-
    // where eps- = sum_i <eps_i>- n_i (x) n_i
    const auto eig = eigenvectors(strain);
    SymmetricTensor<2, 2> strain_minus;
    for (unsigned int i = 0; i < dim; ++i)
    {
        const double ev = macaulay_minus(eig[i].first);
        strain_minus += ev * symmetrize(outer_product(eig[i].second,
                                                      eig[i].second));
    }

    return material.lambda * macaulay_minus(tr_strain) * identity +
           2.0 * material.mu * strain_minus;
}

SymmetricTensor<2, 2>
PhaseField::degraded_stress(const SymmetricTensor<2, 2> &strain,
                            const double damage) const
{
    return degraded_stress(
        strain, damage, material_properties(PhaseFieldVariables::material_1_id));
}

SymmetricTensor<2, 2>
PhaseField::degraded_stress(const SymmetricTensor<2, 2> &strain,
                            const double damage,
                            const MaterialProperties &material) const
{
    const double bounded_damage = std::min(1.0, std::max(0.0, damage));
    const double degradation =
        (1.0 - bounded_damage) * (1.0 - bounded_damage) + k;

    // Spectral split: sigma = g * sigma+ + sigma-
    return degradation * tensile_stress(strain, material) +
           compressive_stress(strain, material);
}

double
PhaseField::elastic_tangent_product(const SymmetricTensor<2, 2> &test_strain,
                                    const SymmetricTensor<2, 2> &trial_strain,
                                    const double damage,
                                    const bool /*tensile_volumetric_active*/) const
{
    return elastic_tangent_product(
        test_strain,
        trial_strain,
        damage,
        false,
        material_properties(PhaseFieldVariables::material_1_id));
}

double
PhaseField::elastic_tangent_product(const SymmetricTensor<2, 2> &test_strain,
                                    const SymmetricTensor<2, 2> &trial_strain,
                                    const double damage,
                                    const bool /*tensile_volumetric_active*/,
                                    const MaterialProperties &material) const
{
    const double bounded_damage = std::min(1.0, std::max(0.0, damage));
    const double degradation =
        (1.0 - bounded_damage) * (1.0 - bounded_damage) + k;

    // Uniform degradation on full isotropic tangent (ref_source approach)
    return degradation * (material.lambda * trace(test_strain) * trace(trial_strain)
                         + 2.0 * material.mu * (test_strain * trial_strain));
}

double PhaseField::H_plus(const SymmetricTensor<2, 2> &strain) const
{
    return tensile_energy_density(strain);
}

double PhaseField::H_plus(const SymmetricTensor<2, 2> &strain,
                          const MaterialProperties &material) const
{
    return tensile_energy_density(strain, material);
}
