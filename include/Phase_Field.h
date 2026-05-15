#ifndef __MAIN_ALL_HEADER_H_INCLUDED__
#define __MAIN_ALL_HEADER_H_INCLUDED__

#include "allheaders.h"

namespace PhaseFieldVariables
{
    extern const double Gc;
    extern const double lo;
    extern const unsigned int material_1_id;
    extern const unsigned int material_2_id;
    extern const double material_1_lambda;
    extern const double material_1_mu;
    extern const double material_1_Gc;
    extern const double material_2_lambda;
    extern const double material_2_mu;
    extern const double material_2_Gc;
    extern const unsigned int bottom_boundary_id;
    extern const unsigned int top_boundary_id;
    extern const unsigned int crack_boundary_id;

    extern const unsigned int n_steps;
    extern const unsigned int max_staggered_iterations_1;
    extern const unsigned int max_staggered_iterations_2;
    extern const double tol_1;
    extern const double tol_2;
    extern const double inc_large;
    extern const double inc_small;
}

class PhaseField
{
public:
    PhaseField( );
    PhaseField(const std::string &output_directory,
               const std::string &force_displacement_file,
               const bool enable_amr);

    void run();
    double lambda, mu; 
    double Gc;
    double lo;
    double k;


private:
    void make_grid();
    void setup_qph();

    void setup_constraints_elastic(double boundary_value);
    void setup_system_elastic(double boundary_value);
    void assemble_system_elastic();
    void solve_linear_system_elastic();

    void setup_constraints_damage();
    void setup_system_damage();
    void assemble_system_damage();
    struct MaterialProperties
    {
        double lambda;
        double mu;
        double Gc;
    };

    MaterialProperties material_properties(const types::material_id material_id) const;
    double H_plus(const SymmetricTensor<2,2> &strain) const;
    double H_plus(const SymmetricTensor<2,2> &strain,
                  const MaterialProperties &material) const;
    double tensile_energy_density(const SymmetricTensor<2,2> &strain) const;
    double tensile_energy_density(const SymmetricTensor<2,2> &strain,
                                  const MaterialProperties &material) const;
    SymmetricTensor<2,2> tensile_stress(const SymmetricTensor<2,2> &strain) const;
    SymmetricTensor<2,2>
    tensile_stress(const SymmetricTensor<2,2> &strain,
                   const MaterialProperties &material) const;
    SymmetricTensor<2,2> compressive_stress(const SymmetricTensor<2,2> &strain) const;
    SymmetricTensor<2,2>
    compressive_stress(const SymmetricTensor<2,2> &strain,
                       const MaterialProperties &material) const;
    SymmetricTensor<2,2>
    degraded_stress(const SymmetricTensor<2,2> &strain,
                    const double damage) const;
    SymmetricTensor<2,2>
    degraded_stress(const SymmetricTensor<2,2> &strain,
                    const double damage,
                    const MaterialProperties &material) const;
    double elastic_tangent_product(const SymmetricTensor<2,2> &test_strain,
                                   const SymmetricTensor<2,2> &trial_strain,
                                   const double damage,
                                   const bool tensile_volumetric_active) const;
    double elastic_tangent_product(const SymmetricTensor<2,2> &test_strain,
                                   const SymmetricTensor<2,2> &trial_strain,
                                   const double damage,
                                   const bool tensile_volumetric_active,
                                   const MaterialProperties &material) const;
    void solve_linear_system_damage();

    bool check_convergence(const double tol);

    void output_results(const unsigned int step) const;
    void update_history_field();
    void load_disp_calculation(double boundary_value);
    void compute_energies(double &elastic_energy,
                          double &fracture_energy) const;
    void refine_grid(double boundary_value);

    MPI_Comm mpi_communicator;
    ConditionalOStream pcout;
    TimerOutput computing_timer;

    parallel::distributed::Triangulation<2> triangulation;

    FESystem<2>    fe_elastic;
    DoFHandler<2>    dof_handler_elastic;
    IndexSet locally_owned_dofs_elastic;
    IndexSet locally_relevant_dofs_elastic;
    AffineConstraints<double> constraints_elastic; 
    LA::MPI::SparseMatrix system_matrix_elastic;
    LA::MPI::Vector locally_relevant_solution_elastic;
    LA::MPI::Vector completely_distributed_solution_elastic_old;
    LA::MPI::Vector completely_distributed_solution_elastic;
    LA::MPI::Vector system_rhs_elastic;
    const QGauss<2> quadrature_formula_elastic;

    FESystem<2>    fe_damage;
    DoFHandler<2>  dof_handler_damage;
    IndexSet locally_owned_dofs_damage;
    IndexSet locally_relevant_dofs_damage;
    AffineConstraints<double> constraints_damage;
    LA::MPI::SparseMatrix system_matrix_damage;
    LA::MPI::Vector locally_relevant_solution_damage;
    LA::MPI::Vector completely_distributed_solution_damage_old;
    LA::MPI::Vector completely_distributed_solution_damage;
    LA::MPI::Vector system_rhs_damage;
    const QGauss<2> quadrature_formula_damage;

    const unsigned int bottom_boundary_id;
    const unsigned int top_boundary_id;
    const unsigned int crack_boundary_id;
    double applied_displacement;
    double force_data;
    unsigned int last_elastic_cg_iterations;
    unsigned int last_damage_cg_iterations;
    double last_elastic_residual;
    double last_damage_residual;

    const std::string output_directory;
    const std::string force_displacement_file;
    const bool enable_amr;
    const unsigned int n_steps;
    const unsigned int max_staggered_iterations_1;
    const unsigned int max_staggered_iterations_2;
    const unsigned int output_interval;
    const unsigned int amr_frequency;
    const unsigned int amr_stop_step;
    const unsigned int amr_max_level;
    const double inc_large;
    const double inc_small;
    const double tol_1;
    const double tol_2;
    const double amr_refine_length;
    Vector<double> ux_disp; 
    Vector<double> uy_disp; 

    class MyQData : public TransferableQuadraturePointData
    {
    public:
        MyQData() = default;
        virtual ~MyQData() = default;

        unsigned int number_of_values() const override
        {
            return 2;  
        }

        void pack_values(std::vector<double> &scalars) const override
        {
            Assert(scalars.size() == 2, ExcInternalError());
            scalars[0] = value_H;
            scalars[1] = value_H_new;
        }

        void unpack_values(const std::vector<double> &scalars) override
        {
            Assert(scalars.size() == 2, ExcInternalError());
            value_H        = scalars[0];
            value_H_new    = scalars[1];
        }

        double value_H;        
        double value_H_new;    
    };

    CellDataStorage<Triangulation<2>::cell_iterator, MyQData>
      quadrature_point_history_field;

};

#endif //__MAIN_ALL_HEADER_H_INCLUDED__
