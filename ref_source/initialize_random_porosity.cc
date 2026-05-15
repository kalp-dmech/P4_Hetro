// Initial random porosity — mirrors the P11 reference pattern.
//
// Responsibilities:
//   1. Place the requested number of pore centers/radii (Gaussian cloud,
//      excluded from the crack boundary and from other pores).
//   2. Write the diffuse pore damage field directly into the damage
//      solution vector (no separate floor field).  After this call
//      `completely_distributed_solution_damage`,
//      `completely_distributed_solution_damage_old`, and
//      `locally_relevant_solution_damage` all carry the pore shape and
//      irreversibility is preserved by the nodal clamp inside run().
//   3. Lock the history variable H = Gc/(2*lo) at quadrature points that
//      lie strongly inside a pore (phi_q > 0.5).  This is the canonical
//      Francfort-Marigo driving force, not an inflated multiple.

#include "Phase_Field.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace
{
    struct Bounds
    {
        Point<2> min;
        Point<2> max;
    };

    std::vector<Point<2>> unpack_points(const std::vector<double> &buffer)
    {
        std::vector<Point<2>> points;
        points.reserve(buffer.size() / 2);

        for (unsigned int i = 0; i + 1 < buffer.size(); i += 2)
            points.emplace_back(buffer[i], buffer[i + 1]);

        return points;
    }
}

double PhaseField::compute_diffuse_phi(const Point<2> &point) const
{
    double phi_value = 0.0;

    for (unsigned int pore_id = 0; pore_id < pore_centers.size(); ++pore_id)
    {
        const double distance = (point - pore_centers[pore_id]).norm();
        const double radius = pore_radii[pore_id];
        double local_phi = 0.0;

        if (distance < radius - lo)
            local_phi = 1.0;
        else if (distance <= radius)
        {
            const double s = (distance - (radius - lo)) / lo;
            local_phi = 0.5 * (1.0 + std::cos(numbers::PI * s));
        }

        phi_value = std::max(phi_value, local_phi);
    }

    return phi_value;
}

void PhaseField::rebuild_porosity_damage_field()
{
    for (const auto &cell : dof_handler_damage.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        for (const auto vertex_number : cell->vertex_indices())
        {
            const Point<2> point = cell->vertex(vertex_number);
            const double phi_value = compute_diffuse_phi(point);

            const auto dof_index = cell->vertex_dof_index(vertex_number, 0);
            if (locally_owned_dofs_damage.is_element(dof_index))
            {
                const double current = completely_distributed_solution_damage[dof_index];
                completely_distributed_solution_damage[dof_index] =
                    std::max(current, std::min(1.0, phi_value));
            }
        }
    }

    completely_distributed_solution_damage.compress(VectorOperation::insert);
}

void PhaseField::enforce_initial_porosity_damage()
{
    constraints_damage.distribute(completely_distributed_solution_damage);

    locally_relevant_solution_damage = completely_distributed_solution_damage;
    locally_relevant_solution_damage.update_ghost_values();
    completely_distributed_solution_damage_old =
        completely_distributed_solution_damage;
}

void PhaseField::lock_initial_porosity_history()
{
    const double H_lock = Gc / (2.0 * lo);

    FEValues<2> fe_values_damage(fe_damage,
                                 quadrature_formula_damage,
                                 update_values);
    std::vector<double> phi_q(quadrature_formula_damage.size());

    for (const auto &cell : dof_handler_damage.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        fe_values_damage.reinit(cell);
        fe_values_damage.get_function_values(locally_relevant_solution_damage,
                                             phi_q);

        const std::vector<std::shared_ptr<MyQData>> local_qph =
            quadrature_point_history_field.get_data(cell);

        for (unsigned int q = 0; q < phi_q.size(); ++q)
            if (phi_q[q] > 0.5)
            {
                local_qph[q]->value_H = H_lock;
                local_qph[q]->value_H_new = H_lock;
            }
    }
}

void PhaseField::initialize_random_porosity()
{
    // Always start from a zero damage field on both the owned and
    // ghosted vectors.  For pore_count == 0 this is all we need.
    completely_distributed_solution_damage = 0.0;
    completely_distributed_solution_damage_old = 0.0;
    locally_relevant_solution_damage = 0.0;

    if (pore_count == 0)
    {
        pcout << "Random porosity disabled because pore_count = 0." << std::endl;
        return;
    }

    AssertThrow(min_pore_size > 0.0,
        ExcMessage("min_pore_size must be > 0."));
    AssertThrow(max_pore_size >= min_pore_size,
        ExcMessage("max_pore_size must be >= min_pore_size."));
    AssertThrow(crack_spread_radius > 0.0,
        ExcMessage("crack_spread_radius must be > 0."));
    AssertThrow(boundary_margin >= 0.0,
        ExcMessage("boundary_margin must be >= 0."));
    AssertThrow(lo > 0.0, ExcMessage("lo/l0 must be > 0."));

    const double min_pore_radius = 0.5 * min_pore_size;
    const double max_pore_radius = 0.5 * max_pore_size;
    Point<2> distribution_center;
    Bounds global_bounds;
    std::vector<Point<2>> crack_points;

    // ----------------------------------------------------------------
    // 1. Placement: pick pore centres and radii (rank 0, then broadcast)
    // ----------------------------------------------------------------
    if (pore_centers.empty())
    {
        pore_centers.clear();
        pore_radii.clear();

        Point<2> local_min;
        Point<2> local_max;
        for (unsigned int d = 0; d < 2; ++d)
        {
            local_min[d] = std::numeric_limits<double>::max();
            local_max[d] = std::numeric_limits<double>::lowest();
        }

        double local_max_cell_diameter = 0.0;

        for (const auto &cell : dof_handler_damage.active_cell_iterators())
            if (cell->is_locally_owned())
            {
                local_max_cell_diameter =
                    std::max(local_max_cell_diameter, cell->diameter());

                for (const auto vertex_number : cell->vertex_indices())
                {
                    const Point<2> &point = cell->vertex(vertex_number);
                    for (unsigned int d = 0; d < 2; ++d)
                    {
                        local_min[d] = std::min(local_min[d], point[d]);
                        local_max[d] = std::max(local_max[d], point[d]);
                    }
                }
            }

        for (unsigned int d = 0; d < 2; ++d)
        {
            global_bounds.min[d] =
                Utilities::MPI::min(local_min[d], mpi_communicator);
            global_bounds.max[d] =
                Utilities::MPI::max(local_max[d], mpi_communicator);
        }

        const double global_max_cell_diameter =
            Utilities::MPI::max(local_max_cell_diameter, mpi_communicator);

        if (global_max_cell_diameter > lo / 3.0)
            pcout << "WARNING: Mesh does not resolve phase-field length scale. max cell diameter h="
                  << global_max_cell_diameter
                  << " while lo/3=" << lo / 3.0
                  << "."
                  << std::endl;

        if (min_pore_size / lo < 1.0)
            pcout << "WARNING: min pore diameter / l0 = "
                  << (min_pore_size / lo)
                  << ". Pores smaller than l0 are not physically meaningful."
                  << std::endl;

        const double margin = boundary_margin + max_pore_radius;
        AssertThrow((global_bounds.max[0] - global_bounds.min[0]) > 2.0 * margin &&
                    (global_bounds.max[1] - global_bounds.min[1]) > 2.0 * margin,
            ExcMessage("Domain is too small for the requested pore sizes."));

        const Point<2> domain_center(
            0.5 * (global_bounds.min[0] + global_bounds.max[0]),
            0.5 * (global_bounds.min[1] + global_bounds.max[1]));
        const double half_width =
            0.5 * (global_bounds.max[0] - global_bounds.min[0]);
        const double half_height =
            0.5 * (global_bounds.max[1] - global_bounds.min[1]);
        const double max_allowed_radius = std::min(half_width, half_height) - margin;
        const double effective_spread_radius =
            std::min(crack_spread_radius, max_allowed_radius);

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 &&
            crack_spread_radius > max_allowed_radius)
            pcout << "Requested crack_spread_radius=" << crack_spread_radius
                  << " is too large for the mesh. Clamping it to "
                  << effective_spread_radius << std::endl;

        distribution_center = domain_center;
        const double pore_spacing_factor = 1.1;

        std::vector<double> local_crack_buffer;
        for (const auto &cell : triangulation.active_cell_iterators())
            if (cell->is_locally_owned())
                for (const auto face_number : cell->face_indices())
                    if (cell->face(face_number)->at_boundary() &&
                        cell->face(face_number)->boundary_id() == 4)
                        for (const auto vertex_number :
                             cell->face(face_number)->vertex_indices())
                        {
                            const Point<2> point =
                                cell->face(face_number)->vertex(vertex_number);
                            local_crack_buffer.push_back(point[0]);
                            local_crack_buffer.push_back(point[1]);
                        }

        const std::vector<std::vector<double>> gathered_crack_buffers =
            Utilities::MPI::gather(mpi_communicator, local_crack_buffer, 0);

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
            for (const auto &buffer : gathered_crack_buffers)
            {
                const std::vector<Point<2>> points = unpack_points(buffer);
                crack_points.insert(crack_points.end(), points.begin(), points.end());
            }

            if (crack_points.empty())
                pcout << "WARNING: no boundary_id=4 crack points found. "
                      << "Pore placement will skip crack avoidance."
                      << std::endl;

            std::seed_seq seed{
                random_seed_base,
                static_cast<unsigned int>(case_number),
                static_cast<unsigned int>(pore_trial_number),
                static_cast<unsigned int>(pore_count)};
            std::mt19937 rng(seed);

            const double sigma = std::max(effective_spread_radius / 3.0, 1e-12);
            std::normal_distribution<double> dist_x(distribution_center[0],
                                                    sigma);
            std::normal_distribution<double> dist_y(distribution_center[1],
                                                    sigma);
            std::uniform_real_distribution<double> radius_dist(min_pore_radius,
                                                               max_pore_radius);

            const unsigned int max_attempts = std::max(10000u, 2000u * pore_count);
            unsigned int attempts = 0;

            while (pore_centers.size() < pore_count && attempts < max_attempts)
            {
                ++attempts;

                const double candidate_radius = radius_dist(rng);
                const Point<2> candidate(dist_x(rng), dist_y(rng));

                if (candidate[0] < (global_bounds.min[0] + boundary_margin + candidate_radius) ||
                    candidate[0] > (global_bounds.max[0] - boundary_margin - candidate_radius) ||
                    candidate[1] < (global_bounds.min[1] + boundary_margin + candidate_radius) ||
                    candidate[1] > (global_bounds.max[1] - boundary_margin - candidate_radius))
                    continue;

                if ((candidate - distribution_center).norm() > effective_spread_radius)
                    continue;

                bool avoids_crack = true;
                for (const auto &crack_point : crack_points)
                    if ((candidate - crack_point).norm() < 2.0 * candidate_radius)
                    {
                        avoids_crack = false;
                        break;
                    }

                if (!avoids_crack)
                    continue;

                bool separated = true;
                for (unsigned int pore_id = 0; pore_id < pore_centers.size(); ++pore_id)
                    if ((candidate - pore_centers[pore_id]).norm() <
                        (pore_spacing_factor * (candidate_radius + pore_radii[pore_id])))
                    {
                        separated = false;
                        break;
                    }

                if (separated)
                {
                    pore_centers.push_back(candidate);
                    pore_radii.push_back(candidate_radius);
                }
            }

            AssertThrow(pore_centers.size() == pore_count,
                ExcMessage("Could not place all pores. Reduce pore_count, "
                           "reduce max_pore_size, or reduce crack_spread_radius."));
        }

        std::vector<double> pore_buffer(3 * pore_count, 0.0);
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
            for (unsigned int pore_id = 0; pore_id < pore_count; ++pore_id)
            {
                pore_buffer[3 * pore_id] = pore_centers[pore_id][0];
                pore_buffer[3 * pore_id + 1] = pore_centers[pore_id][1];
                pore_buffer[3 * pore_id + 2] = pore_radii[pore_id];
            }

            std::ofstream out(output_directory + "/pores.dat");
            AssertThrow(out, ExcMessage("Could not open pores.dat for writing."));
            for (unsigned int pore_id = 0; pore_id < pore_centers.size(); ++pore_id)
                out << pore_centers[pore_id][0] << ' '
                    << pore_centers[pore_id][1] << ' '
                    << pore_radii[pore_id] << '\n';
        }

        pore_buffer = Utilities::MPI::broadcast(mpi_communicator, pore_buffer, 0);

        if (Utilities::MPI::this_mpi_process(mpi_communicator) != 0)
        {
            pore_centers.resize(pore_count);
            pore_radii.resize(pore_count);

            for (unsigned int pore_id = 0; pore_id < pore_count; ++pore_id)
            {
                pore_centers[pore_id] =
                    Point<2>(pore_buffer[3 * pore_id], pore_buffer[3 * pore_id + 1]);
                pore_radii[pore_id] = pore_buffer[3 * pore_id + 2];
            }
        }
    }
    else if (pore_radii.size() != pore_centers.size())
    {
        pore_radii.assign(pore_centers.size(), max_pore_radius);
    }

    // ----------------------------------------------------------------
    // 2. Write the diffuse pore field directly into the damage solution
    //    (vertex DoFs, take max with the current value).
    // ----------------------------------------------------------------
    rebuild_porosity_damage_field();
    enforce_initial_porosity_damage();
    lock_initial_porosity_history();

    pcout << "Pore initialisation: " << pore_centers.size()
          << " pores, H_lock = " << Gc / (2.0 * lo) << std::endl;
}
