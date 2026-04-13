#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <queue>
#include <fstream>
#include <cmath>
#include <set>
#include <limits>
#include <iostream>

#include "ALClasses.h"
#include "IMASH.H"

// ============================================================
// Existing utility / lookup functions
// ============================================================

int FindSubsetIndexByName(IdsNs::IDS& ids, const std::string& name)
{
    auto& subsets = ids._edge_profiles.grid_ggd(0).grid_subset;

    for (int s = subsets.lbound(0); s <= subsets.ubound(0); ++s) {
        if (subsets(s).identifier.name == name) {
            return s;
        }
    }

    throw std::runtime_error("Subset not found: " + name);
}

// NOTE: template definitions should live in the header.
// If you leave them here, you may get linker errors depending on use sites.

// ============================================================
// Small numerical helpers
// ============================================================

double LinearInterp1D(const std::vector<double>& x,
                      const std::vector<double>& y,
                      double xq)
{
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Invalid interpolation arrays");
    }

    if (xq <= x.front()) return y.front();
    if (xq >= x.back())  return y.back();

    auto it = std::upper_bound(x.begin(), x.end(), xq);
    std::size_t i1 = static_cast<std::size_t>(it - x.begin());
    std::size_t i0 = i1 - 1;

    double w = (xq - x[i0]) / (x[i1] - x[i0]);
    return (1.0 - w) * y[i0] + w * y[i1];
}

double WrapAngleDiff(double a, double b)
{
    double d = a - b;
    while (d >  M_PI) d -= 2.0 * M_PI;
    while (d < -M_PI) d += 2.0 * M_PI;
    return d;
}

double ComputeTheta(double R, double Z, double R0, double Z0)
{
    return std::atan2(Z - Z0, R - R0);
}

int FindBracket(const std::vector<double>& x, double xq)
{
    if (x.empty()) return -1;
    if (xq < x.front() || xq > x.back()) return -1;

    auto it = std::upper_bound(x.begin(), x.end(), xq);
    if (it == x.begin()) return 0;
    if (it == x.end())   return static_cast<int>(x.size()) - 2;
    return static_cast<int>((it - x.begin()) - 1);
}

double BilinearInterpEq(const EqGrid2D& g, double Rq, double Zq)
{
    const int i = FindBracket(g.Rvals, Rq);
    const int j = FindBracket(g.Zvals, Zq);

    if (i < 0 || j < 0 || i >= g.nR - 1 || j >= g.nZ - 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double R0 = g.Rvals[i];
    const double R1 = g.Rvals[i + 1];
    const double Z0 = g.Zvals[j];
    const double Z1 = g.Zvals[j + 1];

    const double tx = (Rq - R0) / (R1 - R0);
    const double tz = (Zq - Z0) / (Z1 - Z0);

    auto idx = [&](int ii, int jj) { return jj * g.nR + ii; };

    const double f00 = g.rho[idx(i,     j)];
    const double f10 = g.rho[idx(i + 1, j)];
    const double f01 = g.rho[idx(i,     j + 1)];
    const double f11 = g.rho[idx(i + 1, j + 1)];

    return (1.0 - tx) * (1.0 - tz) * f00
         + tx         * (1.0 - tz) * f10
         + (1.0 - tx) * tz         * f01
         + tx         * tz         * f11;
}

double LocalInterpRhoTheta(double rho_q,
                           double theta_q,
                           const std::vector<EdgeParamPoint>& src,
                           bool use_ne,
                           double sigma_rho,
                           double sigma_theta,
                           double max_drho,
                           double max_dtheta,
                           double min_weight_sum)
{
    double num = 0.0;
    double den = 0.0;

    for (const auto& s : src) {
        const double drho   = s.rho - rho_q;
        const double dtheta = WrapAngleDiff(s.theta, theta_q);

        if (std::abs(drho)   > max_drho)   continue;
        if (std::abs(dtheta) > max_dtheta) continue;

        const double xr = drho   / sigma_rho;
        const double xt = dtheta / sigma_theta;

        const double d2 = xr * xr + xt * xt;
        const double w  = std::exp(-0.5 * d2);

        num += w * (use_ne ? s.ne : s.Te);
        den += w;
    }

    if (den < min_weight_sum) return 0.0;
    return num / den;
}

// ============================================================
// Structured equilibrium-grid helpers
// ============================================================

EqGrid2D ExtractEqGridRho(IdsNs::IDS& ids, double time)
{
    ids._equilibrium.getSlice(time, CLOSEST_SAMPLE);

    auto& ts0 = ids._equilibrium.time_slice(0);
    auto& p2d = ts0.profiles_2d(0);

    auto& R   = p2d.grid.dim1;
    auto& Z   = p2d.grid.dim2;
    auto& psi = p2d.psi;

    const int iR0 = R.lbound(0);
    const int iR1 = R.ubound(0);
    const int iZ0 = Z.lbound(0);
    const int iZ1 = Z.ubound(0);

    const int nR = iR1 - iR0 + 1;
    const int nZ = iZ1 - iZ0 + 1;

    const double psi_axis = ts0.global_quantities.psi_axis;
    const double psi_bnd  = ts0.global_quantities.psi_boundary;

    const double denom = psi_bnd - psi_axis;
    if (std::abs(denom) == 0.0) {
        throw std::runtime_error("ExtractEqGridRho: invalid psi normalization");
    }

    EqGrid2D g;
    g.nR = nR;
    g.nZ = nZ;
    g.Rvals.resize(static_cast<std::size_t>(nR));
    g.Zvals.resize(static_cast<std::size_t>(nZ));
    g.rho.resize(static_cast<std::size_t>(nR * nZ));

    for (int i = 0; i < nR; ++i) {
        g.Rvals[static_cast<std::size_t>(i)] = R(i + iR0);
    }
    for (int j = 0; j < nZ; ++j) {
        g.Zvals[static_cast<std::size_t>(j)] = Z(j + iZ0);
    }

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const int k = j * nR + i;
            const double psi_ij = psi(i + iR0, j + iZ0);
            const double arg = (psi_ij - psi_axis) / denom;
            g.rho[static_cast<std::size_t>(k)] = (arg >= 0.0) ? std::sqrt(arg) : -1.0;
        }
    }

    return g;
}

// ============================================================
// Edge data extraction
// ============================================================

std::vector<EdgeCellPoint> ExtractEdge2D(IdsNs::IDS& ids, double time)
{
    ids._edge_profiles.getSlice(time, CLOSEST_SAMPLE);

    auto& ggd_data  = ids._edge_profiles.ggd(0);
    auto& electrons = ggd_data.electrons;
    auto& edens     = electrons.density;
    auto& etemp     = electrons.temperature;

    auto& grid0   = ids._edge_profiles.grid_ggd(0);
    auto& spaces  = grid0.space;

    if (spaces.extent(0) == 0) {
        throw std::runtime_error("ExtractEdge2D: no spaces found in grid_ggd");
    }

    auto& space0  = spaces(0);
    auto& objects = space0.objects_per_dimension;
    auto& subsets = grid0.grid_subset;

    if (objects.extent(0) < 3) {
        throw std::runtime_error("ExtractEdge2D: objects_per_dimension does not contain nodes/faces/cells");
    }

    auto& nodes_bucket = objects(0).object;
    auto& cells_bucket = objects(2).object;

    constexpr int cells_identifier_index = 5;

    int cell_subset_slot = FindSubsetSlotByIdentifierIndex(subsets, cells_identifier_index);
    if (cell_subset_slot < 0) {
        throw std::runtime_error("ExtractEdge2D: could not find cells subset");
    }

    auto& cell_subset = subsets(cell_subset_slot);

    int ne_entry = FindFieldEntryForSubset(edens, cells_identifier_index);
    int Te_entry = FindFieldEntryForSubset(etemp, cells_identifier_index);

    if (ne_entry < 0) {
        throw std::runtime_error("ExtractEdge2D: could not find density field for cells subset");
    }
    if (Te_entry < 0) {
        throw std::runtime_error("ExtractEdge2D: could not find temperature field for cells subset");
    }

    const int ncells = cell_subset.element.extent(0);

    if (edens(ne_entry).values.extent(0) != ncells) {
        throw std::runtime_error("ExtractEdge2D: density values size does not match cells subset size");
    }
    if (etemp(Te_entry).values.extent(0) != ncells) {
        throw std::runtime_error("ExtractEdge2D: temperature values size does not match cells subset size");
    }

    std::vector<EdgeCellPoint> pts;
    pts.reserve(static_cast<std::size_t>(ncells));

    for (int k = cell_subset.element.lbound(0); k <= cell_subset.element.ubound(0); ++k) {
        int cell_idx = cell_subset.element(k).object(0).index;

        auto& cell = cells_bucket(cell_idx - 1);
        auto& cell_nodes = cell.nodes;

        if (cell_nodes.extent(0) == 0) {
            throw std::runtime_error("ExtractEdge2D: encountered cell with no nodes");
        }

        std::set<int> unique_nodes;
        for (int j = cell_nodes.lbound(0); j <= cell_nodes.ubound(0); ++j) {
            unique_nodes.insert(cell_nodes(j));
        }

        double Rc = 0.0;
        double Zc = 0.0;

        for (int node_idx : unique_nodes) {
            auto& node = nodes_bucket(node_idx - 1);
            auto& geom = node.geometry;

            if (geom.extent(0) != 2) {
                throw std::runtime_error("ExtractEdge2D: expected node.geometry to have extent 2");
            }

            Rc += geom(0);
            Zc += geom(1);
        }

        Rc /= static_cast<double>(unique_nodes.size());
        Zc /= static_cast<double>(unique_nodes.size());

        pts.push_back({Rc, Zc, edens(ne_entry).values(k), etemp(Te_entry).values(k)});
    }

    return pts;
}

void FillEdgeCellPointsFromImas(IdsNs::IDS& ids,
                                double time,
                                std::vector<EdgeCellPoint>& edge_pts)
{
    edge_pts = ExtractEdge2D(ids, time);
}

std::vector<EdgeParamPoint>
BuildEdgeParamPoints(const std::vector<EdgeCellPoint>& edge_pts,
                     const EqGrid2D& eq,
                     double R0, double Z0)
{
    std::vector<EdgeParamPoint> out;
    out.reserve(edge_pts.size());

    for (const auto& p : edge_pts) {
        const double rho_e = BilinearInterpEq(eq, p.R, p.Z);
        if (!std::isfinite(rho_e)) continue;

        EdgeParamPoint q;
        q.R = p.R;
        q.Z = p.Z;
        q.rho = rho_e;
        q.theta = ComputeTheta(p.R, p.Z, R0, Z0);
        q.ne = p.ne;
        q.Te = p.Te;
        out.push_back(q);
    }

    return out;
}

// ============================================================
// Plasma extraction into reusable structures
// ============================================================

void FillCoreProfilesFromImas(IdsNs::IDS& ids,
                              double time,
                              PlasmaData& plasma)
{
    // 1. Core profiles
    ids._core_profiles.getSlice(time, CLOSEST_SAMPLE);

    auto& profiles = ids._core_profiles.profiles_1d;
    if (profiles.extent(0) == 0) {
        throw std::runtime_error("No core_profiles found");
    }

    auto& prof = profiles(0);

    auto& rho_arr = prof.grid.rho_tor_norm;
    auto& ne_arr  = prof.electrons.density;
    auto& Te_arr  = prof.electrons.temperature;

    std::vector<double> rho_1d, ne_1d, Te_1d;
    for (int i = rho_arr.lbound(0); i <= rho_arr.ubound(0); ++i) {
        rho_1d.push_back(rho_arr(i));
        ne_1d.push_back(ne_arr(i));
        Te_1d.push_back(Te_arr(i));
    }

    std::vector<std::size_t> perm(rho_1d.size());
    for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;

    std::sort(perm.begin(), perm.end(),
              [&](std::size_t a, std::size_t b) {
                  return rho_1d[a] < rho_1d[b];
              });

    std::vector<double> rho1d, ne1d, Te1d;
    rho1d.reserve(perm.size());
    ne1d.reserve(perm.size());
    Te1d.reserve(perm.size());

    for (auto k : perm) {
        rho1d.push_back(rho_1d[k]);
        ne1d.push_back(ne_1d[k]);
        Te1d.push_back(Te_1d[k]);
    }

    // 2. Equilibrium
    ids._equilibrium.getSlice(time, CLOSEST_SAMPLE);

    auto& ts0 = ids._equilibrium.time_slice(0);
    auto& p2d = ts0.profiles_2d(0);

    auto& R   = p2d.grid.dim1;
    auto& Z   = p2d.grid.dim2;
    auto& psi = p2d.psi;

    const int iR0 = R.lbound(0);
    const int iR1 = R.ubound(0);
    const int iZ0 = Z.lbound(0);
    const int iZ1 = Z.ubound(0);

    const int nR = iR1 - iR0 + 1;
    const int nZ = iZ1 - iZ0 + 1;

    const double psi_axis = ts0.global_quantities.psi_axis;
    const double psi_bnd  = ts0.global_quantities.psi_boundary;

    const double denom = psi_bnd - psi_axis;
    if (std::abs(denom) == 0.0) {
        throw std::runtime_error("psi normalization invalid");
    }

    plasma.resize(nR, nZ);

    for (int i = 0; i < nR; ++i) {
        plasma.grid.Rvals[static_cast<std::size_t>(i)] = R(i + iR0);
    }
    for (int j = 0; j < nZ; ++j) {
        plasma.grid.Zvals[static_cast<std::size_t>(j)] = Z(j + iZ0);
    }

    // Local flood-fill scratch arrays preserve your original indexing:
    // idx(i,j) = i*nZ + j
    std::vector<double> rho_map(static_cast<std::size_t>(nR * nZ), -1.0);
    std::vector<unsigned char> inside_candidate(static_cast<std::size_t>(nR * nZ), 0);
    std::vector<unsigned char> boundary_connected(static_cast<std::size_t>(nR * nZ), 0);

    auto idx_local = [nZ](int ii, int jj) {
        return ii * nZ + jj;
    };

    const double rho_tol = 1.0e-10;

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const double psi_ij = psi(i + iR0, j + iZ0);
            const double arg = (psi_ij - psi_axis) / denom;

            double rho_ij = -1.0;
            if (arg >= 0.0) {
                rho_ij = std::sqrt(arg);
            }

            rho_map[static_cast<std::size_t>(idx_local(i, j))] = rho_ij;

            if (rho_ij >= 0.0 && rho_ij <= 1.0 + rho_tol) {
                inside_candidate[static_cast<std::size_t>(idx_local(i, j))] = 1;
            }
        }
    }

    std::queue<std::pair<int,int>> q;

    auto try_push = [&](int i, int j)
    {
        if (i < 0 || i >= nR || j < 0 || j >= nZ) return;
        const int k = idx_local(i, j);
        if (!inside_candidate[static_cast<std::size_t>(k)]) return;
        if (boundary_connected[static_cast<std::size_t>(k)]) return;
        boundary_connected[static_cast<std::size_t>(k)] = 1;
        q.push({i, j});
    };

    for (int i = 0; i < nR; ++i) {
        try_push(i, 0);
        try_push(i, nZ - 1);
    }
    for (int j = 0; j < nZ; ++j) {
        try_push(0, j);
        try_push(nR - 1, j);
    }

    while (!q.empty()) {
        auto ij = q.front();
        q.pop();

        const int i = ij.first;
        const int j = ij.second;

        try_push(i - 1, j);
        try_push(i + 1, j);
        try_push(i, j - 1);
        try_push(i, j + 1);
    }

    // Store into reusable structured fields using canonical field indexing
    // handled by PlasmaData / ScalarField2D
    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const double rho_ij = rho_map[static_cast<std::size_t>(idx_local(i, j))];

            plasma.rho(i, j) = rho_ij;

            const bool is_candidate = inside_candidate[static_cast<std::size_t>(idx_local(i, j))] != 0;
            const bool is_boundary_connected = boundary_connected[static_cast<std::size_t>(idx_local(i, j))] != 0;

            double ne_ij = 0.0;
            double Te_ij = 0.0;
            double m_ij  = 0.0;

            if (is_candidate && !is_boundary_connected && rho_ij <= 1.0) {
                ne_ij = LinearInterp1D(rho1d, ne1d, rho_ij);
                Te_ij = LinearInterp1D(rho1d, Te1d, rho_ij);
                m_ij  = 1.0;
            }

            plasma.ne_core(i, j)   = ne_ij;
            plasma.Te_core(i, j)   = Te_ij;
            plasma.core_mask(i, j) = m_ij;
        }
    }
}

void FillBfieldFromImas(IdsNs::IDS& ids,
                        double time,
                        PlasmaData& plasma)
{
    ids._equilibrium.getSlice(time, PREVIOUS_SAMPLE);
    auto& ts0 = ids._equilibrium.time_slice(0);

    const double BT = ids._equilibrium.vacuum_toroidal_field.b0(0);
    const double R0 = ids._equilibrium.vacuum_toroidal_field.r0;

    auto& p2d = ts0.profiles_2d(0);

    auto& Rvec = p2d.grid.dim1;
    auto& Zvec = p2d.grid.dim2;
    auto& BR2D = p2d.b_field_r;
    auto& BZ2D = p2d.b_field_z;

    const int nR = static_cast<int>(Rvec.size());
    const int nZ = static_cast<int>(Zvec.size());

    if (nR <= 0 || nZ <= 0) {
        throw std::runtime_error("FillBfieldFromImas: empty equilibrium grid.");
    }

    if (static_cast<int>(BR2D.size()) != nR * nZ ||
        static_cast<int>(BZ2D.size()) != nR * nZ) {
        throw std::runtime_error("FillBfieldFromImas: inconsistent magnetic field array sizes.");
    }

    if (plasma.empty()) {
        plasma.resize(nR, nZ);
        for (int i = 0; i < nR; ++i) {
            plasma.grid.Rvals[static_cast<std::size_t>(i)] = Rvec(i);
        }
        for (int j = 0; j < nZ; ++j) {
            plasma.grid.Zvals[static_cast<std::size_t>(j)] = Zvec(j);
        }
    } else if (plasma.grid.nR != nR || plasma.grid.nZ != nZ) {
        throw std::runtime_error("FillBfieldFromImas: plasma grid size does not match equilibrium grid.");
    }

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const int k = j * nR + i;
            const double R = Rvec(i);

            plasma.B.br(i, j)   = BR2D(k);
            plasma.B.bz(i, j)   = BZ2D(k);
            plasma.B.bphi(i, j) = BT * R0 / R;
        }
    }
}

void FillEdgeProfilesOnGridFromImas(IdsNs::IDS& ids,
                                    double time,
                                    PlasmaData& plasma,
                                    double sigma_rho,
                                    double sigma_theta)
{
    std::vector<EdgeCellPoint> edge_pts = ExtractEdge2D(ids, time);
    if (edge_pts.empty()) {
        throw std::runtime_error("FillEdgeProfilesOnGridFromImas: no edge points extracted");
    }

    EqGrid2D eq = ExtractEqGridRho(ids, time);

    ids._equilibrium.getSlice(time, CLOSEST_SAMPLE);
    auto& ts0 = ids._equilibrium.time_slice(0);

    const double R_axis = ts0.global_quantities.magnetic_axis.r;
    const double Z_axis = ts0.global_quantities.magnetic_axis.z;

    std::vector<EdgeParamPoint> src = BuildEdgeParamPoints(edge_pts, eq, R_axis, Z_axis);
    if (src.empty()) {
        throw std::runtime_error("FillEdgeProfilesOnGridFromImas: no valid edge param points");
    }

    const double max_drho   = 3.0 * sigma_rho;
    const double max_dtheta = 3.0 * sigma_theta;

    if (plasma.empty()) {
        plasma.resize(eq.nR, eq.nZ);
        plasma.grid.Rvals = eq.Rvals;
        plasma.grid.Zvals = eq.Zvals;
        for (int i = 0; i < eq.nR; ++i) {
            for (int j = 0; j < eq.nZ; ++j) {
                plasma.rho(i, j) = eq.rho[static_cast<std::size_t>(j * eq.nR + i)];
            }
        }
    } else if (plasma.grid.nR != eq.nR || plasma.grid.nZ != eq.nZ) {
        throw std::runtime_error("FillEdgeProfilesOnGridFromImas: plasma grid size mismatch");
    }

    for (int i = 0; i < eq.nR; ++i) {
        for (int j = 0; j < eq.nZ; ++j) {
            const int k = j * eq.nR + i;

            const double Rq   = eq.Rvals[static_cast<std::size_t>(i)];
            const double Zq   = eq.Zvals[static_cast<std::size_t>(j)];
            const double rhoq = eq.rho[static_cast<std::size_t>(k)];
            const double thq  = ComputeTheta(Rq, Zq, R_axis, Z_axis);

            double ne_i = 0.0;
            double Te_i = 0.0;
            double m_i  = 0.0;

            if (rhoq > 1.0) {
                ne_i = LocalInterpRhoTheta(rhoq, thq, src, true,
                                           sigma_rho, sigma_theta,
                                           max_drho, max_dtheta);

                Te_i = LocalInterpRhoTheta(rhoq, thq, src, false,
                                           sigma_rho, sigma_theta,
                                           max_drho, max_dtheta);

                if (ne_i < 0.0) ne_i = 0.0;
                if (Te_i < 0.0) Te_i = 0.0;
                m_i = 1.0;
            }

            plasma.ne_edge(i, j)   = ne_i;
            plasma.Te_edge(i, j)   = Te_i;
            plasma.edge_mask(i, j) = m_i;
        }
    }
}

void MergePlasmaProfiles(PlasmaData& plasma)
{
    if (plasma.empty()) return;

    const int nR = plasma.grid.nR;
    const int nZ = plasma.grid.nZ;

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const double core_m = plasma.core_mask(i, j);
            const double edge_m = plasma.edge_mask(i, j);

            if (core_m > 0.5) {
                plasma.ne(i, j) = plasma.ne_core(i, j);
                plasma.Te(i, j) = plasma.Te_core(i, j);
                plasma.source_mask(i, j) = 1.0;
            } else if (edge_m > 0.5) {
                plasma.ne(i, j) = plasma.ne_edge(i, j);
                plasma.Te(i, j) = plasma.Te_edge(i, j);
                plasma.source_mask(i, j) = 2.0;
            } else {
                plasma.ne(i, j) = 0.0;
                plasma.Te(i, j) = 0.0;
                plasma.source_mask(i, j) = 0.0;
            }
        }
    }
}

PlasmaData ExtractPlasmaData(IdsNs::IDS& ids,
                             double time,
                             bool fill_core,
                             bool fill_bfield,
                             bool fill_edge_on_grid,
                             bool merge_density,
                             bool merge_temperature,
                             double sigma_rho,
                             double sigma_theta)
{
    PlasmaData plasma;

    if (fill_core) {
        FillCoreProfilesFromImas(ids, time, plasma);
    }
    if (fill_bfield) {
        FillBfieldFromImas(ids, time, plasma);
    }
    if (fill_edge_on_grid) {
        FillEdgeProfilesOnGridFromImas(ids, time, plasma, sigma_rho, sigma_theta);
    }
    if (merge_density || merge_temperature) {
        MergePlasmaProfiles(plasma);
    }

    return plasma;
}

// ============================================================
// Wall extraction
// ============================================================

WallData ExtractWall2D(IdsNs::IDS& ids, double time)
{
    ids._wall.getSlice(time, PREVIOUS_SAMPLE);

    WallData wall;

    const int ndesc = static_cast<int>(ids._wall.description_2d.size());

    for (int i1 = 0; i1 < ndesc; ++i1) {
        auto& desc = ids._wall.description_2d(i1);

        // Limiter
        const int nlim = static_cast<int>(desc.limiter.unit.size());
        for (int i2 = 0; i2 < nlim; ++i2) {
            auto& unit = desc.limiter.unit(i2);
            auto& r = unit.outline.r;
            auto& z = unit.outline.z;

            const int n = static_cast<int>(r.size());
            if (n == 0) continue;
            if (static_cast<int>(z.size()) != n) {
                throw std::runtime_error("ExtractWall2D: limiter outline r/z size mismatch.");
            }

            PolylineRZ c("limiter");
            c.closed = false;
            for (int k = 0; k < n; ++k) {
                c.addPoint(r(k), z(k));
            }
            wall.addContour(c);
        }

        // Vessel annular inner/outer
        const int nves = static_cast<int>(desc.vessel.unit.size());
        for (int i2 = 0; i2 < nves; ++i2) {
            auto& unit = desc.vessel.unit(i2);

            auto& rin = unit.annular.outline_inner.r;
            auto& zin = unit.annular.outline_inner.z;
            const int nin = static_cast<int>(rin.size());
            if (nin > 0) {
                if (static_cast<int>(zin.size()) != nin) {
                    throw std::runtime_error("ExtractWall2D: vessel inner outline r/z size mismatch.");
                }
                PolylineRZ c("vessel_inner");
                c.closed = false;
                for (int k = 0; k < nin; ++k) {
                    c.addPoint(rin(k), zin(k));
                }
                wall.addContour(c);
            }

            auto& rout = unit.annular.outline_outer.r;
            auto& zout = unit.annular.outline_outer.z;
            const int nout = static_cast<int>(rout.size());
            if (nout > 0) {
                if (static_cast<int>(zout.size()) != nout) {
                    throw std::runtime_error("ExtractWall2D: vessel outer outline r/z size mismatch.");
                }
                PolylineRZ c("vessel_outer");
                c.closed = false;
                for (int k = 0; k < nout; ++k) {
                    c.addPoint(rout(k), zout(k));
                }
                wall.addContour(c);
            }
        }
    }

    return wall;
}

// ============================================================
// Magnetic-axis helper
// ============================================================

MagneticAxisPoint ExtractMagneticAxis(IdsNs::IDS& ids, double time)
{
    ids._equilibrium.getSlice(time, PREVIOUS_SAMPLE);
    auto& ts0 = ids._equilibrium.time_slice(0);

    MagneticAxisPoint axis;
    axis.R = ts0.global_quantities.magnetic_axis.r;
    axis.Z = ts0.global_quantities.magnetic_axis.z;
    return axis;
}

// ============================================================
// Structure-based debug / Python-friendly writers
// ============================================================

void WriteCore2D(const PlasmaData& plasma,
                 const std::string& filename)
{
    if (plasma.empty()) {
        throw std::runtime_error("WriteCore2D(plasma): plasma is empty");
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("WriteCore2D(plasma): Cannot open output file");
    }

    out << std::setprecision(12);
    out << "# R Z rho ne Te\n";

    const int nR = plasma.grid.nR;
    const int nZ = plasma.grid.nZ;

    for (int i = 0; i < nR; ++i) {
        const double Rv = plasma.grid.Rvals[static_cast<std::size_t>(i)];
        for (int j = 0; j < nZ; ++j) {
            const double Zv = plasma.grid.Zvals[static_cast<std::size_t>(j)];
            out << Rv                  << " "
                << Zv                  << " "
                << plasma.rho(i, j)    << " "
                << plasma.ne_core(i,j) << " "
                << plasma.Te_core(i,j) << "\n";
        }
        out << "\n";
    }
}

void WriteBfield2D(const PlasmaData& plasma,
                   const std::string& filename)
{
    if (plasma.empty()) {
        throw std::runtime_error("WriteBfield2D(plasma): plasma is empty");
    }

    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("WriteBfield2D(plasma): could not open output file " + filename);
    }

    os << "# R Z BR BZ Bphi Bpol Bmag\n";

    const int nR = plasma.grid.nR;
    const int nZ = plasma.grid.nZ;

    for (int i = 0; i < nR; ++i) {
        const double R = plasma.grid.Rvals[static_cast<std::size_t>(i)];
        for (int j = 0; j < nZ; ++j) {
            const double Z    = plasma.grid.Zvals[static_cast<std::size_t>(j)];
            const double BR   = plasma.B.br(i, j);
            const double BZ   = plasma.B.bz(i, j);
            const double Bphi = plasma.B.bphi(i, j);
            const double Bpol = std::sqrt(BR*BR + BZ*BZ);
            const double Bmag = std::sqrt(BR*BR + BZ*BZ + Bphi*Bphi);

            os << R    << " "
               << Z    << " "
               << BR   << " "
               << BZ   << " "
               << Bphi << " "
               << Bpol << " "
               << Bmag << "\n";
        }
    }
}

void WriteWall2D(const WallData& wall,
                 const std::string& filename)
{
    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("WriteWall2D(wall): could not open output file " + filename);
    }

    os << "# type segment R Z\n";

    int segment = 0;
    for (const auto& c : wall.rz_contours) {
        int type = 0;
        if (c.name == "limiter") type = 0;
        else if (c.name == "vessel_inner") type = 1;
        else if (c.name == "vessel_outer") type = 2;

        for (int k = 0; k < c.size(); ++k) {
            os << type << " " << segment << " "
               << c.R[static_cast<std::size_t>(k)] << " "
               << c.Z[static_cast<std::size_t>(k)] << "\n";
        }
        os << "\n";
        ++segment;
    }
}

void WriteMagneticAxis(const MagneticAxisPoint& axis,
                       const std::string& filename)
{
    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("WriteMagneticAxis(axis): could not open output file " + filename);
    }

    os << std::setprecision(12);
    os << axis.R << " " << axis.Z << "\n";
}

void WriteEdge2D(const std::vector<EdgeCellPoint>& edge_pts,
                 const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("WriteEdge2D(edge_pts): Could not open output file: " + filename);
    }

    out << std::setprecision(12);
    out << "# k Rc Zc ne Te\n";

    for (std::size_t k = 0; k < edge_pts.size(); ++k) {
        const auto& p = edge_pts[k];
        out << k    << " "
            << p.R  << " "
            << p.Z  << " "
            << p.ne << " "
            << p.Te << "\n";
    }
}

void WriteEdgeOnGridRhoTheta(const PlasmaData& plasma,
                             const std::string& filename)
{
    if (plasma.empty()) {
        throw std::runtime_error("WriteEdgeOnGridRhoTheta(plasma): plasma is empty");
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("WriteEdgeOnGridRhoTheta(plasma): could not open output file " + filename);
    }

    out << std::setprecision(12);
    out << "# R Z ne_edge_interp Te_edge_interp rho theta\n";

    // theta is not stored in PlasmaData, so write a placeholder NaN here
    const double nanv = std::numeric_limits<double>::quiet_NaN();

    for (int i = 0; i < plasma.grid.nR; ++i) {
        const double Rq = plasma.grid.Rvals[static_cast<std::size_t>(i)];
        for (int j = 0; j < plasma.grid.nZ; ++j) {
            const double Zq = plasma.grid.Zvals[static_cast<std::size_t>(j)];
            out << Rq                  << " "
                << Zq                  << " "
                << plasma.ne_edge(i,j) << " "
                << plasma.Te_edge(i,j) << " "
                << plasma.rho(i,j)     << " "
                << nanv                << "\n";
        }
    }
}

// ============================================================
// Backward-compatible wrappers
// ============================================================

void WriteCore2D(IdsNs::IDS& ids,
                 double time,
                 const std::string& filename)
{
    PlasmaData plasma;
    FillCoreProfilesFromImas(ids, time, plasma);
    WriteCore2D(plasma, filename);
}

void WriteEdge2D(IdsNs::IDS& ids,
                 double time,
                 const std::string& filename)
{
    std::vector<EdgeCellPoint> edge_pts = ExtractEdge2D(ids, time);
    WriteEdge2D(edge_pts, filename);
}

void WriteBfield2D(IdsNs::IDS& ids,
                   double time,
                   std::string filename)
{
    PlasmaData plasma;
    FillBfieldFromImas(ids, time, plasma);
    WriteBfield2D(plasma, filename);
}

void WriteWall2D(IdsNs::IDS& ids,
                 double time,
                 std::string filename)
{
    WallData wall = ExtractWall2D(ids, time);
    WriteWall2D(wall, filename);
}

void WriteMagneticAxis(IdsNs::IDS& ids,
                       double time,
                       const std::string& filename)
{
    MagneticAxisPoint axis = ExtractMagneticAxis(ids, time);
    WriteMagneticAxis(axis, filename);
}

void WriteEdgeOnGridRhoTheta(IdsNs::IDS& ids,
                             double time,
                             const std::string& filename,
                             double sigma_rho,
                             double sigma_theta)
{
    PlasmaData plasma;
    FillCoreProfilesFromImas(ids, time, plasma);
    FillEdgeProfilesOnGridFromImas(ids, time, plasma, sigma_rho, sigma_theta);

    // For compatibility with your old output, compute theta here directly.
    ids._equilibrium.getSlice(time, CLOSEST_SAMPLE);
    auto& ts0 = ids._equilibrium.time_slice(0);
    const double R_axis = ts0.global_quantities.magnetic_axis.r;
    const double Z_axis = ts0.global_quantities.magnetic_axis.z;

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("WriteEdgeOnGridRhoTheta: could not open output file " + filename);
    }

    out << std::setprecision(12);
    out << "# R Z ne_edge_interp Te_edge_interp rho theta\n";

    for (int i = 0; i < plasma.grid.nR; ++i) {
        const double Rq = plasma.grid.Rvals[static_cast<std::size_t>(i)];
        for (int j = 0; j < plasma.grid.nZ; ++j) {
            const double Zq = plasma.grid.Zvals[static_cast<std::size_t>(j)];
            const double rhoq = plasma.rho(i,j);
            const double thq = ComputeTheta(Rq, Zq, R_axis, Z_axis);

            out << Rq                  << " "
                << Zq                  << " "
                << plasma.ne_edge(i,j) << " "
                << plasma.Te_edge(i,j) << " "
                << rhoq                << " "
                << thq                 << "\n";
        }
    }
}