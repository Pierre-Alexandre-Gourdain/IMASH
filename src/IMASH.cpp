#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <queue>
#include <fstream>
#include <cmath>

#include "IMASH.H"

EdgeCut ExtractEdgeCut(IdsNs::IDS& ids, int subset_index)
{
    auto& g      = ids._edge_profiles.ggd_fast(0);
    auto& e      = g.electrons;
    auto& edens  = e.density;
    auto& etemp  = e.temperature;

    auto& grid0   = ids._edge_profiles.grid_ggd(0);
    auto& spaces  = grid0.space;

    if (spaces.extent(0) <= 0) {
        throw std::runtime_error("grid_ggd(0).space is empty");
    }

    auto& space0  = spaces(0);
    auto& objects = space0.objects_per_dimension;

    if (objects.extent(0) < 2) {
        throw std::runtime_error("objects_per_dimension does not contain nodes/faces");
    }

    auto& nodes_bucket = objects(0).object;   // nodes
    auto& faces_bucket = objects(1).object;   // faces

    auto& subsets = grid0.grid_subset;

    if (subset_index < subsets.lbound(0) || subset_index > subsets.ubound(0)) {
        throw std::runtime_error("subset_index out of range");
    }

    auto& subset = subsets(subset_index);

    EdgeCut cut;
    cut.name = subset.identifier.name;

    const int npts = subset.element.extent(0);

    if (edens.extent(0) != npts || etemp.extent(0) != npts) {
        throw std::runtime_error(
            "Field size does not match subset size for subset '" + cut.name + "'"
        );
    }

    cut.points.reserve(npts);

    for (int i = subset.element.lbound(0); i <= subset.element.ubound(0); ++i) {
        auto& el = subset.element(i);

        if (el.object.extent(0) != 1) {
            throw std::runtime_error("Expected one object per subset element");
        }

        int dim      = el.object(0).dimension;
        int face_idx = el.object(0).index;

        if (dim != 2) {
            throw std::runtime_error(
                "Subset '" + cut.name + "' is not face-based (dimension != 2)"
            );
        }

        // IMAS indices appear to be 1-based
        auto& face = faces_bucket(face_idx - 1);

        if (face.nodes.extent(0) != 2) {
            throw std::runtime_error("Expected face to have exactly 2 nodes");
        }

        int n0 = face.nodes(0);
        int n1 = face.nodes(1);

        auto& node0 = nodes_bucket(n0 - 1);
        auto& node1 = nodes_bucket(n1 - 1);

        auto& g0 = node0.geometry;
        auto& g1 = node1.geometry;

        if (g0.extent(0) != 2 || g1.extent(0) != 2) {
            throw std::runtime_error("Expected node.geometry to have extent 2");
        }

        double Rf = 0.5 * (g0(0) + g1(0));
        double Zf = 0.5 * (g0(1) + g1(1));

        EdgeCutPoint p;
        p.R  = Rf;
        p.Z  = Zf;
        p.ne = edens(i).value;
        p.Te = etemp(i).value;

        cut.points.push_back(p);
    }

    return cut;
}

void PrintEdgeCut(const EdgeCut& cut)
{
    std::cout << std::setprecision(10);
    std::cout << "# Edge cut: " << cut.name << "\n";
    std::cout << "# i  R  Z  ne  Te\n";

    for (std::size_t i = 0; i < cut.points.size(); ++i) {
        const auto& p = cut.points[i];
        std::cout << i << "  "
                  << p.R << "  "
                  << p.Z << "  "
                  << p.ne << "  "
                  << p.Te << "\n";
    }
}

void WriteEdgeCut(const EdgeCut& cut, const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    out << std::setprecision(10);
    out << "# Edge cut: " << cut.name << "\n";
    out << "# i R Z ne Te\n";

    for (std::size_t i = 0; i < cut.points.size(); ++i) {
        const auto& p = cut.points[i];
        out << i << " "
            << p.R << " "
            << p.Z << " "
            << p.ne << " "
            << p.Te << "\n";
    }
}

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

template <typename SubsetArray>
int FindSubsetSlotByIdentifierIndex(const SubsetArray& subsets, int identifier_index)
{
    for (int s = subsets.lbound(0); s <= subsets.ubound(0); ++s) {
        if (subsets(s).identifier.index == identifier_index) {
            return s;
        }
    }
    return -1;
}

template <typename FieldArray>
int FindFieldEntryForSubset(const FieldArray& field, int subset_identifier_index)
{
    for (int i = field.lbound(0); i <= field.ubound(0); ++i) {
        if (field(i).grid_subset_index == subset_identifier_index) {
            return i;
        }
    }
    return -1;
}

void WriteEdge2D(IdsNs::IDS& ids, double time, const std::string& filename)
{
	ids._edge_profiles.getSlice(time, CLOSEST_SAMPLE);
	
    auto& ggd_data  = ids._edge_profiles.ggd(0);
    auto& electrons = ggd_data.electrons;
    auto& edens     = electrons.density;
    auto& etemp     = electrons.temperature;

    auto& grid0   = ids._edge_profiles.grid_ggd(0);
    auto& spaces  = grid0.space;

    if (spaces.extent(0) == 0) {
        throw std::runtime_error("No spaces found in grid_ggd");
    }

    auto& space0  = spaces(0);
    auto& objects = space0.objects_per_dimension;
    auto& subsets = grid0.grid_subset;

    if (objects.extent(0) < 3) {
        throw std::runtime_error("objects_per_dimension does not contain nodes/faces/cells");
    }

    auto& nodes_bucket = objects(0).object;   // nodes
    auto& cells_bucket = objects(2).object;   // cells

    constexpr int cells_identifier_index = 5; // from your dataset: "cells"

    int cell_subset_slot = FindSubsetSlotByIdentifierIndex(subsets, cells_identifier_index);
    if (cell_subset_slot < 0) {
        throw std::runtime_error("Could not find cells subset");
    }

    auto& cell_subset = subsets(cell_subset_slot);

    int ne_entry = FindFieldEntryForSubset(edens, cells_identifier_index);
    int Te_entry = FindFieldEntryForSubset(etemp, cells_identifier_index);

    if (ne_entry < 0) {
        throw std::runtime_error("Could not find density field for cells subset");
    }
    if (Te_entry < 0) {
        throw std::runtime_error("Could not find temperature field for cells subset");
    }

    const int ncells = cell_subset.element.extent(0);

    if (edens(ne_entry).values.extent(0) != ncells) {
        throw std::runtime_error("Density values size does not match cells subset size");
    }
    if (etemp(Te_entry).values.extent(0) != ncells) {
        throw std::runtime_error("Temperature values size does not match cells subset size");
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Could not open output file: " + filename);
    }

    out << std::setprecision(12);
    out << "# k Rc Zc ne Te\n";

    for (int k = cell_subset.element.lbound(0); k <= cell_subset.element.ubound(0); ++k) {
        int cell_idx = cell_subset.element(k).object(0).index;

        auto& cell = cells_bucket(cell_idx - 1);
        auto& cell_nodes = cell.nodes;

        if (cell_nodes.extent(0) == 0) {
            throw std::runtime_error("Encountered cell with no nodes");
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
                throw std::runtime_error("Expected node.geometry to have extent 2");
            }

            Rc += geom(0);
            Zc += geom(1);
        }

        Rc /= static_cast<double>(unique_nodes.size());
        Zc /= static_cast<double>(unique_nodes.size());

        double ne = edens(ne_entry).values(k);
        double Te = etemp(Te_entry).values(k);

        out << k << " "
            << Rc << " "
            << Zc << " "
            << ne << " "
            << Te << "\n";
    }
}

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
    std::size_t i1 = it - x.begin();
    std::size_t i0 = i1 - 1;

    double w = (xq - x[i0]) / (x[i1] - x[i0]);
    return (1.0 - w) * y[i0] + w * y[i1];
}

void WriteCore2D(IdsNs::IDS& ids,
                 double time,
                 const std::string& filename)
{
    // ============================================================
    // 1. Core profiles
    // ============================================================
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

    // Sort by rho
    std::vector<std::size_t> perm(rho_1d.size());
    for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;

    std::sort(perm.begin(), perm.end(),
              [&](std::size_t a, std::size_t b) {
                  return rho_1d[a] < rho_1d[b];
              });

    std::vector<double> rho, ne, Te;
    for (auto k : perm) {
        rho.push_back(rho_1d[k]);
        ne.push_back(ne_1d[k]);
        Te.push_back(Te_1d[k]);
    }

    // ============================================================
    // 2. Equilibrium
    // ============================================================
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

    // ============================================================
    // 3. Build rho map and candidate mask
    // ============================================================
    std::vector<double> rho_map(nR * nZ, -1.0);
    std::vector<unsigned char> inside_candidate(nR * nZ, 0);
    std::vector<unsigned char> boundary_connected(nR * nZ, 0);

    auto idx = [nZ](int ii, int jj) {
        return ii * nZ + jj;
    };

    // Small tolerance near rho=1
    const double rho_tol = 1.0e-10;

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const double psi_ij = psi(i + iR0, j + iZ0);
            const double arg = (psi_ij - psi_axis) / denom;

            double rho_ij = -1.0;
            if (arg >= 0.0) {
                rho_ij = std::sqrt(arg);
            }

            rho_map[idx(i, j)] = rho_ij;

            if (rho_ij >= 0.0 && rho_ij <= 1.0 + rho_tol) {
                inside_candidate[idx(i, j)] = 1;
            }
        }
    }

    // ============================================================
    // 4. Flood fill from boundary through candidate region
    //    -> remove all rho<=1 regions touching the box boundary
    // ============================================================
    std::queue<std::pair<int,int>> q;

    auto try_push = [&](int i, int j)
    {
        if (i < 0 || i >= nR || j < 0 || j >= nZ) return;
        const int k = idx(i, j);
        if (!inside_candidate[k]) return;
        if (boundary_connected[k]) return;
        boundary_connected[k] = 1;
        q.push({i, j});
    };

    // Push all boundary cells
    for (int i = 0; i < nR; ++i) {
        try_push(i, 0);
        try_push(i, nZ - 1);
    }
    for (int j = 0; j < nZ; ++j) {
        try_push(0, j);
        try_push(nR - 1, j);
    }

    // 4-neighbor flood fill
    while (!q.empty()) {
        auto [i, j] = q.front();
        q.pop();

        try_push(i - 1, j);
        try_push(i + 1, j);
        try_push(i, j - 1);
        try_push(i, j + 1);
    }

    // ============================================================
    // 5. Output
    // ============================================================
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open output file");
    }

    out << std::setprecision(12);
    out << "# R Z rho ne Te\n";

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
            const double Rv = R(i + iR0);
            const double Zv = Z(j + iZ0);
            const double rho_ij = rho_map[idx(i, j)];

            double ne_ij = 0.0;
            double Te_ij = 0.0;

            const bool is_candidate = inside_candidate[idx(i, j)] != 0;
            const bool is_boundary_connected = boundary_connected[idx(i, j)] != 0;

            // Keep only interior closed rho<=1 regions
            if (is_candidate && !is_boundary_connected && rho_ij <= 1.0) {
                ne_ij = LinearInterp1D(rho, ne, rho_ij);
                Te_ij = LinearInterp1D(rho, Te, rho_ij);
            }

            out << Rv     << " "
                << Zv     << " "
                << rho_ij << " "
                << ne_ij  << " "
                << Te_ij  << "\n";
        }
        out << "\n";
    }
}

void WriteBfield2D(IdsNs::IDS& ids,
                   double time,
                   std::string filename)
{
    ids._equilibrium.getSlice(time, PREVIOUS_SAMPLE);
    auto& ts0 = ids._equilibrium.time_slice(0);
	const double BT = ids._equilibrium.vacuum_toroidal_field.b0(0);
	const double R0 = ids._equilibrium.vacuum_toroidal_field.r0;
    auto& p2d = ts0.profiles_2d(0);

    auto& Rvec = p2d.grid.dim1;
    auto& Zvec = p2d.grid.dim2;

    // Replace these names if your generated IMAS-Cpp bindings differ
    auto& BR2D   = p2d.b_field_r;
    auto& BZ2D   = p2d.b_field_z;

    const int nR = static_cast<int>(Rvec.size());
    const int nZ = static_cast<int>(Zvec.size());

    if (nR <= 0 || nZ <= 0) {
        throw std::runtime_error("WriteBfield2D: empty equilibrium grid.");
    }
auto& psi2D = p2d.psi;

std::cout << "nR = " << Rvec.size() << "\n";
std::cout << "nZ = " << Zvec.size() << "\n";
std::cout << "nR*nZ = " << Rvec.size() * Zvec.size() << "\n";
std::cout << "psi size   = " << psi2D.size() << "\n";
std::cout << "BR size    = " << BR2D.size() << "\n";
std::cout << "BZ size    = " << BZ2D.size() << "\n";
std::cout << "BT size    = " << BT << "\n";

    if (static_cast<int>(BR2D.size())   != nR * nZ ||
	static_cast<int>(BZ2D.size())   != nR * nZ){
        throw std::runtime_error("WriteBfield2D: inconsistent magnetic field array sizes.");
    }

    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("WriteBfield2D: could not open output file " + filename);
    }

    os << "# R Z BR BZ Bphi Bpol Bmag\n";

    for (int i = 0; i < nR; ++i) {
        for (int j = 0; j < nZ; ++j) {
			const int k = j * nR + i;
            const double R    = Rvec(i);
            const double Z    = Zvec(j);
            const double BR   = BR2D(k);
            const double BZ   = BZ2D(k);
            const double Bphi = BT*R0/R;

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

void WriteWall2D(IdsNs::IDS& ids,
                 double time,
                 std::string filename)
{
    ids._wall.getSlice(time, PREVIOUS_SAMPLE);

    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("WriteWall2D: could not open output file " + filename);
    }

    os << "# type segment R Z\n";

    int segment = 0;

    const int ndesc = static_cast<int>(ids._wall.description_2d.size());

    std::cout << "WriteWall2D: number of description_2d entries = " << ndesc << "\n";

    for (int i1 = 0; i1 < ndesc; ++i1) {
        auto& desc = ids._wall.description_2d(i1);

        // ============================================================
        // Limiter outlines
        // ============================================================
        const int nlim = static_cast<int>(desc.limiter.unit.size());
        std::cout << "  description_2d(" << i1 << "): limiter units = " << nlim << "\n";

        for (int i2 = 0; i2 < nlim; ++i2) {
            auto& unit = desc.limiter.unit(i2);
            auto& r = unit.outline.r;
            auto& z = unit.outline.z;

            const int n = static_cast<int>(r.size());
            if (n == 0) continue;
            if (static_cast<int>(z.size()) != n) {
                throw std::runtime_error("WriteWall2D: limiter outline r/z size mismatch.");
            }

            for (int k = 0; k < n; ++k) {
                os << 0 << " " << segment << " " << r(k) << " " << z(k) << "\n";
            }
            os << "\n";
            ++segment;
        }

        // ============================================================
        // Vessel annular outlines
        // ============================================================
        const int nves = static_cast<int>(desc.vessel.unit.size());
        std::cout << "  description_2d(" << i1 << "): vessel units = " << nves << "\n";

        for (int i2 = 0; i2 < nves; ++i2) {
            auto& unit = desc.vessel.unit(i2);

            // ----------------------------
            // Inner outline
            // ----------------------------
            auto& rin = unit.annular.outline_inner.r;
            auto& zin = unit.annular.outline_inner.z;

            const int nin = static_cast<int>(rin.size());
            if (nin > 0) {
                if (static_cast<int>(zin.size()) != nin) {
                    throw std::runtime_error("WriteWall2D: vessel inner outline r/z size mismatch.");
                }

                for (int k = 0; k < nin; ++k) {
                    os << 1 << " " << segment << " " << rin(k) << " " << zin(k) << "\n";
                }
                os << "\n";
                ++segment;
            }

            // ----------------------------
            // Outer outline
            // ----------------------------
            auto& rout = unit.annular.outline_outer.r;
            auto& zout = unit.annular.outline_outer.z;

            const int nout = static_cast<int>(rout.size());
            if (nout > 0) {
                if (static_cast<int>(zout.size()) != nout) {
                    throw std::runtime_error("WriteWall2D: vessel outer outline r/z size mismatch.");
                }

                for (int k = 0; k < nout; ++k) {
                    os << 2 << " " << segment << " " << rout(k) << " " << zout(k) << "\n";
                }
                os << "\n";
                ++segment;
            }
        }
    }

    std::cout << "WriteWall2D: wrote " << segment << " contour segments to " << filename << "\n";
}