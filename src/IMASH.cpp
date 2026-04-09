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
