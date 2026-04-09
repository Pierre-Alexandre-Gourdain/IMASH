// Include the Access Layer
#include "ALClasses.h"
#include "IMASH.H"
#include <iostream>
#include <vector>

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <typeinfo>
#include <cstdlib>
#include <cxxabi.h>


template <typename T>
void print_type(const T&, const char* name)
{
    int status = 0;
    std::unique_ptr<char, void(*)(void*)> demangled(
        abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status),
        std::free
    );

    std::cout << name << " type = "
              << ((status == 0 && demangled) ? demangled.get() : typeid(T).name())
              << "\n";
}

template <typename T>
void print_subset_indices(T& subset, const std::string& name){
    std::cout << "\nSubset: " << name << "\n";
    std::cout << "elements = " << subset.element.extent(0) << "\n";

    for (int i = subset.element.lbound(0); i <= subset.element.ubound(0); ++i) {
        auto& el = subset.element(i);
        auto& objs = el.object;

        std::cout << "  element " << i << ": ";
        for (int j = objs.lbound(0); j <= objs.ubound(0); ++j) {
            std::cout << "(dim=" << objs(j).dimension
                      << ", idx=" << objs(j).index << ") ";
        }
        std::cout << "\n";
    }
}


int main1(int argc, char *argv[]) {
    IdsNs::IDS ids;

    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
    ids._core_profiles.getSlice(60.0, PREVIOUS_SAMPLE);

    auto& profiles = ids._core_profiles.profiles_1d;

    if (profiles.extent(0) == 0) {
        std::cout << "No profiles_1d found\n";
        return 1;
    }

    auto& prof = profiles(0);

    if (prof.ion.extent(0) == 0) {
        std::cout << "No ion species found\n";
        return 1;
    }

    auto& ion0 = prof.ion(0);

    auto& ni = ion0.density;
    auto& Te = prof.electrons.temperature;
    auto& rho = prof.grid.rho_tor_norm;

    std::cout << "ni points  = " << ni.extent(0)  << "\n";
    std::cout << "Te points  = " << Te.extent(0)  << "\n";
    std::cout << "rho points = " << rho.extent(0) << "\n";

    std::ofstream out("core_profiles_t60.dat");
    out << "# rho_tor_norm ni Te\n";

    int i0 = rho.lbound(0);
    int i1 = rho.ubound(0);

    for (int i = i0; i <= i1; ++i) {
        out << rho(i) << " " << ni(i) << " " << Te(i) << "\n";
    }

    out.close();

    std::cout << "Wrote core_profiles_t60.dat\n";
    return 0;
}

int main2(int argc, char *argv[]) {
    IdsNs::IDS ids;

    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
    ids._equilibrium.getSlice(60.0, PREVIOUS_SAMPLE);

    auto& ts0  = ids._equilibrium.time_slice(0);
    auto& p2d  = ts0.profiles_2d(0);

    auto& R    = p2d.grid.dim1;
    auto& Z    = p2d.grid.dim2;
    auto& psi  = p2d.psi;

    std::cout << "R extent   = " << R.extent(0)   << "\n";
    std::cout << "Z extent   = " << Z.extent(0)   << "\n";
    std::cout << "psi rank-2 = "
              << psi.extent(0) << " x " << psi.extent(1) << "\n";

    std::cout << "\nR grid:\n";
    for (int i = R.lbound(0); i <= R.ubound(0); ++i) {
        std::cout << "R(" << i << ") = " << R(i) << "\n";
    }

    std::cout << "\nZ grid:\n";
    for (int j = Z.lbound(0); j <= Z.ubound(0); ++j) {
        std::cout << "Z(" << j << ") = " << Z(j) << "\n";
    }

    std::cout << "\nSample psi values:\n";
    for (int i = psi.lbound(0); i <= std::min(psi.lbound(0) + 2, psi.ubound(0)); ++i) {
        for (int j = psi.lbound(1); j <= std::min(psi.lbound(1) + 2, psi.ubound(1)); ++j) {
            std::cout << "psi(" << i << "," << j << ") = " << psi(i,j) << "\n";
        }
    }

    return 0;
}


int main3() {
    IdsNs::IDS ids;

    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
    ids._edge_profiles.getSlice(60.0, CLOSEST_SAMPLE);

	auto& g      = ids._edge_profiles.ggd_fast(0);
	auto& e      = g.electrons;
	auto& grid0  = ids._edge_profiles.grid_ggd(0);
	auto& spaces = grid0.space;
    auto& edens  = e.density;
    auto& etemp  = e.temperature;
	
	
	std::cout << "edens extent(0) = " << edens.extent(0) << "\n";
	std::cout << "etemp extent(0) = " << etemp.extent(0) << "\n";

	std::cout << "edens bounds = ["
			  << edens.lbound(0) << ", "
			  << edens.ubound(0) << "]\n";

	std::cout << "etemp bounds = ["
			  << etemp.lbound(0) << ", "
			  << etemp.ubound(0) << "]\n";

	std::cout << "spaces extent = " << spaces.extent(0) << "\n";

	auto& space0 = spaces(0);
	auto& objects = space0.objects_per_dimension;

	std::cout << "objects_per_dimension extent = "
			  << objects.extent(0) << "\n";

	print_type(objects, "space0.objects_per_dimension");

	if (objects.extent(0) > 0) print_type(objects(0), "objects(0)");
	if (objects.extent(0) > 1) print_type(objects(1), "objects(1)");
	if (objects.extent(0) > 2) print_type(objects(2), "objects(2)");
	
	auto& obj0 = objects(0);
	auto& obj1 = objects(1);
	auto& obj2 = objects(2);

	print_type(obj0, "obj0");
	print_type(obj1, "obj1");
	print_type(obj2, "obj2");

	auto& o0 = obj0.object;
	auto& o1 = obj1.object;
	auto& o2 = obj2.object;

	std::cout << "obj0.object extent = " << o0.extent(0) << "\n";
	std::cout << "obj1.object extent = " << o1.extent(0) << "\n";
	std::cout << "obj2.object extent = " << o2.extent(0) << "\n";	
	
	auto& subsets = grid0.grid_subset;
	std::cout << "grid_subset extent = " << subsets.extent(0) << "\n";
	print_type(subsets, "grid0.grid_subset");

	if (subsets.extent(0) > 0) print_type(subsets(0), "grid_subset(0)");
	if (subsets.extent(0) > 1) print_type(subsets(1), "grid_subset(1)");
	if (subsets.extent(0) > 2) print_type(subsets(2), "grid_subset(2)");
	
	auto& subset0 = subsets(0);

	std::cout << "grid_subset extent = " << subsets.extent(0) << "\n";

	/* for (int s = subsets.lbound(0); s <= subsets.ubound(0); ++s) {
		auto& subset = subsets(s);

		std::cout << "\n=== subset " << s << " ===\n";
		std::cout << "subset.element extent = "
				  << subset.element.extent(0) << "\n";

		// Identifier fields
		std::cout << "identifier.name        = " << subset.identifier.name << "\n";
		std::cout << "identifier.index       = " << subset.identifier.index << "\n";
		std::cout << "identifier.description = " << subset.identifier.description << "\n";

		if (subset.element.extent(0) > 0) {
			auto& el0 = subset.element(0);
			auto& objs = el0.object;

			std::cout << "subset.element(0).object extent = "
					  << objs.extent(0) << "\n";

			if (objs.extent(0) > 0) {
				std::cout << "first object dimension = "
						  << objs(0).dimension << "\n";
				std::cout << "first object index     = "
						  << objs(0).index << "\n";
			}
		}
	} */
	/* print_subset_indices(subsets(16), "inner_PFR_wall");
	print_subset_indices(subsets(18), "outer_PFR_wall");
	print_subset_indices(subsets(20), "inner_baffle");
	print_subset_indices(subsets(22), "outer_baffle");
	print_subset_indices(subsets(23), "inner_midplane");
	print_subset_indices(subsets(24), "outer_midplane"); */
	auto& subset = subsets(24);

	for (int i = subset.element.lbound(0); i <= subset.element.ubound(0); ++i) {
		auto& el = subset.element(i);
		auto& objs = el.object;

		std::cout << "i = " << i;
		for (int j = objs.lbound(0); j <= objs.ubound(0); ++j) {
			std::cout << "  dim=" << objs(j).dimension
					  << " idx=" << objs(j).index;
		}
		std::cout << "\n";
	}

	auto& faces_bucket = objects(1).object;
	auto& nodes_bucket = objects(0).object;

	std::cout << "faces_bucket extent = " << faces_bucket.extent(0) << "\n";
	std::cout << "nodes_bucket extent = " << nodes_bucket.extent(0) << "\n";

	int idx = 623;   // first outer_midplane face index

	// IMAS subset indices are very likely 1-based, so try idx-1 first
	auto& face = faces_bucket(idx - 1);

	print_type(face, "face");

	// try likely members one by one:
	auto& face_nodes = face.nodes;
	print_type(face_nodes, "face.nodes");
	std::cout << "face.nodes extent = " << face_nodes.extent(0) << "\n";

	int n0 = face_nodes(0);
	int n1 = face_nodes(1);

	std::cout << "face node indices = " << n0 << ", " << n1 << "\n";

	// very likely 1-based, so try -1 first
	auto& node0 = nodes_bucket(n0 - 1);
	auto& node1 = nodes_bucket(n1 - 1);

	print_type(node0, "node0");
	print_type(node1, "node1");
	auto& c0 = node0.geometry;
	auto& c1 = node1.geometry;
	auto& x0 = node0.geometry;
	print_type(x0, "node0.geometry");
	auto& g0 = node0.geometry;
	auto& g1 = node1.geometry;

	print_type(g0, "node0.geometry");
	print_type(g1, "node1.geometry");

	std::cout << "node0.geometry extent = " << g0.extent(0) << "\n";
	std::cout << "node1.geometry extent = " << g1.extent(0) << "\n";

	for (int i = g0.lbound(0); i <= g0.ubound(0); ++i) {
		std::cout << "node0.geometry(" << i << ") = " << g0(i) << "\n";
	}
	for (int i = g1.lbound(0); i <= g1.ubound(0); ++i) {
		std::cout << "node1.geometry(" << i << ") = " << g1(i) << "\n";
	}
	
	print_type(edens(0), "edens(0)");
	print_type(etemp(0), "etemp(0)");

	// first guess
	auto& ne0 = edens(0).value;
	auto& Te0 = etemp(0).value;

	print_type(ne0, "edens(0).value");
	print_type(Te0, "etemp(0).value");	
	return 0;
}

#include <iostream>
#include <iomanip>

int main4() {
    IdsNs::IDS ids;

    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
    ids._edge_profiles.getSlice(60.0, CLOSEST_SAMPLE);

    auto& g      = ids._edge_profiles.ggd_fast(0);
    auto& e      = g.electrons;
    auto& edens  = e.density;
    auto& etemp  = e.temperature;

    auto& grid0   = ids._edge_profiles.grid_ggd(0);
    auto& spaces  = grid0.space;
    auto& space0  = spaces(0);
    auto& objects = space0.objects_per_dimension;

    auto& nodes_bucket = objects(0).object;   // nodes
    auto& faces_bucket = objects(1).object;   // faces

    auto& subsets = grid0.grid_subset;
    auto& subset  = subsets(24);              // outer_midplane

    std::cout << std::setprecision(10);
    std::cout << "# i  R_face  Z_face  ne  Te\n";

    for (int i = 0; i < edens.extent(0); ++i) {
        int face_idx = subset.element(i).object(0).index;

        // subset indices are effectively 1-based
        auto& face = faces_bucket(face_idx - 1);

        int n0 = face.nodes(0);
        int n1 = face.nodes(1);

        auto& node0 = nodes_bucket(n0 - 1);
        auto& node1 = nodes_bucket(n1 - 1);

        auto& g0 = node0.geometry;
        auto& g1 = node1.geometry;

        double Rf = 0.5 * (g0(0) + g1(0));
        double Zf = 0.5 * (g0(1) + g1(1));

        double ne = edens(i).value;
        double Te = etemp(i).value;

        std::cout << i << "  "
                  << Rf << "  "
                  << Zf << "  "
                  << ne << "  "
                  << Te << "\n";
    }

    return 0;
}

int main()
{
    IdsNs::IDS ids;

    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
    ids._edge_profiles.getSlice(60.0, CLOSEST_SAMPLE);

    try {
        int outer_idx = FindSubsetIndexByName(ids, "outer_midplane");
        int inner_idx = FindSubsetIndexByName(ids, "inner_midplane");

        EdgeCut outer_midplane = ExtractEdgeCut(ids, outer_idx);
        EdgeCut inner_midplane = ExtractEdgeCut(ids, inner_idx);

        WriteEdgeCut(outer_midplane, "outer_midplane_edge_profile.dat");
        WriteEdgeCut(inner_midplane, "inner_midplane_edge_profile.dat");

        PrintEdgeCut(outer_midplane);
        PrintEdgeCut(inner_midplane);
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}