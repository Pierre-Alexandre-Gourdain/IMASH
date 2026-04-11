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


/* int main1()
{
    try {
        IdsNs::IDS ids;

        ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
        ids._edge_profiles.getSlice(60.0, CLOSEST_SAMPLE);

        WriteEdgeCells2D(ids, "edge_cells_2d.dat");

        std::cout << "Wrote edge_cells_2d.dat\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
} */

int main()
{
    try {
        IdsNs::IDS ids;

        ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
		
        WriteEdge2D(ids, 60.0, "edge_cells_2d.dat");
        std::cout << "Wrote edge_cells_2d.dat\n";
		
        WriteCore2D(ids, 60.0,"core_on_equilibrium_grid.dat");
        std::cout << "Wrote core_on_equilibrium_grid.dat\n";
		
		WriteBfield2D(ids, 60.0, "bfield.dat");
        std::cout << "Wrote bfield.dat\n";
		
		WriteEdgeOnGridRhoTheta(ids, 60.0, "edge_cells_2d_interpolated.dat",.01);
        std::cout << "Wrote edge_cells_2d_interpolated.dat\n";
		
		/*WriteWall2D(ids, 60.0, "wall_2d.dat"); //needs testing
        std::cout << "Wrote wall_2d.dat\n"; */
		
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}