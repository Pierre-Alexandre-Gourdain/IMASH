// Include the Access Layer
#include "ALClasses.h"
#include "IMASH.H"

#include <iostream>
#include <vector>
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
void print_subset_indices(T& subset, const std::string& name)
{
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

int main()
{
    try {
        const std::string uri = "imas:hdf5?path=/home/pag/iter/134174/117";
        const double time = 60.0;

        IdsNs::IDS ids;
        ids.open(uri, OPEN_PULSE);

        WriteEdge2D(ids, time, "edge_cells_2d.dat");
        std::cout << "Wrote edge_cells_2d.dat\n";

        WriteCore2D(ids, time, "core_on_equilibrium_grid.dat");
        std::cout << "Wrote core_on_equilibrium_grid.dat\n";

        WriteBfield2D(ids, time, "bfield.dat");
        std::cout << "Wrote bfield.dat\n";

        WriteEdgeOnGridRhoTheta(ids, time, "edge_cells_2d_interpolated.dat", 0.01);
        std::cout << "Wrote edge_cells_2d_interpolated.dat\n";

        /*
        WriteWall2D(ids, time, "wall_2d.dat");
        std::cout << "Wrote wall_2d.dat\n";
        */
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}