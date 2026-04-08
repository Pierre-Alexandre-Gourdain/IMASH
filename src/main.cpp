// Include the Access Layer
#include "ALClasses.h"
#include <iostream>
#include <vector>

#include <iostream>
#include <fstream>

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


template <typename T>
void print_type(const T&, const char* name)
{
    std::cout << name << " type = " << typeid(T).name() << "\n";
}

int main() {
    IdsNs::IDS ids;

    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);
    ids._edge_profiles.getSlice(60.0, PREVIOUS_SAMPLE);

    auto& g    = ids._edge_profiles.ggd_fast(0);
    auto& e    = g.electrons;
    auto& ion  = g.ion;
    auto& grid = ids._edge_profiles.grid_ggd;

    auto& ion0   = ion(0);
    auto& edens  = e.density;
    auto& etemp  = e.temperature;
    auto& grid0  = grid(0);
    auto& space  = grid0.space;

    std::cout << "ggd extent         = " << ids._edge_profiles.ggd.extent(0) << "\n";
    std::cout << "grid extent        = " << grid.extent(0) << "\n";
    std::cout << "ion extent         = " << ion.extent(0) << "\n";
    std::cout << "space extent       = " << space.extent(0) << "\n";

    print_type(edens, "e.density");
    print_type(etemp, "e.temperature");
    print_type(ion0,  "ion(0)");
    print_type(space, "grid0.space");

std::cout << "edens extent(0) = " << edens.extent(0) << "\n";
std::cout << "etemp extent(0) = " << etemp.extent(0) << "\n";
    return 0;
}