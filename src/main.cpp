// Include the Access Layer
#include "ALClasses.h"
#include <iostream>
#include <vector>

#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
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