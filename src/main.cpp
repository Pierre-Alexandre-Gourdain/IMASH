#include "ALClasses.h"
#include "IMASH.H"

#include <iostream>
#include <string>
#include <cstdlib>   // std::stod
#include <filesystem>

int main(int argc, char* argv[])
{
    try {
        // --------------------------------------------------
        // Parse arguments
        // --------------------------------------------------
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0]
                      << " <imas_path> <time> [output_dir]\n";
            return 1;
        }

        const std::string imas_path = argv[1];
        const double time = std::stod(argv[2]);

        std::string output_dir = ".";
        if (argc >= 4) {
            output_dir = argv[3];
        }

        std::filesystem::create_directories(output_dir);

        // --------------------------------------------------
        // Open IMAS
        // --------------------------------------------------
        IdsNs::IDS ids;

        std::string uri = "imas:hdf5?path=" + imas_path;
        ids.open(uri, OPEN_PULSE);

        // --------------------------------------------------
        // Core / Edge / B-field outputs
        // --------------------------------------------------
        {
            std::string fname = output_dir + "/edge_cells_2d.dat";
            WriteEdge2D(ids, time, fname);
            std::cout << "Wrote " << fname << "\n";
        }

        {
            std::string fname = output_dir + "/core_on_equilibrium_grid.dat";
            WriteCore2D(ids, time, fname);
            std::cout << "Wrote " << fname << "\n";
        }

        {
            std::string fname = output_dir + "/bfield.dat";
            WriteBfield2D(ids, time, fname);
            std::cout << "Wrote " << fname << "\n";
        }

        {
            std::string fname =
                output_dir + "/edge_cells_2d_interpolated.dat";

            WriteEdgeOnGridRhoTheta(ids, time, fname, 0.01);
            std::cout << "Wrote " << fname << "\n";
        }

        /*
        {
            std::string fname = output_dir + "/wall_2d.dat";
            WriteWall2D(ids, time, fname);
            std::cout << "Wrote " << fname << "\n";
        }
        */

    }
    catch (const std::exception& ex) {
        std::cerr << "Fatal error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}