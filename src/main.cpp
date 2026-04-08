// Include the Access Layer
#include "ALClasses.h"
#include <iostream>

int main(int argc, char *argv[]) {
    // Create a new IDS object
    IdsNs::IDS ids;
    // Open the database entry by providing an IMAS URI
    // ids.open("imas:mdsplus?user=public;pulse=131024;run=41;database=ITER;version=3", OPEN_PULSE);
    ids.open("imas:hdf5?path=/home/pag/iter/134174/117", OPEN_PULSE);

	ids._core_profiles.getSlice(60, PREVIOUS_SAMPLE);
	
    return 0;
}

