#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "ALClasses.h"

struct EdgeCutPoint
{
    double R;
    double Z;
    double ne;
    double Te;
};

struct EdgeCut
{
    std::string name;
    std::vector<EdgeCutPoint> points;
};

EdgeCut ExtractEdgeCut(IdsNs::IDS& ids, int subset_index);

void PrintEdgeCut(const EdgeCut& cut);
void WriteEdgeCut(const EdgeCut& cut, const std::string& filename);
int FindSubsetIndexByName(IdsNs::IDS& ids, const std::string& name);