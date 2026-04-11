#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "ALClasses.h"

struct EdgeCutPoint
{
    double R;
    double Z;
    double ne;
    double Te;
};

struct EdgeCellPoint
{
    double R;
    double Z;
    double ne;
    double Te;
};

struct EqGrid2D
{
    std::vector<double> Rvals;
    std::vector<double> Zvals;
    std::vector<double> rho;   // flattened with k = j*nR + i
    int nR = 0;
    int nZ = 0;
};

struct EdgeParamPoint
{
    double R;
    double Z;
    double rho;
    double theta;
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
void WriteEdge2D(IdsNs::IDS& ids, double time, const std::string& filename);
void WriteCore2D(IdsNs::IDS& ids,double time,const std::string& filename);
void WriteBfield2D(IdsNs::IDS& ids,
                   double time,
                   std::string filename);
void WriteWall2D(IdsNs::IDS& ids,
                 double time,
                 std::string filename);
void WriteMagneticAxis(IdsNs::IDS& ids,
                       double time,
                       const std::string& filename);
void WriteEdgeOnGridRhoTheta(IdsNs::IDS& ids,
                             double time,
                             const std::string& filename,
                             double sigma_rho = 0.1,
                             double sigma_theta = 0.1);