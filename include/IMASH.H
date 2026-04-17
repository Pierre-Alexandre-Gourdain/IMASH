#pragma once

#include <string>
#include <vector>

#include "ALClasses.h"
#include "TokamakData.H"

// ============================================================
// Existing utility / lookup functions
// ============================================================

int FindSubsetIndexByName(IdsNs::IDS& ids, const std::string& name);


template <typename SubsetArray>
int FindSubsetSlotByIdentifierIndex(const SubsetArray& subsets,
                                    int identifier_index)
{
    for (int s = subsets.lbound(0); s <= subsets.ubound(0); ++s) {
        if (subsets(s).identifier.index == identifier_index) {
            return s;
        }
    }
    return -1;
}

template <typename FieldArray>
int FindFieldEntryForSubset(const FieldArray& field,
                            int subset_identifier_index)
{
    for (int i = field.lbound(0); i <= field.ubound(0); ++i) {
        if (field(i).grid_subset_index == subset_identifier_index) {
            return i;
        }
    }
    return -1;
}

// ============================================================
// Small numerical helpers used by the IMAS extraction routines
// ============================================================

double LinearInterp1D(const std::vector<double>& x,
                      const std::vector<double>& y,
                      double xq);

double WrapAngleDiff(double a, double b);

double ComputeTheta(double R, double Z, double R0, double Z0);

int FindBracket(const std::vector<double>& x, double xq);

double BilinearInterpEq(const EqGrid2D& g, double Rq, double Zq);

double LocalInterpRhoTheta(double rho_q,
                           double theta_q,
                           const std::vector<EdgeParamPoint>& src,
                           bool use_ne,
                           double sigma_rho,
                           double sigma_theta,
                           double max_drho,
                           double max_dtheta,
                           double min_weight_sum = 1.0e-12);

// ============================================================
// Edge-cut extraction / debug output
// ============================================================

EdgeCut ExtractEdgeCut(IdsNs::IDS& ids, int subset_index);

void PrintEdgeCut(const EdgeCut& cut);

void WriteEdgeCut(const EdgeCut& cut, const std::string& filename);

// ============================================================
// Structured equilibrium-grid helpers
// These preserve the equilibrium-grid flattening used in IMASH.cpp:
// k = j*nR + i
// ============================================================

EqGrid2D ExtractEqGridRho(IdsNs::IDS& ids, double time);

// ============================================================
// Edge data extraction
// ============================================================

std::vector<EdgeCellPoint> ExtractEdge2D(IdsNs::IDS& ids, double time);

std::vector<EdgeParamPoint>
BuildEdgeParamPoints(const std::vector<EdgeCellPoint>& edge_pts,
                     const EqGrid2D& eq,
                     double R0, double Z0);

// ============================================================
// Plasma extraction into reusable structures
// ============================================================

// Fills plasma.grid, plasma.rho, plasma.ne_core, plasma.Te_core
void FillCoreProfilesFromImas(IdsNs::IDS& ids,
                              double time,
                              PlasmaData& plasma);

// Fills plasma.B on the same equilibrium grid
void FillBfieldFromImas(IdsNs::IDS& ids,
                        double time,
                        PlasmaData& plasma);

// Extract edge values as unstructured cell-center points only
void FillEdgeCellPointsFromImas(IdsNs::IDS& ids,
                                double time,
                                std::vector<EdgeCellPoint>& edge_pts);

// Optional projection of edge data onto the equilibrium grid
void FillEdgeProfilesOnGridFromImas(IdsNs::IDS& ids,
                                    double time,
                                    PlasmaData& plasma,
                                    double sigma_rho = 0.1,
                                    double sigma_theta = 0.1);

// Optional merge step for combined fields
void MergePlasmaProfiles(PlasmaData& plasma);

// Convenience full plasma builder
PlasmaData ExtractPlasmaData(IdsNs::IDS& ids,
                             double time,
                             bool fill_core = true,
                             bool fill_bfield = true,
                             bool fill_edge_on_grid = false,
                             bool merge_density = false,
                             bool merge_temperature = false,
                             double sigma_rho = 0.1,
                             double sigma_theta = 0.1);

// ============================================================
// Wall extraction
// ============================================================

WallData ExtractWall2D(IdsNs::IDS& ids, double time);


MagneticAxisPoint ExtractMagneticAxis(IdsNs::IDS& ids, double time);

// ============================================================
// Debug / Python-friendly writers from data structures
// These preserve your file-based verification workflow
// ============================================================

void WriteCore2D(const PlasmaData& plasma,
                 const std::string& filename);

void WriteBfield2D(const PlasmaData& plasma,
                   const std::string& filename);

void WriteWall2D(const WallData& wall,
                 const std::string& filename);

void WriteMagneticAxis(const MagneticAxisPoint& axis,
                       const std::string& filename);

void WriteEdge2D(const std::vector<EdgeCellPoint>& edge_pts,
                 const std::string& filename);

void WriteEdgeOnGridRhoTheta(const PlasmaData& plasma,
                             const std::string& filename);

// ============================================================
// Backward-compatible wrappers
// These keep the current "read IMAS + write file" workflow alive
// ============================================================

void WriteCore2D(IdsNs::IDS& ids,
                 double time,
                 const std::string& filename);

void WriteEdge2D(IdsNs::IDS& ids,
                 double time,
                 const std::string& filename);

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