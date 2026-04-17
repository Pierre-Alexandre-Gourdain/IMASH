#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <cstddef>
#include <cmath>

// ============================================================
// Basic point types
// ============================================================

struct EdgeCutPoint
{
    double R  = 0.0;
    double Z  = 0.0;
    double ne = 0.0;
    double Te = 0.0;
};

struct EdgeCellPoint
{
    double R  = 0.0;
    double Z  = 0.0;
    double ne = 0.0;
    double Te = 0.0;
};

struct EdgeParamPoint
{
    double R     = 0.0;
    double Z     = 0.0;
    double rho   = 0.0;
    double theta = 0.0;
    double ne    = 0.0;
    double Te    = 0.0;
};

// ============================================================
// Magnetic-axis helper
// ============================================================

struct MagneticAxisPoint
{
    double R = 0.0;
    double Z = 0.0;
};

// ============================================================
// Axisymmetric structured grid in (R,Z)
// ============================================================

class AxisymmetricGrid2D
{
public:
    std::vector<double> Rvals;   // size nR
    std::vector<double> Zvals;   // size nZ

    int nR = 0;
    int nZ = 0;

    AxisymmetricGrid2D() = default;
    AxisymmetricGrid2D(int nR_, int nZ_);

    void resize(int nR_, int nZ_);
    void clear();

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] int size() const noexcept;

    [[nodiscard]] bool validIndex(int iR, int iZ) const noexcept;
    [[nodiscard]] int index(int iR, int iZ) const;
};

// ============================================================
// Scalar field defined on the structured (R,Z) grid
// Flattened with k = iR + nR * iZ
// ============================================================

class ScalarField2D
{
public:
    std::string name;
    std::string units;

    int nR = 0;
    int nZ = 0;
    std::vector<double> values;

    ScalarField2D() = default;
    ScalarField2D(std::string name_, std::string units_, int nR_, int nZ_);

    void resize(int nR_, int nZ_);
    void clear();

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] int size() const noexcept;

    [[nodiscard]] bool validIndex(int iR, int iZ) const noexcept;
    [[nodiscard]] int index(int iR, int iZ) const;

    double& operator()(int iR, int iZ);
    const double& operator()(int iR, int iZ) const;
};

// ============================================================
// Tokamak vector field on structured (R,Z) grid
// Components are (BR, BZ, Bphi)
// ============================================================

class VectorFieldTokamak2D
{
public:
    std::string name  = "B";
    std::string units = "T";

    int nR = 0;
    int nZ = 0;

    std::vector<double> BR;
    std::vector<double> BZ;
    std::vector<double> Bphi;

    VectorFieldTokamak2D() = default;
    VectorFieldTokamak2D(int nR_, int nZ_);

    void resize(int nR_, int nZ_);
    void clear();

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] int size() const noexcept;

    [[nodiscard]] bool validIndex(int iR, int iZ) const noexcept;
    [[nodiscard]] int index(int iR, int iZ) const;

    double& br(int iR, int iZ);
    double& bz(int iR, int iZ);
    double& bphi(int iR, int iZ);

    const double& br(int iR, int iZ) const;
    const double& bz(int iR, int iZ) const;
    const double& bphi(int iR, int iZ) const;
};

// ============================================================
// Legacy-compatible equilibrium-grid container
// ============================================================

struct EqGrid2D
{
    std::vector<double> Rvals;
    std::vector<double> Zvals;
    std::vector<double> rho;   // flattened with k = iR + nR * iZ
    int nR = 0;
    int nZ = 0;
};

// ============================================================
// 2D polyline in (R,Z), useful for wall contours
// ============================================================

class PolylineRZ
{
public:
    std::string name;
    int type = -1;      // e.g. 0=limiter, 1=vessel_inner, 2=vessel_outer
    bool closed = false;

    std::vector<double> R;
    std::vector<double> Z;

    PolylineRZ() = default;
    explicit PolylineRZ(std::string name_);

    void clear();
    void addPoint(double R_, double Z_);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] int size() const noexcept;
};

// ============================================================
// STL / CAD style triangulated surface
// ============================================================

class TriangleSurfaceMesh
{
public:
    std::string name;
    std::string units = "m";

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    std::vector<int> tri_i;
    std::vector<int> tri_j;
    std::vector<int> tri_k;

    TriangleSurfaceMesh() = default;
    explicit TriangleSurfaceMesh(std::string name_);

    void clear();
    void addVertex(double x_, double y_, double z_);
    void addTriangle(int i, int j, int k);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] int numVertices() const noexcept;
    [[nodiscard]] int numTriangles() const noexcept;
};

// ============================================================
// Edge cut container
// ============================================================

struct EdgeCut
{
    std::string name;
    std::vector<EdgeCutPoint> points;

    [[nodiscard]] bool empty() const noexcept { return points.empty(); }
    [[nodiscard]] int size() const noexcept { return static_cast<int>(points.size()); }
};

// ============================================================
// Plasma data on axisymmetric equilibrium grid
// ============================================================

class PlasmaData
{
public:
    AxisymmetricGrid2D grid;

    ScalarField2D rho;

    ScalarField2D ne_core;
    ScalarField2D Te_core;

    ScalarField2D ne_edge;
    ScalarField2D Te_edge;

    ScalarField2D ne;
    ScalarField2D Te;

    VectorFieldTokamak2D B;

    // Optional masks / provenance
    ScalarField2D core_mask;
    ScalarField2D edge_mask;
    ScalarField2D source_mask;

    PlasmaData() = default;
    PlasmaData(int nR, int nZ);

    void resize(int nR, int nZ);
    void clear();

    [[nodiscard]] bool empty() const noexcept;
};

// ============================================================
// Wall geometry can come from IMAS contours and/or STL meshes
// ============================================================

class WallData
{
public:
    std::vector<PolylineRZ> rz_contours;
    std::vector<TriangleSurfaceMesh> surface_meshes;

    void clear();
    void addContour(const PolylineRZ& contour);
    void addSurfaceMesh(const TriangleSurfaceMesh& mesh);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] bool hasRZContours() const noexcept;
    [[nodiscard]] bool hasSurfaceMeshes() const noexcept;
};

// ============================================================
// Diagnostics / launcher geometry
// ============================================================

class WaveguideAntenna
{
public:
    std::string name;

    // Position in tokamak coordinates
    double R   = 0.0;
    double Z   = 0.0;
    double phi = 0.0;

    // Aperture
    double diameter = 0.0;

    // Pointing direction in cylindrical tokamak basis
    double dir_R   = 0.0;
    double dir_Z   = 0.0;
    double dir_phi = 0.0;

    // Optional operating parameters
    double frequency_hz = 0.0;
    std::string polarization;

    [[nodiscard]] bool validPosition() const noexcept;
    [[nodiscard]] bool validDirection() const noexcept;
};

class DiagnosticData
{
public:
    std::vector<WaveguideAntenna> waveguides;

    void clear();
    void addWaveguide(const WaveguideAntenna& wg);

    [[nodiscard]] bool empty() const noexcept;
};

// ============================================================
// Metadata and top-level scene
// ============================================================

class SceneMetadata
{
public:
    std::string imas_uri;
    double time = 0.0;

    std::string length_units      = "m";
    std::string density_units     = "m^-3";
    std::string temperature_units = "eV";
    std::string magnetic_units    = "T";

    bool has_plasma      = false;
    bool has_wall        = false;
    bool has_diagnostics = false;
};

class TokamakScene
{
public:
    PlasmaData plasma;
    WallData wall;
    DiagnosticData diagnostics;
    SceneMetadata meta;

    [[nodiscard]] bool empty() const noexcept;
};

// ============================================================
// Optional runtime / debug options for later loaders/writers
// ============================================================

class SceneOptions
{
public:
    bool load_plasma           = true;
    bool load_wall_from_imas   = true;
    bool load_wall_from_stl    = false;
    bool load_diagnostics      = true;

    bool merge_density         = true;
    bool merge_temperature     = true;

    bool dump_plasma           = false;
    bool dump_wall             = false;
    bool dump_diagnostics      = false;

    std::string dump_directory = ".";
    std::string dump_prefix    = "scene";

    std::vector<std::string> stl_wall_files;
};
