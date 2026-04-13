#include "TokamakData.H"

// ============================================================
// AxisymmetricGrid2D
// ============================================================

AxisymmetricGrid2D::AxisymmetricGrid2D(int nR_, int nZ_)
{
    resize(nR_, nZ_);
}

void AxisymmetricGrid2D::resize(int nR_, int nZ_)
{
    if (nR_ < 0 || nZ_ < 0) {
        throw std::invalid_argument("AxisymmetricGrid2D::resize: negative size");
    }

    nR = nR_;
    nZ = nZ_;

    Rvals.resize(static_cast<std::size_t>(nR));
    Zvals.resize(static_cast<std::size_t>(nZ));
}

void AxisymmetricGrid2D::clear()
{
    Rvals.clear();
    Zvals.clear();
    nR = 0;
    nZ = 0;
}

bool AxisymmetricGrid2D::empty() const noexcept
{
    return (nR == 0 || nZ == 0);
}

int AxisymmetricGrid2D::size() const noexcept
{
    return nR * nZ;
}

bool AxisymmetricGrid2D::validIndex(int iR, int iZ) const noexcept
{
    return (0 <= iR && iR < nR && 0 <= iZ && iZ < nZ);
}

int AxisymmetricGrid2D::index(int iR, int iZ) const
{
    if (!validIndex(iR, iZ)) {
        throw std::out_of_range("AxisymmetricGrid2D::index: invalid (iR,iZ)");
    }
    return iR + nR * iZ;
}

// ============================================================
// ScalarField2D
// ============================================================

ScalarField2D::ScalarField2D(std::string name_,
                             std::string units_,
                             int nR_,
                             int nZ_)
    : name(std::move(name_)),
      units(std::move(units_))
{
    resize(nR_, nZ_);
}

void ScalarField2D::resize(int nR_, int nZ_)
{
    if (nR_ < 0 || nZ_ < 0) {
        throw std::invalid_argument("ScalarField2D::resize: negative size");
    }

    nR = nR_;
    nZ = nZ_;
    values.assign(static_cast<std::size_t>(nR * nZ), 0.0);
}

void ScalarField2D::clear()
{
    values.clear();
    nR = 0;
    nZ = 0;
}

bool ScalarField2D::empty() const noexcept
{
    return values.empty();
}

int ScalarField2D::size() const noexcept
{
    return nR * nZ;
}

bool ScalarField2D::validIndex(int iR, int iZ) const noexcept
{
    return (0 <= iR && iR < nR && 0 <= iZ && iZ < nZ);
}

int ScalarField2D::index(int iR, int iZ) const
{
    if (!validIndex(iR, iZ)) {
        throw std::out_of_range("ScalarField2D::index: invalid (iR,iZ)");
    }
    return iR + nR * iZ;
}

double& ScalarField2D::operator()(int iR, int iZ)
{
    return values.at(static_cast<std::size_t>(index(iR, iZ)));
}

const double& ScalarField2D::operator()(int iR, int iZ) const
{
    return values.at(static_cast<std::size_t>(index(iR, iZ)));
}

// ============================================================
// VectorFieldTokamak2D
// ============================================================

VectorFieldTokamak2D::VectorFieldTokamak2D(int nR_, int nZ_)
{
    resize(nR_, nZ_);
}

void VectorFieldTokamak2D::resize(int nR_, int nZ_)
{
    if (nR_ < 0 || nZ_ < 0) {
        throw std::invalid_argument("VectorFieldTokamak2D::resize: negative size");
    }

    nR = nR_;
    nZ = nZ_;
    const std::size_t n = static_cast<std::size_t>(nR * nZ);

    BR.assign(n, 0.0);
    BZ.assign(n, 0.0);
    Bphi.assign(n, 0.0);
}

void VectorFieldTokamak2D::clear()
{
    BR.clear();
    BZ.clear();
    Bphi.clear();
    nR = 0;
    nZ = 0;
}

bool VectorFieldTokamak2D::empty() const noexcept
{
    return BR.empty() && BZ.empty() && Bphi.empty();
}

int VectorFieldTokamak2D::size() const noexcept
{
    return nR * nZ;
}

bool VectorFieldTokamak2D::validIndex(int iR, int iZ) const noexcept
{
    return (0 <= iR && iR < nR && 0 <= iZ && iZ < nZ);
}

int VectorFieldTokamak2D::index(int iR, int iZ) const
{
    if (!validIndex(iR, iZ)) {
        throw std::out_of_range("VectorFieldTokamak2D::index: invalid (iR,iZ)");
    }
    return iR + nR * iZ;
}

double& VectorFieldTokamak2D::br(int iR, int iZ)
{
    return BR.at(static_cast<std::size_t>(index(iR, iZ)));
}

double& VectorFieldTokamak2D::bz(int iR, int iZ)
{
    return BZ.at(static_cast<std::size_t>(index(iR, iZ)));
}

double& VectorFieldTokamak2D::bphi(int iR, int iZ)
{
    return Bphi.at(static_cast<std::size_t>(index(iR, iZ)));
}

const double& VectorFieldTokamak2D::br(int iR, int iZ) const
{
    return BR.at(static_cast<std::size_t>(index(iR, iZ)));
}

const double& VectorFieldTokamak2D::bz(int iR, int iZ) const
{
    return BZ.at(static_cast<std::size_t>(index(iR, iZ)));
}

const double& VectorFieldTokamak2D::bphi(int iR, int iZ) const
{
    return Bphi.at(static_cast<std::size_t>(index(iR, iZ)));
}

// ============================================================
// PolylineRZ
// ============================================================

PolylineRZ::PolylineRZ(std::string name_)
    : name(std::move(name_))
{
}

void PolylineRZ::clear()
{
    R.clear();
    Z.clear();
    closed = false;
}

void PolylineRZ::addPoint(double R_, double Z_)
{
    R.push_back(R_);
    Z.push_back(Z_);
}

bool PolylineRZ::empty() const noexcept
{
    return R.empty();
}

int PolylineRZ::size() const noexcept
{
    return static_cast<int>(R.size());
}

// ============================================================
// TriangleSurfaceMesh
// ============================================================

TriangleSurfaceMesh::TriangleSurfaceMesh(std::string name_)
    : name(std::move(name_))
{
}

void TriangleSurfaceMesh::clear()
{
    x.clear();
    y.clear();
    z.clear();
    tri_i.clear();
    tri_j.clear();
    tri_k.clear();
}

void TriangleSurfaceMesh::addVertex(double x_, double y_, double z_)
{
    x.push_back(x_);
    y.push_back(y_);
    z.push_back(z_);
}

void TriangleSurfaceMesh::addTriangle(int i, int j, int k)
{
    tri_i.push_back(i);
    tri_j.push_back(j);
    tri_k.push_back(k);
}

bool TriangleSurfaceMesh::empty() const noexcept
{
    return x.empty() || tri_i.empty();
}

int TriangleSurfaceMesh::numVertices() const noexcept
{
    return static_cast<int>(x.size());
}

int TriangleSurfaceMesh::numTriangles() const noexcept
{
    return static_cast<int>(tri_i.size());
}

// ============================================================
// PlasmaData
// ============================================================

PlasmaData::PlasmaData(int nR, int nZ)
{
    resize(nR, nZ);
}

void PlasmaData::resize(int nR, int nZ)
{
    grid.resize(nR, nZ);

    rho       = ScalarField2D("rho",       "",     nR, nZ);

    ne_core   = ScalarField2D("ne_core",   "m^-3", nR, nZ);
    Te_core   = ScalarField2D("Te_core",   "eV",   nR, nZ);

    ne_edge   = ScalarField2D("ne_edge",   "m^-3", nR, nZ);
    Te_edge   = ScalarField2D("Te_edge",   "eV",   nR, nZ);

    ne        = ScalarField2D("ne",        "m^-3", nR, nZ);
    Te        = ScalarField2D("Te",        "eV",   nR, nZ);

    B         = VectorFieldTokamak2D(nR, nZ);

    core_mask   = ScalarField2D("core_mask",   "", nR, nZ);
    edge_mask   = ScalarField2D("edge_mask",   "", nR, nZ);
    source_mask = ScalarField2D("source_mask", "", nR, nZ);
}

void PlasmaData::clear()
{
    grid.clear();

    rho.clear();

    ne_core.clear();
    Te_core.clear();

    ne_edge.clear();
    Te_edge.clear();

    ne.clear();
    Te.clear();

    B.clear();

    core_mask.clear();
    edge_mask.clear();
    source_mask.clear();
}

bool PlasmaData::empty() const noexcept
{
    return grid.empty();
}

// ============================================================
// WallData
// ============================================================

void WallData::clear()
{
    rz_contours.clear();
    surface_meshes.clear();
}

void WallData::addContour(const PolylineRZ& contour)
{
    rz_contours.push_back(contour);
}

void WallData::addSurfaceMesh(const TriangleSurfaceMesh& mesh)
{
    surface_meshes.push_back(mesh);
}

bool WallData::empty() const noexcept
{
    return rz_contours.empty() && surface_meshes.empty();
}

bool WallData::hasRZContours() const noexcept
{
    return !rz_contours.empty();
}

bool WallData::hasSurfaceMeshes() const noexcept
{
    return !surface_meshes.empty();
}

// ============================================================
// WaveguideAntenna / DiagnosticData
// ============================================================

bool WaveguideAntenna::validPosition() const noexcept
{
    return std::isfinite(R) && std::isfinite(Z) && std::isfinite(phi);
}

bool WaveguideAntenna::validDirection() const noexcept
{
    return std::isfinite(dir_R) && std::isfinite(dir_Z) && std::isfinite(dir_phi) &&
           !(dir_R == 0.0 && dir_Z == 0.0 && dir_phi == 0.0);
}

void DiagnosticData::clear()
{
    waveguides.clear();
}

void DiagnosticData::addWaveguide(const WaveguideAntenna& wg)
{
    waveguides.push_back(wg);
}

bool DiagnosticData::empty() const noexcept
{
    return waveguides.empty();
}

// ============================================================
// TokamakScene
// ============================================================

bool TokamakScene::empty() const noexcept
{
    return plasma.empty() && wall.empty() && diagnostics.empty();
}