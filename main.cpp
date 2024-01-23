#include <fmt/core.h>

#include <libmesh/boundary_info.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/implicit_system.h>
#include <libmesh/libmesh.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_base.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/vector_value.h>
#include <libmesh/zero_function.h>

#include <igl/delaunay_triangulation.h>
#include <igl/isolines.h>

void assemble_distance(libMesh::EquationSystems & es, std::string const & name);

int main(int argc, char * argv[])
{
  libMesh::LibMeshInit init(argc, argv);

  libMesh::Mesh mesh{init.comm()};

  uint const n = 20;
  libMesh::MeshTools::Generation::build_square(
      mesh, n, n, 0.0, 1.0, 0.0, 1.0, libMesh::QUAD4);

  mesh.print_info();

  mesh.boundary_info->print_info();

  libMesh::EquationSystems es{mesh};

  libMesh::LinearImplicitSystem & system =
      es.add_system<libMesh::LinearImplicitSystem>("distance");

  system.add_variable("d", libMesh::FIRST, libMesh::LAGRANGE);

  libMesh::ZeroFunction<double> zero;

  system.get_dof_map().add_dirichlet_boundary(
      libMesh::DirichletBoundary{{3}, {0}, &zero});

  system.attach_assemble_function(assemble_distance);

  es.init();

  es.print_info();

  system.solve();

  libMesh::ExodusII_IO io{mesh};
  io.write_timestep("out.e", es, 1, 0.0);

  return 0;
}

void assemble_distance(libMesh::EquationSystems & es, std::string const & name)
{
  using namespace libMesh;

  assert(name == "distance");
  MeshBase const & mesh = es.get_mesh();
  uint const dim = mesh.mesh_dimension();

  ImplicitSystem & system = es.get_system<ImplicitSystem>(name);

  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  std::unique_ptr<FEBase> fe{FEBase::build(dim, fe_type)};
  QGauss qrule{dim, FIFTH};
  fe->attach_quadrature_rule(&qrule);

  std::vector<double> const & JxW = fe->get_JxW();
  std::vector<std::vector<double>> const & phi = fe->get_phi();
  std::vector<std::vector<RealGradient>> const & dphi = fe->get_dphi();
  // std::vector<Point> const & qpoint = fe->get_xyz();

  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

  std::vector<dof_id_type> dof_indices;

  for (auto const elem: mesh.active_local_element_ptr_range())
  {
    dof_map.dof_indices(elem, dof_indices);
    uint const ndofs = dof_indices.size();

    fe->reinit(elem);

    Ke.resize(ndofs, ndofs);
    Fe.resize(ndofs);

    for (uint qp = 0; qp < qrule.n_points(); qp++)
    {
      for (uint i = 0; i < ndofs; i++)
      {
        Fe(i) += JxW[qp] * 1.0 * phi[i][qp];

        for (uint j = 0; j < ndofs; j++)
        {
          Ke(i, j) = JxW[qp] * (dphi[j][qp] * dphi[i][qp]);
        }
      }
    }

    dof_map.heterogeneously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }
}
