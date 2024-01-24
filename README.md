# Turbulence Boundary Layer

Build a structured boundary layer near walls for turbulent simulations.

## Task list

- [ ] define a suitable `struct`/`class` that stores information on the geometry and
  mesh features

- [ ] generate an initial unstructured mesh with `igl`/`cgal`

- [ ] feed the mesh to `libmesh`

- [ ] solve a Poisson problem using `libmesh` with fixed boundary conditions on the
  wall(s)

- [ ] extract isosurfaces of the solution using `igl`

- [ ] use the isosurfaces to define the boundary layer

- [ ] mesh the boundary layer with a structured mesh using `igl`/`cgal`

- [ ] combine the strcutured boundary layer mesh with the unstructured mesh on the
  outside

## External libraries

- `libmesh`

- `igl`

- `cgal` (optional)

### Dependencies

- `eigen`

- `openblas`

- `spack`

### Installation

Using `spack`: use the dedicated environment

```bash
spack env activate spack_env
spack concretize -f
spack install
```

After that, check that all required packages are available

```bash
$> spack find
==> In environment .../tbl/spack_env
==> Root specs
cgal@5.5.2 +eigen  eigen@3.4.0   libmesh@master +eigen+exodusii~mpi~petsc threads=openmp  openblas

==> Installed packages
...
```
