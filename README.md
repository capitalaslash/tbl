# Turbulence Boundary Layer

Build a structured boundary layer near walls for turbulent simulations.

Steps:

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

