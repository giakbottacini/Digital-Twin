from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, lil_matrix
from scipy.sparse.linalg import spsolve

class GwFlowSolver:
    def __init__(self, resolution, field_mean, field_stdev):
        '''
        This class solves the steady state groundwater flow equation (Darcy)
        on a unit square with some simple boundary conditions.
        '''
        # Internatise conductivity field parameters.
        self.field_mean  = field_mean
        self.field_stdev = field_stdev

        # To suppress the output of the model
        set_log_level(LogLevel.ERROR)
        # to restore the output
        # set_log_level(LogLevel.PROGRESS)

        # Head at inflow and outflow.
        self.h_in = 1
        self.h_out = 0

        # Zero flow through boundaries.
        self.q_0 = Constant(0.0)
        
        # Create mesh and define function space
        self.nx = resolution[0]
        self.ny = resolution[1]
        self.mesh = UnitSquareMesh(self.nx, self.ny)

        self.V = FunctionSpace(self.mesh, 'CG', 1)
        self.n = FacetNormal(self.mesh)
        self.d2v = dof_to_vertex_map(self.V)
        
        # Define variational problem
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)
        
        # Define the subdomains.
        sub_domains = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        sub_domains.set_all(0)

        # Sub domain for no-flow (mark whole boundary, inflow and outflow will later be overwritten)
        class Noflow(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary

        # Sub domain for inflow (left)
        class Inflow(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] < DOLFIN_EPS

        # Sub domain for outflow (right)
        class Outflow(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] > 1.0 - DOLFIN_EPS

        # Create instances of the boundaries and mark them with some numbers.
        noflow = Noflow()
        noflow.mark(sub_domains, 0)

        inflow = Inflow()
        inflow.mark(sub_domains, 1)

        outflow = Outflow()
        outflow.mark(sub_domains, 2)

        # Impose the subdomain numbering on external boundary ds.
        self.ds = ds(subdomain_data=sub_domains)

        # Set Dirichlet BC's
        bc_left = DirichletBC(self.V, self.h_in, inflow)
        bc_right = DirichletBC(self.V, self.h_out, outflow)
        self.bcs = [bc_left, bc_right]
        
    def plot_mesh(self):
        
        # This method plots the mesh
        plt.figure(figsize = (10,10))
        plot(self.mesh)
        plt.show()
        
    def set_conductivity(self, random_field = False):
        
        # Set the conductivity
        if np.any(random_field):
            # Make it exponential
            self.conductivity = np.exp(self.field_mean + self.field_stdev*random_field)

        # If no field is given, just set the flow-field to the mean.
        else:
            self.conductivity = np.exp(self.field_mean*np.ones(self.mesh.coordinates().shape[0]))
            
        # Map the random field vector to the domain.
        self.K = Function(self.V)
        self.K.vector().set_local(self.conductivity[self.d2v])

    

    # def set_conductivity(self, parameters, num_zones=9):
    #     """
    #     Imposta il campo di conduttività suddividendo il dominio quadrato (0,1)x(0,1) in 9 zone.
    #     Ogni zona avrà un valore di conduttività aleatorio.

    #     :param random_field: opzionale, campo aleatorio generato (se fornito)
    #     :param num_zones: numero di zone in cui suddividere il dominio (default è 9, che crea una griglia 3x3)
    #     """
    #     # Prendi il numero di punti sulla mesh
    #     num_points = self.mesh.coordinates().shape[0]

    #     # Numero di righe e colonne per la suddivisione 3x3
    #     rows = 3
    #     cols = 3
        
    #     # Ottieni la dimensione della mesh in X e Y
    #     x_min, x_max = self.mesh.coordinates()[:, 0].min(), self.mesh.coordinates()[:, 0].max()
    #     y_min, y_max = self.mesh.coordinates()[:, 1].min(), self.mesh.coordinates()[:, 1].max()
        
    #     # Calcola la dimensione di ogni zona
    #     zone_size_x = (x_max - x_min) / cols
    #     zone_size_y = (y_max - y_min) / rows

    #     self.conductivity = np.zeros(num_points)

    #     # Crea un valore di conduttività casuale per ciascuna delle 9 zone
    #     for i in range(rows):
    #         for j in range(cols):
    #             # Definisci il valore di conduttività per questa zona
    #             zone_value = np.exp(self.field_mean + self.field_stdev*parameters[j+3*i])
                
    #             # Definisci i limiti della zona nel dominio
    #             x_min_zone = x_min + j * zone_size_x
    #             x_max_zone = x_min + (j + 1) * zone_size_x
    #             y_min_zone = y_min + i * zone_size_y
    #             y_max_zone = y_min + (i + 1) * zone_size_y

    #             # print(f"x_min_zone: {x_min_zone}")
    #             # print(f"x_max_zone: {x_max_zone}")
    #             # print(f"y_min_zone: {y_min_zone}")
    #             # print(f"y_max_zone: {y_max_zone}")
                
    #             # Assegna il valore di conduttività a tutti i punti che appartengono a questa zona
    #             # Per ogni punto nella mesh, verifica se rientra nell'area della zona
    #             for idx, coord in enumerate(self.mesh.coordinates()):
    #                 x, y = coord
    #                 if x_min_zone <= x < x_max_zone and y_min_zone <= y < y_max_zone:
    #                     self.conductivity[idx] = zone_value
    #                     # print(f"Zone: {j+3*i}, Index:{idx}, Coord:{x,y}, Valore:{zone_value}")

    #     # Mappa il campo di conduttività sul dominio della mesh
    #     self.K = Function(self.V)
    #     self.K.vector().set_local(self.conductivity[self.d2v])

    
    def solve(self):
        # info(LinearVariationalSolver.default_parameters(), True)
        # Solve the variational problem
        F = inner(grad(self.v), self.K*grad(self.u))*dx - self.v*self.q_0*self.ds(0)
        a, L = lhs(F), rhs(F)
        self.h = Function(self.V)
        # Define the variational problem and solver
        problem = LinearVariationalProblem(a, L, self.h, self.bcs)
        solver = LinearVariationalSolver(problem)

        # Set solver parameters
        prm = solver.parameters
        prm['linear_solver'] = 'gmres'
        prm['preconditioner'] = 'ilu'
        prm['krylov_solver']['absolute_tolerance'] = 1e-12
        prm['krylov_solver']['relative_tolerance'] = 1e-12
        prm['krylov_solver']['maximum_iterations'] = 1000000

        # Solve the variational problem
        solver.solve()
        # solve(a == L, self.h, self.bcs, 
        #       solver_parameters={"linear_solver": "lu"},
        #       form_compiler_parameters={"optimize": True})
    
    def solve_df(self):
        dx, dy = 1.0 / (self.nx), 1.0 / (self.ny)
        
        conductivity_flat = self.K.vector().get_local()

        # Map the conductivity to the finite element vertices
        conductivity_at_vertices = conductivity_flat[self.d2v]

        # Reshape to match the finite difference grid layout
        conductivity = conductivity_at_vertices.reshape((self.ny + 1, self.nx + 1))

        # Number of grid points
        num_points = (self.nx+1) * (self.ny+1)
 
        # Create the coefficient matrix A and the right-hand side vector b
        A = lil_matrix((num_points, num_points))  # Sparse matrix
        b = np.zeros(num_points)  # RHS vector

        # Map 2D indices to 1D
        def idx(i, j):
            return j * (self.nx+1) + i

        # Fill the matrix A and vector b
        for j in range(self.ny+1):
            for i in range(self.nx+1):
                row = idx(i, j)

                if i == 0:  # Left boundary (Dirichlet)
                    A[row, row] = 1
                    b[row] = 1.0  # h_in
                elif i == self.nx:  # Right boundary (Dirichlet)
                    A[row, row] = 1
                    b[row] = 0.0  # h_out
                elif j == 0 or j == self.ny:  # Top and bottom boundaries (Neumann - no flow)
                    A[row, row] = 1  # Simplification for no-flow boundaries
                    b[row] = 0.0
                else:  # Interior points
                    Kx_left = (conductivity[j, i - 1] + conductivity[j, i]) / 2
                    Kx_right = (conductivity[j, i + 1] + conductivity[j, i]) / 2
                    Ky_down = (conductivity[j - 1, i] + conductivity[j, i]) / 2
                    Ky_up = (conductivity[j + 1, i] + conductivity[j, i]) / 2

                    A[row, idx(i - 1, j)] = -Kx_left / dx**2
                    A[row, idx(i + 1, j)] = -Kx_right / dx**2
                    A[row, idx(i, j - 1)] = -Ky_down / dy**2
                    A[row, idx(i, j + 1)] = -Ky_up / dy**2
                    A[row, row] = (Kx_left + Kx_right) / dx**2 + (Ky_down + Ky_up) / dy**2
                    

        # Solve the linear system
        self.h = Function(self.V)
        
        h_flat = spsolve(A.tocsr(), b)

        # Ensure `h_flat` is in the correct order (map grid points to vertices)
        h_vertices = np.zeros_like(self.h.vector().get_local())
        h_vertices[self.d2v] = h_flat  # Map finite difference solution to FEM vertices

        # Assign values to the finite element function
        self.h.vector().set_local(h_vertices)
        self.h.vector().apply("insert")


    def compute_flow(self):
        
        self.Q = VectorFunctionSpace(self.mesh, "CG", 1)
        self.q = project(-self.K*grad(self.h), self.Q)
        
    
    def get_data(self, datapoints):
        
        # Return data from a set of points.
        self.data = np.zeros(len(datapoints))
        for i, datapoint in enumerate(datapoints):
            self.data[i] = self.h(datapoint[0], datapoint[1])
        return self.data

    # from dolfin import vertex_to_dof_map

    # def get_data_mesh(self):
    #     """
    #     Return the solution at all mesh vertices, ordered by vertex index.
    #     """
    #     v2d = vertex_to_dof_map(self.V)
    #     coords = self.mesh.coordinates()               # Vertices in mesh order
    #     h_vec = self.h.vector().get_local()            # Values at DOFs
    #     data_on_vertices = h_vec[v2d]                  # Values reordered to match vertex indices
    #     cells = self.mesh.cells()                      # Cell definitions (by vertex indices)
    #     return coords, cells, data_on_vertices
    

    from dolfin import vertex_to_dof_map, dof_to_vertex_map

    def get_data_mesh(self):
        """
        Return the solution and conductivity at all mesh vertices, ordered by vertex index.
        """
        # Mapping DOFs to vertices
        v2d = vertex_to_dof_map(self.V)

        # Coordinates of mesh vertices
        coords = self.mesh.coordinates()               # (num_vertices, 2)

        # --- Solution values ---
        h_vec = self.h.vector().get_local()            # Solution values at DOFs
        sol_on_vertices = h_vec[v2d]                     # Reordered to match mesh vertex indices

        # --- Conductivity values ---
        k_vec = self.K.vector().get_local()            # Conductivity values at DOFs
        k_on_vertices = k_vec[v2d]                     # Reordered to match mesh vertex indices

        # --- Mesh cell definitions (vertex indices) ---
        cells = self.mesh.cells()                      # Cell connectivity

        return coords, cells, sol_on_vertices, k_on_vertices




    def get_outflow(self):
        return assemble(dot(-self.K*grad(self.h), self.n)*self.ds(2)) # This method works without computing the flow first.

    def plot_solution(self):
        
        # Plot the solution.
        plt.figure(figsize = (12,10))
        p = plot(self.h, cmap = 'magma'); plt.colorbar(p); plt.show()

        plt.colorbar(p)
        plt.title("Soluzione h")

        # Salva l'immagine
        plt.savefig("soluzione.png", dpi=300, bbox_inches='tight')  # puoi cambiare nome o formato

        # Mostra l'immagine e blocca finché non viene chiusa
        plt.show(block=True)
        plt.close()
