import numpy as np
import matplotlib.pyplot as plt

from GwFlow import GwFlowSolver
from random_process import RandomProcess

class Model:
    def __init__(self, resolution, field_mean, field_stdev, mkl, lamb):
        
        '''
        This class is basically a wrapper around GwFlowSolver and RandomProcess.
        It has some functions that makes it simple to access important features.
        '''
        
        # Set up all the parameters needed for the solver and the ransom field.
        self.resolution = resolution
        self.field_mean = field_mean
        self.field_stdev = field_stdev
        self.mkl = mkl
        self.lamb = lamb
        
        # Initialise a solver.
        self.solver = GwFlowSolver(self.resolution, self.field_mean, self.field_stdev)
        self.x = self.solver.mesh.coordinates()[:,0]; self.y = self.solver.mesh.coordinates()[:,1]
        
        # initialise a random process given the solver mesh.
        self.random_process = RandomProcess( self.solver.mesh, self.mkl, self.lamb ) #changed
        
        # Compute the eigenpairs of the covariance matrix in the random process.
        self.random_process.compute_eigenpairs()

        # self.random_process.plot_eigenvalues()

        
        
    def solve(self, parameters = None):

        # self.random_process.plot_eigenvalues()
        
        # Solve the problem, given a vector of modes.
        self.random_process.generate(parameters)
        self.parameters = self.random_process.parameters
        self.solver.set_conductivity(self.random_process.random_field)
        # self.solver.set_conductivity(parameters)
        self.solver.solve()
        
        # self.solver.plot_mesh()
        # self.solver.plot_solution()
        
        
    def solve_df(self, parameters = None):
        
        # Solve the problem, given a vector of modes.
        self.random_process.generate(parameters)
        self.parameters = self.random_process.parameters
        self.solver.set_conductivity(self.random_process.random_field)
        self.solver.solve_df()
        
    def get_solution(self):
        return np.fromiter(map(self.solver.h, self.x, self.y), dtype=float)
        
    def get_data(self, datapoints):
        # Get data from a list of coordinates.
        return self.solver.get_data(datapoints)
    
    def get_data_mesh(self):
        # Get solution on mesh points and coordinates of noeds

        return self.solver.get_data_mesh()
        
    def get_outflow(self):
        return self.solver.get_outflow()
        
    def plot(self, limits = [0,0], transform_field = False):
        
        # This method plots both the random field and the solution.
        
        # First, contruct a random field, given the field parameters.
        if transform_field:
            random_field = np.exp(self.field_mean + self.field_stdev*self.random_process.random_field)
        else:
            random_field = self.field_mean + self.field_stdev*self.random_process.random_field
        
        # Set contour levels.
        if any(limits):
            contour_levels_field = np.linspace(limits[0], limits[1], 100)
        else:
            contour_levels_field = np.linspace(min(random_field), max(random_field), 100)
        
        # Then extract the solution from every node in the solver.
        solution = self.get_solution()
            
        # Set the contour levels.
        contour_levels_solution = np.linspace(min(solution), max(solution), 100)
        
        # Plot field and solution.
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize = (24, 9))
        
        axes[0].set_title('Transmissivity Field', fontdict = {'fontsize': 24})
        axes[0].tick_params(labelsize=16)
        field = axes[0].tricontourf(self.x, 
                                    self.y, 
                                    random_field, 
                                    levels = contour_levels_field, 
                                    cmap = 'plasma');  
        fig.colorbar(field, ax=axes[0])
        
        axes[1].set_title('Solution', fontdict = {'fontsize': 24})
        axes[1].tick_params(labelsize=16)
        solution = axes[1].tricontourf(self.x, 
                                       self.y, 
                                       solution, 
                                       levels = contour_levels_solution, 
                                       cmap = 'plasma');  
        fig.colorbar(solution, ax=axes[1])

    
    def variability_captured(self):
        
        variability = self.random_process.compute_variability()
        print(f"Variability captured with first 8 modes: {variability}")
        self.random_process.plot_variability()


    def plot_multiple_realizations_with_sensors(self, num_realizations=4, transform_field=True):
        """
        Plotta più realizzazioni del campo e della soluzione in un'unica figura 2x4,
        evidenziando i punti sensore.
        """

        x_data = y_data = np.array([0.056, 0.168, 0.280, 0.392, 0.500, 0.612, 0.725, 0.837, 0.949]) 
        sensor_coords = [(x, y) for x in x_data for y in y_data]

        fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(24, 12))

        for i in range(num_realizations):
            self.solve()
            x, y = self.x, self.y

            # Campo (trasformato o no)
            field = self.field_mean + self.field_stdev * self.random_process.random_field
            if transform_field:
                field = np.exp(field)

            solution = self.get_solution()

            levels_field = np.linspace(np.min(field), np.max(field), 100)
            levels_sol = np.linspace(np.min(solution), np.max(solution), 100)

            # Plot campo
            ax_field = axes[0, i]
            cf1 = ax_field.tricontourf(x, y, field, levels=levels_field, cmap='plasma')
            ax_field.set_title(f"Transmissivity field #{i+1}", fontsize=20)
            ax_field.set_aspect('equal')
            fig.colorbar(cf1, ax=ax_field)

            # Plot soluzione
            ax_sol = axes[1, i]
            cf2 = ax_sol.tricontourf(x, y, solution, levels=levels_sol, cmap='viridis')
            ax_sol.set_title(f"Solution #{i+1}", fontsize=20)
            ax_sol.set_aspect('equal')
            for xc, yc in sensor_coords:
                ax_sol.plot(xc, yc, 'ro', markersize=3)
            fig.colorbar(cf2, ax=ax_sol)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig("/root/shared/plot_4.png")
        plt.show()
