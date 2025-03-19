import numpy as np
from numpy.linalg import inv, det

import matplotlib.pyplot as plt

from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh 
from scipy.spatial import distance_matrix


class RandomProcess:
    def __init__(self, dolfin_mesh, mkl, lamb):
        
        '''
        This class sets up a random process on a grid and generates
        a realisation of the process, given parameters or a random vector.
        '''
        
        # Internalise the grid and set number of vertices.
        self.mesh = dolfin_mesh
        self.n_points = self.mesh.num_vertices()
        
        # Save the coordinates in some vectors.
        self.x = self.mesh.coordinates()[:,0]; self.y = self.mesh.coordinates()[:,1]
        
        # Set some random field parameters.
        self.mkl = mkl
        self.lamb = lamb
        
        # Create a snazzy distance-matrix for rapid computation of the covariance matrix.
        dist = distance_matrix(self.mesh.coordinates(), self.mesh.coordinates())
        
        # Compute the covariance between all points in the space.
        self.cov =  np.exp(-0.5*dist**2/self.lamb**2)
    
    def plot_covariance_matrix(self):
        
        # Plot the covariance matrix.
        plt.figure(figsize = (10,8)); plt.imshow(self.cov, cmap = 'binary'); plt.colorbar(); plt.show()
    
    def compute_eigenpairs(self):
        
        # Find eigenvalues and eigenvectors using Arnoldi iteration.
        #eigvals, eigvecs = np.linalg.eigh(self.cov)
        eigvals, eigvecs = eigh(self.cov, eigvals = (self.n_points - self.mkl, self.n_points - 1))
        #eigvals, eigvecs = eigsh(self.cov, self.mkl, which = 'LM')
        
        
        order = np.flip(np.argsort(eigvals))
        self.eigenvalues = eigvals[order]
        self.eigenvectors = eigvecs[:,order]
      
    def generate(self, parameters = None):
        
        # Generate a random field, see
        # Scarth, C., Adhikari, S., Cabral, P. H., Silva, G. H. C., & Prado, A. P. do. (2019). 
        # Random field simulation over curved surfaces: Applications to computational structural mechanics. 
        # Computer Methods in Applied Mechanics and Engineering, 345, 283–301. https://doi.org/10.1016/j.cma.2018.10.026
        
        if parameters is None:
            self.parameters = np.random.normal(size=self.mkl)
            
        else:
            self.parameters = np.array(parameters).flatten()
        
        self.random_field = np.linalg.multi_dot((self.eigenvectors, 
                                                 np.sqrt(np.diag(self.eigenvalues)), 
                                                 self.parameters))
        

    def plot(self, lognormal = False):
        
        # Plot the random field.
        if lognormal:
            random_field = np.exp(self.random_field)
            contour_levels = np.linspace(min(random_field), max(random_field), 20)
        else:
            random_field = self.random_field
            contour_levels = np.linspace(min(random_field), max(random_field), 20)

        plt.figure(figsize = (12,10))
        plt.tricontourf(self.x, self.y, random_field, levels = contour_levels, cmap = 'magma'); 
        plt.colorbar()
        plt.show()

    def plot_same_scale(self, lognormal=False, vmin=None, vmax=None, compare=True):
        # Plot the random field.
        if lognormal:
            random_field = np.exp(self.random_field)
        else:
            random_field = self.random_field

        contour_levels = np.linspace(vmin, vmax, 50)  # Use vmin and vmax for consistent levels
        if not compare:
            plt.figure(figsize=(12, 10))

        plt.tricontourf(self.x, self.y, random_field, levels=contour_levels, cmap='magma', vmin=vmin, vmax=vmax)
        plt.colorbar()
        if not compare:
            plt.show()


    def plot_eigenvalues(self, save_path="/root/shared/eigenvalues.png"):
        plt.figure(figsize=(8, 6))
        plt.plot(self.eigenvalues, marker='o', linestyle='-', color='b')

        # Aggiungi i valori numerici dei primi 6 autovalori con gestione delle sovrapposizioni
        last_y = None  # Per tracciare l'ultima posizione Y usata
        for i in range(min(6, len(self.eigenvalues))):
            y_value = self.eigenvalues[i]
            
            # Se l'etichetta è troppo vicina alla precedente, spostala più in alto
            if last_y is not None and abs(y_value - last_y) < 0.1 * y_value:
                y_value += 0.05 * y_value
            
            plt.text(i, y_value, f"{self.eigenvalues[i]:.2e}",
                    fontsize=9, ha='center', va='bottom', color='black')
            
            last_y = y_value  # Aggiorna la posizione Y dell'ultima etichetta

        plt.xlabel('Index')
        plt.ylabel('Eigenvalue')
        plt.title('Eigenvalues of the Covariance Matrix')
        plt.grid(True)

        # Salva il grafico invece di mostrarlo
        plt.savefig(save_path, bbox_inches='tight')
        # plt.close()  # Chiudi la figura per liberare memoria


