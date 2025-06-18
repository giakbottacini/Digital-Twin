import os
import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from pyDOE2 import lhs as lhcs
from scipy.stats.distributions import norm as norm_dist
from sklearn.decomposition import TruncatedSVD
from tqdm import tqdm

# Environment Setup
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['KERAS_BACKEND'] = 'tensorflow'

# Import your custom model module (assuming correct structure)
from model import Model

def setup_environment():
    """
    Configure the random seed and return the setup parameters.
    """
    np.random.seed(123)
    # n_samples = 128000
    n_samples = 40000
    # resolutions = [(60, 60), (30, 30), (20,20), (15, 15)]
    resolutions = [(60,60)]
    field_mean, field_stdev, lamb_cov = 1, 1, 0.15
    # field_mean, field_stdev, lamb_cov = 1, 1, 0.25  
    mkl_values = [32] 
    # mkl_values = [128]   #Number of modes (dominant patterns) i generate (higher->higher variability captured)
    # x_data = y_data = np.array([0.1, 0.3, 0.5, 0.7, 0.9])  # 5 X 5 grid of sensors
    # x_data = y_data = np.array([0.083, 0.249, 0.415, 0.581, 0.747, 0.913])  # 6 X 6 grid of sensors
    # x_data = y_data = np.array([0.07, 0.211, 0.353, 0.495, 0.637, 0.78, 0.92])  # 7 X 7 grid of sensors
    x_data = y_data = np.array([0.056, 0.168, 0.280, 0.392, 0.500, 0.612, 0.725, 0.837, 0.949])  # 9 X 9 grid of sensors
    # x_data = y_data = np.array([0.25, 0.5, 0.75])  # 3 X 3 grid of sensors
    datapoints = np.array(list(product(x_data, y_data)))  #Sensors's coordinates
    
    return n_samples, resolutions, field_mean, field_stdev, lamb_cov, mkl_values, datapoints

def setup_random_process(solver_high, solver_low):
    """
    Synchronize the random processes between the higher and lower fidelity models
    by matching transmissivity fields across different resolutions.
    """
    coords_high = solver_high.solver.mesh.coordinates()
    coords_low = solver_low.solver.mesh.coordinates()
    
    structured_high = np.array(coords_high).view([('', coords_high.dtype)] * coords_high.shape[1])
    structured_low = np.array(coords_low).view([('', coords_low.dtype)] * coords_low.shape[1])
    
    bool_mask = np.in1d(structured_high, structured_low)
    solver_low.random_process.eigenvalues = solver_high.random_process.eigenvalues
    solver_low.random_process.eigenvectors = solver_high.random_process.eigenvectors[bool_mask]

def solver_data(solver, datapoints, x):
    """
    Run the solver for given input x and return the computed data.
    """
    # solver.variability_captured()  #To see variability captured by first modes

    solver.solve(x)
    # return solver.get_data(datapoints)
    return solver.get_data_mesh()
    

def generate_samples(mkl_values, n_samples):
    # n_regions = 9
    return norm_dist(loc=0, scale=1).ppf(lhcs(mkl_values[0], samples=n_samples))
    # return norm_dist(loc=0, scale=1).ppf(lhcs(n_regions, samples=n_samples))
    

def generate_solver_data(solvers, solver_key, samples, datapoints, path_prefix):
    """
    Generate data for a specific solver and save training/testing sets.
    
    :param solver_key: The key identifying the solver in the `solvers` dictionary.
    :param samples: previously generated samples.
    :param path_prefix: Directory path to save the generated data.
    """
    n_samples = samples.shape[0]
    solver = solvers[solver_key]

    solver.plot_multiple_realizations_with_sensors()
    
    # data = np.zeros((n_samples, len(datapoints)))
    # for i in tqdm(range(n_samples), desc=f"Processing {solver_key} samples"):
    #     data[i, :] = solver_data(solver, datapoints, samples[i, :])

    data = np.zeros((n_samples, 3721))
    k_field = np.zeros((n_samples, 3721))
    for i in tqdm(range(n_samples), desc=f"Processing {solver_key} samples"):
        _, _, data[i, :], _ = solver_data(solver, datapoints, samples[i, :])
        _, _, _, k_field[i, :] = solver_data(solver, datapoints, samples[i, :])

    coords, _, _, _= solver_data(solver, datapoints, samples[0, :])
    np.savetxt(f"{path_prefix}/coordinates.csv", coords, delimiter=",")

    _, cells, _, _= solver_data(solver, datapoints, samples[0, :])
    np.savetxt(f"{path_prefix}/cells.csv", cells, delimiter=",")

    # Genera gli indici mescolati
    indices = np.random.permutation(n_samples)

    # Applica il rimescolamento a samples e data
    samples_shuffled = samples[indices]
    data_shuffled = data[indices]
    k_field_shuffled = k_field[indices]
    
    # Split data into training and testing sets
    split_idx = int(0.8 * n_samples)
    # np.savetxt(f"{path_prefix}/X_train_{solver_key}_60_3.csv", samples[:split_idx], delimiter=",")
    # np.savetxt(f"{path_prefix}/y_train_{solver_key}_60_3.csv", data[:split_idx], delimiter=",")
    # np.savetxt(f"{path_prefix}/X_test_{solver_key}_60_3.csv", samples[split_idx:], delimiter=",")
    # np.savetxt(f"{path_prefix}/y_test_{solver_key}_60_3.csv", data[split_idx:], delimiter=",")

    np.savetxt(f"{path_prefix}/samples_train_{solver_key}_60.csv", samples_shuffled[:split_idx], delimiter=",")
    np.savetxt(f"{path_prefix}/sensorsdata_train_{solver_key}_60.csv", data_shuffled[:split_idx], delimiter=",")
    np.savetxt(f"{path_prefix}/samples_test_{solver_key}_60.csv", samples_shuffled[split_idx:], delimiter=",")
    np.savetxt(f"{path_prefix}/sensorsdata_test_{solver_key}_60.csv", data_shuffled[split_idx:], delimiter=",")

    np.savetxt(f"{path_prefix}/kfield_train_{solver_key}_60.csv", k_field_shuffled[:split_idx], delimiter=",")
    np.savetxt(f"{path_prefix}/kfield_test_{solver_key}_60.csv", k_field_shuffled[split_idx:], delimiter=",")

def project_to_pod_basis(coarse_data, n_components=45):
    """
    Perform Singular Value Decomposition (SVD) and project data onto the POD basis.
    
    :param coarse_data: Data matrix to be projected.
    :param n_components: Number of components to retain in the projection.
    :return: The retained POD basis.
    """
    y_t = coarse_data.T
    U, S, Vh = np.linalg.svd(y_t, full_matrices=False)

    return U[:, :n_components]

def project_and_save_pod(solvers, solver_key, samples, datapoints, path_prefix):
    """
    Project solver data onto POD basis and save the results.
    
    :param solver_key: The key identifying the solver.
    :param samples: Samples.
    :param path_prefix: Directory path to save the results.
    """
    n_samples = samples.shape[0]
    X_train = np.loadtxt(f"{path_prefix}/X_train_{solver_key}_60_3.csv", delimiter=",")
    coarse_sol_train = np.array([solver_data(solvers[solver_key], datapoints, X_train[i, :]) for i in tqdm(range(int(n_samples * 0.9)), desc=f"Processing {solver_key} samples")])
    
    basis = project_to_pod_basis(coarse_sol_train)
    X_train_pod = coarse_sol_train @ basis
    
    X_test = np.loadtxt(f"{path_prefix}/X_test_{solver_key}_60_3.csv", delimiter=",")
    coarse_sol_test = np.array([solver_data(solvers[solver_key], datapoints, X_test[i, :]) for i in tqdm(range(int(n_samples * 0.1)), desc=f"Processing {solver_key} samples")])
    X_test_pod = coarse_sol_test @ basis
    
    # Save projected data and basis
    np.savetxt(f"{path_prefix}/X_train_{solver_key}_pod_60_3.csv", X_train_pod, delimiter=",")
    np.savetxt(f"{path_prefix}/X_test_{solver_key}_pod_60_3.csv", X_test_pod, delimiter=",")
    np.savetxt(f"{path_prefix}/POD_basis_{solver_key}_60_3.csv", basis, delimiter=",")

def print_simulation_parameters(n_samples, resolutions, field_mean, field_stdev, lamb_cov, mkl_values):
    """
    Print the simulation parameters to the screen for logging purposes.
    
    :param n_samples: Number of samples used in the simulation.
    :param resolutions: List of resolutions used in the simulation.
    :param field_mean: Mean value of the field.
    :param field_stdev: Standard deviation of the field.
    :param lamb_cov: Covariance of the lambda parameter.
    :param mkl_values: List of MKL values for the solvers.
    """
    print("Simulation Parameters")
    print("=====================")
    print(f"Number of samples: {n_samples}")
    print(f"Resolutions h_i: {resolutions}")
    print(f"Field Mean: {field_mean}")
    print(f"Field Standard Deviation: {field_stdev}")
    print(f"Lambda Covariance: {lamb_cov}")
    print(f"MKL Values: {mkl_values}")
    print("=====================")



def main():

    print("\nSetup equation solvers \n")

    # Setup environment and solvers
    n_samples, resolutions, field_mean, field_stdev, lamb_cov, mkl_values, datapoints = setup_environment()

    print_simulation_parameters(n_samples, resolutions, field_mean, field_stdev, lamb_cov, mkl_values)

    solvers = {
        "h1": Model(resolutions[0], field_mean, field_stdev, mkl_values[0], lamb_cov)
        # "h2": Model(resolutions[1], field_mean, field_stdev, mkl_values[1], lamb_cov),
        # "h3": Model(resolutions[2], field_mean, field_stdev, mkl_values[2], lamb_cov),
        # "h4": Model(resolutions[3], field_mean, field_stdev, mkl_values[2], lamb_cov)
    }    

    # Setup random processes between solvers
    # setup_random_process(solvers["h1"], solvers["h2"])
    # setup_random_process(solvers["h1"], solvers["h3"])
    # setup_random_process(solvers["h1"], solvers["h4"])

    print("\nGenerate solver data \n")

    # Generate data for all solvers
    samples = generate_samples(mkl_values, n_samples)

    # for key in ["h1", "h2", "h3", "h4"]:
    #     generate_solver_data(solvers, key, samples, datapoints, "/data")

    save_path = '/root/shared/data'

    generate_solver_data(solvers, "h1", samples, datapoints, save_path)

    # print("\nProject solution and save POD \n")
    # datapoints = np.array(list(product(np.linspace(0.0, 1.0, 10), repeat=2)))
    
    # # Perform POD projection and save for all solvers
    # for key in ["h1", "h2", "h3", "h4"]:
    #     project_and_save_pod(solvers, key, samples, datapoints, "/data")

    # project_and_save_pod(solvers, "h1", samples, datapoints, "/data")

if __name__ == "__main__":
    main()
