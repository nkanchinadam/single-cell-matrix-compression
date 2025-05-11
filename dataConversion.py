import os
import gzip
from scipy.io import mmread
from scipy.sparse import csc_matrix, csr_matrix, save_npz
import loompy
import pandas as pd

cellTypes = ["human_kidney", "human_blood", "human_skin"]
folders = ["data/humanKidneyCells", "data/humanBloodCells", "data/humanSkinCells"]

for cell_type, folder in zip(cellTypes, folders):
    input_file = os.path.join(folder, "matrix.mtx.gz")
    
    # Read the market matrix from the gzipped file
    with gzip.open(input_file, 'rt') as f:
        matrix = mmread(f)

    # Convert to dense format
    matrix_dense = matrix.todense()
    dense_df = pd.DataFrame(matrix_dense)

    # Convert to and save CSC format
    matrix_csc = csc_matrix(matrix)
    csc_file_path = os.path.join(folder, f'{cell_type}_matrix_csc.npz')
    save_npz(csc_file_path, matrix_csc)

    # Convert to and save CSR format
    matrix_csr = csr_matrix(matrix)
    csr_file_path = os.path.join(folder, f'{cell_type}_matrix_csr.npz')
    save_npz(csr_file_path, matrix_csr)

    # Convert to and save Loom format
    data = matrix_csc.toarray()  # Loom requires dense data
    row_attrs = {"Gene": [f"Gene_{i}" for i in range(data.shape[0])]}  
    col_attrs = {"Cell": [f"Cell_{j}" for j in range(data.shape[1])]} 
    loom_file_path = os.path.join(folder, f'{cell_type}_matrix.loom')
    loompy.create(loom_file_path, data, row_attrs, col_attrs)