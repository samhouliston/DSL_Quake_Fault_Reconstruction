# Earthquake Fault Network Reconstruction

## Contributions

This repository contains the following contributions:
* Python Implementation of the probabilistic agglomerative fault network clustering algorithm introduced by Kamer et al. (2020) in ```benchmark/src```
* Planar Agglomerative Clustering favoring Planarity (PACP), a version of the algorithm that favors planar clusters via introduction of new hyperparameters
* A synthetic data generator in ```synthetics/synthetic_data_generator.py```

## Running the Kamer algorithm

Change into the correct directory
```
cd benchmark/src
```

Run the algorithm from a python script/notebook on the 3D coordinates ```X``` of the events in your earthquake catalog
```python
from main import run_fault_reconstruction
kernels, labels = run_fault_reconstruction(X, min_sz_cluster=4)
```

## Running PACP

Change into the correct directory
```
cd benchmark/src
```

Run the algorithm from a python script/notebook on the 3D coordinates ```X``` of the events in your earthquake catalog
```python
from main import run_fault_reconstruction
kernels, labels = run_fault_reconstruction(X, min_sz_cluster=10, margin_scale=20, refit_kernels=True)
```
we recommend ```min_sz_cluster >= 10``` for good performance.


## Generating Synthetic Data

### Overview
The script `generate_synthetic_dataset.py` in the `synthetics` folder generates a dataset of a fault network comprising a point-cloud and corresponding cluster labels. 

### Execution
The script can be executed via the command line using the following syntax:

```sh
python generate_synthetic_dataset.py <arguments>
```

The script accepts the following arguments to specify the parameters for dataset generation:

- `-o`, `--output`: Path to the output file.
- `-v`, `--verbose`: Enable verbose mode.
- `--visualise`: Visualise the generated dataset.
- `--n_simple_faults`: Number of simple faults. Default value: 1.
- `--n_bent_faults`: Number of bent faults. Default value: 0.
- `--n_cross_faults`: Number of cross faults. Default value: 0.
- `--n_Y_faults`: Number of Y faults. Default value: 0.
- `--n_parallel_faults`: Number of parallel faults. Default value: 0.
- `--n_structures`: Number of structures. Default value: 0

### Example Command
Here is an example of how to run the script with some arguments:

```sh
python generate_synthetic_dataset.py -o output_path -v --n_simple_faults 3 --n_parallel_faults 3 --visualise
```

### Example Outputs

Example Specific Faults    |  Example Synthetic Dataset
:-------------------------:|:-------------------------:
![image](https://github.com/samhouliston/DSL_Quake_Fault_Reconstruction/assets/93288237/1adcd712-5bd8-46b2-b9be-fa631753af3e)  |  ![image](https://github.com/samhouliston/DSL_Quake_Fault_Reconstruction/assets/93288237/0707998d-b03a-47f9-af3d-f282f6d96f7c)


