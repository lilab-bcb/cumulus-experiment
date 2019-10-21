# Experiment and Benchmark on Cumulus and Pegasus

## Get Raw Data

## Start Docker Container

`
	$ docker run -it --rm -v /path-to-data:/data -v /path-to-output:/output cumulusprod:cumulus-experiment
`

where ``/path-to-data`` is the local directory in which the experiment data are stored, and ``/path-to-output`` is the local directory to which you want to set the experiment output.

## Generate Data for Experiment and Benchmark

First, enter the Conda environment of pegasus by:

`
	$ source activate pegasus-env
`

Then execute
`
	(pegasus-env) $ python generate_data_info.py MantonBM
`
to get all the necessary data for experiment on Manton Bone Marrow dataset.

Execute
`
	(pegasus-env) $ python generate_data_info.py 1M_neurons
`
to get all the necessary data for experiment on Mouse Neuron dataset.

## Experiment on Highly Variable Feature Selection

## Experiment on Batch Correction

## Experiment on k-Nearest-Neighbor

## Experiment on Diffusion Maps

## Experiment on Clustering Algorithms

## Experiment on Visualization Methods

## Benchmark on Analysis Tasks

## Benchmark on Workflows