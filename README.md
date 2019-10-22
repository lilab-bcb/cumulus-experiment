# Experiment and Benchmark on Cumulus and Pegasus

## Software Versions for Benchmark

We benchmark Pegasus, SCANPY, and Seurat with the following versions:

| Software | Version     | Release Date |
|:---------|:-----------:|:------------:|
| Pegasus  | 0.15.0      | 10/02/2019   |
| SCANPY   | 1.4.4.post1 | 07/29/2019   |
| Seurat   | 3.1.0       | 08/20/2019   |

## Get Raw Data

## Start Docker Container

```
$ docker run -it --rm -v /path-to-data:/data -v /path-to-output:/output cumulusprod:cumulus-experiment
```

where ``/path-to-data`` is the local directory in which the experiment data are stored, and ``/path-to-output`` is the local directory to which you want to set the experiment output.

Notice that there are 2 conda environments already installed: ``pegasus-env`` and ``scanpy-env``. You can activate/deactivate either of them with the following commands (taking ``pegasus-env`` as the example):

```
root# source activate pegasus-env
(pegasus-env) root# conda deactivate
$
```

## Generate Data for Experiment and Benchmark

First, enter the Conda environment of pegasus by:

```
root# source activate pegasus-env
```

Then execute

```
(pegasus-env) root# python generate_data_info.py MantonBM
```

to get all the necessary data for experiment on Manton Bone Marrow dataset.

Execute

```
(pegasus-env) root# python generate_data_info.py 1M_neurons
```

to get all the necessary data for experiment on Mouse Neuron dataset.

## Experiment on Highly Variable Feature Selection

First, enter the experiment subfolder:

```
(pegasus-env) root# cd /experiment/highly_variable_features
```

Then execute the experiment:

```
(pegasus-env) root# python gen_result.py
```

When finished, you'll find the figures generated in ``/output``, and CSV files containing markers in the current folder.

## Experiment on Batch Correction

### Ground Truth Cell Types

As MNN and Seurat CCA both fail for the whole Bone Marrow dataset, we use a subset of it, by selecting one channel per donor, for the batch correction benchmark.

First, we need a list of cell types as the ground truth to measure the performance of different batch correction methods. We use the information from the clustering on the dataset not batch-corrected for this purpose.

Enter the folder for generating ground truth information, and execute the script:

```
(pegasus-env) root# cd /experiment/batch_correction/ground
(pegasus-env) root# python gen_ground_h5ad.py
(pegasus-env) root# python gen_celltype.py
```

This will generate a file containing ground truth of cell types (``ground_cell_types.txt``) in folder ``/experiment/batch_correction``.

### Benchmark on Baseline Method

The baseline is running Pegasus clustering without batch correction. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/batch_correction/baseline
(pegasus-env) root# python run_baseline.py
```

When finished, you'll have a result file ``baseline_result.h5ad`` in the current folder.

### Benchmark on Pegasus Batch Correction Method

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/batch_correction/pegasus
(pegasus-env) root# python run_pegasus_batch_correct.py
```

When finished, you'll have a result file ``pegasus_corrected.h5ad``, along with its log file ``pegasus_correct.log`` in the current folder. To calculate its batch correction time, in ``pegasus_correct.log``, add up the time spent on **Estimation on feature statistics per channel** and **Batch correction** steps.

### Benchmark on ComBat

In SCANPY environment, run the following commands:

```
(scanpy-env) root# cd /experiment/batch_correction/combat
(scanpy-env) root# python scanpy_combat.py
```

When finished, you'll have a result file ``scanpy_combat_corrected.h5ad`` in the corrent folder, and you can read the time spent on ComBat from the screen output.

### Benchmark on MNN

In SCANPY environment, run the following commands:

```
(scanpy-env) root# cd /experiment/batch_correction/mnn
(scanpy-env) root# python scanpy_mnn.py
```

When finished, you'll have a result file ``scanpy_mnn_corrected.h5ad`` in the current folder, and you can read the time spent on MNN from the screen output.

### Benchmark on BBKNN

In SCANPY environment, run the following commands:

```
(scanpy-env) root# cd /experiment/batch_correction/bbknn
(scanpy-env) root# python scanpy_bbknn.py
```

When finished, you'll have a result file ``scanpy_bbknn_corrected.h5ad`` in the current folder, and you can read the time spent on BBKNN from the screen output.

### Benchmark on Seurat CCA

Running Seurat batch correction doesn't depend on Python environment. Run the following commands:

```
root# cd /experiment/batch_correction/seurat
root# Rscript seurat_cca.R
```

When finished, you'll have 3 result files: ``matrix.mtx`` for count matrix, ``barcodes.txt`` for cell barcode names, ``genes.txt`` for gene names in the current folder. Besides, there is also a log file ``seurat_cca.log`` in the folder. To get its batch correction time, add up the time spent on **Finding Anchors** and **Integration** together.

### Measure Performance of All Batch Correction Methods

We calculate two measures for each batch correction method result: **kSIM** and **kBET** accept rates. In Pegasus environment, type the following commands:

```
(pegasus-env) root# python measure_result.py baseline
(pegasus-env) root# python measure_result.py pegasus
(pegasus-env) root# python measure_result.py seurat
(pegasus-env) root# python measure_result.py combat
(pegasus-env) root# python measure_result.py bbknn
(pegasus-env) root# python measure_result.py mnn
```

For each of the commands above, you'll see its kSIM and kBET accept rates from the screen output. Besides, the method's result UMAP plot will be generated in ``/output``. When finished, don't forget to update its measures in ``correction_benchmark.txt`` for later plot.

After executing all these commands, run

```
(pegasus-env) root# python measure_result.py plot
```

to generate the measurement plot on batch correction methods in ``/output``.

## Experiment on k-Nearest-Neighbor

### Generate Ground Truth kNN

We use the kNN result by brute force algorithm in _scikit-learn_ as the ground truth. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py brute
```

When finished, you'll have a result file ``baseline_indices.npy`` in the current folder.

### Benchmark on Pegasus kNN Method

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py pegasus
```

When finished, you'll have a result file ``pegasus_indices.npy`` in the current folder, and you can read the time spent on kNN from screen output.

### Benchmark on SCANPY kNN Method

SCANPY uses kNN from _umap-learn_ package. So we directly benchmark this function. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py scanpy
```

When finished, you'll have a result file ``scanpy_indices.npy`` in the current folder, and you can read the time spent on kNN from screen output.

### Benchmark on Seurat kNN Method

Seurat's kNN has two methods: ``nn2`` from _RANN_ package, default method but time-consuming; ``AnnoyNN`` from _RcppAnnoy_ package, not default but more efficient. We decide to choose ``AnnoyNN`` for this benchmark. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# Rscript seurat_knn.R
(pegasus-env) root# python knn_comp.py seurat
```

When finished, you'll have a result file ``seurat_indices_annoy.txt`` in the current folder, and you can read the time spent on kNN from screen output.

## Measure kNN Methods

Given that all the kNN resuls of the 3 methods are calculated, we are ready to measure the recall of them, and generate their time and recall plots. 

First, update ``time_stats.txt`` in ``/experiment/knn_comparison`` by the time information you saw in benchmarks above. Then in Pegasus environment, type:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py plot
```

This will generate the corresponding figures in ``/output``.

## Experiment on Diffusion Maps

## Experiment on Clustering Algorithms

## Experiment on Visualization Methods

## Benchmark on Analysis Tasks

## Benchmark on Workflows