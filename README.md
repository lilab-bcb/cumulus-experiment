# Experiment and Benchmark on Cumulus and Pegasus

## Software Versions for Benchmark

We benchmark Pegasus, SCANPY, and Seurat with the following versions:

| Software | Version     | Release Date |
|:---------|:-----------:|:------------:|
| Pegasus  | 0.15.0      | 10/02/2019   |
| SCANPY   | 1.4.4.post1 | 07/29/2019   |
| Seurat   | 3.1.0       | 08/20/2019   |

Both Pegasus and SCANPY are run on Python 3.7.3. Seurat is run on R 3.6.1.

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

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/diffmap
(pegasus-env) root# python get_diffmap_figures.py
```

When finished, you'll find the figures generated in ``/output``.

## Experiment on Clustering Algorithms

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/clustering
(pegasus-env) root# python algorithm_compare.py
(pegasus-env) root# python spectral_clustering.py
```

When finished, you'll find the figures generated in ``/output``, and AMI results can be read from the screen output. Besides, you can check ``/experiment/pegasus.log`` for the execution time on each of the 4 clustering algorithms, while the time on spectral clustering can be read from the screen output.

## Experiment on Visualization Methods

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/visualization
(pegasus-env) root# python origin_vs_net.py
```

When finished, you'll find the figures generated in ``/output``, and for each visualization method, its kSIM regarding Louvain clustering labels can be seen from screen output. Besides, you can check ``/experiment/pegasus.log`` for the execution time on each visualization method.

## Benchmark on Analysis Tasks

All 3 packages are benchmarked on Bone Marrow dataset, while only Pegasus and SCANPY are benchmarked on Mouse Neuron, because Seurat fails at loading the count matrix step for this big dataset.

### Manton Bone Marrow Dataset

#### Pegasus

To benchmark Pegasus, in Pegasus environment, run the following command:

```
(pegasus-env) root# cd /experiment/overall
(pegasus-env) root# python run_pegasus.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/pegasus.log``.

#### SCANPY

To benchmark SCANPY, in SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/overall
(scanpy-env) root# python run_scanpy.py > scanpy.log
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/scanpy.log``.

#### Seurat

The benchmark on Seurat is a little bit complicated, as it's written in R. 

First, in Pegasus environment, run the following commands to generate a Seurat-compatible h5ad file on the dataset, and convert it into Seurat object format:

```
(pegasus-env) root# cd /experiment/overall
(pegasus-env) root# ./get_seurat_compatible.sh
(pegasus-env) root# Rscript convert_mantonbm_pegasus.R
```

When finished, you'll have a Seurat object in ``/experiment/overall/MantonBM_nonmix.RData`` for benchmarking steps starting from kNN. 

Now in SCANPY environment, run the following command:

```
(scanpy-env) root# Rscript run_seurat3.R
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/seurat.log``.

As Seurat fails in the Batch correction and Leiden clustering steps, I make them as two separate R scripts for users to try themselves.

Run the following command to benchmark on Batch correction using Seurat:

```
root# Rscript seurat_batch_correction.R
```

When terminated with failure, you'll find the time information in its log file ``/experiment/overall/seurat_batch_correction.log``, and error message from screen output.

In SCANPY environment, run the following command to benchmark on Leiden clustereing using Seurat:

```
(scanpy-env) root# Rscript seurat_leiden.R
```

When terminated with failure, you may find information in its log file ``/experiment/overall/seurat_leidenl.log`` and screen output.


### Mouse Neuron Dataset

#### Pegasus

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/overall
(pegasus-env) root# python run_pegasus_1m.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/1m_pegasus.log``.

#### SCANPY

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/overall
(scanpy-env) root# python run_scanpy_1m.py > 1m_scanpy.log
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/1m_scanpy.log``.

#### Seurat

Seurat fails at loading data step. Users can try the following commands in R environment, and check out the error message:

```
> library(Seurat)
> adata <- Read10X("/data/1M_neurons.h5")
```

## Benchmark on Workflows

The benchmark was performed on Google Cloud, with 32 CPUs under the default Haswell platform, and 120 GB memory. The dataset used is Manton Bone Marrow dataset. The analysis tasks performed are:

* Count matrix aggregation
* Highly variable features selection
* Batch correction
* PCA
* Find K nearest neighbors
* Louvain-like clustering
* tSNE-like visualization
* UMAP-like visualization
* Differential expression analysis, and cell type annotation.

Notice that Cumulus is the only one providing **Count matrix aggregation** feature.

### Cumulus

Benchmark on Cumulus is done by running jobs on Terra via Cumulus WDL workflows. And its overall execution time includes all Terra or Google Cloud specific preprocessing and postprocessing phases.

You need to set up your account and workspace on Terra first. 

Then change parameters ``"cumulus.input_file"`` and ``"cumulus.output_name"`` in ``/experiment/cloud/inputs_32.cpu.json`` to your own. Also, you need to specify your own workspace name in ``/experiment/cloud/run_cumulus.sh`` after ``-w`` flag.

After that, in Pegasus environment, run the following command to start the execution:

```
(pegasus-env) root# cd /experiment/cloud
(pegasus-env) root# ./run_cumulus.sh
```

Then check its progress and result using the URL shown on your screen.

### SCANPY

As SCANPY doesn't have a cloud-based interface, its benchmark is performed on a Google Cloud VM. Besides, since we used BBKNN for batch correction, SCANPY doesn't need to find kNN.

Similarly as benchmarks before, in SCANPY environment, run the following commands

```
(scanpy-env) root# cd /experiment/cloud
(scanpy-env) root# python run_scanpy.py > scanpy.log
```

When finished, you'll find execution time for each step in its log file ``/experiment/cloud/scanpy.log``.

### Seurat

Similarly as SCANPY, benchmark on Seurat is performed on a Google Cloud VM, and used as a single-server solution.

Besides, as Seurat would fail for batch correction when using 63 channels as the batches, we instead use batch correction with 8 donors being the batches. Moreover, the batch correction failed when using 10, 15, 20, and 32 threads via R _future_ package. So we simply used 2 threads just to make sure the batch correction step terminates successfully. Then for all the other steps, 32 threads are used whenever possible.

In SCANPY environment, run the following commands

```
(scanpy-env) root# cd /experiment/cloud
(scanpy-env) root# Rscript run_seurat_hvg_batch_correction.R
(scanpy-env) root# Rscript run_seurat_analysis.R
```

When finished, you'll find execution time for each step in log files ``/experiment/cloud/seurat_batch_correction.log`` and ``/experiment/cloud/seurat_analysis.log``.