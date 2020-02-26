# Experiment and Benchmark on Cumulus and Pegasus

This repository stores the code for experiment and benchmark in our BioRxiv paper "[Cumulus: a cloud-based data analysis framework for large-scale single-cell and single-nucleus RNA-seq](https://www.biorxiv.org/content/10.1101/823682v1.abstract)".

To try it out with our software, please pull our [cumulusprod/cumulus-experiment](https://hub.docker.com/repository/docker/cumulusprod/cumulus-experiment) docker image, and run it as a docker container. 

Sections below give instructions on each benchmark, as well as other useful instructions and information.

## Software Versions for Benchmark

The benchmark platform was a single server with 28 CPUs and Ubuntu Linux 18.04 OS. There is also a benchmark on cloud, which used Google Cloud platform, with detailed settings in Section [Benchmark on Workflows](#benchmark-on-workflows).

We benchmark Pegasus, SCANPY, and Seurat with the following versions:

| Software | Version     | Release Date | Language Platform |
|:---------|:-----------:|:------------:|:------------------|
| Pegasus  | 0.15.0      | 10/02/2019   | Python 3.7.3 |
| SCANPY   | 1.4.4.post1 | 07/29/2019   | Python 3.7.3 |
| Seurat   | 3.1.0       | 08/20/2019   | R 3.6.2 |

For versions of software dependencies, please refer to information in our [Dockerfile](https://raw.githubusercontent.com/lilab-bcb/cumulus-experiment/master/docker/Dockerfile).

Non-standard hardware is not required.

## Overview

This benchmark repository consists of the following main parts:

1. **Software usage:** Instructions on both [software installation](#software-installation) and [docker execution](#run-as-a-docker-container) are provided.
2. [Dataset description](#raw-data): Introduce the 3 datasets used in benchmark, and how to retrieve them.
3. [Software demo](#demo): A demo on running Pegasus for downstream analysis.
4. **Reproducing paper results:** Reproduce results on paper.
    * [Precalculation](#precalculation-for-reproducibility): Pegasu PCA and Diffusion map results are compter-dependent. So they need to be loaded in order to generate the same figures in the paper. 
    * [Highly variable feature selection](#experiment-on-highly-variable-feature-selection)
    * [Batch correction](#experiment-on-batch-correction)
    * [Nearest neighbors](#experiment-on-k-nearest-neighbor)
    * **Other steps:** [Diffusion maps](#experiment-on-diffusion-maps), [clustering](#experiment-on-clustering-algorithms), and [visualization](#experiment-on-visualization-methods).
    * [Analysis tasks](#benchmark-on-analysis-tasks): Benchmark on Pegasus, Seurat, and SCANPY over main analysis tasks.
5. **Scalability benchmark**
    * [10X channel](#number-of-10x-channels)
    * [CPU and Data Size](#cpu-and-data-size)
    * [Memory usage](#memory)
6. [Workflow](#benchmark-on-workflows): Benchmark cumulus workflow.

## Software Installation

**Cumulus** is a cloud-based framework for single-cell/single-nucleus RNA-Seq data analysis. It accounts for processing from sequencing output extraction down to mining biological knowledge from gene-count matrix. It's used as [Terra](https://app.terra.bio/) workflows. Its open-source GitHub repository is [here](https://github.com/klarman-cell-observatory/cumulus), and its documentation can be found [here](https://cumulus.readthedocs.io).

**Pegasus** is the analysis module of Cumulus, which is written in Python. Its GitHub repository is [here](https://github.com/klarman-cell-observatory/pegasus), with documentation [here](https://pegasus.readthedocs.io). Pegasus is avalabile on PyPI with package name [pegasuspy](https://pypi.org/project/pegasuspy/).

### Run with Docker

We recommending using Docker to run our docker image for reproducing our benchmark results, as this is a way with **no installation** but just Docker itself (if you don't have it on your computer):

1. Install Docker on your computer following instructions [here](https://docs.docker.com/v17.09/engine/installation/#supported-platforms), if you don't have Docker yet.

2. Sign up with an account for [Docker Hub](https://hub.docker.com/).

3. Sign in for docker on your computer with your Docker Hub account:

```
docker login
```

4. Pull our docker image public on Docker Hub to your computer:

```
docker pull cumulusprod/cumulus-experiment:19.11
```

Then see Section [Run as a Docker Container](#run-as-a-docker-container) for how to run it as a docker container on your computer.

### Local Installation

Otherwise, if you want to try Pegasus on your machine, please following its installation [here](https://pegasus.readthedocs.io/en/0.15.0/installation.html).

## Run as a Docker Container

In terminal of your computer, type the following command to run our docker image as a docker container:

```
docker run -it --rm --name my-experiment -v /path-to-output:/output cumulusprod/cumulus-experiment:19.11
```

where 
* ``/path-to-output`` is the local directory to which you want to set the experiment output;
* ``my-experiment`` is the container name, which can be changed to your preferred name.

Then type the following commands in terminal to fetch the most recent version of experiment scripts:

```
cd /experiment
git pull
```

Notice that there are 2 conda environments already installed: ``pegasus-env`` and ``scanpy-env``. You can activate/deactivate either of them with the following commands (taking ``pegasus-env`` as the example):

```
root# source activate pegasus-env
(pegasus-env) root# conda deactivate
root#
```

where ``root#`` and ``(pegasus-env) root#`` are environment information automatically appearing in terminal. Similarly below.

To detach from this container, press ``Ctrl`` + ``p``, then ``Ctrl`` + ``q``. 

To attach back, type ``docker attach my-experiment``, where ``my-experiment`` is the container name you set in ``docker run`` with ``--name`` option.

To terminate the container, if inside it, type ``exit``; if outside, type ``docker container stop my-experiment``.

## Raw Data

### Human Bone Marrow Dataset

This dataset has 378,000 cells and 33,694 genes before quality control. It consists of 63 channels collected from 8 donors. Donor 6 has 7 channles in use, while each of the other donors provides 8 channels.

This dataset is available at https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79 in **csv**, **loom** and **mtx** formats.

For the experiments, we've provided its **10X h5** and **h5sc** (Cumulus h5) formats in ``/data`` folder in the docker image. Below is their summary:

| File | Description |
|:-----|:------------|
| ``/data/MantonBM_nonmix_10x.h5`` | Bone Marrow dataset in **10X h5** format. Used for SCANPY and Seurat. |
| ``/data/MantonBM_nonmix.h5sc`` | Bone Marrow dataset in **h5sc** format. Used for Pegasus. |
| ``/data/MantonBM_nonmix_tiny.h5sc`` | A subset of ``MantonBM_nonmix.h5sc`` of 8 samples, each from one donor. Used for batch correction benchmark. |

### PBMC Dataset

This data set has 5,025 cells and 33,538 genes, which come from 1 channel. It can be downloaded in command line after running this docker as container:

```
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_filtered_feature_bc_matrix.h5 -O /data/5k_pbmc_v3.h5
```

When finished, the data file is ``/data/5k_pbmc_v3.h5``. It will be used for several benchmark on analysis tasks.

### Mouse Brain Dataset

This dataset has 1,306,127 cells and 27,998 genes before quality control, with 133 channels. It can be downloaded in command line after running this docker as container:

```
wget http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5 -O /data/1M_neurons.h5
```

When finished, the data file is ``/data/1M_neurons.h5``. It will be used for benchmark on analysis tasks.

## Demo

Before reproducing our benchmark and experiment in paper, you can try Pegasus with data provided in our docker as a demo.

### Execution

In the docker container, run the following commands:

```
root# source activate pegasus-env
(pegasus-env) root# pegasus cluster -p 8 --output-filtration-results --plot-filtration-results --correct-batch-effect --diffmap --spectral-leiden --fitsne /data/MantonBM_nonmix.h5sc /output/demo_out
```

This runs clustering with 8 threads: generate Quality-Control (QC) summary as spreadsheet and plots, apply batch correction, compute Diffusion Maps, cluster using Spectral Leiden algorithm, calculate FIt-SNE embedding. Steps like PCA and kNN are done by default. 

Notice that you can add ``--knn-full-speed`` to run kNN with multiple threads. But because this can reduce the reproducibility on kNN result, we choose to do it with single core here. 

For details on options of ``pegasus cluster`` command, please see [here](https://pegasus.readthedocs.io/en/latest/usage.html#pegasus-cluster).

After that, apply Differential Expression (DE) analysis and cell type annotation on clusters:

```
(pegasus-env) root# pegasus de_analysis -p 8 --labels spectral_leiden_labels --t /output/demo_out.h5ad /output/demo_out.de.xlsx
(pegasus-env) root# pegasus annotate_cluster /output/demo_out.h5ad /output/demo_out.anno.txt
```

The DE analysis also uses 8 threads. It applies only Welch's t-test on clusters, and the putative cell type annotation simply uses this test result.

Finally, generete FIt-SNE plots of data regarding cluster labels and channels side-by-side: 

```
(pegasus-env) root# pegasus plot scatter --basis fitsne --attributes spectral_leiden_labels,Channel /output/demo_out.h5ad /output/demo_out.fitsne.pdf
```

### Expected Output

When finished, in folder ``/output``, you can find the following results using the following commands in command-line:

| File | Description |
|:-----|:------------|
| ``demo_out.h5ad`` | Analysis result in **h5ad** format. |
| ``demo_out.filt.xlsx`` | Quality-Control (QC) summary as an Excel sheet. |
| ``demo_out.filt.{UMI, gene, mito}.pdf`` | QC plots regarding UMIs, barcodes, and mitochondrial genes. |
| ``demo_out.de.xlsx`` | Differential Expression (DE) analysis result in Excel format. |
| ``demo_out.anno.txt`` | Cell type annotation on each cluster. |
| ``demo_out.fitsne.pdf`` | FIt-SNE plot on dataset regarding clusters and channels side-by-side. |

For ``demo_out.h5ad``, to load it, in docker container, run

```
(pegasus-env) root# python
>>> import pegasus as pg
>>> adata = pg.read_input("/output/demo_out.h5ad")
```

to load it as an [anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.html) object.

For all the other output files, you can open them in ``/path-to-out`` outside the docker container.

### Expected Runtime

We tested this demo on a MacBook laptop with 2.9 GHz 6-Core Intel i9 CPU (i.e. 12 CPUs to use), 32GB memory, MacOS 10.15.1, and Docker Desktop 2.1.0.4 (with Docker engine 19.03.4). The overall runtime was 18 minutes.

## Reproducing Paper Results

### Precalculation for Reproducibility

As PCA and Diffusion Map results are different among different machines, we use precalculated ones done by our experiment server for reproducint paper results. You can find them in ``/data/precalculated`` folder inside docker image.

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
(pegasus-env) root# python generate_data_info.py 5k_pbmc
```

to get all the necessary data for experiment on PBMC dataset.

Execute

```
(pegasus-env) root# python generate_data_info.py 1M_neurons
```

to get all the necessary data for experiment on Mouse Neuron dataset. 

Notice that Mouse Neuron dataset was only used for runtime benchmark. And since it's a huge dataset, and its processing is memory-consuming, if your computer doesn't have a memory large enough (e.g. 16GB is not enough to hold it due to our test), please consider to only try Bone Marrow and PBMC datasets. All the figures are based on Bone Marrow dataset.

### Experiment on Highly Variable Feature Selection

First, enter the experiment subfolder:

```
(pegasus-env) root# cd /experiment/highly_variable_features
```

Then execute the experiment:

```
(pegasus-env) root# python gen_result.py
```

When finished, you'll find the figures generated in ``/output``, and CSV files containing markers in the current folder. Besides, lists of markers convered by different highly variable feature selection methods are in the following files under folder ``/experiment/highly_variable_features``:

| File | Description |
|------|-------------|
| ``immune_genes.txt`` | Total list of marker genes for comparison. |
| ``pegasus_markers.txt`` | Marker genes covered by Pegasus HVF selection method. |
| ``pegasus_specific.txt`` | Marker genes covered ONLY by Pegasus HVF selection method. |
| ``seurat_markers.txt`` | Marker genes covered by Seurat HVF selection method. |
| ``seurat_specific.txt`` | Marker genes covered ONLY by Seurat HVF selection method. |
| ``common_markers.txt`` | Marker genes covered by both methods. |

### Experiment on Batch Correction

#### Package Summary

The following batch correction methods are compared. The benchmark dataset is the 8-channel subset of bone marrow data: ``/data/MantonBM_nonmix_tiny.h5sc``.

| Method | Package | Version | Release Date |
|--------|---------|---------|--------------|
| L/S adjustment | Pegasus | 0.15.0 | 10/02/2019 |
| ComBat | SCANPY | 1.4.4.post1 | 07/29/2019 |
| MNN | mnnpy | 0.1.9.5 | 02/24/2019 |
| BBKNN | bbknn | 1.3.6 | 08/22/2019 |
| CCA | Seurat | 3.1.0 | 08/20/2019 |

#### Ground Truth Cell Types

As MNN and Seurat CCA both fail for the whole Bone Marrow dataset, we use a subset of it, by selecting one channel per donor, for the batch correction benchmark.

First, we need a list of cell types as the ground truth to measure the performance of different batch correction methods. We use the information from the clustering on the dataset not batch-corrected for this purpose.

Enter the folder for generating ground truth information, and execute the script:

```
(pegasus-env) root# cd /experiment/batch_correction/ground
(pegasus-env) root# python gen_ground_h5ad.py
(pegasus-env) root# python gen_celltype.py
```

This will generate a file containing ground truth of cell types (``ground_cell_types.txt``) in folder ``/experiment/batch_correction``.

Notice that we use PCA coordinates precalculated on our server, so that the clustering result and ground truth are consistent with those shown in the paper.

#### Benchmark on Baseline Method

The baseline is running Pegasus clustering without batch correction. It's already done above, and result is stored as ``/experiment/batch_correction/ground/ground.h5ad``.

#### Benchmark on Pegasus Batch Correction Method

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/batch_correction/pegasus
(pegasus-env) root# python run_pegasus_batch_correct.py
```

When finished, you'll have a result file ``pegasus_corrected.h5ad``, along with its log file ``pegasus_correct.log`` in the current folder. To calculate its batch correction time, in ``pegasus_correct.log``, add up the time spent on **Estimation on feature statistics per channel** and **Batch correction** steps.

#### Benchmark on ComBat

In SCANPY environment, run the following commands:

```
(scanpy-env) root# cd /experiment/batch_correction/combat
(scanpy-env) root# python scanpy_combat.py
```

When finished, you'll have a result file ``scanpy_combat_corrected.h5ad`` in the corrent folder, and you can read the time spent on ComBat from the screen output.

#### Benchmark on MNN

In SCANPY environment, run the following commands:

```
(scanpy-env) root# cd /experiment/batch_correction/mnn
(scanpy-env) root# python scanpy_mnn.py
```

When finished, you'll have a result file ``scanpy_mnn_corrected.h5ad`` in the current folder, and you can read the time spent on MNN from the screen output.

#### Benchmark on BBKNN

In SCANPY environment, run the following commands:

```
(scanpy-env) root# cd /experiment/batch_correction/bbknn
(scanpy-env) root# python scanpy_bbknn.py
```

When finished, you'll have a result file ``scanpy_bbknn_corrected.h5ad`` in the current folder, and you can read the time spent on BBKNN from the screen output.

#### Benchmark on Seurat CCA

Running Seurat batch correction doesn't depend on Python environment. Run the following commands:

```
root# cd /experiment/batch_correction/seurat
root# Rscript seurat_cca.R
```

When finished, you'll have 3 result files: ``matrix.mtx`` for count matrix, ``barcodes.txt`` for cell barcode names, ``genes.txt`` for gene names in the current folder. Besides, there is also a log file ``seurat_cca.log`` in the folder. To get its batch correction time, add up the time spent on **Finding Anchors** and **Integration** together.

#### Measure Performance of All Batch Correction Methods

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

### Experiment on k-Nearest-Neighbor

#### Package Summary

kNN methods of Pegasus, SCANPY, and Seurat are compared. The ground truth of accurate kNN is achieved by bruth force method in ``scikit-learn``. We list the kNN packages that these softwares uses as follow:

| Package | Used By | Version | Release Date |
|---------|---------|---------|--------------|
| scikit-learn | Ground Truth | 0.21.3 | 07/29/2019 |
| hnswlib | Pegasus | 0.3.2.0 | 08/23/2019 |
| umap-learn | SCANPY | 0.3.10 | 08/14/2019 |
| RcppAnnoy | Seurat | 0.0.13 | 09/23/2019 | 

#### Generate Ground Truth kNN

We use the kNN result by brute force algorithm in _scikit-learn_ as the ground truth. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py brute
```

When finished, you'll have a result file ``baseline_indices.npy`` in the current folder.

#### Benchmark on Pegasus kNN Method

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py pegasus
```

When finished, you'll have a result file ``pegasus_indices.npy`` in the current folder, and you can read the time spent on kNN from screen output.

#### Benchmark on SCANPY kNN Method

SCANPY uses kNN from _umap-learn_ package. So we directly benchmark this function. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py scanpy
```

When finished, you'll have a result file ``scanpy_indices.npy`` in the current folder, and you can read the time spent on kNN from screen output.

#### Benchmark on Seurat kNN Method

Seurat's kNN has two methods: ``nn2`` from _RANN_ package, default method but time-consuming; ``AnnoyNN`` from _RcppAnnoy_ package, not default but more efficient. We decide to choose ``AnnoyNN`` for this benchmark. In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# Rscript seurat_knn.R
(pegasus-env) root# python knn_comp.py seurat
```

When finished, you'll have a result file ``seurat_indices_annoy.txt`` in the current folder, and you can read the time spent on kNN from screen output.

#### Measure kNN Methods

Given that all the kNN resuls of the 3 methods are calculated, we are ready to measure the recall of them, and generate their time and recall plots. 

First, update ``time_stats.txt`` in ``/experiment/knn_comparison`` by the time information you saw in benchmarks above. Then in Pegasus environment, type:

```
(pegasus-env) root# cd /experiment/knn_comparison
(pegasus-env) root# python knn_comp.py plot
```

This will generate the corresponding figures in ``/output``.

### Experiment on Diffusion Maps

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/diffmap
(pegasus-env) root# python get_diffmap_figures.py
```

When finished, you'll find the figures generated in ``/output``.

### Experiment on Clustering Algorithms

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/clustering
(pegasus-env) root# python algorithm_compare.py
(pegasus-env) root# python spectral_clustering.py
```

When finished, you'll find the figures generated in ``/output``, and AMI results can be read from the screen output. Besides, you can check ``/experiment/pegasus.log`` for the execution time on each of the 4 clustering algorithms, while the time on spectral clustering can be read from the screen output.

### Experiment on Visualization Methods

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/visualization
(pegasus-env) root# python origin_vs_net.py
```

When finished, you'll find the figures generated in ``/output``, and for each visualization method, its kSIM regarding Louvain clustering labels can be seen from screen output. Besides, you can check ``/experiment/pegasus.log`` for the execution time on each visualization method.

### Benchmark on Analysis Tasks

All 3 packages are benchmarked on PBMC and Bone Marrow datasets, while only Pegasus and SCANPY are benchmarked on Mouse Neuron, because Seurat fails at loading the count matrix step for this big dataset.

#### PBMC Dataset

Only 8 CPUs are used for benchmark this small dataset.

##### Pegasus

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/overall/pbmc
(pegasus-env) root# python run_pegasus_pbmc.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/pbmc/pbmc_pegasus_cpu_8.log``.

##### SCANPY

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/overall/pbmc
(scanpy-env) root# python run_scanpy_pbmc.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/pbmc/pbmc_scanpy_cpu_8.log``.

##### Seurat

First, in Pegasus environment, run the following commands to generate a Seurat-compatible h5ad file on the dataset, and convert it into Seurat object format:

```
(pegasus-env) root# cd /experiment/overall/pbmc
(pegasus-env) root# ./get_seurat_compatible_pbmc.sh
(pegasus-env) root# Rscript convert_pbmc_pegasus.R
```

When finished, you'll have a Seurat object in ``/experiment/overall/pbmc/5k_pbmc_v3.RData`` for benchmarking steps starting from kNN. 

Now in SCANPY environment, run the following command:

```
(scanpy-env) root# Rscript run_seurat_pbmc.R
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/pbmc/pbmc_seurat_cpu_8.log``.

#### Manton Bone Marrow Dataset

All 28 CPUs are used for benchmark on this dataset.

##### Pegasus

To benchmark Pegasus, in Pegasus environment, run the following command:

```
(pegasus-env) root# cd /experiment/overall/MantonBM
(pegasus-env) root# python run_pegasus_mantonbm.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/MantonBM/mantonbm_pegasus_cpu_28.log``.

##### SCANPY

To benchmark SCANPY, in SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/overall/MantonBM
(scanpy-env) root# python run_scanpy_mantonbm.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/MantonBM/mantonbm_scanpy_cpu_28.log``.

##### Seurat

The benchmark on Seurat is a little bit complicated, as it's written in R. 

First, in Pegasus environment, run the following commands to generate a Seurat-compatible h5ad file on the dataset, and convert it into Seurat object format:

```
(pegasus-env) root# cd /experiment/overall/MantonBM
(pegasus-env) root# ./get_seurat_compatible_mantonbm.sh
(pegasus-env) root# Rscript convert_mantonbm_pegasus.R
```

When finished, you'll have a Seurat object in ``/experiment/overall/MantonBM/MantonBM_nonmix.RData`` for benchmarking steps starting from kNN. 

Now in SCANPY environment, run the following command:

```
(scanpy-env) root# Rscript run_seurat_mantonbm.R
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/MantonBM/mantombm_seurat_cpu_28.log``.

As Seurat fails in the Batch correction and Leiden clustering steps, I make them as two separate R scripts for users to try themselves.

Run the following command to benchmark on Batch correction using Seurat:

```
root# Rscript seurat_batch_correction.R
```

When terminated with failure, you'll find the time information in its log file ``/experiment/overall/MantonBM/seurat_batch_correction.log``, and error message from screen output.

In SCANPY environment, run the following command to benchmark on Leiden clustereing using Seurat:

```
(scanpy-env) root# Rscript seurat_leiden.R
```

When terminated with failure, you may find information in its log file ``/experiment/overall/MantonBM/seurat_leiden.log`` and screen output.

#### Mouse Neuron Dataset

All 28 CPUs are used for benchmark on this dataset.

##### Pegasus

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/overall/mouse_neuron
(pegasus-env) root# python run_pegasus_1m.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/mouse_neuron/1m_pegasus_cpu_28.log``.

##### SCANPY

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/overall/mouse_neuron
(scanpy-env) root# python run_scanpy_1m.py
```

When finished, you'll find execution time for each step in its log file ``/experiment/overall/mouse_neuron/1m_scanpy_cpu_28.log``.

Specifically, as the last step, FLE embedding calculation, takes a significantly long time, we put it into a seperate script, and should be run after done with steps above:

```
(scanpy-env) root# python run_scanpy_1m_fle.py
```

When finished, execution time for FLE step will be appended into the same log file above.

##### Seurat

Seurat fails at loading data step. Users can try the following commands in R environment, and check out the error message:

```
> library(Seurat)
> adata <- Read10X("/data/1M_neurons.h5")
```

## Scalability

### Number of 10X channels

The scalability benchmark on cumulus *mkfastq* and *count* steps over Bone Marrow dataset's 63 channels was run on Terra + Google Cloud platform. Runtime per channel is recorded in ``/experiment/scalability/channel/channel_stats.csv``, and total cost is $101.43.

Then in Pegasus environment, running

```
(pegasus-env) root# python plot_channel_stats.py
```

will generate amortized cost per channel in a new CSV file: ``/experiment/scalability/channel/channel_stats_updated.csv``, and 2 plots.

### CPU and Data Size

The scalability on runtime of Pegasus is performed with respect to:

* **Number of threads:** We use 8 CPUs to simulate a normal laptop, and 28 CPUs to simulate a server.
* **Dataset size:** We down-sample Bone Marrow datasets into subsets of sizes 5k, 10k, 25k, 50k, 100k, and 200k cells, along with the whole dataset of 274k cells.

In Pegasus environment, run the following commands

```
(pegasus-env) root# cd /experiment/scalability/runtime
(pegasus-env) root# python down_sample.py sampling
```

to generate subsamples of Bone Marrow dataset of different number of cells.

Then run

```
(pegasus-env) root# python down_sample.py benchmark
```

to start the benchmark on Pegasus over all these subsamples with different number of threads.

When finished, all the ``.log`` files in this folder record the runtime of each step of analysis using Pegasus. The subsample size is inferred from filename.

### Memory

The scalability on memory usage is performed on Pegasus, Seurat, and SCANPY. Only PBMC and Bone Marrow datasets are used, as Seurat fails on loading Mouse Neuron dataset. 

Moreover, since Seurat doesn't have Diffusion maps and FLE features, and it fails at batch correction and Leiden clustering steps on Bone Marrow dataset, we only benchmark 6 analysis tasks here: Highly variable feature selection, PCA, kNN, Louvain clustering, tSNE, and UMAP.

#### PBMC Dataset

Only 8 CPUs are used for this small dataset.

##### Pegasus

In Pegasus environment, run the following commands:

```
(pegasus-env) root# cd /experiment/scalability/memory/pbmc
(pegasus-env) root# mem_monitor.sh > monitor_pegasus_pbmc.log &
(pegasus-env) root# python run_pegasus_pbmc.py
```

When finished, in the same folder, you'll find execution time for each step in its log file ``pbmc_pegasus.log``, and the memory usage during the execution in ``monitor_pegasus_pbmc.log``.

When you are done, you also need to stop logging the memory usage. simply run the following command in your terminal:

```
kill <number>
```

where ``<number>`` is the process ID shown after running ``mem_monitor.sh > monitor_pegasus_pbmc.log &`` command. Similarly as below.

##### SCANPY

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/scalability/memory/pbmc
(scanpy-env) root# mem_monitor.sh > monitor_scanpy_pbmc.log &
(scanpy-env) root# python run_scanpy_pbmc.py
```

When finished, in the same folder, you'll find execution time for each step in its log file ``pbmc_scanpy.log``, and the memory usage during the execution in ``monitor_scanpy_pbmc.log``.

Don't forget to kill the process on logging the memory usage when you are done.

##### Seurat

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/scalability/memory/pbmc
(scanpy-env) root# mem_monitor.sh > monitor_seurat_pbmc.log &
(scanpy-env) root# Rscript run_pegasus_pbmc.py
```

When finished, in the same folder, you'll find execution time for each step in its log file ``pbmc_seurat.log``, and the memory usage during the execution in ``monitor_seurat_pbmc.log``.

Don't forget to kill the process on logging the memory usage when you are done.

#### Bone Marrow Dataset

All 28 CPUs are used for memory benchmark on this dataset.

##### Pegasus

In Pegasus environment, run the following command:

```
(pegasus-env) root# cd /experiment/scalability/memory/MantonBM
(pegasus-env) root# mem_monitor.sh > monitor_pegasus_mantonbm.log &
(pegasus-env) root# python run_pegasus_mantonbm.py
```

When finished, in the same folder, you'll find execution time for each step in its log file ``pbmc_pegasus.log``, and the memory usage during the execution in ``monitor_pegasus_mantonbm.log``.

When you are done, you also need to stop logging the memory usage. simply run the following command in your terminal:

```
kill <number>
```

where ``<number>`` is the process ID shown after running ``mem_monitor.sh > monitor_pegasus_mantonbm.log &`` command. Similarly as below.

##### SCANPY

In SCANPY environment, run the follow commands:

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/scalability/memory/MantonBM
(scanpy-env) root# mem_monitor.sh > monitor_scanpy_mantonbm.log &
(scanpy-env) root# python run_scanpy_mantonbm.py
```

When finished, in the same folder, you'll find execution time for each step in its log file ``mantonbm_scanpy.log``, and the memory usage during the execution in ``monitor_scanpy_mantonbm.log``.

Don't forget to kill the process on logging the memory usage when you are done.

##### Seurat

In SCANPY environment, run the following command:

```
(scanpy-env) root# cd /experiment/scalability/memory/MantonBM
(scanpy-env) root# mem_monitor.sh > monitor_seurat_mantonbm.log &
(scanpy-env) root# Rscript run_pegasus_mantonbm.py
```

When finished, in the same folder, you'll find execution time for each step in its log file ``mantonbm_seurat.log``, and the memory usage during the execution in ``monitor_seurat_mantonbm.log``.

Don't forget to kill the process on logging the memory usage when you are done.

### Benchmark on Workflows

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

#### Cumulus

Benchmark on Cumulus is done by running jobs on Terra via Cumulus WDL workflows. And its overall execution time includes all Terra or Google Cloud specific preprocessing and postprocessing phases.

To run it on Terra, please following [cumulus documentation](https://cumulus.readthedocs.io/en/latest/cumulus.html) and [our tutorial video](https://www.youtube.com/watch?v=zKgo3mf4uRk&t=2s) for this benchmark.

Notice that you need to upload ``/data/MantonBM_nonmix.h5sc`` to the Google bucket of your workspace via [gsutil](https://cloud.google.com/storage/docs/gsutil), change parameters ``"cumulus.input_file"`` and ``"cumulus.output_name"`` in ``/experiment/cloud/inputs_32.cpu.json`` to your own, and upload this JSON file in **cumulus** workflow page in your workspace.

#### SCANPY

As SCANPY doesn't have a cloud-based interface, its benchmark is performed on a Google Cloud VM. Besides, since we used BBKNN for batch correction, SCANPY doesn't need to find kNN.

Similarly as benchmarks before, in SCANPY environment, run the following commands

```
(scanpy-env) root# cd /experiment/cloud
(scanpy-env) root# python run_scanpy.py > scanpy.log
```

When finished, you'll find execution time for each step in its log file ``/experiment/cloud/scanpy.log``.

#### Seurat

Similarly as SCANPY, benchmark on Seurat is performed on a Google Cloud VM, and used as a single-server solution.

Besides, as Seurat would fail for batch correction when using 63 channels as the batches, we instead use batch correction with 8 donors being the batches. Moreover, the batch correction failed when using 10, 15, 20, and 32 threads via R _future_ package. So we simply used 2 threads just to make sure the batch correction step terminates successfully. Then for all the other steps, 32 threads are used whenever possible.

In SCANPY environment, run the following commands

```
(scanpy-env) root# cd /experiment/cloud
(scanpy-env) root# Rscript run_seurat_hvg_batch_correction.R
(scanpy-env) root# Rscript run_seurat_analysis.R
```

When finished, you'll find execution time for each step in log files ``/experiment/cloud/seurat_batch_correction.log`` and ``/experiment/cloud/seurat_analysis.log``.