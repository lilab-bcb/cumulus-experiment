FROM ubuntu:18.04

# Install ssl, curl, ssh, and git
RUN apt-get update -qq && \
	apt-get install -y libssl-dev libcurl4-openssl-dev libssh2-1-dev libhdf5-dev gfortran emacs && \
	apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		git \
		zlib1g-dev \
		libxml2-dev \
		cmake \
		gnupg \
		wget \
		lsb-release \
		libfftw3-dev && \
	apt-get clean

## Install Oracle JDK 12.0.1 ##
ADD jdk-12.0.1_linux-x64_bin.tar.gz /opt/software

## Add script to track resource usage ##
RUN mkdir -p /opt/software/scripts
ADD mem_monitor.sh /opt/software/scripts
RUN chmod a+rx /opt/software/scripts/mem_monitor.sh

ENV PATH=/opt/software/jdk-12.0.1/bin:/opt/software/scripts:$PATH

# Install R 3.6
RUN wget http://security.ubuntu.com/ubuntu/pool/main/t/tzdata/tzdata_2019c-0ubuntu0.18.04_all.deb && \
	dpkg -i tzdata_2019c-0ubuntu0.18.04_all.deb && \
	rm tzdata_2019c-0ubuntu0.18.04_all.deb

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list && \
	apt-get update && \
	apt-get -y install r-base

ADD install_dependencies.R /opt/software/install_dependencies.R
RUN Rscript /opt/software/install_dependencies.R

## Get Seurat 3.1.0 source code ##
RUN wget https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.1.0.tar.gz && \
	tar -zxvf Seurat_3.1.0.tar.gz && \ 
	mv Seurat /opt/software/seurat-3.1.0 && \
	rm Seurat_3.1.0.tar.gz

## INSTALL SEURAT ##
ADD install_seurat.R /opt/software/install_seurat.R
RUN Rscript /opt/software/install_seurat.R

## INSTALL FAST_TSNE for Seurat ##
RUN git clone https://github.com/KlugerLab/FIt-SNE.git /opt/software/FIt-SNE && \
	cd /opt/software/FIt-SNE && \
	g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

ENV PATH=$PATH:/opt/software/FIt-SNE/bin

# Install Python environment
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	bash ./Miniconda3-latest-Linux-x86_64.sh -b

ENV PATH=/root/miniconda3/bin:$PATH
ARG old_path=$PATH

# Install pegasus environment
RUN conda create -n pegasus-env python=3.7.3 -y

RUN conda install -n pegasus-env -y -c anaconda mkl==2019.4 && \
	conda install -n pegasus-env -y -c anaconda mkl_fft==1.0.14 && \
	conda install -n pegasus-env -y -c anaconda mkl_random==1.0.2 && \
	conda install -n pegasus-env -y -c anaconda blas==1.0 && \
	conda install -n pegasus-env -y -c anaconda intel-openmp==2019.4 && \
	conda install -n pegasus-env -y -c anaconda numpy==1.16.4 && \
	conda install -n pegasus-env -y -c conda-forge fftw==3.3.8

ENV PATH=/root/miniconda3/envs/pegasus-env/bin:$old_path
RUN	pip install Cython==0.29.13 && \
	pip install docopt==0.6.2 && \
	pip install fitsne==1.0.1 && \
	pip install h5py==2.9.0 && \
	pip install joblib==0.13.2 && \
	pip install kiwisolver==1.1.0 && \
	pip install leidenalg==0.7.0 && \
	pip install lightgbm==2.2.1 && \
	pip install llvmlite==0.29.0 && \
	pip install matplotlib==3.1.1 && \
	pip install natsort==6.0.0 && \
	pip install numba==0.45.1 && \
	pip install numexpr==2.7.0 && \
	pip install pandas==0.25.1 && \
	pip install pybind11==2.3.0 && \
	pip install pycparser==2.19 && \
	pip install pyparsing==2.4.2 && \
	pip install python-igraph==0.7.1.post6 && \
	pip install scikit-learn==0.21.3 && \
	pip install scikit-misc==0.1.1 && \
	pip install scipy==1.2.2 && \
	pip install six==1.12.0 && \
	pip install statsmodels==0.10.1 && \
	pip install tables==3.5.2 && \
	pip install umap-learn==0.3.10 && \
	pip install anndata==0.6.22.post1 && \
	pip install termcolor==1.1.0 && \
	pip install hnswlib==0.3.2.0 && \
	pip install pegasuspy==0.15.0

# Install scanpy environment
RUN conda create -n scanpy-env python=3.7.3 -y

RUN conda install -n scanpy-env -y -c anaconda mkl==2019.4 && \
	conda install -n scanpy-env -y -c anaconda mkl_fft==1.0.14 && \
	conda install -n scanpy-env -y -c anaconda mkl_random==1.0.2 && \
	conda install -n scanpy-env -y -c anaconda blas==1.0 && \
	conda install -n scanpy-env -y -c anaconda intel-openmp==2019.4 && \
	conda install -n scanpy-env -y -c anaconda numpy==1.16.4

ENV PATH=/root/miniconda3/envs/scanpy-env/bin:$old_path
RUN pip install Cython==0.29.13 && \
	pip install h5py==2.9.0 && \
	pip install joblib==0.13.2 && \
	pip install kiwisolver==1.1.0 && \
	pip install lightgbm==2.2.1 && \
	pip install llvmlite==0.29.0 && \
	pip install matplotlib==3.0.3 && \
	pip install natsort==6.0.0 && \
	pip install networkx==2.3 && \
	pip install numba==0.45.1 && \
	pip install numexpr==2.7.0 && \
	pip install pandas==0.25.1 && \
	pip install patsy==0.5.1 && \
	pip install pybind11==2.3.0 && \
	pip install python-igraph==0.7.1.post6 && \
	pip install scikit-learn==0.21.3 && \
	pip install scikit-misc==0.1.1 && \
	pip install scipy==1.3.1 && \
	pip install six==1.12.0 && \
	pip install statsmodels==0.10.1 && \
	pip install tables==3.5.2 && \
	pip install tqdm==4.32.2 && \
	pip install anndata==0.6.22.post1 && \
	git clone https://github.com/lilab-bcb/bbknn.git /opt/bbknn && cd /opt/bbknn && pip install -e . && \
	pip install scanpy==1.4.4.post1 && \
	pip install MulticoreTSNE==0.1 && \
	pip install louvain==0.6.1 && \
	pip install leidenalg==0.7.0 && \
	pip install umap-learn==0.3.10 && \
	pip install mnnpy==0.1.9.5 && \
	pip install fa2==0.3.5 && \
	pip install termcolor==1.1.0

ENV PATH=$old_path

# Add experiment scripts
RUN git clone https://github.com/lilab-bcb/cumulus-experiment.git /experiment

# Add Data
ADD data /data

WORKDIR /experiment