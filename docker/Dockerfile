FROM ubuntu:18.04

# Install ssl, curl, ssh, and git
RUN apt-get update -qq && \
	apt-get install -y libssl-dev libcurl4-openssl-dev libssh2-1-dev libhdf5-dev gfortran vim emacs && \
	apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		git \
		zlib1g-dev \
		libxml2-dev \
		cmake \
		gnupg \
		wget \
		vim \
		emacs \
		lsb-release \
		libfftw3-dev \
		libboost-iostreams-dev && \
	apt-get clean

# Install Oracle JDK 12
ADD jdk-12.0.1_linux-x64_bin.tar.gz /opt
ENV PATH=/opt/jdk-12.0.1/bin:$PATH

# Install R 3.6
RUN wget http://security.ubuntu.com/ubuntu/pool/main/t/tzdata/tzdata_2019b-0ubuntu0.18.04_all.deb && \
	dpkg -i tzdata_2019b-0ubuntu0.18.04_all.deb && \
	rm tzdata_2019b-0ubuntu0.18.04_all.deb

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list && \
	apt-get update && \
	apt-get -y install r-base

## INSTALL SEURAT ##
ADD install.R /opt/install.R
RUN wget https://github.com/satijalab/seurat/archive/v3.0.2.tar.gz && \
	tar -zxvf v3.0.2.tar.gz && \
	mv seurat-3.0.2/ /opt/seurat-3.0.2 && \
	rm v3.0.2.tar.gz && \
	Rscript /opt/install.R

## INSTALL FAST_TSNE for Seurat ##
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz && \
	tar -zxvf fftw-3.3.8.tar.gz && \
	mv fftw-3.3.8/ /opt/fftw-3.3.8 && \
	rm fftw-3.3.8.tar.gz && \
	cd /opt/fftw-3.3.8 && \
	./configure && \
	make && \
	make install && \
	git clone https://github.com/KlugerLab/FIt-SNE.git /opt/FIt-SNE && \
	cd /opt/FIt-SNE && \
	g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

ENV PATH=$PATH:/opt/FIt-SNE/bin

# Install Python environment
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	bash ./Miniconda3-latest-Linux-x86_64.sh -b

ENV PATH=/root/miniconda3/bin:$PATH
ARG old_path=$PATH

# Install scCloud environment
RUN conda create -n scCloud -y pip

#RUN	conda install -n scCloud -y -c anaconda numpy && \
#	conda install -n scCloud -y -c anaconda cython && \
#	conda install -n scCloud -y -c anaconda pybind11 && \
#	conda install -n scCloud -y -c conda-forge fftw && \
#	conda install -n scCloud -y -c anaconda pytables

RUN conda install -n scCloud -y -c anaconda numpy && \
	conda install -n scCloud -y -c conda-forge fftw

ENV PATH=/root/miniconda3/envs/scCloud/bin:$old_path
RUN	pip install numpy==1.16.3 && \
	pip install Cython==0.29.7 && \
	pip install docopt==0.6.2 && \
	pip install fitsne==1.0.1 && \
	pip install h5py==2.9.0 && \
	pip install joblib==0.13.2 && \
	pip install kiwisolver==1.1.0 && \
	pip install leidenalg==0.7.0 && \
	pip install lightgbm==2.2.1 && \
	pip install llvmlite==0.28.0 && \
	pip install natsort==6.0.0 && \
	pip install numba==0.43.1 && \
	pip install numexpr==2.6.9 && \
	pip install pandas==0.24.2 && \
	pip install pybind11==2.2.4 && \
	pip install pycparser==2.19 && \
	pip install pyparsing==2.4.0 && \
	pip install python-igraph==0.7.1.post6 && \
	pip install scikit-learn==0.21.1 && \
	pip install scikit-misc==0.1.1 && \
	pip install scipy==1.2.1 && \
	pip install six==1.12.0 && \
	pip install statsmodels==0.10.1 && \
	pip install tables==3.5.1 && \
	pip install umap-learn==0.3.9 && \
	git clone https://github.com/klarman-cell-observatory/scCloudPy.git /opt/scCloudPy && cd /opt/scCloudPy && pip install .
	

# Install scanpy environment
RUN conda create -n scanpy -y pip
ENV PATH=/root/miniconda3/envs/scanpy/bin:$old_path
RUN pip install numpy==1.16.3 && \
	pip install Cython==0.29.7 && \
	pip install h5py==2.9.0 && \
	pip install joblib==0.13.2 && \
	pip install kiwisolver==1.1.0 && \
	pip install lightgbm==2.2.1 && \
	pip install llvmlite==0.28.0 && \
	pip install matplotlib==3.0.3 && \
	pip install natsort==6.0.0 && \
	pip install networkx==2.3 && \
	pip install numba==0.43.1 && \
	pip install numexpr==2.6.9 && \
	pip install pandas==0.24.2 && \
	pip install patsy==0.5.1 && \
	pip install pybind11==2.2.4 && \
	pip install python-igraph==0.7.1.post6 && \
	pip install scikit-learn==0.21.2 && \
	pip install scikit-misc==0.1.1 && \
	pip install scipy==1.3.0 && \
	pip install six==1.12.0 && \
	pip install statsmodels==0.10.1 && \
	pip install tables==3.5.1 && \
	pip install tqdm==4.32.2 && \
	pip install anndata==0.6.22.post1 && \
	git clone https://github.com/lilab-cbb/bbknn.git /opt/bbknn-1.3.4 && cd /opt/bbknn-1.3.4 && pip install -e . && \
	pip install scanpy==1.4.4 && \
	pip install MulticoreTSNE==0.1 && \
	pip install louvain==0.6.1 && \
	pip install leidenalg==0.7.0 && \
	pip install umap-learn==0.3.9 && \
	pip install mnnpy==0.1.9.5 && \
	pip install fa2==0.3.5

ENV PATH=$old_path