# Base Ubuntu
FROM ubuntu:20.04

# Se crean las variantes de trabajo del ambiente
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV R_LIBS_SITE=/usr/local/lib/R/site-library

# Crea la libreria de R
RUN mkdir -p /usr/local/lib/R/site-library && \
    chmod -R 777 /usr/local/lib/R/site-library && \
    chown -R root:staff /usr/local/lib/R && \
    chmod -R g+w /usr/local/lib/R && \
    chmod -R 777 /usr/local/lib/R/site-library

# Actualiza e instala dependencias
RUN apt-get update && apt-get install -y \
    software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update && apt-get install -y \
    libstdc++6 \
    dirmngr \
    apt-transport-https \
    ca-certificates \
    wget \
    build-essential \
    gfortran \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libgsl-dev \
    libssl-dev \
    libz-dev \
    git \
    unzip \
    python3 \
    python3-pip \
    python3-dev \
    default-jdk \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libboost-all-dev \
    cmake \
    samtools \
    bedtools \
    fastqc \
    parallel \
    pigz \
    curl \
    rsync \
    zip \
    bzip2 \
    ncbi-blast+ \
    htop \
    tree \
    screen \
    gawk \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libsqlite3-dev \
    gnupg \
    pandoc \
    pandoc-citeproc \
    gcc \
    g++ \
    make \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Instalación Nextflow
RUN curl -s https://get.nextflow.io | bash
RUN mv nextflow /usr/local/bin/

# Instalacion de R y su configuracion 
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    gpg --dearmor > /usr/share/keyrings/cran-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" > \
    /etc/apt/sources.list.d/cran.list && \
    apt-get update && apt-get install -y r-base r-base-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    chmod -R 777 /usr/local/lib/R/site-library

# Permisos de R
RUN chmod -R 777 /usr/local/lib/R/site-library && \
    chown -R root:staff /usr/local/lib/R && \
    chmod -R g+w /usr/local/lib/R

# Instalacion de paquetes de R
RUN R -e 'install.packages(c("BiocManager", "rmarkdown", "flexdashboard", "plotly", "tidyverse", \
    "DT", "pheatmap", "RColorBrewer", "ggrepel", "scales", "kableExtra", "Matrix", "Rcpp", \
    "RcppArmadillo", "RcppEigen"), repos="https://cloud.r-project.org", dependencies=TRUE)'

# Instalación DESeq2
RUN R -e 'BiocManager::install("DESeq2", update=FALSE, ask=FALSE, force=TRUE); \
    library(DESeq2)'

# Instalación Bioconductor y sus paquetes
RUN R -e 'BiocManager::install(version = "3.20"); \
    BiocManager::install(c("GenomicFeatures", "GenomicAlignments", "Biobase", "S4Vectors", \
    "IRanges", "tximport", "EnhancedVolcano", "edgeR", "limma", "fgsea", \
    "GSEABase", "GSVA", "rtracklayer", "AnnotationDbi", "org.Hs.eg.db"), \
    update=FALSE, ask=FALSE)'

# Instalación de paquetes CRAN
RUN R -e 'install.packages(c("tidyverse", "ggplot2", "viridis", "RColorBrewer", \
    "ggridges", "igraph", "optparse", "flexdashboard", "plotly", "DT", \
    "pheatmap", "ggrepel", "scales", "kableExtra", "data.table", "reshape2", \
    "gridExtra", "ComplexHeatmap", "UpSetR", "corrplot", "VennDiagram", \
    "rmarkdown", "gprofiler2"), repos="https://cloud.r-project.org", dependencies=TRUE)'

# Instalación STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && \
    make STAR && \
    mv STAR /usr/local/bin/ && \
    cd / && \
    rm -rf STAR-2.7.10b 2.7.10b.tar.gz

# Instalación Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz && \
    tar xf salmon-1.9.0_linux_x86_64.tar.gz && \
    mv salmon-1.9.0_linux_x86_64/bin/* /usr/local/bin/ && \
    mv salmon-1.9.0_linux_x86_64/lib/* /usr/local/lib/ && \
    rm -rf salmon-1.9.0_linux_x86_64*

# Install Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/ && \
    rm Trimmomatic-0.39.zip

# Instalación StringTie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar xf stringtie-2.2.1.Linux_x86_64.tar.gz && \
    mv stringtie-2.2.1.Linux_x86_64/stringtie /usr/local/bin/ && \
    rm -rf stringtie-2.2.1.Linux_x86_64*

# Instalación gffcompare
RUN wget https://github.com/gpertea/gffcompare/releases/download/v0.12.6/gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    tar xf gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    mv gffcompare-0.12.6.Linux_x86_64/gffcompare /usr/local/bin/ && \
    rm -rf gffcompare-0.12.6.Linux_x86_64*

# Instalación gffread
RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar xf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin/ && \
    rm -rf gffread-0.12.7.Linux_x86_64*

# Instalación Miniconda 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.11.0-2-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Configuracion e instalacion de conda en el ambiente
ENV PATH="/opt/conda/bin:$PATH"
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n rnaseq_env python=3.9 && \
    echo "source activate rnaseq_env" > ~/.bashrc

# Instalación de los paquetes de conda en el ambiente
RUN conda Instalación -n rnaseq_env -y -c bioconda -c conda-forge \
    multiqc=1.14 \
    fastqc=0.11.9 \
    suppa=2.3 \
    deeptools \
    pandas=1.5.3 \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    biopython \
    scikit-learn \
    statsmodels \
    && conda clean -a

# Instalación rMATS v4.3.0
RUN git clone https://github.com/Xinglab/rmats-turbo.git /opt/rmats && \
    cd /opt/rmats && \
    git checkout v4.3.0 && \
    /opt/conda/envs/rnaseq_env/bin/pip Instalación numpy cython && \
    ./build_rmats && \
    chmod +x /opt/rmats/rmats.py && \
    ln -s /opt/rmats/rmats.py /usr/local/bin/rmats.py

# Instalación de los paquetes de python en conda
RUN /opt/conda/envs/rnaseq_env/bin/pip Instalación \
    pyarrow \
    openpyxl

# Configuracion del ambiente
ENV PATH=$PATH:/opt/GSEA_4.3.3:/opt/FastQC:/opt/rmats/rMATS-turbo-Linux-UCS4
ENV TRIMMOMATIC_HOME=/opt/Trimmomatic-0.39
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64

# Configuracion de Tximport
COPY tximport.R /usr/local/bin/
RUN chmod +x /usr/local/bin/tximport.R

# Verfica que conda este en el path
ENV PATH="/opt/conda/envs/rnaseq_env/bin:$PATH"
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate rnaseq_env" >> ~/.bashrc

# Crea el link de multiqc
RUN ln -s /opt/conda/envs/rnaseq_env/bin/multiqc /usr/local/bin/multiqc

# Indica el directorio de trabajo
WORKDIR /data

# Verifica las instalaciones
RUN R -e 'library(DESeq2); packageVersion("DESeq2")' && \
    fastqc --version && \
    bedtools --version && \
    samtools --version && \
    gffread --version && \
    STAR --version && \
    salmon --version && \
    multiqc --version && \
    gffcompare --version && \
    stringtie --version && \
    python3 -c "import scipy; print('scipy version:', scipy.__version__)" && \
    python3 -c "import pandas; print('pandas version:', pandas.__version__)" && \
    bamCoverage --version && \
    rmats.py --version && \
    suppa.py -h && \
    [ -f /usr/local/bin/tximport.R ] && echo "tximport.R installed successfully"

CMD ["/bin/bash"]
