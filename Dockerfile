# syntax=docker/dockerfile:1
FROM debian:11.1

RUN apt-get update && apt-get install -y \
    bison \
    cmake \
    curl \
    flex \
    g++ \
    libaec-dev \
    libboost-all-dev \
    libc6-dev \
    libcunit1-dev \
    libcurl4-openssl-dev \
    libfreetype6-dev \
    libgeos-dev \
    libgsl-dev \
    libopenmpi-dev \
    libproj-dev \
    libsnappy-dev \
    libudunits2-dev \
    libzip-dev \
    nlohmann-json3-dev \
    pkg-config \
    proj-bin \
    python3-dev \
    python3-pip \
    python3-tk \
    && rm -rf /var/lib/apt/lists/*

# Build hdf5
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /tmp
ENV HDF5_USE_FILE_LOCKING=FALSE
RUN curl -SL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.bz2 \
    | tar -xjC /tmp && \
    cd hdf5-1.12.1 && \
    ./configure --enable-parallel --enable-build-mode=production --enable-hl --enable-tools --prefix=/usr/local &&\
    make -j4 &&\
    make install && \
    rm -r /tmp/hdf5-1.12.1

# Build pnetcdf
WORKDIR /tmp
ENV CC=mpicc
RUN curl -SL http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/pnetcdf-1.12.2.tar.gz \
    | tar -xzC /tmp && \
    cd pnetcdf-1.12.2 && \
    ./configure --disable-fortran && \
    make -j4 && \
    make install && \
    rm -r /tmp/pnetcdf-1.12.2

# Build netcdf
WORKDIR /tmp
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
RUN curl -SL https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.8.1.tar.gz \
    | tar -xzC /tmp && \
    cd netcdf-c-4.8.1 && \
    CC=mpicc ./configure --disable-shared --enable-parallel-tests --enable-utilities --prefix=/usr/local && \
    make -j4 && \
    make install && \
    rm -r /tmp/netcdf-c-4.8.1

# install codipack
WORKDIR /tmp
ENV CC=gcc
RUN curl -SL https://github.com/SciCompKL/CoDiPack/archive/refs/tags/v2.0.2.tar.gz \
    | tar -xzC /tmp && \
    mv /tmp/CoDiPack-2.0.2/include/* /usr/include/. && \
    rm -r /tmp/CoDiPack-2.0.2

# install doxygen
WORKDIR /tmp
RUN curl -SL https://www.doxygen.nl/files/doxygen-1.9.6.src.tar.gz \
    | tar -xzC /tmp && \
    cd doxygen-1.9.6 && \
    mkdir build && \
    cd build && \
    cmake -j4 -G "Unix Makefiles" .. && \
    make -j 4 && \
    make install && \
    rm -r /tmp/doxygen-1.9.6


# install all python dependencies
WORKDIR /app
COPY requirements.txt .
RUN python3 --version && pip3 install -r requirements.txt