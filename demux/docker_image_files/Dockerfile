# Use a base Ubuntu image
FROM ubuntu:22.04


# Update package lists and install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    xorg-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    gfortran \
    libreadline-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libharfbuzz-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgc-dev \
    python3-pip \
    less \
    grep \
    vim 

WORKDIR /tmp/

# Copying the demux repo and sample data

COPY ./16S-demux /16S-demux
COPY ./test_fastq_data /test_fastq_data
COPY ./docker_image_files/requirements.txt ./requirements.txt

# Install Python requirements
RUN pip3 install -r ./requirements.txt 

WORKDIR /16S-demux
