FROM ubuntu:jammy-20221130

# GHCR related label
LABEL org.opencontainers.image.source=https://github.com/id-bioinfo/TIPars
LABEL org.opencontainers.image.description="TIPars - Taxa Insertion by Parsimony on Ubuntu Jammy"
LABEL org.opencontainers.image.licenses=LGPL-2.1

# Copy all files to /tipars/
COPY ./ /tipars/

# Update apt and apt-get
RUN apt-get update
RUN apt update
RUN apt upgrade -y

# Install java JRE
RUN apt install default-jre -y

# Setting timezone for installing python without manual user input
ENV TZ=Asia/Hong_Kong
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# Install deadsnakes for python installation
RUN apt install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa -y
# Install python3.11
RUN apt install python3.11 -y
# Install pip3
RUN apt install python3-pip -y
# Update python3 command to use python3.11 instead of 3.10
RUN update-alternatives --install /usr/bin/python3 python /usr/bin/python3.11 10

# Install Perl and Perl libraries 
RUN apt install perl -y
RUN apt-get install libtime-hires-perl -y
RUN apt install libdata-dump-perl -y

# Install ete3 w/ pip
RUN python3.11 -m pip install six numpy
RUN python3.11 -m pip install ete3

# Install pastml w/ pip
RUN python3.11 -m pip install pastml

# Install gcc and OpenMP
RUN apt install gcc -y
RUN apt-get install libomp-dev -y

WORKDIR /home

# Default CMD: Tipars toy test
CMD ["/tipars/tipars", "-t", "/tipars/Benchmark_datasets/NDV/NDV_tree.nwk", "-s", "/tipars/Benchmark_datasets/NDV/NDV_taxa.fas", "-a", "/tipars/Benchmark_datasets/NDV/NDV_anc.fas", "-q", "/tipars/Benchmark_datasets/NDV/NDV_query.fas", "-o", "/tipars/Benchmark_datasets/NDV/tipars.tree"]