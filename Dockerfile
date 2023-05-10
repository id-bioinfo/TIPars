FROM centos:7
WORKDIR /code
COPY . /code/

# update yum
RUN yum update -y

# install wget
RUN yum install -y wget

# install perl
RUN yum install -y perl 
RUN rpm -qa
RUN yum install -y 'perl(Data::Dumper)'

# install conda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda env create -f ete3-py2-env.yml 
RUN conda create --name pastml-py3 python=3.6
# RUN source activate ete3-py2 && conda install -c etetoolkit ete3 ete_toolchain && source deactivate
RUN source activate pastml-py3 && pip install -r pastml-py3-req.txt && conda deactivate

# install java 11
RUN yum install -y java-11-openjdk-devel
RUN yum install -y java-11-openjdk
RUN java -version

# install nextflow
RUN yum install -y which
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
RUN ./nextflow self-update