FROM ubuntu:16.04

RUN apt-get -y update 
RUN apt-get install -y git python3 python3-pip gcc psutils cmake wget make
RUN apt-get install -y gfortran libblas-dev liblapack-dev python3-tk
RUN pip3 install cython matplotlib pytest scipy psutil pyvtk ipywidgets

RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /usr/local
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR /usr/local/fidimag/bin
RUN bash install-fftw.sh
RUN bash install-sundials.sh

ENV PYTHONPATH="/usr/local/fidimag:$PYTHONPATH"
ENV LD_LIBRARY_PATH="/usr/local/fidimag/local/lib:$LD_LIBRARY_PATH"

WORKDIR /usr/local/fidimag
RUN python3 setup.py build_ext --inplace

ENV OMP_NUM_THREADS=2
WORKDIR /io
