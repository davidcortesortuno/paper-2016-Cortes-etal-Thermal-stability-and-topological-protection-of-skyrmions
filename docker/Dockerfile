FROM ubuntu:16.04

RUN apt-get -y update 
RUN apt-get install -y git python3 python3-pip gcc psutils cmake wget make
RUN apt-get install -y gfortran libblas-dev liblapack-dev python3-tk sudo fonts-lato
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
RUN python3 -c "import matplotlib"
# Headless Matplotlib:
ENV MPLBACKEND=Agg

ENV OMP_NUM_THREADS=2
WORKDIR /io

# Set up user so that we do not run as root
# (not necessary if we run docker with --user)
# -s for shell, -m to create home directory, -G groups
#
# RUN useradd -m -s /bin/bash -G sudo fidimag && \
#     echo "fidimag:docker" | chpasswd && \
#     echo "fidimag ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
# # -R recursive
# RUN chown -R fidimag:fidimag /usr/local/fidimag

