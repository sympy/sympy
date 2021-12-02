FROM continuumio/anaconda3

WORKDIR /root

RUN apt-get update \
    && apt-get install -y libc6-i386 libc6 linux-headers-amd64 git make zip graphviz chromium-browser texlive-xetex texlive-fonts-recommended texlive-latex-extra librsvg2-bin docbook2x latexmk libatlas-dev libatlas-base-dev liblapack-dev gfortran gcc dbus \
    && apt-get -y clean

RUN conda config --add channels conda-forge
RUN conda config --set always_yes yes
RUN conda install python=3.7 rever requests requests-oauthlib \
    && conda clean --all
RUN conda info
RUN /opt/conda/bin/pip install xonda

# Make matplotlib tests work
# https://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
ENV MATPLOTLIBRC=/root/matplotlibrc
RUN mkdir -p $MATPLOTLIBRC
RUN echo "backend : Agg" > $MATPLOTLIBRC/matplotlibrc

RUN git clone git://github.com/sympy/sympy.git

WORKDIR /root/sympy/release

ENTRYPOINT ["./pull_and_run_rever.sh"]
