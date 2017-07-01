FROM continuumio/anaconda3

WORKDIR /home

RUN apt-get update \
    && apt-get install -y git make zip graphviz inkscape texlive-xetex texlive-fonts-recommended texlive-latex-extra librsvg2-bin docbook2x \
    && apt-get -y clean

RUN conda config --add channels conda-forge
RUN conda config --set always_yes yes
RUN conda install python=3.6 rever requests requests-oauthlib \
    && conda clean --all
RUN conda info
RUN /opt/conda/bin/pip install xonda

RUN git clone git://github.com/sympy/sympy.git

WORKDIR /home/sympy/release

ENTRYPOINT ["/opt/conda/bin/rever"]
