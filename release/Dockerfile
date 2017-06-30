FROM sympy-release-base

WORKDIR /home/sympy/release

# TODO: How can we keep from hard-coding this?
RUN git checkout 1.1
RUN git pull

ENTRYPOINT ["/opt/conda/envs/sympy-release/bin/rever"]
