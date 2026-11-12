FROM continuumio/miniconda3

COPY beyondcell_environment.yml /tmp/beyondcell_environment.yml

RUN conda env create -f /tmp/beyondcell_environment.yml && \
    conda clean -afy && \
    conda config --system --set auto_activate_base false

RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate beyondcell" >> /root/.bashrc

ENV PATH=/opt/conda/envs/beyondcell/bin:$PATH

WORKDIR /workspace

CMD ["/bin/bash"]
