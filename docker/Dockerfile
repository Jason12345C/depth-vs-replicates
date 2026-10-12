# Pin this base tag too (check hub.docker.com/r/mambaorg/micromamba/tags)
FROM mambaorg/micromamba:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
