FROM lifebitai/phewas@sha256:359d5028fb441b6a014137b64c5eaca1a0409c8a24ed1ff3e3b7d90dc07631cc
# See OG Dockerfile here: https://github.com/lifebit-ai/phewas/blob/8352b4c8c26ca086862eba0b4cd2d5c694060d71/docker/Dockerfile

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y procps && apt-get clean -y