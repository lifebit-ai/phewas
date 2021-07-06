FROM quay.io/lifebitai/phewas:latest
# Copy additonal scripts
COPY bin/* /opt/bin/
RUN chmod +x /opt/bin/*
ENV PATH="$PATH:/opt/bin/"
