FROM sahilgar/bigmechisi

ARG INDRA_USER
ARG INDRA_BRANCH

ARG INDRA_DB_USER
ARG INDRA_DB_BRANCH

ENV DIRPATH /root
ENV PYTHONPATH "${DIRPATH}/indra_db:${PYTHONPATH}"
WORKDIR $DIRPATH

RUN apt-get update && \
    apt-get install -y git wget python3 python3-pip

# Install libpq5
RUN apt-get update && \
    apt-get install -y libpq5 libpq-dev

# Make sure pip is up-to-date
RUN pip3 install -U pip

# Install indra
RUN git clone https://github.com/$INDRA_USER/indra.git && \
    cd indra && \
    git checkout $INDRA_BRACH && \
    pip3 install . && \
    cd ..

# Update the protmapper
RUN pip3 install -U git+https://github.com/indralab/protmapper.git

# Install psycopg2
RUN git clone https://github.com/psycopg/psycopg2.git && \
    cd psycopg2 && \
    python3 setup.py build && \
    python3 setup.py install

# Install pgcopy
RUN git clone https://github.com/pagreene/pgcopy.git && \
    cd pgcopy && \
    python3 setup.py install

# Install adeft
RUN pip3 install adeft
RUN python3 -m adeft.download

# Install indra_db
RUN git clone https://github.com/$INDRA_DB_USER/indra_db.git && \
    cd indra_db && \
    git checkout $INDRA_DB_BRANCH && \
    pip3 install . && \
    cd ..

# Check the installation (and make sure any files are downloaded)
RUN python3 -m indra_db.reading.read_db_aws --help
