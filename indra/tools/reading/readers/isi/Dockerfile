FROM sahilgar/bigmechisi

ARG INDRA_USER
ARG INDRA_BRANCH
RUN echo $INDRA_USER $INDRA_BRANCH

ARG INDRA_DB_USER
ARG INDRA_DB_BRANCH
RUN echo $INDRA_DB_USER $INDRA_DB_BRANCH

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
RUN echo "Cloning" https://github.com/$INDRA_USER/indra.git && \
    git clone https://github.com/$INDRA_USER/indra.git && \
    cd indra && \
    echo "================================================" && \
    git fetch --all && \
    git branch --list --all && \
    echo "Changing to branch" $INDRA_BRANCH &&\
    git checkout $INDRA_BRANCH && \
    git status && \
    git show && \
    echo "================================================" && \
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
RUN echo "Cloning" https://github.com/$INDRA_DB_USER/indra_db.git && \
    git clone https://github.com/$INDRA_DB_USER/indra_db.git && \
    cd indra_db && \
    echo "================================================" && \
    git fetch --all && \
    git branch --list --all && \
    echo "Changing to branch" $INDRA_DB_BRANCH &&\
    git checkout $INDRA_DB_BRANCH && \
    git status && \
    git show && \
    echo "================================================" && \
    pip3 install . && \
    cd ..

# Check the installation (and make sure any files are downloaded)
RUN python3 -m indra_db.reading.read_db_aws --help

# Set the in-isi environment variable
ENV IN_ISI_DOCKER=true
