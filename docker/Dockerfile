FROM ubuntu:focal

ARG BUILD_BRANCH
ARG READING_BRANCH

# Set working folder
ENV DIRPATH /sw
WORKDIR $DIRPATH

RUN apt-get update && \
    # Install Java
    apt-get install -y openjdk-8-jdk && \
    # jnius-indra requires cython which requires gcc
    apt-get install -y git wget zip unzip bzip2 gcc graphviz graphviz-dev \
        pkg-config python3 python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python

# Set default character encoding
# See http://stackoverflow.com/questions/27931668/encoding-problems-when-running-an-app-in-docker-python-java-ruby-with-u/27931669
# See http://stackoverflow.com/questions/39760663/docker-ubuntu-bin-sh-1-locale-gen-not-found
RUN apt-get install -y locales && \
    locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8  #

# Set environment variables
ENV BNGPATH=$DIRPATH/BioNetGen-2.4.0
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

# Install INDRA and dependencies
RUN git clone https://github.com/sorgerlab/indra.git && \
    cd indra && \
    git checkout $BUILD_BRANCH && \
    git branch && \
    mkdir /root/.config && \
    mkdir /root/.config/indra && \
    echo "[indra]" > /root/.config/indra/config.ini && \
    # Install Python dependencies
    pip install --upgrade pip && \
    # Install cython first for pyjnius
    pip install cython && \
    pip install -e .[all] && \
    pip uninstall -y enum34 && \
    # Pre-build the bio ontology
    python -m indra.ontology.bio build && \
    # Download Adeft models
    python -m adeft.download && \
    # Download protmapper resources
    python -m protmapper.resources && \
    # Install BioNetGen
    cd $DIRPATH && \
    wget "https://github.com/RuleWorld/bionetgen/releases/download/BioNetGen-2.4.0/BioNetGen-2.4.0-Linux.tgz" \
        -O bionetgen.tar.gz -nv && \
    tar xzf bionetgen.tar.gz && \
    # Install things related to API deployment
    pip install gunicorn

# Install indra_reading
RUN git clone https://github.com/indralab/indra_reading.git && \
    cd indra_reading && \
    git checkout $READING_BRANCH && \
    echo $READING_BRANCH && \
    pip install -e .

# Add and set up reading systems
# ------------------------------
# SPARSER
ENV SPARSERPATH=$DIRPATH/sparser
ADD r3.core $SPARSERPATH/r3.core
ADD save-semantics.sh $SPARSERPATH/save-semantics.sh
ADD version.txt $SPARSERPATH/version.txt
RUN chmod +x $SPARSERPATH/save-semantics.sh && \
    chmod +x $SPARSERPATH/r3.core

# REACH
# Default character encoding for Java in Docker is not UTF-8, which
# leads to problems with REACH; so we set option
# See https://github.com/docker-library/openjdk/issues/32
ENV JAVA_TOOL_OPTIONS -Dfile.encoding=UTF8
ENV REACHDIR=$DIRPATH/reach
ENV REACHPATH=$REACHDIR/reach-1.6.3-e48717.jar
ENV REACH_VERSION=1.6.3-e48717
ADD reach-1.6.3-e48717.jar $REACHPATH

# MTI
ADD mti_jars.zip $DIRPATH
RUN mkdir $DIRPATH/mti_jars && \
    unzip $DIRPATH/mti_jars.zip -d $DIRPATH/mti_jars/
