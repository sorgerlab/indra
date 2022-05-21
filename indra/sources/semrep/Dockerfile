FROM ubuntu:latest

ARG UMLS_API_KEY

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

WORKDIR /sw

RUN ln -sf /bin/bash /bin/sh && \
    apt-get update && \
    # Install Java
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y git bzip2 gcc pkg-config python3 python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    pip install umls_downloader

RUN cd /sw && \
    umls_downloader custom --url https://data.lhncbc.nlm.nih.gov/umls-restricted/ii/tools/SemRep_SemMedDB_SKR/public_semrep_v1.8_Linux.tar.bz2 \
        --output public_semrep_v1.8_Linux.tar.bz2 && \
    umls_downloader custom --url https://data.lhncbc.nlm.nih.gov/umls-restricted/ii/tools/SemRep_SemMedDB_SKR/public_semrep_v1.8_lex_db_18_Linux.tar.bz2 \
        --output public_semrep_v1.8_lex_db_18_Linux.tar.bz2 && \
    umls_downloader custom --url https://data.lhncbc.nlm.nih.gov/umls-restricted/ii/tools/SemRep_SemMedDB_SKR/public_semrep_v1.8_hier_18_Linux.tar.bz2 \
        --output public_semrep_v1.8_hier_18_Linux.tar.bz2 && \
    umls_downloader custom --url https://data.lhncbc.nlm.nih.gov/umls-restricted/ii/tools/MetaMap/download/public_mm_linux_main_2016v2.tar.bz2 \
        --output public_mm_linux_main_2016v2.tar.bz2

RUN cd /sw && \
    tar xvjf public_semrep_v1.8_Linux.tar.bz2 && \
    tar xvjf public_semrep_v1.8_lex_db_18_Linux.tar.bz2 && \
    tar xvjf public_semrep_v1.8_hier_18_Linux.tar.bz2 && \
    tar xvjf public_mm_linux_main_2016v2.tar.bz2

RUN cd /sw/public_mm && \
    ./bin/install.sh && \
    cd /sw/public_semrep && \
    ./bin/install.sh \
    cd /sw/public_semrep/lib && \
    # This symlink is missing and is required for the system to work
    ln -s libpcre.so.0.0.1 libpcre.so.1

