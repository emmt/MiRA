FROM ferreol/yeti:latest

MAINTAINER Ferréol Soulez <ferreol.soulez@univ-lyon1.fr>

RUN echo http://nl.alpinelinux.org/alpine/edge/testing >> /etc/apk/repositories && \
    apk add --update --no-cache --virtual .build-deps  bash  make wget build-base zip  &&\
    apk add --update  --no-cache fftw-dev  &&\
    mkdir -p /opt && cd /opt && \
    wget https://github.com/emmt/OptimPackLegacy/archive/master.zip && \
    unzip master.zip && rm master.zip && cd OptimPackLegacy-master/ && \
    yorick/configure && make && make install && \
    cd /opt && rm -rf OptimPackLegacy-master && \
    wget https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.4.1.tar.gz && \
    tar -xvzf nfft-3.4.1.tar.gz && rm nfft-3.4.1.tar.gz && cd nfft-3.4.1/ && \
    ./configure --enable-openmp  --enable-all && make && make install && \
    cd .. && rm -rf  nfft-3.4.1 && \
    wget https://github.com/emmt/ynfft/archive/master.zip && \
    unzip master.zip && rm master.zip && cd ynfft-master/ && \
    ./configure  --cflags='-I/usr/local/include -fopenmp -Ofast -march=native -mfpmath=sse -pedantic -pipe -std=c99' && \
    make && make install && cd /opt && rm -rf ynfft-master && \
    wget -q  https://github.com/emmt/MiRA/archive/master.zip && \ 
    unzip master.zip && rm master.zip && cd MiRA-master/ && \ 
    ./configure --prefix=/usr/local && make install && \
    cd /opt && rm -rf /opt/MiRA-master && \
    apk del .build-deps  && \
    rm -rf /var/cache/apk/* && rm -rf /usr/local/share/doc

ENTRYPOINT ["/usr/local/bin/ymira"]
