FROM quay.io/pypa/manylinux1_x86_64:latest

#Taken from https://github.com/skvark/opencv-python/blob/master/docker/Dockerfile_x86_64
RUN curl -O -L https://cmake.org/files/v3.9/cmake-3.9.0.tar.gz && \
    tar -xf cmake-3.9.0.tar.gz && \
    cd cmake-3.9.0 && \
    #manylinux1 provides curl-devel equivalent and libcurl statically linked
    # against the same newer OpenSSL as other source-built tools
    # (1.0.2s as of this writing)
    yum -y install zlib-devel && \
    #configure does a bootstrap make under the hood
    export MAKEFLAGS=-j$(getconf _NPROCESSORS_ONLN) && \
    ./configure --system-curl && \
    make && \
    make install && \
    cd .. && \
    rm -rf cmake-3.9.0*

# https://trac.ffmpeg.org/wiki/CompilationGuide/Centos#GettheDependencies
# manylinux provides the toolchain and git; we provide cmake
RUN yum install -y  lapack-devel  blas-devel


ENV PATH "$HOME/bin:$PATH"
