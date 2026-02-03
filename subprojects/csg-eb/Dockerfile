# Docker Image for running csg-eb tests

FROM ubuntu:20.04

LABEL maintainer="Mark Meredith <mark.meredith@netl.doe.gov>"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update \
  && apt-get -qq -y install --no-install-recommends \
  build-essential=12.* \
  git=1:2.* \
  file=1:5.* \
  openscad=2019.* \
  pkg-config=0.29.* \
  python3-pip=20.* \
  python3-venv=3.8.* \
  libgmp-dev \
  libmpfr-dev \
  libboost-all-dev \
  && apt-get clean && rm -rf /var/lib/apt/lists/* \

  && pip3 install --no-cache-dir \
  cmake==3.* \
  ninja==1.* \

  && git clone -c feature.manyFiles=true https://github.com/spack/spack.git /usr/local/spack \

  && PATH=/usr/local/spack/bin:$PATH

RUN git clone --depth 1 --branch v2.13.7 https://github.com/catchorg/Catch2
WORKDIR /Catch2
RUN cmake -S . -B build -DCATCH_BUILD_TESTING=OFF
RUN cmake --build build --target install
WORKDIR /

RUN git clone --depth 1 --branch v5.3 https://github.com/CGAL/cgal
WORKDIR /cgal
RUN cmake -S . -B build
RUN cmake --build build --target install
WORKDIR /

RUN git clone --branch 3.2.2 https://github.com/taocpp/PEGTL
WORKDIR /PEGTL
RUN cmake -S . -B build -DPEGTL_BUILD_TESTS=OFF -DPEGTL_BUILD_EXAMPLES=OFF
RUN cmake --build build --target install
WORKDIR /

RUN  useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
