FROM ubuntu:22.04
COPY ./bin/* /usr/bin/
ENV PATH="/usr/bin/:${PATH}"

RUN apt-get update && apt-get install -y \
    zlib1g \
    libc6 \
    libcurl4 \
    libssl3 \
    libnghttp2-14 \
    libidn2-0 \
    librtmp1 \
    libssh-4 \
    libpsl5 \
    libkrb5-3 \
    libldap-2.5-0 \
    libzstd1 \
    libbrotli1 \
    libunistring2 \
    libgnutls30 \
    libhogweed6 \
    libnettle8 \
    libgmp10 \
    libk5crypto3 \
    libcom-err2 \
    libkrb5support0 \
    libsasl2-2 \
    libp11-kit0 \
    libtasn1-6 \
    libkeyutils1 \
    libffi8 \
    libgssapi-krb5-2 \
    libkrb5-3

