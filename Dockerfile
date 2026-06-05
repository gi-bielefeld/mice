FROM alpine:3.22 AS builder

RUN apk add --no-cache tar gzip wget cargo git alpine-sdk zlib-dev openmp-dev libgomp
WORKDIR /opt/software

RUN git clone https://github.com/gi-bielefeld/mice.git
RUN cargo install --path ./mice --root .

RUN wget "https://github.com/algbio/ggcat/releases/download/v2.0.0/ggcat-x86_64-unknown-linux-musl.tar.gz"; \
tar xvzf ggcat-x86_64-unknown-linux-musl.tar.gz

RUN git clone https://github.com/lucaparmigiani/gfa2gff
RUN cd gfa2gff && make

FROM alpine:3.22

RUN apk add --no-cache bash gcompat libgcc zlib openmp libgomp
WORKDIR /opt/software

COPY --from=builder /opt/software/bin/mice /opt/software/mice
COPY --from=builder /opt/software/ggcat /opt/software/ggcat
COPY --from=builder /opt/software/gfa2gff/gfa2gff /opt/software/gfa2gff

ENV PATH="/opt/software:${PATH}"