
FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y gcc g++ make

COPY src /app/src

WORKDIR /app/src

RUN rm -f /app/src/freec && make && chmod +x /app/src/freec

ARG USER_ID=1000
ARG GROUP_ID=1000

RUN groupadd --gid $GROUP_ID -r user && useradd -r --uid $USER_ID --gid $GROUP_ID user
USER user

ENTRYPOINT ["./freec", "-conf"]
CMD ["/app/data/test/config_BL.txt"]