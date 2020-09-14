
# ---- Base Node ----
FROM ubuntu AS base
# install the core dependencies
RUN apt-get update -y 
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get install -y fastqc
RUN apt-get install -y samtools
RUN apt-get install -y bwa
RUN apt-get install -y freebayes

# set working directory
WORKDIR /app
CMD ["fastqc", "samtools", "freebayes", "bwa"]
