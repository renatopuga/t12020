
FROM gitpod/workspace-full

RUN sudo apt-get update -y 
RUN sudo DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN sudo apt-get install -y fastqc
RUN sudo apt-get install -y samtools
RUN sudo apt-get install -y bwa
RUN sudo apt-get install -y freebayes

# set working directory
WORKDIR /app
CMD ["fastqc", "samtools", "freebayes", "bwa"]
