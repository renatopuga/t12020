FROM gitpod/workspace-full

# Install custom tools, runtimes, etc.
# For example "bastet", a command-line tetris clone:
# RUN brew install bastet
#
# More information: https://www.gitpod.io/docs/config-docker/

RUN sudo apt-get update -y \
  && sudo DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata \ 
  && sudo apt-get install -y fastqc samtools bwa freebayes

RUN sudo apt-get install -y cpanm \
    && sudo cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib) \
    && sudo cpanm install DBI Spread Module::Build Try::Tiny DBD::mysql

RUN git clone https://github.com/Ensembl/ensembl-vep.git \
    && cd ensembl-vep \ 
    && git pull \
    && git checkout release/101 \
    && perl INSTALL.pl --NO_HTSLIB --AUTO ap --PLUGINS gnomADc

RUN sudo mkdir /opt/vep/ \
    && sudo mkdir /opt/vep/.vep \
    && sudo chmod 777 /opt/vep -R 
  
CMD ["fastqc", "samtools", "freebayes", "bwa"]
