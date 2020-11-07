# pipeline todos os comandos

#----- Estrutura de Diretórios

# voltar para workspace
cd /workspace/t12020

# criar o diretorio no home (caso ele ainda nao exista)
mkdir bioinfo/
mkdir bioinfo/app
mkdir bioinfo/resultados
mkdir bioinfo/reference
mkdir bioinfo/data
mkdir bioinfo/data/fastq

#----- Diretório de Resultados

# voltar para workspace
cd /workspace/t12020

# entrar no diretorio resultados
cd bioinfo/resultados/

# mkdir para criar o diretorio das amostras 003, 017 e 019
mkdir 003 017 019

#---- (gdown.pl) Script para download de arquivos no Google Drive 

# voltar para workspace
cd /workspace/t12020

# clonar repositorio
git clone https://github.com/circulosmeos/gdown.pl.git

#----- Download dos Arquivos .FASTQ (no Google Drive) direto do terminal

# voltar para workspace
cd /workspace/t12020

# rodar o script ./gdown.pl download
./gdown.pl/gdown.pl https://drive.google.com/open?id=1LlBQ2BNI_vF-E_4ntrEexJQ2-Wea2e0o 003.fastq.gz
./gdown.pl/gdown.pl https://drive.google.com/open?id=1KKmk90gUk0174MvATZzWRwyqRSHJNZDv 017.fastq.gz
./gdown.pl/gdown.pl https://drive.google.com/open?id=11jtliN2G0vTs50-z79QxXLWkStPeuh2E 019.fastq.gz

# movendo arquivos fastq para outro diretorio
mv 003.fastq.gz 017.fastq.gz 019.fastq.gz  bioinfo/data/fastq

#---- Download dos Cromossomos (chr13 e chr17) da UCSC
# volta para workspace
cd /workspace/t12020

# entrar em referencias
cd bioinfo/reference

# wget para fazer download do chr13
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz

# wget para fazer download do chr17 
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz

# cat para concatenar os arquivo do cromossomo 13 e 17 em hg19.fa
zcat chr13.fa.gz chr17.fa.gz > hg19.fa

# rm para deletar os arquivos chr13.fa e chr17.fa
rm chr13.fa.gz chr17.fa.gz

#----- Annovar: instalação e download dos bancos (exac e clinvar) 

# install dabases annovar
# NOTA: entreno no site do ANNOVAR com seu e-mail e salve o arquivo no diretorio: bionfo/app/

cd /workspace/t12020

# download do arquivo annovar.latest.tar.gz
./gdown.pl/gdown.pl https://drive.google.com/file/d/1XVnRT0GUKFuQifvoROgHMqoy4aV3PCGX/view?usp=sharing annovar.tar.gz

# mover annovar
mv annovar.tar.gz bioinfo/app

# entrar no diretorio
cd bioinfo/app/

# descompactar 
tar -zxvf annovar.tar.gz

# entrar no diretorio annovar
cd annovar

# baixar as bases: clinvar e exac03
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20200316  humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/

#----- BWA index

# voltar para workspace
cd /workspace/t12020

# entrar em reference
cd bioinfo/reference

# bwa index para gerar o index da referencia hg19.fa
bwa index hg19.fa 

# voltar para workspace
cd /workspace/t12020

#----- FastQC para gerar relatorio de qualidade dos arquivo FASTQ
# opcao (-o ./) diz para salvar o resultado no mesmo arquivo em que o comando esta sendo rodado.

# voltar para workspace
cd /workspace/t12020

# rodar fastqc e salvar o resultado de cada amostra em seu diretorio
fastqc -o bioinfo/resultados/003/ bioinfo/data/fastq/003.fastq.gz


#----- cutadapt: instalação e run

# instalar cutadapt usando pip
python3 -m pip install --user --upgrade cutadapt

# voltar para workspace
cd /workspace/t12020

# cutadapt para remover sequencias de tamanho menores 100pb e maiores que 220
cutadapt --minimum-length 100 --maximum-length 220 -q 15  -o bioinfo/resultados/003/003.cutadapt.fastq  bioinfo/data/fastq/003.fastq.gz


#------ BWA-mem: maximal exact matches

# voltar para workspace
cd /workspace/t12020

# rodar bwa para alinhar as sequencias contra o genoma de referencia
bwa mem -R '@RG\tID:003\tSM:003_NGSA\tLB:XPTO\tPL:Ion'  bioinfo/reference/hg19.fa bioinfo/resultados/003/003.cutadapt.fastq > bioinfo/resultados/003/003.sam


#------ Samtools: fixmate e sort

# voltar para workspace
cd /workspace/t12020

# fixmate: converter sam to bam
samtools fixmate bioinfo/resultados/003/003.sam bioinfo/resultados/003/003.bam

# sort: ordenar pela coordenadas genomicas
samtools sort -O bam -o bioinfo/resultados/003/003_sort.bam -T /tmp/ bioinfo/resultados/003/003.bam

# index: gerar o index (.bai) do arquivo bam
samtools index bioinfo/resultados/003/003_sort.bam 

#------ Samtools: freebayes

# voltar para workspace
cd /workspace/t12020

# freebayes: chamar variantes
# NOTA: ler sobre os parametros no github
freebayes -f bioinfo/reference/hg19.fa -F 0.01 -C 1 --pooled-continuous bioinfo/resultados/003/003_sort.bam > bioinfo/resultados/003/003.vcf


#------ Anotacao: 1/2 - converter VCF para AVINPUT

# voltar para workspace
cd /workspace/t12020

perl bioinfo/app/annovar/convert2annovar.pl -format vcf4 bioinfo/resultados/003/003.vcf > bioinfo/resultados/003/003.avinput

# comando head para listar as 10 primeiras linhas do arquivo 003.avinput
head bioinfo/resultados/003/003.avinput

#------ Anotacao: 2/2 - Adicionar informacoes do refeseq, exac e clinvar

# voltar para workspace
cd /workspace/t12020

perl bioinfo/app/annovar/table_annovar.pl bioinfo/resultados/003/003.avinput bioinfo/app/annovar/humandb/ -buildver hg19 -out bioinfo/resultados/003/003 -remove -protocol refGene,exac03,clinvar_20200316 -operation g,f,f -nastring . 

# filtrar o resultado (_multianno) com o comando grep.
# Pegar o cabecalho e apenas linhas do arquivo com a palavra exonic
grep "^Chr\|exonic" bioinfo/resultados/003/003.hg19_multianno.txt  | head
