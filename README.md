[![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/renatopuga/t12020)

# t12020 - Meu Primeiro Pipeline de NGS

# Olá Ambiente de Bioinformática
***Linux:*** O Linux é um sistema operacional, assim como o Windows da Microsoft e o Mac OS da Apple. Ele foi criado pelo finlandês Linus Torvalds, e o nome é a mistura do nome do criador com Unix, um antigo sistema operacional da empresa de mesmo nome. [O que é Linux?](https://www.techtudo.com.br/artigos/noticia/2011/12/o-que-e-linux.html).

***Pipelines:*** Primeiro, o pipeline não é um termo de bioinformática, é na verdade um termo de ciência da computação, envolve encadeamento de processos / threads / funções etc. Em resumo, o resultado de um processo é a entrada de um novo processo na sequência. [A review of bioinformatic pipeline frameworks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429012/)

# Qual ambiente vamos utilizar?
Vamos utilizar uma instância na Amazon Cloud com Sistema Operacional Linux Ubuntu. Os programas para que possamos criar nosso pipeline de análises de sequenciamento de nova geração já foram instalados e serão apenas listados nesse documento. O objetivo é melhorar a experiência de usuários iniciantes em "digitar comandos em um terminal linux" e ser um roteiro para que todos possam se localizar.

# Programas Utilizados
Listamos os programas previamente instalados em nosso ambiente para executar o pipeline de chamada de varinates:

* annovar [Download](http://annovar.openbioinformatics.org/)
* bwa [Download](http://bio-bwa.sourceforge.net/)
* FastQC [Download](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* freebayes [Download](https://github.com/ekg/freebayes)
* samtools [Download](http://samtools.sourceforge.net/)


# gitbash

* Para Windows - Terminal
	https://github.com/git-for-windows/git/releases/download/v2.28.0.windows.1/PortableGit-2.28.0-32-bit.7z.exe
	
* Para Macbook: 
	Buscar no Spotlight (buscar por Terminal - dois cliques)
	
* Para Linux: 
	Procure por Terminal


# Primeiros Comandos

```bash
# ls para listar arquivos e diretorios
ls

# pwd para informar o diretorio atual
pwd
```

Fonte: [Comandos Básicos do Terminal Linux](http://swcarpentry.github.io/shell-novice/)

# Estrutura de Diretórios
Estrutura de diretórios: sequências, programas e arquivos de referência.

* bioinfo
	* bioinfo/app
	* bioinfo/data
		* bioinfo/data/fastq
	* bioinfo/reference
  * bioinfo/resultados


# Criando a Estrutura de Diretórios

```bash

# volta para workspace
cd /workspace/t12020

# criar o diretorio no home (caso ele ainda nao exista)
mkdir bioinfo/
mkdir bioinfo/app
mkdir bioinfo/resultados
mkdir bioinfo/reference
mkdir bioinfo/data
mkdir bioinfo/data/fastq

```


# Diretório de Resultados

```bash
# volta para workspace
cd /workspace/t12020

# cd para entrar no diretorio resultados
cd bioinfo/resultados/

# mkdir para criar o diretorio das amostras 003, 017 e 019
mkdir 003 017 019

# volta para workspace
cd /workspace/t12020

```


# Download dos Arquivos .FASTQ (no Google Drive) direto do terminal

```bash
# volta para workspace
cd /workspace/t12020

# download do script para copiar os dados 
git clone https://github.com/circulosmeos/gdown.pl.git

# rodar o script ./gdown.pl download
./gdown.pl/gdown.pl https://drive.google.com/open?id=1LlBQ2BNI_vF-E_4ntrEexJQ2-Wea2e0o 003.fastq.gz
./gdown.pl/gdown.pl https://drive.google.com/open?id=1KKmk90gUk0174MvATZzWRwyqRSHJNZDv 017.fastq.gz
./gdown.pl/gdown.pl https://drive.google.com/open?id=11jtliN2G0vTs50-z79QxXLWkStPeuh2E 019.fastq.gz

# movendo arquivos fastq para outro diretorio
mv 003.fastq.gz 017.fastq.gz 019.fastq.gz  bioinfo/data/fastq

# voltar para o home
cd /workspace/t12020
```


# Instalação

No terminal do Linux, vamos instalar alguns pacotes: (o resto vem instalado).

``` bash
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

# voltar para workspace
cd /workspace/t12020
```

# Download das Referências
Nesta estapa vamos utilizar dois cromossomos humanos (chr13 e chr17), nossos genes de interesse são: BRCA1 e BRCA2. [NCBI BRCA 1 and 2](https://www.ncbi.nlm.nih.gov/books/NBK470239/)

Acessar o site: [Sequence and Annotation Downloads](http://hgdownload.cse.ucsc.edu/downloads.html)


```bash
# voltar para workspace
cd /workspace/t12020

cd bioinfo/reference

# wget para fazer download do chr13
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz

# wget para fazer download do chr17 
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz

# cat para concatenar os arquivo do cromossomo 13 e 17 em hg19.fa
zcat chr13.fa.gz chr17.fa.gz > hg19.fa

# rm para deletar os arquivos chr13.fa e chr17.fa
rm chr13.fa.gz chr17.fa.gz

# voltar para workspace
cd /workspace/t12020
```

# BWA: index reference
Agora, é preciso indexar o arquivo FASTA da referência. Todos os programas de alinhmaneto tem esta etapa que otimiza o processo de alinhamento das nossas sequências na referência.

***NOTA:*** A etapa de index é feita apenas uma vez para cada arquivo de referência:

```bash
# voltar para workspace
cd /workspace/t12020

cd bioinfo/reference

# bwa index para gerar o index da referencia hg19.fa

# entrar no diretorio reference
cd bioinfo/reference

# indexar o arquivo hg19.fa
bwa index hg19.fa 

# voltar para workspace
cd /workspace/t12020
```

# FastQC: Relatório de Controle de Qualidade
Gerar relatório de controle de qualidade com FastQC (Tempo ~10s):

```bash
# fastqc para gerar relatorio de qualidade dos arquivo FASTQ
# opcao (-o ./) diz para salvar o resultado no mesmo arquivo em que o comando esta sendo rodado.

# voltar para workspace
cd /workspace/t12020

# rodar fastqc e salvar o resultado de cada amostra em seu diretorio
fastqc -o bioinfo/resultados/003/ bioinfo/data/fastq/003.fastq.gz 
```

FastQC 003 resultado

### Tarefa 01: Repetir o processo para as amostras 017 e 019

~10 min


#cutadapt: localizar e remover adaptadres
O Cutadapt localiza e remove sequências de adaptadores, primers, caudas poly-A e outros tipos de sequência indesejada. Aqui vamos utilizar as funções para "trimar" sequências pequenas e maiores do que o esperado. (Tempo ~3s):

> Gitpod Note: Execute esse comando no terminal do gitpod caso o cutadapt não esteja instalado.

```bash
python3 -m pip install --user --upgrade cutadapt
```

```bash

 -m LEN[:LEN2], --minimum-length=LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0
 -M LEN[:LEN2], --maximum-length=LEN[:LEN2]
                        Discard reads longer than LEN. Default: no limit
 
# voltar para workspace
cd /workspace/t12020

# cutadapt para remover sequencias de tamanho menores 100pb e maiores que 220
cutadapt --minimum-length 100 --maximum-length 220 -q 15  -o bioinfo/resultados/003/003.cutadapt.fastq  bioinfo/data/fastq/003.fastq.gz
```

003.cutadapt.fastq resultado

### Tarefa 02: Repetir o processo para as amostras 017 e 019

# BWA-mem: maximal exact matches
Alinha sequencias de tamanho 70bp-1Mbp com o algoritmo BWA-MEM. Em resumo o algoritmo trabalha com "alinhamento por sementes" com maximal exact matches (MEMs) e então estendendo sementes com o algoritmo Smith-Waterman (SW). Link. Tempo (~60s):

```bash
# voltar para workspace
cd /workspace/t12020

# rodar bwa para alinhar as sequencias contra o genoma de referencia
bwa mem -R '@RG\tID:003\tSM:003_NGSA\tLB:XPTO\tPL:Ion'  bioinfo/reference/hg19.fa bioinfo/resultados/003/003.cutadapt.fastq > bioinfo/resultados/003/003.sam
```

BWA-mem 003 resultado

### Tarefa 03: Repetir o processo para as amostras 017 e 019

# samtools: fixmate, sort e index

## samtools fixmate

Preencha coordenadas de posicionamento, posicione FLAGs relacionadas a partir a alinhamentos classificados por nome. Tempo (~5s):

```bash
# voltar para workspace
cd /workspace/t12020

samtools fixmate bioinfo/resultados/003/003.sam bioinfo/resultados/003/003.bam
```

SAMTOOLS fixmate 003 resultado

### Tarefa 04: Repetir o processo para as amostras 017 e 019

Tempo: ~10min

## samtools sort

O ***samtools sort*** vai ordenar de nome para ordem de coordenadas. Tempo (5s):

```bash
# voltar para workspace
cd /workspace/t12020

time samtools sort -O bam -o bioinfo/resultados/003/003_sort.bam -T /tmp/ bioinfo/resultados/003/003.bam 
```

SAMTOOLS sort 003 resultado

### Tarefa 05: Repetir o processo para as amostras 017 e 019

Tempo: ~1min

## samtools index

O ***samtools index*** cria um index (.BAI) do arquivo binário (.BAM):

```bash
# voltar para workspace
cd /workspace/t12020

samtools index bioinfo/resultados/003/003_sort.bam 
```

SAMTOOLS index 003 resultado

### Tarefa 06: Repetir o processo para as amostras 017 e 019

Tempo: ~1min

# freebayes: chamador de variantes
O FreeBayes é um detector variante genético Bayesiano projetado para encontrar pequenos polimorfismos, especificamente SNPs (polimorfismos de nucleotídeo único), indels (inserções e deleções), MNPs (polimorfismos de múltiplos nucleotídeos) e eventos complexos (eventos compostos de inserção e substituição) menores que os comprimento de um alinhamento de seqüenciamento de leitura curta. Link. Tempo (~6min):

```bash
# voltar para workspace
cd /workspace/t12020

freebayes -f bioinfo/reference/hg19.fa -F 0.01 -C 1 --pooled-continuous bioinfo/resultados/003/003_sort.bam > bioinfo/resultados/003/003.vcf
```

### Parâmetros


```bash
 -C --min-alternate-count N
                   Require at least this count of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position.  default: 2

 --pooled-continuous
                   Output all alleles which pass input filters, regardles of
                   genotyping outcome or model.
```


FREEBAYES call variant 003 resultado

### Tarefa 07: Repetir o processo para as amostras 017 e 019

Tempo: ~10min

# annovar: anotador de variantes
ANNOVAR éma ferramenta eficiente para anotar funcionalmente variantes genéticas detectadas a partir de diversos genomas (incluindo o genoma humano hg18, hg19, hg38, bem como mouse, verme, mosca, levedura e muitos outros). [Link]. Aqui vamos converter .VCF para .avinput. Tempo (~5s):


```bash 
# voltar para workspace
cd /workspace/t12020

perl bioinfo/app/annovar/convert2annovar.pl -format vcf4 bioinfo/resultados/003/003.vcf > bioinfo/resultados/003/003.avinput
``` 

***output: mensagens na tela***
 
```bash
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=refGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile resultados/003/003.refGene -exonsort resultados/003/003.avinput /bioinfo/app/annovar/humandb/>
NOTICE: Output files were written to resultados/003/003.refGene.variant_function, resultados/003/003.refGene.exonic_variant_function
Error: cannot read from --queryfile (resultados/003/003.avinput): No such file or directory
Error running system command: <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile resultados/003/003.refGene -exonsort resultados/003/003.avinput /bioinfo/app/annovar/humandb/>
```

**Arquivo .avinput**


```bash
# voltar para workspace
cd /workspace/t12020

# comando head para listar as 10 primeiras linhas do arquivo 003.avinput
head resultados/003/003.avinput

chr13	19127992	19127992	G	A	hom	11.313	1
chr13	19128097	19128097	A	G	hom	19.0318	1
chr13	19128113	19128113	T	A	hom	15.1147	1
chr13	19650636	19650636	T	CAA	hom	12.2464	1
chr13	19650650	19650650	T	A	hom	2.53014	1
chr13	19650671	19650671	G	C	hom	19.0318	1
chr13	19650676	19650683	GCCTGAGC	ACACGACT	hom	11.313	1
chr13	19650695	19650703	TGTATGGAT	CATACAGAG	hom	14.1494	1
chr13	19650738	19650746	TCCTTCACG	CCCTGGACA	hom	16.0868	1
chr13	19650770	19650770	G	A	hom	19.0318	1
```
 
ANNOVAR convert2annovar 003 resultado

### Tarefa 08: Repetir o processo para as amostras 017 e 019

Anotar as variantes chamadas utilizando algumas bases de dados públicas: Tempo (~5s).

```
# voltar para workspace
cd /workspace/t12020

perl bioinfo/app/annovar/table_annovar.pl bioinfo/resultados/003/003.avinput bioinfo/app/annovar/humandb/ -buildver hg19 -out bioinfo/resultados/003/003 -remove -protocol refGene,exac03,clinvar_20200316 -operation g,f,f -nastring . 
```

***output:***

```bash
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=refGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile resultados/003/003.refGene -exonsort resultados/003/003.avinput /bioinfo/app/annovar/humandb/>
NOTICE: Output files were written to resultados/003/003.refGene.variant_function, resultados/003/003.refGene.exonic_variant_function
NOTICE: Reading gene annotation from /bioinfo/app/annovar/humandb/hg19_refGene.txt ... Done with 63481 transcripts (including 15216 without coding sequence annotation) for 27720 unique genes
NOTICE: Processing next batch with 2312 unique variants in 2312 input lines
NOTICE: Reading FASTA sequences from /bioinfo/app/annovar/humandb/hg19_refGeneMrna.fa ... Done with 11 sequences
WARNING: A total of 402 sequences will be ignored due to lack of correct ORF annotation
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=exac03
NOTICE: Finished reading 8 column headers for '-dbtype exac03'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype exac03 -buildver hg19 -outfile resultados/003/003 resultados/003/003.avinput /bioinfo/app/annovar/humandb/ -otherinfo>
NOTICE: the --dbtype exac03 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to resultados/003/003.hg19_exac03_dropped, other variants are written to resultados/003/003.hg19_exac03_filtered
NOTICE: Processing next batch with 2312 unique variants in 2312 input lines
NOTICE: Database index loaded. Total number of bins is 749886 and the number of bins to be scanned is 70
NOTICE: Scanning filter database /bioinfo/app/annovar/humandb/hg19_exac03.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=clinvar_20190114
NOTICE: Finished reading 5 column headers for '-dbtype clinvar_20190114'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype clinvar_20190114 -buildver hg19 -outfile resultados/003/003 resultados/003/003.avinput /bioinfo/app/annovar/humandb/ -otherinfo>
NOTICE: the --dbtype clinvar_20190114 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to resultados/003/003.hg19_clinvar_20190114_dropped, other variants are written to resultados/003/003.hg19_clinvar_20190114_filtered
NOTICE: Processing next batch with 2312 unique variants in 2312 input lines
NOTICE: Database index loaded. Total number of bins is 45640 and the number of bins to be scanned is 54
NOTICE: Scanning filter database /bioinfo/app/annovar/humandb/hg19_clinvar_20190114.txt...Done
-----------------------------------------------------------------
NOTICE: Multianno output file is written to resultados/003/003.hg19_multianno.txt
```

**Arquivo .hg19_multiano.txt**

```bash

# voltar para workspace
cd /workspace/t12020

# filtro com o comando grep para buscar apenas o cabeçalho e variantes do tipo exonic
grep "^Chr\|exonic" resultados/003/003.hg19_multianno.txt  | head

Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	ExAC_ALL	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	CLNALLELEID	CLNDN	CLNDISDB	CLNREVSTAT	CLNSIG
chr13	32906565	32906565	-	A	exonic	BRCA2	.	frameshift insertion	BRCA2:NM_000059:exon10:c.951dupA:p.T317fs	8.321e-06	0	0	0.0001	0	0	0	0	67538	Hereditary_breast_and_ovarian_cancer_syndrome|Familial_cancer_of_breast|Hereditary_cancer-predisposing_syndrome|Breast-ovarian_cancer,_familial_2|not_provided	MeSH:D061325,MedGen:C0677776,Orphanet:ORPHA145|MedGen:C0006142,OMIM:114480,Orphanet:ORPHA227535,SNOMED_CT:254843006|MedGen:C0027672,SNOMED_CT:699346009|MedGen:C2675520,OMIM:612555|MedGen:CN517202	reviewed_by_expert_panel	Pathogenic
chr13	32906729	32906729	A	C	exonic	BRCA2	.	nonsynonymous SNV	BRCA2:NM_000059:exon10:c.A1114C:p.N372H	0.2779	0.1249	0.3049	0.2728	0.2331	0.2818	0.2677	0.3558	24368	Hereditary_breast_and_ovarian_cancer_syndrome|Familial_cancer_of_breast|Fanconi_anemia|Hereditary_cancer-predisposing_syndrome|Ductal_breast_carcinoma|Breast-ovarian_cancer,_familial_2|not_specified|not_provided	MeSH:D061325,MedGen:C0677776,Orphanet:ORPHA145|MedGen:C0006142,OMIM:114480,Orphanet:ORPHA227535,SNOMED_CT:254843006|MedGen:C0015625,Orphanet:ORPHA84,SNOMED_CT:30575002|MedGen:C0027672,SNOMED_CT:699346009|MedGen:C1527349|MedGen:C2675520,OMIM:612555|MedGen:CN169374|MedGen:CN517202	reviewed_by_expert_panel	Benign
chr13	32907215	32907215	-	A	exonic	BRCA2	.	frameshift insertion	BRCA2:NM_000059:exon10:c.1601dupA:p.E534fs	.	.
chr13	32907303	32907303	G	-	exonic	BRCA2	.	frameshift deletion	BRCA2:NM_000059:exon10:c.1688delG:p.W563fs	.	234658	Hereditary_cancer-predisposing_syndrome|Breast-ovarian_cancer,_familial_2	MedGen:C0027672,SNOMED_CT:699346009|MedGen:C2675520,OMIM:612555	reviewed_by_expert_panel	Pathogenic
chr13	32907421	32907421	A	-	exonic	BRCA2	.	frameshift deletion	BRCA2:NM_000059:exon10:c.1806delA:p.G602fs	.	46319	Hereditary_breast_and_ovarian_cancer_syndrome|Hereditary_cancer-predisposing_syndrome|Breast-ovarian_cancer,_familial_2|not_provided	MeSH:D061325,MedGen:C0677776,Orphanet:ORPHA145|MedGen:C0027672,SNOMED_CT:699346009|MedGen:C2675520,OMIM:612555|MedGen:CN517202	reviewed_by_expert_panel	Pathogenic
chr13	32910430	32910430	C	T	exonic	BRCA2	.	synonymous SNV	BRCA2:NM_000059:exon11:c.C1938T:p.S646S	0.0009	0.0002	0.0008	0	0	0.0015	0	0	65898	Hereditary_breast_and_ovarian_cancer_syndrome|Familial_cancer_of_breast|Fanconi_anemia|Hereditary_cancer-predisposing_syndrome|Breast-ovarian_cancer,_familial_2|not_specified|not_provided	MeSH:D061325,MedGen:C0677776,Orphanet:ORPHA145|MedGen:C0006142,OMIM:114480,Orphanet:ORPHA227535,SNOMED_CT:254843006|MedGen:C0015625,Orphanet:ORPHA84,SNOMED_CT:30575002|MedGen:C0027672,SNOMED_CT:699346009|MedGen:C2675520,OMIM:612555|MedGen:CN169374|MedGen:CN517202	reviewed_by_expert_panel	Benign
chr13	32910661	32910661	-	A	exonic	BRCA2	.	frameshift insertion	BRCA2:NM_000059:exon11:c.2170dupA:p.S723fs	.	46332	Neoplasm_of_the_breast|Hereditary_breast_and_ovarian_cancer_syndrome|Breast-ovarian_cancer,_familial_2	Human_Phenotype_Ontology:HP:0100013,MeSH:D001943,MedGen:C1458155,Orphanet:ORPHA180250,SNOMED_CT:126926005|MeSH:D061325,MedGen:C0677776,Orphanet:ORPHA145|MedGen:C2675520,OMIM:612555	reviewed_by_expert_panel	Pathogenic
chr13	32911321	32911334	TAAAAAAGATTTGG	AAAAAAAGATTTTGGT	exonic	BRCA2	.	frameshift substitution	BRCA2:NM_000059:exon11:c.2829_2842AAAAAAAGATTTTGGT	.	.	.	.	.	.	.	.	.	.	.	.	.
chr13	32911443	32911443	A	-	exonic	BRCA2	.	frameshift deletion	BRCA2:NM_000059:exon11:c.2951delA:p.E984fs	.	248944	Hereditary_cancer-predisposing_syndrome|Breast-ovarian_cancer,_familial_2	MedGen:C0027672,SNOMED_CT:699346009|MedGen:C2675520,OMIM:612555	reviewed_by_expert_panel	Pathogenic
```

Utilize o commando: 

```bash
less -SN resultados/003/003.hg19_multianno.txt
``` 

para visualizar o arquivo de anotação. Para sair do comando less pressione a tecla `q`.

### Tarefa 09: Repetir o processo para as amostras 017 e 019
