## Análises de ancestralidade global e local para caracterização de zonas híbridas


### Backgroud teórico

**Zona híbrida**:

Zonas híbridas são regiões geográficas onde duas linhagens diferentes apresentam uma zona de contato com reprodução ativa, frequentemente com prole fértil. Essa miscigenação pode ter uma dinâmica muito diversa dependendo de vários aspéctos das populações originais e de demais interações bióticas e abióticas. Da mesma forma, as consequências desse processo também pode variar muito, podendo afetar positivamente e negativamente as populações ancestrais, dependendo da dinâmica e das forças seletivas.

**Mistura genômica**:

Em um cenário onde indivíduos de duas linhagens diferentes (vamos chamar aqui de linhagens ancestrais) se cruzam e geram um indivíduo híbrido F1 (primeira geração). Considerando que esses ancestrais são diploides com cariótipo similar, cada par de cromossomos deste indivíduo híbrido terá um cromossomo de cada linhagem acestral, ou seja, ele vai apresentar 50% (ou 0,5) do meterial genômica de cada ancestral. Na medida que esse híbrido F1 se reproduz com outro F1 (ou outro indivíduo de uma das linhagens ancestrais), essa proporção de 0,5/0,5 vai começar a diluir, gerando indivíduos com uma certa predominância de um ou outro ancestral (por exemplo 0,6/0,4). Essa proporção, inclusive ao longo de próximas gerações de mistura, vai variar respondendo a forças seletivas atuando sobre esses híbridos e sobre a dinâmica da miscigenação entre as duas linhagens. Dessa forma, podemos usar esse tipo de informação para caracterizar aspectos evolutivos sobre a zona híbrida.

**Ancestralidade global x ancestralidade local**:

Com genomas completos, temos algumas possibilidades de caracterizar de forma bem refinada a dinâmica genômica de miscigenação. O que chamaremos aqui de ancestralidade global é a proporção geral, em determinado indivíduo de uma população, do genoma vindo de cada ancestral. Por exemplo, o indivíduo F1 comentado anteriormente terá uma proporção de ancestralidade global de 0,5/0,5. Esse tipo de informação pode ser inferida mesmo com marcadores moleculares clássicos. Mas, com genomas completos, ela se torna mais robusta, mais realista, e sem potenciais viéses de introgressão diferencial. Por outro lado a ancestralidade local só pode ser inferida a partir de genomas completos, visto que ela consegue apresentar, ao longo do genoma de determinado indivíduo, quais segmentos são originais de uma ou outra linhagem ancestral. Com esse tipo de abordagem é possível encontrar regiões adaptativas, datar há quantas gerações a miscigenação está ocorrendo e observar de forma bastante detalhada a dinâmica de mistura genômica.


### Exemplo prático

**Dados usados na demonstração**:

Cromossomo A1, sequenciados à baixa cobertura, de ~60 indivíduos amostrados em um transecto abrangendo a zona híbrida entre *L. guttulus* e *L. geoffroyi*. 

**Pipeline de preparação dos dados (não detalhado)**:

- *Processamento inicial*: Filtragem das reads ([fastp](https://github.com/OpenGene/fastp))> mapeamento contra a referência ([BWA](https://bio-bwa.sourceforge.net/bwa.shtml), [Samtools](http://www.htslib.org/doc/samtools.html)) >  filtros de mapeamento (Samtools).
- *Genotype likelihoods*: mapeamento todos indivíduos > maf e genótipo “estimado” ([ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD))
- *SNP calling*: mapeamento todos indivíduos > haplotype calling ([GATK](https://gatk.broadinstitute.org/hc/en-us)) > filtragem cobertura, qualidade, SNP (vcftools, bcftools)
Um exemplo simples de como fazer esses processos pode ser encontrado [aqui](https://github.com/ffertrindade/EvolGenomics/tree/dev/working/3-2%20Population%20Genomics).

**Análise de ancestralidade global:**

Software:
- NGSadmix [software](http://www.popgen.dk/software/index.php/NgsAdmix) e [artigo](https://academic.oup.com/genetics/article/195/3/693/5935455).

Arquivos iniciais e scripts:

- *leopardus_63ind_chrA1_minInd21.beagle.gz*: arquivo com os genotype likelihoods estimados para cada indivíduo.
- *leopardus_63ind_chrA1_minInd21.mafs.gz*: arquivo com as frequências alélicas da população.
- *plotNGSadmix.R*: script para plotar resultado da análise do NGSadmix.
- *leopardus_63ind_chrA1_minInd21.popinfo*: arquivo tabular com informações das amostras.

Comandos:

- Rodar NGSadmix
```
NGSadmix -likes leopardus_63ind_chrA1_minInd21.beagle.gz -K 2 -P 3 -minMaf 0.05 -o leopardus_63ind_chrA1_minInd21_k2
```
- Plotar NGSadmix no R Studio

Arquivos gerados:

- *leopardus_63ind_chrA1_minInd21_k2.qopt*: proporção por indivíduo considerando o número de linhagens ancestrais modeladas.
- *leopardus_63ind_chrA1_minInd21_k2.fopt.gz*: frequências por sítio calculadas considerando o número de linhagens ancestrais modeladas.
- *leopardus_63ind_chrA1_minInd21_k2.log*: informações da análise.

**Análise de ancestralidade local:**

Software:

- Ancestry_HMM [software](https://github.com/russcd/Ancestry_HMM) e [artigo](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006529).

Arquivos iniciais e scripts:

- *leopardus_63ind.chrA1.filtered.vcf*: arquivo VCF filtrado com todos os indivíduos do transecto (híbridos e não híbridos).
- *popNames.txt*: arquivo de texto com indicação de quais indivíduos pertencem a cada população (ex. parental 1, parental 2 e híbridos).
- [*vcf2aHMM.py*](https://github.com/russcd/Ancestry_HMM/tree/master/scripts): script para criar arquivo de entrada para rodar o program Ancestry_HMM.
- *samplesPloidy.txt*: lista dos indivíduos a serem avaliados quanto à ancestralidade local (no caso, os híbridos, mas é possível testar com alguns indivíduos das linhagens ancestrais) seguido da ploidia; será usado como arquivo de entrada do Ancestry_HMM.
- *chr.txt*: lista dos cromossomos sendo usados.
- *getAvrLAI.py*: script para preparar arquivos para plotar resultado de posterior e calcular as proporções médias de ancestralidade local.
- *plotLocalAncestry.R*: script para plotar ancestralidade local do genoma.

Comandos:

- Criando painel:
```
python3 vcf2aHMM.py --vcf leopardus_63ind.chrA1.filtered.vcf --pop popNames.txt --rate 1.9 --minGT 0.2 > leopardus_63ind.popNames_parentals.minGT02.panel
```
```
usage: vcf2aHMM.py [-h] [--vcf VCF] [--map MAP] [--pop POP] [--rate RATE] [--dist DIST] [--minDif MINDIF] [--minGT MINGT] [--minDP MINDP] [-geno]

Read in files and paremeters

optional arguments:
  -h, --help       show this help message and exit
  --vcf VCF        name of the vcf
  --map MAP        name of the recombination map file (units in cM/Mb)
  --pop POP        name of the population identity file
  --rate RATE      mean recombination rate estimate in cM/Mb
  --dist DIST      minimum distance between sites in bp
  --minDif MINDIF  minimum frequency difference in parental genotypes
  --minGT MINGT    minimum mean rate of genotype calls each parent population (and admixed with -geno flag
  --minDP MINDP    minimum mean depth for admixed population (if no -geno flag)
  -geno            use admixed genotypes instead of read counts
```

- Rodando Ancestry_HMM:
```
ancestry_hmm -i leopardus_63ind.popNames_parentals.minGT02.panel -s samplesPloidy.txt -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 -1 0.5 --ne 10000 -e 1e-3 1> leopardus_63ind.popNames_parentals.minGT02.out 2> leopardus_63ind.popNames_parentals.minGT02.err
```
```
ancestry_hmm usage:

        required:
                -i [string]             input file name
                -s [string]             sample id and ploidy file
                -a [int] [float] [float] ...
                        number of ancestral populations and ancestry proportion attributable to each
                -p [int] [int] [float]
                        ancestry pulse with format, ancestral population, time,
                        and proportion of final ancestry from this pulse
                        negative time or proportions indicate that parameters are to be estimated

        optional:
                --help                  print this help statement
                --ne [int]              effective population size of the admixed population
                -g                      samples are specified with genotypes rather than read counts
                --precision [int]       modify float and double precision to int
                -v                      viterbi decoding
                -b [int] [int]          number of bootstraps and bootstrap block size in number of SNPs
                --tmax [int]            maximum time of an admixture pulse
                --tmin [int]            minimum time of an admixture pulse
                --tolerance [float]     distance in lnL units to just convergence
                -e [float]              error rates
                -E                      site specific error rates are included
                --fix                   ancestral allele frequencies are certain

        optional and relevant only for multiple pulse models:
                --output-ancestry       output ancestry posteriors rather than pulses
                -r [int]                number of random restarts during nelder-mead optimization
                --pmax [int]            maximum proportion ancestry in an admixture pulse
                --pmin [int]            minimum proportion ancestry in an admixture pulse
```

- Plotar resultados:
```
python getAvrLAI.py bLge-251 ~/aHMM/genotropics/parentals.minGT02/ chr.file
Rscript plotLocalAncestry.R bLge-251 ~/aHMM/genotropics/parentals.minGT02/
```

Arquivos gerados*:*

- *leopardus_63ind.popNames_parentals.minGT02.panel*: painel gerado pelo script vcf2aHMM.py, usando os arquivos samples.vcf e popNames.txt, queue será usado como arquivo de entrada do Ancestry_HMM.
- *sample***.posterior*: prob. posterior para determinada posição do genoma, do indivíduo sendo analisado, ser de cada linhagem ancestral.
- *leopardus_63ind.popNames_parentals.minGT02.err*: detalhamento da corrida do Ancestry_HMM.
- *leopardus_63ind.popNames_parentals.minGT02.out*: resultado da datação do Ancestry_HMM.
