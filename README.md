# 🧬 Genomic Data Extractor & Protein Clustering Pipeline

![Python](https://img.shields.io/badge/Python-3.11%2B-blue?style=for-the-badge&logo=python)
![Bioinformatics](https://img.shields.io/badge/NCBI-GenBank%20%7C%20RefSeq-green?style=for-the-badge)
![MMseqs2](https://img.shields.io/badge/Engine-MMseqs2%20(AVX2)-orange?style=for-the-badge)
![WSL2](https://img.shields.io/badge/Environment-Windows%20%2B%20WSL2%20(Ubuntu)-purple?style=for-the-badge&logo=linux)
![Status](https://img.shields.io/badge/TCC-Conclu%C3%ADdo%20(3.835%20Lotes)-success?style=for-the-badge)

> **Arquitetura de Engenharia de Dados, Automação e Clusterização de Sequências Proteicas em Larga Escala via NCBI, Python e MMseqs2.**

Este repositório contém o código-fonte, os scripts de automação e a engenharia de dados desenvolvidos como parte do **Trabalho de Conclusão de Curso (TCC)** em **Bacharelado em Engenharia de Software** da **Universidade Positivo (Curitiba - PR)**.

---

## 📌 1. Visão Geral e Problema Resolvido

Na bioinformática contemporânea, o volume de dados biológicos depositados em repositórios públicos como o **NCBI (*National Center for Biotechnology Information*)** cresce exponencialmente, atingindo a escala de Terabytes e Petabytes. A extração manual ou o processamento ingênuo dessas bases em máquinas locais esbarra em cinco gargalos críticos de engenharia:

1. **Timeouts de Rede e Limitação de Taxa (API):** Requisições síncronas genéricas travam conexões HTTP ao tentar baixar sequências genômicas inteiras não curadas.
2. **Complexidade Algorítmica em Larga Escala:** Ferramentas clássicas de agrupamento baseadas em *Greedy Incremental Clustering* (como o **CD-HIT**) tornam-se lentas diante de dezenas de milhões de sequências.
3. **Gargalo de I/O entre Sistemas Operacionais (Protocolo 9P):** Executar binários nativos de Linux (via WSL2) lendo e gravando bancos temporários diretamente na partição NTFS do Windows (`/mnt/c/`) degrada o desempenho de entrada e saída em até $100\times$.
4. **Vazamento de Memória RAM (`vmmem`):** O acúmulo de cache de página no subsistema Linux e a retenção de dicionários extensos em Python causam estouro de memória (*Out of Memory* - OOM).
5. **Explosão Volumétrica e Bloqueio em Nuvem (`Errno 28`):** Gravar a cadeia polipeptídica completa dentro de um único arquivo CSV sincronizado com o Google Drive gera mais de $2,5\text{ TB}$ de texto e bloqueia o descritor de arquivos durante o *upload*.

**Este projeto soluciona todos esses gargalos** através de uma pipeline otimizada que reduziu uma carga de trabalho estimada em mais de 60 dias para apenas **72 horas (3 dias)** em um computador pessoal.

---

## 🏗️ 2. Evolução Arquitetural da Pipeline

### 🔹 Prototipagem Exploratória (API REST & Conteúdo GC)
Antes da consolidação deste repositório, foram realizados estudos preliminares de extração pontual de genes específicos (ex.: *Homo sapiens BRCA1*) diretamente do banco `nucleotide` do NCBI utilizando **Biopython (`Bio.Entrez` e `Bio.SeqIO`)** com filtro de curadoria `src db refseq[Properties]` e cálculo termodinâmico de **Conteúdo GC** ($100 \times \frac{G + C}{\text{Tamanho\_bp}}$), além de ensaios laboratoriais comparando **CD-HIT** e **MMseqs2** (*Many-against-Many sequence searching*).

### 🔹 Entrega 1 (Commit `97bfb03`) — Estruturação Base em Lote
A primeira entrega versionada neste repositório estruturou o fluxo automatizado em três módulos (`baixar_gz.py`, `extrair_gz.py` e `agrupar_gz_ncbi.py`) para coleta via FTP do NCBI e agrupamento inicial por similaridade.

### 🔹 Versão 2.0 (Atual) — Pipeline Massiva FTP + WSL2 + MMseqs2 + Google Drive
Focada na extração paralelizada e no agrupamento por homologia de milhões de proteínas do GenBank em ambiente híbrido **Windows + Linux (WSL2/Ubuntu)**:

1. **Coleta Massiva Paralelizada:** Download via servidor FTP do NCBI (`/ncbi-asn1/protein_fasta/`) utilizando `concurrent.futures` com até **15 fluxos simultâneos (*threads*)**, armazenando os **3.893 arquivos `.fsa_aa.gz`** diretamente em um volume de $5\text{ TB}$ no Google Drive (`G:/Meu Drive/base_genes_tcc/arquivos_gz`).
2. **Isolamento de I/O no Sistema de Arquivos Nativo Linux (`/tmp`):** Para eliminar a lentidão da ponte 9P (`/mnt/c/`), o Python extrai um arquivo `.gz` por vez e instrui o WSL2 a criar o banco (`mmseqs createdb`), executar o agrupamento em cascata (`mmseqs cluster --min-seq-id 0.5 -c 0.8`) e converter o resultado (`mmseqs createtsv`) integralmente dentro do disco virtual nativo (`ext4`) do Linux (`/tmp`).
3. **Reciclagem Ativa de Memória RAM:** Combinação do arquivo `.wslconfig` (`autoMemoryReclaim=dropcache`) com destruição explícita de objetos em Python (`del dicionario_sequencias` e `gc.collect()`) a cada arquivo processado.
4. **Classificação Computacional de Aminoácidos Essenciais:** Verificação matemática em conjuntos (`set.issubset`) da presença dos **9 aminoácidos essenciais** (`H, I, L, K, M, F, T, W, V`), atribuindo a *flag* computacional `"Completa"` ou `"Incompleta"`.
5. **Desacoplamento da Cadeia Polipeptídica e Fatiamento em Nuvem:** Remoção da string da sequência completa do CSV final (reduzindo o volume projetado de $2,5\text{ TB}$ para poucos Gigabytes) e gravação em lotes individuais no Google Drive para evitar o erro `OSError: [Errno 28] No space left on device`.
6. **Tolerância a Falhas, Idempotência e Consolidação (`glob`):** Captura automática de arquivos corrompidos (`gzip.BadGzipFile`, `EOFError` e `zlib.error`), registro de *checkpoint* em `historico_arquivos_processados.txt` e rotina final de mesclagem com validação matemática de integridade.

---

## 📂 3. Estrutura do Repositório

```text
genomic-data-extractor/
├── baixar_dados_brutos_ncbi/
│   ├── baixar_gz.py                        # Coleta paralelizada (15 threads) via servidor FTP do NCBI
│   └── extrair_gz.py                       # Módulo utilitário de descompressão de arquivos .gz em lote
├── agrupar_gz_ncbi.py                      # Script Principal v2.0 (Orquestrador WSL2, MMseqs2, Classificação e Mesclagem)
├── historico_arquivos_processados.txt      # Controle de idempotência (checkpoint em tempo real dos arquivos concluídos)
├── log.txt                                 # Auditoria do sistema, rastreamento de arquivos corrompidos e progresso
└── README.md                               # Documentação técnica e acadêmica do projeto
