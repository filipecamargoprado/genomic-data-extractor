# 🧬 Genomic Data Extractor & Protein Clustering Pipeline

![Python](https://img.shields.io/badge/Python-3.11%2B-blue?style=for-the-badge&logo=python)
![Bioinformatics](https://img.shields.io/badge/NCBI-GenBank%20%7C%20RefSeq-green?style=for-the-badge)
![MMseqs2](https://img.shields.io/badge/Engine-MMseqs2%20(AVX2)-orange?style=for-the-badge)
![WSL2](https://img.shields.io/badge/Environment-Windows%20%2B%20WSL2%20(Ubuntu)-purple?style=for-the-badge&logo=linux)
![Status](https://img.shields.io/badge/TCC-Conclu%C3%ADdo%20(3.835%20Lotes)-success?style=for-the-badge)

> **Arquitetura de Engenharia de Dados, Automação e Clusterização de Sequências Genômicas e Proteicas em Larga Escala via NCBI, Python e MMseqs2.**

Este repositório contém o código-fonte, os scripts de automação e a engenharia de dados desenvolvidos como parte do **Trabalho de Conclusão de Curso (TCC)** em **Bacharelado em Engenharia de Software** da **Universidade Positivo (Curitiba - PR)**.

---

## 📌 1. Visão Geral e Problema Resolvido

Na bioinformática contemporânea, o volume de dados biológicos depositados em repositórios públicos como o **NCBI (*National Center for Biotechnology Information*)** cresce exponencialmente, atingindo a escala de Terabytes e Petabytes. A extração manual ou o processamento ingênuo dessas bases em máquinas locais esbarra em cinco gargalos críticos de engenharia:

1. **Timeouts de Rede e Limitação de Taxa (API):** Requisições genéricas travam conexões HTTP ao tentar baixar sequências genômicas inteiras não curadas.
2. **Complexidade Algorítmica Quadrática:** Ferramentas clássicas de agrupamento baseadas em *Greedy Incremental Clustering* (como o **CD-HIT**) tornam-se inviáveis diante de dezenas de milhões de sequências.
3. **Gargalo de I/O entre Sistemas Operacionais (Protocolo 9P):** Executar binários nativos de Linux (via WSL2) lendo e gravando diretamente na partição NTFS do Windows (`/mnt/c/`) degrada a performance em até $100\times$.
4. **Vazamento de Memória RAM (`vmmem`):** O acúmulo de cache de página no subsistema Linux e a retenção de dicionários gigantes em Python causam estouro de memória (*Out of Memory* - OOM).
5. **Explosão Volumétrica e Bloqueio em Nuvem (`Errno 28`):** Gravar a cadeia completa de aminoácidos dentro de um único arquivo CSV sincronizado com o Google Drive gera mais de $2,5\text{ TB}$ de texto e bloqueia o descritor de arquivos.

**Este projeto soluciona todos esses gargalos** através de uma pipeline evoluída de protótipos exploratórios para uma arquitetura escalável de extração e agrupamento em lote, reduzindo uma carga de trabalho estimada em mais de 60 dias para apenas **72 horas (3 dias)** em um computador pessoal.

---

## 🏗️ 2. Evolução Arquitetural da Pipeline

### 🔹 Prototipagem Inicial (Estudo Preliminar) — Extração via API REST & Conteúdo GC
Antes da estruturação deste repositório, realizou-se um estudo exploratório voltado à extração pontual de genes específicos (ex.: *Homo sapiens BRCA1*) diretamente do banco `nucleotide` do NCBI utilizando **Biopython (`Bio.Entrez` e `Bio.SeqIO`)**:
- **Filtro Inteligente RefSeq:** Injeção da diretiva `src db refseq[Properties]` na *query* para ignorar fragmentos genômicos pesados ($\approx 81.000\text{ bp}$) e baixar apenas sequências de referência curadas ($\approx 7.000\text{ bp}$), eliminando *timeouts*.
- **Stream Parsing em Memória:** Uso de `io.StringIO` e gerenciadores de contexto (`with`) para leitura segura sem geração de lixo em disco.
- **Cálculo Termodinâmico (Conteúdo GC):** Filtragem de sequências com comprimento $> 500\text{ bp}$ e cálculo automático do percentual de Guanina e Citosina:

$$GC\% = 100 \times \frac{\text{contagem}(G) + \text{contagem}(C)}{\text{Tamanho\_bp}}$$

### 🔹 Ensaios Comparativos de Agrupamento (`CD-HIT` vs. `MMseqs2`)
Durante a transição para o processamento massivo das divisões do **GenBank** (**3.893 arquivos compactados `.fsa_aa.gz`**), testou-se a ferramenta **CD-HIT** (`-M 0 -T 0`). Contudo, sua heurística de contagem de *k-mers* apresentou escalabilidade quase quadrática. A arquitetura migrou definitivamente para o **MMseqs2** (*Many-against-Many sequence searching*), executado com vetorização **AVX2** em ambiente Linux.

### 🔹 Arquitetura Versionada (Entrega 1 `97bfb03` $\rightarrow$ Versão Atual) — Pipeline FTP + WSL2 + MMseqs2 + Google Drive
A partir do commit inicial (`97bfb03`), consolidou-se a pipeline massiva em ambiente híbrido **Windows + Linux (WSL2/Ubuntu)**, posteriormente refatorada para máxima eficiência de memória e disco:
1. **Coleta Massiva Paralelizada:** Download via servidor FTP do NCBI (`/ncbi-asn1/protein_fasta/`) utilizando `concurrent.futures` com até **15 fluxos simultâneos (*threads*)**, armazenando os 3.893 arquivos `.gz` do GenBank diretamente em um volume de $5\text{ TB}$ no Google Drive (`G:/Meu Drive/base_genes_tcc/arquivos_gz`).
2. **Isolamento de I/O no Disco Virtual Linux (`/tmp`):** Para eliminar a lentidão do protocolo 9P (`/mnt/c/`), o Python extrai um arquivo `.gz` por vez e instrui o WSL2 a criar o banco (`mmseqs createdb`), executar o agrupamento (`mmseqs cluster --min-seq-id 0.5 -c 0.8`) e gerar o arquivo tabular (`mmseqs createtsv`) integralmente dentro da partição nativa `ext4` (`/tmp`) do Linux.
3. **Reciclagem Ativa de Memória RAM:** Combinação do arquivo `.wslconfig` (`autoMemoryReclaim=dropcache`) com destruição explícita de objetos em Python (`del dicionario_sequencias` e `gc.collect()`) a cada arquivo processado.
4. **Classificação Computacional Nutricional:** Verificação matemática em conjuntos (`set.issubset`) da presença dos **9 aminoácidos essenciais** (`H, I, L, K, M, F, T, W, V`), classificando cada proteína como `"Completa"` ou `"Incompleta"`.
5. **Desacoplamento da Cadeia Polipeptídica e Fatiamento em Nuvem:** Remoção da string da sequência completa do CSV final (reduzindo o volume projetado de $2,5\text{ TB}$ para poucos Gigabytes) e gravação em lotes CSV parciais no Google Drive para evitar o erro `OSError: [Errno 28] No space left on device`.
6. **Tolerância a Falhas, Idempotência e Mesclagem (`glob`):** Captura automática de arquivos corrompidos (`gzip.BadGzipFile`, `EOFError` e `zlib.error`), registro de estado em `historico_arquivos_processados.txt` e rotina final de unificação automática com validação matemática de integridade.

---

## 📂 3. Estrutura do Repositório

    genomic-data-extractor/
    ├── baixar_dados_brutos_ncbi/
    │   ├── baixar_gz.py                        # Coleta paralelizada (15 threads) via servidor FTP do NCBI
    │   └── extrair_gz.py                       # Módulo utilitário de descompressão de arquivos .gz em lote
    ├── agrupar_gz_ncbi.py                      # Script Principal (Orquestrador WSL2, MMseqs2, Classificação e Mesclagem)
    ├── rodar_automacao.bat                     # Inicializador autônomo em lote para Windows (execução sem IDE)
    ├── historico_arquivos_processados.txt      # Controle de idempotência (checkpoint em tempo real dos arquivos concluídos)
    ├── log.txt                                 # Auditoria do sistema, rastreamento de arquivos corrompidos e progresso
    └── README.md                               # Documentação técnica e acadêmica do projeto

*(Nota: Os diretórios de saída dos arquivos compactados e dos lotes CSV parciais são gerados dinamicamente na unidade de armazenamento em nuvem configurada, como `G:/Meu Drive/base_genes_tcc/`, mantendo o repositório leve e focado no código-fonte e nos logs de rastreabilidade).*

---

## 💻 4. Requisitos de Hardware e Ambiente

Para reproduzir o processamento da base completa sem gargalos, a máquina deve atender aos seguintes requisitos:

- **Processador (CPU):** Mínimo de 4 núcleos / 8 *threads* (Recomendado: 8 núcleos / 16 *threads*) com suporte obrigatório ao conjunto de instruções vetoriais **AVX2** e virtualização **Intel VT-x / AMD-V** ativa na BIOS.
- **Memória RAM:** Mínimo de 8 GB a 16 GB (Recomendado: 32 GB) com configuração de limpeza de cache no `.wslconfig`.
- **Armazenamento Local (Disco `C:`):** SSD NVMe M.2 com pelo menos **15 GB a 20 GB livres** para área de manobra (*staging* em `C:/temp_base_genes`).
- **Armazenamento em Nuvem / Externo (Unidade `G:`):** Google Drive Desktop (ou HD Externo) com pelo menos **500 GB a 5 TB disponíveis**.
- **Sistema Operacional:** Windows 10/11 integrado ao **WSL2 (Ubuntu)** e **Python 3.11+**.

---

## 🚀 5. Guia de Instalação e Execução

### Passo 1: Clonar o Repositório

No terminal do Windows (PowerShell ou Prompt de Comando), execute:

    git clone https://github.com/filipecamargoprado/genomic-data-extractor.git
    cd genomic-data-extractor

### Passo 2: Configurar o WSL2 (Ubuntu) e o Motor MMseqs2

1. Certifique-se de que o WSL2 e o Ubuntu estão instalados no Windows:

       wsl --install -d Ubuntu

2. Dentro do terminal do **Ubuntu (WSL2)**, baixe e instale o binário vetorizado **AVX2** do MMseqs2:

       wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
       tar xvfz mmseqs-linux-avx2.tar.gz
       sudo mv mmseqs/bin/mmseqs /usr/local/bin/
       rm -rf mmseqs mmseqs-linux-avx2.tar.gz

3. *(Recomendado)* Crie o arquivo `C:\Users\<SeuUsuario>\.wslconfig` no Windows para evitar retenção de memória RAM pelo processo `vmmem`:

       [wsl2]
       memory=8GB
       autoMemoryReclaim=dropcache

### Passo 3: Executar a Pipeline

Verifique os caminhos das variáveis `download_path` e `g_drive` no topo do arquivo `agrupar_gz_ncbi.py` e dê um **duplo-clique** no arquivo:

    rodar_automacao.bat

*(Ou execute diretamente pelo terminal com `python agrupar_gz_ncbi.py`).*

---

## 📊 6. Esquema de Dados de Saída (Data Schema)

O arquivo consolidado (`base_proteinas_agrupadas_final.csv`) é gerado com o delimitador `;` (ponto e vírgula), pronto para ingestão em bancos relacionais ou *frameworks* de Inteligência Artificial:

| Nome da Coluna | Tipo de Dado | Descrição Técnica |
| :--- | :--- | :--- |
| `Arquivo_Origem` | `VARCHAR` | Nome do lote compactado original do NCBI (ex.: `gbvrl18.fsa_aa.gz`). |
| `Grupo_Similaridade` | `VARCHAR` | Identificador do *cluster* baseado na proteína centroide (ex.: `Grupo_CAL6443558.1`). |
| `ID_da_Proteina` | `VARCHAR` | Identificador único (*Accession Number*) da proteína membro avaliada. |
| `Tipo_Nutricional` | `VARCHAR` | Classificação computacional (`Completa` ou `Incompleta`) baseada nos 9 aminoácidos essenciais. |

> 💡 **Recuperação da Cadeia Polipeptídica Sob Demanda:** A sequência primária completa de aminoácidos foi desacoplada do CSV final para evitar a geração de $2,5\text{ TB}$ de texto redundante. Para recuperar a sequência bruta de qualquer proteína de interesse, basta cruzar os campos `Arquivo_Origem` e `ID_da_Proteina` diretamente no arquivo `.gz` correspondente.

---

## 📈 7. Resultados Alcançados e Integração com IA

- **Validação Estatística Inicial:** **11.322 proteínas** processadas e agrupadas em **7.339 *clusters* de similaridade** (redução de $35,18\%$ de redundância com parâmetros $Id \ge 50\%$ e $Cov \ge 80\%$).
- **Escala Completa de Produção (~72 horas):** **3.893 arquivos `.fsa_aa.gz`** varridos na base bruta do GenBank, resultando em **3.835 lotes íntegros** processados, classificados e mesclados automaticamente (com isolamento de 58 arquivos corrompidos na origem).
- **Integração com Inteligência Artificial (Etapa Seguinte):** Dos 3.835 lotes consolidados, **3.428 pares `.gz`/CSV** compuseram a base histórica elegível utilizada no módulo de Aprendizado de Máquina do projeto (repositório [`tcc-proteinas-ia`](https://github.com/KatsuMouley/tcc-proteinas-ia)), servindo como infraestrutura para geração de *embeddings* proteicos (ESM-2), particionamento livre de vazamento por homologia e treinamento de redes neurais (MLP) em total conformidade com a **LGPD (*Privacy by Design*)**.

---

## 👥 8. Autores e Orientação

Projeto desenvolvido pelos acadêmicos do curso de **Bacharelado em Engenharia de Software** da **Universidade Positivo (UP)**:

- **Caio Tsutiya Marcidele**
- **Filipe Camargo Prado**
- **Guilherme Fernando Nascimento**
- **Lohanna Bieti Vieira**
- **Vitor Luiz Mougenot**
- **Yuri Congeca Terato Ramos**

**Orientador:** Prof. Aryel Marlus Repula de Oliveira

---
*Curitiba, Paraná — Brasil (2026)*
