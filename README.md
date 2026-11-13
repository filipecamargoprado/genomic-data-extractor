# 🧬 Genomic Data Extractor & Protein Clustering Pipeline

![Python](https://img.shields.io/badge/Python-3.11%2B-blue?style=for-the-badge&logo=python)
![Bioinformatics](https://img.shields.io/badge/NCBI-GenBank%20%7C%20RefSeq-green?style=for-the-badge)
![MMseqs2](https://img.shields.io/badge/Engine-MMseqs2%20(AVX2)-orange?style=for-the-badge)
![WSL2](https://img.shields.io/badge/Environment-Windows%20%2B%20WSL2%20(Ubuntu)-purple?style=for-the-badge&logo=linux)
![Status](https://img.shields.io/badge/TCC%20I-Conclu%C3%ADdo%20(3.835%20Lotes)-success?style=for-the-badge)

> **Arquitetura de Engenharia de Dados, Automação e Clusterização de Sequências Genômicas e Proteicas em Larga Escala via NCBI, Python e MMseqs2.**

Este repositório contém o código-fonte, os scripts de automação e a engenharia de dados desenvolvidos como parte do **Trabalho de Conclusão de Curso (TCC)** em Engenharia de Software da **Universidade Positivo (Curitiba - PR)**.

---

## 📌 1. Visão Geral e Problema Resolvido

Na bioinformática contemporânea, o volume de dados biológicos depositados em repositórios públicos como o **NCBI (*National Center for Biotechnology Information*)** cresce exponencialmente, atingindo a escala de Terabytes e Petabytes. A extração manual ou o processamento ingênuo dessas bases em máquinas locais esbarra em cinco gargalos críticos de engenharia:

1. **Timeouts de Rede e Limitação de Taxa (API):** Requisições genéricas travam conexões HTTP ao tentar baixar sequências genômicas inteiras não curadas.
2. **Complexidade Algorítmica Quadrática:** Ferramentas clássicas de agrupamento baseadas em *Greedy Incremental Clustering* (como o **CD-HIT**) tornam-se inviáveis diante de dezenas de milhões de sequências.
3. **Gargalo de I/O entre Sistemas Operacionais (Protocolo 9P):** Executar binários nativos de Linux (via WSL2) lendo e gravando diretamente na partição NTFS do Windows (`/mnt/c/`) degrada a performance em até $100\times$.
4. **Vazamento de Memória RAM (`vmmem`):** O acúmulo de cache de página no subsistema Linux e a retenção de dicionários gigantes em Python causam estouro de memória (*Out of Memory* - OOM).
5. **Explosão Volumétrica e Bloqueio em Nuvem (`Errno 28`):** Gravar a cadeia completa de aminoácidos dentro de um único arquivo CSV sincronizado com o Google Drive gera mais de $2,5\text{ TB}$ de texto e bloqueia o descritor de arquivos.

**Este projeto soluciona todos esses gargalos** através de uma pipeline evoluída em duas gerações arquiteturais (**Versão 1.0** e **Versão 2.0**), reduzindo uma carga de trabalho estimada em mais de 60 dias para apenas **72 horas (3 dias)** em um computador pessoal.

---

## 🏗️ 2. Arquitetura do Projeto (Evolução v1.0 $\rightarrow$ v2.0)

### 🔹 Geração 1 (`tratamento_genes.py`) — Extração Curada via API REST & Conteúdo GC
Focada na extração pontual de genes específicos (ex.: *Homo sapiens BRCA1*) diretamente do banco `nucleotide` do NCBI utilizando **Biopython (`Bio.Entrez` e `Bio.SeqIO`)**:
* **Filtro Inteligente RefSeq:** Injeção da diretiva `src db refseq[Properties]` na *query* para ignorar fragmentos genômicos pesados ($\approx 81.000\text{ bp}$) e baixar apenas sequências de referência curadas ($\approx 7.000\text{ bp}$), eliminando *timeouts*.
* **Stream Parsing em Memória:** Uso de `io.StringIO` e gerenciadores de contexto (`with`) para leitura segura sem geração de lixo em disco.
* **Cálculo Termodinâmico (Conteúdo GC):** Filtragem de sequências com comprimento $> 500\text{ bp}$ e cálculo automático do percentual de Guanina e Citosina:
  $$GC\% = 100 \times \frac{\text{contagem}(G) + \text{contagem}(C)}{\text{Tamanho\_BP}}$$
* **Exportação Tabular:** Geração automatizada de planilhas `.xlsx` via `pandas` e `openpyxl`.

### 🔹 Ensaios Intermediários (`CD-HIT` vs. `MMseqs2`)
Durante a transição para o processamento de **3.893 arquivos compactados (`.fsa_aa.gz`)**, testou-se a ferramenta **CD-HIT** (`-M 0 -T 0`). Contudo, sua heurística de contagem de *k-mers* apresentou escalabilidade quase quadrática. A arquitetura migrou definitivamente para o **MMseqs2** (*Many-against-Many sequence searching*), cujo módulo *Linclust* opera em **tempo linear $O(N)$** com vetorização **AVX2**.

### 🔹 Geração 2 (`agrupar_gz_ncbi.py`) — Pipeline Massiva FTP + WSL2 + MMseqs2 + Google Drive
Focada na extração paralelizada e no agrupamento por homologia de milhões de proteínas em ambiente híbrido **Windows + Linux (WSL2/Ubuntu)**:
1. **Coleta Massiva Paralelizada:** Download via FTP do NCBI utilizando `concurrent.futures` com até **15 fluxos simultâneos (*threads*)**, armazenando os 3.893 arquivos `.gz` diretamente em um volume de $5\text{ TB}$ no Google Drive (`G:/Meu Drive/base_genes_tcc/arquivos_gz`).
2. **Isolamento de I/O em Memória Virtual Linux (`/tmp`):** Para eliminar a lentidão do protocolo 9P (`/mnt/c/`), o Python extrai um arquivo `.gz` por vez e instrui o WSL2 a criar o banco (`mmseqs createdb`), agrupar (`mmseqs cluster --min-seq-id 0.5 -c 0.8`) e gerar o arquivo tabular (`mmseqs createtsv`) integralmente dentro da partição nativa `/tmp` do Linux.
3. **Reciclagem Ativa de Memória RAM:** Combinação do arquivo `.wslconfig` (`autoMemoryReclaim=dropcache`) com destruição explícita de objetos em Python (`del dicionario_sequencias` e `gc.collect()`) a cada arquivo processado.
4. **Classificação Computacional Nutricional:** Verificação matemática em conjuntos (`set.issubset`) da presença dos **9 aminoácidos essenciais** (`H, I, L, K, M, F, T, W, V`), classificando cada proteína como `"Completa"` ou `"Incompleta"`.
5. **Desacoplamento da Fita Bruta e Fatiamento em Nuvem:** Remoção da string da sequência completa do CSV final (reduzindo o volume projetado de $2,5\text{ TB}$ para poucos Gigabytes) e gravação em lotes individuais na pasta `csv_resultados` para evitar o erro `OSError: [Errno 28] No space left on device`.
6. **Tolerância a Falhas, Idempotência e Mesclagem (`glob`):** Captura automática de arquivos corrompidos (`gzip.BadGzipFile`, `EOFError` e `zlib.error`), registro de estado em `historico_arquivos_processados.txt` e rotina final de unificação automática com validação matemática de integridade.

---

## 📂 3. Estrutura do Repositório

```text
genomic-data-extractor/
├── baixar_dados_brutos_ncbi/              # Scripts de coleta paralelizada (15 threads) via servidor FTP do NCBI
├── csv_resultados/                        # Diretório de saída para os lotes CSV parciais processados
├── agrupar_gz_ncbi.py                     # Script Principal v2.0 (Orquestrador WSL2, MMseqs2, Classificação e Mesclagem)
├── extrair_gz.py                          # Módulo utilitário de descompressão de arquivos .gz em lote
├── tratamento_genes.py                    # Script Principal v1.0 (Extração via API Entrez, RefSeq, Conteúdo GC e Excel)
├── rodar_automacao.bat                    # Inicializador autônomo em lote para Windows (execução sem IDE)
├── historico_arquivos_processados.txt     # Controle de idempotência (checkpoint em tempo real dos arquivos concluídos)
├── log.txt                                # Auditoria do sistema, rastreamento de arquivos corrompidos e progresso
└── README.md                              # Documentação técnica e acadêmica do projeto
```

---

## 💻 4. Requisitos de Hardware e Ambiente

Para reproduzir o processamento da base completa da Versão 2.0 sem travamentos, a máquina deve atender aos seguintes requisitos:

* **Processador (CPU):** Mínimo de 4 núcleos / 8 *threads* (Recomendado: 8 núcleos / 16 *threads*) com suporte obrigatório ao conjunto de instruções vetoriais **AVX2** e virtualização **Intel VT-x / AMD-V** ativa na BIOS.
* **Memória RAM:** Mínimo de 8 GB (Recomendado: 16 GB a 32 GB) com configuração de limpeza de cache no `.wslconfig`.
* **Armazenamento Local (Disco `C:`):** SSD NVMe M.2 com pelo menos **15 GB a 20 GB livres** para área de manobra (*staging* em `C:/temp_base_genes`).
* **Armazenamento em Nuvem / Externo (Unidade `G:`):** Google Drive Desktop (ou HD Externo) com pelo menos **500 GB a 5 TB disponíveis**.
* **Sistema Operacional:** Windows 10/11 integrado ao **WSL2 (Ubuntu)** e **Python 3.11+**.

---

## 🚀 5. Guia de Instalação e Execução

### Passo 1: Clonar o Repositório e Preparar o Ambiente Python

No terminal do Windows (PowerShell ou Prompt de Comando), execute:

```bash
git clone https://github.com/filipecamargoprado/genomic-data-extractor.git
cd genomic-data-extractor

# Criar e ativar o ambiente virtual (opcional para v1.0)
python -m venv venv
venv\Scripts\activate

# Instalar dependências da Versão 1.0
pip install pandas biopython openpyxl
```

### Passo 2: Configurar o WSL2 (Ubuntu) e o Motor MMseqs2 (Versão 2.0)

1. Certifique-se de que o WSL2 e o Ubuntu estão instalados no Windows:
   ```powershell
   wsl --install -d Ubuntu
   ```
2. Dentro do terminal do **Ubuntu (WSL2)**, baixe e instale o binário vetorizado **AVX2** do MMseqs2:
   ```bash
   wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
   tar xvfz mmseqs-linux-avx2.tar.gz
   sudo mv mmseqs/bin/mmseqs /usr/local/bin/
   rm -rf mmseqs mmseqs-linux-avx2.tar.gz
   ```
3. *(Recomendado)* Crie o arquivo `C:\Users\<SeuUsuario>\.wslconfig` no Windows para evitar retenção de memória RAM pelo processo `vmmem`:
   ```ini
   [wsl2]
   memory=8GB
   autoMemoryReclaim=dropcache
   ```

### Passo 3: Executar a Pipeline

* **Para rodar a Versão 1.0 (API Entrez + Excel):**
  Abra o arquivo `tratamento_genes.py`, configure seu e-mail no cabeçalho `Entrez.email` e execute:
  ```bash
  python tratamento_genes.py
  ```

* **Para rodar a Versão 2.0 (Clusterização Massiva com MMseqs2):**
  Verifique os caminhos das variáveis `download_path` e `g_drive` no topo do arquivo `agrupar_gz_ncbi.py` e dê um **duplo-clique** no arquivo:
  ```cmd
  rodar_automacao.bat
  ```
  *(Ou execute diretamente pelo terminal com `python agrupar_gz_ncbi.py`).*

---

## 📊 6. Esquema de Dados de Saída (Data Schemas)

### Tabela de Saída da Versão 1.0 (`dados_geneticos.xlsx`)

| Nome da Coluna | Tipo de Dado | Descrição Técnica |
| :--- | :--- | :--- |
| `ID_Banco` | `VARCHAR` | Chave Primária única (*Accession Number*) atribuída pelo NCBI ao gene. |
| `Resumo` | `TEXT` | Cabeçalho de metadados contendo espécie, nome do gene e tipo molecular. |
| `Tamanho_BP` | `INTEGER` | Comprimento exato da sequência em Pares de Bases (*Base Pairs*). |
| `Conteudo_GC_%` | `FLOAT` | Percentual de nucleotídeos Guanina (`G`) e Citosina (`C`) na sequência. |
| `Sequencia_Completa` | `LONG TEXT` | Cadeia contínua de nucleotídeos (`A, T, C, G`) higienizada. |

### Tabela Consolidada da Versão 2.0 (`base_proteinas_agrupadas_final.csv`)

Separada pelo delimitador `;` (ponto e vírgula), pronta para ingestão em bancos relacionais ou *frameworks* de Inteligência Artificial:

| Nome da Coluna | Tipo de Dado | Descrição Técnica |
| :--- | :--- | :--- |
| `Arquivo_Origem` | `VARCHAR` | Nome do lote compactado original do NCBI (ex.: `gbvrl18.fsa_aa.gz`). |
| `Grupo_Similaridade` | `VARCHAR` | Identificador do *cluster* baseado na proteína centroide (ex.: `Grupo_CAL6443558.1`). |
| `ID_da_Proteina` | `VARCHAR` | Identificador único (*Accession Number*) da proteína membro avaliada. |
| `Tipo_Nutricional` | `VARCHAR` | Classificação computacional (`Completa` ou `Incompleta`) baseada nos 9 aminoácidos essenciais. |

> 💡 **Recuperação de Sequência Sob Demanda:** A cadeia primária completa de aminoácidos foi desacoplada do CSV final para evitar a geração de $2,5\text{ TB}$ de texto redundante. Para recuperar a sequência bruta de qualquer proteína de interesse, basta cruzar os campos `Arquivo_Origem` e `ID_da_Proteina` diretamente no repositório `.gz`.

---

## 📈 7. Resultados Alcançados e Alinhamento com IA

* **Validação Estatística Inicial:** **11.322 proteínas** processadas e agrupadas em **7.339 *clusters* de similaridade** (redução de $35,18\%$ de redundância com parâmetros $Id \ge 50\%$ e $Cov \ge 80\%$).
* **Escala Completa de Produção:** **3.893 arquivos `.fsa_aa.gz`** varridos na base bruta, resultando em **3.835 lotes íntegros** processados, classificados e mesclados automaticamente (com isolamento cirúrgico de 58 arquivos corrompidos na origem).
* **Preparação para Inteligência Artificial (TCC II):** A base estruturada resolve a camada de Engenharia de Dados para aplicação de algoritmos de Aprendizado de Máquina (*Random Forest*, *XGBoost* e Redes Neurais), mitigação de desbalanceamento de classes (via *SMOTE* / *Undersampling*), análise de Variantes de Significado Incerto (VUS), modelagem estrutural (*AlphaFold2*) e total conformidade com a **LGPD (*Privacy by Design*)** por operar exclusivamente sobre dados públicos anonimizados.

---

## 👥 8. Autores e Orientação

Projeto desenvolvido pelos acadêmicos do curso de **Bacharelado em Engenharia da Computação / Engenharia de Software** da **Universidade Positivo (UP)**:

* **Caio Tsutiya Marcidele**
* **Filipe Camargo Prado**
* **Guilherme Fernando Nascimento**
* **Lohanna Bieti Vieira**
* **Vitor Luiz Mougenot**
* **Yuri Congeca Terato Ramos**

**Orientador:** Prof. Aryel Marlus Repula de Oliveira

---
*Curitiba, Paraná — Brasil (2026)*
