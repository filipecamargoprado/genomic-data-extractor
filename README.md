# 🧬 Genomic Data Extractor

### Automação, Engenharia de Dados e Processamento de Sequências Genômicas via API do NCBI

Este repositório contém a arquitetura de engenharia de dados e os algoritmos desenvolvidos como parte integrante do **Trabalho de Conclusão de Curso (TCC) em Engenharia de Software**. O ecossistema foi projetado para atuar como um pipeline automatizado de extração, tratamento, filtragem e estruturação de dados biológicos complexos, consumindo diretamente as APIs do **NCBI** (*National Center for Biotechnology Information*).

---

## 🎯 Contexto e Objetivo do Projeto

Na bioinformática contemporânea, o volume de dados genômicos cresce em escala exponencial (atingindo a casa dos Petabytes). O acesso manual a essas informações ou o download de arquivos genômicos brutos sem critérios rígidos de subamostragem geram sérios problemas operacionais, tais como:

1. **Estouro de Memória (OOM):** Sequências genômicas completas de mamíferos possuem bilhões de caracteres.
2. **Timeouts de Rede e Bloqueios:** Requisições síncronas pesadas travam as conexões com os servidores globais do NCBI, resultando em rejeição de tráfego.
3. **Massa de Dados Não Estruturada:** O formato padrão de intercâmbio biológico (FASTA) é puramente textual e semiestruturado, exigindo parsers complexos para torná-lo legível para algoritmos de Inteligência Artificial.

**Este projeto resolve esses gargalos.** Através da construção de scripts robustos em Python 3.11+, o sistema executa buscas cirúrgicas por identificadores biomédicos e aplica filtros de curadoria rigorosa em nível de servidor. O resultado final é a conversão de strings biológicas massivas in situ em matrizes tabulares limpas, prontas para alimentar pipelines de Data Science e modelos preditivos de Machine Learning.

---

## ✨ Arquitetura e Funcionalidades do Pipeline

* **Extração Otimizada (Entrez E-utils API):** Integração via protocolo HTTP utilizando o módulo `Bio.Entrez` para realizar consultas indexadas e buscar chaves primárias globais (*Accession Numbers*).
* **Filtro de Curadoria Inteligente (RefSeq):** Injeção de propriedades lógicas diretamente nas queries de busca (`[Gene Name]`, `[Organism]`, `src db refseq`). Isso força o servidor do NCBI a ignorar fragmentos parciais, clones clonados ou sequências redundantes, baixando exclusivamente as **Sequências de Referência Oficiais** revisadas por curadores humanos.
* **Stream Parsing de Arquivos FASTA:** Utilização do `Bio.SeqIO` com buffers de memória dinâmicos (`StringIO`) para realizar o processamento em tempo real do arquivo texto. O parser quebra o arquivo delimitado por cabeçalhos `>` e constrói objetos estruturados em tempo de execução.
* **Cálculo Estatístico de Engenharia Genética:** Algoritmização da métrica de **Conteúdo GC** (Guanina + Citosina em relação ao total de bases). Esta variável estatística determina a estabilidade térmica das fitas de DNA e serve como uma *feature* crucial para a identificação de regiões promotoras (Ilhas CpG) por modelos de IA.
* **Sanitização de Dados:** Sanitização automática das strings de DNA para caixa alta (`.upper()`), eliminação de espaçamentos em branco que corrompem a leitura de vetores e limitação de comprimento de descrições textuais para evitar quebra de células em sistemas de armazenamento tabular.

---

## 🛠️ Tecnologias e Dependências do Ecossistema

* **Runtime Principal:** Python 3.11+
* **Framework de Bioinformática:** `Biopython` (Módulos de rede `Entrez` e de parsing de sequências `SeqIO`)
* **Processamento de Matrizes:** `pandas` (Abstração e manipulação de DataFrames de alta performance)
* **Engine de Persistência Tabular:** `openpyxl` (Motor de escrita e serialização para planilhas em formato OpenXML)

---

## 📦 Estrutura do Repositório Local

```text
genomic-data-extractor/
│
├── baixar_dados_brutos_ncbi/       # Diretório destinado ao armazenamento de arquivos FASTA intermediários
├── agrupar_gz_ncbi.py              # Script utilitário para manipulação e junção de grandes volumes de dados compactados
├── extrair_gz.py                   # Módulo de descompressão de dados em lote vindos de repositórios biológicos
├── tratamento_genes.py             # Script principal (Pipeline de Automação, Processamento e Geração do Excel)
├── historico_arquivos_processados.txt # Arquivo de log de controle para idempotência de arquivos tratados
├── log.txt                         # Relatórios de rastreabilidade de erros e logs de execução do terminal
└── README.md                       # Documentação técnica detalhada do projeto
```

---

## 🚀 Guia de Implantação e Execução

### 1. Clonar o Repositório e Preparar o Ambiente

Abra o terminal do seu sistema operacional (o sistema foi validado em ambiente Windows 10/11 e Linux) e execute os comandos abaixo para clonar o repositório e isolar o ambiente de execução:

```bash
# Clonar o repositório institucional do GitHub
git clone [https://github.com/filipecamargoprado/genomic-data-extractor.git](https://github.com/filipecamargoprado/genomic-data-extractor.git)

# Navegar para o diretório raiz do projeto
cd genomic-data-extractor

# Criar o ambiente virtual isolado (Virtual Environment)
python -m venv venv

# Ativar o ambiente virtual
# No Windows (Prompt de Comando):
venv\Scripts\activate
# No Linux/Mac:
source venv/bin/activate
```

### 2. Instalar as Dependências do Projeto

Instale os pacotes de software com as versões compatíveis de produção através do gerenciador de pacotes do Python:

```bash
pip install pandas biopython openpyxl
```

### 3. Configurar os Parâmetros da API e Executar

Antes de executar, certifique-se de alterar as boas práticas de cabeçalho HTTP no arquivo de execução principal (`tratamento_genes.py`). O NCBI exige o envio de um e-mail válido para identificação de bots e prevenção de Denial of Service (DoS):

```python
# Dentro do script tratamento_genes.py, altere para suas credenciais de uso:
Entrez.email = "seu_email@exemplo.com"
Entrez.tool = "AutomacaoGeral_Bot"
```

Para disparar o pipeline completo de extração do gene alvo (ex: *Homo sapiens BRCA1*), execute o script chamando o interpretador ativo:

```bash
python tratamento_genes.py
```

---

## 📊 Arquitetura de Dados de Saída (Data Schema)

O pipeline exportará as sequências limpas diretamente para o caminho configurado no script (Ex: `C:\_TESTE\dados_geneticos.xlsx`). O arquivo resultante possui um esquema de dados estrito composto por colunas estruturadas, prontas para consumo estatístico imediato:

| Nome da Coluna | Tipo de Dado | Descrição Técnica no Banco de Dados |
| :--- | :--- | :--- |
| **`ID_Banco`** | `VARCHAR / String` | Chave Primária Global única (*Accession Number*) atribuída pelo NCBI ao registro biológico. |
| **`Resumo`** | `TEXT / String` | Metadados limpos contendo a taxonomia completa da espécie, nome do gene e tipo molecular. |
| **`Tamanho_BP`** | `INTEGER / Int` | Comprimento exato da sequência molecular mensurado em **Pares de Bases (Base Pairs)**. |
| **`Conteudo_GC_%`** | `FLOAT / Decimal` | Métrica percentual arredondada que expressa a densidade dos nucleotídeos Guanina e Citosina na amostra. |
| **`Sequencia_Completa`** | `LONG TEXT / String` | A cadeia polimórfica bruta de nucleotídeos (letras A, T, C, G) unificada em formato de string contínua. |

---

## 🔬 Alinhamento com a Pesquisa Científica de TCC

Este pipeline resolve a camada inicial de **Engenharia de Dados**, fornecendo a base de dados tratada para os tópicos mais profundos discutidos no escopo teórico do Trabalho de Conclusão de Curso:

* **Variantes VUS (Variant of Uncertain Significance):** A estruturação das sequências e o cálculo do conteúdo GC em formato tabular são cruciais para cruzar dados do GenBank com bases clínicas mutacionais (como o *ClinVar*), auxiliando médicos a desvendarem se uma mutação de significado incerto é de fato patogênica ou benigna.

* **Mitigação do Desbalanceamento de Classes via IA:** Modelos de inteligência artificial aplicados à genética lidam com bancos altamente desbalanceados (mutações benignas ocorrem em abundância na população, enquanto as causadoras de doenças raras são escassas). A geração deste arquivo `.xlsx` limpo permite a aplicação imediata de técnicas de reamostragem estatística (**SMOTE** ou **Undersampling**) antes do treinamento de algoritmos de classificação como **Random Forest** ou **XGBoost**.

* **Modelagem de Estruturas Tridimensionais (AlphaFold2):** A correta extração das sequências primárias de nucleotídeos (representadas na última coluna da planilha) é o dado bruto de entrada exigido por redes neurais profundas de predição estrutural de proteínas, permitindo antever o dobramento tridimensional das macromoléculas para o desenvolvimento de novos fármacos na **Medicina de Precisão**.

* **Governança de Dados e Conformidade (LGPD):** Embora dados genéticos sejam enquadrados como **Dados Pessoais Sensíveis** pela Lei Geral de Proteção de Dados (LGPD), este software adota uma postura *Privacy by Design*. Ao extrair unicamente informações de bases públicas globais e realizar a separação total entre sequências biológicas de referência e identificadores de prontuários clínicos privados, o ecossistema atua em total conformidade jurídica, garantindo a anonimização sistêmica dos dados.

---

> **Desenvolvido por:** Filipe Camargo Prado  
> **Curso:** Engenharia de Software  
> **Foco de Pesquisa:** Engenharia de Dados, Automação de Processos, Bioinformática e Integração de Sistemas Computacionais com Bases de Dados Biológicas Globais.
