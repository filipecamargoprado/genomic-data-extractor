import os
import gzip
import shutil
import subprocess
import csv
import logging
from datetime import datetime, timedelta

# ==========================================
# 1. VARIÁVEIS DE CAMINHOS
# ==========================================

download_path = "C:/Users/filip/Downloads/codigos_gerar_genes"
# Arquivo que apenas anota quem já foi feito (para o script pular)
ARQUIVO_LOG = os.path.join(download_path, "historico_arquivos_processados.txt")
# Novo arquivo de log geral do sistema (onde tudo fica registrado)
ARQUIVO_LOG_SISTEMA = os.path.join(download_path, "log.txt")

NOME_AUTOMACAO = "Agrupamento MMseqs2"

# Chamamos o motor do Linux
EXECUTAVEL_LINUX = "wsl"

g_drive = "G:/Meu Drive/base_genes_tcc"
DIR_GZ = os.path.join(g_drive, "arquivos_gz")
ARQUIVO_CSV_FINAL = os.path.join(g_drive, "base_proteinas_agrupadas_final.csv")

# Pasta local temporária para não pesar o Drive
DIR_TEMP_LOCAL = "C:/temp_base_genes"

os.makedirs(DIR_TEMP_LOCAL, exist_ok=True)

# Lista dos 9 aminoácidos essenciais (código de 1 letra)
AMINOACIDOS_ESSENCIAIS = {'H', 'I', 'L', 'K', 'M', 'F', 'T', 'W', 'V'}

# ==========================================
# 2. SISTEMA DE LOGS E CONTROLE
# ==========================================
def setup_logging(log_path, dias_manter=5):
    if os.path.exists(log_path):
        try:
            mod_time = os.path.getmtime(log_path)
            agora = datetime.now()
            idade = agora - datetime.fromtimestamp(mod_time)
            
            if idade > timedelta(days=dias_manter):
                os.remove(log_path)
                print(f"[SISTEMA] Log antigo de {idade.days} dias removido com sucesso.")
        except Exception as e:
            print(f"[AVISO] Não foi possível remover log antigo: {e}")

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, encoding='utf-8', mode='a'),
            logging.StreamHandler()
        ]
    )
    logging.info(f"--- INICIANDO {NOME_AUTOMACAO} ---")

def carregar_arquivos_ja_processados():
    if not os.path.exists(ARQUIVO_LOG):
        return set()
    with open(ARQUIVO_LOG, 'r', encoding='utf-8') as log:
        return set(linha.strip() for linha in log)

def registrar_sucesso(nome_arquivo):
    with open(ARQUIVO_LOG, 'a', encoding='utf-8') as log:
        log.write(nome_arquivo + '\n')

def descompactar_gz(caminho_gz, caminho_fasta):
    nome_arquivo = os.path.basename(caminho_gz)
    logging.info(f"Extraindo: {nome_arquivo}")
    try:
        with gzip.open(caminho_gz, 'rb') as f_in:
            with open(caminho_fasta, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return True
    except (gzip.BadGzipFile, EOFError) as e:
        logging.error(f"ERRO CRITICO: O arquivo {nome_arquivo} esta corrompido ou incompleto. Pulando...")
        if os.path.exists(caminho_fasta):
            os.remove(caminho_fasta)
        return False

def rodar_mmseqs(caminho_fasta, pasta_saida):
    nome_arquivo_apenas = os.path.basename(caminho_fasta)
    logging.info(f"Agrupando proteinas de: {nome_arquivo_apenas} com MMseqs2...")
    
    # MMseqs2 exige uma pasta temporária própria para guardar os arquivos em banco de dados
    pasta_db_tmp = os.path.join(pasta_saida, "mmseqs_tmp")
    os.makedirs(pasta_db_tmp, exist_ok=True)
    
    nome_db = nome_arquivo_apenas + "_db"
    nome_db_clu = nome_arquivo_apenas + "_clu"
    nome_saida_tsv = nome_arquivo_apenas + "_mmseqs_out.tsv"
    
    # 1. Converte o FASTA em banco de dados interno do MMseqs2
    cmd_createdb = [EXECUTAVEL_LINUX, "mmseqs", "createdb", nome_arquivo_apenas, nome_db]
    
    # 2. Faz o clustering de alta sensibilidade com min_seq_id = 0.5 (50%)
    cmd_cluster = [EXECUTAVEL_LINUX, "mmseqs", "cluster", nome_db, nome_db_clu, "mmseqs_tmp", "--min-seq-id", "0.5", "-c", "0.8"]
    
    # 3. Transforma o resultado num arquivo TSV legível e bonito (cluster_id, sequence_id)
    cmd_createtsv = [EXECUTAVEL_LINUX, "mmseqs", "createtsv", nome_db, nome_db, nome_db_clu, nome_saida_tsv]

    try:
        subprocess.run(cmd_createdb, check=True, cwd=pasta_saida, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(cmd_cluster, check=True, cwd=pasta_saida, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(cmd_createtsv, check=True, cwd=pasta_saida, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        logging.error(f"Erro ao agrupar o arquivo {nome_arquivo_apenas} no MMseqs2.")
        return False

def extrair_dicionario_fasta(caminho_fasta):
    seqs = {}
    with open(caminho_fasta, 'r', encoding='utf-8') as f:
        header = ''
        seq = []
        for line in f:
            linha = line.strip()
            if linha.startswith('>'):
                if header:
                    seqs[header] = "".join(seq)
                header = linha.split()[0].replace('>', '') 
                seq = []
            else:
                seq.append(linha)
        if header:
            seqs[header] = "".join(seq)
    return seqs

def verificar_tipo_nutricional(sequencia):
    if not sequencia or sequencia == "Sequência não encontrada":
        return "N/A"
    
    letras_presentes = set(sequencia.upper())
    if AMINOACIDOS_ESSENCIAIS.issubset(letras_presentes):
        return "Completa"
    else:
        return "Incompleta"

def inicializar_csv():
    if not os.path.exists(ARQUIVO_CSV_FINAL):
        with open(ARQUIVO_CSV_FINAL, 'w', newline='', encoding='utf-8') as csv_out:
            escritor = csv.writer(csv_out, delimiter=';')
            escritor.writerow(['Arquivo_Origem', 'Grupo_Similaridade', 'ID_da_Proteina', 'Tipo_Nutricional', 'Sequencia_Completa_da_Proteina'])

def salvar_resultados_no_csv(pasta_temp, arquivo_origem, caminho_fasta_original):
    dicionario_sequencias = extrair_dicionario_fasta(caminho_fasta_original)
    
    with open(ARQUIVO_CSV_FINAL, 'a', newline='', encoding='utf-8') as csv_out:
        escritor = csv.writer(csv_out, delimiter=';')
        
        for item in os.listdir(pasta_temp):
            # O MMseqs gera arquivos tabulares .tsv maravilhosos (sem sujeira e confusão como o clstr)
            if item.endswith('_mmseqs_out.tsv'):
                caminho_completo = os.path.join(pasta_temp, item)
                
                with open(caminho_completo, 'r', encoding='utf-8') as f_in:
                    for linha in f_in:
                        linha_limpa = linha.strip()
                        if not linha_limpa:
                            continue
                        
                        # O TSV vem como: Representante (Cluster) [tab] Proteína_Membro
                        partes = linha_limpa.split('\t')
                        if len(partes) >= 2:
                            representante = partes[0]
                            membro = partes[1]
                            
                            cluster_atual = f"Grupo_{representante}"
                            id_proteina = membro
                                
                            seq_completa = dicionario_sequencias.get(id_proteina, "Sequência não encontrada")
                            tipo_nutri = verificar_tipo_nutricional(seq_completa)
                            
                            escritor.writerow([arquivo_origem, cluster_atual, id_proteina, tipo_nutri, seq_completa])

# ==========================================
# 3. FLUXO PRINCIPAL DE EXECUÇÃO
# ==========================================
def iniciar_automacao():
    setup_logging(ARQUIVO_LOG_SISTEMA)
    
    processados = carregar_arquivos_ja_processados()
    arquivos_gz = [f for f in os.listdir(DIR_GZ) if f.endswith('.gz')]

    logging.info(f"Total de arquivos .gz encontrados: {len(arquivos_gz)}")
    logging.info(f"Arquivos já processados anteriormente: {len(processados)}")
    
    inicializar_csv()

    for arquivo in arquivos_gz:
        if arquivo in processados:
            logging.info(f"Pulando {arquivo} (Já processado).")
            continue
            
        caminho_gz = os.path.join(DIR_GZ, arquivo)
        
        nome_base = arquivo.replace('.gz', '')
        if '.' in nome_base:
            nome_fas = nome_base.split('.')[0] + ".fas"
        else:
            nome_fas = nome_base + ".fas"
            
        caminho_fasta_local = os.path.join(DIR_TEMP_LOCAL, nome_fas)
        
        descompactado_com_sucesso = descompactar_gz(caminho_gz, caminho_fasta_local)
        
        if not descompactado_com_sucesso:
            continue
        
        sucesso = rodar_mmseqs(caminho_fasta_local, DIR_TEMP_LOCAL)
        
        if sucesso:
            salvar_resultados_no_csv(DIR_TEMP_LOCAL, arquivo, caminho_fasta_local)
            registrar_sucesso(arquivo)
            logging.info(f"Sucesso! {arquivo} consolidado no CSV final.")
            
        for item in os.listdir(DIR_TEMP_LOCAL):
            caminho_item = os.path.join(DIR_TEMP_LOCAL, item)
            try:
                if os.path.isdir(caminho_item):
                    shutil.rmtree(caminho_item)
                else:
                    os.remove(caminho_item)
            except Exception:
                pass

if __name__ == "__main__":
    iniciar_automacao()