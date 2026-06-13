import os
import gzip
import sys
import pandas as pd
from Bio import SeqIO
from io import TextIOWrapper

# --- CONFIGURAÇÕES DE DIRETÓRIO ---

download_path = "C:/Users/Filipe Camargo Prado/Downloads/codigos_gerar_genes"
CAMINHO_HISTORICO = os.path.join(download_path, "historico_processamento.txt")

g_drive = "G:\Meu Drive\base_genes_tcc"
PASTA_FONTE = os.path.join(g_drive, "arquivos_gz")
ARQUIVO_CSV_FINAL = os.path.join(g_drive, "base_proteinas_consolidada.csv")

def processar_arquivos_ftp_para_csv():
    """
    Varre a pasta local, lê os arquivos .gz sem descompactar no disco
    e extrai os dados das proteínas direto para um CSV unificado com controle de progresso.
    """
    arquivos_gz = [f for f in os.listdir(PASTA_FONTE) if f.endswith('.gz')]
    
    if not arquivos_gz:
        print(f"⚠️ Nenhum arquivo .gz encontrado na pasta: {PASTA_FONTE}")
        return

    # Carrega o histórico de arquivos que já foram processados com sucesso em rodadas anteriores
    processados_anteriormente = set()
    if os.path.exists(CAMINHO_HISTORICO):
        with open(CAMINHO_HISTORICO, "r", encoding="utf-8") as h:
            processados_anteriormente = set(linha.strip() for list_item in h for linha in [list_item] if linha.strip())

    print(f"📚 Encontrados {len(arquivos_gz)} arquivos no total. {len(processados_anteriormente)} já foram consolidados.")
    
    total_proteinas_processadas = 0
    
    for idx, nome_arquivo in enumerate(arquivos_gz, start=1):
        # TRAVA DE SEGURANÇA: Se o arquivo já foi processado com sucesso antes, pula para não duplicar no CSV
        if nome_arquivo in processados_anteriormente:
            continue
            
        caminho_arquivo_gz = os.path.join(PASTA_FONTE, nome_arquivo)
        print(f"\n📦 [{idx}/{len(arquivos_gz)}] Lendo arquivo: {nome_arquivo}")
        
        registros_lote = []
        contador_local = 0
        
        try:
            with gzip.open(caminho_arquivo_gz, 'rb') as f_bytes:
                with TextIOWrapper(f_bytes, encoding='utf-8') as f_texto:
                    
                    for registro in SeqIO.parse(f_texto, "fasta"):
                        contador_local += 1
                        total_proteinas_processadas += 1
                        
                        if total_proteinas_processadas % 1000 == 0:
                            sys.stdout.write(f"\r🔄 Processados nesta rodada: {total_proteinas_processadas} proteínas...")
                            sys.stdout.flush()
                        
                        seq = str(registro.seq).upper()
                        gc_val = 100 * float(seq.count("G") + seq.count("C")) / len(seq) if len(seq) > 0 else 0
                        
                        registros_lote.append({
                            "ID_NCBI": registro.id,
                            "Definicao": registro.description[:100],
                            "Tamanho_AA": len(seq),
                            "GC_Percent": f"{gc_val:.2f}%",
                            "Sequencia": seq
                        })
                        
                        if len(registros_lote) >= 50000:
                            salvar_lote_no_csv(registros_lote)
                            registros_lote = []
            
            if registros_lote:
                salvar_lote_no_csv(registros_lote)
                
            print(f"\n   ✅ Concluído! {contador_local} proteínas extraídas deste arquivo.")
            
            # SALVA NO HISTÓRICO: Arquivo lido com sucesso de ponta a ponta, marca como feito
            with open(CAMINHO_HISTORICO, "a", encoding="utf-8") as h:
                h.write(f"{nome_arquivo}\n")
            
        except Exception as e:
            print(f"\n❌ Erro crítico ao processar o arquivo {nome_arquivo}: {e}")
            if "decompressing" in str(e) or "invalid" in str(e):
                print(f"🗑️ Deletando arquivo corrompido para re-download: {nome_arquivo}")
                try:
                    os.remove(caminho_arquivo_gz)
                except:
                    pass
            continue

    print(f"\n🎉 SUCESSO ABSOLUTO! Todas as bases disponíveis foram consolidadas.")
    print(f"💾 Arquivo final gerado em: {ARQUIVO_CSV_FINAL}")

def salvar_lote_no_csv(lista_dados):
    """
    Salva uma lista de dicionários diretamente no arquivo CSV no modo Append.
    """
    df_lote = pd.DataFrame(lista_dados)
    existe_csv = os.path.exists(ARQUIVO_CSV_FINAL)
    
    # Grava no modo 'a' (append) garantindo escrita ultra rápida
    df_lote.to_csv(ARQUIVO_CSV_FINAL, mode='a', index=False, header=not existe_csv, encoding='utf-8')

if __name__ == "__main__":
    # Garante que o terminal trate corretamente caracteres especiais na saída visual
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding='utf-8')
        
    print("🚀 Iniciando a esteira de processamento de Big Data do TCC...")
    processar_arquivos_ftp_para_csv()