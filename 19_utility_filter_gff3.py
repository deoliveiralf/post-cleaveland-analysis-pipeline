#!/usr/bin/env python3
# filter_gff3.py - Filtra o GFF3 baseado nos transcritos de interesse

import pandas as pd
import os
import re
import sys
import argparse
import logging
from collections import defaultdict
from datetime import datetime

def setup_logging(output_dir):
    """Configura o sistema de logging"""
    log_file = os.path.join(output_dir, f'filter_gff3_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def load_target_transcripts(mirna_target_file):
    """Carrega a lista de transcritos alvo"""
    logger = logging.getLogger(__name__)
    try:
        df = pd.read_csv(mirna_target_file, sep='\t')
        transcripts = set(df['Target ID'].dropna().unique())
        logger.info(f"Total de transcritos alvo únicos: {len(transcripts)}")
        return transcripts
    except Exception as e:
        logger.error(f"Erro ao carregar arquivo miRNA-target: {e}")
        raise

def filter_gff3(input_gff3, target_transcripts, output_gff3):
    """Filtra o arquivo GFF3 mantendo apenas os transcritos de interesse"""
    logger = logging.getLogger(__name__)
    mrna_ids = set()
    transcripts_found = set()
    
    try:
        with open(input_gff3, 'r') as infile, open(output_gff3, 'w') as outfile:
            # Primeira passagem: identificar mRNAs de interesse
            for line in infile:
                if line.startswith('#') or not line.strip():
                    continue
                
                cols = line.strip().split('\t')
                if len(cols) < 9 or cols[2] != "mRNA":
                    continue
                
                attrs = cols[8]
                tx_id = re.search(r'ID=([^;]+)', attrs).group(1)
                tx_name = re.search(r'Name=([^;]+)', attrs).group(1) if 'Name=' in attrs else tx_id
                
                if tx_name in target_transcripts or tx_id in target_transcripts:
                    mrna_ids.add(tx_id)
                    transcripts_found.add(tx_name)
                    outfile.write(line)
            
            # Segunda passagem: escrever features dos mRNAs selecionados
            infile.seek(0)
            for line in infile:
                if line.startswith('#') or not line.strip():
                    continue
                
                cols = line.strip().split('\t')
                if len(cols) < 9:
                    continue
                
                # Escreve mRNAs já escritos na primeira passagem
                if cols[2] == "mRNA":
                    tx_id = re.search(r'ID=([^;]+)', cols[8]).group(1)
                    if tx_id in mrna_ids:
                        continue
                
                # Escreve features filhos dos mRNAs selecionados
                if 'Parent=' in cols[8]:
                    parents = re.findall(r'Parent=([^;]+)', cols[8])
                    if any(parent in mrna_ids for parent in parents):
                        outfile.write(line)
        
        logger.info(f"Transcritos encontrados no GFF3: {len(transcripts_found)}/{len(target_transcripts)}")
        logger.info(f"GFF3 filtrado salvo em: {output_gff3}")
        
        # Verifica se todos os transcritos foram encontrados
        missing = target_transcripts - transcripts_found
        if missing:
            logger.warning(f"Transcritos não encontrados no GFF3: {len(missing)}")
            with open(os.path.join(os.path.dirname(output_gff3), 'missing_transcripts.txt'), 'w') as f:
                f.write("\n".join(sorted(missing)))
        
        return output_gff3
        
    except Exception as e:
        logger.error(f"Erro durante o filtro do GFF3: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Filtra arquivo GFF3 baseado nos transcritos alvo")
    parser.add_argument('input_gff3', help='Arquivo GFF3 original')
    parser.add_argument('mirna_target_table', help='Tabela de pares miRNA-target')
    parser.add_argument('-o', '--output', default='filtered_annotation.gff3', help='Arquivo GFF3 filtrado de saída')
    parser.add_argument('--output-dir', default='.', help='Diretório de saída')
    
    args = parser.parse_args()
    
    # Configuração inicial
    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logging(args.output_dir)
    output_gff3 = os.path.join(args.output_dir, args.output)
    
    try:
        logger.info("Iniciando filtro do GFF3")
        
        # Etapa 1: Carrega transcritos alvo
        target_transcripts = load_target_transcripts(args.mirna_target_table)
        
        # Etapa 2: Filtra GFF3
        filter_gff3(args.input_gff3, target_transcripts, output_gff3)
        
        logger.info("Processo concluído com sucesso!")
        
    except Exception as e:
        logger.error(f"Erro fatal: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
