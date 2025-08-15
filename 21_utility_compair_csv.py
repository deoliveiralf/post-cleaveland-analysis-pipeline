#!/usr/bin/env python3
# comparador_csv.py - Compara resultados de múltiplas execuções do cleavage_annotator

import pandas as pd
from pathlib import Path
import argparse

def comparar_csvs(arquivos, output_comparacao):
    """
    Compara múltiplos arquivos CSV gerados pelo cleavage_annotator
    e gera um relatório das diferenças.
    """
    # Carrega todos os arquivos
    dfs = {}
    for arquivo in arquivos:
        nome = Path(arquivo).stem
        dfs[nome] = pd.read_csv(arquivo)
    
    # Verifica se temos dados para comparar
    if len(dfs) < 2:
        print("Erro: Necessário pelo menos 2 arquivos para comparação")
        return
    
    # Cria um relatório de comparação
    with open(output_comparacao, 'w') as f:
        # 1. Comparação básica
        f.write("=== COMPARAÇÃO BÁSICA ===\n")
        for nome, df in dfs.items():
            f.write(f"{nome}: {len(df)} linhas\n")
        
        # 2. Diferenças nas linhas
        f.write("\n=== LINHAS ÚNICAS EM CADDA ARQUIVO ===\n")
        todos = pd.concat(dfs.values()).drop_duplicates(keep=False)
        
        for nome, df in dfs.items():
            unicas = todos[todos.isin(df.to_dict('list')).all(axis=1)]
            f.write(f"\n{nome} tem {len(unicas)} linhas únicas:\n")
            f.write(unicas.to_string() + "\n")
        
        # 3. Comparação estatística
        f.write("\n=== ESTATÍSTICAS COMPARATIVAS ===\n")
        stats = {}
        for nome, df in dfs.items():
            stats[nome] = {
                'Total': len(df),
                'miRNAs Únicos': df['miRNA'].nunique(),
                'Transcripts Únicos': df['Transcript'].nunique(),
                'Regiões Distintas': df['Region'].value_counts().to_dict()
            }
        
        for nome, dados in stats.items():
            f.write(f"\nEstatísticas - {nome}:\n")
            for k, v in dados.items():
                f.write(f"{k}: {v}\n")
        
        # 4. Exemplo de diferenças (mostra as primeiras 5)
        f.write("\n=== EXEMPLOS DE DIFERENÇAS (5 primeiras) ===\n")
        base = next(iter(dfs.values()))
        outros = [df for nome, df in dfs.items() if nome != next(iter(dfs))]
        
        for i, outro in enumerate(outros, 1):
            diff = pd.concat([base, outro]).drop_duplicates(keep=False)
            f.write(f"\nDiferenças entre {next(iter(dfs))} e {list(dfs.keys())[i]}:\n")
            f.write(diff.head(5).to_string() + "\n")
    
    print(f"Relatório de comparação salvo em: {output_comparacao}")

def main():
    parser = argparse.ArgumentParser(description='Compara múltiplos arquivos CSV do cleavage_annotator')
    parser.add_argument('arquivos', nargs='+', help='Arquivos CSV para comparar')
    parser.add_argument('-o', '--output', default='comparacao.txt', help='Arquivo de saída para o relatório')
    args = parser.parse_args()

    comparar_csvs(args.arquivos, args.output)

if __name__ == "__main__":
    main()
