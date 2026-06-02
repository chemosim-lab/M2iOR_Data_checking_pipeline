import argparse
from pathlib import Path
from scripts.process_excel import process_excel
from scripts.fetch_uniprot import fetch_sequences
from scripts.align_sequences import find_mutations
import pandas as pd

LARAVEL_DATA_PATH = Path("../laravel/database/seeders/data")

def process_group(excel_path: Path) -> None:
    group_name = excel_path.stem           # ex: "group_001"
    output_path = LARAVEL_DATA_PATH / group_name
    output_path.mkdir(parents=True, exist_ok=True)

    # 1. Lire et nettoyer l'Excel
    df = process_excel(excel_path)

    # 2. Récupérer les séquences de référence via UniProt
    df = fetch_sequences(df)

    # 3. Calculer les mutations
    df_experiments, df_mutations = find_mutations(df)

    # 4. Extraire les récepteurs uniques
    df_receptors = df[['uniprot_accession']].drop_duplicates()

    # 5. Écrire les CSV directement dans Laravel
    df_receptors.to_csv(   output_path / "receptors.csv",    index=False)
    df_experiments.to_csv( output_path / "experiments.csv",  index=False)
    df_mutations.to_csv(   output_path / "mutations.csv",    index=False)

    print(f"✓ {group_name} → {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file",  help="Traiter un seul fichier Excel")
    parser.add_argument("--all",   action="store_true", help="Traiter tous les Excel")
    args = parser.parse_args()

    input_path = Path("input")

    if args.file:
        process_group(input_path / args.file)
    elif args.all:
        for excel_file in sorted(input_path.glob("*.xlsx")):
            process_group(excel_file)
