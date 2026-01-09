from rdkit import Chem  # type: ignore[import-untyped]
from rdkit.Chem import Draw  # type: ignore[import-untyped]
from backend.chemistry import get_drugs_by_protein

print("Enter Protein (e.g., VEGFR or PDGFR):")
protein = input().strip().upper()

# Get drugs for the protein using chemistry.py function
drugs = get_drugs_by_protein(protein)

if drugs:
    for drug in drugs:
        name = drug["name"]
        smiles = drug["smiles"]

        # 1. Print the text info for both drugs
        print(f"\nDrug: {name}")
        print(f"SMILES: {smiles}")

        # 2. Generate the molecule object
        mol = Chem.MolFromSmiles(smiles)

        if mol:
            # 3. Create the image and open it in your default viewer
            img = Draw.MolToImage(mol, size=(400, 400))
            img.show()  # This replaces the need for IPython's display()
        else:
            print(f"Error: Could not process SMILES for {name}")
else:
    print("No drugs found for this protein.")
