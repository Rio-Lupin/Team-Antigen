from rdkit import Chem
from rdkit.Chem import Draw

# Dictionary containing both drugs and their SMILES

drug_data = {
    "Sunitinib": "CCN(CC)CCNC(=O)c1c(C)[nH]c(\C=C2/C(=O)Nc3ccc(F)cc23)c1C",
    "Cabozantinib": "COc1cc2c(cc1OC)nccc2Oc3ccc(cc3)NC(=O)C4(CC4)C(=O)Nc5ccc(cc5)F"
}

print("Enter Protein (e.g., VEGFR or PDGFR):")
protein = input().strip().upper()

if protein in ["VEGFR", "PDGFR"]:
    for name, smiles in drug_data.items():
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
    