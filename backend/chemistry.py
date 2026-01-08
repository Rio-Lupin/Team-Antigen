from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Contrib.SA_Score import sascorer

# took from demo
def smiles_to_svg(smiles: str, width: int = 400, height: int = 400) -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES")

    Chem.rdCoordGen.AddCoords(mol)
    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    # set drawing options on drawer.getOptions()
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()

def get_synthesis_score(smiles: str):
    mol = Chem.MolFromSmiles(smiles) # Convert SMILES to RDKit Molecule object
    if mol:
        score = sascorer.calculateScore(mol) # Calculate the synthetic accessibility score
        return round(score, 2)
    return None

def calculate_drug_likeness(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        score = QED.qed(mol) # Calculate the QED score
        return round(score, 3)
    return None


drug_list=["Sunitinib", "Cabozantinib"]


print("Enter Protein:")
protein = input()

if protein == "VEGFR" or protein == "vegfr":
        print( drug_list[0], drug_list[1])

elif protein == "PDGFR" or protein == "pdgfr":
        print( drug_list[0], drug_list[1])
    
else:
    print("No drugs found for this protein")
