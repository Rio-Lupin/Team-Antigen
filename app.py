from flask import Flask, render_template, jsonify
from chemistry import get_ui_metadata, validate_and_optimize_compound

app = Flask(__name__)


#Protein Data
proteins = [
        {
        "id": "VEGF Receptor 2",
        "drugs": [
            {"name": "Sunitinib-", "compound": "placer1"},
            {"name": "Cabozantinib","compound": "placer2"}]
        }
        ,
        {
          "id": "PDGF Receptor beta" ,
          "drugs": [
            {"name": "Sunitinib-", "compound": "placer3"},
            {"name": "Cabozantinib","compound": "placer4"}]
        },
        {
        "id": "mTOR",
        "drugs": [{"name": "Everolimus", "compound": "placer5"}, {"name": "Temsirolimus","compound": "placer1"}]
        }
]

COMPOUND_DATABASE = {
    "placer1": "CN1CC(C2=C(C1)C=CC(=C2)F)C3=C(NC(=C3)C)C=O", # Sunitinib-like
    "placer2": "COC1=CC2=C(C=C1OCC3CCN(CC3)C4=CC=C(C=C4)F)N=CN=C2NC5=CC(=C(C=C5)Cl)F",
    "placer3": "Cc1c(C[NH+]2CCCC2)sc(NC(=O)Nc2ccc(Cl)cc2)c1",
    # ... add others
}
        

@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")

#Getting the protein - Nyosha 
@app.route('/api/proteins', methods=['GET'])
def protein_and_drugs_jsonify():
    return jsonify(proteins)
    
@app.route('/api/molecule-metadata/<compound_id>', methods=['GET'])
def get_molecule_ui_data(compound_id):
    """Provides the SVG, clickable atom indices, and coordinates."""
    smiles = COMPOUND_DATABASE.get(compound_id)
    if not smiles:
        return jsonify({"error": "Compound not found"}), 404
    
    try:
        # This calls the get_ui_metadata function from your backend module
        metadata = get_ui_metadata(smiles)
        return jsonify(metadata)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/api/modify-drug', methods=['POST'])
def modify_drug():
    """Receives the swap request and returns new scores + new SVG."""
    data = request.json
    
    # Expected payload: { "smiles": "...", "atom_idx": 5, "new_group": "Methyl" }
    old_smiles = data.get("smiles")
    atom_idx = int(data.get("atom_idx"))
    new_group = data.get("new_group")
    
    if not all([old_smiles, atom_idx is not None, new_group]):
        return jsonify({"error": "Missing parameters"}), 400
        
    # Call the main validation function
    result = validate_and_optimize_compound(old_smiles, atom_idx, new_group)
    
    # This result already contains stability_score, drug_likeness, and molecule_svg
    return jsonify(result)