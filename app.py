from flask import Flask, render_template, jsonify, request  # type: ignore[import-untyped]
from backend.chemistry import (
    get_ui_metadata,
    validate_and_optimize_compound,
    get_drugs_by_protein,
    get_smiles_by_drug_name,
    get_all_drugs,
)

app = Flask(__name__)


# Protein Data
proteins = [
    {
        "id": "VEGFR",
        "drugs": [
            {"name": "Sunitinib", "compound": "placer1"},
            {"name": "Cabozantinib", "compound": "placer2"},
        ],
    },
    {
        "id": "PDGFR",
        "drugs": [
            {"name": "Sunitinib", "compound": "placer1"},
            {"name": "Cabozantinib", "compound": "placer2"},
        ],
    },
    {
        "id": "mTOR",
        "drugs": [
            {"name": "Everolimus", "compound": "placer3"},
            {"name": "Temsirolimus", "compound": "placer4"},
        ],
    },
]

COMPOUND_DATABASE = {
    "placer1": "CN1CC(C2=C(C1)C=CC(=C2)F)C3=C(NC(=C3)C)C=O",  # Sunitinib-like
    "placer2": "COC1=CC2=C(C=C1OCC3CCN(CC3)C4=CC=C(C=C4)F)N=CN=C2NC5=CC(=C(C=C5)Cl)F",
    "placer3": "Cc1c(C[NH+]2CCCC2)sc(NC(=O)Nc2ccc(Cl)cc2)c1",
    "placer4": "COC1CC(CCC1OC(=O)C(CO)(CO)C)CC(C1OC(=O)C2CCCCN2C(=O)C(=O)C2(O)OC(CCC2C)CC(OC)C(=CC=CC=CC(CC(C(=O)C(C(C(=CC(C(=O)C1)C)C)O)OC)C)C)C)C"
}


@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")


# Getting the protein - Nyosha
@app.route("/api/proteins", methods=["GET"])
def protein_and_drugs_jsonify():
    return jsonify(proteins)


@app.route("/api/molecule-metadata/<compound_id>", methods=["GET"])
def get_molecule_ui_data(compound_id):
    """Provides the SVG, clickable atom indices, and coordinates."""
    # First try COMPOUND_DATABASE (for compound IDs like "placer1")
    smiles = COMPOUND_DATABASE.get(compound_id)

    # If not found, try drug name lookup (from test-1.py integration)
    if not smiles:
        smiles = get_smiles_by_drug_name(compound_id)

    if not smiles:
        return jsonify({"error": "Compound not found"}), 404

    try:
        # calls the get_ui_metadata function from backend module
        metadata = get_ui_metadata(smiles)
        return jsonify(metadata)
    except Exception as e:
        app.logger.exception(
            "Error while generating UI metadata for compound_id %s with SMILES %s",
            compound_id,
            smiles,
        )
        return (
            jsonify(
                {"error": "Internal server error while generating molecule metadata"}
            ),
            500,
        )


@app.route("/api/drugs-by-protein/<protein_name>", methods=["GET"])
def get_drugs_for_protein(protein_name):
    """Get all drugs associated with a protein name.

    Returns list of drugs with their names and SMILES strings.
    """
    try:
        drugs = get_drugs_by_protein(protein_name)
        if not drugs:
            return (
                jsonify(
                    {
                        "error": f"No drugs found for protein: {protein_name}",
                        "drugs": [],
                    }
                ),
                404,
            )

        return jsonify({"protein": protein_name, "drugs": drugs})
    except Exception as e:
        app.logger.exception("Error getting drugs for protein %s", protein_name)
        return jsonify({"error": "Internal server error"}), 500


@app.route("/api/drugs", methods=["GET"])
def get_all_drugs_endpoint():
    """Get all available drugs in the database."""
    try:
        drugs = get_all_drugs()
        return jsonify({"drugs": drugs})
    except Exception as e:
        app.logger.exception("Error getting all drugs")
        return jsonify({"error": "Internal server error"}), 500


@app.route("/api/molecule-metadata-by-smiles", methods=["GET"])
def get_molecule_metadata_by_smiles():
    """Get molecule metadata by providing SMILES directly as query parameter.

    Query params:
        smiles: SMILES string of the molecule
    """
    smiles = request.args.get("smiles")
    if not smiles:
        return jsonify({"error": "Missing 'smiles' query parameter"}), 400

    try:
        metadata = get_ui_metadata(smiles)
        return jsonify(metadata)
    except Exception as e:
        app.logger.exception("Error generating metadata for SMILES: %s", smiles)
        return jsonify({"error": f"Error generating molecule metadata: {str(e)}"}), 500


@app.route("/api/modify-drug", methods=["POST"])
def modify_drug():
    """Receives the swap request and returns new scores + new SVG."""
    data = request.json or {}

    # Expected payload: { "smiles": "...", "atom_idx": 5, "new_group": "Methyl" }
    old_smiles = data.get("smiles")
    raw_atom_idx = data.get("atom_idx")
    new_group = data.get("new_group")

    # Validate required parameters
    if not old_smiles or raw_atom_idx is None or not new_group:
        return (
            jsonify(
                {"error": "Missing required parameters: smiles, atom_idx, new_group"}
            ),
            400,
        )

    # Validate and parse atom_idx
    try:
        atom_idx = int(raw_atom_idx)
    except (TypeError, ValueError):
        return jsonify({"error": "Invalid atom_idx: must be an integer"}), 400
        
    # main validation function
    try:
        result = validate_and_optimize_compound(old_smiles, atom_idx, new_group)
        return jsonify(result)
    except ValueError as e:
        app.logger.warning("Validation error for SMILES %s: %s", old_smiles, str(e))
        return jsonify({"error": f"Validation failed: {str(e)}"}), 400
    except Exception as e:
        app.logger.exception(
            "Error validating compound with SMILES %s, atom_idx %d, new_group %s",
            old_smiles,
            atom_idx,
            new_group,
        )
        return (
            jsonify({"error": "Internal server error while processing compound"}),
            500,
        )
