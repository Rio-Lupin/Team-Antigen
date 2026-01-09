# Team Antigen - Architecture Documentation

## Table of Contents
1. [Overall Architecture](#overall-architecture)
2. [How Backend Components Work Together](#how-backend-components-work-together)
3. [API Endpoints Reference](#api-endpoints-reference)
4. [Data Flow Examples](#data-flow-examples)

---

## Overall Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                      USER INTERFACE                         │
│  (Web Browser - HTML/CSS/JavaScript)                       │
└───────────────────────┬─────────────────────────────────────┘
                        │ HTTP Requests (GET/POST)
                        ▼
┌─────────────────────────────────────────────────────────────┐
│                    Flask API Server                         │
│                      (app.py)                               │
│  • REST API Endpoints                                       │
│  • Request/Response Handling                                │
│  • Error Handling                                           │
└───────────────────────┬─────────────────────────────────────┘
                        │ Function Calls
                        ▼
┌─────────────────────────────────────────────────────────────┐
│              Chemistry Backend Module                       │
│              (backend/chemistry.py)                         │
│  • Drug Data Storage (drug_data, PROTEIN_TO_DRUGS)          │
│  • Molecule Processing (RDKit)                             │
│  • SVG Generation                                           │
│  • Drug Property Calculations                               │
│  • Molecule Modification                                    │
└───────────────────────┬─────────────────────────────────────┘
                        │
                        │ (Used by)
                        ▼
┌─────────────────────────────────────────────────────────────┐
│              CLI Testing Tool                               │
│              (backend/test-1.py)                            │
│  • Interactive Command-Line Tool                             │
│  • Displays Molecule Images                                 │
│  • Uses chemistry.py functions                              │
└─────────────────────────────────────────────────────────────┘
```

---

## How Backend Components Work Together

### 1. **chemistry.py** - The Core Module
**Purpose:** Central data storage and chemistry functions

**Key Components:**
- **Data Storage:**
  - `drug_data`: Dictionary mapping drug names → SMILES strings
  - `PROTEIN_TO_DRUGS`: Dictionary mapping protein names → list of drug names
  - `COMMON_R_GROUPS`: Dictionary of available R-group substitutions

- **Main Functions:**
  - `get_drugs_by_protein(protein_name)` → Returns list of drugs for a protein
  - `get_smiles_by_drug_name(drug_name)` → Returns SMILES string for a drug
  - `get_all_drugs()` → Returns all drugs in database
  - `get_modifiable_atoms(smiles)` → Returns list of atom indices that can be modified
  - `get_ui_metadata(smiles)` → Returns SVG, clickable atoms, coordinates, R-group options
  - `smiles_to_svg(smiles, ...)` → Generates SVG with clickable atom overlays
  - `validate_and_optimize_compound(smiles, atom_idx, new_group)` → Modifies molecule and validates
  - `apply_swap(smiles, atom_idx, new_group_smiles)` → Performs R-group substitution
  - `calculate_stability_score(mol)` → Calculates SAScore
  - `calculate_drug_likeness(mol)` → Calculates QED score
  - `check_lipinski_rules(mol)` → Validates Lipinski's Rule of Five
  - `check_structural_alerts(mol)` → Detects PAINS alerts

**Example:**
```python
from backend.chemistry import get_drugs_by_protein, get_ui_metadata

# Get drugs for a protein
drugs = get_drugs_by_protein("VEGFR")
# Returns: [{"name": "Sunitinib", "smiles": "CCN(CC)..."}, ...]

# Get molecule metadata with SVG
metadata = get_ui_metadata("CCN(CC)CCNC(=O)...")
# Returns: {
#   "molecule_svg": "<svg>...</svg>",
#   "modifiable_atoms": [0, 1, 3, ...],
#   "r_group_options": {"Methyl": "C", "Fluorine": "F", ...},
#   ...
# }
```

### 2. **test-1.py** - CLI Testing Tool
**Purpose:** Standalone command-line tool for testing/displaying molecules

**How it works:**
1. Asks user for protein input
2. Calls `get_drugs_by_protein()` from `chemistry.py`
3. Displays molecule images using RDKit's `Draw.MolToImage()`

**Key Point:** `test-1.py` does NOT duplicate data - it imports from `chemistry.py`

**Usage:**
```bash
python backend/test-1.py
# Enter: VEGFR
# Shows molecule images for Sunitinib and Cabozantinib
```

### 3. **app.py** - Flask API Server
**Purpose:** REST API that connects frontend to backend

**How it works:**
1. Receives HTTP requests from frontend
2. Calls functions from `chemistry.py`
3. Returns JSON responses

**Key Point:** All data comes from `chemistry.py` - no duplicate storage

---

## API Endpoints Reference

### 1. Get All Proteins
```
GET /api/proteins
```
**Response:**
```json
[
  {
    "id": "VEGFR",
    "drugs": [
      {"name": "Sunitinib", "compound": "placer1"},
      {"name": "Cabozantinib", "compound": "placer2"}
    ]
  },
  ...
]
```

### 2. Get Drugs by Protein
```
GET /api/drugs-by-protein/<protein_name>
```
**Example:** `GET /api/drugs-by-protein/VEGFR`

**Response:**
```json
{
  "protein": "VEGFR",
  "drugs": [
    {
      "name": "Sunitinib",
      "smiles": "CCN(CC)CCNC(=O)c1c(C)[nH]c(...)c1C"
    },
    {
      "name": "Cabozantinib",
      "smiles": "COc1cc2c(cc1OC)nccc2Oc3ccc(...)F"
    }
  ]
}
```

### 3. Get Molecule Metadata (by compound ID or drug name)
```
GET /api/molecule-metadata/<compound_id>
```
**Example:** 
- `GET /api/molecule-metadata/placer1` (compound ID)
- `GET /api/molecule-metadata/Sunitinib` (drug name)

**Response:**
```json
{
  "molecule_svg": "<svg>...</svg>",
  "modifiable_atoms": [0, 1, 3, 4, 5, ...],
  "atom_info": {
    "0": {"symbol": "C", "degree": 2, "is_in_ring": false, "is_aromatic": false},
    ...
  },
  "atom_coordinates": {
    "0": {"x": 44.76, "y": 124.83},
    ...
  },
  "r_group_options": {
    "Fluorine": "F",
    "Chlorine": "Cl",
    "Methyl": "C",
    ...
  },
  "stability_score": 4.2,
  "drug_likeness": 0.75,
  "lipinski_compliance": {
    "complies": true,
    "violations": 0,
    "molecular_weight": 398.5,
    "logP": 3.2,
    "h_donors": 2,
    "h_acceptors": 5
  },
  "structural_alerts": {
    "has_alert": false,
    "alert_reason": "None"
  }
}
```

### 4. Get Molecule Metadata (by SMILES)
```
GET /api/molecule-metadata-by-smiles?smiles=<SMILES_STRING>
```
**Example:** `GET /api/molecule-metadata-by-smiles?smiles=CCN(CC)CCNC(=O)...`

**Response:** Same as endpoint #3

### 5. Get All Drugs
```
GET /api/drugs
```
**Response:**
```json
{
  "drugs": {
    "Sunitinib": "CCN(CC)CCNC(=O)...",
    "Cabozantinib": "COc1cc2c(cc1OC)...",
    ...
  }
}
```

### 6. Modify Drug (Apply R-group substitution)
```
POST /api/modify-drug
Content-Type: application/json

{
  "smiles": "CCN(CC)CCNC(=O)...",
  "atom_idx": 5,
  "new_group": "Methyl"
}
```

**Response:**
```json
{
  "success": true,
  "new_smiles": "CCN(CC)CCNC(=O)...",
  "molecule_svg": "<svg>...</svg>",
  "stability_score": 4.2,
  "drug_likeness": 0.75,
  "lipinski_compliance": {
    "complies": true,
    "violations": 0,
    "molecular_weight": 398.5,
    "logP": 3.2,
    "h_donors": 2,
    "h_acceptors": 5
  },
  "structural_alerts": {
    "has_alert": false,
    "alert_reason": "None"
  },
  "is_real_compound": false,
  "compound_name": null
}
```

**Available R-groups:** `"Methyl"`, `"Fluorine"`, `"Chlorine"`, `"Bromine"`, `"Iodine"`, `"Hydroxyl"`, `"Amino"`, `"Phenyl"`, `"Nitro"`, `"Cyano"`, `"Trifluoromethyl"`, `"Methoxy"`, `"Ethoxy"`, `"Aminomethyl"`, `"Carboxyl"`, `"Aldehyde"`, `"Tert-butyl"`, `"Ethyl"`

(See `COMMON_R_GROUPS` in `chemistry.py` for full list)

**Error Response:**
```json
{
  "success": false,
  "error": "Failed to apply swap. Invalid atom index or incompatible substitution.",
  "new_smiles": null,
  ...
}
```

---

## Data Flow Examples

### Example 1: User Views Drugs for a Protein

```
1. Frontend: User types "VEGFR" in search field
   ↓
2. Frontend: fetch('/api/drugs-by-protein/VEGFR')
   ↓
3. app.py: get_drugs_for_protein("VEGFR")
   ↓
4. chemistry.py: get_drugs_by_protein("VEGFR")
   ↓
5. chemistry.py: Looks up PROTEIN_TO_DRUGS["VEGFR"]
   ↓
6. chemistry.py: Retrieves SMILES from drug_data for each drug
   ↓
7. chemistry.py: Returns [{"name": "Sunitinib", "smiles": "..."}, ...]
   ↓
8. app.py: Returns JSON response
   ↓
9. Frontend: Displays drug cards with molecule SVGs
```

### Example 2: User Clicks on a Molecule Atom and Applies R-group Substitution

```
1. Frontend: User clicks atom #5 in SVG (clickable circle with data-atom-id="5")
   ↓
2. Frontend: showRGroupDropdown() - displays dropdown with R-group options
   ↓
3. Frontend: User selects "Methyl" from dropdown
   ↓
4. Frontend: fetch('/api/modify-drug', {
     method: 'POST',
     body: JSON.stringify({
       smiles: currentSmiles,
       atom_idx: 5,
       new_group: "Methyl"
     })
   })
   ↓
5. app.py: modify_drug() receives request
   ↓
6. app.py: Validates required parameters (smiles, atom_idx, new_group)
   ↓
7. chemistry.py: validate_and_optimize_compound(smiles, 5, "Methyl")
   ↓
8. chemistry.py: Resolves "Methyl" → "C" from COMMON_R_GROUPS
   ↓
9. chemistry.py: apply_swap() modifies molecule at atom index 5
   ↓
10. chemistry.py: Validates new SMILES string
    ↓
11. chemistry.py: Calculates properties:
    - calculate_stability_score() → SAScore
    - calculate_drug_likeness() → QED
    - check_lipinski_rules() → Lipinski compliance
    - check_structural_alerts() → PAINS detection
    ↓
12. chemistry.py: get_modifiable_atoms() for new molecule
    ↓
13. chemistry.py: smiles_to_svg() generates SVG with:
    - Clickable atom overlays
    - Highlighted modified atom
    - Correct coordinates using RDKit GetDrawCoords()
    ↓
14. chemistry.py: Returns result dictionary with:
    - success: true
    - new_smiles
    - molecule_svg
    - All calculated properties
    ↓
15. app.py: Returns JSON response
    ↓
16. Frontend: Updates molecule display with new SVG
    ↓
17. Frontend: Updates drug name to "DrugName-test"
    ↓
18. Frontend: Updates SMILES string
    ↓
19. Frontend: Updates analytics panel with new scores
    ↓
20. Frontend: Re-attaches click handlers to new SVG for further modifications
```

### Example 3: Molecule Visualization with Clickable Atoms

```
1. Frontend: Requests molecule metadata
   GET /api/molecule-metadata/Sunitinib
   ↓
2. app.py: get_molecule_ui_data("Sunitinib")
   ↓
3. chemistry.py: get_ui_metadata(smiles)
   ↓
4. chemistry.py: get_modifiable_atoms(smiles) → [0, 1, 3, 4, 5, ...]
   ↓
5. chemistry.py: smiles_to_svg(smiles, modifiable_atoms=[0, 1, 3, ...])
   ↓
6. RDKit: Generates 2D coordinates
   ↓
7. RDKit: Draws molecule to SVG
   ↓
8. chemistry.py: For each modifiable atom:
    - Calls drawer.GetDrawCoords(atom_idx) to get screen coordinates
    - Adds clickable circle overlay at exact position
    - Sets data-atom-id attribute for JavaScript interaction
   ↓
9. chemistry.py: Returns complete SVG with embedded clickable circles
   ↓
10. Frontend: Inserts SVG into DOM
    ↓
11. Frontend: attachAtomClickHandlers() finds all .clickable-atom elements
    ↓
12. Frontend: Each circle has click handler that shows R-group dropdown
```

---

## File Structure Summary

```
Team-Antigen/
├── app.py                    # Flask API server (entry point)
├── backend/
│   ├── chemistry.py         # Core chemistry module (data + functions)
│   └── test-1.py            # CLI testing tool (optional)
├── templates/
│   └── base.html            # Frontend HTML template
├── requirements.txt         # Python dependencies
├── README.md               # User documentation
└── ARCHITECTURE.md         # This file
```

---

## Key Technical Details

### Coordinate System for Clickable Atoms

- Uses RDKit's `GetDrawCoords()` method to get exact screen coordinates after drawing
- Coordinates are stored as `(x, y)` tuples in viewBox coordinate system
- Clickable circles are positioned using `cx` and `cy` attributes in SVG
- Ensures perfect alignment between rendered atoms and clickable overlays

### Molecule Modification

- Uses RDKit's `ReplaceSubstructs()` for R-group substitution
- Validates chemical valency to ensure modifications are chemically valid
- Automatically recalculates all properties after modification
- Preserves modifiable atom information in modified molecules

### Property Calculations

- **SAScore**: Uses RDKit's synthetic accessibility scoring
- **QED**: Uses RDKit's Quantitative Estimate of Drug-likeness
- **Lipinski's Rule of Five**: Calculates MW, LogP, HBD, HBA and checks compliance
- **PAINS**: Uses RDKit's FilterCatalog to detect problematic substructures

### Error Handling

- All API endpoints return consistent error format: `{"error": "message"}`
- Invalid SMILES strings return appropriate error messages
- Failed modifications return `success: false` with error description
- Frontend handles errors gracefully with user-friendly messages

---

## Dependencies

### Python Packages
- `flask` - Web framework
- `rdkit` - Cheminformatics toolkit
- `rdkit-pypi` - RDKit Python bindings (if using PyPI installation)

### Optional
- `chembl-webresource-client` - For ChEMBL database queries (optional feature)

---

## Development Notes

- `test-1.py` is for testing/development only - not used by the web app
- `COMPOUND_DATABASE` in `app.py` is legacy - prefer using drug names from `chemistry.py`
- Molecule modifications preserve valency (chemical correctness)
- All validation scores are calculated automatically after modifications
- SVG coordinates use RDKit's internal coordinate system for accuracy
