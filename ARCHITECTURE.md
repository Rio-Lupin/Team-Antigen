# Team Antigen - Architecture & Frontend Guide

## 📋 Table of Contents
1. [Overall Architecture](#overall-architecture)
2. [How Backend Components Work Together](#how-backend-components-work-together)
3. [API Endpoints Reference](#api-endpoints-reference)
4. [Frontend Developer Guide](#frontend-developer-guide)
5. [Data Flow Examples](#data-flow-examples)

---

## 🏗️ Overall Architecture

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

- **Main Functions:**
  - `get_drugs_by_protein(protein_name)` → Returns list of drugs for a protein
  - `get_smiles_by_drug_name(drug_name)` → Returns SMILES string for a drug
  - `get_all_drugs()` → Returns all drugs in database
  - `get_ui_metadata(smiles)` → Returns SVG, clickable atoms, coordinates
  - `validate_and_optimize_compound(smiles, atom_idx, new_group)` → Modifies molecule and validates

**Example:**
```python
from backend.chemistry import get_drugs_by_protein

drugs = get_drugs_by_protein("VEGFR")
# Returns: [{"name": "Sunitinib", "smiles": "CCN(CC)..."}, ...]
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

##  API Endpoints Reference

### 1. Get All Proteins
```
GET /api/proteins
```
**Response:**
```json
[
  {
    "id": "VEGF Receptor 2",
    "drugs": [
      {"name": "Sunitinib-", "compound": "placer1"},
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
  "modifiable_atoms": [5, 12, 18, 23],
  "atom_coordinates": {
    "5": {"x": 123.4, "y": 456.7},
    "12": {"x": 234.5, "y": 567.8},
    ...
  },
  "smiles": "CCN(CC)CCNC(=O)..."
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
    "Cabozantinib": "COc1cc2c(cc1OC)..."
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
  "new_smiles": "CCN(CC)CCNC(=O)...",
  "molecule_svg": "<svg>...</svg>",
  "stability_score": 4.2,
  "drug_likeness": 0.75,
  "lipinski": {
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

**Available R-groups:** `"Methyl"`, `"Fluorine"`, `"Chlorine"`, `"Hydroxyl"`, `"Amino"`, `"Phenyl"`, etc.
(See `COMMON_R_GROUPS` in `chemistry.py` for full list)

---

##  Frontend Developer Guide

### What You Need to Know

#### 1. **Displaying Molecules**

The API returns SVG strings that can be directly inserted into HTML:

```javascript
// Fetch molecule metadata
fetch('/api/molecule-metadata/Sunitinib')
  .then(res => res.json())
  .then(data => {
    // Insert SVG directly into DOM
    document.getElementById('molecule-container').innerHTML = data.molecule_svg;
    
    // Attach click handlers to modifiable atoms
    attachAtomClickHandlers(data.modifiable_atoms, data.atom_coordinates);
  });
```

#### 2. **Clickable Atoms**

The SVG includes clickable overlay circles for modifiable atoms. Each circle has:
- `class="clickable-atom"`
- `data-atom-idx="5"` (the atom index)
- Transparent fill with hover effects

**JavaScript Example:**
```javascript
function attachAtomClickHandlers(modifiableAtoms, coordinates) {
  modifiableAtoms.forEach(atomIdx => {
    const circle = document.querySelector(`[data-atom-idx="${atomIdx}"]`);
    if (circle) {
      circle.addEventListener('click', () => {
        showRGroupModal(atomIdx, coordinates[atomIdx]);
      });
    }
  });
}
```

#### 3. **R-group Selection Modal**

When user clicks an atom, show a modal with R-group options:

```javascript
function showRGroupModal(atomIdx, coordinates) {
  const rGroups = [
    "Methyl", "Fluorine", "Chlorine", "Hydroxyl", 
    "Amino", "Phenyl", "Nitro", "Cyano"
  ];
  
  // Show modal with dropdown/buttons for R-groups
  // When user selects one, call modifyDrug()
}
```

#### 4. **Modifying Molecules**

Send modification request to backend:

```javascript
async function modifyDrug(smiles, atomIdx, newGroup) {
  const response = await fetch('/api/modify-drug', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      smiles: smiles,
      atom_idx: atomIdx,
      new_group: newGroup
    })
  });
  
  const result = await response.json();
  
  // Update UI with:
  // - result.molecule_svg (new molecule image)
  // - result.stability_score
  // - result.drug_likeness
  // - result.lipinski (compliance info)
  // - result.structural_alerts
}
```

#### 5. **Displaying Validation Results**

Show the validation scores in a results panel:

```javascript
function displayValidationResults(result) {
  // Stability Score (lower is better, scale 1-10)
  document.getElementById('stability').textContent = result.stability_score;
  
  // Drug-likeness (0-1, higher is better)
  document.getElementById('drug-likeness').textContent = 
    (result.drug_likeness * 100).toFixed(1) + '%';
  
  // Lipinski's Rule of Five
  const lipinski = result.lipinski;
  document.getElementById('lipinski-status').textContent = 
    lipinski.complies ? '✓ Compliant' : '✗ Violations: ' + lipinski.violations;
  document.getElementById('molecular-weight').textContent = 
    lipinski.molecular_weight + ' Da';
  document.getElementById('logp').textContent = lipinski.logP;
  
  // Structural Alerts
  if (result.structural_alerts.has_alert) {
    document.getElementById('alerts').textContent = 
      '⚠️ ' + result.structural_alerts.alert_reason;
  }
}
```

#### 6. **Complete Workflow Example**

```javascript
// 1. User selects a protein
async function loadDrugsForProtein(proteinName) {
  const response = await fetch(`/api/drugs-by-protein/${proteinName}`);
  const data = await response.json();
  
  // Display list of drugs
  data.drugs.forEach(drug => {
    displayDrugCard(drug.name, drug.smiles);
  });
}

// 2. User clicks on a drug
async function loadMolecule(drugName) {
  const response = await fetch(`/api/molecule-metadata/${drugName}`);
  const data = await response.json();
  
  // Display SVG
  document.getElementById('molecule').innerHTML = data.molecule_svg;
  
  // Make atoms clickable
  attachAtomClickHandlers(data.modifiable_atoms, data.atom_coordinates);
  
  // Store current SMILES for modification
  window.currentSmiles = data.smiles;
}

// 3. User clicks an atom and selects R-group
async function applyModification(atomIdx, rGroup) {
  const response = await fetch('/api/modify-drug', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      smiles: window.currentSmiles,
      atom_idx: atomIdx,
      new_group: rGroup
    })
  });
  
  const result = await response.json();
  
  // Update molecule display
  document.getElementById('molecule').innerHTML = result.molecule_svg;
  
  // Show validation results
  displayValidationResults(result);
  
  // Update current SMILES
  window.currentSmiles = result.new_smiles;
}
```

---

##  Data Flow Examples

### Example 1: User Views Drugs for a Protein

```
1. Frontend: User types "VEGFR" and clicks search
   ↓
2. Frontend: fetch('/api/drugs-by-protein/VEGFR')
   ↓
3. app.py: get_drugs_for_protein("VEGFR")
   ↓
4. chemistry.py: get_drugs_by_protein("VEGFR")
   ↓
5. chemistry.py: Looks up PROTEIN_TO_DRUGS["VEGFR"]
   ↓
6. chemistry.py: Returns [{"name": "Sunitinib", "smiles": "..."}, ...]
   ↓
7. app.py: Returns JSON response
   ↓
8. Frontend: Displays drug cards
```

### Example 2: User Clicks on a Molecule Atom

```
1. Frontend: User clicks atom #5 in SVG
   ↓
2. Frontend: showRGroupModal(5, coordinates[5])
   ↓
3. Frontend: User selects "Methyl" from dropdown
   ↓
4. Frontend: fetch('/api/modify-drug', {smiles, atom_idx: 5, new_group: "Methyl"})
   ↓
5. app.py: modify_drug() receives request
   ↓
6. chemistry.py: validate_and_optimize_compound(smiles, 5, "Methyl")
   ↓
7. chemistry.py: apply_swap() modifies molecule
   ↓
8. chemistry.py: Calculates scores (stability, drug-likeness, Lipinski, PAINS)
   ↓
9. chemistry.py: Generates new SVG
   ↓
10. app.py: Returns JSON with new_smiles, molecule_svg, scores
    ↓
11. Frontend: Updates molecule display and shows validation results
```

---

## File Structure Summary

```
Team-Antigen/
├── app.py                    # Flask API server
├── backend/
│   ├── chemistry.py         # Core chemistry module (data + functions)
│   └── test-1.py            # CLI testing tool
├── templates/
│   └── base.html            # Frontend HTML/JS
└── requirements.txt         # Python dependencies
```

---

##  Key Points for Frontend Dev

1. **All molecule data comes from `chemistry.py`** - no hardcoding needed
2. **SVG strings are ready to use** - just insert into DOM with `innerHTML`
3. **Clickable atoms are already in SVG** - just attach event listeners
4. **R-group names match exactly** - use names from `COMMON_R_GROUPS` dictionary
5. **Atom indices are 0-based** - first atom is index 0
6. **All API responses are JSON** - easy to parse and use
7. **Error handling** - check for `error` field in responses

---

## 📝 Notes

- `test-1.py` is for testing only - not used by the web app
- `COMPOUND_DATABASE` in `app.py` is legacy - prefer using drug names from `chemistry.py`
- Molecule modifications preserve valency (chemical correctness)
- All validation scores are calculated automatically
