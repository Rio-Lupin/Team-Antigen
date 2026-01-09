"""
Chemoinformatics Module for Drug Discovery Workflow
Provides functions for molecule modification, validation, and visualization.
"""

from typing import Dict, List, Optional, Tuple, Any
import logging
import re
from rdkit import Chem  # type: ignore[import-untyped]
from rdkit.Chem import QED, AllChem, Descriptors, rdMolDescriptors, FilterCatalog  # type: ignore[import-untyped]
from rdkit.Chem import rdCoordGen  # type: ignore[import-untyped]
from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore[import-untyped]
from rdkit import rdBase  # type: ignore[import-untyped]

# Try multiple import paths for SA_Score (varies by installation method)
# Note: SA_Score is in RDKit's Contrib directory which may not be in Python path
SA_SCORE_AVAILABLE = False
sascorer = None

try:
    from rdkit.Contrib.SA_Score import sascorer  # type: ignore[import-untyped]

    SA_SCORE_AVAILABLE = True
except ImportError:
    try:
        # Alternative import path for Conda installations
        import sys
        import os

        contrib_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "share", "RDKit", "Contrib"
        )
        if os.path.exists(contrib_path):
            sys.path.insert(0, contrib_path)
            from SA_Score import sascorer  # type: ignore[import-untyped]

            SA_SCORE_AVAILABLE = True
    except ImportError:
        logging.warning("SA_Score not available. Stability scoring will be disabled.")
        SA_SCORE_AVAILABLE = False
        sascorer = None

try:
    from chembl_webresource_client.new_client import new_client  # type: ignore[import-untyped]

    CHEMBL_AVAILABLE = True
except ImportError:
    logging.warning(
        "chembl_webresource_client not available. ChEMBL checks will be disabled."
    )
    CHEMBL_AVAILABLE = False

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Drug data: Maps drug names to their SMILES strings
# Integrated from test-1.py
drug_data = {
    "Sunitinib": "CCN(CC)CCNC(=O)c1c(C)[nH]c(\\C=C2/C(=O)Nc3ccc(F)cc23)c1C",
    "Cabozantinib": "COc1cc2c(cc1OC)nccc2Oc3ccc(cc3)NC(=O)C4(CC4)C(=O)Nc5ccc(cc5)F",
    "Everolimus": "CC1CCC2CC(C(=CC=CC=CC(CC(C(=O)C(C(C(=CC(C(=O)CC(OC(=O)C3CCCCN3C(=O)C(=O)C1(O2)O)C(C)CC4CCC(C(C4)OC)OCCO)C)C)O)OC)C)C)C)OC",
    "Temsirolimus": "COC1CC(CCC1OC(=O)C(CO)(CO)C)CC(C1OC(=O)C2CCCCN2C(=O)C(=O)C2(O)OC(CCC2C)CC(OC)C(=CC=CC=CC(CC(C(=O)C(C(C(=CC(C(=O)C1)C)C)O)OC)C)C)C)C",
}

# Protein to drug mapping
# Maps protein names (case-insensitive) to lists of drug names
PROTEIN_TO_DRUGS = {
    "VEGFR": ["Sunitinib", "Cabozantinib"],
    "PDGFR": ["Sunitinib", "Cabozantinib"],
    "VEGF RECEPTOR 2": ["Sunitinib", "Cabozantinib"],
    "PDGF RECEPTOR BETA": ["Sunitinib", "Cabozantinib"],
    "MTOR": ["Everolimus", "Temsirolimus"],
}

# Common R-group substitution options for the dropdown
# Format: {display_name: SMILES_string}
# Note: These are written to be attached via a single bond
# The attachment point is implicit (where the group connects)
COMMON_R_GROUPS = {
    "Fluorine": "F",
    "Chlorine": "Cl",
    "Bromine": "Br",
    "Iodine": "I",
    "Methyl": "C",
    "Ethyl": "CC",
    "Hydroxyl": "O",
    "Amino": "N",
    "Nitro": "[N+](=O)[O-]",
    "Cyano": "C#N",
    "Trifluoromethyl": "C(F)(F)F",
    "Methoxy": "OC",
    "Ethoxy": "OCC",
    "Aminomethyl": "CN",
    "Carboxyl": "C(=O)O",
    "Aldehyde": "C=O",
    "Phenyl": "c1ccccc1",  # Aromatic ring
    "Tert-butyl": "CC(C)(C)C",
}


def get_drugs_by_protein(protein_name: str) -> List[Dict[str, str]]:
    """
    Get list of drugs associated with a protein.

    Args:
        protein_name: Name of the protein (case-insensitive)

    Returns:
        List of dictionaries with 'name' and 'smiles' keys for each drug
        Returns empty list if protein not found
    """
    protein_upper = protein_name.strip().upper()

    # Try exact match first
    drug_names = PROTEIN_TO_DRUGS.get(protein_upper, [])

    # If no exact match, try partial matching
    if not drug_names:
        for key, drugs in PROTEIN_TO_DRUGS.items():
            if protein_upper in key or key in protein_upper:
                drug_names = drugs
                break

    # Build result list with drug names and SMILES
    result = []
    for drug_name in drug_names:
        smiles = drug_data.get(drug_name)
        if smiles:
            result.append({"name": drug_name, "smiles": smiles})

    return result


def get_smiles_by_drug_name(drug_name: str) -> Optional[str]:
    """
    Get SMILES string for a drug by name.

    Args:
        drug_name: Name of the drug

    Returns:
        SMILES string if found, None otherwise
    """
    return drug_data.get(drug_name)


def get_all_drugs() -> Dict[str, str]:
    """
    Get all drugs in the database.

    Returns:
        Dictionary mapping drug names to SMILES strings
    """
    return drug_data.copy()


def get_modifiable_atoms(smiles: str) -> List[int]:
    """
    Identifies atoms that are safe to modify (peripheral/R-groups, not core ring atoms).

    Args:
        smiles: SMILES string of the molecule

    Returns:
        List of atom indices that can be safely modified

    Raises:
        ValueError: If the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    modifiable_indices = []

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        symbol = atom.GetSymbol()

        # Skip if atom is in a ring (core structure)
        if atom.IsInRing():
            # Allow halogens in rings (common substitution sites)
            if symbol not in ["F", "Cl", "Br", "I"]:
                continue

        # Allow terminal atoms (only one neighbor)
        if atom.GetDegree() == 1:
            modifiable_indices.append(atom_idx)
            continue

        # Allow halogens regardless of connectivity
        if symbol in ["F", "Cl", "Br", "I"]:
            modifiable_indices.append(atom_idx)
            continue

        # Allow carbon/nitrogen/oxygen with low connectivity (peripheral groups)
        if symbol in ["C", "N", "O"] and atom.GetDegree() <= 2:
            # Check if it's not part of an aromatic ring system
            if not atom.GetIsAromatic():
                modifiable_indices.append(atom_idx)

    return modifiable_indices


def apply_swap(base_smiles: str, atom_idx: int, new_group_smiles: str) -> Optional[str]:
    """
    Replaces an atom in the molecule with a new group, handling valency correctly.

    Uses RDKit reactions to ensure proper bond handling and valency.

    Args:
        base_smiles: Original SMILES string
        atom_idx: Index of the atom to replace (0-based)
        new_group_smiles: SMILES string of the replacement group

    Returns:
        New SMILES string if successful, None otherwise
    """
    try:
        # Parse the base molecule
        mol = Chem.MolFromSmiles(base_smiles)
        if mol is None:
            logger.error(f"Invalid base SMILES: {base_smiles}")
            return None

        # Validate atom index
        if atom_idx < 0 or atom_idx >= mol.GetNumAtoms():
            logger.error(
                f"Invalid atom index: {atom_idx} (molecule has {mol.GetNumAtoms()} atoms)"
            )
            return None

        atom_to_replace = mol.GetAtomWithIdx(atom_idx)
        neighbors = list(atom_to_replace.GetNeighbors())

        # If it's a terminal atom (single neighbor), we can do a simple replacement
        if len(neighbors) == 1:
            return _replace_terminal_atom(mol, atom_idx, new_group_smiles)
        else:
            # For atoms with multiple neighbors, we need a more sophisticated approach
            return _replace_multi_neighbor_atom(mol, atom_idx, new_group_smiles)

    except Exception as e:
        logger.error(f"Error in apply_swap: {str(e)}")
        return None


def _replace_terminal_atom(
    mol: Chem.Mol, atom_idx: int, new_group_smiles: str
) -> Optional[str]:
    """Replace a terminal atom (only one neighbor) with a new group."""
    try:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = list(atom.GetNeighbors())
        if len(neighbors) != 1:
            return None

        neighbor = neighbors[0]
        neighbor_idx = neighbor.GetIdx()
        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), atom_idx)
        bond_order = bond.GetBondTypeAsDouble()

        # Parse the new group
        new_group_mol = Chem.MolFromSmiles(new_group_smiles)
        if new_group_mol is None:
            logger.error(f"Invalid new group SMILES: {new_group_smiles}")
            return None

        # For simple atoms like F, Cl, etc., use direct replacement
        if new_group_mol.GetNumAtoms() == 1:
            new_symbol = new_group_mol.GetAtomWithIdx(0).GetSymbol()
            editable_mol = Chem.RWMol(mol)
            new_atom = Chem.Atom(new_symbol)
            editable_mol.ReplaceAtom(atom_idx, new_atom)
            final_mol = editable_mol.GetMol()

            try:
                Chem.SanitizeMol(final_mol)
                return Chem.MolToSmiles(final_mol)
            except:
                return None

        # For multi-atom groups, use reaction-based replacement
        # Strategy: Cut the bond, then attach the new group
        # Reaction pattern: [*:1]-[*:2] -> [*:1]-[NewGroup]
        try:
            # First, try using a reaction with proper attachment point
            # The new group needs to have an attachment point marked
            # For now, we'll try attaching at the first atom of the new group

            # Create a modified new group SMILES with attachment point
            # We assume the new group attaches via its first atom
            reaction_smarts = f"[*:1]-[*:2]>>[*:1]-{new_group_smiles}"
            rxn = AllChem.ReactionFromSmarts(reaction_smarts)

            if rxn is not None:
                # Prepare molecule with mapping
                mol_with_map = Chem.RWMol(mol)
                mol_with_map.GetAtomWithIdx(neighbor_idx).SetAtomMapNum(1)
                mol_with_map.GetAtomWithIdx(atom_idx).SetAtomMapNum(2)
                mol_with_map = mol_with_map.GetMol()

                # Run reaction
                products = rxn.RunReactants((mol_with_map,))

                if products and len(products) > 0:
                    product = products[0][0]
                    # Remove atom map numbers
                    for atom in product.GetAtoms():
                        atom.SetAtomMapNum(0)
                    Chem.SanitizeMol(product)
                    return Chem.MolToSmiles(product)
        except Exception as e:
            logger.debug(
                f"Reaction-based replacement failed: {str(e)}, trying alternative method"
            )

        # Alternative: Manual bond cutting and group attachment
        try:
            editable_mol = Chem.RWMol(mol)
            stored_neighbor_idx = neighbor_idx

            # Parse the new group
            new_group_mol = Chem.MolFromSmiles(new_group_smiles)
            if new_group_mol:
                # Add the new group to the molecule
                new_group_atoms = []
                for atom in new_group_mol.GetAtoms():
                    new_atom_idx = editable_mol.AddAtom(atom)
                    new_group_atoms.append(new_atom_idx)

                # Add bonds from the new group
                for bond in new_group_mol.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()
                    editable_mol.AddBond(
                        new_group_atoms[begin_idx],
                        new_group_atoms[end_idx],
                        bond.GetBondType(),
                    )

                # Connect the new group to the neighbor (attach at first atom)
                editable_mol.AddBond(
                    stored_neighbor_idx, new_group_atoms[0], Chem.BondType.SINGLE
                )

                # Remove the old atom
                editable_mol.RemoveAtom(atom_idx)

                final_mol = editable_mol.GetMol()
                Chem.SanitizeMol(final_mol)
                return Chem.MolToSmiles(final_mol)
        except Exception as e:
            logger.debug(f"Manual replacement failed: {str(e)}")

        # Final fallback
        return _direct_replacement_fallback(mol, atom_idx, new_group_smiles)

    except Exception as e:
        logger.error(f"Error in _replace_terminal_atom: {str(e)}")
        return None


def _replace_multi_neighbor_atom(
    mol: Chem.Mol, atom_idx: int, new_group_smiles: str
) -> Optional[str]:
    """Replace an atom with multiple neighbors (more complex case)."""
    # For hackathon purposes, we'll be more conservative here
    # Only allow simple atom replacements for multi-neighbor atoms
    try:
        new_group_mol = Chem.MolFromSmiles(new_group_smiles)
        if new_group_mol is None or new_group_mol.GetNumAtoms() != 1:
            logger.warning(
                f"Cannot replace multi-neighbor atom with complex group: {new_group_smiles}"
            )
            return None

        new_symbol = new_group_mol.GetAtomWithIdx(0).GetSymbol()
        editable_mol = Chem.RWMol(mol)
        new_atom = Chem.Atom(new_symbol)
        editable_mol.ReplaceAtom(atom_idx, new_atom)
        final_mol = editable_mol.GetMol()

        try:
            Chem.SanitizeMol(final_mol)
            return Chem.MolToSmiles(final_mol)
        except Exception as e:
            logger.warning(f"Sanitization failed: {str(e)}")
            return None

    except Exception as e:
        logger.error(f"Error in _replace_multi_neighbor_atom: {str(e)}")
        return None


def _direct_replacement_fallback(
    mol: Chem.Mol, atom_idx: int, new_group_smiles: str
) -> Optional[str]:
    """Fallback method using direct atom type replacement."""
    try:
        new_group_mol = Chem.MolFromSmiles(new_group_smiles)
        if new_group_mol is None:
            return None

        if new_group_mol.GetNumAtoms() == 1:
            editable_mol = Chem.RWMol(mol)
            new_symbol = new_group_mol.GetAtomWithIdx(0).GetSymbol()
            new_atom = Chem.Atom(new_symbol)
            editable_mol.ReplaceAtom(atom_idx, new_atom)
            final_mol = editable_mol.GetMol()

            try:
                Chem.SanitizeMol(final_mol)
                return Chem.MolToSmiles(final_mol)
            except:
                return None

        return None
    except:
        return None


def calculate_stability_score(mol: Chem.Mol) -> Optional[float]:
    """
    Calculate the synthetic accessibility score (SAScore) for a molecule.
    Lower scores indicate easier synthesis.

    Args:
        mol: RDKit molecule object

    Returns:
        SAScore (typically 1-10, lower is better) or None if invalid
    """
    if not SA_SCORE_AVAILABLE or sascorer is None:
        logger.warning("SA_Score not available. Cannot calculate stability score.")
        return None
    try:
        score = sascorer.calculateScore(mol)
        return round(score, 2)
    except Exception as e:
        logger.error(f"Error calculating stability score: {str(e)}")
        return None


def calculate_drug_likeness(mol: Chem.Mol) -> Optional[float]:
    """
    Calculate the Quantitative Estimate of Drug-likeness (QED) score.
    Higher scores (closer to 1.0) indicate better drug-likeness.

    Args:
        mol: RDKit molecule object

    Returns:
        QED score (0-1, higher is better) or None if invalid
    """
    try:
        score = QED.qed(mol)
        return round(score, 3)
    except Exception as e:
        logger.error(f"Error calculating drug-likeness: {str(e)}")
        return None


def check_lipinski_rules(mol: Chem.Mol) -> Dict[str, Any]:
    """Checks Lipinski's Rule of Five compliance.

    Args:
        mol: RDKit molecule object
    Returns:
        Dictionary with compliance status and property values

    """
    try:
        mw = Descriptors.MolWt(mol)  # Molecular weight, not exact mass
        logp = Descriptors.MolLogP(mol)  # Standard LogP calculation
        h_donors = Descriptors.NumHDonors(mol)  # Correct function name
        h_acceptors = Descriptors.NumHAcceptors(mol)  # Correct function name
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)  # Number of rotatable bonds

        violations = sum([mw > 500, logp > 5, h_donors > 5, h_acceptors > 10, rotatable_bonds > 10])

        complies = violations == 0

        return {
            "complies": complies,
            "violations": violations,
            "molecular_weight": round(mw, 2),
            "logP": round(logp, 2),
            "h_donors": h_donors,
            "h_acceptors": h_acceptors,
            "rotatable_bonds": rotatable_bonds,
        }
    except Exception as e:
        logger.error(f"Lipinski check failed: {e}")
        return {}


def check_structural_alerts(mol: Chem.Mol) -> Dict[str, Any]:
    """Checks Structural Alerts (PAINS).

    Args:
        mol: RDKit molecule object

    Returns:
        Dictionary with alert status and reason
    """
    try:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        has_alert = catalog.HasMatch(mol)
        alert_reason = (
            catalog.GetFirstMatch(mol).GetDescription() if has_alert else "None"
        )

        return {"has_alert": has_alert, "alert_reason": alert_reason}
    except Exception as e:
        logger.error(f"PAINS alert: {str(e)}")
        return {}


def verify_if_real_compound(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Check if a compound exists in the ChEMBL database.

    Args:
        smiles: SMILES string to search for

    Returns:
        Tuple of (is_real_compound: bool, compound_name: Optional[str])

    Raises:
        RuntimeError: If ChEMBL client is not available
    """
    if not CHEMBL_AVAILABLE:
        raise RuntimeError("ChEMBL webresource client is not available")

    try:
        molecule = new_client.molecule
        # Canonicalize the SMILES first for better matching
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, None

        canon_smiles = Chem.MolToSmiles(mol)

        # Search ChEMBL database
        res = molecule.filter(
            molecule_structures__canonical_smiles__flex=canon_smiles
        ).only(["pref_name", "molecule_chembl_id"])

        if len(res) > 0:
            # Get the preferred name, or fallback to ChEMBL ID
            pref_name = res[0].get("pref_name")
            if not pref_name:
                pref_name = res[0].get("molecule_chembl_id", "Unknown")
            return True, pref_name

        return False, None

    except Exception as e:
        logger.error(f"Error checking ChEMBL database: {str(e)}")
        # Return False on error to avoid breaking the workflow
        return False, None


def smiles_to_svg(
    smiles: str,
    width: int = 400,
    height: int = 400,
    highlight_atom_idx: Optional[int] = None,
    modifiable_atoms: Optional[List[int]] = None,
) -> str:
    """
    Generate an SVG representation of a molecule.

    Args:
        smiles: SMILES string of the molecule
        width: SVG width in pixels
        height: SVG height in pixels
        highlight_atom_idx: Optional atom index to highlight (for modified atoms)
        modifiable_atoms: Optional list of atom indices that should be clickable

    Returns:
        SVG string (not bytes) for easy web display

    Raises:
        ValueError: If the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    try:
        # Generate 2D coordinates
        rdCoordGen.AddCoords(mol)
    except:
        # Fallback to default coordinate generation
        AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    # Get drawer options to understand coordinate system
    drawer_opts = drawer.drawOptions()

    # Configure highlighting if atom index is provided
    highlight_atoms = []
    highlight_bonds = []
    highlight_radii = {}
    highlight_colors = {}

    if highlight_atom_idx is not None and 0 <= highlight_atom_idx < mol.GetNumAtoms():
        highlight_atoms = [highlight_atom_idx]
        # Highlight bonds connected to this atom
        atom = mol.GetAtomWithIdx(highlight_atom_idx)
        for neighbor in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond:
                highlight_bonds.append(bond.GetIdx())
        highlight_radii[highlight_atom_idx] = 0.5
        # Use a bright color for highlighting (yellow)
        highlight_colors[highlight_atom_idx] = (1.0, 1.0, 0.0)

    # Draw the molecule
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=highlight_bonds,
        highlightAtomRadii=highlight_radii,
        highlightAtomColors=highlight_colors,
    )
    drawer.FinishDrawing()

    svg_text = drawer.GetDrawingText()

    # Get actual drawn coordinates using RDKit's GetDrawCoords
    # This gives us the exact screen positions after drawing (most accurate method)
    draw_coords_dict = {}
    if modifiable_atoms:
        try:
            # GetDrawCoords returns Point2D objects with x and y attributes
            # These are the actual screen coordinates RDKit uses
            for atom_idx in modifiable_atoms:
                if 0 <= atom_idx < mol.GetNumAtoms():
                    try:
                        pos = drawer.GetDrawCoords(atom_idx)
                        # Point2D has x and y attributes
                        draw_coords_dict[atom_idx] = (float(pos.x), float(pos.y))
                    except Exception as e:
                        logger.debug(
                            f"Could not get draw coords for atom {atom_idx}: {str(e)}"
                        )
                        pass
        except Exception as e:
            logger.debug(f"Could not use GetDrawCoords: {str(e)}")

    # Add clickable overlays for modifiable atoms
    if modifiable_atoms:
        svg_text = _add_clickable_atoms_to_svg(
            svg_text,
            mol,
            modifiable_atoms,
            width,
            height,
            drawer,
            drawer_opts,
            draw_coords_dict,
        )

    return svg_text


def _add_clickable_atoms_to_svg(
    svg_text: str,
    mol: Chem.Mol,
    modifiable_atoms: List[int],
    width: int,
    height: int,
    drawer=None,
    drawer_opts=None,
    draw_coords_dict: Optional[Dict[int, Tuple[float, float]]] = None,
) -> str:
    """
    Add clickable overlay circles to SVG for modifiable atoms.
    Inserts transparent circles with data attributes for JavaScript interaction.

    Extracts actual atom positions from SVG text elements that RDKit draws,
    which provides the most accurate coordinates matching the rendered atoms.
    """
    try:
        clickable_circles = []
        conf = mol.GetConformer()

        # First, try to use RDKit's GetDrawCoords if available (most accurate)
        if draw_coords_dict:
            for atom_idx in modifiable_atoms:
                if atom_idx in draw_coords_dict:
                    svg_x_pos, svg_y_pos = draw_coords_dict[atom_idx]
                    clickable_circles.append(
                        f'<circle cx="{svg_x_pos:.2f}" cy="{svg_y_pos:.2f}" r="15" '
                        f'fill="rgba(43, 140, 255, 0.2)" stroke="#2b8cff" stroke-width="2" '
                        f'cursor="pointer" data-atom-id="{atom_idx}" class="clickable-atom" '
                        f'style="pointer-events: all;"/>'
                    )

        # Fallback: Calculate from molecule coordinates if SVG parsing didn't work
        if not clickable_circles:
            conf = mol.GetConformer()
            coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            if not coords:
                return svg_text

            # Parse SVG viewBox (RDKit may or may not include it)
            viewbox_match = re.search(r'viewBox="([^"]*)"', svg_text)
            if viewbox_match:
                viewbox_parts = viewbox_match.group(1).split()
                if len(viewbox_parts) >= 4:
                    vb_x = float(viewbox_parts[0])
                    vb_y = float(viewbox_parts[1])
                    vb_width = float(viewbox_parts[2])
                    vb_height = float(viewbox_parts[3])
                else:
                    vb_x, vb_y, vb_width, vb_height = 0, 0, width, height
            else:
                # No viewBox - RDKit uses pixel coordinates directly
                vb_x, vb_y, vb_width, vb_height = 0, 0, width, height

            # Get molecule bounds
            min_x = min(c.x for c in coords)
            max_x = max(c.x for c in coords)
            min_y = min(c.y for c in coords)
            max_y = max(c.y for c in coords)
            mol_width = max_x - min_x if max_x != min_x else 1
            mol_height = max_y - min_y if max_y != min_y else 1

            # RDKit drawer uses padding (default 0.05 or 5% margin on each side)
            # Get actual padding from drawer options if available
            padding = 0.05
            if drawer_opts is not None:
                try:
                    padding = getattr(drawer_opts, "padding", 0.05)
                except:
                    pass
            elif drawer is not None:
                try:
                    opts = drawer.drawOptions()
                    padding = getattr(opts, "padding", 0.05)
                except:
                    pass

            # Calculate available drawing area
            if vb_width > 0:
                available_width = vb_width * (1 - 2 * padding)
                available_height = vb_height * (1 - 2 * padding)
                draw_center_x = vb_x + vb_width / 2
                draw_center_y = vb_y + vb_height / 2
            else:
                # No viewBox - use pixel dimensions
                available_width = width * (1 - 2 * padding)
                available_height = height * (1 - 2 * padding)
                draw_center_x = width / 2
                draw_center_y = height / 2

            # Calculate scale (RDKit uses the smaller scale to maintain aspect ratio)
            scale_x = available_width / mol_width if mol_width > 0 else 1
            scale_y = available_height / mol_height if mol_height > 0 else 1
            scale = min(scale_x, scale_y)

            # RDKit centers the molecule in the drawing area
            mol_center_x = (min_x + max_x) / 2
            mol_center_y = (min_y + max_y) / 2

            # Calculate offset to center molecule
            # RDKit centers the molecule in the drawing area
            # SVG and RDKit both use Y-axis increasing downward (no flip needed)
            offset_x = draw_center_x - mol_center_x * scale
            offset_y = draw_center_y - mol_center_y * scale

            # Create clickable circles
            for atom_idx in modifiable_atoms:
                if 0 <= atom_idx < len(coords):
                    pos = coords[atom_idx]
                    svg_x_pos = offset_x + pos.x * scale
                    svg_y_pos = offset_y + pos.y * scale

                    clickable_circles.append(
                        f'<circle cx="{svg_x_pos:.2f}" cy="{svg_y_pos:.2f}" r="15" '
                        f'fill="rgba(43, 140, 255, 0.2)" stroke="#2b8cff" stroke-width="2" '
                        f'cursor="pointer" data-atom-id="{atom_idx}" class="clickable-atom" '
                        f'style="pointer-events: all;"/>'
                    )

        # Insert the clickable circles before closing </svg> tag
        if clickable_circles:
            circles_svg = "\n".join(clickable_circles)
            svg_text = svg_text.replace("</svg>", f"{circles_svg}\n</svg>")

        return svg_text
    except Exception as e:
        logger.warning(f"Could not add clickable atoms to SVG: {str(e)}")
        return svg_text


def get_ui_metadata(smiles: str) -> Dict[str, Any]:
    """
    Generate UI metadata for the frontend.
    Returns information about which atoms are clickable and what dropdown options are available.

    Args:
        smiles: SMILES string of the molecule

    Returns:
        Dictionary containing:
            - modifiable_atoms: List of atom indices that can be modified
            - atom_info: Dictionary mapping atom indices to their properties
            - atom_coordinates: Dictionary mapping atom indices to their (x, y) coordinates
            - r_group_options: Available R-group substitution options
            - molecule_svg: SVG representation of the molecule with clickable overlays
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Generate coordinates if not present
        try:
            rdCoordGen.AddCoords(mol)
        except:
            AllChem.Compute2DCoords(mol)

        modifiable_atoms = get_modifiable_atoms(smiles)

        # Get atom information for modifiable atoms
        atom_info = {}
        atom_coordinates = {}
        conf = mol.GetConformer()

        for idx in modifiable_atoms:
            atom = mol.GetAtomWithIdx(idx)
            atom_info[idx] = {
                "symbol": atom.GetSymbol(),
                "degree": atom.GetDegree(),
                "is_in_ring": atom.IsInRing(),
                "is_aromatic": atom.GetIsAromatic(),
            }
            # Get coordinates
            pos = conf.GetAtomPosition(idx)
            atom_coordinates[idx] = {
                "x": float(pos.x),
                "y": float(pos.y),
            }

        # Generate SVG with clickable overlays
        molecule_svg = smiles_to_svg(smiles, modifiable_atoms=modifiable_atoms)

        return {
            "modifiable_atoms": modifiable_atoms,
            "atom_info": atom_info,
            "atom_coordinates": atom_coordinates,
            "r_group_options": COMMON_R_GROUPS,
            "molecule_svg": molecule_svg,
        }

    except Exception as e:
        logger.error(f"Error generating UI metadata: {str(e)}")
        raise


def _find_modified_atom_in_new_molecule(
    new_smiles: str, new_group_smiles: str, original_atom_idx: int
) -> Optional[int]:
    """
    Attempt to find the modified atom in the new molecule.
    This is a heuristic approach for highlighting purposes.

    Args:
        new_smiles: SMILES of the modified molecule
        new_group_smiles: SMILES of the replacement group
        original_atom_idx: Original atom index (for reference)

    Returns:
        Atom index in new molecule if found, None otherwise
    """
    try:
        new_mol = Chem.MolFromSmiles(new_smiles)
        if new_mol is None:
            return None

        # Parse the new group to get its symbol/pattern
        group_mol = Chem.MolFromSmiles(new_group_smiles)
        if group_mol is None:
            return None

        # If it's a single atom, search for that atom type in terminal positions
        if group_mol.GetNumAtoms() == 1:
            target_symbol = group_mol.GetAtomWithIdx(0).GetSymbol()
            # Look for terminal atoms with matching symbol
            for atom in new_mol.GetAtoms():
                if atom.GetSymbol() == target_symbol and atom.GetDegree() == 1:
                    return atom.GetIdx()

        # Fallback: return None (no highlighting)
        return None
    except:
        return None


def validate_and_optimize_compound(
    old_smiles: str, atom_idx: int, new_group: str
) -> Dict[str, Any]:
    """
    Main function for Function 2: Modify a drug and validate the result.

    This function:
    1. Applies the swap to create a new compound
    2. Validates the compound (stability, drug-likeness, ChEMBL check)
    3. Returns visualization with highlighted modified atom

    Args:
        old_smiles: Original SMILES string
        atom_idx: Index of atom to modify
        new_group: SMILES string of the replacement group (key from COMMON_R_GROUPS)

    Returns:
        Dictionary containing:
            - success: bool indicating if the swap was successful
            - new_smiles: New SMILES string (if successful)
            - stability_score: SAScore
            - drug_likeness: QED score
            - is_real_compound: bool
            - compound_name: Name if found in ChEMBL
            - molecule_svg: SVG with highlighted atom
            - error: Error message if failed
    """
    result = {
        "success": False,
        "new_smiles": None,
        "stability_score": None,
        "drug_likeness": None,
        "lipinski_compliance": None,
        "structural_alerts": None,
        "is_real_compound": False,
        "compound_name": None,
        "molecule_svg": None,
        "error": None,
    }

    try:
        # Resolve the R-group SMILES if a key was provided
        new_group_smiles = COMMON_R_GROUPS.get(new_group, new_group)

        # Apply the swap
        new_smiles = apply_swap(old_smiles, atom_idx, new_group_smiles)

        if new_smiles is None:
            result["error"] = (
                "Failed to apply swap. Invalid atom index or incompatible substitution."
            )
            return result

        result["new_smiles"] = new_smiles
        result["success"] = True

        # Validate the new compound
        mol = Chem.MolFromSmiles(new_smiles)
        if mol is None:
            result["error"] = "Resulting SMILES is invalid."
            result["success"] = False
            return result
        result["stability_score"] = calculate_stability_score(mol)
        result["drug_likeness"] = calculate_drug_likeness(mol)
        result["lipinski_compliance"] = check_lipinski_rules(mol)
        result["structural_alerts"] = check_structural_alerts(mol)

        # Check ChEMBL database
        try:
            is_real, compound_name = verify_if_real_compound(new_smiles)
            result["is_real_compound"] = is_real
            result["compound_name"] = compound_name
        except RuntimeError:
            # ChEMBL not available, but continue with other validations
            logger.warning("ChEMBL check skipped (not available)")

        # Generate visualization with highlighted modified atom
        try:
            # Get modifiable atoms for the new molecule
            modifiable_atoms = get_modifiable_atoms(new_smiles)

            # Try to find the modified atom in the new molecule
            highlight_idx = _find_modified_atom_in_new_molecule(
                new_smiles, new_group_smiles, atom_idx
            )
            result["molecule_svg"] = smiles_to_svg(
                new_smiles,
                highlight_atom_idx=highlight_idx,
                modifiable_atoms=modifiable_atoms,
            )
        except Exception as e:
            logger.warning(f"Could not generate SVG: {str(e)}")
            # Try without highlighting as fallback
            try:
                modifiable_atoms = get_modifiable_atoms(new_smiles)
                result["molecule_svg"] = smiles_to_svg(
                    new_smiles,
                    highlight_atom_idx=None,
                    modifiable_atoms=modifiable_atoms,
                )
            except:
                pass

        return result

    except Exception as e:
        logger.error(f"Error in validate_and_optimize_compound: {str(e)}")
        result["error"] = str(e)
        return result


# Backward compatibility aliases
get_synthesis_score = calculate_stability_score
