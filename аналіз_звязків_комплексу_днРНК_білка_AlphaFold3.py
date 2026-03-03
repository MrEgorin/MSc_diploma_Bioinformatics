from pymol import cmd
import os
import math

# --- Parameters ---
pdb_file = r"C:\docs\imbg_crispr\adar_anrassf1_complex\fold_2025_10_27_13_48\complex_A_A.pdb"
output_dir = r"C:\docs\imbg_crispr\adar_anrassf1_complex\fold_2025_10_27_13_48"
outfile = os.path.join(output_dir, "chain_contacts_full.txt")
cutoff = 4.5  # Å, distance cutoff for contacts

# --- Function to classify interaction type ---
def classify_interaction(atom1, atom2, distance):
    """
    Classify interaction type based on atom types and distance.
    Returns: Interaction type as a string.
    """
    # Get atom elements (e.g., C, N, O) from atom names
    elem1 = atom1.name[0] if atom1.name else "X"
    elem2 = atom2.name[0] if atom2.name else "X"
    
    # Hydrogen bond: N or O atoms within 2.5–3.5 Å
    if (elem1 in ['N', 'O'] and elem2 in ['N', 'O']) and 2.5 <= distance <= 3.5:
        return "Hydrogen Bond"
    # Hydrophobic: C-C contacts within 3.5–4.5 Å
    elif (elem1 == 'C' and elem2 == 'C') and 3.5 <= distance <= cutoff:
        return "Hydrophobic"
    # Other contacts within cutoff are van der Waals
    else:
        return "Van der Waals"

# --- Initialize PyMOL ---
cmd.reinitialize()

# --- Load PDB file ---
if not os.path.exists(pdb_file):
    print(f"❌ PDB file not found: {pdb_file}")
    exit(1)
cmd.load(pdb_file, "complex")

# --- Find available chains ---
chains = sorted(list(set(a.chain for a in cmd.get_model("all").atom if a.chain)))
if len(chains) < 2:
    print("❌ Found fewer than two chains in the PDB file!")
    exit(1)

# Assume first chain is protein, second is RNA (adjust if needed)
chain1, chain2 = chains[:2]
print(f"🔍 Analyzing contacts between chain {chain1} (protein) and chain {chain2} (RNA)...")

# --- Create selections ---
sel1 = f"chain {chain1}"
sel2 = f"chain {chain2}"

# --- Find all atomic pairs between chains within cutoff ---
cmd.select("sel1", sel1)
cmd.select("sel2", sel2)
pairs = cmd.find_pairs("sel1", "sel2", cutoff=cutoff)

# --- Write contacts to file ---
os.makedirs(output_dir, exist_ok=True)
with open(outfile, "w", encoding="utf-8") as f:
    # Write header
    f.write("chain1\tresn1\tresi1\tatom1\tchain2\tresn2\tresi2\tatom2\tdistance(Å)\tinteraction_type\n")
    
    # Process each pair
    contact_count = 0
    for (a1, a2) in pairs:
        # Get atom objects
        atom1 = cmd.get_model(a1).atom[0]
        atom2 = cmd.get_model(a2).atom[0]
        
        # Calculate distance
        distance = cmd.get_distance(a1, a2)
        
        # Classify interaction
        interaction_type = classify_interaction(atom1, atom2, distance)
        
        # Write to file
        f.write(f"{atom1.chain}\t{atom1.resn}\t{atom1.resi}\t{atom1.name}\t"
                f"{atom2.chain}\t{atom2.resn}\t{atom2.resi}\t{atom2.name}\t"
                f"{distance:.2f}\t{interaction_type}\n")
        contact_count += 1

# --- Clean up PyMOL selections ---
cmd.delete("sel1")
cmd.delete("sel2")

# --- Print summary ---
print(f"Saved {contact_count} contacts to file: {outfile}")
print(f"Output format: chain1, residue_name1, residue_id1, atom1, chain2, residue_name2, residue_id2, atom2, distance(Å), interaction_type")
