import os
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Superimposer, Select
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Функція для завантаження структури з PDB файлу
def load_structure(pdb_file):
    if not os.path.exists(pdb_file):
        print(f"Файл {pdb_file} не існує!")
        return None
    parser = PDBParser(QUIET=True)
    return parser.get_structure(os.path.basename(pdb_file).split('.')[0], pdb_file)

# Функція для вибору атомів backbone (для РНК - фосфор P)
def get_backbone_atoms(struct, chain_id='A'):
    backbone_atoms = []
    for model in struct:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if 'P' in residue:
                    backbone_atoms.append(residue['P'])
    return backbone_atoms

# Функція для обрізання атомів до мінімальної довжини
def trim_atoms(ref_atoms, mob_atoms):
    min_len = min(len(ref_atoms), len(mob_atoms))
    return ref_atoms[:min_len], mob_atoms[:min_len]

# Функція для обчислення RMSD за допомогою Superimposer
def compute_rmsd_and_align(ref_struct, mobile_struct, ref_atoms, mobile_atoms):
    ref_atoms, mob_atoms = trim_atoms(ref_atoms, mobile_atoms)
    if len(ref_atoms) != len(mob_atoms):
        print(f"Кількість атомів не співпадає після обрізання: {len(ref_atoms)} vs {len(mob_atoms)}")
        return None
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)
    sup.apply(mobile_struct.get_atoms())
    return sup.rms

# Функція для попарного вирівнювання (pairwise)
def pairwise_alignment(files):
    print("=== Попарне структурне вирівнювання ===")
    rmsd_pairs = {}
    for i, f1 in enumerate(files):
        for j, f2 in enumerate(files[i+1:]):
            ref = load_structure(f1)
            mob = load_structure(f2)
            if ref is None or mob is None:
                continue
            ref_atoms = get_backbone_atoms(ref)
            mob_atoms = get_backbone_atoms(mob)
            rmsd = compute_rmsd_and_align(ref, mob, ref_atoms, mob_atoms)
            if rmsd is not None:
                key = (os.path.basename(f1), os.path.basename(f2))
                rmsd_pairs[key] = rmsd
                print(f"RMSD між {key[0]} і {key[1]}: {rmsd:.2f} Å")
    return rmsd_pairs

# Функція для all-to-one вирівнювання відносно reference
def all_to_one_alignment(files, ref_file, base_output_dir="align", threshold_bad=5.0):
    ref_struct = load_structure(ref_file)
    if ref_struct is None:
        return
    ref_atoms = get_backbone_atoms(ref_struct)
    
    aligned_structs = {}
    rmsd_values = {}
    
    print(f"\n=== Вирівнювання всіх до {os.path.basename(ref_file)} ===")
    for f in files:
        if f == ref_file:
            continue
        mob_struct = load_structure(f)
        if mob_struct is None:
            continue
        mob_atoms = get_backbone_atoms(mob_struct)
        rmsd = compute_rmsd_and_align(ref_struct, mob_struct, ref_atoms, mob_atoms)
        if rmsd is not None:
            base_name = os.path.basename(f).split('.')[0]
            suffix = "_bad" if rmsd > threshold_bad else ""
            output_dir = os.path.join(os.path.dirname(ref_file), base_output_dir, f"aligned_to_{os.path.basename(ref_file).split('.')[0]}")
            os.makedirs(output_dir, exist_ok=True)
            output_pdb = os.path.join(output_dir, f"{base_name}{suffix}_aligned.pdb")
            
            # Збереження вирівняної структури
            io = PDBIO()
            class AllSelect(Select):
                def accept_residue(self, residue):
                    return 1  # Зберегти всі
            io.set_structure(mob_struct)
            io.save(output_pdb, AllSelect())
            
            aligned_structs[f] = mob_struct
            rmsd_values[f] = rmsd
            print(f"RMSD {os.path.basename(f)} до {os.path.basename(ref_file)}: {rmsd:.2f} Å (збережено як {output_pdb})")
    
    # Статистика
    if rmsd_values:
        mean_rmsd = np.mean(list(rmsd_values.values()))
        std_rmsd = np.std(list(rmsd_values.values()))
        print(f"Середній RMSD: {mean_rmsd:.2f} ± {std_rmsd:.2f} Å")
        max_rmsd = max(rmsd_values.values())
        print(f"Максимальний RMSD: {max_rmsd:.2f} Å")
    
    return aligned_structs, rmsd_values

# Функція для обчислення пер-резіду RMSD (для консервативних регіонів)
def per_residue_rmsd(ref_chain, aligned_chains):
    n_res = len(list(ref_chain))
    per_res_rmsd = np.zeros(n_res)
    n_structs = len(aligned_chains)
    
    for i, res in enumerate(ref_chain):
        if 'P' not in res:
            continue
        ref_coord = res['P'].get_coord()
        res_rmsd = []
        for chain in aligned_chains.values():
            if i < len(list(chain)):
                mob_res = list(chain)[i]
                if 'P' in mob_res:
                    mob_coord = mob_res['P'].get_coord()
                    dist = np.linalg.norm(ref_coord - mob_coord)
                    res_rmsd.append(dist)
        if res_rmsd:
            per_res_rmsd[i] = np.mean(res_rmsd)
    
    return per_res_rmsd

# Функція для знаходження консервативних регіонів (низький RMSD, e.g. < 2Å)
def find_conservative_regions(ref_struct, aligned_structs, rmsd_per_res, threshold=2.0, min_length=5):
    ref_chain = list(ref_struct.get_chains())[0]
    aligned_chains = {f: list(s.get_chains())[0] for f, s in aligned_structs.items()}
    
    conservative_regions = []
    in_region = False
    start = 0
    
    for i, rmsd in enumerate(rmsd_per_res):
        if rmsd < threshold:
            if not in_region:
                start = i
                in_region = True
        else:
            if in_region and (i - start) >= min_length:
                conservative_regions.append((start, i-1))
            in_region = False
    
    if in_region and (len(list(ref_chain)) - start) >= min_length:
        conservative_regions.append((start, len(list(ref_chain))-1))
    
    return conservative_regions, aligned_chains

# Функція для збереження консервативного регіону в PDB і FASTA
def save_conservative_region(ref_struct, aligned_structs, region, reg_id, base_output_dir="align"):
    start, end = region
    ref_chain = list(ref_struct.get_chains())[0]
    
    # PDB для reference
    io = PDBIO()
    class RegionSelect(Select):
        def accept_residue(self, residue):
            res_id = residue.get_id()[1]
            return 1 if start <= res_id - 1 <= end else 0
    ref_copy = ref_struct.copy()
    output_dir = os.path.join(os.path.dirname(ref_struct.id + '.pdb'), base_output_dir, "conservative")
    os.makedirs(output_dir, exist_ok=True)
    io.set_structure(ref_copy)
    io.save(os.path.join(output_dir, f"{reg_id}_ref.pdb"), RegionSelect())
    
    # FASTA для всіх
    records = []
    structs_dict = {os.path.basename(ref_struct.id + '.pdb'): ref_struct}
    structs_dict.update({os.path.basename(k): v for k, v in aligned_structs.items()})
    for f_name, s in structs_dict.items():
        chain = list(s.get_chains())[0]
        seq = ''.join([res.get_resname().strip() for res in chain if start <= list(chain).index(res) <= end and res.get_resname().strip() in 'AUGC'])
        records.append(SeqRecord(Seq(seq), id=f"{f_name.split('.')[0]}_{reg_id}"))
    
    fasta_file = os.path.join(output_dir, f"{reg_id}.fasta")
    SeqIO.write(records, fasta_file, "fasta")
    print(f"Збережено регіон {reg_id} (позиції {start+1}-{end+1}) в {fasta_file} і PDB")

# Головна функція
def main():
    pdb_dir = r"C:\docs\imbg_crispr\output_molprobity"
    files = [
        r"C:\docs\imbg_crispr\output_molprobity\alphafoldFH.pdb",
        r"C:\docs\imbg_crispr\output_molprobity\farfar2_S_000001_059_1FH.pdb",
        r"C:\docs\imbg_crispr\output_molprobity\nufold_5lnc_rank_1FH.pdb",
        r"C:\docs\imbg_crispr\output_molprobity\rosetta_fold_model1_cleanFH.pdb"
    ]
    print(f"Файли для обробки: {[os.path.basename(f) for f in files]}")
    
    # 1. Попарне вирівнювання
    pairwise_rmsd = pairwise_alignment(files)
    
    # 2. All-to-one до FARFAR2
    farfar_ref = r"C:\docs\imbg_crispr\output_molprobity\farfar2_S_000001_059_1FH.pdb"
    aligned_to_farfar, rmsd_to_farfar = all_to_one_alignment(files, farfar_ref)
    
    # 3. All-to-one до AlphaFold
    alphafold_ref = r"C:\docs\imbg_crispr\output_molprobity\alphafoldFH.pdb"
    if os.path.exists(alphafold_ref):
        aligned_to_alpha, rmsd_to_alpha = all_to_one_alignment(files, alphafold_ref)
        final_aligned = aligned_to_alpha
    else:
        print("AlphaFold PDB не знайдено, використовуємо FARFAR2")
        final_aligned = aligned_to_farfar
    
    # 4. Консервативні регіони
    ref_struct = load_structure(farfar_ref)
    if ref_struct and final_aligned:
        ref_chain = list(ref_struct.get_chains())[0]
        aligned_chains = {f: list(s.get_chains())[0] for f, s in final_aligned.items()}
        rmsd_per_res = per_residue_rmsd(ref_chain, aligned_chains)
        regions, _ = find_conservative_regions(ref_struct, final_aligned, rmsd_per_res)
        for i, reg in enumerate(regions):
            save_conservative_region(ref_struct, final_aligned, reg, f"region_{i+1}")
        print(f"Знайдено {len(regions)} консервативних регіонів")
    else:
        print("Не вдалося обчислити консервативні регіони")

if __name__ == "__main__":
    main()
