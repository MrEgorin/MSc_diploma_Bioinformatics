from Bio.PDB import PDBParser, NeighborSearch
import sys

# Шлях до PDB-файлу 
pdb_file = r"C:\docs\4sem\diplom\cluster1_1.pdb"

# Парсер PDB
parser = PDBParser(QUIET=True)
structure = parser.get_structure('complex', pdb_file)

# Функція для визначення, чи атом з білка (амінокислота) чи РНК (нуклеотид)
def is_protein_residue(res):
    return res.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def is_rna_residue(res):
    return res.get_resname() in ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']  # Для РНК/ДНК

# Збираємо всі атоми
all_atoms = [atom for model in structure for chain in model for res in chain for atom in res]

# NeighborSearch для швидкого пошуку сусідів (cutoff 3.5 Å для H-bonds)
ns = NeighborSearch(all_atoms)

# Потенційні донори/акцептори для H-bonds
protein_hbond_atoms = ['N', 'O', 'S']  # Атоми в білку (N-H, O-H тощо)
rna_hbond_atoms = ['O', 'N', 'P']      # Атоми в РНК (O в цукрі/фосфаті, N/O в базі)

# Збираємо взаємодії
interactions = []

for model in structure:
    for chain_prot in model:
        for res_prot in chain_prot:
            if not is_protein_residue(res_prot):
                continue
            prot_id = res_prot.id[1]  # Номер амінокислоти
            prot_name = res_prot.resname
            for atom_prot in res_prot:
                if atom_prot.element not in protein_hbond_atoms:
                    continue
                # Шукаємо сусідів в РНК
                neighbors = ns.search(atom_prot.coord, 3.5, 'A')  # Атоми в радіусі 3.5 Å
                for atom_rna in neighbors:
                    res_rna = atom_rna.get_parent()
                    if not is_rna_residue(res_rna):
                        continue
                    if atom_rna.element not in rna_hbond_atoms:
                        continue
                    rna_id = res_rna.id[1]  # Номер нуклеотиду
                    rna_name = res_rna.resname
                    distance = atom_prot - atom_rna  # Відстань
                    interactions.append({
                        'prot_res': prot_name + str(prot_id),
                        'prot_atom': atom_prot.get_name(),
                        'rna_res': rna_name + str(rna_id),
                        'rna_atom': atom_rna.get_name(),
                        'type': 'H-bond',
                        'distance': round(distance, 2)
                    })

# Сортуємо за відстанню або номером (тут за номером білка)
interactions.sort(key=lambda x: (int(x['prot_res'][3:]), x['distance']))

# Вивід у табличному форматі
print(f"{'Амінокислота':<15} {'Атом білка':<10} {'Нуклеотида':<10} {'Атом РНК':<10} {'Тип взаємодії':<15} {'Відстань Å':<10}")
print('-' * 80)
for inter in interactions:
    print(f"{inter['prot_res']:<15} {inter['prot_atom']:<10} {inter['rna_res']:<10} {inter['rna_atom']:<10} {inter['type']:<15} {inter['distance']:<10}")

# Якщо взаємодій немає, вивід помилки
if not interactions:
    print("Не знайдено взаємодій. Перевірте PDB-файл, ланцюги або cutoff.")

# Збереження результату в файл у потрібній папці
output_file = r"C:\docs\4sem\diplom\result_complex_11_hd.txt"

with open(output_file, "w", encoding="utf-8") as f:
    # Заголовок
    f.write(f"{'Амінокислота':<15} {'Атом білка':<10} {'Нуклеотид':<10} {'Атом РНК':<10} {'Тип взаємодії':<15} {'Відстань Å':<10}\n")
    f.write("-" * 80 + "\n")
    
    # Всі взаємодії
    for inter in interactions:
        line = f"{inter['prot_res']:<15} {inter['prot_atom']:<10} {inter['rna_res']:<10} {inter['rna_atom']:<10} {inter['type']:<15} {inter['distance']:<10}\n"
        f.write(line)

print(f"\nРезультати збережено у файл: {output_file}")

