#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Простий аналіз взаємодій РНК-білок через PyMOL
Показує: які нуклеотиди з якими амінокислотами взаємодіють
"""

from pymol import cmd
import os

# ============================================================
# НАЛАШТУВАННЯ
# ============================================================
WORK_DIR = r"C:\docs\imbg_crispr\adar_anrassf1_complex\mprdock"
PDB_FILE = os.path.join(WORK_DIR, "dock_top1.pdb")
OUTPUT_FILE = os.path.join(WORK_DIR, "interactions_simple.txt")

# Відстані для різних типів взаємодій
DIST_HBOND = 3.5      # Водневі зв'язки
DIST_SALT = 4.0       # Солеві містки
DIST_HYDROPHOBIC = 4.5  # Гідрофобні контакти
DIST_GENERAL = 4.0    # Загальні контакти

# ============================================================
# ЗАВАНТАЖЕННЯ СТРУКТУРИ
# ============================================================
print(f"Завантажуємо: {PDB_FILE}")
cmd.load(PDB_FILE, "complex")

# ============================================================
# ВИЗНАЧЕННЯ СЕЛЕКЦІЙ
# ============================================================
# Білок
cmd.select("protein", "polymer.protein")

# РНК (різні формати: A/C/G/U або RA/RC/RG/RU)
cmd.select("rna", "resn A+C+G+U+RA+RC+RG+RU+DA+DC+DG+DT")

# Перевірка
n_prot = cmd.count_atoms("protein")
n_rna = cmd.count_atoms("rna")
print(f"Білок: {n_prot} атомів")
print(f"РНК: {n_rna} атомів")

if n_prot == 0 or n_rna == 0:
    print("ПОМИЛКА: Не знайдено білок або РНК!")
    exit(1)

# ============================================================
# ФУНКЦІЯ: Знаходження взаємодій
# ============================================================
def find_interactions(sel1, sel2, distance, int_type):
    """
    Знаходить атоми з sel1, що близькі до sel2
    Повертає список: [(залишок1, атом1, залишок2, атом2, відстань)]
    """
    interactions = []
    
    # Отримуємо всі атоми з першої селекції
    atoms1 = []
    cmd.iterate(sel1, "atoms1.append((resi, resn, name, chain))", space={'atoms1': atoms1})
    
    # Для кожного атома з sel1 шукаємо близькі в sel2
    for i, (resi1, resn1, name1, chain1) in enumerate(atoms1):
        # Створюємо тимчасову селекцію для цього атома
        temp_sel = f"temp_atom_{i}"
        cmd.select(temp_sel, f"{sel1} and resi {resi1} and chain {chain1} and name {name1}")
        
        # Шукаємо близькі атоми в sel2
        nearby_sel = f"nearby_{i}"
        cmd.select(nearby_sel, f"{sel2} within {distance} of {temp_sel}")
        
        # Отримуємо інформацію про близькі атоми
        atoms2 = []
        cmd.iterate(nearby_sel, "atoms2.append((resi, resn, name, chain))", space={'atoms2': atoms2})
        
        # Для кожної пари обчислюємо точну відстань
        for resi2, resn2, name2, chain2 in atoms2:
            # Створюємо селекції для обох атомів
            atom1_sel = f"{sel1} and resi {resi1} and chain {chain1} and name {name1}"
            atom2_sel = f"{sel2} and resi {resi2} and chain {chain2} and name {name2}"
            
            # Обчислюємо відстань
            dist = cmd.distance("temp_dist", atom1_sel, atom2_sel)
            cmd.delete("temp_dist")
            
            if dist <= distance:
                interactions.append({
                    'prot_res': f"{resn1}{resi1}",
                    'prot_atom': name1,
                    'rna_res': f"{resn2}{resi2}",
                    'rna_atom': name2,
                    'distance': round(dist, 2),
                    'type': int_type
                })
        
        # Видаляємо тимчасові селекції
        cmd.delete(temp_sel)
        cmd.delete(nearby_sel)
    
    return interactions

# ============================================================
# ПОШУК РІЗНИХ ТИПІВ ВЗАЄМОДІЙ
# ============================================================
all_interactions = []

print("\n" + "="*60)
print("АНАЛІЗ ВЗАЄМОДІЙ")
print("="*60)

# ------------------------------------------------------------
# 1. ВОДНЕВІ ЗВ'ЯЗКИ (O, N атоми)
# ------------------------------------------------------------
print("\n1. Шукаємо водневі зв'язки...")
cmd.select("donors_prot", "protein and (elem N or elem O)")
cmd.select("acceptors_rna", "rna and (elem N or elem O)")

hbonds = find_interactions("donors_prot", "acceptors_rna", DIST_HBOND, "H-bond")
all_interactions.extend(hbonds)
print(f"   Знайдено: {len(hbonds)} водневих зв'язків")

# ------------------------------------------------------------
# 2. СОЛЕНІ МІСТКИ (позитивні з негативними)
# ------------------------------------------------------------
print("\n2. Шукаємо солені містки...")
# Позитивні: LYS (NZ), ARG (NH1, NH2), HIS (ND1, NE2)
cmd.select("positive", "protein and (resn LYS and name NZ) or (resn ARG and name NH*) or (resn HIS and name ND1+NE2)")
# Негативні: фосфати РНК
cmd.select("negative", "rna and name OP1+OP2+O5'+O3'")

salt = find_interactions("positive", "negative", DIST_SALT, "Salt-bridge")
all_interactions.extend(salt)
print(f"   Знайдено: {len(salt)} солених містків")

# ------------------------------------------------------------
# 3. ГІДРОФОБНІ КОНТАКТИ (вуглеці)
# ------------------------------------------------------------
print("\n3. Шукаємо гідрофобні контакти...")
# Гідрофобні залишки білка
cmd.select("hydrophobic_prot", "protein and resn ALA+VAL+LEU+ILE+PHE+TRP+MET+PRO and elem C")
# Вуглеці основ РНК
cmd.select("hydrophobic_rna", "rna and (name C2+C4+C5+C6+C8)")

hydro = find_interactions("hydrophobic_prot", "hydrophobic_rna", DIST_HYDROPHOBIC, "Hydrophobic")
all_interactions.extend(hydro)
print(f"   Знайдено: {len(hydro)} гідрофобних контактів")

# ------------------------------------------------------------
# 4. ПІ-СТЕКІНГ (ароматичні кільця)
# ------------------------------------------------------------
print("\n4. Шукаємо π-стекінг...")
# Ароматичні залишки
cmd.select("aromatic_prot", "protein and resn PHE+TYR+TRP+HIS and elem C")
cmd.select("aromatic_rna", "rna and (name C2+C4+C5+C6+C8)")

pi_stack = find_interactions("aromatic_prot", "aromatic_rna", 4.5, "Pi-stacking")
all_interactions.extend(pi_stack)
print(f"   Знайдено: {len(pi_stack)} π-стекінгів")

# ============================================================
# ГРУПУВАННЯ ПО ПАРАХ ЗАЛИШКІВ
# ============================================================
from collections import defaultdict

grouped = defaultdict(list)
for inter in all_interactions:
    key = (inter['prot_res'], inter['rna_res'])
    grouped[key].append(inter)

# Сортуємо по найкоротшій відстані
sorted_pairs = sorted(grouped.items(), 
                     key=lambda x: min(i['distance'] for i in x[1]))

# ============================================================
# ЗБЕРЕЖЕННЯ РЕЗУЛЬТАТІВ
# ============================================================
print(f"\n" + "="*60)
print(f"ЗБЕРЕЖЕННЯ РЕЗУЛЬТАТІВ")
print("="*60)

with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
    f.write("Амінокислота\tАтом_білка\tНуклеотид\tАтом_РНК\tТип_взаємодії\tВідстань_Å\n")
    f.write("-" * 80 + "\n")
    
    for (prot_res, rna_res), interactions in sorted_pairs:
        for inter in interactions:
            f.write(f"{inter['prot_res']}\t"
                   f"{inter['prot_atom']}\t"
                   f"{inter['rna_res']}\t"
                   f"{inter['rna_atom']}\t"
                   f"{inter['type']}\t"
                   f"{inter['distance']}\n")

print(f"✓ Результати збережено: {OUTPUT_FILE}")
print(f"✓ Всього взаємодій: {len(all_interactions)}")
print(f"✓ Унікальних пар: {len(grouped)}")

# ============================================================
# ВИВЕДЕННЯ ТОП-20 НА ЕКРАН
# ============================================================
print(f"\n" + "="*60)
print("ТОП-20 ВЗАЄМОДІЙ (найкоротші відстані)")
print("="*60)

for i, ((prot_res, rna_res), interactions) in enumerate(sorted_pairs[:20], 1):
    min_dist = min(inter['distance'] for inter in interactions)
    types = set(inter['type'] for inter in interactions)
    
    print(f"\n{i}. {prot_res} ↔ {rna_res} [{min_dist:.2f} Å]")
    print(f"   Типи: {', '.join(types)}")
    
    for inter in interactions:
        print(f"   • {inter['prot_atom']} → {inter['rna_atom']} "
              f"({inter['type']}, {inter['distance']} Å)")

# ============================================================
# СТАТИСТИКА
# ============================================================
print(f"\n" + "="*60)
print("СТАТИСТИКА")
print("="*60)

type_counts = defaultdict(int)
for inter in all_interactions:
    type_counts[inter['type']] += 1

for int_type, count in sorted(type_counts.items()):
    print(f"{int_type:15s}: {count:4d}")

print(f"\n{'ЗАГАЛОМ':15s}: {len(all_interactions):4d}")

# Найактивніші залишки
prot_residues = defaultdict(int)
rna_residues = defaultdict(int)

for inter in all_interactions:
    prot_residues[inter['prot_res']] += 1
    rna_residues[inter['rna_res']] += 1

print(f"\n" + "="*60)
print("ТОП-10 НАЙАКТИВНІШИХ АМІНОКИСЛОТ")
print("="*60)
for res, count in sorted(prot_residues.items(), key=lambda x: x[1], reverse=True)[:10]:
    print(f"{res:10s}: {count:3d} взаємодій")

print(f"\n" + "="*60)
print("ТОП-10 НАЙАКТИВНІШИХ НУКЛЕОТИДІВ")
print("="*60)
for res, count in sorted(rna_residues.items(), key=lambda x: x[1], reverse=True)[:10]:
    print(f"{res:10s}: {count:3d} взаємодій")

print(f"\n✓ Аналіз завершено!")
print(f"✓ Детальні результати у файлі: {OUTPUT_FILE}")

# Очищаємо
cmd.delete("all")


