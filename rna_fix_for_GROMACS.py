#!/usr/bin/env python3
"""
Автоматичне виправлення відсутнього O5' у першому нуклеотиді (A 1)
Вставляє O5' після C5' з фіктивними координатами
Перенумеровує всі наступні атоми
Зберігає: rna_fixed.pdb
"""

import os
import sys
from pathlib import Path

# === Налаштування ===
INPUT_PDB = Path("C:/docs/imbg_crispr/docking/fixrna/rna.pdb")
OUTPUT_PDB = Path("C:/docs/imbg_crispr/docking/fixrna/rna_fixed.pdb")

# Перевірка існування файлу
if not INPUT_PDB.exists():
    print(f"Помилка: файл не знайдено: {INPUT_PDB}")
    sys.exit(1)

# Читання PDB
with open(INPUT_PDB, 'r') as f:
    lines = f.readlines()

# Шукаємо C5' у першому залишку (A 1)
c5p_line = None
c5p_index = None

for i, line in enumerate(lines):
    if line.startswith("ATOM") or line.startswith("HETATM"):
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21]
        res_num = int(line[22:26].strip())
        
        if atom_name == "C5'" and res_name == "A" and res_num == 1:
            c5p_line = line
            c5p_index = i
            break

if c5p_line is None:
    print("Помилка: не знайдено C5' у залишку A 1")
    sys.exit(1)

print(f"Знайдено C5' на рядку {c5p_index + 1}")

# Парсимо координати C5'
x = float(c5p_line[30:38])
y = float(c5p_line[38:46])
z = float(c5p_line[46:54])

# Фіктивні координати для O5': на 0.1 нм ліворуч від C5' (вздовж -X)
o5p_x = x - 1.0  # 1.0 Å = 0.1 нм
o5p_y = y
o5p_z = z

# Формуємо новий рядок O5'
new_o5p_line = (
    f"ATOM  {c5p_index + 1:>5d}  O5' A A   1"
    f"{o5p_x:12.3f}{o5p_y:8.3f}{o5p_z:8.3f}  1.00  0.00           O  \n"
)

# Вставляємо O5' після C5'
lines.insert(c5p_index + 1, new_o5p_line)

# Перенумеровуємо всі наступні атоми
for i in range(c5p_index + 2, len(lines)):
    line = lines[i]
    if line.startswith("ATOM") or line.startswith("HETATM"):
        # Отримуємо старий номер
        old_num = int(line[6:11])
        # Новий номер
        new_num = old_num + 1
        # Замінюємо
        new_line = line[:6] + f"{new_num:5d}" + line[11:]
        lines[i] = new_line

# Зберігаємо
with open(OUTPUT_PDB, 'w') as f:
    f.writelines(lines)


