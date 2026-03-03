import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(r"C:\docs\imbg_crispr\docking\md\productive_md")
print(f"Робоча папка: {os.getcwd()}")

# Перевіряємо, чи є файли
required_files = [f + ".xvg" for f in ["temperature","potential","pressure","rmsd","rmsf","rgyr","hbond","sasa"]]
missing = [f for f in required_files if not os.path.exists(f)]
if missing:
    print("Не знайдено файли:", missing)
    input("Натисни Enter, щоб закрити...")
    exit()

files = ["temperature","potential","pressure","rmsd","rmsf","rgyr","hbond","sasa"]

input("Натисни Enter для української версії...")

titles_ua = {
    "temperature": "Температура",
    "potential":   "Потенційна енергія",
    "pressure":    "Тиск",
    "rmsd":        "RMSD остова РНК",
    "rmsf":        "RMSF по нуклеотидах",
    "rgyr":        "Радіус гірації",
    "hbond":       "Кількість водневих зв’язків РНК–ліганд",
    "sasa":        "Доступна поверхня розчинника (SASA)"
}
units_ua = {
    "temperature": "К", "potential": "кДж/моль", "pressure": "бар",
    "rmsd": "нм", "rmsf": "нм", "rgyr": "нм", "hbond": "шт", "sasa": "нм²"
}

for f in files:
    data = np.loadtxt(f + ".xvg", comments=["@", "#"])
    t, y = data[:, 0], data[:, 1]

    plt.figure(figsize=(11, 6))
    plt.plot(t, y, color="#1f77b4", linewidth=2.5)
    plt.title(titles_ua[f], fontsize=20, pad=20)
    plt.xlabel("Час (нс)" if f in ["rmsd"] else "Час (пс)", fontsize=16)
    plt.ylabel(units_ua[f], fontsize=16)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(4)
    plt.close()


