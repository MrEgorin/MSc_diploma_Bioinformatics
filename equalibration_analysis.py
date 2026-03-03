!pip install gromacs
!pip install -q mdtraj matplotlib numpy 2>/dev/null
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import mdtraj as md
import os

# -------------------------------------------------
def read_xvg(f):
    if not os.path.exists(f): return None, None
    data = np.loadtxt(f, comments=["#","@","&"])
    return data[:,0]/1000, data[:,1]      # ps → ns

def статистика(час, значення, останні_ns=50.0):
    if час is None or len(час)<100: return "немає даних"
    відсіч = час[-1] - останні_ns
    маска = час >= max(відсіч, час[-1]*0.7)
    return f"{np.mean(значення[маска]):.2f} ± {np.std(значення[маска]):.2f}"

# ------------------------------------------------- RMSD + Rg
top   = glob("*.gro") + glob("*.pdb") + glob("*.tpr")
traj  = glob("*.xtc") + glob("*.trr") + glob("*.dcd")

if top and traj:
    t = md.load(traj[0], top=top[0])
    sel = t.top.select("protein or resname LIG DRG UNK MOL")
    if len(sel)==0:
        sel = t.top.select("not water and not resname NA CL K ION SOL WAT")
    t = t.atom_slice(sel)

    rmsd = md.rmsd(t, t[0], 0) * 10          # Å
    rgyr = md.compute_rg(t) * 10             # Å
    час_traj = t.time / 1e6                  # ns

    rmsd_stat = статистика(час_traj, rmsd)
    rgyr_stat = статистика(час_traj, rgyr)

    print(f"RMSD (останні нс): {rmsd_stat} Å")
    print(f"Радіус гірації (останні нс): {rgyr_stat} Å")
else:
    rmsd = rgyr = час_traj = None
    rmsd_stat = rgyr_stat = "немає траєкторії"
    print("Увага: .gro/.pdb + .xtc/.trr не знайдено → RMSD/Rg пропущено")

# ------------------------------------------------- Малюємо
plt.style.use('seaborn-v0_8-darkgrid')
fig = plt.figure(figsize=(22, 14))
gs  = fig.add_gridspec(3, 4, hspace=0.32, wspace=0.38)

всі_файли = {
    "Температура (K)"       : glob("*temperature*.xvg"),
    "Тиск (бар)"            : glob("*pressure*.xvg"),
    "Густина (кг/м³)"       : glob("*density*.xvg"),
    "Потенційна енергія"    : glob("*potential*.xvg"),
    "Об'єм боксу (нм³)"     : glob("*volume*.xvg")
}

for idx, (назва, файли) in enumerate(всі_файли.items()):
    ax = fig.add_subplot(gs[idx//3, idx%3])
    for f in файли:
        t, y = read_xvg(f)
        if t is None: continue
        ax.plot(t, y, lw=1.7, label=os.path.basename(f).replace(".xvg",""))
        стат = статистика(t, y)
        ax.text(0.97, 0.93, f": {стат}", transform=ax.transAxes,
                ha='right', va='top', fontsize=11,
                bbox=dict(facecolor='white', alpha=0.8))
    ax.set_title(назва, fontsize=14)
    ax.set_xlabel("Час (нс)")
    if "Температура" in назва: ax.axhline(300, color='red', ls='--', alpha=0.7)
    if "Тиск"        in назва: ax.axhline(1,   color='red', ls='--', alpha=0.7)
    if "Густина"     in назва: ax.axhline(1000,color='red', ls='--', alpha=0.7)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.4)

# RMSD і Rg внизу
if rmsd is not None:
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(час_traj, rmsd, color='#1f77b4', lw=2.2)
    ax.set_title(f"RMSD (Å)\n{rmsd_stat}")
    ax.set_xlabel("Час (нс)")
    ax.grid(alpha=0.4)

    ax = fig.add_subplot(gs[2, 1])
    ax.plot(час_traj, rgyr, color='#ff7f0e', lw=2.2)
    ax.set_title(f"Радіус гірації (Å)\n{rgyr_stat}")
    ax.set_xlabel("Час (нс)")
    ax.grid(alpha=0.4)

plt.suptitle("ПОВНА ПЕРЕВІРКА РІВНОВАГИ МОЛЕКУЛЯРНОЇ ДИНАМІКИ", fontsize=19, y=0.98)
plt.show()
