# Transport de charges 1D instationnaire dans un matériau diélectrique

Simulation 1D (géométrie plane, anode à gauche / cathode à droite) du
transport bipolaire de charges (électrons et trous, libres et piégés)
dans un diélectrique d'épaisseur L, par :

- **volumes finis** sur un **maillage non uniforme** raffiné aux électrodes,
- **schéma MUSCL** (pente limitée minmod) pour l'advection, flux centré
  pour la diffusion,
- **Poisson** résolue à chaque pas par l'**algorithme de Thomas**,
- **RK2 (Heun)** en temps avec **pas de temps adaptatif**,
- injections d'électrodes **Richardson–Schottky** (dépendantes du signe
  du champ E), piégeage / dépiégeage / recombinaison.

## Fichiers

| Fichier | Rôle |
|---|---|
| `transport_charges.py` | code de simulation complet (commenté, sous-programmes) |
| `generer_pdf_explication.py` | génère `explication.pdf` (documentation complète du modèle et du schéma) |
| `explication.pdf` | documentation : modèle, équations, conditions aux limites, méthode numérique, utilisation |
| `resultats/` | sorties : `resultats.npz` + figures PNG |

## Utilisation

```bash
pip install numpy matplotlib
python3 transport_charges.py          # simulation par défaut (Va = 1000 V, 50 s)
python3 generer_pdf_explication.py    # (re)génère la documentation PDF
```

Personnalisation (potentiel d'anode constant / rampe / sinus, durée,
paramètres physiques et numériques) : tout est dans la classe
`Parametres` de `transport_charges.py`.

```python
import transport_charges as tc
p = tc.Parametres()
p.va_type = 'sinus'          # ou 'constante' / 'rampe'
p.va0, p.va_amp, p.va_freq = 500., 300., 0.05
p.tmax = 100.
res = tc.simuler(p)
tc.tracer(res, p)
```

## Sorties

Toutes les grandeurs en fonction du temps **et** de l'espace :
densités de charges libres et piégées, potentiel V(x,t), champ E(x,t),
charges totales intégrées, potentiel d'anode, courants aux électrodes
(`resultats.npz` + `profils_charges.png`, `potentiel_champ.png`,
`evolutions_temporelles.png`, `maillage.png`).
