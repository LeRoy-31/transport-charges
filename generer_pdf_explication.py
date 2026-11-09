# -*- coding: utf-8 -*-
"""
Génère le document PDF "explication.pdf" qui explique en détail :
  - le modèle physique (équations de Poisson, transport, piégeage),
  - les conditions aux limites et les injections d'électrodes,
  - la méthode numérique (volumes finis, schéma MUSCL, RK2, maillage
    non uniforme, pas de temps adaptatif),
  - l'utilisation du code et le contenu des sorties.

Le PDF est composé de pages A4 générées avec matplotlib (mathtext pour
les équations). Lancer simplement :  python3 generer_pdf_explication.py
"""

import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TITRE = "Transport de charges dans un matériau diélectrique"
SOUS_TITRE = ("Modèle 1D instationnaire bipolaire (électrons / trous, "
              "libres / piégés)\nVolumes finis — schéma MUSCL — "
              "documentation et explication")


def nouvelle_page():
    """Crée une page A4 blanche."""
    fig = plt.figure(figsize=(8.27, 11.69))
    return fig


def bloc_texte(fig, x, y, lignes, taille=10.5, interligne=1.45):
    """Écrit une liste de (texte, style) sur la page.

    'style' : 'titre' | 'texte' | 'eq' | 'puce' | 'code'
    Les textes peuvent contenir des sauts de ligne : chaque ligne
    physique occupe une ligne verticale (évite les chevauchements).
    """
    for texte, style in lignes:
        lignes_txt = texte.split("\n")
        if style == "titre":
            y -= 0.022
            fig.text(x, y, lignes_txt[0], fontsize=13, fontweight="bold",
                     color="#1a3d6d")
            y -= 0.006
        elif style == "eq":
            for ln in lignes_txt:
                y -= 0.030
                fig.text(x + 0.03, y, ln, fontsize=12, va="top")
            y -= 0.006
        elif style == "code":
            for ln in lignes_txt:
                y -= 0.014
                fig.text(x + 0.03, y, ln if ln else " ", fontsize=9.5,
                         va="top", family="monospace")
        elif style == "puce":
            # retour à la ligne manuel (largeur ~92 caractères)
            lignes_txt = textwrap.fill(texte, 92).split("\n")
            for j, ln in enumerate(lignes_txt):
                y -= interligne * 0.0125
                if j == 0:
                    fig.text(x, y, u"\u2022", fontsize=10.5)
                fig.text(x + 0.012, y, ln, fontsize=taille, va="top")
        else:
            lignes_txt = textwrap.fill(texte, 100).split("\n")
            for ln in lignes_txt:
                y -= interligne * 0.0125
                fig.text(x, y, ln, fontsize=taille, va="top")
    return y


def page_titre(pdf):
    fig = nouvelle_page()
    fig.text(0.5, 0.62, TITRE, ha="center", fontsize=22,
             fontweight="bold", color="#1a3d6d")
    fig.text(0.5, 0.54, SOUS_TITRE, ha="center", fontsize=13)
    fig.text(0.5, 0.40,
             "Contenu :\n"
             u"  1. Modèle physique et hypothèses\n"
             u"  2. Équations du modèle\n"
             u"  3. Conditions aux limites et injections\n"
             u"  4. Méthode numérique (volumes finis, MUSCL, RK2)\n"
             u"  5. Maillage non uniforme et pas de temps adaptatif\n"
             u"  6. Utilisation du code et sorties\n"
             u"  7. Structure du programme",
             ha="center", fontsize=12, va="top")
    fig.text(0.5, 0.06, "Généré automatiquement — projet transport de charges",
             ha="center", fontsize=8, color="gray")
    pdf.savefig(fig)
    plt.close(fig)


def page_modele(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96, "1. Modèle physique et hypothèses",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        ("Le matériau diélectrique (épaisseur L = 200 µm) est placé entre "
         "deux électrodes planes :", "texte"),
        (u"l'ANODE à gauche (x = 0), portée au potentiel Va(t) — constant, "
         u"rampe ou sinusoïde —", "puce"),
        (u"la CATHODE à droite (x = L), au potentiel nul.", "puce"),
        ("", "texte"),
        ("Quatre populations de charges sont suivies, exprimées en densité "
         "de charge [C/m³] :", "texte"),
        (u"n_eµ : électrons LIBRES (négatifs)      n_et : électrons PIÉGÉS",
         "puce"),
        (u"n_hµ : trous LIBRES (positifs)          n_ht : trous PIÉGÉS",
         "puce"),
        ("", "texte"),
        ("À t = 0 le matériau est vierge : aucune charge (libre ou piégée) "
         "n'est présente.", "texte"),
        ("", "texte"),
        ("Le modèle couple :", "texte"),
        (u"l'équation de POISSON, qui relie le potentiel et le champ "
         u"électrique à la densité de charge nette ;", "puce"),
        (u"des équations d'ADVECTION / DIFFUSION instationnaires pour les "
         u"charges libres (dérive dans le champ + diffusion thermique), "
         u"comportant des termes sources non constants en temps et en "
         u"espace : piégeage, dépiégeage thermique et recombinaison ;",
         "puce"),
        (u"des équations de RÉACTION (locales en espace) pour les charges "
         u"piégées.", "puce"),
        ("", "texte"),
        ("Constantes et paramètres (valeurs par défaut du code) :", "titre"),
        (u"kB = 1,381×10⁻²³ J/K ;  q = 1,6×10⁻¹⁹ C ;  A = 1,2×10⁶ A m⁻²K⁻² ;"
         u"  ν = 6,2×10¹² s⁻¹", "puce"),
        (u"ε₀ = 8,85×10⁻¹² F/m ;  εr = 2,3 ;  ε = ε₀εr ;  T = 293 K ;"
         u"  L = 200 µm", "puce"),
        (u"w_e = w_h = 1,12 eV ;  w_me = w_mh = 0,65 eV ;"
         u"  a_e = a_h = 3,8×10⁻⁹ m", "puce"),
        (u"S₁ = S₂ = S₃ = 1 ;  B_e = B_h = 0,1 s⁻¹ ;  D_e = D_h = 1 ;"
         u"  n₀et = n₀ht = 100 C/m³ ;  w_tre = w_trh = 0,6 eV", "puce"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


def page_equations(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96, "2. Équations du modèle",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        ("2.1 Équation de Poisson", "titre"),
        (r"$\rho(x,t) = -n_{e\mu} - n_{et} + n_{h\mu} + n_{ht}$", "eq"),
        (r"$\dfrac{\partial^2 V}{\partial x^2} = -\dfrac{\rho}{\varepsilon_0"
         r"\varepsilon_r}\qquad E(x,t) = -\dfrac{\partial V}{\partial x}$",
         "eq"),
        ("(le signe est celui de la forme physique  div E = ρ/ε  avec "
         "E = −grad V)", "texte"),
        ("", "texte"),
        ("2.2 Courants de dérive — diffusion", "titre"),
        (r"$j_e = -\,n_{e\mu}\,\mu_e E - D_e\,\partial_x n_{e\mu}$"
         r"$\qquad j_h = +\,n_{h\mu}\,\mu_h E - D_h\,\partial_x n_{h\mu}$",
         "eq"),
        (u"les TROUS dérivent dans le sens de E, les ÉLECTRONS à "
         u"contre-champ ;", "texte"),
        (r"$D_{e,h} = \dfrac{k_BT}{q}\,\mu_{e,h}$   (relation d'Einstein)",
         "eq"),
        ("", "texte"),
        ("2.3 Mobilité dépendante du champ (avec limite à champ nul)", "titre"),
        (r"$\mu_{e,h}(E) = \dfrac{2\nu a_{e,h}}{|E|}\,"
         r"e^{-q w_{m e,h}/k_B T}\,\sinh\!\left(\dfrac{q E a_{e,h}}"
         r"{2 k_B T}\right)$", "eq"),
        (r"quand $E \to 0$ : $\ \mu \to \nu a^2 q\,/(k_B T)\,"
         r"e^{-q w_m/k_B T}$ (valeur finie, utilisée sous un seuil pour "
         r"éviter la division par zéro)", "texte"),
        ("", "texte"),
        ("2.4 Conservation des charges libres (advection / diffusion + "
         "sources)", "titre"),
        (r"$\partial_t n_{e\mu} + \partial_x j_e = -S_1\dfrac{\mu_e}"
         r"{\varepsilon}n_{ht}n_{e\mu} - S_3\dfrac{\mu_e+\mu_h}"
         r"{\varepsilon}n_{h\mu}n_{e\mu} - B_e n_{e\mu}\!\left(1-"
         r"\dfrac{n_{et}}{n_{0et}}\right) + \nu n_{et}e^{-q w_{tre}/k_BT}$",
         "eq"),
        (r"$\partial_t n_{h\mu} + \partial_x j_h = -S_2\dfrac{\mu_h}"
         r"{\varepsilon}n_{et}n_{h\mu} - S_3\dfrac{\mu_e+\mu_h}"
         r"{\varepsilon}n_{h\mu}n_{e\mu} - B_h n_{h\mu}\!\left(1-"
         r"\dfrac{n_{ht}}{n_{0ht}}\right) + \nu n_{ht}e^{-q w_{trh}/k_BT}$",
         "eq"),
        ("", "texte"),
        ("2.5 Évolution des charges piégées (réaction)", "titre"),
        (r"$\partial_t n_{et} = -S_1\dfrac{\mu_e}{\varepsilon}n_{ht}"
         r"n_{e\mu} + B_e n_{e\mu}\!\left(1-\dfrac{n_{et}}{n_{0et}}\right)"
         r" - \nu n_{et}e^{-q w_{tre}/k_BT}$", "eq"),
        (r"$\partial_t n_{ht} = -S_2\dfrac{\mu_h}{\varepsilon}n_{et}"
         r"n_{h\mu} + B_h n_{h\mu}\!\left(1-\dfrac{n_{ht}}{n_{0ht}}\right)"
         r" - \nu n_{ht}e^{-q w_{trh}/k_BT}$", "eq"),
        ("Les termes S₁/S₂/S₃ modélisent les recombinaisons (libre↔piégé et "
         "libre↔libre), B le piégeage dans les pièges de capacité n₀, et "
         "ν·exp(−qw_tr/k_BT) le dépiégeage thermique.", "texte"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


def page_cl(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96, "3. Conditions aux limites et injections",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        ("3.1 Potentiel", "titre"),
        (u"V(0, t) = Va(t) à l'anode : Va constant (= va0), rampe "
         u"(= pente·t) ou sinusoïde (va0 + amplitude·sin(2πf t)) — "
         u"choisi par p.va_type ;", "puce"),
        (u"V(L, t) = 0 à la cathode.", "puce"),
        ("", "texte"),
        ("3.2 Injection d'électrodes (Richardson — Schottky)", "titre"),
        (u"Les TROUS sont injectés à l'ANODE et les ÉLECTRONS à la "
         u"CATHODE, uniquement si le champ local E a le bon signe "
         u"(sinon l'électrode est bloquante) :", "texte"),
        (r"$j(w,E) = \dfrac{A T^2}{q}\,e^{-q w/k_B T}\left[\exp\!\left("
         r"\dfrac{q}{k_B T}\sqrt{\dfrac{qE}{4\pi\varepsilon}}\right)"
         r"-1\right]$", "eq"),
        (u"avec w = w_h à l'anode et w = w_e à la cathode ; j = 0 si "
         u"E ≤ 0.", "texte"),
        ("", "texte"),
        (u"L'injection est ensuite LIMITÉE à θ·n₀·Δx_bord par pas de temps "
         u"(θ = 0,05) : pendant le tout début de la simulation, le champ "
         u"n'est pas encore écranté et le courant de Schottky serait "
         u"immense ; cette limite stabilise le transitoire d'écran de "
         u"charge. Dès que l'accumulation aux électrodes écrante le champ "
         u"(E → 0), l'injection décroît naturellement (rétroaction de "
         u"Poisson).", "texte"),
        ("", "texte"),
        ("3.3 Porteurs atteignant l'électrode opposée", "titre"),
        (u"Un porteur qui atteint l'électrode opposée en est EXTRAIT "
         u"(flux sortant upwind) ; l'électrode est bloquante pour "
         u"l'entrée de l'autre signe. Aucun flux diffusif n'échange "
         u"avec l'extérieur.", "texte"),
        ("", "texte"),
        ("3.4 Courant mesuré aux électrodes", "titre"),
        (u"J = j_h − j_e (les densités sont des densités de charge, "
         u"les flux sont directement en A/m²). Le code enregistre le "
         u"courant total à l'anode et à la cathode en fonction du temps.",
         "texte"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


def page_numerique(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96,
             "4. Méthode numérique — volumes finis et schéma MUSCL",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        ("4.1 Formulation volumes finis", "titre"),
        (u"Pour chaque volume de contrôle i (largeur Δxᵢ), l'équation de "
         u"conservation  ∂ₜn + ∂ₓφ = S  est intégrée :", "texte"),
        (r"$\dfrac{n_i^{k+1}-n_i^k}{\Delta t} = S_i - \dfrac{\phi_{i+1/2}"
         r"- \phi_{i-1/2}}{\Delta x_i}$", "eq"),
        (u"Le flux φ est la somme d'un flux ADVECTIF (dérive) et d'un flux "
         u"DIFFUSIF.", "texte"),
        ("", "texte"),
        ("4.2 Schéma MUSCL (advection)", "titre"),
        (u"Dans chaque cellule, on reconstruit une pente limitée "
         u"(limiteur MINMOD) :", "texte"),
        (r"$s_i = \mathrm{minmod}\!\left(\dfrac{n_{i+1}-n_i}{\Delta x_i},"
         r"\ \dfrac{n_i-n_{i-1}}{\Delta x_i}\right)$", "eq"),
        (u"Les états reconstruits de part et d'autre de la face sont "
         u"alors :", "texte"),
        (r"$n_L = n_i + s_i\Delta x_i/2$"
         r"$\qquad n_R = n_{i+1} - s_{i+1}\Delta x_{i+1}/2$", "eq"),
        (u"et le flux advectif est UPWIND sur ces états : "
         u"φ = u⁺n_L + u⁻n_R (u vitesse de dérive, u⁺ = max(u,0)). "
         u"Le limiteur rend le schéma du 2ᵉ ordre là où la solution est "
         u"lisse et du 1ᵉʳ ordre près des fronts, sans oscillations "
         u"parasites.", "texte"),
        ("", "texte"),
        ("4.3 Diffusion", "titre"),
        (u"Le flux diffusif aux faces internes est centré : "
         u"φ_d = −D_face (n_R − n_L)/distance entre centres.", "texte"),
        ("", "texte"),
        ("4.4 Intégration temporelle RK2 (Heun)", "titre"),
        (u"À chaque pas : (1) résolution de Poisson → champ E ; "
         u"(2) calcul des flux MUSCL et des sources ; (3) étage 1 de RK2 ; "
         u"(4) résolution de Poisson sur l'état intermédiaire ; "
         u"(5) étage 2 ; (6) moyenne. Les densités sont contraintes "
         u"positives et les pièges bornés par leur capacité n₀.", "texte"),
        ("", "texte"),
        ("4.5 Résolution de Poisson", "titre"),
        (u"La discrétisation volumes finis de V'' = −ρ/ε conduit à un "
         u"système TRIDIAGONAL en le potentiel aux faces, résolu par "
         u"l'ALGORITHME DE THOMAS (coût linéaire en N). Le champ est "
         u"ensuite calculé aux centres des cellules : "
         u"E_i = −(V_{i+1} − V_i)/Δx_i.", "texte"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


def page_maillage_dt(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96,
             "5. Maillage non uniforme et pas de temps adaptatif",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        ("5.1 Maillage raffiné aux électrodes", "titre"),
        (u"Les gradients les plus forts (écran de charge, injection) sont "
         u"près des électrodes : le maillage y est donc plus fin. On "
         u"utilise un mapping TANH du paramètre uniforme f ∈ [0,1] :", "texte"),
        (r"$x(f) = \dfrac{L}{2}\left[1 + \dfrac{\tanh\!\left(g(2f-1)\right)}"
         r"{\tanh g}\right]$", "eq"),
        (u"avec g = p.etirement (défaut 4). Les N cellules sont fines aux "
         u"deux bords (≈ L·g·sech²g /(N·tanh g) ≈ 0,1 µm) et grossissent "
         u"régulièrement vers le centre (≈ 8 µm).", "texte"),
        ("", "texte"),
        ("5.2 Pas de temps adaptatif", "titre"),
        (u"Le pas de temps est recalculé à chaque itération comme le "
         u"minimum de :", "texte"),
        (u"CFL advection/diffusion :  Δt ≤ cfl·Δxᵢ / (|u| + 2D/Δxᵢ) ;",
         "puce"),
        (u"limite de flux entrant dans chaque cellule (raideur de "
         u"l'injection) ;", "puce"),
        (u"raideur des sources : Δt ≤ cfl·n/|S| (piégeage, dépiégeage, "
         u"recombinaison) ;", "puce"),
        (u"bornes : fin de simulation et instants d'enregistrement des "
         u"instantanés.", "puce"),
        ("", "texte"),
        (u"Au démarrage, le transitoire d'écran de charge (≈ 10⁻¹⁶ s) "
         u"impose des pas minuscules ; dès que le champ est écranté, "
         u"le pas croît vers la seconde. Une simulation de 50 s demande "
         u"typiquement ~15 000 itérations.", "texte"),
        ("", "texte"),
        ("5.3 Stabilité et robustesse", "titre"),
        (u"positivité des densités libres (clamp n ≥ 0) ;", "puce"),
        (u"pièges bornés : 0 ≤ n_et ≤ n₀et, 0 ≤ n_ht ≤ n₀ht ;", "puce"),
        (u"limiteur minmod (TVB) : pas d'oscillations numériques sur les "
         u"fronts de charge ;", "puce"),
        (u"injection limitée en charge par pas (cf. §3.2).", "puce"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


def page_utilisation(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96, "6. Utilisation du code et sorties",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        ("6.1 Lancement", "titre"),
        (u"python3 transport_charges.py", "code"),
        (u"Tous les paramètres sont regroupés dans la classe Parametres "
         u"(constantes physiques, géométrie, coefficients S/B/n₀, "
         u"conditions aux limites, paramètres numériques). Exemple de "
         u"personnalisation :", "texte"),
        (u"import transport_charges as tc\n"
         u"p = tc.Parametres()\n"
         u"p.va_type = 'sinus'      # ou 'constante' / 'rampe'\n"
         u"p.va0, p.va_amp, p.va_freq = 500., 300., 0.05\n"
         u"p.tmax = 100.            # durée [s]\n"
         u"res = tc.simuler(p)      # simulation\n"
         u"tc.tracer(res, p)        # figures", "code"),
        ("", "texte"),
        ("6.2 Sorties (dossier resultats/)", "titre"),
        (u"resultats.npz : maillage (xc, xf), instantanés à ~25 instants "
         u"des 4 densités, du potentiel V(x,t) et du champ E(x,t), "
         u"historiques temporels (charges totales intégrées, potentiel "
         u"d'anode, courants aux électrodes) ;", "puce"),
        (u"profils_charges.png : n_eµ, n_hµ, n_et, n_ht en fonction de x "
         u"à différents instants ;", "puce"),
        (u"potentiel_champ.png : V(x,t) et E(x,t) ;", "puce"),
        (u"evolutions_temporelles.png : charges totales, Va(t), courants "
         u"d'électrodes ;", "puce"),
        (u"maillage.png : répartition des faces du maillage.", "puce"),
        ("", "texte"),
        ("6.3 Comportement physique attendu", "titre"),
        (u"Au tout début, le champ est uniforme (E ≈ Va/L) : l'injection "
         u"est forte ; les charges s'accumulent devant les électrodes et "
         u"ÉCRANTENT le champ local, ce qui auto-limite l'injection "
         u"(rétroaction Poisson ↔ transport).", "puce"),
        (u"Les porteurs dérivent lentement (mobilité ~10⁻¹⁴ m²/V/s), se "
         u"piègent (B), se dépiègent thermiquement (ν·e^(−qw_tr/kT)) et "
         u"recombinent (S₁, S₂, S₃) ; les pièges se remplissent "
         u"progressivement, ce qui est le mécanisme de vieillissement "
         u"modélisé.", "puce"),
        (u"Compte tenu des constantes fournies, le régime est dominé par "
         u"l'accumulation près des électrodes sur des durées de l'ordre "
         u"de la seconde à la minute.", "puce"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


def page_structure(pdf):
    fig = nouvelle_page()
    fig.text(0.07, 0.96, "7. Structure du programme",
             fontsize=15, fontweight="bold", color="#1a3d6d")
    y = 0.90
    y = bloc_texte(fig, 0.07, y, [
        (u"transport_charges.py — sous-programmes :", "titre"),
        (u"Parametres : classe regroupant toutes les données d'entrée ;",
         "puce"),
        (u"potentiel_anode(t, p) : Va(t) — constante / rampe / sinus ;",
         "puce"),
        (u"maillage(p) : maillage tanh non uniforme (xc, xf, dx) ;", "puce"),
        (u"resoudre_poisson(rho, xf, dx, Va, p) : tridiagonale de Thomas "
         u"→ V aux faces, E aux centres ;", "puce"),
        (u"mobilite(E, a, w_m, p) : µ(E) avec limite analytique à champ "
         u"faible/nul ;", "puce"),
        (u"courant_injection(E_face, w, p) : courant de Richardson–"
         u"Schottky (avec garde de signe) ;", "puce"),
        (u"minmod(a, b) : limiteur de pente ;", "puce"),
        (u"flux_interieurs(n, u, D, dx, xf) : flux MUSCL advectif + "
         u"diffusion centrée aux faces internes ;", "puce"),
        (u"termes_sources(...) : piégeage / dépiégeage / recombinaison "
         u"des 4 équations ;", "puce"),
        (u"etat_derivees(...) : Poisson → E, mobilités, flux, conditions "
         u"aux limites d'injection, bilan de volumes finis ;", "puce"),
        (u"pas_temps(...) : pas adaptatif (CFL, flux, sources) ;", "puce"),
        (u"simuler(p) : boucle RK2, enregistrements npz + instantanés ;",
         "puce"),
        (u"tracer(res, p) : figures PNG.", "puce"),
        ("", "texte"),
        (u"Remarque sur les conventions de signe", "titre"),
        (u"Avec l'anode à gauche (Va > 0), le champ E = −∂V/∂x pointe "
         u"vers +x (de l'anode vers la cathode). Les équations du cahier "
         u"des charges écrivent les deux courants avec le même signe "
         u"(+µE) : c'est la convention de la géométrie cathode-à-gauche. "
         u"Pour la géométrie demandée (anode à gauche), le code utilise "
         u"la forme physiquement cohérente : j_h = +n_hµ_h E − D∂ₓn_h "
         u"pour les trous et j_e = −n_eµ_e E − D∂ₓn_e pour les électrons, "
         u"de sorte que chaque porteur dérive bien de son électrode "
         u"d'injection vers l'électrode opposée. Toutes les autres "
         u"équations (sources, piégeage, injection) sont implémentées "
         u"littéralement.", "texte"),
    ])
    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages("explication.pdf") as pdf:
        page_titre(pdf)
        page_modele(pdf)
        page_equations(pdf)
        page_cl(pdf)
        page_numerique(pdf)
        page_maillage_dt(pdf)
        page_utilisation(pdf)
        page_structure(pdf)
    print("explication.pdf généré.")
