# -*- coding: utf-8 -*-
"""
=============================================================================
 TRANSPORT DE CHARGES 1D INSTATIONNAIRE DANS UN MATERIAU DIELECTRIQUE
=============================================================================

Modèle bipolaire (électrons / trous, libres / piégés) type "Le Roy" :
  - équation de Poisson pour le potentiel V(x,t) et le champ E(x,t),
  - équations d'advection/diffusion (dérive - diffusion) pour les charges
    libres, avec termes sources de piégeage / dépiégeage / recombinaison,
  - équations de réaction (ODS locales) pour les charges piégées.

Discrétisation :
  - volumes finis sur un maillage 1D NON UNIFORME (raffiné aux électrodes),
  - schéma MUSCL (reconstruction linéaire à pente limitée "minmod")
    pour la partie advective, flux centré pour la diffusion,
  - intégration temporelle RK2 (Heun / SSP), pas de temps adaptatif,
  - Poisson résolu à chaque étage par l'algorithme de Thomas (tridiagonal).

Géométrie :  ANODE à gauche (x = 0, potentiel Va(t) imposé,
             configurable : constante / rampe / sinusoïde),
             CATHODE à droite (x = L, potentiel nul).

Conventions de signe (important, voir le PDF d'explication) :
  E = -dV/dx. Avec Va > 0 le champ pointe vers +x (de l'anode vers la
  cathode). Les TROUS dérivent dans le sens de E (flux  +mu_h*E*n),
  les ELECTRONS dérivent à contre-champ (flux  -mu_e*E*n), ce qui est
  la forme physiquement cohérente des équations du cahier des charges
  pour la géométrie anode-à-gauche / cathode-à-droite.

  Les densités n sont des densités de CHARGE en C/m^3 (positives en
  magnitude). La densité volumique de charge nette vaut :
      rho = -n_emu - n_et + n_hmu + n_ht
  et l'équation de Poisson s'écrit  d2V/dx2 = -rho/eps.

Injections aux électrodes (courant de Richardson - Schottky) :
  - TROUS injectés à l'ANODE (x=0)  si E(0,t) > 0,
  - ELECTRONS injectés à la CATHODE (x=L) si E(L,t) > 0,
  - injection nulle sinon (dépendance au signe du champ demandée).

Sorties : toutes les densités (libres et piégées), le potentiel et le
champ électrique, en fonction du temps ET de l'espace (fichier .npz +
figures PNG), courants aux électrodes et charges totales en fonction
du temps.

Auteur : généré pour le projet "transport de charges dans les diélectriques".
=============================================================================
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# 1. PARAMETRES ET CONSTANTES PHYSIQUES
# ============================================================================

class Parametres:
    """Toutes les constantes et données d'entrée du modèle."""

    def __init__(self):
        # --- constantes physiques -----------------------------------------
        self.kB   = 1.381e-23      # constante de Boltzmann [J/K]
        self.q    = 1.6e-19        # charge élémentaire [C]
        self.A    = 1.2e6          # constante de Richardson [A/m^2/K^2]
        self.eps0 = 8.85e-12       # permittivité du vide [F/m]
        self.epsR = 2.3            # permittivité relative du matériau [-]
        self.eps  = self.eps0 * self.epsR  # permittivité du matériau [F/m]
        self.nu   = 6.2e12         # fréquence de saut / phonon [s^-1]

        # --- géométrie et température --------------------------------------
        self.L    = 200e-6         # épaisseur du matériau [m]
        self.TK   = 20.0 + 273.0   # température [K]

        # --- barrières et pièges -------------------------------------------
        self.w_e  = 1.12           # barrière d'injection électrons [eV]
        self.w_h  = 1.12           # barrière d'injection trous [eV]
        self.w_me = 0.65           # barrière de mobilité électrons [eV]
        self.w_mh = 0.65           # barrière de mobilité trous [eV]
        self.a_e  = 3.8e-9         # distance inter-sites électrons [m]
        self.a_h  = 3.8e-9         # distance inter-sites trous [m]
        self.w_tre = 0.6           # barrière de dépiégeage électrons [eV]
        self.w_trh = 0.6           # barrière de dépiégeage trous [eV]

        # --- coefficients du modèle -----------------------------------------
        self.S1 = 1.0              # recombinaison libre e / piégé h
        self.S2 = 1.0              # recombinaison libre h / piégé e
        self.S3 = 1.0              # recombinaison libre e / libre h
        self.Be = 0.1              # capture électrons [s^-1]
        self.Bh = 0.1              # capture trous [s^-1]
        self.n0et = 100.0          # capacité des pièges à électrons [C/m^3]
        self.n0ht = 100.0          # capacité des pièges à trous [C/m^3]

        # --- conditions aux limites en potentiel ----------------------------
        # type : 'constante' | 'rampe' | 'sinus'
        self.va_type = 'constante'
        self.va0      = 1000.0     # amplitude / valeur [V]
        self.va_ramp  = 10.0       # vitesse de rampe [V/s]
        self.va_amp   = 500.0      # amplitude du sinus [V]
        self.va_freq  = 0.1        # fréquence du sinus [Hz]

        # --- paramètres numériques ------------------------------------------
        self.N        = 200        # nombre de volumes de contrôle
        self.etirement = 4.0       # paramètre g du mapping tanh (raffinement
                                   # aux électrodes ; plus g est grand, plus
                                   # le maillage est fin près des électrodes)
        self.tmax     = 50.0       # durée de la simulation [s]
        self.cfl_adv  = 0.4        # CFL advection / diffusion
        self.cfl_src  = 0.2        # CFL termes sources
        self.theta_inj = 0.05      # fraction de n0 injectée par pas (limite)
        self.n_snap   = 25         # nombre d'instantanés enregistrés


def potentiel_anode(t, p):
    """Potentiel imposé à l'anode Va(t) [V].

    Trois formes disponibles, choisies via p.va_type :
      - 'constante' : Va = va0,
      - 'rampe'     : Va = va_ramp * t (montée linéaire),
      - 'sinus'     : Va = va0 + va_amp * sin(2*pi*f*t).
    """
    if p.va_type == 'constante':
        return p.va0
    elif p.va_type == 'rampe':
        return p.va_ramp * t
    elif p.va_type == 'sinus':
        return p.va0 + p.va_amp * np.sin(2.0 * np.pi * p.va_freq * t)
    else:
        raise ValueError("Type de potentiel d'anode inconnu : %s" % p.va_type)


# ============================================================================
# 2. MAILLAGE 1D NON UNIFORME (raffiné aux deux électrodes)
# ============================================================================

def maillage(p):
    """Construit le maillage non uniforme (exactement N cellules).

    Mapping tanh : x(f) = L/2 * (1 + tanh(g*(2f-1))/tanh(g)), f dans [0,1].
    Les mailles sont fines aux deux électrodes et s'épaississent
    régulièrement vers le centre du domaine ; g = p.etirement contrôle
    le raffinement (dx au bord ~ L*g*sech^2(g)/(N*tanh(g))).

    Returns
    -------
    xc : (N,)   centres des cellules
    xf : (N+1,) faces (xf[0]=0 anode, xf[N]=L cathode)
    dx : (N,)   largeurs des cellules
    """
    N, g = p.N, p.etirement
    f = np.linspace(0.0, 1.0, N + 1)
    xf = p.L * 0.5 * (1.0 + np.tanh(g * (2.0 * f - 1.0)) / np.tanh(g))
    xf[-1] = p.L
    dx = np.diff(xf)
    xc = 0.5 * (xf[:-1] + xf[1:])
    return xc, xf, dx


# ============================================================================
# 3. EQUATION DE POISSON (algorithme de Thomas, tridiagonal)
# ============================================================================

def resoudre_poisson(rho, xf, dx, Va, p):
    """Résout  d2V/dx2 = -rho/eps  avec V(0)=Va, V(L)=0 (volumes finis).

    Inconnues : potentiel aux N+1 faces. Pour la cellule i :
        (V[i+1]-V[i])/h[i] - (V[i]-V[i-1])/h[i-1] = -rho[i]*dx[i]/eps
    (intégration de V'' sur la cellule ; h = distances entre faces).

    Returns
    -------
    V : (N+1,) potentiel aux faces
    E : (N,)   champ électrique aux centres des cellules, E = -dV/dx
    """
    N = len(dx)
    h = np.diff(xf)                       # (N,) distances face->face
    rhs = -rho * dx / p.eps               # second membre des cellules

    # coefficients tridiagonaux pour les cellules internes 1..N-1
    a = np.zeros(N - 1)                   # sous-diagonale
    b = np.zeros(N - 1)                   # diagonale
    c = np.zeros(N - 1)                   # sur-diagonale
    d = np.zeros(N - 1)                   # second membre
    for i in range(1, N):
        k = i - 1
        a[k] = -1.0 / h[i - 1]
        b[k] = 1.0 / h[i - 1] + 1.0 / h[i]
        c[k] = -1.0 / h[i]
        d[k] = rhs[i]

    # conditions aux limites de Dirichlet
    d[0] -= a[0] * Va
    d[-1] -= c[-1] * 0.0

    # ---- algorithme de Thomas (élimination avant / substitution arrière)
    cp = np.copy(c)
    dp = np.copy(d)
    for k in range(1, N - 1):
        m = a[k] / b[k - 1]
        b[k] -= m * cp[k - 1]
        dp[k] -= m * dp[k - 1]
    Vi = np.zeros(N - 1)
    Vi[-1] = dp[-1] / b[-1]
    for k in range(N - 3, -1, -1):
        Vi[k] = (dp[k] - cp[k] * Vi[k + 1]) / b[k]

    V = np.concatenate(([Va], Vi, [0.0]))
    # champ électrique aux CENTRES des cellules (convention décalée) :
    # E_i = -(V_{i+1} - V_i)/dx_i ; aux électrodes on utilise la maille de bord
    E = -np.diff(V) / h
    return V, E


# ============================================================================
# 4. MOBILITE DEPENDANTE DU CHAMP (avec limite des faibles champs)
# ============================================================================

def mobilite(E, a, w_m, p):
    """Mobilité  mu(E) = (2*nu*a/|E|) * exp(-q*w_m/(kB*T)) * sinh(q*E*a/(2*kB*T)).

    Pour E -> 0, sinh(xE)/E -> x, donc la mobilité tend vers la valeur
    finie  mu0 = nu*a^2*q/(kB*T) * exp(-q*w_m/(kB*T)).
    On utilise cette limite analytique sous un seuil |E| < E_min afin
    d'éviter toute division par zéro (champ faible ou nul).
    """
    E_min = 1.0 * p.kB * p.TK / (p.q * a)          # seuil : argument ~ 0.5
    mu0 = p.nu * a * a * p.q / (p.kB * p.TK) * np.exp(-p.q * w_m / (p.kB * p.TK))
    sgn = np.where(np.asarray(E) >= 0.0, 1.0, -1.0)
    Ebig = np.where(np.abs(E) >= E_min, E, sgn * E_min)
    mu = np.where(
        np.abs(E) >= E_min,
        (2.0 * p.nu * a / Ebig)
        * np.exp(-p.q * w_m / (p.kB * p.TK))
        * np.sinh(p.q * Ebig * a / (2.0 * p.kB * p.TK)),
        mu0,
    )
    return mu


# ============================================================================
# 5. COURANT D'INJECTION (Richardson - Schottky)
# ============================================================================

def courant_injection(E_face, w, p):
    """Courant d'injection à une électrode [A/m^2] (grandeur signée >= 0).

        j = (A*T^2/q) * exp(-q*w/(kB*T))
            * ( exp( (q/(kB*T)) * sqrt(q*E/(4*pi*eps)) ) - 1 )

    L'injection n'a lieu QUE si le champ local pousse les charges vers
    l'intérieur du matériau (E > 0 au sens de la géométrie choisie) ;
    sinon j = 0 (dépendance au signe du champ).
    """
    E = max(E_face, 0.0)
    pref = p.A * p.TK ** 2 / p.q * np.exp(-p.q * w / (p.kB * p.TK))
    arg = (p.q / (p.kB * p.TK)) * np.sqrt(p.q * E / (4.0 * np.pi * p.eps))
    return pref * (np.expm1(arg))


# ============================================================================
# 6. FLUX MUSCL (advection limitee + diffusion centree)
# ============================================================================

def minmod(a, b):
    """Limiteur de pente minmod : renvoie a si a,b de même signe (|a| min)."""
    same = (a * b) > 0.0
    return np.where(same, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def flux_interieurs(n, u, D, dx, xf):
    """Flux aux faces internes pour l'équation d'advection/diffusion.

    Reconstruction MUSCL : dans la cellule i, pente limitée
        s_i = minmod(n[i+1]-n[i], n[i]-n[i-1]) / dx[i]
    états reconstruits de part et d'autre de la face f :
        n_L = n[i]   + s_i   * dx[i]/2      (côté gauche)
        n_R = n[i+1] - s_{i+1}* dx[i+1]/2   (côté droite)
    Flux advectif upwind sur les états reconstruits :
        phi = u^+ * n_L + u^- * n_R
    Flux diffusif centré :
        phi -= D_face * (n_R - n_L) / distance_entre_centres

    Parameters
    ----------
    n : (N,)  densité dans les cellules
    u : (N,)  vitesse de dérive aux CENTRES (u = +/- mu*E interpolé)
    D : (N,)  coefficient de diffusion aux centres

    Returns
    -------
    phi : (N+1,) flux aux faces (partie interne ; faces de bord = 0,
          remplies ensuite par les conditions d'injection).
    """
    N = len(n)
    phi = np.zeros(N + 1)

    # pentes limitées (minmod) ; conditions de bord : pente nulle
    dm = np.zeros(N)
    dp = np.zeros(N)
    dm[1:] = n[1:] - n[:-1]
    dp[:-1] = n[1:] - n[:-1]
    s = minmod(dp, dm) / dx                      # pente par cellule

    # états reconstruits de part et d'autre de chaque face interne
    nL = n[:-1] + 0.5 * s[:-1] * dx[:-1]          # (N-1,) face 1..N-1
    nR = n[1:] - 0.5 * s[1:] * dx[1:]            # (N-1,)

    # vitesse et diffusion aux faces internes (moyenne harmonique simple)
    uf = 0.5 * (u[:-1] + u[1:])
    Df = 0.5 * (D[:-1] + D[1:])
    dist = 0.5 * (dx[:-1] + dx[1:])              # centre -> centre

    up = np.maximum(uf, 0.0)
    um = np.minimum(uf, 0.0)
    phi[1:N] = up * nL + um * nR - Df * (nR - nL) / dist
    return phi


# ============================================================================
# 7. TERMES SOURCES (piégeage / dépiégeage / recombinaison)
# ============================================================================

def termes_sources(n_emu, n_hmu, n_et, n_ht, mu_e, mu_h, p):
    """Termes sources des 4 équations de transport (formes du cahier des charges).

    Rappel : les équations s'écrivent  dn/dt = S - dphi/dx ; ici on ne
    retourne que S. Les pièges sont limités à leur capacité n0et/n0ht
    (facteur 1 - n_t/n0), les densités sont en C/m^3.
    """
    # coefficients thermiquement activés
    kT = p.kB * p.TK
    detr_e = p.nu * np.exp(-p.q * p.w_tre / kT)   # dépiégeage électrons [s^-1]
    detr_h = p.nu * np.exp(-p.q * p.w_trh / kT)   # dépiégeage trous [s^-1]

    # recombinaisons (coefficients S du modèle, divisés par eps)
    rec_e_ht = p.S1 * mu_e / p.eps * n_ht * n_emu        # libre e / piégé h
    rec_h_et = p.S2 * mu_h / p.eps * n_et * n_hmu        # libre h / piégé e
    rec_e_h  = p.S3 * (mu_e + mu_h) / p.eps * n_hmu * n_emu  # libre / libre

    # piégeage / dépiégeage
    trap_e = p.Be * n_emu * (1.0 - n_et / p.n0et)
    trap_h = p.Bh * n_hmu * (1.0 - n_ht / p.n0ht)
    detrap_e = detr_e * n_et
    detrap_h = detr_h * n_ht

    # ---- équations --------------------------------------------------------
    S_emu = -rec_e_ht - rec_e_h - trap_e + detrap_e   # électrons libres
    S_hmu = -rec_h_et - rec_e_h - trap_h + detrap_h  # trous libres
    S_et  = -rec_e_ht + trap_e - detrap_e             # électrons piégés
    S_ht  = -rec_h_et + trap_h - detrap_h             # trous piégés
    return S_emu, S_hmu, S_et, S_ht


# ============================================================================
# 8. UN PAS DE TEMPS COMPLETE (transport + sources, RK2)
# ============================================================================

def etat_derivees(n_emu, n_hmu, n_et, n_ht, xf, dx, t, p, dt=None):
    """Calcule E, les flux aux faces et les sources pour un état donné.

    Si dt est fourni, le flux d'injection est LIMITÉ de façon à ne pas
    injecter plus de theta*n0*piège par pas dans la maille de bord :
        j_eff = min(j_inj, theta*n0*dx_bord/dt)
    (stabilise la phase transitoire d'écran de charge sans réduire dt).

    Retourne (dn_emu, dn_hmu, dn_et, dn_ht, E_centres, phi_h, phi_e).
    """
    N = len(dx)
    Va = potentiel_anode(t, p)

    # --- Poisson -> champ aux faces ----------------------------------------
    rho = -n_emu - n_et + n_hmu + n_ht
    V, E = resoudre_poisson(rho, xf, dx, Va, p)

    # --- mobilités et diffusion aux centres (interpolation des faces) -----
    E_c = E                                       # champ aux centres (N)
    mu_e = mobilite(E_c, p.a_e, p.w_me, p)
    mu_h = mobilite(E_c, p.a_h, p.w_mh, p)
    DIFF_e = p.kB * p.TK / p.q * mu_e               # relation d'Einstein
    DIFF_h = p.kB * p.TK / p.q * mu_h
    De = DIFF_e
    Dh = DIFF_h

    # --- vitesses de dérive aux centres ------------------------------------    # trous : +mu*E ; électrons : -mu*E (dérive à contre-champ)
    u_h = mu_h * E_c
    u_e = -mu_e * E_c

    # --- flux internes (MUSCL) ----------------------------------------------
    phi_h = flux_interieurs(n_hmu, u_h, Dh, dx, xf)
    phi_e = flux_interieurs(n_emu, u_e, De, dx, xf)

    # --- conditions aux limites : injections Richards- Schottky ------------
    # anode (x=0) : trous injectés si E(0) > 0 ; cathode (x=L) :
    # électrons injectés si E(L) > 0 ; sinon l'électrode est bloquante.
    # Le porteur atteignant l'électrode opposée en est extrait (sortie
    # upwind) ; aucune diffusion n'échange avec l'extérieur.
    if E_c[0] > 0.0:
        jh = courant_injection(E_c[0], p.w_h, p)          # flux entrant (+x)
        if dt is not None:
            jh = min(jh, p.theta_inj * p.n0ht * dx[0] / dt)
        phi_h[0] = jh
    else:
        phi_h[0] = min(u_h[0] * n_hmu[0], 0.0)            # extraction seule
    if E_c[-1] > 0.0:
        je = courant_injection(E_c[-1], p.w_e, p)         # entrée en -x
        if dt is not None:
            je = min(je, p.theta_inj * p.n0et * dx[-1] / dt)
        phi_e[N] = -je
    else:
        phi_e[N] = max(u_e[-1] * n_emu[-1], 0.0)          # extraction seule
    phi_e[0] = min(u_e[0] * n_emu[0], 0.0)                  # anode bloquante e
    phi_h[N] = max(u_h[-1] * n_hmu[-1], 0.0)               # cathode bloquante h

    # --- bilan de volume fini : dn/dt = S - (phi[i+1]-phi[i])/dx -----------
    S_emu, S_hmu, S_et, S_ht = termes_sources(
        n_emu, n_hmu, n_et, n_ht, mu_e, mu_h, p)
    div = lambda phi: (phi[1:] - phi[:-1]) / dx
    dn_emu = S_emu - div(phi_e)
    dn_hmu = S_hmu - div(phi_h)
    dn_et = S_et
    dn_ht = S_ht

    return dn_emu, dn_hmu, dn_et, dn_ht, E, phi_h, phi_e


def pas_temps(n_emu, n_hmu, n_et, n_ht, dn, dx, E, phi_h, phi_e, p):
    """Pas de temps adaptatif : min(CFL advection/diffusion, CFL sources,
    limite d'injection, limite de flux entrant)."""
    N = len(dx)
    E_c = E
    mu_e = mobilite(E_c, p.a_e, p.w_me, p)
    mu_h = mobilite(E_c, p.a_h, p.w_mh, p)
    u = np.maximum(np.abs(mu_e * E_c), np.abs(mu_h * E_c))
    D = p.kB * p.TK / p.q * np.maximum(mu_e, mu_h)

    # 1) advection + diffusion explicites
    rate = u + 2.0 * D / dx
    dt_adv = p.cfl_adv * np.min(dx / rate)

    # 2) flux entrants (raideur de l'injection d'électrode)
    inflow = (np.maximum(phi_h[:-1], 0.0) + np.maximum(-phi_e[1:], 0.0)
              + np.maximum(-phi_h[1:], 0.0) + np.maximum(phi_e[:-1], 0.0))
    dt_flux = np.min((n_hmu + n_emu + 1e-3 * (p.n0et + p.n0ht)) * dx
                     / np.maximum(inflow, 1e-30))

    # 3) raideur des sources (piégeage / dépiégeage / recombinaison)
    dt_src = np.inf
    for n, dnk in ((n_emu, dn[0]), (n_hmu, dn[1]),
                   (n_et, dn[2]), (n_ht, dn[3])):
        actif = np.abs(dnk) > 0.0
        if np.any(actif):
            with np.errstate(over='ignore', divide='ignore'):
                dt_src = min(dt_src, p.cfl_src * np.min(
                    (n[actif] + 1e-3 * p.n0et) / np.abs(dnk[actif])))

    return min(dt_adv, p.cfl_adv * dt_flux, dt_src)


# ============================================================================
# 9. SIMULATION COMPLETE
# ============================================================================

def simuler(p, dossier="resultats"):
    """Boucle en temps RK2. Retourne un dictionnaire de résultats."""
    import os
    os.makedirs(dossier, exist_ok=True)

    xc, xf, dx = maillage(p)
    N = len(dx)

    # ---- conditions initiales : matériau vierge (aucune charge) ------------
    n_emu = np.zeros(N)
    n_hmu = np.zeros(N)
    n_et = np.zeros(N)
    n_ht = np.zeros(N)

    # ---- structures d'enregistrement ---------------------------------------
    t_snap = np.linspace(0.0, p.tmax, p.n_snap)
    i_snap = 1
    snaps = {'t': [], 'n_emu': [], 'n_hmu': [], 'n_et': [], 'n_ht': [],
             'V': [], 'E': []}
    hist = {'t': [0.0], 'Q_emu': [0.0], 'Q_hmu': [0.0], 'Q_et': [0.0],
            'Q_ht': [0.0], 'Va': [potentiel_anode(0.0, p)],
            'j_anode': [0.0], 'j_cathode': [0.0]}

    def enregistrer(t, E, V):
        snaps['t'].append(t)
        for key, arr in (('n_emu', n_emu), ('n_hmu', n_hmu),
                        ('n_et', n_et), ('n_ht', n_ht)):
            snaps[key].append(arr.copy())
        snaps['V'].append(V.copy())
        snaps['E'].append(E.copy())

    # état initial
    rho = np.zeros(N)
    V0, E0 = resoudre_poisson(rho, xf, dx, potentiel_anode(0.0, p), p)
    enregistrer(0.0, E0, V0)

    t = 0.0
    n_iter = 0
    dt_prev = None
    while t < p.tmax:
        # ---- étage 1 de RK2 (injection limitée par le pas précédent) ----------
        d1 = etat_derivees(n_emu, n_hmu, n_et, n_ht, xf, dx, t, p, dt=dt_prev)
        dt = pas_temps(n_emu, n_hmu, n_et, n_ht, d1, dx, d1[4], d1[5], d1[6], p)
        dt = min(dt, p.tmax - t)
        if i_snap < len(t_snap) and t + dt > t_snap[i_snap]:
            dt = t_snap[i_snap] - t
        n1 = (n_emu + dt * d1[0], n_hmu + dt * d1[1],
              n_et + dt * d1[2], n_ht + dt * d1[3])

        # ---- étage 2 de RK2 (injection limitée par le pas courant) ------------
        d2 = etat_derivees(n1[0], n1[1], n1[2], n1[3], xf, dx, t + dt, p, dt=dt)
        dt_prev = dt
        n_emu = n_emu + 0.5 * dt * (d1[0] + d2[0])
        n_hmu = n_hmu + 0.5 * dt * (d1[1] + d2[1])
        n_et = np.clip(n_et + 0.5 * dt * (d1[2] + d2[2]), 0.0, p.n0et)
        n_ht = np.clip(n_ht + 0.5 * dt * (d1[3] + d2[3]), 0.0, p.n0ht)
        n_emu = np.maximum(n_emu, 0.0)   # positivité des densités
        n_hmu = np.maximum(n_hmu, 0.0)

        t += dt
        n_iter += 1

        # ---- sorties temporelles -------------------------------------------
        E = d2[4]
        j_ano = d2[5][0] - d2[6][0]      # courant total à l'anode
        j_cath = d2[5][N] - d2[6][N]     # courant total à la cathode
        hist['t'].append(t)
        hist['Q_emu'].append(np.sum(n_emu * dx))
        hist['Q_hmu'].append(np.sum(n_hmu * dx))
        hist['Q_et'].append(np.sum(n_et * dx))
        hist['Q_ht'].append(np.sum(n_ht * dx))
        hist['Va'].append(potentiel_anode(t, p))
        hist['j_anode'].append(j_ano)
        hist['j_cathode'].append(j_cath)

        # ---- instantanés -----------------------------------------------------
        if i_snap < len(t_snap) and t >= t_snap[i_snap] - 1e-12:
            # potentiel et champ recalculés pour l'instantané
            rho = -n_emu - n_et + n_hmu + n_ht
            V, E = resoudre_poisson(rho, xf, dx, potentiel_anode(t, p), p)
            enregistrer(t, E, V)
            i_snap += 1

        if n_iter % 5000 == 0:
            print("  t = %.3e s  (pas %d, dt = %.3e s)" % (t, n_iter, dt))

    print("Simulation terminée : %d itérations, t = %.3e s" % (n_iter, t))

    # ---- assemblage des résultats ------------------------------------------
    res = {'xc': xc, 'xf': xf, 'dx': dx, 'snaps': snaps, 'hist': hist,
           'n_iter': n_iter}
    np.savez(os.path.join(dossier, "resultats.npz"),
             xc=xc, xf=xf,
             t_snap=np.array(snaps['t']),
             **{"n_emu_snap": np.array(snaps['n_emu']),
                "n_hmu_snap": np.array(snaps['n_hmu']),
                "n_et_snap": np.array(snaps['n_et']),
                "n_ht_snap": np.array(snaps['n_ht']),
                "V_snap": np.array(snaps['V']),
                "E_snap": np.array(snaps['E']),
                "t_hist": np.array(hist['t']),
                "Q_emu": np.array(hist['Q_emu']),
                "Q_hmu": np.array(hist['Q_hmu']),
                "Q_et": np.array(hist['Q_et']),
                "Q_ht": np.array(hist['Q_ht']),
                "Va_hist": np.array(hist['Va']),
                "j_anode": np.array(hist['j_anode']),
                "j_cathode": np.array(hist['j_cathode'])})
    return res


# ============================================================================
# 10. FIGURES
# ============================================================================

def tracer(res, p, dossier="resultats"):
    import os
    os.makedirs(dossier, exist_ok=True)
    xc = res['xc']
    snaps = res['snaps']
    hist = res['hist']
    xmm = xc * 1e6

    # --- profils de charges à différents instants ---------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    labels = [('n_emu', "Électrons libres  $n_{e\\mu}$"),
              ('n_hmu', "Trous libres  $n_{h\\mu}$"),
              ('n_et', "Électrons piégés  $n_{et}$"),
              ('n_ht', "Trous piégés  $n_{ht}$")]
    for ax, (key, lab) in zip(axes.ravel(), labels):
        for k in range(len(snaps['t'])):
            ax.plot(xmm, snaps[key][k], lw=1.2,
                    label="t = %.2g s" % snaps['t'][k])
        ax.set_title(lab)
        ax.set_ylabel("densité de charge [C/m$^3$]")
        ax.grid(alpha=0.3)
    axes[1, 0].set_xlabel("x [µm]")
    axes[1, 1].set_xlabel("x [µm]")
    fig.suptitle("Profils de charges dans le diélectrique")
    fig.tight_layout()
    fig.savefig(os.path.join(dossier, "profils_charges.png"), dpi=130)
    plt.close(fig)

    # --- champ et potentiel --------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    xf = res['xf'] * 1e6
    for k in range(len(snaps['t'])):
        axes[0].plot(xf, snaps['V'][k], lw=1.2,
                     label="t = %.2g s" % snaps['t'][k])
        axes[1].plot(xmm, snaps['E'][k] * 1e-6, lw=1.2,
                     label="t = %.2g s" % snaps['t'][k])
    axes[0].set_title("Potentiel V(x,t)")
    axes[0].set_ylabel("V [volt]")
    axes[1].set_title("Champ électrique E(x,t)")
    axes[1].set_ylabel("E [MV/m]")
    for ax in axes:
        ax.set_xlabel("x [µm]")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(dossier, "potentiel_champ.png"), dpi=130)
    plt.close(fig)

    # --- évolutions temporelles globales --------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    th = np.array(hist['t'])
    axes[0].plot(th, np.array(hist['Q_emu']), label="$Q_{e\\mu}$")
    axes[0].plot(th, np.array(hist['Q_hmu']), label="$Q_{h\\mu}$")
    axes[0].plot(th, np.array(hist['Q_et']), label="$Q_{et}$")
    axes[0].plot(th, np.array(hist['Q_ht']), label="$Q_{ht}$")
    axes[0].set_title("Charges totales (intégrées sur l'épaisseur)")
    axes[0].set_ylabel("charge par unité de surface [C/m$^2$]")
    axes[0].legend()
    axes[1].plot(th, np.array(hist['Va']))
    axes[1].set_title("Potentiel d'anode $V_a(t)$")
    axes[1].set_ylabel("V [volt]")
    axes[2].semilogy(th[1:], np.abs(np.array(hist['j_anode'])[1:]) + 1e-30,
                     label="anode")
    axes[2].semilogy(th[1:], np.abs(np.array(hist['j_cathode'])[1:]) + 1e-30,
                     label="cathode")
    axes[2].set_title("Courants aux électrodes")
    axes[2].set_ylabel("|j| [A/m$^2$]")
    axes[2].legend()
    for ax in axes:
        ax.set_xlabel("temps [s]")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(dossier, "evolutions_temporelles.png"), dpi=130)
    plt.close(fig)

    # --- maillage -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.eventplot(res['xf'] * 1e6, lineoffsets=0.5, linelengths=0.4)
    ax.set_title("Maillage non uniforme (raffiné aux électrodes)")
    ax.set_xlabel("x [µm]")
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(dossier, "maillage.png"), dpi=130)
    plt.close(fig)


# ============================================================================
# 11. PROGRAMME PRINCIPAL
# ============================================================================

if __name__ == "__main__":
    p = Parametres()
    p.va_type = 'constante'   # choisir : 'constante' | 'rampe' | 'sinus'
    p.va0 = 1000.0             # potentiel d'anode [V]
    p.tmax = 50.0              # durée [s]

    print("Transport de charges 1D dans un diélectrique")
    print("  L = %.0f µm, N = %d cellules, Va(%s) = %g V"
          % (p.L * 1e6, p.N, p.va_type, p.va0))
    res = simuler(p)
    tracer(res, p)
    print("Résultats écrits dans 'resultats/' (npz + figures PNG).")
