"""Physical constants (CODATA 2018), element symbols, isotope masses and covalent radii."""

import math

HARTREE_J = 4.3597447222071e-18
BOHR_M = 5.29177210903e-11
BOHR_TO_ANG = 0.529177210903
ANG_TO_BOHR = 1.0 / BOHR_TO_ANG
AMU_KG = 1.66053906660e-27
C_CM_S = 2.99792458e10

# sqrt(Eh / (bohr^2 u)) / (2 pi c)  ->  cm^-1   (= 5140.4871...)
FREQ_FACTOR = math.sqrt(HARTREE_J / (BOHR_M**2 * AMU_KG)) / (2 * math.pi * C_CM_S)
# IR intensity (km/mol) per (d mu / d Q)^2 in e^2/u (Gaussian's convention)
IR_KMMOL_PER_E2_AMU = 974.8801

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe").split()

# Mass of the most abundant isotope (u), the Gaussian default. Used only when the input file carries no
# masses; every run records which masses were used and where they came from.
ISOTOPE_MASS = {
    "H": 1.00782503223, "He": 4.00260325413, "Li": 7.0160034366, "Be": 9.012183065, "B": 11.00930536,
    "C": 12.0, "N": 14.00307400443, "O": 15.99491461957, "F": 18.99840316273, "Ne": 19.9924401762,
    "Na": 22.989769282, "Mg": 23.985041697, "Al": 26.98153853, "Si": 27.97692653465, "P": 30.97376199842,
    "S": 31.9720711744, "Cl": 34.968852682, "Ar": 39.9623831237, "K": 38.9637064864, "Ca": 39.962590863,
    "Mn": 54.93804391, "Fe": 55.93493633, "Co": 58.93319429, "Ni": 57.93534241, "Cu": 62.92959772,
    "Zn": 63.92914201, "Se": 79.9165218, "Br": 78.9183376, "I": 126.9044719,
}

# Covalent radii (angstrom), Cordero et al., Dalton Trans. 2008, 2832.
COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03,
    "Ca": 1.76, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Se": 1.20,
    "Br": 1.20, "I": 1.39,
}
