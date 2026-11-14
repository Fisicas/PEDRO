# References

Every paper and program PEDRO's methods, documentation, validation and comparisons rely on. Each entry with a DOI was
checked against its Crossref record on 13 November 2026 (`docs/references_check/`); entries without a DOI are books,
software or reports that have none.

## Potential and total energy distributions; VEDA

1. E. B. Wilson Jr., J. C. Decius, P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational
   Spectra*, McGraw-Hill, New York, 1955. *GF method, internal coordinates and the B matrix (theory §1–3).*
2. Y. Morino, K. Kuchitsu, "A note on the classification of normal vibrations of molecules", *J. Chem. Phys.* **20**,
   1809–1810 (1952). [10.1063/1.1700299](https://doi.org/10.1063/1.1700299) *The diagonal potential energy distribution.*
3. G. Keresztury, G. Jalsovszky, "An alternative calculation of the vibrational potential energy distribution",
   *J. Mol. Struct.* **10**, 304–305 (1971). [10.1016/0022-2860(71)87090-4](https://doi.org/10.1016/0022-2860(71)87090-4)
   *The total energy distribution T = L ∘ L⁻ᵀ used by PEDRO.*
4. M. H. Jamróz, "Vibrational Energy Distribution Analysis (VEDA): Scopes and limitations", *Spectrochim. Acta A*
   **114**, 220–230 (2013). [10.1016/j.saa.2013.05.096](https://doi.org/10.1016/j.saa.2013.05.096) *VEDA 4; the
   methanol and 2-fluoropropane coordinate sets and PED values PEDRO is compared with.*
5. M. H. Jamróz, J. Cz. Dobrowolski, R. Brzozowski, "Vibrational modes of 2,6-, 2,7-, and 2,3-diisopropylnaphthalene.
   A DFT study", *J. Mol. Struct.* **787**, 172–183 (2006).
   [10.1016/j.molstruc.2005.10.044](https://doi.org/10.1016/j.molstruc.2005.10.044) *PED methodology cited in VEDA's
   documentation.*
6. M. H. Jamróz, *VEDA 4* (program, version 4.99e/f), with its example files and usage notes. *Reference program; its
   own example inputs and outputs are compared with PEDRO (not redistributed here).*

## Internal coordinates, numerics and molecular graphs

7. P. Pulay, G. Fogarasi, F. Pang, J. E. Boggs, "Systematic ab initio gradient calculation of molecular geometries,
   force constants, and dipole moment derivatives", *J. Am. Chem. Soc.* **101**, 2550–2560 (1979).
   [10.1021/ja00504a009](https://doi.org/10.1021/ja00504a009) *Local symmetry coordinates.*
8. G. Fogarasi, X. Zhou, P. W. Taylor, P. Pulay, "The calculation of ab initio molecular geometries: efficient
   optimization by natural internal coordinates and empirical correction by offset forces", *J. Am. Chem. Soc.*
   **114**, 8191–8201 (1992). [10.1021/ja00047a032](https://doi.org/10.1021/ja00047a032) *Natural internal coordinates
   (automatic sets).*
9. A. Blondel, M. Karplus, "New formulation for derivatives of torsion angles and improper torsion angles in
   molecular mechanics: Elimination of singularities", *J. Comput. Chem.* **17**, 1132–1141 (1996).
   [10.1002/(SICI)1096-987X(19960715)17:9<1132::AID-JCC5>3.0.CO;2-T](https://doi.org/10.1002/(SICI)1096-987X(19960715)17:9%3C1132::AID-JCC5%3E3.0.CO;2-T)
   *Analytic torsion B-matrix rows.*
10. R. Sinkhorn, "A relationship between arbitrary positive matrices and doubly stochastic matrices", *Ann. Math.
    Statist.* **35**, 876–879 (1964). [10.1214/aoms/1177703591](https://doi.org/10.1214/aoms/1177703591)
    *Balancing used for the VEDA-compatible PED.*
11. B. Cordero et al., "Covalent radii revisited", *Dalton Trans.* 2832–2838 (2008).
    [10.1039/b801115j](https://doi.org/10.1039/b801115j) *Bond detection.*
12. W. Kabsch, "A solution for the best rotation to relate two sets of vectors", *Acta Cryst. A* **32**, 922–923
    (1976). [10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873) *RMSD used to choose among graph
    isomorphisms.*
13. B. Weisfeiler, A. A. Leman, "The reduction of a graph to canonical form and the algebra which appears therein",
    *Nauchno-Tekhnicheskaya Informatsia, Ser. 2*, **9**, 12–16 (1968). *Color refinement in the isomorphism search.*
14. H. W. Kuhn, "The Hungarian method for the assignment problem", *Naval Res. Logist. Q.* **2**, 83–97 (1955).
    [10.1002/nav.3800020109](https://doi.org/10.1002/nav.3800020109) *Mode matching and GVPT2 mode mapping.*
15. E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor, "CODATA recommended values of the fundamental physical
    constants: 2018", *Rev. Mod. Phys.* **93**, 025010 (2021).
    [10.1103/RevModPhys.93.025010](https://doi.org/10.1103/RevModPhys.93.025010) *Physical constants.*
16. J. W. Ochterski, *Vibrational Analysis in Gaussian*, Gaussian, Inc. (1999). *Gaussian's conventions for normal
    modes and projections.*

## Partial Hessians and molecules in large models

17. J. D. Head, "Computation of vibrational frequencies for adsorbates on surfaces", *Int. J. Quantum Chem.* **65**,
    827–838 (1997). [10.1002/(SICI)1097-461X(1997)65:5<827::AID-QUA47>3.3.CO;2-L](https://doi.org/10.1002/(sici)1097-461x(1997)65:5%3C827::aid-qua47%3E3.3.co;2-l)
    *Partial Hessian vibrational analysis.*
18. H. Li, J. H. Jensen, "Partial Hessian vibrational analysis: the localization of the molecular vibrational energy
    and entropy", *Theor. Chem. Acc.* **107**, 211–219 (2002). [10.1007/s00214-001-0317-7](https://doi.org/10.1007/s00214-001-0317-7)
19. A. Ghysels, T. Verstraelen, K. Hemelsoet, M. Waroquier, V. Van Speybroeck, "TAMkin: A versatile package for
    vibrational analysis and chemical kinetics", *J. Chem. Inf. Model.* **50**, 1736–1750 (2010).
    [10.1021/ci100099g](https://doi.org/10.1021/ci100099g)
20. A. Ghysels et al., "Efficient calculation of QM/MM frequencies with the mobile block Hessian", *J. Chem. Theory
    Comput.* **7**, 496–514 (2011). [10.1021/ct100473f](https://doi.org/10.1021/ct100473f)
21. J. Vester, J. M. H. Olsen, "Assessing the partial Hessian approximation in QM/MM-based vibrational analysis",
    *J. Chem. Theory Comput.* **20**, 9533–9546 (2024). [10.1021/acs.jctc.4c00882](https://doi.org/10.1021/acs.jctc.4c00882)
    *Advice not to project rigid-body motion from partial Hessians (PEDRO case P).*
22. Y. Tao, C. Tian, N. Verma, W. Zou, C. Wang, D. Cremer, E. Kraka, "Recovering intrinsic fragmental vibrations using
    the generalized subsystem vibrational analysis", *J. Chem. Theory Comput.* **14**, 2558–2569 (2018).
    [10.1021/acs.jctc.7b01171](https://doi.org/10.1021/acs.jctc.7b01171) *GSVA.*
23. Y. Tao, W. Zou, S. Nanayakkara, M. Freindorf, E. Kraka, "A revised formulation of the generalized subsystem
    vibrational analysis (GSVA)", *Theor. Chem. Acc.* **140**, 31 (2021).
    [10.1007/s00214-021-02727-y](https://doi.org/10.1007/s00214-021-02727-y) *rev-GSVA, used as a cross-check.*
24. R. T. Moura Jr., M. Quintano, J. J. Antonio, M. Freindorf, E. Kraka, "Automatic generation of local vibrational
    mode parameters: from small to large molecules and QM/MM systems", *J. Phys. Chem. A* **126**, 9313–9331 (2022).
    [10.1021/acs.jpca.2c07871](https://doi.org/10.1021/acs.jpca.2c07871) *LModeAGen; closest prior QM/MM analysis.*
25. W. Zou et al., *UniMoVib* (program), https://github.com/zorkzou/UniMoVib. *Used for the GSVA and harmonic
    cross-checks.*
26. K. Schwinn, N. Ferré, M. Huix-Rotllant, "Efficient analytic second derivative of electrostatic embedding QM/MM
    energy: normal mode analysis of plant cryptochrome", *J. Chem. Theory Comput.* **16**, 3816–3824 (2020).
    [10.1021/acs.jctc.9b01145](https://doi.org/10.1021/acs.jctc.9b01145)
27. M. Huix-Rotllant, K. Schwinn, N. Ferré, "Infrared spectroscopy from electrostatic embedding QM/MM: local normal
    mode analysis of infrared spectra of arabidopsis thaliana plant cryptochrome", *Phys. Chem. Chem. Phys.* **23**,
    1666–1674 (2021). [10.1039/D0CP06070D](https://doi.org/10.1039/d0cp06070d)
28. M. Reiher, J. Neugebauer, "A mode-selective quantum chemical method for tracking molecular vibrations applied to
    functionalized carbon nanotubes", *J. Chem. Phys.* **118**, 1634–1641 (2003).
    [10.1063/1.1523908](https://doi.org/10.1063/1.1523908)
29. K. Kiewisch, S. Luber, J. Neugebauer, M. Reiher, "Intensity tracking for vibrational spectra of large
    molecules", *Chimia* **63**, 270–274 (2009). [10.2533/chimia.2009.270](https://doi.org/10.2533/chimia.2009.270)
30. K. Oenen, D. F. Dinu, K. R. Liedl, "Determining internal coordinate sets for optimal representation of molecular
    vibration", *J. Chem. Phys.* **160**, 014104 (2024). [10.1063/5.0180657](https://doi.org/10.1063/5.0180657) *NOMODECO.*

## Applications that motivated PEDRO

31. D. P. N. Le, G. Hastings, S. Gozem, "Assigning flavin's difference-FTIR spectral bands in solution: frequency and
    intensity shifts in flavin's 1-electron and 2-electron reduced states", *Phys. Chem. Chem. Phys.* **27**,
    24024–24039 (2025). [10.1039/D5CP02306H](https://doi.org/10.1039/d5cp02306h)
32. D. P. N. Le, G. Hastings, S. Gozem, "QM/MM calculated triplet–singlet infrared difference spectra of LOV2
    photoreceptor domains and FMN in solution", *ACS Omega* **11**, 51851–51863 (2026).
    [10.1021/acsomega.6c06386](https://doi.org/10.1021/acsomega.6c06386) *Reports that VEDA could not be applied to the
    QM/MM calculations.*

## Anharmonic analysis (GVPT2)

33. V. Barone, "Anharmonic vibrational properties by a fully automated second-order perturbative approach",
    *J. Chem. Phys.* **122**, 014108 (2005). [10.1063/1.1824881](https://doi.org/10.1063/1.1824881)
34. J. Bloino, V. Barone, "A second-order perturbation theory route to vibrational averages and transition
    properties of molecules: General formulation and application to infrared and vibrational circular dichroism
    spectroscopies", *J. Chem. Phys.* **136**, 124108 (2012). [10.1063/1.3695210](https://doi.org/10.1063/1.3695210)
35. J. M. L. Martin, T. J. Lee, P. R. Taylor, J.-P. François, "The anharmonic force field of ethylene, C₂H₄, by means
    of accurate ab initio calculations", *J. Chem. Phys.* **103**, 2589–2602 (1995).
    [10.1063/1.469681](https://doi.org/10.1063/1.469681) *Fermi-resonance test.*

## Programs and methods used for the example and validation calculations

36. M. J. Frisch et al., *Gaussian 16*, Revision A.03, Gaussian, Inc., Wallingford CT (2016).
37. F. Neese, "Software update: The ORCA program system—Version 6.0", *WIREs Comput. Mol. Sci.* **15**, e70019 (2025).
    [10.1002/wcms.70019](https://doi.org/10.1002/wcms.70019)
38. A. D. Becke, "Density-functional thermochemistry. III. The role of exact exchange", *J. Chem. Phys.* **98**,
    5648–5652 (1993). [10.1063/1.464913](https://doi.org/10.1063/1.464913)
39. C. Lee, W. Yang, R. G. Parr, "Development of the Colle-Salvetti correlation-energy formula into a functional of
    the electron density", *Phys. Rev. B* **37**, 785–789 (1988). [10.1103/PhysRevB.37.785](https://doi.org/10.1103/PhysRevB.37.785)
40. T. H. Dunning Jr., "Gaussian basis sets for use in correlated molecular calculations. I. The atoms boron through
    neon and hydrogen", *J. Chem. Phys.* **90**, 1007–1023 (1989). [10.1063/1.456153](https://doi.org/10.1063/1.456153)
41. R. A. Kendall, T. H. Dunning Jr., R. J. Harrison, "Electron affinities of the first-row atoms revisited.
    Systematic basis sets and wave functions", *J. Chem. Phys.* **96**, 6796–6806 (1992).
    [10.1063/1.462569](https://doi.org/10.1063/1.462569)
42. F. Weigend, R. Ahlrichs, "Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence
    quality for H to Rn: Design and assessment of accuracy", *Phys. Chem. Chem. Phys.* **7**, 3297–3305 (2005).
    [10.1039/B508541A](https://doi.org/10.1039/b508541a)
43. S. Grimme, S. Ehrlich, L. Goerigk, "Effect of the damping function in dispersion corrected density functional
    theory", *J. Comput. Chem.* **32**, 1456–1465 (2011). [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
44. A. K. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard III, W. M. Skiff, "UFF, a full periodic table force field
    for molecular mechanics and molecular dynamics simulations", *J. Am. Chem. Soc.* **114**, 10024–10035 (1992).
    [10.1021/ja00051a040](https://doi.org/10.1021/ja00051a040)
45. S. Dapprich, I. Komáromi, K. S. Byun, K. Morokuma, M. J. Frisch, "A new ONIOM implementation in Gaussian98. Part I.
    The calculation of energies, gradients, vibrational frequencies and electric field derivatives", *J. Mol.
    Struct. (THEOCHEM)* **461–462**, 1–21 (1999). [10.1016/S0166-1280(98)00475-8](https://doi.org/10.1016/s0166-1280(98)00475-8)
