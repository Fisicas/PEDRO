"""Look up every reference PEDRO's documentation relies on in Crossref and record what the registry says.

    python crossref_lookup_v01.py

For entries with a DOI the DOI is resolved; otherwise the citation string is searched (query.bibliographic) and the
top hit is recorded with its score, for a human to accept or reject. Writes crossref_results_v01.json (refuses to
overwrite). Only public bibliographic metadata is sent.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "crossref_results_v01.json"
UA = {"User-Agent": "pedro-reference-check/0.1 (bibliographic metadata lookup)"}

REFS = [  # (key, DOI or None, citation string used when no DOI is known)
    ("wilson1955", None, "Wilson Decius Cross Molecular Vibrations McGraw-Hill 1955"),
    ("morino1952", None, "Morino Kuchitsu 1952 J. Chem. Phys. 20 1809 potential energy distribution"),
    ("keresztury1971", None, "Keresztury Jalsovszky 1971 J. Mol. Struct. 10 304 alternative calculation of vibrational potential energy distribution"),
    ("jamroz2013", "10.1016/j.saa.2013.05.096", ""),
    ("jamroz2006", None, "Jamroz Dobrowolski 2006 Journal of Molecular Structure 787 172"),
    ("pulay1979", None, "Pulay Fogarasi Pang Boggs 1979 J. Am. Chem. Soc. 101 2550 systematic ab initio gradient calculation"),
    ("fogarasi1992", None, "Fogarasi Zhou Taylor Pulay 1992 J. Am. Chem. Soc. 114 8191 natural internal coordinates"),
    ("blondel1996", None, "Blondel Karplus 1996 J. Comput. Chem. 17 1132 new formulation for derivatives of torsion angles"),
    ("sinkhorn1964", None, "Sinkhorn 1964 Annals of Mathematical Statistics 35 876 relationship between arbitrary positive matrices and doubly stochastic matrices"),
    ("head1997", None, "Head 1997 International Journal of Quantum Chemistry 65 827 computation of vibrational frequencies for adsorbates on surfaces"),
    ("li2002", None, "Li Jensen 2002 Theoretical Chemistry Accounts 107 211 partial Hessian vibrational analysis"),
    ("cordero2008", None, "Cordero Gomez Platero-Prats 2008 Dalton Transactions 2832 covalent radii revisited"),
    ("tiesinga2021", None, "Tiesinga Mohr Newell Taylor 2021 Reviews of Modern Physics 93 025010 CODATA recommended values 2018"),
    ("barone2005", None, "Barone 2005 J. Chem. Phys. 122 014108 anharmonic vibrational properties by a fully automated second-order perturbative approach"),
    ("bloino2012", None, "Bloino Barone 2012 J. Chem. Phys. 136 124108 second-order perturbation theory approach for computing vibrational spectra"),
    ("martin1995", None, "Martin Lee Taylor Francois 1995 J. Chem. Phys. 103 2589 anharmonic force field of ethylene"),
    ("kabsch1976", None, "Kabsch 1976 Acta Crystallographica A32 922 solution for the best rotation to relate two sets of vectors"),
    ("kuhn1955", None, "Kuhn 1955 Naval Research Logistics Quarterly Hungarian method for the assignment problem"),
    ("weisfeiler1968", None, "Weisfeiler Leman 1968 reduction of a graph to a canonical form"),
    ("tao2018", None, "Tao Tian Verma Zou Wang Cremer Kraka 2018 Recovering intrinsic fragmental vibrations using the generalized subsystem vibrational analysis J. Chem. Theory Comput."),
    ("zou2021revgsva", None, "Zou Tao Freindorf Cremer Kraka 2020 local vibrational mode analysis revised generalized subsystem vibrational analysis"),
    ("moura2023", "10.1021/acs.jpca.2c07871", ""),
    ("schwinn2020", None, "Schwinn Ferre Huix-Rotllant 2020 Efficient analytic second derivative of electrostatic embedding QM/MM energy normal mode analysis of plant cryptochrome"),
    ("schwinn2021", None, "Schwinn Ferre Huix-Rotllant 2021 Physical Chemistry Chemical Physics QM/MM normal mode localized vibrational analysis"),
    ("ghysels2010tamkin", None, "Ghysels Verstraelen Hemelsoet Waroquier Van Speybroeck 2010 TAMkin toolkit for normal mode analysis J. Chem. Inf. Model."),
    ("ghysels2010qmmm", None, "Ghysels Woodcock Larkin Miller Shao Kong Van Neck Van Speybroeck Waroquier Brooks 2011 efficient calculation of QM/MM frequencies with the mobile block Hessian J. Chem. Theory Comput."),
    ("vester2024", "10.1021/acs.jctc.4c00882", ""),
    ("oenen2024", "10.1063/5.0180657", ""),
    ("neugebauer2009", "10.2533/chimia.2009.270", ""),
    ("reiher2003", None, "Reiher Neugebauer 2003 a mode-selective quantum chemical method for tracking molecular vibrations applied to functionalized carbon nanotubes J. Chem. Phys."),
    ("le2026", "10.1021/acsomega.6c06386", ""),
    ("flavin2025", "10.1039/d5cp02306h", ""),
    ("unimovib", None, "Zou Tao Freindorf Makos Verma Cremer Kraka UniMoVib unified interface molecular harmonic vibrational frequency calculations"),
    ("becke1993", None, "Becke 1993 J. Chem. Phys. 98 5648 density-functional thermochemistry III role of exact exchange"),
    ("lee1988", None, "Lee Yang Parr 1988 Phys. Rev. B 37 785 development of the Colle-Salvetti correlation-energy formula"),
    ("dunning1989", None, "Dunning 1989 J. Chem. Phys. 90 1007 Gaussian basis sets for use in correlated molecular calculations atoms boron through neon and hydrogen"),
    ("kendall1992", None, "Kendall Dunning Harrison 1992 J. Chem. Phys. 96 6796 electron affinities of the first-row atoms revisited"),
    ("weigend2005", None, "Weigend Ahlrichs 2005 Phys. Chem. Chem. Phys. 7 3297 balanced basis sets of split valence triple zeta valence and quadruple zeta valence quality"),
    ("grimme2011", None, "Grimme Ehrlich Goerigk 2011 J. Comput. Chem. 32 1456 effect of the damping function in dispersion corrected density functional theory"),
    ("neese2025orca6", None, "Neese 2025 software update the ORCA program system version 6.0 WIREs Computational Molecular Science"),
    ("rappe1992", None, "Rappe Casewit Colwell Goddard Skiff 1992 UFF a full periodic table force field J. Am. Chem. Soc. 114 10024"),
    ("dapprich1999", None, "Dapprich Komaromi Byun Morokuma Frisch 1999 new ONIOM implementation in Gaussian98 J. Mol. Struct. THEOCHEM"),
]


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


def summarise(m):
    au = m.get("author", [])
    return {"doi": m.get("DOI"), "title": (m.get("title") or [""])[0], "container": (m.get("container-title") or [""])[0],
            "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0], "volume": m.get("volume"),
            "page": m.get("page") or m.get("article-number"), "first_author": (au[0].get("family") if au else None),
            "authors": [f"{a.get('given', '')} {a.get('family', '')}".strip() for a in au]}


def main():
    if OUT.exists():
        sys.exit("refusing to overwrite " + OUT.name)
    res = {}
    for key, doi, cit in REFS:
        try:
            if doi:
                m = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
                res[key] = {"mode": "doi", **summarise(m)}
            else:
                q = urllib.parse.urlencode({"query.bibliographic": cit, "rows": 1})
                items = get("https://api.crossref.org/works?" + q)["message"]["items"]
                res[key] = {"mode": "search", "query": cit, "score": items[0].get("score") if items else None,
                            **(summarise(items[0]) if items else {})}
        except Exception as e:  # noqa: BLE001 - record and continue
            res[key] = {"error": str(e)}
        r = res[key]
        print(f"{key:20s} {r.get('mode', '?'):6s} {r.get('year')} {r.get('first_author')} | {r.get('title', r.get('error'))[:70]} | {r.get('container', '')[:30]} {r.get('volume')} {r.get('page')}", flush=True)
        time.sleep(0.3)
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
