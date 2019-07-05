#!/usr/bin/env python
# encoding: utf-8

name = "n-Heptane"
shortDesc = u"/gss_gpfs_scratch/westgroup/Importer/RMG-models/n-Heptane/nc7_ver3.1_mech.txt"
longDesc = u"""
n-Heptane, Detailed Mechanism, Version 3.1

This mechanism represents an updated release of the previous version 
available on the website (Version 3.0). The mechanism is based on the 
previously developed and very successful mechanism of Curran et al. 
1998 [1]. Version 3.1 fixes some bugs found in Version 3.0. This 
detailed chemical kinetic mechanism has been developed and validated 
by comparison to experiments in shock tubes and rapid compression 
machines. Over the series of experiments numerically investigated, 
the initial pressure ranged from 3 to 50 atm, the temperature from 
650 to 1200 K, and equivalence ratios from 0.3 to 1.0. The mechanism 
performs well at both low and high temperature and over a broad 
pressure range important for internal combustion engines.

References for Mechanism

Mehl M., W.J. Pitz, C.K. Westbrook, H.J. Curran, "Kinetic Modeling of Gasoline Surrogate Components and Mixtures Under Engine Conditions", Proceedings of the Combustion Institute 33:193-200 (2011).

M. Mehl, W. J. Pitz, M. Sj?berg and J. E. Dec, "Detailed kinetic modeling of low-temperature heat release for PRF fuels in an HCCI engine," SAE 2009 International Powertrains, Fuels and Lubricants Meeting, SAE Paper No. 2009-01-1806, Florence, Italy, 2009. Available at www.sae.org.

Other reference

[1] Curran, H. J., P. Gaffuri, W. J. Pitz, and C. K. Westbrook, "A Comprehensive Modeling Study of n-Heptane Oxidation" Combustion and Flame 114:149-177 (1998).


Downloaded from https://www-pls.llnl.gov/?url=science_and_technology-chemistry-combustion-n_heptane_version_3 in September 2013

Abstract
========================
Real fuels are complex mixtures of thousands of hydrocarbon compounds 
including linear and branched paraffins, naphthenes, olefins and 
aromatics. It is generally agreed that their behavior can be effectively 
reproduced by simpler fuel surrogates containing a limited number of components.

In this work, an improved version of the kinetic model by the 
authors is used to analyze the combustion behavior of several 
components relevant to gasoline surrogate formulation. Particular 
attention is devoted to linear and branched saturated hydrocarbons 
(PRF mixtures), olefins (1-hexene) and aromatics (toluene). Model 
predictions for pure components, binary mixtures and multi-component 
gasoline surrogates are compared with recent experimental information 
collected in rapid compression machine, shock tube and jet stirred 
reactors covering a wide range of conditions pertinent to internal 
combustion engines (3?50 atm, 650?1200 K, stoichiometric fuel/air mixtures). 
Simulation results are discussed focusing attention on the mixing effects of the fuel components.
"""
entry(
    index = 1,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.547e+15, 'cm^3/(mol*s)'),
        n = -0.406,
        Ea = (16600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> O + OH""",
)

entry(
    index = 2,
    label = "O + H2 <=> H + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(50800, 'cm^3/(mol*s)'), n=2.67, Ea=(6292, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O + H2 <=> H + OH""",
)

entry(
    index = 3,
    label = "OH + H2 <=> H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.16e+08, 'cm^3/(mol*s)'),
        n = 1.51,
        Ea = (3430, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + H2 <=> H + H2O""",
)

entry(
    index = 4,
    label = "O + H2O <=> OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.97e+06, 'cm^3/(mol*s)'),
        n = 2.02,
        Ea = (13400, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + H2O <=> OH + OH""",
)

entry(
    index = 5,
    label = "H2 <=> H + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (4.577e+19, 'cm^3/(mol*s)'),
            n = -1.4,
            Ea = (104400, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'[H][H]': 2.5, '[C]=O': 1.9, 'O=C=O': 3.8, 'O': 12},
    ),
    shortDesc = u"""The chemkin file reaction is H2 <=> H + H""",
)

entry(
    index = 6,
    label = "O2 <=> O + O",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (4.42e+17, 'cm^3/(mol*s)'),
            n = -0.634,
            Ea = (118900, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 3.8, 'CC': 3, 'O': 12, '[H][H]': 2.5, '[He]': 0.83, '[C]=O': 1.9, '[Ar]': 0.83},
    ),
    shortDesc = u"""The chemkin file reaction is O2 <=> O + O""",
)

entry(
    index = 7,
    label = "OH <=> O + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (9.78e+17, 'cm^3/(mol*s)'),
            n = -0.743,
            Ea = (102100, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 12, '[H][H]': 2.5, '[He]': 0.75, '[C]=O': 1.5, '[Ar]': 0.75},
    ),
    shortDesc = u"""The chemkin file reaction is OH <=> O + H""",
)

entry(
    index = 8,
    label = "H2O <=> H + OH",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.907e+23, 'cm^3/(mol*s)'),
            n = -1.83,
            Ea = (118500, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'CC': 3, 'O': 12, '[H][H]': 0.73, '[He]': 0.38, '[Ar]': 0.38},
    ),
    shortDesc = u"""The chemkin file reaction is H2O <=> H + OH""",
)

entry(
    index = 9,
    label = "H + O2 <=> HO2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.475e+12, 'cm^3/(mol*s)'), n=0.6, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.482e+16, 'cm^6/(mol^2*s)'),
            n = -0.411,
            Ea = (-1115, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 3.8, 'CC': 3, 'O': 14, '[H][H]': 1.3, '[He]': 0.67, '[C]=O': 1.9, '[Ar]': 0.67},
    ),
    shortDesc = u"""The chemkin file reaction is H + O2 <=> HO2""",
)

entry(
    index = 10,
    label = "HO2 + H <=> H2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.66e+13, 'cm^3/(mol*s)'), n=0, Ea=(823, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> H2 + O2""",
)

entry(
    index = 11,
    label = "HO2 + H <=> OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.079e+13, 'cm^3/(mol*s)'), n=0, Ea=(295, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + H <=> OH + OH""",
)

entry(
    index = 12,
    label = "HO2 + O <=> OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.25e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + O <=> OH + O2""",
)

entry(
    index = 13,
    label = "HO2 + OH <=> H2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.973e+10, 'cm^3/(mol*s)'),
        n = 0.962,
        Ea = (-328.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HO2 + OH <=> H2O + O2""",
)

entry(
    index = 14,
    label = "H2O2 + O2 <=> HO2 + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (1.136e+16, 'cm^3/(mol*s)'),
                n = -0.347,
                Ea = (49730, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(
                A = (2.141e+13, 'cm^3/(mol*s)'),
                n = -0.347,
                Ea = (37280, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H2O2 + O2 <=> HO2 + HO2""",
)

entry(
    index = 15,
    label = "H2O2 <=> OH + OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.951e+14, 's^-1'), n=0, Ea=(48430, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.202e+17, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (45500, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (1e-30, 'K'),
        T1 = (1e+30, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 3.8, 'CC': 3, 'O': 12, '[H][H]': 2.5, '[He]': 0.64, '[C]=O': 1.9, '[Ar]': 0.64},
    ),
    shortDesc = u"""The chemkin file reaction is H2O2 <=> OH + OH""",
)

entry(
    index = 16,
    label = "H2O2 + H <=> H2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(3970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + H <=> H2O + OH""",
)

entry(
    index = 17,
    label = "H2O2 + H <=> H2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.15e+10, 'cm^3/(mol*s)'), n=1, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + H <=> H2 + HO2""",
)

entry(
    index = 18,
    label = "H2O2 + O <=> OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.55e+06, 'cm^3/(mol*s)'), n=2, Ea=(3970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + O <=> OH + HO2""",
)

entry(
    index = 19,
    label = "H2O2 + OH <=> H2O + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(427.2, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.7e+18, 'cm^3/(mol*s)'), n=0, Ea=(29410, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is H2O2 + OH <=> H2O + HO2""",
)

entry(
    index = 20,
    label = "CO + O <=> CO2",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(1.8e+10, 'cm^3/(mol*s)'), n=0, Ea=(2384, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.35e+24, 'cm^6/(mol^2*s)'),
            n = -2.788,
            Ea = (4191, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 3.5, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.5, '[O][O]': 6, '[C]=O': 1.5, '[Ar]': 0.5},
    ),
    shortDesc = u"""The chemkin file reaction is CO + O <=> CO2""",
)

entry(
    index = 21,
    label = "CO + O2 <=> CO2 + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.05e+12, 'cm^3/(mol*s)'), n=0, Ea=(42540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO + O2 <=> CO2 + O""",
)

entry(
    index = 22,
    label = "CO + OH <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (223000, 'cm^3/(mol*s)'),
        n = 1.89,
        Ea = (-1158, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CO + OH <=> CO2 + H""",
)

entry(
    index = 23,
    label = "CO + HO2 <=> CO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(23000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CO + HO2 <=> CO2 + OH""",
)

entry(
    index = 24,
    label = "HCO <=> H + CO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (4.75e+11, 'cm^3/(mol*s)'),
            n = 0.66,
            Ea = (14870, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 12, '[H][H]': 2, '[C]=O': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is HCO <=> H + CO""",
)

entry(
    index = 25,
    label = "HCO + O2 <=> CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.58e+12, 'cm^3/(mol*s)'), n=0, Ea=(410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O2 <=> CO + HO2""",
)

entry(
    index = 26,
    label = "HCO + H <=> CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.34e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + H <=> CO + H2""",
)

entry(
    index = 27,
    label = "HCO + O <=> CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.02e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO + OH""",
)

entry(
    index = 28,
    label = "HCO + O <=> CO2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + O <=> CO2 + H""",
)

entry(
    index = 29,
    label = "HCO + OH <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.02e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + OH <=> CO + H2O""",
)

entry(
    index = 30,
    label = "HCO + CH3 <=> CH4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.65e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + CH3 <=> CH4 + CO""",
)

entry(
    index = 31,
    label = "HCO + HO2 <=> CH2O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.499e+14, 'cm^3/(mol*s)'),
        n = -0.061,
        Ea = (13920, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + HO2 <=> CH2O + O2""",
)

entry(
    index = 32,
    label = "HCO + HO2 => CO2 + H + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HO2 => CO2 + H + OH""",
)

entry(
    index = 33,
    label = "O2CHO <=> HCO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.959e+15, 's^-1'), n=-1.126, Ea=(41000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2CHO <=> HCO + O2""",
)

entry(
    index = 34,
    label = "CH2O + O2CHO <=> HCO + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + O2CHO <=> HCO + HO2CHO""",
)

entry(
    index = 35,
    label = "HO2CHO <=> OCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.01e+14, 's^-1'), n=0, Ea=(40150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2CHO <=> OCHO + OH""",
)

entry(
    index = 36,
    label = "OCHO <=> H + CO2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (5.318e+14, 'cm^3/(mol*s)'),
            n = -0.353,
            Ea = (17580, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is OCHO <=> H + CO2""",
)

entry(
    index = 37,
    label = "CH2O + CO <=> HCO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.186e+13, 'cm^3/(mol*s)'),
        n = 0.37,
        Ea = (73040, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + CO <=> HCO + HCO""",
)

entry(
    index = 38,
    label = "HCO + HCO => H2 + CO + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + HCO => H2 + CO + CO""",
)

entry(
    index = 39,
    label = "HCO + H <=> CH2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.09e+12, 'cm^3/(mol*s)'),
            n = 0.48,
            Ea = (-260, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.35e+24, 'cm^6/(mol^2*s)'),
            n = -2.57,
            Ea = (1425, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7824,
        T3 = (271, 'K'),
        T1 = (2755, 'K'),
        T2 = (6570, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is HCO + H <=> CH2O""",
)

entry(
    index = 40,
    label = "CO + H2 <=> CH2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (4.3e+07, 'cm^3/(mol*s)'),
            n = 1.5,
            Ea = (79600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (5.07e+27, 'cm^6/(mol^2*s)'),
            n = -3.42,
            Ea = (84348, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.932,
        T3 = (197, 'K'),
        T1 = (1540, 'K'),
        T2 = (10300, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CO + H2 <=> CH2O""",
)

entry(
    index = 41,
    label = "CH2O + OH <=> HCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.82e+07, 'cm^3/(mol*s)'),
        n = 1.63,
        Ea = (-1055, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + OH <=> HCO + H2O""",
)

entry(
    index = 42,
    label = "CH2O + H <=> HCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.74e+07, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (2740, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> HCO + H2""",
)

entry(
    index = 43,
    label = "CH2O + O <=> HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.26e+09, 'cm^3/(mol*s)'),
        n = 1.15,
        Ea = (2260, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + O <=> HCO + OH""",
)

entry(
    index = 44,
    label = "CH2O + CH3 <=> HCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(38.3, 'cm^3/(mol*s)'), n=3.36, Ea=(4312, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3 <=> HCO + CH4""",
)

entry(
    index = 45,
    label = "CH2O + HO2 <=> HCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.0071, 'cm^3/(mol*s)'),
        n = 4.517,
        Ea = (6580, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + HO2 <=> HCO + H2O2""",
)

entry(
    index = 46,
    label = "HOCH2O <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.056e+21, 's^-1'), n=-2.336, Ea=(25730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2O <=> CH2O + OH""",
)

entry(
    index = 47,
    label = "HOCH2O <=> HOCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2O <=> HOCHO + H""",
)

entry(
    index = 48,
    label = "HOCHO <=> CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.45e+12, 's^-1'), n=0, Ea=(60470, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO <=> CO + H2O""",
)

entry(
    index = 49,
    label = "HOCHO <=> CO2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.95e+09, 's^-1'), n=0, Ea=(48520, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO <=> CO2 + H2""",
)

entry(
    index = 50,
    label = "HOCHO <=> HCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.471e+22, 's^-1'), n=-1.542, Ea=(110700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO <=> HCO + OH""",
)

entry(
    index = 51,
    label = "HOCHO + O2 <=> OCHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.101e+12, 'cm^3/(mol*s)'),
        n = -0.308,
        Ea = (59880, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + O2 <=> OCHO + HO2""",
)

entry(
    index = 52,
    label = "HOCHO + OH => H2O + CO2 + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.62e+06, 'cm^3/(mol*s)'),
        n = 2.06,
        Ea = (916, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + OH => H2O + CO2 + H""",
)

entry(
    index = 53,
    label = "HOCHO + OH => H2O + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.85e+07, 'cm^3/(mol*s)'),
        n = 1.51,
        Ea = (-962, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + OH => H2O + CO + OH""",
)

entry(
    index = 54,
    label = "HOCHO + H => H2 + CO2 + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (4.24e+06, 'cm^3/(mol*s)'),
        n = 2.1,
        Ea = (4868, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + H => H2 + CO2 + H""",
)

entry(
    index = 55,
    label = "HOCHO + H => H2 + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (6.03e+13, 'cm^3/(mol*s)'),
        n = -0.35,
        Ea = (2988, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + H => H2 + CO + OH""",
)

entry(
    index = 56,
    label = "HOCHO + CH3 => CH4 + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.9e-07, 'cm^3/(mol*s)'), n=5.8, Ea=(2200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + CH3 => CH4 + CO + OH""",
)

entry(
    index = 57,
    label = "HOCHO + HO2 <=> OCHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.549e+12, 'cm^3/(mol*s)'),
        n = 0.04,
        Ea = (34470, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + HO2 <=> OCHO + H2O2""",
)

entry(
    index = 58,
    label = "HOCHO + HO2 => H2O2 + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCHO + HO2 => H2O2 + CO + OH""",
)

entry(
    index = 59,
    label = "HOCHO + O => CO + OH + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.77e+18, 'cm^3/(mol*s)'),
        n = -1.9,
        Ea = (2975, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + O => CO + OH + OH""",
)

entry(
    index = 60,
    label = "HOCHO + HCO <=> CH2O + OCHO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.584e+11, 'cm^3/(mol*s)'),
        n = 0.04,
        Ea = (26750, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HOCHO + HCO <=> CH2O + OCHO""",
)

entry(
    index = 61,
    label = "CH3O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.8e+13, 's^-1'), n=0, Ea=(26170, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.867e+25, 'cm^3/(mol*s)'),
            n = -3,
            Ea = (24307, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.9,
        T3 = (2500, 'K'),
        T1 = (1300, 'K'),
        T2 = (1e+99, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C]=O': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH3O <=> CH2O + H""",
)

entry(
    index = 62,
    label = "CH3O + O2 <=> CH2O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.38e-19, 'cm^3/(mol*s)'),
        n = 9.5,
        Ea = (-5501, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3O + O2 <=> CH2O + HO2""",
)

entry(
    index = 63,
    label = "CH2O + CH3O <=> CH3OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.62e+11, 'cm^3/(mol*s)'), n=0, Ea=(2294, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3O <=> CH3OH + HCO""",
)

entry(
    index = 64,
    label = "CH4 + CH3O <=> CH3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(611.9, 'cm^3/(mol*s)'), n=2.867, Ea=(8248, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH3O <=> CH3 + CH3OH""",
)

entry(
    index = 65,
    label = "CH3O + CH3 <=> CH2O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH3 <=> CH2O + CH4""",
)

entry(
    index = 66,
    label = "CH3O + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + H <=> CH2O + H2""",
)

entry(
    index = 67,
    label = "CH3O + HO2 <=> CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + HO2 <=> CH2O + H2O2""",
)

entry(
    index = 68,
    label = "CH2O + H <=> CH2OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.4e+11, 'cm^3/(mol*s)'),
            n = 0.454,
            Ea = (3600, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.27e+32, 'cm^6/(mol^2*s)'),
            n = -4.82,
            Ea = (6530, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.7187,
        T3 = (103, 'K'),
        T1 = (1291, 'K'),
        T2 = (4160, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[C]=O': 1.5},
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + H <=> CH2OH""",
)

entry(
    index = 69,
    label = "CH2OH + O2 <=> CH2O + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.51e+15, 'cm^3/(mol*s)'), n=-1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2.41e+14, 'cm^3/(mol*s)'), n=0, Ea=(5017, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2OH + O2 <=> CH2O + HO2""",
)

entry(
    index = 70,
    label = "CH2OH + H <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + H <=> CH2O + H2""",
)

entry(
    index = 71,
    label = "CH2OH + HO2 <=> CH2O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HO2 <=> CH2O + H2O2""",
)

entry(
    index = 72,
    label = "CH2OH + HCO <=> CH2O + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HCO <=> CH2O + CH2O""",
)

entry(
    index = 73,
    label = "CH2OH + CH3O <=> CH2O + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH3O <=> CH2O + CH3OH""",
)

entry(
    index = 74,
    label = "CH2OH + CH2O <=> CH3OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18780, 'cm^3/(mol*s)'), n=2.722, Ea=(4208, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + CH2O <=> CH3OH + HCO""",
)

entry(
    index = 75,
    label = "OH + CH2OH <=> H2O + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH + CH2OH <=> H2O + CH2O""",
)

entry(
    index = 76,
    label = "O + CH2OH <=> OH + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O + CH2OH <=> OH + CH2O""",
)

entry(
    index = 77,
    label = "CH2O + CH3OH <=> CH2OH + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.498e+12, 'cm^3/(mol*s)'),
        n = 0.659,
        Ea = (68460, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3OH <=> CH2OH + CH2OH""",
)

entry(
    index = 78,
    label = "CH2OH + HO2 <=> HOCH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OH + HO2 <=> HOCH2O + OH""",
)

entry(
    index = 79,
    label = "OCH2O2H <=> CH2O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.278e+18, 's^-1'), n=-1.8, Ea=(10460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCH2O2H <=> CH2O + HO2""",
)

entry(
    index = 80,
    label = "OCH2O2H <=> HOCH2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(8600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCH2O2H <=> HOCH2O2""",
)

entry(
    index = 81,
    label = "HOCH2O2 + HO2 <=> HOCH2O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2O2 + HO2 <=> HOCH2O2H + O2""",
)

entry(
    index = 82,
    label = "HOCH2O2H <=> HOCH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.023e+21, 's^-1'), n=-1.92, Ea=(42490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2O2H <=> HOCH2O + OH""",
)

entry(
    index = 83,
    label = "CH3OH <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.9e+16, 's^-1'), n=0, Ea=(91730, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.95e+44, 'cm^3/(mol*s)'),
            n = -7.35,
            Ea = (95460, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.414,
        T3 = (279, 'K'),
        T1 = (5459, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH <=> CH3 + OH""",
)

entry(
    index = 84,
    label = "CH3OH <=> CH2OH + H",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.69e+16, 's^-1'), n=-0.08, Ea=(98940, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.34e+40, 'cm^3/(mol*s)'),
            n = -6.33,
            Ea = (103100, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.773,
        T3 = (693, 'K'),
        T1 = (5333, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH <=> CH2OH + H""",
)

entry(
    index = 85,
    label = "CH3OH + H <=> CH3O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(6095, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH3O + H2""",
)

entry(
    index = 86,
    label = "CH3OH + H <=> CH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.204e+06, 'cm^3/(mol*s)'),
        n = 2.4,
        Ea = (2583, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + H <=> CH2OH + H2""",
)

entry(
    index = 87,
    label = "CH3OH + O <=> CH2OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(388000, 'cm^3/(mol*s)'), n=2.5, Ea=(3080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O <=> CH2OH + OH""",
)

entry(
    index = 88,
    label = "CH3OH + OH <=> CH3O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(513000, 'cm^3/(mol*s)'), n=2.13, Ea=(2450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + OH <=> CH3O + H2O""",
)

entry(
    index = 89,
    label = "CH3OH + OH <=> CH2OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.44e+06, 'cm^3/(mol*s)'), n=2, Ea=(-839, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + OH <=> CH2OH + H2O""",
)

entry(
    index = 90,
    label = "CH3OH + O2 <=> CH2OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.05e+13, 'cm^3/(mol*s)'), n=0, Ea=(44900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + O2 <=> CH2OH + HO2""",
)

entry(
    index = 91,
    label = "CH3OH + HO2 <=> CH2OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(10800, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + HO2 <=> CH2OH + H2O2""",
)

entry(
    index = 92,
    label = "CH3OH + CH3 <=> CH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(31.9, 'cm^3/(mol*s)'), n=3.17, Ea=(7172, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3 <=> CH2OH + CH4""",
)

entry(
    index = 93,
    label = "CH3O + CH3OH <=> CH2OH + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(4074, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + CH3OH <=> CH2OH + CH3OH""",
)

entry(
    index = 94,
    label = "CH3OH + CH2O <=> CH3O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.981e+12, 'cm^3/(mol*s)'),
        n = 0.452,
        Ea = (81490, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH2O <=> CH3O + CH3O""",
)

entry(
    index = 95,
    label = "CH3 + H <=> CH4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.27e+16, 'cm^3/(mol*s)'),
            n = -0.6,
            Ea = (383, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (2.477e+33, 'cm^6/(mol^2*s)'),
            n = -4.76,
            Ea = (2444, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.783,
        T3 = (74, 'K'),
        T1 = (2940, 'K'),
        T2 = (6960, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + H <=> CH4""",
)

entry(
    index = 96,
    label = "CH4 + H <=> CH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(614000, 'cm^3/(mol*s)'), n=2.5, Ea=(9587, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + H <=> CH3 + H2""",
)

entry(
    index = 97,
    label = "CH4 + OH <=> CH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(58300, 'cm^3/(mol*s)'), n=2.6, Ea=(2190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + OH <=> CH3 + H2O""",
)

entry(
    index = 98,
    label = "CH4 + O <=> CH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.02e+09, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (8600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH4 + O <=> CH3 + OH""",
)

entry(
    index = 99,
    label = "CH4 + HO2 <=> CH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(11.3, 'cm^3/(mol*s)'), n=3.74, Ea=(21010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + HO2 <=> CH3 + H2O2""",
)

entry(
    index = 100,
    label = "CH4 + CH2 <=> CH3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.46e+06, 'cm^3/(mol*s)'), n=2, Ea=(8270, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH2 <=> CH3 + CH3""",
)

entry(
    index = 101,
    label = "CH3 + OH <=> CH2O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+09, 'cm^3/(mol*s)'), n=0.5, Ea=(-1755, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2O + H2""",
)

entry(
    index = 102,
    label = "CH3 + OH <=> CH2(S) + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.508e+17, 'cm^3/(mol*s)'),
        n = -1.34,
        Ea = (1417, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2(S) + H2O""",
)

entry(
    index = 103,
    label = "CH3 + OH <=> CH3O + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.943e+07, 'cm^3/(mol*s)'),
        n = 1.343,
        Ea = (11200, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH3O + H""",
)

entry(
    index = 104,
    label = "CH3 + OH <=> CH2OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.09e+07, 'cm^3/(mol*s)'),
        n = 1.596,
        Ea = (4506, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2OH + H""",
)

entry(
    index = 105,
    label = "CH3 + OH <=> CH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+07, 'cm^3/(mol*s)'), n=1.6, Ea=(5420, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + OH <=> CH2 + H2O""",
)

entry(
    index = 106,
    label = "CH3 + HO2 <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1e+12, 'cm^3/(mol*s)'),
        n = 0.269,
        Ea = (-687.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + HO2 <=> CH3O + OH""",
)

entry(
    index = 107,
    label = "CH3 + HO2 <=> CH4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (116000, 'cm^3/(mol*s)'),
        n = 2.23,
        Ea = (-3022, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + HO2 <=> CH4 + O2""",
)

entry(
    index = 108,
    label = "CH3 + O <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.54e+13, 'cm^3/(mol*s)'),
        n = 0.05,
        Ea = (-136, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + O <=> CH2O + H""",
)

entry(
    index = 109,
    label = "CH3 + O2 <=> CH3O + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.546e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (28320, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH3O + O""",
)

entry(
    index = 110,
    label = "CH3 + O2 <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.641, 'cm^3/(mol*s)'), n=3.283, Ea=(8105, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH2O + OH""",
)

entry(
    index = 111,
    label = "CH3 + O2 <=> CH3O2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.812e+09, 'cm^3/(mol*s)'), n=0.9, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(6.85e+24, 'cm^6/(mol^2*s)'), n=-3, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.6,
        T3 = (1000, 'K'),
        T1 = (70, 'K'),
        T2 = (1700, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + O2 <=> CH3O2""",
)

entry(
    index = 112,
    label = "CH3O2 + CH2O <=> CH3O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + CH2O <=> CH3O2H + HCO""",
)

entry(
    index = 113,
    label = "CH4 + CH3O2 <=> CH3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+11, 'cm^3/(mol*s)'), n=0, Ea=(18480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH3O2 <=> CH3 + CH3O2H""",
)

entry(
    index = 114,
    label = "CH3OH + CH3O2 <=> CH2OH + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+12, 'cm^3/(mol*s)'), n=0, Ea=(13710, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + CH3O2 <=> CH2OH + CH3O2H""",
)

entry(
    index = 115,
    label = "CH3O2 + CH3 <=> CH3O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.08e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1411, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + CH3 <=> CH3O + CH3O""",
)

entry(
    index = 116,
    label = "CH3O2 + HO2 <=> CH3O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.47e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1570, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + HO2 <=> CH3O2H + O2""",
)

entry(
    index = 117,
    label = "CH3O2 + CH3O2 => CH2O + CH3OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (3.11e+14, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (-1051, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3O2 + CH3O2 => CH2O + CH3OH + O2""",
)

entry(
    index = 118,
    label = "CH3O2 + CH3O2 => O2 + CH3O + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3O2 + CH3O2 => O2 + CH3O + CH3O""",
)

entry(
    index = 119,
    label = "CH3O2 + H <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + H <=> CH3O + OH""",
)

entry(
    index = 120,
    label = "CH3O2 + O <=> CH3O + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + O <=> CH3O + O2""",
)

entry(
    index = 121,
    label = "CH3O2 + OH <=> CH3OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + OH <=> CH3OH + O2""",
)

entry(
    index = 122,
    label = "CH3O2H <=> CH3O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.31e+14, 's^-1'), n=0, Ea=(42300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2H <=> CH3O + OH""",
)

entry(
    index = 123,
    label = "CH2(S) <=> CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) <=> CH2""",
)

entry(
    index = 124,
    label = "CH2(S) + CH4 <=> CH3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 'cm^3/(mol*s)'), n=0, Ea=(-570, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + CH4 <=> CH3 + CH3""",
)

entry(
    index = 125,
    label = "CH2(S) + O2 => CO + OH + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + O2 => CO + OH + H""",
)

entry(
    index = 126,
    label = "CH2(S) + H2 <=> CH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H2 <=> CH3 + H""",
)

entry(
    index = 127,
    label = "CH2(S) + H <=> CH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H <=> CH2 + H""",
)

entry(
    index = 128,
    label = "CH2(S) + H <=> CH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + H <=> CH + H2""",
)

entry(
    index = 129,
    label = "CH2(S) + O => CO + H + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + O => CO + H + H""",
)

entry(
    index = 130,
    label = "CH2(S) + OH <=> CH2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + OH <=> CH2O + H""",
)

entry(
    index = 131,
    label = "CH2(S) + CO2 <=> CH2O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + CO2 <=> CH2O + CO""",
)

entry(
    index = 132,
    label = "CH2 + H <=> CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.5e+16, 'cm^3/(mol*s)'), n=-0.8, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.2e+27, 'cm^6/(mol^2*s)'),
            n = -3.14,
            Ea = (1230, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.68,
        T3 = (78, 'K'),
        T1 = (1995, 'K'),
        T2 = (5590, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H <=> CH3""",
)

entry(
    index = 133,
    label = "CH2 + O2 <=> CH2O + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 <=> CH2O + O""",
)

entry(
    index = 134,
    label = "CH2 + O2 => CO2 + H + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 => CO2 + H + H""",
)

entry(
    index = 135,
    label = "CH2 + O2 => CO + OH + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O2 => CO + OH + H""",
)

entry(
    index = 136,
    label = "CH2 + O => CO + H + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + O => CO + H + H""",
)

entry(
    index = 137,
    label = "CH2 + H <=> CH + H2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+18, 'cm^3/(mol*s)'), n=-1.56, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(
                A = (2.7e+11, 'cm^3/(mol*s)'),
                n = 0.67,
                Ea = (25700, 'cal/mol'),
                T0 = (1, 'K'),
            ),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + H <=> CH + H2""",
)

entry(
    index = 138,
    label = "CH2 + OH <=> CH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+07, 'cm^3/(mol*s)'), n=2, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2 + OH <=> CH + H2O""",
)

entry(
    index = 139,
    label = "CH + O2 <=> HCO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O2 <=> HCO + O""",
)

entry(
    index = 140,
    label = "C + OH <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + OH <=> CO + H""",
)

entry(
    index = 141,
    label = "C + O2 <=> CO + O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C + O2 <=> CO + O""",
)

entry(
    index = 142,
    label = "CH + H <=> C + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H <=> C + H2""",
)

entry(
    index = 143,
    label = "CH + O <=> CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + O <=> CO + H""",
)

entry(
    index = 144,
    label = "CH + OH <=> HCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + OH <=> HCO + H""",
)

entry(
    index = 145,
    label = "CH + H2O <=> H + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.713e+13, 'cm^3/(mol*s)'), n=0, Ea=(-755, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + H2O <=> H + CH2O""",
)

entry(
    index = 146,
    label = "CH + CO2 <=> HCO + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(685, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + CO2 <=> HCO + CO""",
)

entry(
    index = 147,
    label = "CH3 + CH3 <=> C2H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (9.214e+16, 'cm^3/(mol*s)'),
            n = -1.17,
            Ea = (635.8, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.135e+36, 'cm^6/(mol^2*s)'),
            n = -5.246,
            Ea = (1705, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.405,
        T3 = (1120, 'K'),
        T1 = (69.6, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH3 + CH3 <=> C2H6""",
)

entry(
    index = 148,
    label = "C2H5 + H <=> C2H6",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (5.21e+17, 'cm^3/(mol*s)'),
            n = -0.99,
            Ea = (1580, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.99e+41, 'cm^6/(mol^2*s)'),
            n = -7.08,
            Ea = (6685, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.842,
        T3 = (125, 'K'),
        T1 = (2219, 'K'),
        T2 = (6882, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + H <=> C2H6""",
)

entry(
    index = 149,
    label = "C2H6 + H <=> C2H5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.15e+08, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (7530, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + H <=> C2H5 + H2""",
)

entry(
    index = 150,
    label = "C2H6 + O <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.55e+06, 'cm^3/(mol*s)'),
        n = 2.4,
        Ea = (5830, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H6 + O <=> C2H5 + OH""",
)

entry(
    index = 151,
    label = "C2H6 + OH <=> C2H5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.48e+07, 'cm^3/(mol*s)'), n=1.9, Ea=(950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + OH <=> C2H5 + H2O""",
)

entry(
    index = 152,
    label = "C2H6 + O2 <=> C2H5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+13, 'cm^3/(mol*s)'), n=0, Ea=(51870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + O2 <=> C2H5 + HO2""",
)

entry(
    index = 153,
    label = "C2H6 + CH3 <=> C2H5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51e-07, 'cm^3/(mol*s)'), n=6, Ea=(6047, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3 <=> C2H5 + CH4""",
)

entry(
    index = 154,
    label = "C2H6 + HO2 <=> C2H5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(34.6, 'cm^3/(mol*s)'), n=3.61, Ea=(16920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + HO2 <=> C2H5 + H2O2""",
)

entry(
    index = 155,
    label = "C2H6 + CH3O2 <=> C2H5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19.4, 'cm^3/(mol*s)'), n=3.64, Ea=(17100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3O2 <=> C2H5 + CH3O2H""",
)

entry(
    index = 156,
    label = "C2H6 + CH3O <=> C2H5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+11, 'cm^3/(mol*s)'), n=0, Ea=(7090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3O <=> C2H5 + CH3OH""",
)

entry(
    index = 157,
    label = "C2H6 + CH <=> C2H5 + CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(-260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH <=> C2H5 + CH2""",
)

entry(
    index = 158,
    label = "CH2(S) + C2H6 <=> CH3 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + C2H6 <=> CH3 + C2H5""",
)

entry(
    index = 159,
    label = "C2H4 + H <=> C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.081e+12, 'cm^3/(mol*s)'),
            n = 0.454,
            Ea = (1822, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.2e+42, 'cm^6/(mol^2*s)'),
            n = -7.62,
            Ea = (6970, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.975,
        T3 = (210, 'K'),
        T1 = (984, 'K'),
        T2 = (4374, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H5""",
)

entry(
    index = 160,
    label = "H2 + CH3O2 <=> H + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + CH3O2 <=> H + CH3O2H""",
)

entry(
    index = 161,
    label = "H2 + C2H5O2 <=> H + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+14, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + C2H5O2 <=> H + C2H5O2H""",
)

entry(
    index = 162,
    label = "C2H5 + C2H3 <=> C2H4 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.859e+11, 'cm^3/(mol*s)'),
        n = 0.11,
        Ea = (-4300, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + C2H3 <=> C2H4 + C2H4""",
)

entry(
    index = 163,
    label = "CH3 + C2H5 <=> CH4 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(11800, 'cm^3/(mol*s)'), n=2.45, Ea=(-2921, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C2H5 <=> CH4 + C2H4""",
)

entry(
    index = 164,
    label = "C2H5 + H <=> CH3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.69e+13, 'cm^3/(mol*s)'), n=0, Ea=(220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + H <=> CH3 + CH3""",
)

entry(
    index = 165,
    label = "C2H5 + H <=> C2H4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + H <=> C2H4 + H2""",
)

entry(
    index = 166,
    label = "C2H5 + O <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O <=> CH3CHO + H""",
)

entry(
    index = 167,
    label = "C2H5 + HO2 <=> C2H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + HO2 <=> C2H5O + OH""",
)

entry(
    index = 168,
    label = "CH3O2 + C2H5 <=> CH3O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + C2H5 <=> CH3O + C2H5O""",
)

entry(
    index = 169,
    label = "C2H5O + O2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.28e+10, 'cm^3/(mol*s)'), n=0, Ea=(1097, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O + O2 <=> CH3CHO + HO2""",
)

entry(
    index = 170,
    label = "C2H5O <=> CH3 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.321e+20, 's^-1'), n=-2.018, Ea=(20750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O <=> CH3 + CH2O""",
)

entry(
    index = 171,
    label = "C2H5O <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.428e+15, 's^-1'), n=-0.687, Ea=(22230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O <=> CH3CHO + H""",
)

entry(
    index = 172,
    label = "C2H5O2 <=> C2H5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.312e+62, 's^-1'), n=-14.784, Ea=(49180, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 <=> C2H5 + O2""",
)

entry(
    index = 173,
    label = "C2H5O2 + CH2O <=> C2H5O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 + CH2O <=> C2H5O2H + HCO""",
)

entry(
    index = 174,
    label = "CH4 + C2H5O2 <=> CH3 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+11, 'cm^3/(mol*s)'), n=0, Ea=(18480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + C2H5O2 <=> CH3 + C2H5O2H""",
)

entry(
    index = 175,
    label = "CH3OH + C2H5O2 <=> CH2OH + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+12, 'cm^3/(mol*s)'), n=0, Ea=(13710, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + C2H5O2 <=> CH2OH + C2H5O2H""",
)

entry(
    index = 176,
    label = "C2H5O2 + HO2 <=> C2H5O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 + HO2 <=> C2H5O2H + O2""",
)

entry(
    index = 177,
    label = "C2H6 + C2H5O2 <=> C2H5 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.6, 'cm^3/(mol*s)'), n=3.76, Ea=(17200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + C2H5O2 <=> C2H5 + C2H5O2H""",
)

entry(
    index = 178,
    label = "C2H5O2H <=> C2H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.31e+14, 's^-1'), n=0, Ea=(42300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2H <=> C2H5O + OH""",
)

entry(
    index = 179,
    label = "C2H4O2H <=> C2H5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.374e+47, 's^-1'), n=-12.115, Ea=(31020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O2H <=> C2H5 + O2""",
)

entry(
    index = 180,
    label = "C2H5 + O2 <=> C2H4 + HO2",
    degeneracy = 1,
    duplicate = True,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(
                A = (7.561e+14, 'cm^3/(mol*s)'),
                n = -1.01,
                Ea = (4749, 'cal/mol'),
                T0 = (1, 'K'),
            ),
            Arrhenius(A=(0.4, 'cm^3/(mol*s)'), n=3.88, Ea=(13620, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> C2H4 + HO2""",
)

entry(
    index = 181,
    label = "C2H5 + O2 <=> C2H4O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.626e+11, 'cm^3/(mol*s)'),
        n = -0.31,
        Ea = (6150, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> C2H4O1-2 + OH""",
)

entry(
    index = 182,
    label = "C2H5 + O2 <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(826.5, 'cm^3/(mol*s)'), n=2.41, Ea=(5285, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + O2 <=> CH3CHO + OH""",
)

entry(
    index = 183,
    label = "C2H5O2 <=> C2H4O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.276e+39, 's^-1'), n=-8.479, Ea=(45170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 <=> C2H4O2H""",
)

entry(
    index = 184,
    label = "C2H5O2 <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.52e+41, 's^-1'), n=-10.2, Ea=(43710, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 <=> CH3CHO + OH""",
)

entry(
    index = 185,
    label = "C2H5O2 <=> C2H4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.815e+38, 's^-1'), n=-8.45, Ea=(37890, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 <=> C2H4 + HO2""",
)

entry(
    index = 186,
    label = "C2H5O2 <=> C2H4O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+43, 's^-1'), n=-10.46, Ea=(45580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 <=> C2H4O1-2 + OH""",
)

entry(
    index = 187,
    label = "C2H4O2H <=> C2H4O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.848e+30, 's^-1'), n=-6.08, Ea=(20660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O2H <=> C2H4O1-2 + OH""",
)

entry(
    index = 188,
    label = "C2H4O2H <=> C2H4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+34, 's^-1'), n=-7.25, Ea=(23250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O2H <=> C2H4 + HO2""",
)

entry(
    index = 189,
    label = "C2H4O2H <=> CH3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.188e+34, 's^-1'), n=-9.02, Ea=(29210, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O2H <=> CH3CHO + OH""",
)

entry(
    index = 190,
    label = "C2H4O1-2 <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.63e+13, 's^-1'), n=0, Ea=(57200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 <=> CH3 + HCO""",
)

entry(
    index = 191,
    label = "C2H4O1-2 <=> CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.407e+12, 's^-1'), n=0, Ea=(53800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 <=> CH3CHO""",
)

entry(
    index = 192,
    label = "C2H4O1-2 + OH <=> C2H3O1-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.78e+13, 'cm^3/(mol*s)'), n=0, Ea=(3610, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + OH <=> C2H3O1-2 + H2O""",
)

entry(
    index = 193,
    label = "C2H4O1-2 + H <=> C2H3O1-2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(9680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + H <=> C2H3O1-2 + H2""",
)

entry(
    index = 194,
    label = "C2H4O1-2 + HO2 <=> C2H3O1-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + HO2 <=> C2H3O1-2 + H2O2""",
)

entry(
    index = 195,
    label = "C2H4O1-2 + CH3O2 <=> C2H3O1-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + CH3O2 <=> C2H3O1-2 + CH3O2H""",
)

entry(
    index = 196,
    label = "C2H4O1-2 + C2H5O2 <=> C2H3O1-2 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + C2H5O2 <=> C2H3O1-2 + C2H5O2H""",
)

entry(
    index = 197,
    label = "C2H4O1-2 + CH3 <=> C2H3O1-2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07e+12, 'cm^3/(mol*s)'), n=0, Ea=(11830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + CH3 <=> C2H3O1-2 + CH4""",
)

entry(
    index = 198,
    label = "C2H4O1-2 + CH3O <=> C2H3O1-2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(6750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4O1-2 + CH3O <=> C2H3O1-2 + CH3OH""",
)

entry(
    index = 199,
    label = "C2H3O1-2 <=> CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.5e+14, 's^-1'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3O1-2 <=> CH3CO""",
)

entry(
    index = 200,
    label = "C2H3O1-2 <=> CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3O1-2 <=> CH2CHO""",
)

entry(
    index = 201,
    label = "CH3CHO <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.687e+20, 's^-1'), n=-1.342, Ea=(86950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO <=> CH3 + HCO""",
)

entry(
    index = 202,
    label = "CH3CHO + H <=> CH3CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.37e+13, 'cm^3/(mol*s)'), n=0, Ea=(3642, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + H <=> CH3CO + H2""",
)

entry(
    index = 203,
    label = "CH3CHO + O <=> CH3CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.94e+12, 'cm^3/(mol*s)'), n=0, Ea=(1868, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O <=> CH3CO + OH""",
)

entry(
    index = 204,
    label = "CH3CHO + OH <=> CH3CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.37e+12, 'cm^3/(mol*s)'), n=0, Ea=(-619, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH3CO + H2O""",
)

entry(
    index = 205,
    label = "CH3CHO + O2 <=> CH3CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(39150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + O2 <=> CH3CO + HO2""",
)

entry(
    index = 206,
    label = "CH3CHO + CH3 <=> CH3CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.000708, 'cm^3/(mol*s)'),
        n = 4.58,
        Ea = (1966, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3 <=> CH3CO + CH4""",
)

entry(
    index = 207,
    label = "CH3CHO + HO2 <=> CH3CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + HO2 <=> CH3CO + H2O2""",
)

entry(
    index = 208,
    label = "CH3O2 + CH3CHO <=> CH3O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + CH3CHO <=> CH3O2H + CH3CO""",
)

entry(
    index = 209,
    label = "CH3CHO + CH3CO3 <=> CH3CO + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + CH3CO3 <=> CH3CO + CH3CO3H""",
)

entry(
    index = 210,
    label = "CH3CHO + OH <=> CH3 + HOCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+15, 'cm^3/(mol*s)'), n=-1.076, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH3 + HOCHO""",
)

entry(
    index = 211,
    label = "CH3CHO + OH <=> CH2CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(172000, 'cm^3/(mol*s)'), n=2.4, Ea=(815, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHO + OH <=> CH2CHO + H2O""",
)

entry(
    index = 212,
    label = "CH3CO <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Lindemann(
        arrheniusHigh = Arrhenius(A=(3e+12, 's^-1'), n=0, Ea=(16720, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(1.2e+15, 'cm^3/(mol*s)'), n=0, Ea=(12518, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO <=> CH3 + CO""",
)

entry(
    index = 213,
    label = "CH3CO + H <=> CH2CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + H <=> CH2CO + H2""",
)

entry(
    index = 214,
    label = "CH3CO + O <=> CH2CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + O <=> CH2CO + OH""",
)

entry(
    index = 215,
    label = "CH3CO + CH3 <=> CH2CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO + CH3 <=> CH2CO + CH4""",
)

entry(
    index = 216,
    label = "CH3CO3 <=> CH3CO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.863e+19, 's^-1'), n=-1.949, Ea=(38530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO3 <=> CH3CO + O2""",
)

entry(
    index = 217,
    label = "CH3CO3 + HO2 <=> CH3CO3H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO3 + HO2 <=> CH3CO3H + O2""",
)

entry(
    index = 218,
    label = "H2O2 + CH3CO3 <=> HO2 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+12, 'cm^3/(mol*s)'), n=0, Ea=(9936, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + CH3CO3 <=> HO2 + CH3CO3H""",
)

entry(
    index = 219,
    label = "CH4 + CH3CO3 <=> CH3 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+11, 'cm^3/(mol*s)'), n=0, Ea=(18480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + CH3CO3 <=> CH3 + CH3CO3H""",
)

entry(
    index = 220,
    label = "CH2O + CH3CO3 <=> HCO + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3CO3 <=> HCO + CH3CO3H""",
)

entry(
    index = 221,
    label = "C2H6 + CH3CO3 <=> C2H5 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + CH3CO3 <=> C2H5 + CH3CO3H""",
)

entry(
    index = 222,
    label = "CH3CO3H <=> CH3CO2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.01e+14, 's^-1'), n=0, Ea=(40150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CO3H <=> CH3CO2 + OH""",
)

entry(
    index = 223,
    label = "CH3CO2 <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(4.4e+15, 'cm^3/(mol*s)'), n=0, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3CO2 <=> CH3 + CO2""",
)

entry(
    index = 224,
    label = "CH2CHO <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.071e+15, 's^-1'), n=-0.342, Ea=(50600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHO <=> CH2CO + H""",
)

entry(
    index = 225,
    label = "CH2CHO + O2 => CH2O + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (8.95e+13, 'cm^3/(mol*s)'),
        n = -0.6,
        Ea = (10120, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH2CHO + O2 => CH2O + CO + OH""",
)

entry(
    index = 226,
    label = "CH2 + CO <=> CH2CO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (2.69e+33, 'cm^6/(mol^2*s)'),
            n = -5.11,
            Ea = (7095, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5907,
        T3 = (275, 'K'),
        T1 = (1226, 'K'),
        T2 = (5185, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is CH2 + CO <=> CH2CO""",
)

entry(
    index = 227,
    label = "CH2CO + H <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(3400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> CH3 + CO""",
)

entry(
    index = 228,
    label = "CH2CO + H <=> HCCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+14, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + H <=> HCCO + H2""",
)

entry(
    index = 229,
    label = "CH2CO + O <=> CH2 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+12, 'cm^3/(mol*s)'), n=0, Ea=(1350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + O <=> CH2 + CO2""",
)

entry(
    index = 230,
    label = "CH2CO + O <=> HCCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + O <=> HCCO + OH""",
)

entry(
    index = 231,
    label = "CH2CO + OH <=> HCCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> HCCO + H2O""",
)

entry(
    index = 232,
    label = "CH2CO + OH <=> CH2OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CO + OH <=> CH2OH + CO""",
)

entry(
    index = 233,
    label = "CH2(S) + CH2CO <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + CH2CO <=> C2H4 + CO""",
)

entry(
    index = 234,
    label = "HCCO + OH => H2 + CO + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + OH => H2 + CO + CO""",
)

entry(
    index = 235,
    label = "H + HCCO <=> CH2(S) + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + HCCO <=> CH2(S) + CO""",
)

entry(
    index = 236,
    label = "HCCO + O => H + CO + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O => H + CO + CO""",
)

entry(
    index = 237,
    label = "HCCO + O2 => OH + CO + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4.2e+10, 'cm^3/(mol*s)'), n=0, Ea=(850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCCO + O2 => OH + CO + CO""",
)

entry(
    index = 238,
    label = "HCCO <=> CH + CO",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(6.5e+15, 'cm^3/(mol*s)'), n=0, Ea=(58820, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is HCCO <=> CH + CO""",
)

entry(
    index = 239,
    label = "CH + CH2O <=> H + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.46e+13, 'cm^3/(mol*s)'), n=0, Ea=(-515, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + CH2O <=> H + CH2CO""",
)

entry(
    index = 240,
    label = "CH + HCCO <=> CO + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + HCCO <=> CO + C2H2""",
)

entry(
    index = 241,
    label = "C2H3 + H <=> C2H4",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(
            A = (1.36e+14, 'cm^3/(mol*s)'),
            n = 0.173,
            Ea = (660, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        arrheniusLow = Arrhenius(
            A = (1.4e+30, 'cm^6/(mol^2*s)'),
            n = -3.86,
            Ea = (3320, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.782,
        T3 = (207.5, 'K'),
        T1 = (2663, 'K'),
        T2 = (6095, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H4""",
)

entry(
    index = 242,
    label = "C2H4 <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(8e+12, 's^-1'), n=0.44, Ea=(88770, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (1.58e+51, 'cm^3/(mol*s)'),
            n = -9.3,
            Ea = (97800, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.735,
        T3 = (180, 'K'),
        T1 = (1035, 'K'),
        T2 = (5417, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 <=> C2H2 + H2""",
)

entry(
    index = 243,
    label = "C2H4 + H <=> C2H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.07e+07, 'cm^3/(mol*s)'),
        n = 1.93,
        Ea = (12950, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + H <=> C2H3 + H2""",
)

entry(
    index = 244,
    label = "C2H4 + O <=> CH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.564e+06, 'cm^3/(mol*s)'),
        n = 1.88,
        Ea = (183, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH3 + HCO""",
)

entry(
    index = 245,
    label = "C2H4 + O <=> CH2CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.986e+06, 'cm^3/(mol*s)'),
        n = 1.88,
        Ea = (183, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H4 + O <=> CH2CHO + H""",
)

entry(
    index = 246,
    label = "C2H4 + OH <=> C2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+06, 'cm^3/(mol*s)'), n=2, Ea=(2500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + OH <=> C2H3 + H2O""",
)

entry(
    index = 247,
    label = "C2H4 + CH3 <=> C2H3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.62, 'cm^3/(mol*s)'), n=3.7, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3 <=> C2H3 + CH4""",
)

entry(
    index = 248,
    label = "C2H4 + O2 <=> C2H3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(58200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + O2 <=> C2H3 + HO2""",
)

entry(
    index = 249,
    label = "C2H4 + CH3O <=> C2H3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(6750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3O <=> C2H3 + CH3OH""",
)

entry(
    index = 250,
    label = "C2H4 + CH3O2 <=> C2H3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.23e+12, 'cm^3/(mol*s)'), n=0, Ea=(17190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3O2 <=> C2H3 + CH3O2H""",
)

entry(
    index = 251,
    label = "C2H4 + C2H5O2 <=> C2H3 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.23e+12, 'cm^3/(mol*s)'), n=0, Ea=(17190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + C2H5O2 <=> C2H3 + C2H5O2H""",
)

entry(
    index = 252,
    label = "C2H4 + CH3CO3 <=> C2H3 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3CO3 <=> C2H3 + CH3CO3H""",
)

entry(
    index = 253,
    label = "C2H4 + CH3O2 <=> C2H4O1-2 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.82e+12, 'cm^3/(mol*s)'), n=0, Ea=(17110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + CH3O2 <=> C2H4O1-2 + CH3O""",
)

entry(
    index = 254,
    label = "C2H4 + C2H5O2 <=> C2H4O1-2 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.82e+12, 'cm^3/(mol*s)'), n=0, Ea=(17110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + C2H5O2 <=> C2H4O1-2 + C2H5O""",
)

entry(
    index = 255,
    label = "C2H4 + HO2 <=> C2H4O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.23e+12, 'cm^3/(mol*s)'), n=0, Ea=(17190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + HO2 <=> C2H4O1-2 + OH""",
)

entry(
    index = 256,
    label = "CH + CH4 <=> C2H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH + CH4 <=> C2H4 + H""",
)

entry(
    index = 257,
    label = "CH2(S) + CH3 <=> C2H4 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2(S) + CH3 <=> C2H4 + H""",
)

entry(
    index = 258,
    label = "C2H2 + H <=> C2H3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(2400, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.8e+40, 'cm^6/(mol^2*s)'),
            n = -7.27,
            Ea = (7220, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.751,
        T3 = (98.5, 'K'),
        T1 = (1302, 'K'),
        T2 = (4167, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + H <=> C2H3""",
)

entry(
    index = 259,
    label = "C2H3 + O2 <=> C2H2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.12e-06, 'cm^3/(mol*s)'), n=6, Ea=(9484, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> C2H2 + HO2""",
)

entry(
    index = 260,
    label = "C2H3 + O2 <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.5e+28, 'cm^3/(mol*s)'),
        n = -5.312,
        Ea = (6500, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2O + HCO""",
)

entry(
    index = 261,
    label = "C2H3 + O2 <=> CH2CHO + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.5e+14, 'cm^3/(mol*s)'),
        n = -0.611,
        Ea = (5260, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3 + O2 <=> CH2CHO + O""",
)

entry(
    index = 262,
    label = "CH3 + C2H3 <=> CH4 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.92e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C2H3 <=> CH4 + C2H2""",
)

entry(
    index = 263,
    label = "C2H3 + H <=> C2H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + H <=> C2H2 + H2""",
)

entry(
    index = 264,
    label = "C2H3 + OH <=> C2H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + OH <=> C2H2 + H2O""",
)

entry(
    index = 265,
    label = "C2H + H <=> C2H2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+17, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.75e+33, 'cm^6/(mol^2*s)'),
            n = -4.8,
            Ea = (1900, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.646,
        T3 = (132, 'K'),
        T1 = (1315, 'K'),
        T2 = (5566, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C2H + H <=> C2H2""",
)

entry(
    index = 266,
    label = "C2H2 + O2 <=> HCCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+08, 'cm^3/(mol*s)'), n=1.5, Ea=(30100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O2 <=> HCCO + OH""",
)

entry(
    index = 267,
    label = "O + C2H2 <=> C2H + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.6e+19, 'cm^3/(mol*s)'),
        n = -1.4,
        Ea = (28950, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + C2H2 <=> C2H + OH""",
)

entry(
    index = 268,
    label = "C2H2 + O <=> CH2 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.94e+06, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> CH2 + CO""",
)

entry(
    index = 269,
    label = "C2H2 + O <=> HCCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.35e+07, 'cm^3/(mol*s)'), n=2, Ea=(1900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + O <=> HCCO + H""",
)

entry(
    index = 270,
    label = "C2H2 + OH <=> C2H + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.37e+07, 'cm^3/(mol*s)'), n=2, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> C2H + H2O""",
)

entry(
    index = 271,
    label = "C2H2 + OH <=> CH2CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.236e+13, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (12000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CH2CO + H""",
)

entry(
    index = 272,
    label = "C2H2 + OH <=> CH3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.000483, 'cm^3/(mol*s)'), n=4, Ea=(-2000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H2 + OH <=> CH3 + CO""",
)

entry(
    index = 273,
    label = "OH + C2H2 <=> H + HCCOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(504000, 'cm^3/(mol*s)'), n=2.3, Ea=(13500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OH + C2H2 <=> H + HCCOH""",
)

entry(
    index = 274,
    label = "H + HCCOH <=> H + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + HCCOH <=> H + CH2CO""",
)

entry(
    index = 275,
    label = "C2H5OH <=> CH2OH + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2e+23, 's^-1'), n=-1.68, Ea=(96400, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.11e+85, 'cm^3/(mol*s)'),
            n = -18.84,
            Ea = (113100, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (550, 'K'),
        T1 = (825, 'K'),
        T2 = (6100, 'K'),
        efficiencies = {'[H][H]': 2, '[C]=O': 2, 'O=C=O': 3, 'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5OH <=> CH2OH + CH3""",
)

entry(
    index = 276,
    label = "C2H5OH <=> C2H5 + OH",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.4e+23, 's^-1'), n=-1.62, Ea=(99540, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.11e+85, 'cm^3/(mol*s)'),
            n = -18.8,
            Ea = (118770, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.5,
        T3 = (650, 'K'),
        T1 = (800, 'K'),
        T2 = (1e+15, 'K'),
        efficiencies = {'[H][H]': 2, '[C]=O': 2, 'O=C=O': 3, 'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5OH <=> C2H5 + OH""",
)

entry(
    index = 277,
    label = "C2H5OH <=> C2H4 + H2O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(132000, 's^-1'), n=2.52, Ea=(60660, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.09e+55, 'cm^3/(mol*s)'),
            n = -10.92,
            Ea = (62644, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.897,
        T3 = (1e+10, 'K'),
        T1 = (1, 'K'),
        T2 = (5e+09, 'K'),
        efficiencies = {'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5OH <=> C2H4 + H2O""",
)

entry(
    index = 278,
    label = "C2H5OH <=> CH3CHO + H2",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.24e+11, 's^-1'), n=0.095, Ea=(91010, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.46e+87, 'cm^3/(mol*s)'),
            n = -19.42,
            Ea = (115580, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.9,
        T3 = (900, 'K'),
        T1 = (1100, 'K'),
        T2 = (3500, 'K'),
        efficiencies = {'O': 5},
    ),
    shortDesc = u"""The chemkin file reaction is C2H5OH <=> CH3CHO + H2""",
)

entry(
    index = 279,
    label = "C2H5OH + O2 <=> PC2H4OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(52800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + O2 <=> PC2H4OH + HO2""",
)

entry(
    index = 280,
    label = "C2H5OH + O2 <=> SC2H4OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 'cm^3/(mol*s)'), n=0, Ea=(50150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + O2 <=> SC2H4OH + HO2""",
)

entry(
    index = 281,
    label = "C2H5OH + OH <=> PC2H4OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+11, 'cm^3/(mol*s)'), n=0.4, Ea=(717, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + OH <=> PC2H4OH + H2O""",
)

entry(
    index = 282,
    label = "C2H5OH + OH <=> SC2H4OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.56e+10, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (-380, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5OH + OH <=> SC2H4OH + H2O""",
)

entry(
    index = 283,
    label = "C2H5OH + OH <=> C2H5O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+10, 'cm^3/(mol*s)'), n=0.8, Ea=(2534, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + OH <=> C2H5O + H2O""",
)

entry(
    index = 284,
    label = "C2H5OH + H <=> PC2H4OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1880, 'cm^3/(mol*s)'), n=3.2, Ea=(7150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + H <=> PC2H4OH + H2""",
)

entry(
    index = 285,
    label = "C2H5OH + H <=> SC2H4OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(179000, 'cm^3/(mol*s)'), n=2.53, Ea=(3420, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + H <=> SC2H4OH + H2""",
)

entry(
    index = 286,
    label = "C2H5OH + H <=> C2H5O + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(53600, 'cm^3/(mol*s)'), n=2.53, Ea=(4405, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + H <=> C2H5O + H2""",
)

entry(
    index = 287,
    label = "C2H5OH + HO2 <=> PC2H4OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23790, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + HO2 <=> PC2H4OH + H2O2""",
)

entry(
    index = 288,
    label = "C2H5OH + HO2 <=> SC2H4OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + HO2 <=> SC2H4OH + H2O2""",
)

entry(
    index = 289,
    label = "C2H5OH + HO2 <=> C2H5O + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + HO2 <=> C2H5O + H2O2""",
)

entry(
    index = 290,
    label = "C2H5OH + CH3O2 <=> PC2H4OH + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(12300, 'cm^3/(mol*s)'), n=2.55, Ea=(15750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + CH3O2 <=> PC2H4OH + CH3O2H""",
)

entry(
    index = 291,
    label = "C2H5OH + CH3O2 <=> SC2H4OH + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8200, 'cm^3/(mol*s)'), n=2.55, Ea=(10750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + CH3O2 <=> SC2H4OH + CH3O2H""",
)

entry(
    index = 292,
    label = "C2H5OH + CH3O2 <=> C2H5O + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + CH3O2 <=> C2H5O + CH3O2H""",
)

entry(
    index = 293,
    label = "C2H5OH + O <=> PC2H4OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(969, 'cm^3/(mol*s)'), n=3.23, Ea=(4658, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + O <=> PC2H4OH + OH""",
)

entry(
    index = 294,
    label = "C2H5OH + O <=> SC2H4OH + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(145000, 'cm^3/(mol*s)'), n=2.47, Ea=(876, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + O <=> SC2H4OH + OH""",
)

entry(
    index = 295,
    label = "C2H5OH + O <=> C2H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00146, 'cm^3/(mol*s)'),
        n = 4.73,
        Ea = (1727, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5OH + O <=> C2H5O + OH""",
)

entry(
    index = 296,
    label = "C2H5OH + CH3 <=> PC2H4OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(330, 'cm^3/(mol*s)'), n=3.3, Ea=(12290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + CH3 <=> PC2H4OH + CH4""",
)

entry(
    index = 297,
    label = "C2H5OH + CH3 <=> SC2H4OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19.93, 'cm^3/(mol*s)'), n=3.37, Ea=(7634, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + CH3 <=> SC2H4OH + CH4""",
)

entry(
    index = 298,
    label = "C2H5OH + CH3 <=> C2H5O + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.035, 'cm^3/(mol*s)'), n=3.57, Ea=(7721, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + CH3 <=> C2H5O + CH4""",
)

entry(
    index = 299,
    label = "C2H5OH + C2H5 <=> PC2H4OH + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(13400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + C2H5 <=> PC2H4OH + C2H6""",
)

entry(
    index = 300,
    label = "C2H5OH + C2H5 <=> SC2H4OH + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5OH + C2H5 <=> SC2H4OH + C2H6""",
)

entry(
    index = 301,
    label = "PC2H4OH <=> C2H4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.047e+25, 's^-1'), n=-3.99, Ea=(30390, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC2H4OH <=> C2H4 + OH""",
)

entry(
    index = 302,
    label = "SC2H4OH <=> CH3CHO + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is SC2H4OH <=> CH3CHO + H""",
)

entry(
    index = 303,
    label = "O2C2H4OH <=> PC2H4OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.9e+16, 's^-1'), n=-1, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C2H4OH <=> PC2H4OH + O2""",
)

entry(
    index = 304,
    label = "O2C2H4OH => OH + CH2O + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(18900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C2H4OH => OH + CH2O + CH2O""",
)

entry(
    index = 305,
    label = "SC2H4OH + O2 <=> CH3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.81e+06, 'cm^3/(mol*s)'), n=2, Ea=(1641, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC2H4OH + O2 <=> CH3CHO + HO2""",
)

entry(
    index = 306,
    label = "CH3COCH3 <=> CH3CO + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.108e+21, 's^-1'), n=-1.57, Ea=(84680, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (7.013e+89, 'cm^3/(mol*s)'),
            n = -20.38,
            Ea = (107150, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.863,
        T3 = (1e+10, 'K'),
        T1 = (416.4, 'K'),
        T2 = (3.29e+09, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 <=> CH3CO + CH3""",
)

entry(
    index = 307,
    label = "CH3COCH3 + OH <=> CH3COCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(125000, 'cm^3/(mol*s)'), n=2.483, Ea=(445, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + OH <=> CH3COCH2 + H2O""",
)

entry(
    index = 308,
    label = "CH3COCH3 + H <=> CH3COCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(980000, 'cm^3/(mol*s)'), n=2.43, Ea=(5160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + H <=> CH3COCH2 + H2""",
)

entry(
    index = 309,
    label = "CH3COCH3 + O <=> CH3COCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.13e+11, 'cm^3/(mol*s)'),
        n = 0.211,
        Ea = (4890, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + O <=> CH3COCH2 + OH""",
)

entry(
    index = 310,
    label = "CH3COCH3 + CH3 <=> CH3COCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.96e+11, 'cm^3/(mol*s)'), n=0, Ea=(9784, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + CH3 <=> CH3COCH2 + CH4""",
)

entry(
    index = 311,
    label = "CH3COCH3 + CH3O <=> CH3COCH2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.34e+11, 'cm^3/(mol*s)'), n=0, Ea=(6460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + CH3O <=> CH3COCH2 + CH3OH""",
)

entry(
    index = 312,
    label = "CH3COCH3 + O2 <=> CH3COCH2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+13, 'cm^3/(mol*s)'), n=0, Ea=(48500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + O2 <=> CH3COCH2 + HO2""",
)

entry(
    index = 313,
    label = "CH3COCH3 + HO2 <=> CH3COCH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + HO2 <=> CH3COCH2 + H2O2""",
)

entry(
    index = 314,
    label = "CH3COCH3 + CH3O2 <=> CH3COCH2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + CH3O2 <=> CH3COCH2 + CH3O2H""",
)

entry(
    index = 315,
    label = "CH3COCH2 <=> CH2CO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(31000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH2 <=> CH2CO + CH3""",
)

entry(
    index = 316,
    label = "CH3COCH2O2 <=> CH3COCH2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.021e+15, 's^-1'), n=-0.956, Ea=(24460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH2O2 <=> CH3COCH2 + O2""",
)

entry(
    index = 317,
    label = "CH3COCH3 + CH3COCH2O2 <=> CH3COCH2 + CH3COCH2O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH3 + CH3COCH2O2 <=> CH3COCH2 + CH3COCH2O2H""",
)

entry(
    index = 318,
    label = "CH2O + CH3COCH2O2 <=> HCO + CH3COCH2O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.288e+11, 'cm^3/(mol*s)'), n=0, Ea=(9000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O + CH3COCH2O2 <=> HCO + CH3COCH2O2H""",
)

entry(
    index = 319,
    label = "HO2 + CH3COCH2O2 => CH3COCH2O2H + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + CH3COCH2O2 => CH3COCH2O2H + O2""",
)

entry(
    index = 320,
    label = "CH3COCH2O2H <=> CH3COCH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH2O2H <=> CH3COCH2O + OH""",
)

entry(
    index = 321,
    label = "CH3COCH2O <=> CH3CO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.732e+20, 's^-1'), n=-2.176, Ea=(17260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3COCH2O <=> CH3CO + CH2O""",
)

entry(
    index = 322,
    label = "C2H3CHO <=> C2H3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.003e+24, 's^-1'), n=-2.135, Ea=(103400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO <=> C2H3 + HCO""",
)

entry(
    index = 323,
    label = "C2H3CHO + H <=> C2H3CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.34e+13, 'cm^3/(mol*s)'), n=0, Ea=(3300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + H <=> C2H3CO + H2""",
)

entry(
    index = 324,
    label = "C2H3CHO + O <=> C2H3CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.94e+12, 'cm^3/(mol*s)'), n=0, Ea=(1868, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + O <=> C2H3CO + OH""",
)

entry(
    index = 325,
    label = "C2H3CHO + OH <=> C2H3CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.24e+06, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (-962, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + OH <=> C2H3CO + H2O""",
)

entry(
    index = 326,
    label = "C2H3CHO + O2 <=> C2H3CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.005e+13, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (40700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + O2 <=> C2H3CO + HO2""",
)

entry(
    index = 327,
    label = "C2H3CHO + HO2 <=> C2H3CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + HO2 <=> C2H3CO + H2O2""",
)

entry(
    index = 328,
    label = "C2H3CHO + CH3 <=> C2H3CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.608e+06, 'cm^3/(mol*s)'),
        n = 1.78,
        Ea = (5911, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + CH3 <=> C2H3CO + CH4""",
)

entry(
    index = 329,
    label = "C2H3CHO + C2H3 <=> C2H3CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.74e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + C2H3 <=> C2H3CO + C2H4""",
)

entry(
    index = 330,
    label = "C2H3CHO + CH3O <=> C2H3CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(3300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + CH3O <=> C2H3CO + CH3OH""",
)

entry(
    index = 331,
    label = "C2H3CHO + CH3O2 <=> C2H3CO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + CH3O2 <=> C2H3CO + CH3O2H""",
)

entry(
    index = 332,
    label = "C2H3CO <=> C2H3 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.37e+21, 's^-1'), n=-2.179, Ea=(39410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CO <=> C2H3 + CO""",
)

entry(
    index = 333,
    label = "C2H5CHO <=> C2H5 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.496e+27, 's^-1'), n=-3.205, Ea=(87040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO <=> C2H5 + HCO""",
)

entry(
    index = 334,
    label = "C2H5CHO + H <=> C2H5CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(4200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + H <=> C2H5CO + H2""",
)

entry(
    index = 335,
    label = "C2H5CHO + O <=> C2H5CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(1790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + O <=> C2H5CO + OH""",
)

entry(
    index = 336,
    label = "C2H5CHO + OH <=> C2H5CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + OH <=> C2H5CO + H2O""",
)

entry(
    index = 337,
    label = "C2H5CHO + CH3 <=> C2H5CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.608e+06, 'cm^3/(mol*s)'),
        n = 1.78,
        Ea = (5911, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + CH3 <=> C2H5CO + CH4""",
)

entry(
    index = 338,
    label = "C2H5CHO + HO2 <=> C2H5CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + HO2 <=> C2H5CO + H2O2""",
)

entry(
    index = 339,
    label = "C2H5CHO + CH3O <=> C2H5CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(3300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + CH3O <=> C2H5CO + CH3OH""",
)

entry(
    index = 340,
    label = "C2H5CHO + CH3O2 <=> C2H5CO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + CH3O2 <=> C2H5CO + CH3O2H""",
)

entry(
    index = 341,
    label = "C2H5CHO + C2H5 <=> C2H5CO + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + C2H5 <=> C2H5CO + C2H6""",
)

entry(
    index = 342,
    label = "C2H5CHO + C2H5O <=> C2H5CO + C2H5OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.026e+11, 'cm^3/(mol*s)'), n=0, Ea=(3300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + C2H5O <=> C2H5CO + C2H5OH""",
)

entry(
    index = 343,
    label = "C2H5CHO + C2H5O2 <=> C2H5CO + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + C2H5O2 <=> C2H5CO + C2H5O2H""",
)

entry(
    index = 344,
    label = "C2H5CHO + O2 <=> C2H5CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.005e+13, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (40700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + O2 <=> C2H5CO + HO2""",
)

entry(
    index = 345,
    label = "C2H5CHO + CH3CO3 <=> C2H5CO + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + CH3CO3 <=> C2H5CO + CH3CO3H""",
)

entry(
    index = 346,
    label = "C2H5CHO + C2H3 <=> C2H5CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + C2H3 <=> C2H5CO + C2H4""",
)

entry(
    index = 347,
    label = "C2H5CO <=> C2H5 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.46e+23, 's^-1'), n=-3.208, Ea=(17550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CO <=> C2H5 + CO""",
)

entry(
    index = 348,
    label = "CH3OCH3 <=> CH3 + CH3O",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(7.25e+21, 's^-1'), n=-0.94, Ea=(80250, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.5e+60, 'cm^3/(mol*s)'),
            n = -11.56,
            Ea = (101000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.183,
        T3 = (1.3, 'K'),
        T1 = (13000, 'K'),
        T2 = (6.71e+09, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 <=> CH3 + CH3O""",
)

entry(
    index = 349,
    label = "CH3OCH3 + OH <=> CH3OCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.324e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-651.7, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + OH <=> CH3OCH2 + H2O""",
)

entry(
    index = 350,
    label = "CH3OCH3 + H <=> CH3OCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.721e+06, 'cm^3/(mol*s)'),
        n = 2.09,
        Ea = (3384, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + H <=> CH3OCH2 + H2""",
)

entry(
    index = 351,
    label = "CH3OCH3 + O <=> CH3OCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.75e+08, 'cm^3/(mol*s)'),
        n = 1.36,
        Ea = (2250, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + O <=> CH3OCH2 + OH""",
)

entry(
    index = 352,
    label = "CH3OCH3 + HO2 <=> CH3OCH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.68e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + HO2 <=> CH3OCH2 + H2O2""",
)

entry(
    index = 353,
    label = "CH3OCH3 + CH3O2 <=> CH3OCH2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.68e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + CH3O2 <=> CH3OCH2 + CH3O2H""",
)

entry(
    index = 354,
    label = "CH3OCH3 + CH3 <=> CH3OCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.445e-06, 'cm^3/(mol*s)'),
        n = 5.73,
        Ea = (5700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + CH3 <=> CH3OCH2 + CH4""",
)

entry(
    index = 355,
    label = "CH3OCH3 + O2 <=> CH3OCH2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(44910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + O2 <=> CH3OCH2 + HO2""",
)

entry(
    index = 356,
    label = "CH3OCH3 + CH3O <=> CH3OCH2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.02e+11, 'cm^3/(mol*s)'), n=0, Ea=(4074, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + CH3O <=> CH3OCH2 + CH3OH""",
)

entry(
    index = 357,
    label = "CH3OCH3 + CH3OCH2O2 <=> CH3OCH2 + CH3OCH2O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + CH3OCH2O2 <=> CH3OCH2 + CH3OCH2O2H""",
)

entry(
    index = 358,
    label = "CH3OCH3 + O2CHO <=> CH3OCH2 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(44250, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + O2CHO <=> CH3OCH2 + HO2CHO""",
)

entry(
    index = 359,
    label = "CH3OCH3 + OCHO <=> CH3OCH2 + HOCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH3 + OCHO <=> CH3OCH2 + HOCHO""",
)

entry(
    index = 360,
    label = "CH3OCH2 <=> CH2O + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.6e+13, 's^-1'), n=0, Ea=(25500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2 <=> CH2O + CH3""",
)

entry(
    index = 361,
    label = "CH3OCH2 + CH3O <=> CH3OCH3 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2 + CH3O <=> CH3OCH3 + CH2O""",
)

entry(
    index = 362,
    label = "CH3OCH2 + CH2O <=> CH3OCH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5490, 'cm^3/(mol*s)'), n=2.8, Ea=(5862, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2 + CH2O <=> CH3OCH3 + HCO""",
)

entry(
    index = 363,
    label = "CH3OCH2 + CH3CHO <=> CH3OCH3 + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.26e+12, 'cm^3/(mol*s)'), n=0, Ea=(8499, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2 + CH3CHO <=> CH3OCH3 + CH3CO""",
)

entry(
    index = 364,
    label = "CH3OCH2O2 <=> CH3OCH2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.439e+19, 's^-1'), n=-1.594, Ea=(36240, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O2 <=> CH3OCH2 + O2""",
)

entry(
    index = 365,
    label = "CH3OCH2O2 + CH2O <=> CH3OCH2O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O2 + CH2O <=> CH3OCH2O2H + HCO""",
)

entry(
    index = 366,
    label = "CH3OCH2O2 + CH3CHO <=> CH3OCH2O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O2 + CH3CHO <=> CH3OCH2O2H + CH3CO""",
)

entry(
    index = 367,
    label = "CH3OCH2O2 + CH3OCH2O2 => O2 + CH3OCH2O + CH3OCH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.21e+23, 'cm^3/(mol*s)'), n=-4.5, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O2 + CH3OCH2O2 => O2 + CH3OCH2O + CH3OCH2O""",
)

entry(
    index = 368,
    label = "CH3OCH2O2H <=> CH3OCH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.106e+22, 's^-1'), n=-2.124, Ea=(43830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O2H <=> CH3OCH2O + OH""",
)

entry(
    index = 369,
    label = "CH3OCH2O <=> CH3O + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.384e+19, 's^-1'), n=-2.014, Ea=(25190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O <=> CH3O + CH2O""",
)

entry(
    index = 370,
    label = "CH3OCH2O + O2 <=> CH3OCHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O + O2 <=> CH3OCHO + HO2""",
)

entry(
    index = 371,
    label = "CH3OCH2O <=> CH3OCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.06e+12, 's^-1'), n=0.056, Ea=(8218, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O <=> CH3OCHO + H""",
)

entry(
    index = 372,
    label = "CH3OCH2O2 <=> CH2OCH2O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+10, 's^-1'), n=0, Ea=(21580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCH2O2 <=> CH2OCH2O2H""",
)

entry(
    index = 373,
    label = "CH2OCH2O2H => OH + CH2O + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(20760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCH2O2H => OH + CH2O + CH2O""",
)

entry(
    index = 374,
    label = "O2CH2OCH2O2H <=> CH2OCH2O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.924e+19, 's^-1'), n=-1.622, Ea=(36270, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2CH2OCH2O2H <=> CH2OCH2O2H + O2""",
)

entry(
    index = 375,
    label = "O2CH2OCH2O2H <=> HO2CH2OCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+10, 's^-1'), n=0, Ea=(18580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2CH2OCH2O2H <=> HO2CH2OCHO + OH""",
)

entry(
    index = 376,
    label = "HO2CH2OCHO <=> OCH2OCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+16, 's^-1'), n=0, Ea=(40500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2CH2OCHO <=> OCH2OCHO + OH""",
)

entry(
    index = 377,
    label = "OCH2OCHO <=> CH2O + OCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.898e+19, 's^-1'), n=-2.201, Ea=(31850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCH2OCHO <=> CH2O + OCHO""",
)

entry(
    index = 378,
    label = "OCH2OCHO <=> HOCH2OCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is OCH2OCHO <=> HOCH2OCO""",
)

entry(
    index = 379,
    label = "HOCH2OCO <=> HOCH2O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.238e+19, 's^-1'), n=-2.021, Ea=(19690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2OCO <=> HOCH2O + CO""",
)

entry(
    index = 380,
    label = "HOCH2OCO <=> CH2OH + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.413e+17, 's^-1'), n=-1.574, Ea=(22120, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOCH2OCO <=> CH2OH + CO2""",
)

entry(
    index = 381,
    label = "CH3OCHO <=> CH2OCHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.241e+19, 's^-1'), n=-1.15, Ea=(102500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO <=> CH2OCHO + H""",
)

entry(
    index = 382,
    label = "CH3OCHO <=> CH3OCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.325e+19, 's^-1'), n=-1, Ea=(100100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO <=> CH3OCO + H""",
)

entry(
    index = 383,
    label = "CH3OCHO <=> CH3OH + CO",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(62500, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.143e+60, 'cm^3/(mol*s)'),
            n = -12.07,
            Ea = (75400, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.78,
        T3 = (8.28e+09, 'K'),
        T1 = (438.9, 'K'),
        T2 = (6.7e+08, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCHO <=> CH3OH + CO""",
)

entry(
    index = 384,
    label = "CH3OCHO <=> CH3O + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.374e+16, 's^-1'), n=-0.01, Ea=(97090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO <=> CH3O + HCO""",
)

entry(
    index = 385,
    label = "CH3OCHO <=> CH3 + OCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.214e+17, 's^-1'), n=-0.53, Ea=(79970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO <=> CH3 + OCHO""",
)

entry(
    index = 386,
    label = "CH3OCHO + O2 <=> CH3OCO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(49700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + O2 <=> CH3OCO + HO2""",
)

entry(
    index = 387,
    label = "CH3OCHO + O2 <=> CH2OCHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.05e+13, 'cm^3/(mol*s)'), n=0, Ea=(52000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + O2 <=> CH2OCHO + HO2""",
)

entry(
    index = 388,
    label = "CH3OCHO + OH <=> CH3OCO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.58e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(934, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + OH <=> CH3OCO + H2O""",
)

entry(
    index = 389,
    label = "CH3OCHO + OH <=> CH2OCHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + OH <=> CH2OCHO + H2O""",
)

entry(
    index = 390,
    label = "CH3OCHO + HO2 <=> CH3OCO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + HO2 <=> CH3OCO + H2O2""",
)

entry(
    index = 391,
    label = "CH3OCHO + HO2 <=> CH2OCHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + HO2 <=> CH2OCHO + H2O2""",
)

entry(
    index = 392,
    label = "CH3OCHO + O <=> CH3OCO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(275500, 'cm^3/(mol*s)'), n=2.45, Ea=(2830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + O <=> CH3OCO + OH""",
)

entry(
    index = 393,
    label = "CH3OCHO + O <=> CH2OCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(980000, 'cm^3/(mol*s)'), n=2.43, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + O <=> CH2OCHO + OH""",
)

entry(
    index = 394,
    label = "CH3OCHO + H <=> CH3OCO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(650000, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + H <=> CH3OCO + H2""",
)

entry(
    index = 395,
    label = "CH3OCHO + H <=> CH2OCHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665000, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + H <=> CH2OCHO + H2""",
)

entry(
    index = 396,
    label = "CH3OCHO + CH3 <=> CH3OCO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.755, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + CH3 <=> CH3OCO + CH4""",
)

entry(
    index = 397,
    label = "CH3OCHO + CH3 <=> CH2OCHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.452, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + CH3 <=> CH2OCHO + CH4""",
)

entry(
    index = 398,
    label = "CH3OCHO + CH3O <=> CH3OCO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.48e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + CH3O <=> CH3OCO + CH3OH""",
)

entry(
    index = 399,
    label = "CH3OCHO + CH3O <=> CH2OCHO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + CH3O <=> CH2OCHO + CH3OH""",
)

entry(
    index = 400,
    label = "CH3OCHO + CH3O2 <=> CH3OCO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + CH3O2 <=> CH3OCO + CH3O2H""",
)

entry(
    index = 401,
    label = "CH3OCHO + CH3O2 <=> CH2OCHO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + CH3O2 <=> CH2OCHO + CH3O2H""",
)

entry(
    index = 402,
    label = "CH3OCHO + HCO <=> CH3OCO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.4e+06, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (17010, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + HCO <=> CH3OCO + CH2O""",
)

entry(
    index = 403,
    label = "CH3OCHO + HCO <=> CH2OCHO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(102500, 'cm^3/(mol*s)'), n=2.5, Ea=(18430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCHO + HCO <=> CH2OCHO + CH2O""",
)

entry(
    index = 404,
    label = "CH2OCHO <=> CH3OCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.62e+11, 's^-1'), n=-0.03, Ea=(38180, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCHO <=> CH3OCO""",
)

entry(
    index = 405,
    label = "CH3OCO <=> CH3 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.59e+14, 's^-1'), n=-0.172, Ea=(16010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCO <=> CH3 + CO2""",
)

entry(
    index = 406,
    label = "CH3OCO <=> CH3O + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.431e+15, 's^-1'), n=-0.041, Ea=(23770, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OCO <=> CH3O + CO""",
)

entry(
    index = 407,
    label = "CH2OCHO <=> CH2O + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.661e+12, 's^-1'), n=0.12, Ea=(27440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2OCHO <=> CH2O + HCO""",
)

entry(
    index = 408,
    label = "C3H8 <=> CH3 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.29e+37, 's^-1'), n=-5.84, Ea=(97380, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (5.64e+74, 'cm^3/(mol*s)'),
            n = -15.74,
            Ea = (98714, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.31,
        T3 = (50, 'K'),
        T1 = (3000, 'K'),
        T2 = (9000, 'K'),
        efficiencies = {'C': 2, 'O=C=O': 2, 'CC': 3, 'O': 6, '[H][H]': 2, '[He]': 0.7, '[C]=O': 1.5, '[Ar]': 0.7},
    ),
    shortDesc = u"""The chemkin file reaction is C3H8 <=> CH3 + C2H5""",
)

entry(
    index = 409,
    label = "C3H8 <=> NC3H7 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+17, 's^-1'), n=-0.357, Ea=(101200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 <=> NC3H7 + H""",
)

entry(
    index = 410,
    label = "C3H8 <=> IC3H7 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.377e+18, 's^-1'), n=-0.671, Ea=(98680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 <=> IC3H7 + H""",
)

entry(
    index = 411,
    label = "C3H8 + O2 <=> IC3H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(49640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O2 <=> IC3H7 + HO2""",
)

entry(
    index = 412,
    label = "C3H8 + O2 <=> NC3H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O2 <=> NC3H7 + HO2""",
)

entry(
    index = 413,
    label = "H + C3H8 <=> H2 + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + C3H8 <=> H2 + IC3H7""",
)

entry(
    index = 414,
    label = "H + C3H8 <=> H2 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.33e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6756, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + C3H8 <=> H2 + NC3H7""",
)

entry(
    index = 415,
    label = "C3H8 + O <=> IC3H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(549000, 'cm^3/(mol*s)'), n=2.5, Ea=(3140, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O <=> IC3H7 + OH""",
)

entry(
    index = 416,
    label = "C3H8 + O <=> NC3H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.71e+06, 'cm^3/(mol*s)'),
        n = 2.4,
        Ea = (5505, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H8 + O <=> NC3H7 + OH""",
)

entry(
    index = 417,
    label = "C3H8 + OH <=> NC3H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H8 + OH <=> NC3H7 + H2O""",
)

entry(
    index = 418,
    label = "C3H8 + OH <=> IC3H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H8 + OH <=> IC3H7 + H2O""",
)

entry(
    index = 419,
    label = "C3H8 + HO2 <=> IC3H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(58800, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + HO2 <=> IC3H7 + H2O2""",
)

entry(
    index = 420,
    label = "C3H8 + HO2 <=> NC3H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(81000, 'cm^3/(mol*s)'), n=2.5, Ea=(16690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + HO2 <=> NC3H7 + H2O2""",
)

entry(
    index = 421,
    label = "CH3 + C3H8 <=> CH4 + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(64000, 'cm^3/(mol*s)'), n=2.17, Ea=(7520, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C3H8 <=> CH4 + IC3H7""",
)

entry(
    index = 422,
    label = "CH3 + C3H8 <=> CH4 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.904, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C3H8 <=> CH4 + NC3H7""",
)

entry(
    index = 423,
    label = "IC3H7 + C3H8 <=> NC3H7 + C3H8",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+10, 'cm^3/(mol*s)'), n=0, Ea=(12900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + C3H8 <=> NC3H7 + C3H8""",
)

entry(
    index = 424,
    label = "C2H3 + C3H8 <=> C2H4 + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C3H8 <=> C2H4 + IC3H7""",
)

entry(
    index = 425,
    label = "C2H3 + C3H8 <=> C2H4 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C3H8 <=> C2H4 + NC3H7""",
)

entry(
    index = 426,
    label = "C2H5 + C3H8 <=> C2H6 + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C3H8 <=> C2H6 + IC3H7""",
)

entry(
    index = 427,
    label = "C2H5 + C3H8 <=> C2H6 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C3H8 <=> C2H6 + NC3H7""",
)

entry(
    index = 428,
    label = "C3H8 + C3H5-A <=> NC3H7 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(20500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + C3H5-A <=> NC3H7 + C3H6""",
)

entry(
    index = 429,
    label = "C3H8 + C3H5-A <=> IC3H7 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(16200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + C3H5-A <=> IC3H7 + C3H6""",
)

entry(
    index = 430,
    label = "C3H8 + CH3O <=> NC3H7 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + CH3O <=> NC3H7 + CH3OH""",
)

entry(
    index = 431,
    label = "C3H8 + CH3O <=> IC3H7 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + CH3O <=> IC3H7 + CH3OH""",
)

entry(
    index = 432,
    label = "CH3O2 + C3H8 <=> CH3O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(81000, 'cm^3/(mol*s)'), n=2.5, Ea=(16690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + C3H8 <=> CH3O2H + NC3H7""",
)

entry(
    index = 433,
    label = "CH3O2 + C3H8 <=> CH3O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(58800, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + C3H8 <=> CH3O2H + IC3H7""",
)

entry(
    index = 434,
    label = "C2H5O2 + C3H8 <=> C2H5O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(81000, 'cm^3/(mol*s)'), n=2.5, Ea=(16690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 + C3H8 <=> C2H5O2H + NC3H7""",
)

entry(
    index = 435,
    label = "C2H5O2 + C3H8 <=> C2H5O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(58800, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 + C3H8 <=> C2H5O2H + IC3H7""",
)

entry(
    index = 436,
    label = "NC3H7O2 + C3H8 <=> NC3H7O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C3H8 <=> NC3H7O2H + NC3H7""",
)

entry(
    index = 437,
    label = "NC3H7O2 + C3H8 <=> NC3H7O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C3H8 <=> NC3H7O2H + IC3H7""",
)

entry(
    index = 438,
    label = "IC3H7O2 + C3H8 <=> IC3H7O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C3H8 <=> IC3H7O2H + NC3H7""",
)

entry(
    index = 439,
    label = "IC3H7O2 + C3H8 <=> IC3H7O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C3H8 <=> IC3H7O2H + IC3H7""",
)

entry(
    index = 440,
    label = "C3H8 + CH3CO3 <=> IC3H7 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + CH3CO3 <=> IC3H7 + CH3CO3H""",
)

entry(
    index = 441,
    label = "C3H8 + CH3CO3 <=> NC3H7 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + CH3CO3 <=> NC3H7 + CH3CO3H""",
)

entry(
    index = 442,
    label = "C3H8 + O2CHO <=> NC3H7 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(55200, 'cm^3/(mol*s)'), n=2.55, Ea=(16480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O2CHO <=> NC3H7 + HO2CHO""",
)

entry(
    index = 443,
    label = "C3H8 + O2CHO <=> IC3H7 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14750, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H8 + O2CHO <=> IC3H7 + HO2CHO""",
)

entry(
    index = 444,
    label = "IC3H7 <=> H + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.919e+13, 's^-1'), n=-0.025, Ea=(37690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 <=> H + C3H6""",
)

entry(
    index = 445,
    label = "IC3H7 + H <=> C2H5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + H <=> C2H5 + CH3""",
)

entry(
    index = 446,
    label = "IC3H7 + O2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e-19, 'cm^3/(mol*s)'), n=0, Ea=(5020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + O2 <=> C3H6 + HO2""",
)

entry(
    index = 447,
    label = "IC3H7 + OH <=> C3H6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + OH <=> C3H6 + H2O""",
)

entry(
    index = 448,
    label = "IC3H7 + O <=> CH3COCH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.818e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + O <=> CH3COCH3 + H""",
)

entry(
    index = 449,
    label = "IC3H7 + O <=> CH3CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.818e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + O <=> CH3CHO + CH3""",
)

entry(
    index = 450,
    label = "NC3H7 <=> CH3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.97e+40, 's^-1'), n=-8.6, Ea=(41430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7 <=> CH3 + C2H4""",
)

entry(
    index = 451,
    label = "NC3H7 <=> H + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.78e+39, 's^-1'), n=-8.1, Ea=(46580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7 <=> H + C3H6""",
)

entry(
    index = 452,
    label = "NC3H7 + O2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-19, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7 + O2 <=> C3H6 + HO2""",
)

entry(
    index = 453,
    label = "C2H5CHO + NC3H7 <=> C2H5CO + C3H8",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + NC3H7 <=> C2H5CO + C3H8""",
)

entry(
    index = 454,
    label = "C2H5CHO + IC3H7 <=> C2H5CO + C3H8",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + IC3H7 <=> C2H5CO + C3H8""",
)

entry(
    index = 455,
    label = "C2H5CHO + C3H5-A <=> C2H5CO + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHO + C3H5-A <=> C2H5CO + C3H6""",
)

entry(
    index = 456,
    label = "C3H6 <=> C2H3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.73e+62, 's^-1'), n=-13.28, Ea=(123200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 <=> C2H3 + CH3""",
)

entry(
    index = 457,
    label = "C3H6 <=> C3H5-A + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.01e+61, 's^-1'), n=-13.26, Ea=(118500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 <=> C3H5-A + H""",
)

entry(
    index = 458,
    label = "C3H6 <=> C3H5-S + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.71e+69, 's^-1'), n=-16.09, Ea=(140000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 <=> C3H5-S + H""",
)

entry(
    index = 459,
    label = "C3H6 <=> C3H5-T + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.62e+71, 's^-1'), n=-16.58, Ea=(139300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 <=> C3H5-T + H""",
)

entry(
    index = 460,
    label = "C3H6 + O <=> C2H5 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.58e+07, 'cm^3/(mol*s)'),
        n = 1.76,
        Ea = (-1216, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> C2H5 + HCO""",
)

entry(
    index = 461,
    label = "C3H6 + O => CH2CO + CH3 + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+07, 'cm^3/(mol*s)'), n=1.76, Ea=(76, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O => CH2CO + CH3 + H""",
)

entry(
    index = 462,
    label = "C3H6 + O => CH3CHCO + H + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+07, 'cm^3/(mol*s)'), n=1.76, Ea=(76, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O => CH3CHCO + H + H""",
)

entry(
    index = 463,
    label = "C3H6 + O <=> C3H5-A + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.24e+11, 'cm^3/(mol*s)'),
        n = 0.7,
        Ea = (5884, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> C3H5-A + OH""",
)

entry(
    index = 464,
    label = "C3H6 + O <=> C3H5-S + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0.7, Ea=(8959, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> C3H5-S + OH""",
)

entry(
    index = 465,
    label = "C3H6 + O <=> C3H5-T + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.03e+10, 'cm^3/(mol*s)'),
        n = 0.7,
        Ea = (7632, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H6 + O <=> C3H5-T + OH""",
)

entry(
    index = 466,
    label = "C3H6 + OH <=> C3H5-A + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + OH <=> C3H5-A + H2O""",
)

entry(
    index = 467,
    label = "C3H6 + OH <=> C3H5-S + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.11e+06, 'cm^3/(mol*s)'), n=2, Ea=(2778, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + OH <=> C3H5-S + H2O""",
)

entry(
    index = 468,
    label = "C3H6 + OH <=> C3H5-T + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.11e+06, 'cm^3/(mol*s)'), n=2, Ea=(1451, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + OH <=> C3H5-T + H2O""",
)

entry(
    index = 469,
    label = "C3H6 + HO2 <=> C3H5-A + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27000, 'cm^3/(mol*s)'), n=2.5, Ea=(12340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + HO2 <=> C3H5-A + H2O2""",
)

entry(
    index = 470,
    label = "C3H6 + HO2 <=> C3H5-S + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18000, 'cm^3/(mol*s)'), n=2.5, Ea=(27620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + HO2 <=> C3H5-S + H2O2""",
)

entry(
    index = 471,
    label = "C3H6 + HO2 <=> C3H5-T + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9000, 'cm^3/(mol*s)'), n=2.5, Ea=(23590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + HO2 <=> C3H5-T + H2O2""",
)

entry(
    index = 472,
    label = "C3H6 + H <=> C3H5-A + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(173000, 'cm^3/(mol*s)'), n=2.5, Ea=(2492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> C3H5-A + H2""",
)

entry(
    index = 473,
    label = "C3H6 + H <=> C3H5-S + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(804000, 'cm^3/(mol*s)'), n=2.5, Ea=(12280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> C3H5-S + H2""",
)

entry(
    index = 474,
    label = "C3H6 + H <=> C3H5-T + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(405000, 'cm^3/(mol*s)'), n=2.5, Ea=(9794, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> C3H5-T + H2""",
)

entry(
    index = 475,
    label = "C3H6 + H <=> C2H4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.3e+13, 'cm^3/(mol*s)'), n=0, Ea=(2547, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + H <=> C2H4 + CH3""",
)

entry(
    index = 476,
    label = "C3H6 + O2 <=> C3H5-A + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(39900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O2 <=> C3H5-A + HO2""",
)

entry(
    index = 477,
    label = "C3H6 + O2 <=> C3H5-S + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(62900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O2 <=> C3H5-S + HO2""",
)

entry(
    index = 478,
    label = "C3H6 + O2 <=> C3H5-T + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(60700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + O2 <=> C3H5-T + HO2""",
)

entry(
    index = 479,
    label = "C3H6 + CH3 <=> C3H5-A + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.21, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> C3H5-A + CH4""",
)

entry(
    index = 480,
    label = "C3H6 + CH3 <=> C3H5-S + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.348, 'cm^3/(mol*s)'), n=3.5, Ea=(12850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> C3H5-S + CH4""",
)

entry(
    index = 481,
    label = "C3H6 + CH3 <=> C3H5-T + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.84, 'cm^3/(mol*s)'), n=3.5, Ea=(11660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3 <=> C3H5-T + CH4""",
)

entry(
    index = 482,
    label = "C3H6 + C2H5 <=> C3H5-A + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(9800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + C2H5 <=> C3H5-A + C2H6""",
)

entry(
    index = 483,
    label = "C3H6 + CH3CO3 <=> C3H5-A + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3CO3 <=> C3H5-A + CH3CO3H""",
)

entry(
    index = 484,
    label = "C3H6 + CH3O2 <=> C3H5-A + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + CH3O2 <=> C3H5-A + CH3O2H""",
)

entry(
    index = 485,
    label = "C3H6 + HO2 <=> C3H6O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.29e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + HO2 <=> C3H6O1-2 + OH""",
)

entry(
    index = 486,
    label = "C3H6 + C2H5O2 <=> C3H5-A + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + C2H5O2 <=> C3H5-A + C2H5O2H""",
)

entry(
    index = 487,
    label = "C3H6 + NC3H7O2 <=> C3H5-A + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + NC3H7O2 <=> C3H5-A + NC3H7O2H""",
)

entry(
    index = 488,
    label = "C3H6 + IC3H7O2 <=> C3H5-A + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + IC3H7O2 <=> C3H5-A + IC3H7O2H""",
)

entry(
    index = 489,
    label = "C3H6OH <=> C3H6 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.343e+15, 's^-1'), n=-0.805, Ea=(27900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OH <=> C3H6 + OH""",
)

entry(
    index = 490,
    label = "HOC3H6O2 <=> C3H6OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.873e+19, 's^-1'), n=-1.897, Ea=(34290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOC3H6O2 <=> C3H6OH + O2""",
)

entry(
    index = 491,
    label = "HOC3H6O2 => CH3CHO + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(18900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HOC3H6O2 => CH3CHO + CH2O + OH""",
)

entry(
    index = 492,
    label = "C3H5-A <=> C2H2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.397e+48, 's^-1'), n=-9.9, Ea=(82080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A <=> C2H2 + CH3""",
)

entry(
    index = 493,
    label = "C3H5-A <=> C3H4-A + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.194e+13, 's^-1'), n=0.216, Ea=(61930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A <=> C3H4-A + H""",
)

entry(
    index = 494,
    label = "C3H5-A + HO2 <=> C3H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + HO2 <=> C3H5O + OH""",
)

entry(
    index = 495,
    label = "C3H5-A + CH3O2 <=> C3H5O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + CH3O2 <=> C3H5O + CH3O""",
)

entry(
    index = 496,
    label = "C3H5-A + H <=> C3H4-A + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1232, 'cm^3/(mol*s)'), n=3.035, Ea=(2582, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + H <=> C3H4-A + H2""",
)

entry(
    index = 497,
    label = "C3H5-A + CH3 <=> C3H4-A + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + CH3 <=> C3H4-A + CH4""",
)

entry(
    index = 498,
    label = "C3H5-A + C2H5 <=> C2H6 + C3H4-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + C2H5 <=> C2H6 + C3H4-A""",
)

entry(
    index = 499,
    label = "C3H5-A + C2H5 <=> C2H4 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + C2H5 <=> C2H4 + C3H6""",
)

entry(
    index = 500,
    label = "C3H5-A + C2H3 <=> C2H4 + C3H4-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + C2H3 <=> C2H4 + C3H4-A""",
)

entry(
    index = 501,
    label = "C3H4-A + C3H6 <=> C3H5-A + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.749e+08, 'cm^3/(mol*s)'),
        n = 0.734,
        Ea = (28700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H4-A + C3H6 <=> C3H5-A + C3H5-A""",
)

entry(
    index = 502,
    label = "C3H5-A + O2 <=> C3H4-A + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.18e+21, 'cm^3/(mol*s)'),
        n = -2.85,
        Ea = (30760, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-A + O2 <=> C3H4-A + HO2""",
)

entry(
    index = 503,
    label = "C3H5-A + O2 <=> CH2CHO + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.14e+15, 'cm^3/(mol*s)'),
        n = -1.21,
        Ea = (21050, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-A + O2 <=> CH2CHO + CH2O""",
)

entry(
    index = 504,
    label = "C3H5-A + O2 <=> C2H3CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.47e+13, 'cm^3/(mol*s)'),
        n = -0.44,
        Ea = (23020, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-A + O2 <=> C2H3CHO + OH""",
)

entry(
    index = 505,
    label = "C3H5-A + O2 => C2H2 + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (9.72e+29, 'cm^3/(mol*s)'),
        n = -5.71,
        Ea = (21450, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-A + O2 => C2H2 + CH2O + OH""",
)

entry(
    index = 506,
    label = "C3H5-S <=> C2H2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.598e+39, 's^-1'), n=-8.17, Ea=(42030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-S <=> C2H2 + CH3""",
)

entry(
    index = 507,
    label = "C3H5-S <=> C3H4-P + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.187e+15, 's^-1'), n=-0.79, Ea=(37480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-S <=> C3H4-P + H""",
)

entry(
    index = 508,
    label = "C3H5-S + O2 <=> CH3CHO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.335e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-S + O2 <=> CH3CHO + HCO""",
)

entry(
    index = 509,
    label = "C3H5-S + H <=> C3H4-A + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.333e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-S + H <=> C3H4-A + H2""",
)

entry(
    index = 510,
    label = "C3H5-S + CH3 <=> C3H4-A + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-S + CH3 <=> C3H4-A + CH4""",
)

entry(
    index = 511,
    label = "C3H5-T <=> C2H2 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.163e+40, 's^-1'), n=-8.31, Ea=(45110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-T <=> C2H2 + CH3""",
)

entry(
    index = 512,
    label = "C3H5-T <=> C3H4-A + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.508e+14, 's^-1'), n=-0.44, Ea=(40890, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-T <=> C3H4-A + H""",
)

entry(
    index = 513,
    label = "C3H5-T <=> C3H4-P + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.075e+15, 's^-1'), n=-0.6, Ea=(38490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-T <=> C3H4-P + H""",
)

entry(
    index = 514,
    label = "C3H5-T + O2 <=> C3H4-A + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.89e+30, 'cm^3/(mol*s)'),
        n = -5.59,
        Ea = (15540, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-T + O2 <=> C3H4-A + HO2""",
)

entry(
    index = 515,
    label = "C3H5-T + O2 <=> CH3COCH2 + O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.81e+17, 'cm^3/(mol*s)'),
        n = -1.36,
        Ea = (5580, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-T + O2 <=> CH3COCH2 + O""",
)

entry(
    index = 516,
    label = "C3H5-T + O2 <=> CH2O + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.71e+25, 'cm^3/(mol*s)'),
        n = -3.96,
        Ea = (7043, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5-T + O2 <=> CH2O + CH3CO""",
)

entry(
    index = 517,
    label = "C3H5-T + H <=> C3H4-P + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.333e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-T + H <=> C3H4-P + H2""",
)

entry(
    index = 518,
    label = "C3H5-T + CH3 <=> C3H4-P + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-T + CH3 <=> C3H4-P + CH4""",
)

entry(
    index = 519,
    label = "C3H4-A <=> C3H3 + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.143e+17, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (70000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C3H4-A <=> C3H3 + H""",
)

entry(
    index = 520,
    label = "C3H4-A <=> C3H4-P",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.202e+15, 's^-1'), n=0, Ea=(92400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A <=> C3H4-P""",
)

entry(
    index = 521,
    label = "C3H4-A + O2 <=> C3H3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(39160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + O2 <=> C3H3 + HO2""",
)

entry(
    index = 522,
    label = "C3H4-A + HO2 <=> CH2CO + CH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + HO2 <=> CH2CO + CH2 + OH""",
)

entry(
    index = 523,
    label = "C3H4-A + OH <=> CH2CO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+12, 'cm^3/(mol*s)'), n=0, Ea=(-397, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + OH <=> CH2CO + CH3""",
)

entry(
    index = 524,
    label = "C3H4-A + OH <=> C3H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + OH <=> C3H3 + H2O""",
)

entry(
    index = 525,
    label = "C3H4-A + O <=> C2H4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(1600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + O <=> C2H4 + CO""",
)

entry(
    index = 526,
    label = "C3H4-A + O <=> C2H2 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.003, 'cm^3/(mol*s)'), n=4.61, Ea=(-4243, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + O <=> C2H2 + CH2O""",
)

entry(
    index = 527,
    label = "C3H4-A + H <=> C3H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+07, 'cm^3/(mol*s)'), n=2, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + H <=> C3H3 + H2""",
)

entry(
    index = 528,
    label = "C3H4-A + CH3 <=> C3H3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0367, 'cm^3/(mol*s)'), n=4.01, Ea=(6830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + CH3 <=> C3H3 + CH4""",
)

entry(
    index = 529,
    label = "C3H4-A + C3H5-A <=> C3H3 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + C3H5-A <=> C3H3 + C3H6""",
)

entry(
    index = 530,
    label = "C3H4-A + C2H <=> C3H3 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + C2H <=> C3H3 + C2H2""",
)

entry(
    index = 531,
    label = "C3H4-P <=> C3H3 + H",
    degeneracy = 1,
    kinetics = ThirdBody(
        arrheniusLow = Arrhenius(
            A = (1.143e+17, 'cm^3/(mol*s)'),
            n = 0,
            Ea = (70000, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C3H4-P <=> C3H3 + H""",
)

entry(
    index = 532,
    label = "C3H4-P <=> C2H + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+16, 's^-1'), n=0, Ea=(100000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P <=> C2H + CH3""",
)

entry(
    index = 533,
    label = "C3H4-P + O2 <=> HCCO + OH + CH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=1.5, Ea=(30100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + O2 <=> HCCO + OH + CH2""",
)

entry(
    index = 534,
    label = "C3H4-P + O2 <=> C3H3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + O2 <=> C3H3 + HO2""",
)

entry(
    index = 535,
    label = "C3H4-P + HO2 <=> C2H4 + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + HO2 <=> C2H4 + CO + OH""",
)

entry(
    index = 536,
    label = "C3H4-P + OH <=> C3H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+07, 'cm^3/(mol*s)'), n=2, Ea=(1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + OH <=> C3H3 + H2O""",
)

entry(
    index = 537,
    label = "C3H4-P + OH <=> CH2CO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.0005, 'cm^3/(mol*s)'), n=4.5, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + OH <=> CH2CO + CH3""",
)

entry(
    index = 538,
    label = "C3H4-P + O <=> C2H3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(2010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + O <=> C2H3 + HCO""",
)

entry(
    index = 539,
    label = "C3H4-P + O <=> HCCO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.6e+08, 'cm^3/(mol*s)'), n=1, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + O <=> HCCO + CH3""",
)

entry(
    index = 540,
    label = "C3H4-P + O <=> HCCO + CH2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e-19, 'cm^3/(mol*s)'), n=0, Ea=(2010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + O <=> HCCO + CH2 + H""",
)

entry(
    index = 541,
    label = "C3H4-P + O <=> C3H3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.65e+08, 'cm^3/(mol*s)'),
        n = 1.5,
        Ea = (8600, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H4-P + O <=> C3H3 + OH""",
)

entry(
    index = 542,
    label = "C3H4-P + H <=> C3H3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+07, 'cm^3/(mol*s)'), n=2, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + H <=> C3H3 + H2""",
)

entry(
    index = 543,
    label = "C3H4-P + CH3 <=> C3H3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5, 'cm^3/(mol*s)'), n=3.5, Ea=(5600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + CH3 <=> C3H3 + CH4""",
)

entry(
    index = 544,
    label = "C3H4-P + C2H <=> C3H3 + C2H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + C2H <=> C3H3 + C2H2""",
)

entry(
    index = 545,
    label = "C3H4-P + C2H3 <=> C3H3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + C2H3 <=> C3H3 + C2H4""",
)

entry(
    index = 546,
    label = "C3H4-P + C3H5-A <=> C3H3 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-P + C3H5-A <=> C3H3 + C3H6""",
)

entry(
    index = 547,
    label = "C3H3 + O <=> CH2O + C2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + O <=> CH2O + C2H""",
)

entry(
    index = 548,
    label = "C3H3 + OH <=> C3H2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + OH <=> C3H2 + H2O""",
)

entry(
    index = 549,
    label = "C3H3 + O2 <=> CH2CO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+10, 'cm^3/(mol*s)'), n=0, Ea=(2870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + O2 <=> CH2CO + HCO""",
)

entry(
    index = 550,
    label = "C3H3 + CH3 <=> C2H5 + C2H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.299e+15, 'cm^3/(mol*s)'),
        n = -0.79,
        Ea = (45630, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H3 + CH3 <=> C2H5 + C2H""",
)

entry(
    index = 551,
    label = "C3H2 + O2 <=> HCO + HCCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H2 + O2 <=> HCO + HCCO""",
)

entry(
    index = 552,
    label = "C3H4-A + HO2 <=> C2H4 + CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + HO2 <=> C2H4 + CO + OH""",
)

entry(
    index = 553,
    label = "C3H4-A + HO2 <=> C3H3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(14000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H4-A + HO2 <=> C3H3 + H2O2""",
)

entry(
    index = 554,
    label = "C2H2 + CH3 <=> C3H4-P + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.229e+08, 'cm^3/(mol*s)'),
        n = 1.143,
        Ea = (12090, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> C3H4-P + H""",
)

entry(
    index = 555,
    label = "C2H2 + CH3 <=> C3H4-A + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.74e+19, 'cm^3/(mol*s)'),
        n = -2.08,
        Ea = (31590, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H2 + CH3 <=> C3H4-A + H""",
)

entry(
    index = 556,
    label = "C3H3 + H <=> C3H2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H3 + H <=> C3H2 + H2""",
)

entry(
    index = 557,
    label = "C3H2 + OH <=> C2H2 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H2 + OH <=> C2H2 + HCO""",
)

entry(
    index = 558,
    label = "C3H2 + O2 => HCCO + CO + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H2 + O2 => HCCO + CO + H""",
)

entry(
    index = 559,
    label = "CH3CHCO + OH => C2H5 + CO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.73e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCO + OH => C2H5 + CO2""",
)

entry(
    index = 560,
    label = "CH3CHCO + OH => SC2H4OH + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCO + OH => SC2H4OH + CO""",
)

entry(
    index = 561,
    label = "CH3CHCO + H => C2H5 + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(1459, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCO + H => C2H5 + CO""",
)

entry(
    index = 562,
    label = "CH3CHCO + O => CH3CHO + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCO + O => CH3CHO + CO""",
)

entry(
    index = 563,
    label = "NC3H7 + HO2 <=> NC3H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7 + HO2 <=> NC3H7O + OH""",
)

entry(
    index = 564,
    label = "IC3H7 + HO2 <=> IC3H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7 + HO2 <=> IC3H7O + OH""",
)

entry(
    index = 565,
    label = "CH3O2 + NC3H7 <=> CH3O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + NC3H7 <=> CH3O + NC3H7O""",
)

entry(
    index = 566,
    label = "CH3O2 + IC3H7 <=> CH3O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + IC3H7 <=> CH3O + IC3H7O""",
)

entry(
    index = 567,
    label = "NC3H7O2 <=> NC3H7 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+20, 's^-1'), n=-1.616, Ea=(35960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 <=> NC3H7 + O2""",
)

entry(
    index = 568,
    label = "IC3H7O2 <=> IC3H7 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.132e+22, 's^-1'), n=-2.167, Ea=(38160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 <=> IC3H7 + O2""",
)

entry(
    index = 569,
    label = "NC3H7O2 + CH2O <=> NC3H7O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + CH2O <=> NC3H7O2H + HCO""",
)

entry(
    index = 570,
    label = "NC3H7O2 + CH3CHO <=> NC3H7O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + CH3CHO <=> NC3H7O2H + CH3CO""",
)

entry(
    index = 571,
    label = "IC3H7O2 + CH2O <=> IC3H7O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + CH2O <=> IC3H7O2H + HCO""",
)

entry(
    index = 572,
    label = "IC3H7O2 + CH3CHO <=> IC3H7O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + CH3CHO <=> IC3H7O2H + CH3CO""",
)

entry(
    index = 573,
    label = "NC3H7O2 + HO2 <=> NC3H7O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + HO2 <=> NC3H7O2H + O2""",
)

entry(
    index = 574,
    label = "IC3H7O2 + HO2 <=> IC3H7O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + HO2 <=> IC3H7O2H + O2""",
)

entry(
    index = 575,
    label = "C2H4 + NC3H7O2 <=> C2H3 + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + NC3H7O2 <=> C2H3 + NC3H7O2H""",
)

entry(
    index = 576,
    label = "C2H4 + IC3H7O2 <=> C2H3 + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + IC3H7O2 <=> C2H3 + IC3H7O2H""",
)

entry(
    index = 577,
    label = "CH3OH + NC3H7O2 <=> CH2OH + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + NC3H7O2 <=> CH2OH + NC3H7O2H""",
)

entry(
    index = 578,
    label = "CH3OH + IC3H7O2 <=> CH2OH + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + IC3H7O2 <=> CH2OH + IC3H7O2H""",
)

entry(
    index = 579,
    label = "C2H3CHO + NC3H7O2 <=> C2H3CO + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + NC3H7O2 <=> C2H3CO + NC3H7O2H""",
)

entry(
    index = 580,
    label = "C2H3CHO + IC3H7O2 <=> C2H3CO + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + IC3H7O2 <=> C2H3CO + IC3H7O2H""",
)

entry(
    index = 581,
    label = "CH4 + NC3H7O2 <=> CH3 + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(24640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + NC3H7O2 <=> CH3 + NC3H7O2H""",
)

entry(
    index = 582,
    label = "CH4 + IC3H7O2 <=> CH3 + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(24640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + IC3H7O2 <=> CH3 + IC3H7O2H""",
)

entry(
    index = 583,
    label = "NC3H7O2 + CH3O2 => NC3H7O + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + CH3O2 => NC3H7O + CH3O + O2""",
)

entry(
    index = 584,
    label = "IC3H7O2 + CH3O2 => IC3H7O + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + CH3O2 => IC3H7O + CH3O + O2""",
)

entry(
    index = 585,
    label = "H2 + NC3H7O2 <=> H + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + NC3H7O2 <=> H + NC3H7O2H""",
)

entry(
    index = 586,
    label = "H2 + IC3H7O2 <=> H + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + IC3H7O2 <=> H + IC3H7O2H""",
)

entry(
    index = 587,
    label = "IC3H7O2 + C2H6 <=> IC3H7O2H + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C2H6 <=> IC3H7O2H + C2H5""",
)

entry(
    index = 588,
    label = "NC3H7O2 + C2H6 <=> NC3H7O2H + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C2H6 <=> NC3H7O2H + C2H5""",
)

entry(
    index = 589,
    label = "IC3H7O2 + C2H5CHO <=> IC3H7O2H + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C2H5CHO <=> IC3H7O2H + C2H5CO""",
)

entry(
    index = 590,
    label = "NC3H7O2 + C2H5CHO <=> NC3H7O2H + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C2H5CHO <=> NC3H7O2H + C2H5CO""",
)

entry(
    index = 591,
    label = "IC3H7O2 + CH3CO3 => IC3H7O + CH3CO2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + CH3CO3 => IC3H7O + CH3CO2 + O2""",
)

entry(
    index = 592,
    label = "NC3H7O2 + CH3CO3 => NC3H7O + CH3CO2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + CH3CO3 => NC3H7O + CH3CO2 + O2""",
)

entry(
    index = 593,
    label = "IC3H7O2 + C2H5O2 => IC3H7O + C2H5O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C2H5O2 => IC3H7O + C2H5O + O2""",
)

entry(
    index = 594,
    label = "NC3H7O2 + C2H5O2 => NC3H7O + C2H5O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C2H5O2 => NC3H7O + C2H5O + O2""",
)

entry(
    index = 595,
    label = "IC3H7O2 + IC3H7O2 => O2 + IC3H7O + IC3H7O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + IC3H7O2 => O2 + IC3H7O + IC3H7O""",
)

entry(
    index = 596,
    label = "NC3H7O2 + NC3H7O2 => O2 + NC3H7O + NC3H7O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + NC3H7O2 => O2 + NC3H7O + NC3H7O""",
)

entry(
    index = 597,
    label = "IC3H7O2 + NC3H7O2 => IC3H7O + NC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + NC3H7O2 => IC3H7O + NC3H7O + O2""",
)

entry(
    index = 598,
    label = "IC3H7O2 + CH3 <=> IC3H7O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + CH3 <=> IC3H7O + CH3O""",
)

entry(
    index = 599,
    label = "IC3H7O2 + C2H5 <=> IC3H7O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C2H5 <=> IC3H7O + C2H5O""",
)

entry(
    index = 600,
    label = "IC3H7O2 + IC3H7 <=> IC3H7O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + IC3H7 <=> IC3H7O + IC3H7O""",
)

entry(
    index = 601,
    label = "IC3H7O2 + NC3H7 <=> IC3H7O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + NC3H7 <=> IC3H7O + NC3H7O""",
)

entry(
    index = 602,
    label = "IC3H7O2 + C3H5-A <=> IC3H7O + C3H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C3H5-A <=> IC3H7O + C3H5O""",
)

entry(
    index = 603,
    label = "NC3H7O2 + CH3 <=> NC3H7O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + CH3 <=> NC3H7O + CH3O""",
)

entry(
    index = 604,
    label = "NC3H7O2 + C2H5 <=> NC3H7O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C2H5 <=> NC3H7O + C2H5O""",
)

entry(
    index = 605,
    label = "NC3H7O2 + IC3H7 <=> NC3H7O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + IC3H7 <=> NC3H7O + IC3H7O""",
)

entry(
    index = 606,
    label = "NC3H7O2 + NC3H7 <=> NC3H7O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + NC3H7 <=> NC3H7O + NC3H7O""",
)

entry(
    index = 607,
    label = "NC3H7O2 + C3H5-A <=> NC3H7O + C3H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C3H5-A <=> NC3H7O + C3H5O""",
)

entry(
    index = 608,
    label = "NC3H7O2H <=> NC3H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2H <=> NC3H7O + OH""",
)

entry(
    index = 609,
    label = "IC3H7O2H <=> IC3H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.45e+15, 's^-1'), n=0, Ea=(42600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2H <=> IC3H7O + OH""",
)

entry(
    index = 610,
    label = "NC3H7O <=> C2H5 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.716e+21, 's^-1'), n=-2.449, Ea=(15700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O <=> C2H5 + CH2O""",
)

entry(
    index = 611,
    label = "NC3H7O <=> C2H5CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.899e+10, 's^-1'), n=0.746, Ea=(19800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O <=> C2H5CHO + H""",
)

entry(
    index = 612,
    label = "IC3H7O <=> CH3 + CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.328e+19, 's^-1'), n=-1.696, Ea=(17140, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O <=> CH3 + CH3CHO""",
)

entry(
    index = 613,
    label = "IC3H7O <=> CH3COCH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.663e+14, 's^-1'), n=-0.483, Ea=(20080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O <=> CH3COCH3 + H""",
)

entry(
    index = 614,
    label = "IC3H7O + O2 <=> CH3COCH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.09e+09, 'cm^3/(mol*s)'), n=0, Ea=(390, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O + O2 <=> CH3COCH3 + HO2""",
)

entry(
    index = 615,
    label = "NC3H7O2 <=> C3H6OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 <=> C3H6OOH1-2""",
)

entry(
    index = 616,
    label = "NC3H7O2 <=> C3H6OOH1-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.125e+11, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 <=> C3H6OOH1-3""",
)

entry(
    index = 617,
    label = "IC3H7O2 <=> C3H6OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.8e+12, 's^-1'), n=0, Ea=(29400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 <=> C3H6OOH2-1""",
)

entry(
    index = 618,
    label = "IC3H7O2 <=> C3H6OOH2-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.23e+35, 's^-1'), n=-6.96, Ea=(48880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 <=> C3H6OOH2-2""",
)

entry(
    index = 619,
    label = "C3H6OOH1-2 <=> C3H6O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-2 <=> C3H6O1-2 + OH""",
)

entry(
    index = 620,
    label = "C3H6OOH1-3 <=> C3H6O1-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-3 <=> C3H6O1-3 + OH""",
)

entry(
    index = 621,
    label = "C3H6OOH2-1 <=> C3H6O1-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-1 <=> C3H6O1-2 + OH""",
)

entry(
    index = 622,
    label = "C3H6OOH1-2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.834e+15, 's^-1'), n=-1.3, Ea=(15950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-2 <=> C3H6 + HO2""",
)

entry(
    index = 623,
    label = "C3H6OOH2-1 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.239e+18, 's^-1'), n=-2, Ea=(18970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-1 <=> C3H6 + HO2""",
)

entry(
    index = 624,
    label = "C3H6OOH1-3 => OH + CH2O + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.035e+15, 's^-1'), n=-0.79, Ea=(27400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-3 => OH + CH2O + C2H4""",
)

entry(
    index = 625,
    label = "C3H6OOH2-1 <=> C2H3OOH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.54e+27, 's^-1'), n=-5.14, Ea=(38320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-1 <=> C2H3OOH + CH3""",
)

entry(
    index = 626,
    label = "C3H6OOH1-2 => C2H4 + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.31e+33, 's^-1'), n=-7.01, Ea=(48120, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-2 => C2H4 + CH2O + OH""",
)

entry(
    index = 627,
    label = "C3H6OOH2-2 <=> CH3COCH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+14, 's^-1'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-2 <=> CH3COCH3 + OH""",
)

entry(
    index = 628,
    label = "C3H6OOH1-2O2 <=> C3H6OOH1-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.391e+25, 's^-1'), n=-2.945, Ea=(40100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-2O2 <=> C3H6OOH1-2 + O2""",
)

entry(
    index = 629,
    label = "C3H6OOH1-3O2 <=> C3H6OOH1-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.853e+20, 's^-1'), n=-1.626, Ea=(35690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-3O2 <=> C3H6OOH1-3 + O2""",
)

entry(
    index = 630,
    label = "C3H6OOH2-1O2 <=> C3H6OOH2-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.227e+22, 's^-1'), n=-2.244, Ea=(37820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-1O2 <=> C3H6OOH2-1 + O2""",
)

entry(
    index = 631,
    label = "C3H6OOH1-2O2 <=> C3KET12 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(26400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-2O2 <=> C3KET12 + OH""",
)

entry(
    index = 632,
    label = "C3H6OOH1-3O2 <=> C3KET13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(21400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-3O2 <=> C3KET13 + OH""",
)

entry(
    index = 633,
    label = "C3H6OOH2-1O2 <=> C3KET21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-1O2 <=> C3KET21 + OH""",
)

entry(
    index = 634,
    label = "C3H6OOH2-1O2 <=> C3H51-2,3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.125e+11, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH2-1O2 <=> C3H51-2,3OOH""",
)

entry(
    index = 635,
    label = "C3H6OOH1-2O2 <=> C3H51-2,3OOH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+11, 's^-1'), n=0, Ea=(29400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6OOH1-2O2 <=> C3H51-2,3OOH""",
)

entry(
    index = 636,
    label = "C3H51-2,3OOH <=> AC3H5OOH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.56e+13, 's^-1'), n=-0.49, Ea=(17770, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H51-2,3OOH <=> AC3H5OOH + HO2""",
)

entry(
    index = 637,
    label = "C3H52-1,3OOH <=> C3H6OOH1-3O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.255e+12, 's^-1'), n=-0.36, Ea=(13940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H52-1,3OOH <=> C3H6OOH1-3O2""",
)

entry(
    index = 638,
    label = "C3H52-1,3OOH <=> AC3H5OOH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+14, 's^-1'), n=-0.63, Ea=(17250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H52-1,3OOH <=> AC3H5OOH + HO2""",
)

entry(
    index = 639,
    label = "C3KET12 => CH3CHO + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.45e+15, 's^-1'), n=0, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3KET12 => CH3CHO + HCO + OH""",
)

entry(
    index = 640,
    label = "C3KET13 => CH2O + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3KET13 => CH2O + CH2CHO + OH""",
)

entry(
    index = 641,
    label = "C3KET21 => CH2O + CH3CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3KET21 => CH2O + CH3CO + OH""",
)

entry(
    index = 642,
    label = "AC3H5OOH <=> C3H5O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.88e+19, 's^-1'), n=-1.46, Ea=(45370, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5OOH <=> C3H5O + OH""",
)

entry(
    index = 643,
    label = "C3H5O <=> C2H3CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(29100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5O <=> C2H3CHO + H""",
)

entry(
    index = 644,
    label = "C3H5O <=> C2H3 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.464e+20, 's^-1'), n=-1.968, Ea=(35090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5O <=> C2H3 + CH2O""",
)

entry(
    index = 645,
    label = "C3H5O + O2 <=> C2H3CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5O + O2 <=> C2H3CHO + HO2""",
)

entry(
    index = 646,
    label = "C2H3OOH <=> CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.4e+14, 's^-1'), n=0, Ea=(43000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3OOH <=> CH2CHO + OH""",
)

entry(
    index = 647,
    label = "C3H6O1-2 <=> C2H4 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+14, 's^-1'), n=0, Ea=(60000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 <=> C2H4 + CH2O""",
)

entry(
    index = 648,
    label = "C3H6O1-2 + OH => CH2O + C2H3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 + OH => CH2O + C2H3 + H2O""",
)

entry(
    index = 649,
    label = "C3H6O1-2 + H => CH2O + C2H3 + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.63e+07, 'cm^3/(mol*s)'), n=2, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 + H => CH2O + C2H3 + H2""",
)

entry(
    index = 650,
    label = "C3H6O1-2 + O => CH2O + C2H3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.43e+13, 'cm^3/(mol*s)'), n=0, Ea=(5200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 + O => CH2O + C2H3 + OH""",
)

entry(
    index = 651,
    label = "C3H6O1-2 + HO2 => CH2O + C2H3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 + HO2 => CH2O + C2H3 + H2O2""",
)

entry(
    index = 652,
    label = "C3H6O1-2 + CH3O2 => CH2O + C2H3 + CH3O2H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 + CH3O2 => CH2O + C2H3 + CH3O2H""",
)

entry(
    index = 653,
    label = "C3H6O1-2 + CH3 => CH2O + C2H3 + CH4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-2 + CH3 => CH2O + C2H3 + CH4""",
)

entry(
    index = 654,
    label = "C3H6O1-3 <=> C2H4 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+14, 's^-1'), n=0, Ea=(60000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 <=> C2H4 + CH2O""",
)

entry(
    index = 655,
    label = "C3H6O1-3 + OH => CH2O + C2H3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 + OH => CH2O + C2H3 + H2O""",
)

entry(
    index = 656,
    label = "C3H6O1-3 + O => CH2O + C2H3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.43e+13, 'cm^3/(mol*s)'), n=0, Ea=(5200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 + O => CH2O + C2H3 + OH""",
)

entry(
    index = 657,
    label = "C3H6O1-3 + H => CH2O + C2H3 + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.63e+07, 'cm^3/(mol*s)'), n=2, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 + H => CH2O + C2H3 + H2""",
)

entry(
    index = 658,
    label = "C3H6O1-3 + CH3O2 => CH2O + C2H3 + CH3O2H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 + CH3O2 => CH2O + C2H3 + CH3O2H""",
)

entry(
    index = 659,
    label = "C3H6O1-3 + HO2 => CH2O + C2H3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 + HO2 => CH2O + C2H3 + H2O2""",
)

entry(
    index = 660,
    label = "C3H6O1-3 + CH3 => CH2O + C2H3 + CH4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6O1-3 + CH3 => CH2O + C2H3 + CH4""",
)

entry(
    index = 661,
    label = "IC3H7O2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.015e+43, 's^-1'), n=-9.409, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 <=> C3H6 + HO2""",
)

entry(
    index = 662,
    label = "NC3H7O2 <=> C3H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.112, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 <=> C3H6 + HO2""",
)

entry(
    index = 663,
    label = "C4H10 <=> C2H5 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(2.72e+15, 's^-1'), n=0, Ea=(75610, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(4.72e+18, 'cm^3/(mol*s)'), n=0, Ea=(49576, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.72,
        T3 = (1500, 'K'),
        T1 = (1e-10, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C4H10 <=> C2H5 + C2H5""",
)

entry(
    index = 664,
    label = "C4H10 <=> NC3H7 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.28e+14, 's^-1'), n=0, Ea=(69900, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(5.34e+17, 'cm^3/(mol*s)'), n=0, Ea=(42959, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.72,
        T3 = (1500, 'K'),
        T1 = (1e-10, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is C4H10 <=> NC3H7 + CH3""",
)

entry(
    index = 665,
    label = "C4H10 <=> PC4H9 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.344e+17, 's^-1'), n=-0.356, Ea=(101200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 <=> PC4H9 + H""",
)

entry(
    index = 666,
    label = "C4H10 <=> SC4H9 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.975e+18, 's^-1'), n=-0.694, Ea=(98720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 <=> SC4H9 + H""",
)

entry(
    index = 667,
    label = "C4H10 + O2 <=> PC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(52340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O2 <=> PC4H9 + HO2""",
)

entry(
    index = 668,
    label = "C4H10 + O2 <=> SC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(49800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O2 <=> SC4H9 + HO2""",
)

entry(
    index = 669,
    label = "C4H10 + C3H5-A <=> PC4H9 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(20500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C3H5-A <=> PC4H9 + C3H6""",
)

entry(
    index = 670,
    label = "C4H10 + C3H5-A <=> SC4H9 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.16e+11, 'cm^3/(mol*s)'), n=0, Ea=(16400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C3H5-A <=> SC4H9 + C3H6""",
)

entry(
    index = 671,
    label = "C4H10 + C2H5 <=> PC4H9 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.58e+11, 'cm^3/(mol*s)'), n=0, Ea=(12300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C2H5 <=> PC4H9 + C2H6""",
)

entry(
    index = 672,
    label = "C4H10 + C2H5 <=> SC4H9 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C2H5 <=> SC4H9 + C2H6""",
)

entry(
    index = 673,
    label = "C4H10 + C2H3 <=> PC4H9 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(18000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C2H3 <=> PC4H9 + C2H4""",
)

entry(
    index = 674,
    label = "C4H10 + C2H3 <=> SC4H9 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C2H3 <=> SC4H9 + C2H4""",
)

entry(
    index = 675,
    label = "C4H10 + CH3 <=> PC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.904, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3 <=> PC4H9 + CH4""",
)

entry(
    index = 676,
    label = "C4H10 + CH3 <=> SC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.02, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3 <=> SC4H9 + CH4""",
)

entry(
    index = 677,
    label = "C4H10 + H <=> PC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(188000, 'cm^3/(mol*s)'), n=2.75, Ea=(6280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + H <=> PC4H9 + H2""",
)

entry(
    index = 678,
    label = "C4H10 + H <=> SC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + H <=> SC4H9 + H2""",
)

entry(
    index = 679,
    label = "C4H10 + OH <=> PC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H10 + OH <=> PC4H9 + H2O""",
)

entry(
    index = 680,
    label = "C4H10 + OH <=> SC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.34e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H10 + OH <=> SC4H9 + H2O""",
)

entry(
    index = 681,
    label = "C4H10 + O <=> PC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+14, 'cm^3/(mol*s)'), n=0, Ea=(7850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O <=> PC4H9 + OH""",
)

entry(
    index = 682,
    label = "C4H10 + O <=> SC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.62e+13, 'cm^3/(mol*s)'), n=0, Ea=(5200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O <=> SC4H9 + OH""",
)

entry(
    index = 683,
    label = "C4H10 + HO2 <=> PC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40.8, 'cm^3/(mol*s)'), n=3.59, Ea=(17160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + HO2 <=> PC4H9 + H2O2""",
)

entry(
    index = 684,
    label = "C4H10 + HO2 <=> SC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(126.4, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + HO2 <=> SC4H9 + H2O2""",
)

entry(
    index = 685,
    label = "C4H10 + CH3O <=> PC4H9 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3O <=> PC4H9 + CH3OH""",
)

entry(
    index = 686,
    label = "C4H10 + CH3O <=> SC4H9 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3O <=> SC4H9 + CH3OH""",
)

entry(
    index = 687,
    label = "C4H10 + C2H5O <=> PC4H9 + C2H5OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C2H5O <=> PC4H9 + C2H5OH""",
)

entry(
    index = 688,
    label = "C4H10 + C2H5O <=> SC4H9 + C2H5OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + C2H5O <=> SC4H9 + C2H5OH""",
)

entry(
    index = 689,
    label = "C4H10 + PC4H9 <=> SC4H9 + C4H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + PC4H9 <=> SC4H9 + C4H10""",
)

entry(
    index = 690,
    label = "C4H10 + CH3CO3 <=> PC4H9 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3CO3 <=> PC4H9 + CH3CO3H""",
)

entry(
    index = 691,
    label = "C4H10 + CH3CO3 <=> SC4H9 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + CH3CO3 <=> SC4H9 + CH3CO3H""",
)

entry(
    index = 692,
    label = "C4H10 + O2CHO <=> PC4H9 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.68e+13, 'cm^3/(mol*s)'), n=0, Ea=(20440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O2CHO <=> PC4H9 + HO2CHO""",
)

entry(
    index = 693,
    label = "C4H10 + O2CHO <=> SC4H9 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H10 + O2CHO <=> SC4H9 + HO2CHO""",
)

entry(
    index = 694,
    label = "CH3O2 + C4H10 <=> CH3O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.386, 'cm^3/(mol*s)'), n=3.97, Ea=(18280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + C4H10 <=> CH3O2H + PC4H9""",
)

entry(
    index = 695,
    label = "CH3O2 + C4H10 <=> CH3O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20.37, 'cm^3/(mol*s)'), n=3.58, Ea=(14810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + C4H10 <=> CH3O2H + SC4H9""",
)

entry(
    index = 696,
    label = "C2H5O2 + C4H10 <=> C2H5O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40.8, 'cm^3/(mol*s)'), n=3.59, Ea=(17160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 + C4H10 <=> C2H5O2H + PC4H9""",
)

entry(
    index = 697,
    label = "C2H5O2 + C4H10 <=> C2H5O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(126.4, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5O2 + C4H10 <=> C2H5O2H + SC4H9""",
)

entry(
    index = 698,
    label = "NC3H7O2 + C4H10 <=> NC3H7O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C4H10 <=> NC3H7O2H + PC4H9""",
)

entry(
    index = 699,
    label = "NC3H7O2 + C4H10 <=> NC3H7O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C4H10 <=> NC3H7O2H + SC4H9""",
)

entry(
    index = 700,
    label = "IC3H7O2 + C4H10 <=> IC3H7O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C4H10 <=> IC3H7O2H + PC4H9""",
)

entry(
    index = 701,
    label = "IC3H7O2 + C4H10 <=> IC3H7O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C4H10 <=> IC3H7O2H + SC4H9""",
)

entry(
    index = 702,
    label = "PC4H9O2 + C3H8 <=> PC4H9O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C3H8 <=> PC4H9O2H + NC3H7""",
)

entry(
    index = 703,
    label = "PC4H9O2 + C3H8 <=> PC4H9O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C3H8 <=> PC4H9O2H + IC3H7""",
)

entry(
    index = 704,
    label = "PC4H9O2 + C4H10 <=> PC4H9O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C4H10 <=> PC4H9O2H + PC4H9""",
)

entry(
    index = 705,
    label = "PC4H9O2 + C4H10 <=> PC4H9O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C4H10 <=> PC4H9O2H + SC4H9""",
)

entry(
    index = 706,
    label = "SC4H9O2 + C3H8 <=> SC4H9O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C3H8 <=> SC4H9O2H + NC3H7""",
)

entry(
    index = 707,
    label = "SC4H9O2 + C3H8 <=> SC4H9O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C3H8 <=> SC4H9O2H + IC3H7""",
)

entry(
    index = 708,
    label = "SC4H9O2 + C4H10 <=> SC4H9O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C4H10 <=> SC4H9O2H + PC4H9""",
)

entry(
    index = 709,
    label = "SC4H9O2 + C4H10 <=> SC4H9O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C4H10 <=> SC4H9O2H + SC4H9""",
)

entry(
    index = 710,
    label = "PC4H9 <=> C2H5 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.504e+12, 's^-1'), n=0.463, Ea=(29470, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9 <=> C2H5 + C2H4""",
)

entry(
    index = 711,
    label = "SC4H9 <=> C3H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.803e+10, 's^-1'), n=1.044, Ea=(30350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9 <=> C3H6 + CH3""",
)

entry(
    index = 712,
    label = "PC4H9 <=> C4H8-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.622e+12, 's^-1'), n=0.253, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9 <=> C4H8-1 + H""",
)

entry(
    index = 713,
    label = "SC4H9 <=> C4H8-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.844e+11, 's^-1'), n=0.337, Ea=(35520, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9 <=> C4H8-2 + H""",
)

entry(
    index = 714,
    label = "SC4H9 <=> C4H8-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.027e+11, 's^-1'), n=0.591, Ea=(36820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9 <=> C4H8-1 + H""",
)

entry(
    index = 715,
    label = "PC4H9 + O2 <=> C4H8-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e-18, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9 + O2 <=> C4H8-1 + HO2""",
)

entry(
    index = 716,
    label = "SC4H9 + O2 <=> C4H8-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e-18, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9 + O2 <=> C4H8-1 + HO2""",
)

entry(
    index = 717,
    label = "SC4H9 + O2 <=> C4H8-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e-18, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9 + O2 <=> C4H8-2 + HO2""",
)

entry(
    index = 718,
    label = "C4H8-1 <=> C3H5-A + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.081e+19, 's^-1'), n=-1.256, Ea=(76510, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 <=> C3H5-A + CH3""",
)

entry(
    index = 719,
    label = "C4H8-1 <=> C2H3 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.882e+23, 's^-1'), n=-1.99, Ea=(101600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 <=> C2H3 + C2H5""",
)

entry(
    index = 720,
    label = "C4H8-1 <=> H + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.724e+14, 's^-1'), n=-0.111, Ea=(85200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 <=> H + C4H71-3""",
)

entry(
    index = 721,
    label = "C4H8-1 + O2 <=> C4H71-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(37190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + O2 <=> C4H71-3 + HO2""",
)

entry(
    index = 722,
    label = "C4H8-1 + H <=> C4H71-1 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(781000, 'cm^3/(mol*s)'), n=2.5, Ea=(12290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + H <=> C4H71-1 + H2""",
)

entry(
    index = 723,
    label = "C4H8-1 + H <=> C4H71-2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(390000, 'cm^3/(mol*s)'), n=2.5, Ea=(5821, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + H <=> C4H71-2 + H2""",
)

entry(
    index = 724,
    label = "C4H8-1 + H <=> C4H71-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + H <=> C4H71-3 + H2""",
)

entry(
    index = 725,
    label = "C4H8-1 + H <=> C4H71-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665100, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + H <=> C4H71-4 + H2""",
)

entry(
    index = 726,
    label = "C4H8-1 + OH <=> C4H71-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.14e+06, 'cm^3/(mol*s)'), n=2, Ea=(2778, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + OH <=> C4H71-1 + H2O""",
)

entry(
    index = 727,
    label = "C4H8-1 + OH <=> C4H71-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.22e+06, 'cm^3/(mol*s)'), n=2, Ea=(1451, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + OH <=> C4H71-2 + H2O""",
)

entry(
    index = 728,
    label = "C4H8-1 + OH <=> C4H71-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + OH <=> C4H71-3 + H2O""",
)

entry(
    index = 729,
    label = "C4H8-1 + OH <=> C4H71-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + OH <=> C4H71-4 + H2O""",
)

entry(
    index = 730,
    label = "C4H8-1 + CH3 <=> C4H71-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3 <=> C4H71-3 + CH4""",
)

entry(
    index = 731,
    label = "C4H8-1 + CH3 <=> C4H71-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.452, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3 <=> C4H71-4 + CH4""",
)

entry(
    index = 732,
    label = "C4H8-1 + HO2 <=> C4H71-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + HO2 <=> C4H71-3 + H2O2""",
)

entry(
    index = 733,
    label = "C4H8-1 + HO2 <=> C4H71-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2380, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + HO2 <=> C4H71-4 + H2O2""",
)

entry(
    index = 734,
    label = "C4H8-1 + CH3O2 <=> C4H71-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3O2 <=> C4H71-3 + CH3O2H""",
)

entry(
    index = 735,
    label = "C4H8-1 + CH3O2 <=> C4H71-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2380, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3O2 <=> C4H71-4 + CH3O2H""",
)

entry(
    index = 736,
    label = "C4H8-1 + CH3O <=> C4H71-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3O <=> C4H71-3 + CH3OH""",
)

entry(
    index = 737,
    label = "C4H8-1 + CH3O <=> C4H71-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3O <=> C4H71-4 + CH3OH""",
)

entry(
    index = 738,
    label = "C4H8-1 + CH3CO3 <=> C4H71-3 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(8000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3CO3 <=> C4H71-3 + CH3CO3H""",
)

entry(
    index = 739,
    label = "C4H8-1 + C3H5-A <=> C4H71-3 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.9e+10, 'cm^3/(mol*s)'), n=0, Ea=(12400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + C3H5-A <=> C4H71-3 + C3H6""",
)

entry(
    index = 740,
    label = "C4H8-1 + C4H6 <=> C4H71-3 + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.35e+12, 'cm^3/(mol*s)'), n=0, Ea=(46720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + C4H6 <=> C4H71-3 + C4H71-3""",
)

entry(
    index = 741,
    label = "C4H8-1 + C2H5O2 <=> C4H71-3 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + C2H5O2 <=> C4H71-3 + C2H5O2H""",
)

entry(
    index = 742,
    label = "C4H8-1 + NC3H7O2 <=> C4H71-3 + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + NC3H7O2 <=> C4H71-3 + NC3H7O2H""",
)

entry(
    index = 743,
    label = "C4H8-1 + IC3H7O2 <=> C4H71-3 + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + IC3H7O2 <=> C4H71-3 + IC3H7O2H""",
)

entry(
    index = 744,
    label = "C4H8-1 + PC4H9O2 <=> C4H71-3 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + PC4H9O2 <=> C4H71-3 + PC4H9O2H""",
)

entry(
    index = 745,
    label = "C4H8-1 + SC4H9O2 <=> C4H71-3 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + SC4H9O2 <=> C4H71-3 + SC4H9O2H""",
)

entry(
    index = 746,
    label = "C4H8-1 + CH3O2 => C4H8O1-2 + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + CH3O2 => C4H8O1-2 + CH3O""",
)

entry(
    index = 747,
    label = "C4H8-2 <=> H + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.337e+14, 's^-1'), n=0.143, Ea=(87890, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 <=> H + C4H71-3""",
)

entry(
    index = 748,
    label = "C4H8-2 + O2 <=> C4H71-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(39390, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + O2 <=> C4H71-3 + HO2""",
)

entry(
    index = 749,
    label = "C4H8-2 + H <=> C4H71-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(346000, 'cm^3/(mol*s)'), n=2.5, Ea=(2492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + H <=> C4H71-3 + H2""",
)

entry(
    index = 750,
    label = "C4H8-2 + OH <=> C4H71-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.24e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + OH <=> C4H71-3 + H2O""",
)

entry(
    index = 751,
    label = "C4H8-2 + CH3 <=> C4H71-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.42, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + CH3 <=> C4H71-3 + CH4""",
)

entry(
    index = 752,
    label = "C4H8-2 + HO2 <=> C4H71-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19280, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + HO2 <=> C4H71-3 + H2O2""",
)

entry(
    index = 753,
    label = "C4H8-2 + CH3O2 <=> C4H71-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19280, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + CH3O2 <=> C4H71-3 + CH3O2H""",
)

entry(
    index = 754,
    label = "C4H8-2 + CH3O <=> C4H71-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(18, 'cm^3/(mol*s)'), n=2.95, Ea=(11990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + CH3O <=> C4H71-3 + CH3OH""",
)

entry(
    index = 755,
    label = "C4H8-2 + C2H5O2 <=> C4H71-3 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + C2H5O2 <=> C4H71-3 + C2H5O2H""",
)

entry(
    index = 756,
    label = "C4H8-2 + NC3H7O2 <=> C4H71-3 + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + NC3H7O2 <=> C4H71-3 + NC3H7O2H""",
)

entry(
    index = 757,
    label = "C4H8-2 + IC3H7O2 <=> C4H71-3 + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + IC3H7O2 <=> C4H71-3 + IC3H7O2H""",
)

entry(
    index = 758,
    label = "C4H8-2 + PC4H9O2 <=> C4H71-3 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + PC4H9O2 <=> C4H71-3 + PC4H9O2H""",
)

entry(
    index = 759,
    label = "C4H8-2 + SC4H9O2 <=> C4H71-3 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + SC4H9O2 <=> C4H71-3 + SC4H9O2H""",
)

entry(
    index = 760,
    label = "C4H8-1 + HO2 => C4H8O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(14340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-1 + HO2 => C4H8O1-2 + OH""",
)

entry(
    index = 761,
    label = "C4H8-2 + HO2 => C4H8O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.62e+11, 'cm^3/(mol*s)'), n=0, Ea=(12310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + HO2 => C4H8O2-3 + OH""",
)

entry(
    index = 762,
    label = "C4H8-2 + CH3O2 => C4H8O2-3 + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.62e+11, 'cm^3/(mol*s)'), n=0, Ea=(12310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8-2 + CH3O2 => C4H8O2-3 + CH3O""",
)

entry(
    index = 763,
    label = "PC4H8OH <=> C4H8-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.079e+16, 's^-1'), n=-0.699, Ea=(28090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H8OH <=> C4H8-1 + OH""",
)

entry(
    index = 764,
    label = "SC4H8OH <=> C4H8-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.379e+17, 's^-1'), n=-1.253, Ea=(29920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H8OH <=> C4H8-2 + OH""",
)

entry(
    index = 765,
    label = "C4H8OH-1O2 <=> PC4H8OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.753e+20, 's^-1'), n=-1.944, Ea=(35520, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OH-1O2 <=> PC4H8OH + O2""",
)

entry(
    index = 766,
    label = "C4H8OH-2O2 <=> SC4H8OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.686e+20, 's^-1'), n=-1.968, Ea=(35510, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OH-2O2 <=> SC4H8OH + O2""",
)

entry(
    index = 767,
    label = "C4H8OH-1O2 => C2H5CHO + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OH-1O2 => C2H5CHO + CH2O + OH""",
)

entry(
    index = 768,
    label = "C4H8OH-2O2 => OH + CH3CHO + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(25000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OH-2O2 => OH + CH3CHO + CH3CHO""",
)

entry(
    index = 769,
    label = "C4H71-1 <=> C2H2 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.066e+15, 's^-1'), n=-0.56, Ea=(30320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-1 <=> C2H2 + C2H5""",
)

entry(
    index = 770,
    label = "C4H71-2 <=> C3H4-A + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.589e+14, 's^-1'), n=-0.71, Ea=(31260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-2 <=> C3H4-A + CH3""",
)

entry(
    index = 771,
    label = "C4H71-4 <=> C2H4 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.766e+12, 's^-1'), n=-0.22, Ea=(36290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-4 <=> C2H4 + C2H3""",
)

entry(
    index = 772,
    label = "C4H72-2 <=> C3H4-P + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.326e+10, 's^-1'), n=0.52, Ea=(30020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H72-2 <=> C3H4-P + CH3""",
)

entry(
    index = 773,
    label = "C4H71-3 <=> C4H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.2e+14, 's^-1'), n=0, Ea=(49300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 <=> C4H6 + H""",
)

entry(
    index = 774,
    label = "C4H71-3 + C2H5 <=> C4H8-1 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.59e+12, 'cm^3/(mol*s)'), n=0, Ea=(-131, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + C2H5 <=> C4H8-1 + C2H4""",
)

entry(
    index = 775,
    label = "C4H71-3 + CH3O <=> C4H8-1 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.41e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + CH3O <=> C4H8-1 + CH2O""",
)

entry(
    index = 776,
    label = "C4H71-3 + O <=> C2H3CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + O <=> C2H3CHO + CH3""",
)

entry(
    index = 777,
    label = "C4H71-3 + HO2 <=> C4H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + HO2 <=> C4H7O + OH""",
)

entry(
    index = 778,
    label = "C4H71-3 + CH3O2 <=> C4H7O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + CH3O2 <=> C4H7O + CH3O""",
)

entry(
    index = 779,
    label = "C3H5-A + C4H71-3 <=> C3H6 + C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.31e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + C4H71-3 <=> C3H6 + C4H6""",
)

entry(
    index = 780,
    label = "C4H71-3 + O2 <=> C4H6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+09, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + O2 <=> C4H6 + HO2""",
)

entry(
    index = 781,
    label = "H + C4H71-3 <=> C4H6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.16e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H + C4H71-3 <=> C4H6 + H2""",
)

entry(
    index = 782,
    label = "C2H5 + C4H71-3 <=> C4H6 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C4H71-3 <=> C4H6 + C2H6""",
)

entry(
    index = 783,
    label = "C2H3 + C4H71-3 <=> C2H4 + C4H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C4H71-3 <=> C2H4 + C4H6""",
)

entry(
    index = 784,
    label = "C4H71-3 + C2H5O2 <=> C4H7O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + C2H5O2 <=> C4H7O + C2H5O""",
)

entry(
    index = 785,
    label = "IC3H7O2 + C4H71-3 <=> IC3H7O + C4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + C4H71-3 <=> IC3H7O + C4H7O""",
)

entry(
    index = 786,
    label = "NC3H7O2 + C4H71-3 <=> NC3H7O + C4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + C4H71-3 <=> NC3H7O + C4H7O""",
)

entry(
    index = 787,
    label = "C4H7O <=> CH3CHO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+14, 's^-1'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7O <=> CH3CHO + C2H3""",
)

entry(
    index = 788,
    label = "C4H7O <=> C2H3CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+14, 's^-1'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7O <=> C2H3CHO + CH3""",
)

entry(
    index = 789,
    label = "C4H6 <=> C2H3 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.027e+19, 's^-1'), n=-1, Ea=(98150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 <=> C2H3 + C2H3""",
)

entry(
    index = 790,
    label = "C4H6 + OH <=> C2H5 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + OH <=> C2H5 + CH2CO""",
)

entry(
    index = 791,
    label = "C4H6 + OH <=> CH2O + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + OH <=> CH2O + C3H5-A""",
)

entry(
    index = 792,
    label = "C4H6 + OH <=> C2H3 + CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + OH <=> C2H3 + CH3CHO""",
)

entry(
    index = 793,
    label = "C4H6 + O <=> C2H4 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + O <=> C2H4 + CH2CO""",
)

entry(
    index = 794,
    label = "C4H6 + O <=> CH2O + C3H4-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6 + O <=> CH2O + C3H4-A""",
)

entry(
    index = 795,
    label = "C2H3 + C2H4 <=> C4H6 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+11, 'cm^3/(mol*s)'), n=0, Ea=(7300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C2H4 <=> C4H6 + H""",
)

entry(
    index = 796,
    label = "C4H8O1-2 + OH => CH2O + C3H5-A + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-2 + OH => CH2O + C3H5-A + H2O""",
)

entry(
    index = 797,
    label = "C4H8O1-2 + H => CH2O + C3H5-A + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-2 + H => CH2O + C3H5-A + H2""",
)

entry(
    index = 798,
    label = "C4H8O1-2 + O => CH2O + C3H5-A + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-2 + O => CH2O + C3H5-A + OH""",
)

entry(
    index = 799,
    label = "C4H8O1-2 + HO2 => CH2O + C3H5-A + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-2 + HO2 => CH2O + C3H5-A + H2O2""",
)

entry(
    index = 800,
    label = "C4H8O1-2 + CH3O2 => CH2O + C3H5-A + CH3O2H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-2 + CH3O2 => CH2O + C3H5-A + CH3O2H""",
)

entry(
    index = 801,
    label = "C4H8O1-2 + CH3 => CH2O + C3H5-A + CH4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-2 + CH3 => CH2O + C3H5-A + CH4""",
)

entry(
    index = 802,
    label = "C4H8O1-3 + OH => CH2O + C3H5-A + H2O",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H8O1-3 + OH => CH2O + C3H5-A + H2O""",
)

entry(
    index = 803,
    label = "C4H8O1-3 + H => CH2O + C3H5-A + H2",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(3.51e+07, 'cm^3/(mol*s)'), n=2, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H8O1-3 + H => CH2O + C3H5-A + H2""",
)

entry(
    index = 804,
    label = "C4H8O1-3 + O => CH2O + C3H5-A + OH",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1.124e+14, 'cm^3/(mol*s)'), n=0, Ea=(5200, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H8O1-3 + O => CH2O + C3H5-A + OH""",
)

entry(
    index = 805,
    label = "C4H8O1-3 + HO2 => CH2O + C3H5-A + H2O2",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H8O1-3 + HO2 => CH2O + C3H5-A + H2O2""",
)

entry(
    index = 806,
    label = "C4H8O1-3 + CH3O2 => CH2O + C3H5-A + CH3O2H",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H8O1-3 + CH3O2 => CH2O + C3H5-A + CH3O2H""",
)

entry(
    index = 807,
    label = "C4H8O1-3 + CH3 => CH2O + C3H5-A + CH4",
    degeneracy = 1,
    duplicate = True,
    reversible = False,
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
            Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
        ],
    ),
    shortDesc = u"""The chemkin file reaction is C4H8O1-3 + CH3 => CH2O + C3H5-A + CH4""",
)

entry(
    index = 808,
    label = "C4H8O1-4 + OH => CH2O + C3H5-A + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-4 + OH => CH2O + C3H5-A + H2O""",
)

entry(
    index = 809,
    label = "C4H8O1-4 + H => CH2O + C3H5-A + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-4 + H => CH2O + C3H5-A + H2""",
)

entry(
    index = 810,
    label = "C4H8O1-4 + O => CH2O + C3H5-A + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-4 + O => CH2O + C3H5-A + OH""",
)

entry(
    index = 811,
    label = "C4H8O1-4 + HO2 => CH2O + C3H5-A + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-4 + HO2 => CH2O + C3H5-A + H2O2""",
)

entry(
    index = 812,
    label = "C4H8O1-4 + CH3O2 => CH2O + C3H5-A + CH3O2H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-4 + CH3O2 => CH2O + C3H5-A + CH3O2H""",
)

entry(
    index = 813,
    label = "C4H8O1-4 + CH3 => CH2O + C3H5-A + CH4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O1-4 + CH3 => CH2O + C3H5-A + CH4""",
)

entry(
    index = 814,
    label = "C4H8O2-3 + OH => CH2O + C3H5-A + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O2-3 + OH => CH2O + C3H5-A + H2O""",
)

entry(
    index = 815,
    label = "C4H8O2-3 + H => CH2O + C3H5-A + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O2-3 + H => CH2O + C3H5-A + H2""",
)

entry(
    index = 816,
    label = "C4H8O2-3 + O => CH2O + C3H5-A + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O2-3 + O => CH2O + C3H5-A + OH""",
)

entry(
    index = 817,
    label = "C4H8O2-3 + HO2 => CH2O + C3H5-A + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O2-3 + HO2 => CH2O + C3H5-A + H2O2""",
)

entry(
    index = 818,
    label = "C4H8O2-3 + CH3O2 => CH2O + C3H5-A + CH3O2H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O2-3 + CH3O2 => CH2O + C3H5-A + CH3O2H""",
)

entry(
    index = 819,
    label = "C4H8O2-3 + CH3 => CH2O + C3H5-A + CH4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8O2-3 + CH3 => CH2O + C3H5-A + CH4""",
)

entry(
    index = 820,
    label = "PC4H9O2 <=> PC4H9 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.849e+20, 's^-1'), n=-1.642, Ea=(35930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 <=> PC4H9 + O2""",
)

entry(
    index = 821,
    label = "SC4H9O2 <=> SC4H9 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.329e+22, 's^-1'), n=-2.216, Ea=(38160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 <=> SC4H9 + O2""",
)

entry(
    index = 822,
    label = "SC4H9O2 + CH2O <=> SC4H9O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + CH2O <=> SC4H9O2H + HCO""",
)

entry(
    index = 823,
    label = "SC4H9O2 + CH3CHO <=> SC4H9O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + CH3CHO <=> SC4H9O2H + CH3CO""",
)

entry(
    index = 824,
    label = "SC4H9O2 + HO2 <=> SC4H9O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + HO2 <=> SC4H9O2H + O2""",
)

entry(
    index = 825,
    label = "IC3H7O2 + PC4H9 <=> IC3H7O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + PC4H9 <=> IC3H7O + PC4H9O""",
)

entry(
    index = 826,
    label = "IC3H7O2 + SC4H9 <=> IC3H7O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + SC4H9 <=> IC3H7O + SC4H9O""",
)

entry(
    index = 827,
    label = "NC3H7O2 + PC4H9 <=> NC3H7O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + PC4H9 <=> NC3H7O + PC4H9O""",
)

entry(
    index = 828,
    label = "NC3H7O2 + SC4H9 <=> NC3H7O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + SC4H9 <=> NC3H7O + SC4H9O""",
)

entry(
    index = 829,
    label = "SC4H9O2 + SC4H9O2 => O2 + SC4H9O + SC4H9O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + SC4H9O2 => O2 + SC4H9O + SC4H9O""",
)

entry(
    index = 830,
    label = "SC4H9O2 + NC3H7O2 => SC4H9O + NC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + NC3H7O2 => SC4H9O + NC3H7O + O2""",
)

entry(
    index = 831,
    label = "SC4H9O2 + IC3H7O2 => SC4H9O + IC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + IC3H7O2 => SC4H9O + IC3H7O + O2""",
)

entry(
    index = 832,
    label = "SC4H9O2 + C2H5O2 => SC4H9O + C2H5O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C2H5O2 => SC4H9O + C2H5O + O2""",
)

entry(
    index = 833,
    label = "SC4H9O2 + CH3O2 => SC4H9O + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + CH3O2 => SC4H9O + CH3O + O2""",
)

entry(
    index = 834,
    label = "SC4H9O2 + CH3CO3 => SC4H9O + CH3CO2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + CH3CO3 => SC4H9O + CH3CO2 + O2""",
)

entry(
    index = 835,
    label = "PC4H9O2 + HO2 => PC4H9O + OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e-14, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + HO2 => PC4H9O + OH + O2""",
)

entry(
    index = 836,
    label = "SC4H9O2 + HO2 => SC4H9O + OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e-14, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + HO2 => SC4H9O + OH + O2""",
)

entry(
    index = 837,
    label = "H2 + PC4H9O2 <=> H + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + PC4H9O2 <=> H + PC4H9O2H""",
)

entry(
    index = 838,
    label = "H2 + SC4H9O2 <=> H + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + SC4H9O2 <=> H + SC4H9O2H""",
)

entry(
    index = 839,
    label = "C2H6 + PC4H9O2 <=> C2H5 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + PC4H9O2 <=> C2H5 + PC4H9O2H""",
)

entry(
    index = 840,
    label = "C2H6 + SC4H9O2 <=> C2H5 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H6 + SC4H9O2 <=> C2H5 + SC4H9O2H""",
)

entry(
    index = 841,
    label = "PC4H9O2 + C2H5CHO <=> PC4H9O2H + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C2H5CHO <=> PC4H9O2H + C2H5CO""",
)

entry(
    index = 842,
    label = "SC4H9O2 + C2H5CHO <=> SC4H9O2H + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C2H5CHO <=> SC4H9O2H + C2H5CO""",
)

entry(
    index = 843,
    label = "SC4H9O2 + CH3 <=> SC4H9O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + CH3 <=> SC4H9O + CH3O""",
)

entry(
    index = 844,
    label = "SC4H9O2 + C2H5 <=> SC4H9O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C2H5 <=> SC4H9O + C2H5O""",
)

entry(
    index = 845,
    label = "SC4H9O2 + IC3H7 <=> SC4H9O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + IC3H7 <=> SC4H9O + IC3H7O""",
)

entry(
    index = 846,
    label = "SC4H9O2 + NC3H7 <=> SC4H9O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + NC3H7 <=> SC4H9O + NC3H7O""",
)

entry(
    index = 847,
    label = "SC4H9O2 + PC4H9 <=> SC4H9O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + PC4H9 <=> SC4H9O + PC4H9O""",
)

entry(
    index = 848,
    label = "SC4H9O2 + SC4H9 <=> SC4H9O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + SC4H9 <=> SC4H9O + SC4H9O""",
)

entry(
    index = 849,
    label = "SC4H9O2 + C3H5-A <=> SC4H9O + C3H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + C3H5-A <=> SC4H9O + C3H5O""",
)

entry(
    index = 850,
    label = "PC4H9O2 + CH2O <=> PC4H9O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + CH2O <=> PC4H9O2H + HCO""",
)

entry(
    index = 851,
    label = "PC4H9O2 + CH3CHO <=> PC4H9O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + CH3CHO <=> PC4H9O2H + CH3CO""",
)

entry(
    index = 852,
    label = "PC4H9O2 + HO2 <=> PC4H9O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + HO2 <=> PC4H9O2H + O2""",
)

entry(
    index = 853,
    label = "C3H6 + PC4H9O2 <=> C3H5-A + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + PC4H9O2 <=> C3H5-A + PC4H9O2H""",
)

entry(
    index = 854,
    label = "C3H6 + SC4H9O2 <=> C3H5-A + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6 + SC4H9O2 <=> C3H5-A + SC4H9O2H""",
)

entry(
    index = 855,
    label = "C2H4 + PC4H9O2 <=> C2H3 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + PC4H9O2 <=> C2H3 + PC4H9O2H""",
)

entry(
    index = 856,
    label = "C2H4 + SC4H9O2 <=> C2H3 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + SC4H9O2 <=> C2H3 + SC4H9O2H""",
)

entry(
    index = 857,
    label = "CH3OH + PC4H9O2 <=> CH2OH + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + PC4H9O2 <=> CH2OH + PC4H9O2H""",
)

entry(
    index = 858,
    label = "CH3OH + SC4H9O2 <=> CH2OH + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3OH + SC4H9O2 <=> CH2OH + SC4H9O2H""",
)

entry(
    index = 859,
    label = "C2H3CHO + PC4H9O2 <=> C2H3CO + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + PC4H9O2 <=> C2H3CO + PC4H9O2H""",
)

entry(
    index = 860,
    label = "C2H3CHO + SC4H9O2 <=> C2H3CO + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHO + SC4H9O2 <=> C2H3CO + SC4H9O2H""",
)

entry(
    index = 861,
    label = "CH4 + PC4H9O2 <=> CH3 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(24640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + PC4H9O2 <=> CH3 + PC4H9O2H""",
)

entry(
    index = 862,
    label = "CH4 + SC4H9O2 <=> CH3 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(24640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH4 + SC4H9O2 <=> CH3 + SC4H9O2H""",
)

entry(
    index = 863,
    label = "C4H71-3 + PC4H9O2 <=> C4H7O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + PC4H9O2 <=> C4H7O + PC4H9O""",
)

entry(
    index = 864,
    label = "C4H71-3 + SC4H9O2 <=> C4H7O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + SC4H9O2 <=> C4H7O + SC4H9O""",
)

entry(
    index = 865,
    label = "H2O2 + PC4H9O2 <=> HO2 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + PC4H9O2 <=> HO2 + PC4H9O2H""",
)

entry(
    index = 866,
    label = "H2O2 + SC4H9O2 <=> HO2 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + SC4H9O2 <=> HO2 + SC4H9O2H""",
)

entry(
    index = 867,
    label = "PC4H9O2 + PC4H9O2 => O2 + PC4H9O + PC4H9O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + PC4H9O2 => O2 + PC4H9O + PC4H9O""",
)

entry(
    index = 868,
    label = "PC4H9O2 + SC4H9O2 => PC4H9O + SC4H9O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + SC4H9O2 => PC4H9O + SC4H9O + O2""",
)

entry(
    index = 869,
    label = "PC4H9O2 + NC3H7O2 => PC4H9O + NC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + NC3H7O2 => PC4H9O + NC3H7O + O2""",
)

entry(
    index = 870,
    label = "PC4H9O2 + IC3H7O2 => PC4H9O + IC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + IC3H7O2 => PC4H9O + IC3H7O + O2""",
)

entry(
    index = 871,
    label = "PC4H9O2 + C2H5O2 => PC4H9O + C2H5O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C2H5O2 => PC4H9O + C2H5O + O2""",
)

entry(
    index = 872,
    label = "PC4H9O2 + CH3O2 => PC4H9O + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + CH3O2 => PC4H9O + CH3O + O2""",
)

entry(
    index = 873,
    label = "PC4H9O2 + CH3CO3 => PC4H9O + CH3CO2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + CH3CO3 => PC4H9O + CH3CO2 + O2""",
)

entry(
    index = 874,
    label = "PC4H9O2 + CH3 <=> PC4H9O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + CH3 <=> PC4H9O + CH3O""",
)

entry(
    index = 875,
    label = "PC4H9O2 + C2H5 <=> PC4H9O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C2H5 <=> PC4H9O + C2H5O""",
)

entry(
    index = 876,
    label = "PC4H9O2 + IC3H7 <=> PC4H9O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + IC3H7 <=> PC4H9O + IC3H7O""",
)

entry(
    index = 877,
    label = "PC4H9O2 + NC3H7 <=> PC4H9O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + NC3H7 <=> PC4H9O + NC3H7O""",
)

entry(
    index = 878,
    label = "PC4H9O2 + PC4H9 <=> PC4H9O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + PC4H9 <=> PC4H9O + PC4H9O""",
)

entry(
    index = 879,
    label = "PC4H9O2 + SC4H9 <=> PC4H9O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + SC4H9 <=> PC4H9O + SC4H9O""",
)

entry(
    index = 880,
    label = "PC4H9O2 + C3H5-A <=> PC4H9O + C3H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + C3H5-A <=> PC4H9O + C3H5O""",
)

entry(
    index = 881,
    label = "PC4H9 + HO2 <=> PC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9 + HO2 <=> PC4H9O + OH""",
)

entry(
    index = 882,
    label = "SC4H9 + HO2 <=> SC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9 + HO2 <=> SC4H9O + OH""",
)

entry(
    index = 883,
    label = "CH3O2 + PC4H9 <=> CH3O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + PC4H9 <=> CH3O + PC4H9O""",
)

entry(
    index = 884,
    label = "CH3O2 + SC4H9 <=> CH3O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + SC4H9 <=> CH3O + SC4H9O""",
)

entry(
    index = 885,
    label = "PC4H9O2H <=> PC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2H <=> PC4H9O + OH""",
)

entry(
    index = 886,
    label = "SC4H9O2H <=> SC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.45e+15, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2H <=> SC4H9O + OH""",
)

entry(
    index = 887,
    label = "PC4H9O <=> NC3H7 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.558e+21, 's^-1'), n=-2.444, Ea=(15230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O <=> NC3H7 + CH2O""",
)

entry(
    index = 888,
    label = "SC4H9O <=> CH3 + C2H5CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.381e+16, 's^-1'), n=-0.893, Ea=(15200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O <=> CH3 + C2H5CHO""",
)

entry(
    index = 889,
    label = "SC4H9O <=> C2H5 + CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.489e+22, 's^-1'), n=-2.757, Ea=(12650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O <=> C2H5 + CH3CHO""",
)

entry(
    index = 890,
    label = "PC4H9O2 <=> C4H8OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 <=> C4H8OOH1-2""",
)

entry(
    index = 891,
    label = "PC4H9O2 <=> C4H8OOH1-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 <=> C4H8OOH1-3""",
)

entry(
    index = 892,
    label = "PC4H9O2 <=> C4H8OOH1-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.688e+09, 's^-1'), n=0, Ea=(22350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 <=> C4H8OOH1-4""",
)

entry(
    index = 893,
    label = "SC4H9O2 <=> C4H8OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 <=> C4H8OOH2-1""",
)

entry(
    index = 894,
    label = "SC4H9O2 <=> C4H8OOH2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 <=> C4H8OOH2-3""",
)

entry(
    index = 895,
    label = "SC4H9O2 <=> C4H8OOH2-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+10, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 <=> C4H8OOH2-4""",
)

entry(
    index = 896,
    label = "PC4H9O2 <=> C4H8-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 <=> C4H8-1 + HO2""",
)

entry(
    index = 897,
    label = "SC4H9O2 <=> C4H8-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.075e+42, 's^-1'), n=-9.41, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 <=> C4H8-1 + HO2""",
)

entry(
    index = 898,
    label = "SC4H9O2 <=> C4H8-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 <=> C4H8-2 + HO2""",
)

entry(
    index = 899,
    label = "C4H8OOH1-2 <=> C4H8-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.828e+16, 's^-1'), n=-1.488, Ea=(16260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-2 <=> C4H8-1 + HO2""",
)

entry(
    index = 900,
    label = "C4H8OOH2-1 <=> C4H8-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.851e+20, 's^-1'), n=-2.574, Ea=(21180, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-1 <=> C4H8-1 + HO2""",
)

entry(
    index = 901,
    label = "C4H8OOH2-3 <=> C4H8-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.224e+19, 's^-1'), n=-2.513, Ea=(21020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-3 <=> C4H8-2 + HO2""",
)

entry(
    index = 902,
    label = "C4H8OOH1-2 => C4H8O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-2 => C4H8O1-2 + OH""",
)

entry(
    index = 903,
    label = "C4H8OOH1-3 => C4H8O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-3 => C4H8O1-3 + OH""",
)

entry(
    index = 904,
    label = "C4H8OOH1-4 => C4H8O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-4 => C4H8O1-4 + OH""",
)

entry(
    index = 905,
    label = "C4H8OOH2-1 => C4H8O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-1 => C4H8O1-2 + OH""",
)

entry(
    index = 906,
    label = "C4H8OOH2-3 => C4H8O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-3 => C4H8O2-3 + OH""",
)

entry(
    index = 907,
    label = "C4H8OOH2-4 => C4H8O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-4 => C4H8O1-3 + OH""",
)

entry(
    index = 908,
    label = "C4H8OOH1-1 <=> NC3H7CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+14, 's^-1'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-1 <=> NC3H7CHO + OH""",
)

entry(
    index = 909,
    label = "C4H8OOH2-2 <=> C2H5COCH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+14, 's^-1'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-2 <=> C2H5COCH3 + OH""",
)

entry(
    index = 910,
    label = "C4H8OOH1-3 => OH + CH2O + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6.635e+13, 's^-1'), n=-0.16, Ea=(29900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-3 => OH + CH2O + C3H6""",
)

entry(
    index = 911,
    label = "C4H8OOH2-4 => OH + CH3CHO + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.945e+18, 's^-1'), n=-1.63, Ea=(26790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-4 => OH + CH3CHO + C2H4""",
)

entry(
    index = 912,
    label = "C4H8OOH1-2O2 <=> C4H8OOH1-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.586e+24, 's^-1'), n=-2.709, Ea=(39860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-2O2 <=> C4H8OOH1-2 + O2""",
)

entry(
    index = 913,
    label = "C4H8OOH1-3O2 <=> C4H8OOH1-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.601e+22, 's^-1'), n=-2.234, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-3O2 <=> C4H8OOH1-3 + O2""",
)

entry(
    index = 914,
    label = "C4H8OOH1-4O2 <=> C4H8OOH1-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.569e+20, 's^-1'), n=-1.611, Ea=(35680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-4O2 <=> C4H8OOH1-4 + O2""",
)

entry(
    index = 915,
    label = "C4H8OOH2-1O2 <=> C4H8OOH2-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.821e+20, 's^-1'), n=-1.622, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-1O2 <=> C4H8OOH2-1 + O2""",
)

entry(
    index = 916,
    label = "C4H8OOH2-3O2 <=> C4H8OOH2-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.516e+22, 's^-1'), n=-2.218, Ea=(37880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-3O2 <=> C4H8OOH2-3 + O2""",
)

entry(
    index = 917,
    label = "C4H8OOH2-4O2 <=> C4H8OOH2-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.821e+20, 's^-1'), n=-1.622, Ea=(35700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-4O2 <=> C4H8OOH2-4 + O2""",
)

entry(
    index = 918,
    label = "C4H8OOH1-2O2 <=> NC4KET12 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-2O2 <=> NC4KET12 + OH""",
)

entry(
    index = 919,
    label = "C4H8OOH1-3O2 <=> NC4KET13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(21400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-3O2 <=> NC4KET13 + OH""",
)

entry(
    index = 920,
    label = "C4H8OOH1-4O2 <=> NC4KET14 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(19350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH1-4O2 <=> NC4KET14 + OH""",
)

entry(
    index = 921,
    label = "C4H8OOH2-1O2 <=> NC4KET21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-1O2 <=> NC4KET21 + OH""",
)

entry(
    index = 922,
    label = "C4H8OOH2-3O2 <=> NC4KET23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-3O2 <=> NC4KET23 + OH""",
)

entry(
    index = 923,
    label = "C4H8OOH2-4O2 <=> NC4KET24 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8OOH2-4O2 <=> NC4KET24 + OH""",
)

entry(
    index = 924,
    label = "NC4KET12 => C2H5CHO + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4KET12 => C2H5CHO + HCO + OH""",
)

entry(
    index = 925,
    label = "NC4KET13 => CH3CHO + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4KET13 => CH3CHO + CH2CHO + OH""",
)

entry(
    index = 926,
    label = "NC4KET14 => CH2CH2CHO + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4KET14 => CH2CH2CHO + CH2O + OH""",
)

entry(
    index = 927,
    label = "NC4KET21 => CH2O + C2H5CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4KET21 => CH2O + C2H5CO + OH""",
)

entry(
    index = 928,
    label = "NC4KET23 => CH3CHO + CH3CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4KET23 => CH3CHO + CH3CO + OH""",
)

entry(
    index = 929,
    label = "NC4KET24 => CH2O + CH3COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4KET24 => CH2O + CH3COCH2 + OH""",
)

entry(
    index = 930,
    label = "C2H5COCH3 + OH <=> CH2CH2COCH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.55e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + OH <=> CH2CH2COCH3 + H2O""",
)

entry(
    index = 931,
    label = "C2H5COCH3 + OH <=> CH3CHCOCH3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + OH <=> CH3CHCOCH3 + H2O""",
)

entry(
    index = 932,
    label = "C2H5COCH3 + OH <=> C2H5COCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + OH <=> C2H5COCH2 + H2O""",
)

entry(
    index = 933,
    label = "C2H5COCH3 + HO2 <=> CH2CH2COCH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + HO2 <=> CH2CH2COCH3 + H2O2""",
)

entry(
    index = 934,
    label = "C2H5COCH3 + HO2 <=> CH3CHCOCH3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + HO2 <=> CH3CHCOCH3 + H2O2""",
)

entry(
    index = 935,
    label = "C2H5COCH3 + HO2 <=> C2H5COCH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(14690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + HO2 <=> C2H5COCH2 + H2O2""",
)

entry(
    index = 936,
    label = "C2H5COCH3 + O <=> CH2CH2COCH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.25e+13, 'cm^3/(mol*s)'), n=0, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + O <=> CH2CH2COCH3 + OH""",
)

entry(
    index = 937,
    label = "C2H5COCH3 + O <=> CH3CHCOCH3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.07e+13, 'cm^3/(mol*s)'), n=0, Ea=(3400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + O <=> CH3CHCOCH3 + OH""",
)

entry(
    index = 938,
    label = "C2H5COCH3 + O <=> C2H5COCH2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(5962, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + O <=> C2H5COCH2 + OH""",
)

entry(
    index = 939,
    label = "C2H5COCH3 + H <=> CH2CH2COCH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.16e+06, 'cm^3/(mol*s)'), n=2, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + H <=> CH2CH2COCH3 + H2""",
)

entry(
    index = 940,
    label = "C2H5COCH3 + H <=> CH3CHCOCH3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.46e+06, 'cm^3/(mol*s)'), n=2, Ea=(3200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + H <=> CH3CHCOCH3 + H2""",
)

entry(
    index = 941,
    label = "C2H5COCH3 + H <=> C2H5COCH2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(6357, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + H <=> C2H5COCH2 + H2""",
)

entry(
    index = 942,
    label = "C2H5COCH3 + O2 <=> CH2CH2COCH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.05e+13, 'cm^3/(mol*s)'), n=0, Ea=(51310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + O2 <=> CH2CH2COCH3 + HO2""",
)

entry(
    index = 943,
    label = "C2H5COCH3 + O2 <=> CH3CHCOCH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(41970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + O2 <=> CH3CHCOCH3 + HO2""",
)

entry(
    index = 944,
    label = "C2H5COCH3 + O2 <=> C2H5COCH2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.05e+13, 'cm^3/(mol*s)'), n=0, Ea=(49150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + O2 <=> C2H5COCH2 + HO2""",
)

entry(
    index = 945,
    label = "C2H5COCH3 + CH3 <=> CH2CH2COCH3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(31.9, 'cm^3/(mol*s)'), n=3.17, Ea=(7172, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3 <=> CH2CH2COCH3 + CH4""",
)

entry(
    index = 946,
    label = "C2H5COCH3 + CH3 <=> CH3CHCOCH3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.74, 'cm^3/(mol*s)'), n=3.46, Ea=(3680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3 <=> CH3CHCOCH3 + CH4""",
)

entry(
    index = 947,
    label = "C2H5COCH3 + CH3 <=> C2H5COCH2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.62e+11, 'cm^3/(mol*s)'), n=0, Ea=(9630, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3 <=> C2H5COCH2 + CH4""",
)

entry(
    index = 948,
    label = "C2H5COCH3 + CH3O <=> CH2CH2COCH3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3O <=> CH2CH2COCH3 + CH3OH""",
)

entry(
    index = 949,
    label = "C2H5COCH3 + CH3O <=> CH3CHCOCH3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(2771, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3O <=> CH3CHCOCH3 + CH3OH""",
)

entry(
    index = 950,
    label = "C2H5COCH3 + CH3O <=> C2H5COCH2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(4660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3O <=> C2H5COCH2 + CH3OH""",
)

entry(
    index = 951,
    label = "C2H5COCH3 + CH3O2 <=> CH2CH2COCH3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3O2 <=> CH2CH2COCH3 + CH3O2H""",
)

entry(
    index = 952,
    label = "C2H5COCH3 + CH3O2 <=> CH3CHCOCH3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3O2 <=> CH3CHCOCH3 + CH3O2H""",
)

entry(
    index = 953,
    label = "C2H5COCH3 + CH3O2 <=> C2H5COCH2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(17580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + CH3O2 <=> C2H5COCH2 + CH3O2H""",
)

entry(
    index = 954,
    label = "C2H5COCH3 + C2H3 <=> CH2CH2COCH3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + C2H3 <=> CH2CH2COCH3 + C2H4""",
)

entry(
    index = 955,
    label = "C2H5COCH3 + C2H3 <=> CH3CHCOCH3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(3400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + C2H3 <=> CH3CHCOCH3 + C2H4""",
)

entry(
    index = 956,
    label = "C2H5COCH3 + C2H3 <=> C2H5COCH2 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.15e+10, 'cm^3/(mol*s)'), n=0, Ea=(4278, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + C2H3 <=> C2H5COCH2 + C2H4""",
)

entry(
    index = 957,
    label = "C2H5COCH3 + C2H5 <=> CH2CH2COCH3 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(13400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + C2H5 <=> CH2CH2COCH3 + C2H6""",
)

entry(
    index = 958,
    label = "C2H5COCH3 + C2H5 <=> CH3CHCOCH3 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+10, 'cm^3/(mol*s)'), n=0, Ea=(8600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + C2H5 <=> CH3CHCOCH3 + C2H6""",
)

entry(
    index = 959,
    label = "C2H5COCH3 + C2H5 <=> C2H5COCH2 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(11600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH3 + C2H5 <=> C2H5COCH2 + C2H6""",
)

entry(
    index = 960,
    label = "CH3CHOOCOCH3 <=> CH3CHCOCH3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.372e+17, 's^-1'), n=-1.69, Ea=(28460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOOCOCH3 <=> CH3CHCOCH3 + O2""",
)

entry(
    index = 961,
    label = "CH3CHOOCOCH3 <=> CH2CHOOHCOCH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.9e+12, 's^-1'), n=0, Ea=(29700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHOOCOCH3 <=> CH2CHOOHCOCH3""",
)

entry(
    index = 962,
    label = "CH2CHOOHCOCH3 <=> C2H3COCH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.026e+19, 's^-1'), n=-2.35, Ea=(14130, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CHOOHCOCH3 <=> C2H3COCH3 + HO2""",
)

entry(
    index = 963,
    label = "CH2CH2CHO <=> C2H4 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.127e+13, 's^-1'), n=-0.52, Ea=(24590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2CHO <=> C2H4 + HCO""",
)

entry(
    index = 964,
    label = "CH2CH2COCH3 <=> C2H4 + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(18000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CH2COCH3 <=> C2H4 + CH3CO""",
)

entry(
    index = 965,
    label = "C2H5COCH2 <=> CH2CO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 's^-1'), n=0, Ea=(35000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COCH2 <=> CH2CO + C2H5""",
)

entry(
    index = 966,
    label = "CH3CHCOCH3 <=> C2H3COCH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.417e+16, 's^-1'), n=-0.82, Ea=(41770, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCOCH3 <=> C2H3COCH3 + H""",
)

entry(
    index = 967,
    label = "CH3CHCOCH3 <=> CH3CHCO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.406e+15, 's^-1'), n=-0.44, Ea=(38340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3CHCOCH3 <=> CH3CHCO + CH3""",
)

entry(
    index = 968,
    label = "NC3H7CHO + O2 <=> NC3H7CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(120000, 'cm^3/(mol*s)'), n=2.5, Ea=(37560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + O2 <=> NC3H7CO + HO2""",
)

entry(
    index = 969,
    label = "NC3H7CHO + OH <=> NC3H7CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+06, 'cm^3/(mol*s)'), n=1.8, Ea=(-1300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + OH <=> NC3H7CO + H2O""",
)

entry(
    index = 970,
    label = "NC3H7CHO + H <=> NC3H7CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.14e+09, 'cm^3/(mol*s)'),
        n = 1.12,
        Ea = (2320, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + H <=> NC3H7CO + H2""",
)

entry(
    index = 971,
    label = "NC3H7CHO + O <=> NC3H7CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.94e+12, 'cm^3/(mol*s)'), n=0, Ea=(1868, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + O <=> NC3H7CO + OH""",
)

entry(
    index = 972,
    label = "NC3H7CHO + HO2 <=> NC3H7CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40900, 'cm^3/(mol*s)'), n=2.5, Ea=(10200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + HO2 <=> NC3H7CO + H2O2""",
)

entry(
    index = 973,
    label = "NC3H7CHO + CH3 <=> NC3H7CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (0.00289, 'cm^3/(mol*s)'),
        n = 4.62,
        Ea = (3210, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + CH3 <=> NC3H7CO + CH4""",
)

entry(
    index = 974,
    label = "NC3H7CHO + CH3O <=> NC3H7CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(3300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + CH3O <=> NC3H7CO + CH3OH""",
)

entry(
    index = 975,
    label = "NC3H7CHO + CH3O2 <=> NC3H7CO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40900, 'cm^3/(mol*s)'), n=2.5, Ea=(10200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + CH3O2 <=> NC3H7CO + CH3O2H""",
)

entry(
    index = 976,
    label = "NC3H7CHO + OH <=> C3H6CHO-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.28e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + OH <=> C3H6CHO-1 + H2O""",
)

entry(
    index = 977,
    label = "NC3H7CHO + OH <=> C3H6CHO-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.68e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + OH <=> C3H6CHO-2 + H2O""",
)

entry(
    index = 978,
    label = "NC3H7CHO + OH <=> C3H6CHO-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(552, 'cm^3/(mol*s)'), n=3.12, Ea=(-1176, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + OH <=> C3H6CHO-3 + H2O""",
)

entry(
    index = 979,
    label = "NC3H7CHO + HO2 <=> C3H6CHO-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23790, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + HO2 <=> C3H6CHO-1 + H2O2""",
)

entry(
    index = 980,
    label = "NC3H7CHO + HO2 <=> C3H6CHO-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + HO2 <=> C3H6CHO-2 + H2O2""",
)

entry(
    index = 981,
    label = "NC3H7CHO + HO2 <=> C3H6CHO-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.44e+12, 'cm^3/(mol*s)'),
        n = 0.05,
        Ea = (17880, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + HO2 <=> C3H6CHO-3 + H2O2""",
)

entry(
    index = 982,
    label = "NC3H7CHO + CH3O2 <=> C3H6CHO-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23790, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + CH3O2 <=> C3H6CHO-1 + CH3O2H""",
)

entry(
    index = 983,
    label = "NC3H7CHO + CH3O2 <=> C3H6CHO-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + CH3O2 <=> C3H6CHO-2 + CH3O2H""",
)

entry(
    index = 984,
    label = "NC3H7CHO + CH3O2 <=> C3H6CHO-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.44e+12, 'cm^3/(mol*s)'),
        n = 0.05,
        Ea = (17880, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7CHO + CH3O2 <=> C3H6CHO-3 + CH3O2H""",
)

entry(
    index = 985,
    label = "NC3H7CO <=> NC3H7 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(9600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7CO <=> NC3H7 + CO""",
)

entry(
    index = 986,
    label = "C3H6CHO-1 <=> C2H4 + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.4e+11, 's^-1'), n=0, Ea=(21970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6CHO-1 <=> C2H4 + CH2CHO""",
)

entry(
    index = 987,
    label = "C3H6CHO-3 <=> C2H5CHCO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.431e+15, 's^-1'), n=-0.6, Ea=(40400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6CHO-3 <=> C2H5CHCO + H""",
)

entry(
    index = 988,
    label = "C3H6CHO-3 <=> C2H3CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.174e+14, 's^-1'), n=-0.39, Ea=(29900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6CHO-3 <=> C2H3CHO + CH3""",
)

entry(
    index = 989,
    label = "C3H6CHO-2 <=> SC3H5CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.947e+12, 's^-1'), n=-0.15, Ea=(31300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6CHO-2 <=> SC3H5CHO + H""",
)

entry(
    index = 990,
    label = "C3H6CHO-2 <=> C3H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.249e+12, 's^-1'), n=-0.18, Ea=(21900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6CHO-2 <=> C3H6 + HCO""",
)

entry(
    index = 991,
    label = "C2H5CHCO + OH => NC3H7 + CO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.73e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHCO + OH => NC3H7 + CO2""",
)

entry(
    index = 992,
    label = "C2H5CHCO + H => NC3H7 + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(1459, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHCO + H => NC3H7 + CO""",
)

entry(
    index = 993,
    label = "C2H5CHCO + O => C3H6 + CO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(-437, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5CHCO + O => C3H6 + CO2""",
)

entry(
    index = 994,
    label = "SC3H5CHO + OH <=> SC3H5CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is SC3H5CHO + OH <=> SC3H5CO + H2O""",
)

entry(
    index = 995,
    label = "SC3H5CO <=> C3H5-S + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.6e+15, 's^-1'), n=0, Ea=(23000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC3H5CO <=> C3H5-S + CO""",
)

entry(
    index = 996,
    label = "SC3H5CHO + HO2 <=> SC3H5CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC3H5CHO + HO2 <=> SC3H5CO + H2O2""",
)

entry(
    index = 997,
    label = "SC3H5CHO + CH3 <=> SC3H5CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+12, 'cm^3/(mol*s)'), n=0, Ea=(8700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC3H5CHO + CH3 <=> SC3H5CO + CH4""",
)

entry(
    index = 998,
    label = "SC3H5CHO + O <=> SC3H5CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.18e+12, 'cm^3/(mol*s)'), n=0, Ea=(1389, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC3H5CHO + O <=> SC3H5CO + OH""",
)

entry(
    index = 999,
    label = "SC3H5CHO + O2 <=> SC3H5CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC3H5CHO + O2 <=> SC3H5CO + HO2""",
)

entry(
    index = 1000,
    label = "SC3H5CHO + H <=> SC3H5CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC3H5CHO + H <=> SC3H5CO + H2""",
)

entry(
    index = 1001,
    label = "C2H3COCH3 + OH => CH3CHO + CH3CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3COCH3 + OH => CH3CHO + CH3CO""",
)

entry(
    index = 1002,
    label = "C2H3COCH3 + OH => CH2CO + C2H3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3COCH3 + OH => CH2CO + C2H3 + H2O""",
)

entry(
    index = 1003,
    label = "C2H3COCH3 + HO2 => CH2CHO + CH3CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6.03e+09, 'cm^3/(mol*s)'), n=0, Ea=(7949, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3COCH3 + HO2 => CH2CHO + CH3CO + OH""",
)

entry(
    index = 1004,
    label = "C2H3COCH3 + HO2 => CH2CO + C2H3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3COCH3 + HO2 => CH2CO + C2H3 + H2O2""",
)

entry(
    index = 1005,
    label = "C2H3COCH3 + CH3O2 => CH2CHO + CH3CO + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.97e+11, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3COCH3 + CH3O2 => CH2CHO + CH3CO + CH3O""",
)

entry(
    index = 1006,
    label = "C2H3COCH3 + CH3O2 => CH2CO + C2H3 + CH3O2H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(17580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3COCH3 + CH3O2 => CH2CO + C2H3 + CH3O2H""",
)

entry(
    index = 1007,
    label = "IC4H10 <=> CH3 + IC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.83e+16, 's^-1'), n=0, Ea=(79900, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(A=(2.41e+19, 'cm^3/(mol*s)'), n=0, Ea=(52576, 'cal/mol'), T0=(1, 'K')),
        alpha = 0.25,
        T3 = (750, 'K'),
        T1 = (1e-10, 'K'),
        T2 = (1e+10, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is IC4H10 <=> CH3 + IC3H7""",
)

entry(
    index = 1008,
    label = "IC4H10 <=> TC4H9 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.51e+98, 's^-1'), n=-23.81, Ea=(145300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 <=> TC4H9 + H""",
)

entry(
    index = 1009,
    label = "IC4H10 <=> IC4H9 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.85e+95, 's^-1'), n=-23.11, Ea=(147600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 <=> IC4H9 + H""",
)

entry(
    index = 1010,
    label = "IC4H10 + H <=> TC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(602000, 'cm^3/(mol*s)'), n=2.4, Ea=(2583, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + H <=> TC4H9 + H2""",
)

entry(
    index = 1011,
    label = "IC4H10 + H <=> IC4H9 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.81e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6756, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H10 + H <=> IC4H9 + H2""",
)

entry(
    index = 1012,
    label = "IC4H10 + CH3 <=> TC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.904, 'cm^3/(mol*s)'), n=3.46, Ea=(4598, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3 <=> TC4H9 + CH4""",
)

entry(
    index = 1013,
    label = "IC4H10 + CH3 <=> IC4H9 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.36, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3 <=> IC4H9 + CH4""",
)

entry(
    index = 1014,
    label = "IC4H10 + OH <=> TC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (29250, 'cm^3/(mol*s)'),
        n = 2.531,
        Ea = (-1659, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H10 + OH <=> TC4H9 + H2O""",
)

entry(
    index = 1015,
    label = "IC4H10 + OH <=> IC4H9 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (66540, 'cm^3/(mol*s)'),
        n = 2.665,
        Ea = (-168.9, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H10 + OH <=> IC4H9 + H2O""",
)

entry(
    index = 1016,
    label = "IC4H10 + C2H5 <=> IC4H9 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51e+12, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + C2H5 <=> IC4H9 + C2H6""",
)

entry(
    index = 1017,
    label = "IC4H10 + C2H5 <=> TC4H9 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(7900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + C2H5 <=> TC4H9 + C2H6""",
)

entry(
    index = 1018,
    label = "IC4H10 + HO2 <=> IC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(61.2, 'cm^3/(mol*s)'), n=3.59, Ea=(17160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + HO2 <=> IC4H9 + H2O2""",
)

entry(
    index = 1019,
    label = "IC4H10 + HO2 <=> TC4H9 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(433.2, 'cm^3/(mol*s)'), n=3.01, Ea=(12090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + HO2 <=> TC4H9 + H2O2""",
)

entry(
    index = 1020,
    label = "IC4H10 + O <=> TC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (196800, 'cm^3/(mol*s)'),
        n = 2.402,
        Ea = (1150, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H10 + O <=> TC4H9 + OH""",
)

entry(
    index = 1021,
    label = "IC4H10 + O <=> IC4H9 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.046e+07, 'cm^3/(mol*s)'),
        n = 2.034,
        Ea = (5136, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H10 + O <=> IC4H9 + OH""",
)

entry(
    index = 1022,
    label = "IC4H10 + CH3O <=> IC4H9 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3O <=> IC4H9 + CH3OH""",
)

entry(
    index = 1023,
    label = "IC4H10 + CH3O <=> TC4H9 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.9e+10, 'cm^3/(mol*s)'), n=0, Ea=(2800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3O <=> TC4H9 + CH3OH""",
)

entry(
    index = 1024,
    label = "IC4H10 + O2 <=> IC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + O2 <=> IC4H9 + HO2""",
)

entry(
    index = 1025,
    label = "IC4H10 + O2 <=> TC4H9 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(48200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + O2 <=> TC4H9 + HO2""",
)

entry(
    index = 1026,
    label = "IC4H10 + CH3O2 <=> IC4H9 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.079, 'cm^3/(mol*s)'), n=3.97, Ea=(18280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3O2 <=> IC4H9 + CH3O2H""",
)

entry(
    index = 1027,
    label = "IC4H10 + C2H5O2 <=> IC4H9 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + C2H5O2 <=> IC4H9 + C2H5O2H""",
)

entry(
    index = 1028,
    label = "IC4H10 + CH3CO3 <=> IC4H9 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3CO3 <=> IC4H9 + CH3CO3H""",
)

entry(
    index = 1029,
    label = "IC4H10 + NC3H7O2 <=> IC4H9 + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + NC3H7O2 <=> IC4H9 + NC3H7O2H""",
)

entry(
    index = 1030,
    label = "IC4H10 + IC3H7O2 <=> IC4H9 + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + IC3H7O2 <=> IC4H9 + IC3H7O2H""",
)

entry(
    index = 1031,
    label = "IC4H10 + IC4H9O2 <=> IC4H9 + IC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + IC4H9O2 <=> IC4H9 + IC4H9O2H""",
)

entry(
    index = 1032,
    label = "IC4H10 + TC4H9O2 <=> IC4H9 + TC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.55e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + TC4H9O2 <=> IC4H9 + TC4H9O2H""",
)

entry(
    index = 1033,
    label = "IC4H10 + O2CHO <=> IC4H9 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.52e+13, 'cm^3/(mol*s)'), n=0, Ea=(20440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + O2CHO <=> IC4H9 + HO2CHO""",
)

entry(
    index = 1034,
    label = "IC4H10 + O2CHO <=> TC4H9 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + O2CHO <=> TC4H9 + HO2CHO""",
)

entry(
    index = 1035,
    label = "IC4H10 + SC4H9O2 <=> IC4H9 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.25e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + SC4H9O2 <=> IC4H9 + SC4H9O2H""",
)

entry(
    index = 1036,
    label = "IC4H10 + SC4H9O2 <=> TC4H9 + SC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + SC4H9O2 <=> TC4H9 + SC4H9O2H""",
)

entry(
    index = 1037,
    label = "IC4H10 + PC4H9O2 <=> IC4H9 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.25e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + PC4H9O2 <=> IC4H9 + PC4H9O2H""",
)

entry(
    index = 1038,
    label = "IC4H10 + PC4H9O2 <=> TC4H9 + PC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + PC4H9O2 <=> TC4H9 + PC4H9O2H""",
)

entry(
    index = 1039,
    label = "IC4H10 + CH3O2 <=> TC4H9 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(136.6, 'cm^3/(mol*s)'), n=3.12, Ea=(13190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3O2 <=> TC4H9 + CH3O2H""",
)

entry(
    index = 1040,
    label = "IC4H10 + C2H5O2 <=> TC4H9 + C2H5O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + C2H5O2 <=> TC4H9 + C2H5O2H""",
)

entry(
    index = 1041,
    label = "IC4H10 + CH3CO3 <=> TC4H9 + CH3CO3H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + CH3CO3 <=> TC4H9 + CH3CO3H""",
)

entry(
    index = 1042,
    label = "IC4H10 + NC3H7O2 <=> TC4H9 + NC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + NC3H7O2 <=> TC4H9 + NC3H7O2H""",
)

entry(
    index = 1043,
    label = "IC4H10 + IC3H7O2 <=> TC4H9 + IC3H7O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + IC3H7O2 <=> TC4H9 + IC3H7O2H""",
)

entry(
    index = 1044,
    label = "IC4H10 + IC4H9O2 <=> TC4H9 + IC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + IC4H9O2 <=> TC4H9 + IC4H9O2H""",
)

entry(
    index = 1045,
    label = "IC4H10 + TC4H9O2 <=> TC4H9 + TC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(16000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + TC4H9O2 <=> TC4H9 + TC4H9O2H""",
)

entry(
    index = 1046,
    label = "IC4H10 + IC4H9 <=> TC4H9 + IC4H10",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(7900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H10 + IC4H9 <=> TC4H9 + IC4H10""",
)

entry(
    index = 1047,
    label = "IC4H9 + HO2 <=> IC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9 + HO2 <=> IC4H9O + OH""",
)

entry(
    index = 1048,
    label = "TC4H9 + HO2 <=> TC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9 + HO2 <=> TC4H9O + OH""",
)

entry(
    index = 1049,
    label = "CH3O2 + IC4H9 <=> CH3O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + IC4H9 <=> CH3O + IC4H9O""",
)

entry(
    index = 1050,
    label = "CH3O2 + TC4H9 <=> CH3O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + TC4H9 <=> CH3O + TC4H9O""",
)

entry(
    index = 1051,
    label = "IC4H9 <=> IC4H8 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.371e+13, 's^-1'), n=0.124, Ea=(33660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9 <=> IC4H8 + H""",
)

entry(
    index = 1052,
    label = "IC4H9 <=> C3H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.504e+11, 's^-1'), n=0.773, Ea=(30700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9 <=> C3H6 + CH3""",
)

entry(
    index = 1053,
    label = "TC4H9 <=> H + IC4H8",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.128e+12, 's^-1'), n=0.703, Ea=(36560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9 <=> H + IC4H8""",
)

entry(
    index = 1054,
    label = "TC4H9 + O2 <=> IC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.837, 'cm^3/(mol*s)'), n=3.59, Ea=(11960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9 + O2 <=> IC4H8 + HO2""",
)

entry(
    index = 1055,
    label = "IC4H9 + O2 <=> IC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07, 'cm^3/(mol*s)'), n=3.71, Ea=(9322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9 + O2 <=> IC4H8 + HO2""",
)

entry(
    index = 1056,
    label = "NC3H7O2 + IC4H9 <=> NC3H7O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + IC4H9 <=> NC3H7O + IC4H9O""",
)

entry(
    index = 1057,
    label = "NC3H7O2 + TC4H9 <=> NC3H7O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + TC4H9 <=> NC3H7O + TC4H9O""",
)

entry(
    index = 1058,
    label = "NC3H7O2 + IC4H7 <=> NC3H7O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + IC4H7 <=> NC3H7O + IC4H7O""",
)

entry(
    index = 1059,
    label = "SC4H9O2 + IC4H9 <=> SC4H9O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + IC4H9 <=> SC4H9O + IC4H9O""",
)

entry(
    index = 1060,
    label = "SC4H9O2 + TC4H9 <=> SC4H9O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + TC4H9 <=> SC4H9O + TC4H9O""",
)

entry(
    index = 1061,
    label = "PC4H9O2 + IC4H9 <=> PC4H9O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + IC4H9 <=> PC4H9O + IC4H9O""",
)

entry(
    index = 1062,
    label = "PC4H9O2 + TC4H9 <=> PC4H9O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + TC4H9 <=> PC4H9O + TC4H9O""",
)

entry(
    index = 1063,
    label = "PC4H9O2 + IC4H7 <=> PC4H9O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + IC4H7 <=> PC4H9O + IC4H7O""",
)

entry(
    index = 1064,
    label = "SC4H9O2 + IC4H7 <=> SC4H9O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + IC4H7 <=> SC4H9O + IC4H7O""",
)

entry(
    index = 1065,
    label = "IC4H9O2 <=> IC4H9 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.64e+19, 's^-1'), n=-1.575, Ea=(36080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 <=> IC4H9 + O2""",
)

entry(
    index = 1066,
    label = "TC4H9O2 <=> TC4H9 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.331e+24, 's^-1'), n=-2.472, Ea=(37870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 <=> TC4H9 + O2""",
)

entry(
    index = 1067,
    label = "IC4H9O2 + C4H10 <=> IC4H9O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C4H10 <=> IC4H9O2H + SC4H9""",
)

entry(
    index = 1068,
    label = "TC4H9O2 + C4H10 <=> TC4H9O2H + SC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C4H10 <=> TC4H9O2H + SC4H9""",
)

entry(
    index = 1069,
    label = "IC4H9O2 + C4H10 <=> IC4H9O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C4H10 <=> IC4H9O2H + PC4H9""",
)

entry(
    index = 1070,
    label = "TC4H9O2 + C4H10 <=> TC4H9O2H + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C4H10 <=> TC4H9O2H + PC4H9""",
)

entry(
    index = 1071,
    label = "IC3H7O2 + IC4H9 <=> IC3H7O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + IC4H9 <=> IC3H7O + IC4H9O""",
)

entry(
    index = 1072,
    label = "IC3H7O2 + TC4H9 <=> IC3H7O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + TC4H9 <=> IC3H7O + TC4H9O""",
)

entry(
    index = 1073,
    label = "IC3H7O2 + IC4H7 <=> IC3H7O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + IC4H7 <=> IC3H7O + IC4H7O""",
)

entry(
    index = 1074,
    label = "IC4H9O2 + C3H6 <=> IC4H9O2H + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C3H6 <=> IC4H9O2H + C3H5-A""",
)

entry(
    index = 1075,
    label = "TC4H9O2 + C3H6 <=> TC4H9O2H + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.24e+11, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C3H6 <=> TC4H9O2H + C3H5-A""",
)

entry(
    index = 1076,
    label = "IC4H9O2 + IC4H8 <=> IC4H9O2H + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + IC4H8 <=> IC4H9O2H + IC4H7""",
)

entry(
    index = 1077,
    label = "TC4H9O2 + IC4H8 <=> TC4H9O2H + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + IC4H8 <=> TC4H9O2H + IC4H7""",
)

entry(
    index = 1078,
    label = "PC4H9O2 + IC4H8 <=> PC4H9O2H + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9O2 + IC4H8 <=> PC4H9O2H + IC4H7""",
)

entry(
    index = 1079,
    label = "SC4H9O2 + IC4H8 <=> SC4H9O2H + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC4H9O2 + IC4H8 <=> SC4H9O2H + IC4H7""",
)

entry(
    index = 1080,
    label = "IC3H7O2 + IC4H8 <=> IC3H7O2H + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7O2 + IC4H8 <=> IC3H7O2H + IC4H7""",
)

entry(
    index = 1081,
    label = "NC3H7O2 + IC4H8 <=> NC3H7O2H + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7O2 + IC4H8 <=> NC3H7O2H + IC4H7""",
)

entry(
    index = 1082,
    label = "IC4H9O2 + C4H8-1 <=> IC4H9O2H + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C4H8-1 <=> IC4H9O2H + C4H71-3""",
)

entry(
    index = 1083,
    label = "TC4H9O2 + C4H8-1 <=> TC4H9O2H + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C4H8-1 <=> TC4H9O2H + C4H71-3""",
)

entry(
    index = 1084,
    label = "IC4H9O2 + C4H8-2 <=> IC4H9O2H + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C4H8-2 <=> IC4H9O2H + C4H71-3""",
)

entry(
    index = 1085,
    label = "TC4H9O2 + C4H8-2 <=> TC4H9O2H + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(14900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C4H8-2 <=> TC4H9O2H + C4H71-3""",
)

entry(
    index = 1086,
    label = "C2H4 + TC4H9O2 <=> C2H3 + TC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+11, 'cm^3/(mol*s)'), n=0, Ea=(17110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H4 + TC4H9O2 <=> C2H3 + TC4H9O2H""",
)

entry(
    index = 1087,
    label = "TC4H9O2 + CH4 <=> TC4H9O2H + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH4 <=> TC4H9O2H + CH3""",
)

entry(
    index = 1088,
    label = "H2 + TC4H9O2 <=> H + TC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + TC4H9O2 <=> H + TC4H9O2H""",
)

entry(
    index = 1089,
    label = "TC4H9O2 + C2H6 <=> TC4H9O2H + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H6 <=> TC4H9O2H + C2H5""",
)

entry(
    index = 1090,
    label = "TC4H9O2 + C3H8 <=> TC4H9O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C3H8 <=> TC4H9O2H + IC3H7""",
)

entry(
    index = 1091,
    label = "TC4H9O2 + C3H8 <=> TC4H9O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C3H8 <=> TC4H9O2H + NC3H7""",
)

entry(
    index = 1092,
    label = "TC4H9O2 + CH3OH <=> TC4H9O2H + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH3OH <=> TC4H9O2H + CH2OH""",
)

entry(
    index = 1093,
    label = "TC4H9O2 + C2H5OH <=> TC4H9O2H + PC2H4OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H5OH <=> TC4H9O2H + PC2H4OH""",
)

entry(
    index = 1094,
    label = "TC4H9O2 + C2H5OH <=> TC4H9O2H + SC2H4OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H5OH <=> TC4H9O2H + SC2H4OH""",
)

entry(
    index = 1095,
    label = "IC4H9O2 + CH3CHO <=> IC4H9O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH3CHO <=> IC4H9O2H + CH3CO""",
)

entry(
    index = 1096,
    label = "TC4H9O2 + CH3CHO <=> TC4H9O2H + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH3CHO <=> TC4H9O2H + CH3CO""",
)

entry(
    index = 1097,
    label = "IC4H9O2 + C2H3CHO <=> IC4H9O2H + C2H3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H3CHO <=> IC4H9O2H + C2H3CO""",
)

entry(
    index = 1098,
    label = "TC4H9O2 + C2H3CHO <=> TC4H9O2H + C2H3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H3CHO <=> TC4H9O2H + C2H3CO""",
)

entry(
    index = 1099,
    label = "IC4H9O2 + C2H5CHO <=> IC4H9O2H + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H5CHO <=> IC4H9O2H + C2H5CO""",
)

entry(
    index = 1100,
    label = "TC4H9O2 + C2H5CHO <=> TC4H9O2H + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H5CHO <=> TC4H9O2H + C2H5CO""",
)

entry(
    index = 1101,
    label = "IC4H9O2 + HO2 <=> IC4H9O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + HO2 <=> IC4H9O2H + O2""",
)

entry(
    index = 1102,
    label = "TC4H9O2 + HO2 <=> TC4H9O2H + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + HO2 <=> TC4H9O2H + O2""",
)

entry(
    index = 1103,
    label = "IC4H9O2 + H2O2 <=> IC4H9O2H + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + H2O2 <=> IC4H9O2H + HO2""",
)

entry(
    index = 1104,
    label = "TC4H9O2 + H2O2 <=> TC4H9O2H + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + H2O2 <=> TC4H9O2H + HO2""",
)

entry(
    index = 1105,
    label = "IC4H9O2 + CH2O <=> IC4H9O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(9000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH2O <=> IC4H9O2H + HCO""",
)

entry(
    index = 1106,
    label = "TC4H9O2 + CH2O <=> TC4H9O2H + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(9000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH2O <=> TC4H9O2H + HCO""",
)

entry(
    index = 1107,
    label = "IC4H9O2 + CH3O2 => IC4H9O + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH3O2 => IC4H9O + CH3O + O2""",
)

entry(
    index = 1108,
    label = "TC4H9O2 + CH3O2 => TC4H9O + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH3O2 => TC4H9O + CH3O + O2""",
)

entry(
    index = 1109,
    label = "IC4H9O2 + C2H5O2 => IC4H9O + C2H5O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H5O2 => IC4H9O + C2H5O + O2""",
)

entry(
    index = 1110,
    label = "TC4H9O2 + C2H5O2 => TC4H9O + C2H5O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H5O2 => TC4H9O + C2H5O + O2""",
)

entry(
    index = 1111,
    label = "IC4H9O2 + CH3CO3 => IC4H9O + CH3CO2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH3CO3 => IC4H9O + CH3CO2 + O2""",
)

entry(
    index = 1112,
    label = "TC4H9O2 + CH3CO3 => TC4H9O + CH3CO2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH3CO3 => TC4H9O + CH3CO2 + O2""",
)

entry(
    index = 1113,
    label = "IC4H9O2 + IC4H9O2 => O2 + IC4H9O + IC4H9O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + IC4H9O2 => O2 + IC4H9O + IC4H9O""",
)

entry(
    index = 1114,
    label = "IC4H9O2 + TC4H9O2 => IC4H9O + TC4H9O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + TC4H9O2 => IC4H9O + TC4H9O + O2""",
)

entry(
    index = 1115,
    label = "TC4H9O2 + TC4H9O2 => O2 + TC4H9O + TC4H9O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + TC4H9O2 => O2 + TC4H9O + TC4H9O""",
)

entry(
    index = 1116,
    label = "IC4H9O2 + PC4H9O2 => IC4H9O + PC4H9O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + PC4H9O2 => IC4H9O + PC4H9O + O2""",
)

entry(
    index = 1117,
    label = "TC4H9O2 + PC4H9O2 => TC4H9O + PC4H9O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + PC4H9O2 => TC4H9O + PC4H9O + O2""",
)

entry(
    index = 1118,
    label = "IC4H9O2 + SC4H9O2 => IC4H9O + SC4H9O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + SC4H9O2 => IC4H9O + SC4H9O + O2""",
)

entry(
    index = 1119,
    label = "TC4H9O2 + SC4H9O2 => TC4H9O + SC4H9O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + SC4H9O2 => TC4H9O + SC4H9O + O2""",
)

entry(
    index = 1120,
    label = "IC4H9O2 + NC3H7O2 => IC4H9O + NC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + NC3H7O2 => IC4H9O + NC3H7O + O2""",
)

entry(
    index = 1121,
    label = "TC4H9O2 + NC3H7O2 => TC4H9O + NC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + NC3H7O2 => TC4H9O + NC3H7O + O2""",
)

entry(
    index = 1122,
    label = "IC4H9O2 + IC3H7O2 => IC4H9O + IC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + IC3H7O2 => IC4H9O + IC3H7O + O2""",
)

entry(
    index = 1123,
    label = "TC4H9O2 + IC3H7O2 => TC4H9O + IC3H7O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + IC3H7O2 => TC4H9O + IC3H7O + O2""",
)

entry(
    index = 1124,
    label = "IC4H9O2 + HO2 => IC4H9O + OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + HO2 => IC4H9O + OH + O2""",
)

entry(
    index = 1125,
    label = "TC4H9O2 + HO2 => TC4H9O + OH + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + HO2 => TC4H9O + OH + O2""",
)

entry(
    index = 1126,
    label = "IC4H9O2 + CH3 <=> IC4H9O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH3 <=> IC4H9O + CH3O""",
)

entry(
    index = 1127,
    label = "IC4H9O2 + C2H5 <=> IC4H9O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H5 <=> IC4H9O + C2H5O""",
)

entry(
    index = 1128,
    label = "IC4H9O2 + IC3H7 <=> IC4H9O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + IC3H7 <=> IC4H9O + IC3H7O""",
)

entry(
    index = 1129,
    label = "IC4H9O2 + NC3H7 <=> IC4H9O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + NC3H7 <=> IC4H9O + NC3H7O""",
)

entry(
    index = 1130,
    label = "IC4H9O2 + PC4H9 <=> IC4H9O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + PC4H9 <=> IC4H9O + PC4H9O""",
)

entry(
    index = 1131,
    label = "IC4H9O2 + SC4H9 <=> IC4H9O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + SC4H9 <=> IC4H9O + SC4H9O""",
)

entry(
    index = 1132,
    label = "IC4H9O2 + IC4H9 <=> IC4H9O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + IC4H9 <=> IC4H9O + IC4H9O""",
)

entry(
    index = 1133,
    label = "IC4H9O2 + TC4H9 <=> IC4H9O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + TC4H9 <=> IC4H9O + TC4H9O""",
)

entry(
    index = 1134,
    label = "IC4H9O2 + C3H5-A <=> IC4H9O + C3H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C3H5-A <=> IC4H9O + C3H5O""",
)

entry(
    index = 1135,
    label = "IC4H9O2 + C4H71-3 <=> IC4H9O + C4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C4H71-3 <=> IC4H9O + C4H7O""",
)

entry(
    index = 1136,
    label = "IC4H9O2 + IC4H7 <=> IC4H9O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + IC4H7 <=> IC4H9O + IC4H7O""",
)

entry(
    index = 1137,
    label = "TC4H9O2 + CH3 <=> TC4H9O + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + CH3 <=> TC4H9O + CH3O""",
)

entry(
    index = 1138,
    label = "TC4H9O2 + C2H5 <=> TC4H9O + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C2H5 <=> TC4H9O + C2H5O""",
)

entry(
    index = 1139,
    label = "TC4H9O2 + IC3H7 <=> TC4H9O + IC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + IC3H7 <=> TC4H9O + IC3H7O""",
)

entry(
    index = 1140,
    label = "TC4H9O2 + NC3H7 <=> TC4H9O + NC3H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + NC3H7 <=> TC4H9O + NC3H7O""",
)

entry(
    index = 1141,
    label = "TC4H9O2 + PC4H9 <=> TC4H9O + PC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + PC4H9 <=> TC4H9O + PC4H9O""",
)

entry(
    index = 1142,
    label = "TC4H9O2 + SC4H9 <=> TC4H9O + SC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + SC4H9 <=> TC4H9O + SC4H9O""",
)

entry(
    index = 1143,
    label = "TC4H9O2 + IC4H9 <=> TC4H9O + IC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + IC4H9 <=> TC4H9O + IC4H9O""",
)

entry(
    index = 1144,
    label = "TC4H9O2 + TC4H9 <=> TC4H9O + TC4H9O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + TC4H9 <=> TC4H9O + TC4H9O""",
)

entry(
    index = 1145,
    label = "TC4H9O2 + C3H5-A <=> TC4H9O + C3H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C3H5-A <=> TC4H9O + C3H5O""",
)

entry(
    index = 1146,
    label = "TC4H9O2 + C4H71-3 <=> TC4H9O + C4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + C4H71-3 <=> TC4H9O + C4H7O""",
)

entry(
    index = 1147,
    label = "TC4H9O2 + IC4H7 <=> TC4H9O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 + IC4H7 <=> TC4H9O + IC4H7O""",
)

entry(
    index = 1148,
    label = "IC4H9O2 + C2H4 <=> IC4H9O2H + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H4 <=> IC4H9O2H + C2H3""",
)

entry(
    index = 1149,
    label = "IC4H9O2 + CH4 <=> IC4H9O2H + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.13e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH4 <=> IC4H9O2H + CH3""",
)

entry(
    index = 1150,
    label = "H2 + IC4H9O2 <=> H + IC4H9O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+13, 'cm^3/(mol*s)'), n=0, Ea=(26030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2 + IC4H9O2 <=> H + IC4H9O2H""",
)

entry(
    index = 1151,
    label = "IC4H9O2 + C2H6 <=> IC4H9O2H + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H6 <=> IC4H9O2H + C2H5""",
)

entry(
    index = 1152,
    label = "IC4H9O2 + C3H8 <=> IC4H9O2H + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C3H8 <=> IC4H9O2H + IC3H7""",
)

entry(
    index = 1153,
    label = "IC4H9O2 + C3H8 <=> IC4H9O2H + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+13, 'cm^3/(mol*s)'), n=0, Ea=(20460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C3H8 <=> IC4H9O2H + NC3H7""",
)

entry(
    index = 1154,
    label = "IC4H9O2 + CH3OH <=> IC4H9O2H + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + CH3OH <=> IC4H9O2H + CH2OH""",
)

entry(
    index = 1155,
    label = "IC4H9O2 + C2H5OH <=> IC4H9O2H + PC2H4OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(19360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H5OH <=> IC4H9O2H + PC2H4OH""",
)

entry(
    index = 1156,
    label = "IC4H9O2 + C2H5OH <=> IC4H9O2H + SC2H4OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 + C2H5OH <=> IC4H9O2H + SC2H4OH""",
)

entry(
    index = 1157,
    label = "IC4H9O2H <=> IC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2H <=> IC4H9O + OH""",
)

entry(
    index = 1158,
    label = "TC4H9O2H <=> TC4H9O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.95e+15, 's^-1'), n=0, Ea=(42540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2H <=> TC4H9O + OH""",
)

entry(
    index = 1159,
    label = "IC4H9O + HO2 <=> IC3H7CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O + HO2 <=> IC3H7CHO + H2O2""",
)

entry(
    index = 1160,
    label = "IC4H9O + OH <=> IC3H7CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O + OH <=> IC3H7CHO + H2O""",
)

entry(
    index = 1161,
    label = "IC4H9O + CH3 <=> IC3H7CHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O + CH3 <=> IC3H7CHO + CH4""",
)

entry(
    index = 1162,
    label = "IC4H9O + O <=> IC3H7CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O + O <=> IC3H7CHO + OH""",
)

entry(
    index = 1163,
    label = "IC4H9O + H <=> IC3H7CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O + H <=> IC3H7CHO + H2""",
)

entry(
    index = 1164,
    label = "IC4H9O <=> IC3H7CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+14, 's^-1'), n=0, Ea=(21500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O <=> IC3H7CHO + H""",
)

entry(
    index = 1165,
    label = "IC4H9O <=> CH2O + IC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+14, 's^-1'), n=0, Ea=(17500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O <=> CH2O + IC3H7""",
)

entry(
    index = 1166,
    label = "TC4H9O <=> CH3COCH3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.558e+22, 's^-1'), n=-2.548, Ea=(18650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O <=> CH3COCH3 + CH3""",
)

entry(
    index = 1167,
    label = "IC4H9O + O2 <=> IC3H7CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.93e+11, 'cm^3/(mol*s)'), n=0, Ea=(1660, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O + O2 <=> IC3H7CHO + HO2""",
)

entry(
    index = 1168,
    label = "TC4H9O + O2 <=> IC4H8O + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(4700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O + O2 <=> IC4H8O + HO2""",
)

entry(
    index = 1169,
    label = "IC4H8O <=> IC3H7CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.18e+13, 's^-1'), n=0, Ea=(52720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O <=> IC3H7CHO""",
)

entry(
    index = 1170,
    label = "IC4H8O + OH <=> IC3H6CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O + OH <=> IC3H6CHO + H2O""",
)

entry(
    index = 1171,
    label = "IC4H8O + H <=> IC3H6CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O + H <=> IC3H6CHO + H2""",
)

entry(
    index = 1172,
    label = "IC4H8O + HO2 <=> IC3H6CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O + HO2 <=> IC3H6CHO + H2O2""",
)

entry(
    index = 1173,
    label = "IC4H8O + CH3O2 <=> IC3H6CHO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(19000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O + CH3O2 <=> IC3H6CHO + CH3O2H""",
)

entry(
    index = 1174,
    label = "IC4H8O + CH3 <=> IC3H6CHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O + CH3 <=> IC3H6CHO + CH4""",
)

entry(
    index = 1175,
    label = "IC4H8O + O <=> IC3H6CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O + O <=> IC3H6CHO + OH""",
)

entry(
    index = 1176,
    label = "IC3H7CHO <=> TC3H6CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.304e+18, 's^-1'), n=-0.91, Ea=(92000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO <=> TC3H6CHO + H""",
)

entry(
    index = 1177,
    label = "IC3H7CHO <=> IC3H7 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.129e+17, 's^-1'), n=-0.03, Ea=(79760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO <=> IC3H7 + HCO""",
)

entry(
    index = 1178,
    label = "IC3H7CHO + HO2 <=> IC3H7CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + HO2 <=> IC3H7CO + H2O2""",
)

entry(
    index = 1179,
    label = "IC3H7CHO + HO2 <=> TC3H6CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+10, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + HO2 <=> TC3H6CHO + H2O2""",
)

entry(
    index = 1180,
    label = "IC3H7CHO + CH3 <=> IC3H7CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+12, 'cm^3/(mol*s)'), n=0, Ea=(8700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + CH3 <=> IC3H7CO + CH4""",
)

entry(
    index = 1181,
    label = "IC3H7CHO + O <=> IC3H7CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.18e+12, 'cm^3/(mol*s)'), n=0, Ea=(1389, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + O <=> IC3H7CO + OH""",
)

entry(
    index = 1182,
    label = "IC3H7CHO + O2 <=> IC3H7CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + O2 <=> IC3H7CO + HO2""",
)

entry(
    index = 1183,
    label = "IC3H7CHO + OH <=> IC3H7CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + OH <=> IC3H7CO + H2O""",
)

entry(
    index = 1184,
    label = "IC3H7CHO + OH <=> TC3H6CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.684e+12, 'cm^3/(mol*s)'), n=0, Ea=(-781, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + OH <=> TC3H6CHO + H2O""",
)

entry(
    index = 1185,
    label = "IC3H7CHO + H <=> IC3H7CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + H <=> IC3H7CO + H2""",
)

entry(
    index = 1186,
    label = "IC3H7CHO + OH <=> IC3H6CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + OH <=> IC3H6CHO + H2O""",
)

entry(
    index = 1187,
    label = "IC3H7CHO + HO2 <=> IC3H6CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27400, 'cm^3/(mol*s)'), n=2.55, Ea=(15500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + HO2 <=> IC3H6CHO + H2O2""",
)

entry(
    index = 1188,
    label = "IC3H7CHO + CH3O2 <=> IC3H6CHO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CHO + CH3O2 <=> IC3H6CHO + CH3O2H""",
)

entry(
    index = 1189,
    label = "IC3H7CO <=> IC3H7 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.869e+20, 's^-1'), n=-2.194, Ea=(14970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H7CO <=> IC3H7 + CO""",
)

entry(
    index = 1190,
    label = "IC3H6CHO <=> C3H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.031e+15, 's^-1'), n=-0.62, Ea=(23170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H6CHO <=> C3H6 + HCO""",
)

entry(
    index = 1191,
    label = "IC3H6CHO <=> C2H3CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.425e+13, 's^-1'), n=-0.27, Ea=(22470, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H6CHO <=> C2H3CHO + CH3""",
)

entry(
    index = 1192,
    label = "IC4H8OH <=> IC4H8 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.232e+14, 's^-1'), n=-0.562, Ea=(28050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OH <=> IC4H8 + OH""",
)

entry(
    index = 1193,
    label = "IO2C4H8OH <=> IC4H8OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.917e+21, 's^-1'), n=-2.347, Ea=(35790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IO2C4H8OH <=> IC4H8OH + O2""",
)

entry(
    index = 1194,
    label = "IO2C4H8OH => CH3COCH3 + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(18900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IO2C4H8OH => CH3COCH3 + CH2O + OH""",
)

entry(
    index = 1195,
    label = "IC4H9O2 <=> IC4H8O2H-I",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 <=> IC4H8O2H-I""",
)

entry(
    index = 1196,
    label = "TC4H9O2 <=> TC4H8O2H-I",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+11, 's^-1'), n=0, Ea=(34500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 <=> TC4H8O2H-I""",
)

entry(
    index = 1197,
    label = "IC4H9O2 <=> IC4H8O2H-T",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(29200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 <=> IC4H8O2H-T""",
)

entry(
    index = 1198,
    label = "IC4H9O2 <=> IC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.265e+35, 's^-1'), n=-7.22, Ea=(39490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H9O2 <=> IC4H8 + HO2""",
)

entry(
    index = 1199,
    label = "TC4H9O2 <=> IC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.612e+42, 's^-1'), n=-9.41, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H9O2 <=> IC4H8 + HO2""",
)

entry(
    index = 1200,
    label = "IC4H8OOH-IO2 <=> IC4H8O2H-I + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.44e+20, 's^-1'), n=-1.627, Ea=(35690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-IO2 <=> IC4H8O2H-I + O2""",
)

entry(
    index = 1201,
    label = "TC4H8OOH-IO2 <=> TC4H8O2H-I + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.167e+22, 's^-1'), n=-2.257, Ea=(37800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H8OOH-IO2 <=> TC4H8O2H-I + O2""",
)

entry(
    index = 1202,
    label = "IC4H8OOH-TO2 <=> IC4H8O2H-T + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.266e+27, 's^-1'), n=-3.233, Ea=(39640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-TO2 <=> IC4H8O2H-T + O2""",
)

entry(
    index = 1203,
    label = "IC4H8OOH-IO2 <=> IC4KETII + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 's^-1'), n=0, Ea=(21400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-IO2 <=> IC4KETII + OH""",
)

entry(
    index = 1204,
    label = "IC4H8OOH-TO2 <=> IC4KETIT + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 's^-1'), n=0, Ea=(31500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-TO2 <=> IC4KETIT + OH""",
)

entry(
    index = 1205,
    label = "TC4H8OOH-IO2 <=> TIC4H7Q2-I",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H8OOH-IO2 <=> TIC4H7Q2-I""",
)

entry(
    index = 1206,
    label = "TIC4H7Q2-I <=> IC4H7OOH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.94e+20, 's^-1'), n=-2.19, Ea=(22590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TIC4H7Q2-I <=> IC4H7OOH + HO2""",
)

entry(
    index = 1207,
    label = "IC4H8OOH-IO2 <=> IIC4H7Q2-I",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+10, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-IO2 <=> IIC4H7Q2-I""",
)

entry(
    index = 1208,
    label = "IC4H8OOH-IO2 <=> IIC4H7Q2-T",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(29200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-IO2 <=> IIC4H7Q2-T""",
)

entry(
    index = 1209,
    label = "IC4H8OOH-TO2 <=> TIC4H7Q2-I",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(34500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OOH-TO2 <=> TIC4H7Q2-I""",
)

entry(
    index = 1210,
    label = "IIC4H7Q2-I <=> AC3H5OOH + CH2O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.631e+19, 's^-1'), n=-1.74, Ea=(38310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IIC4H7Q2-I <=> AC3H5OOH + CH2O2H""",
)

entry(
    index = 1211,
    label = "IIC4H7Q2-T <=> IC4H7OOH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.076e+17, 's^-1'), n=-1.56, Ea=(18390, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IIC4H7Q2-T <=> IC4H7OOH + HO2""",
)

entry(
    index = 1212,
    label = "CH2O2H <=> CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+14, 's^-1'), n=0, Ea=(1500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2O2H <=> CH2O + OH""",
)

entry(
    index = 1213,
    label = "IC4KETII => CH2O + C2H5CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4KETII => CH2O + C2H5CO + OH""",
)

entry(
    index = 1214,
    label = "IC4KETIT => CH3COCH3 + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.5e+15, 's^-1'), n=0, Ea=(42540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4KETIT => CH3COCH3 + HCO + OH""",
)

entry(
    index = 1215,
    label = "TC4H8O2H-I <=> IC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.073e+20, 's^-1'), n=-2.085, Ea=(19390, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H8O2H-I <=> IC4H8 + HO2""",
)

entry(
    index = 1216,
    label = "IC4H8O2H-T <=> IC4H8 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.526e+16, 's^-1'), n=-1.109, Ea=(17560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O2H-T <=> IC4H8 + HO2""",
)

entry(
    index = 1217,
    label = "IC4H8O2H-I => CC4H8O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(19500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O2H-I => CC4H8O + OH""",
)

entry(
    index = 1218,
    label = "IC4H8O2H-T => IC4H8O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.38e+12, 's^-1'), n=0, Ea=(14800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O2H-T => IC4H8O + OH""",
)

entry(
    index = 1219,
    label = "TC4H8O2H-I => IC4H8O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 's^-1'), n=0, Ea=(17000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H8O2H-I => IC4H8O + OH""",
)

entry(
    index = 1220,
    label = "IC4H8O2H-I => OH + CH2O + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.451e+15, 's^-1'), n=-0.68, Ea=(29170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8O2H-I => OH + CH2O + C3H6""",
)

entry(
    index = 1221,
    label = "IC4H8 <=> C3H5-T + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.92e+66, 's^-1'), n=-14.22, Ea=(128100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 <=> C3H5-T + CH3""",
)

entry(
    index = 1222,
    label = "IC4H8 <=> IC4H7 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.07e+55, 's^-1'), n=-11.49, Ea=(114300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 <=> IC4H7 + H""",
)

entry(
    index = 1223,
    label = "IC4H8 + H <=> C3H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.68e+33, 'cm^3/(mol*s)'),
        n = -5.72,
        Ea = (20000, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H8 + H <=> C3H6 + CH3""",
)

entry(
    index = 1224,
    label = "IC4H8 + H <=> IC4H7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(340000, 'cm^3/(mol*s)'), n=2.5, Ea=(2492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + H <=> IC4H7 + H2""",
)

entry(
    index = 1225,
    label = "IC4H8 + O => CH2CO + CH3 + CH3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.33e+07, 'cm^3/(mol*s)'), n=1.76, Ea=(76, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + O => CH2CO + CH3 + CH3""",
)

entry(
    index = 1226,
    label = "IC4H8 + O => IC3H6CO + H + H",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.66e+07, 'cm^3/(mol*s)'), n=1.76, Ea=(76, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + O => IC3H6CO + H + H""",
)

entry(
    index = 1227,
    label = "IC4H8 + O <=> IC4H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.206e+11, 'cm^3/(mol*s)'),
        n = 0.7,
        Ea = (7633, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H8 + O <=> IC4H7 + OH""",
)

entry(
    index = 1228,
    label = "IC4H8 + CH3 <=> IC4H7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.42, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + CH3 <=> IC4H7 + CH4""",
)

entry(
    index = 1229,
    label = "IC4H8 + HO2 <=> IC4H7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19280, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + HO2 <=> IC4H7 + H2O2""",
)

entry(
    index = 1230,
    label = "IC4H8 + O2CHO <=> IC4H7 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19280, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + O2CHO <=> IC4H7 + HO2CHO""",
)

entry(
    index = 1231,
    label = "IC4H8 + O2 <=> IC4H7 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(39900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + O2 <=> IC4H7 + HO2""",
)

entry(
    index = 1232,
    label = "IC4H8 + C3H5-A <=> IC4H7 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(20500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + C3H5-A <=> IC4H7 + C3H6""",
)

entry(
    index = 1233,
    label = "IC4H8 + C3H5-S <=> IC4H7 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(20500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + C3H5-S <=> IC4H7 + C3H6""",
)

entry(
    index = 1234,
    label = "IC4H8 + C3H5-T <=> IC4H7 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(20500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + C3H5-T <=> IC4H7 + C3H6""",
)

entry(
    index = 1235,
    label = "IC4H8 + OH <=> IC4H7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.2e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + OH <=> IC4H7 + H2O""",
)

entry(
    index = 1236,
    label = "IC4H8 + O <=> IC3H7 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.58e+07, 'cm^3/(mol*s)'),
        n = 1.76,
        Ea = (-1216, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H8 + O <=> IC3H7 + HCO""",
)

entry(
    index = 1237,
    label = "IC4H8 + CH3O2 <=> IC4H7 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(19280, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + CH3O2 <=> IC4H7 + CH3O2H""",
)

entry(
    index = 1238,
    label = "IC4H8 + HO2 <=> IC4H8O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.29e+12, 'cm^3/(mol*s)'), n=0, Ea=(13340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + HO2 <=> IC4H8O + OH""",
)

entry(
    index = 1239,
    label = "IC4H7 + O2 <=> IC3H5CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.47e+13, 'cm^3/(mol*s)'),
        n = -0.45,
        Ea = (23020, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H7 + O2 <=> IC3H5CHO + OH""",
)

entry(
    index = 1240,
    label = "IC4H7 + O2 <=> CH3COCH2 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.14e+15, 'cm^3/(mol*s)'),
        n = -1.21,
        Ea = (21050, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H7 + O2 <=> CH3COCH2 + CH2O""",
)

entry(
    index = 1241,
    label = "IC4H7 + O2 => C3H4-A + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (7.29e+29, 'cm^3/(mol*s)'),
        n = -5.71,
        Ea = (21450, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H7 + O2 => C3H4-A + CH2O + OH""",
)

entry(
    index = 1242,
    label = "IC4H7 + O <=> IC3H5CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7 + O <=> IC3H5CHO + H""",
)

entry(
    index = 1243,
    label = "IC4H7 <=> C3H4-A + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.23e+47, 's^-1'), n=-9.74, Ea=(74260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7 <=> C3H4-A + CH3""",
)

entry(
    index = 1244,
    label = "CH3O2 + IC4H7 <=> CH3O + IC4H7O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + IC4H7 <=> CH3O + IC4H7O""",
)

entry(
    index = 1245,
    label = "IC4H7 + HO2 <=> IC4H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7 + HO2 <=> IC4H7O + OH""",
)

entry(
    index = 1246,
    label = "IC4H7O <=> C3H5-T + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.925e+21, 's^-1'), n=-2.391, Ea=(35590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O <=> C3H5-T + CH2O""",
)

entry(
    index = 1247,
    label = "IC4H7O <=> IC4H6OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.391e+11, 's^-1'), n=0, Ea=(15600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O <=> IC4H6OH""",
)

entry(
    index = 1248,
    label = "IC4H7O <=> IC3H5CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 's^-1'), n=0, Ea=(29100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O <=> IC3H5CHO + H""",
)

entry(
    index = 1249,
    label = "IC4H6OH + H2 <=> IC4H7OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(21600, 'cm^3/(mol*s)'), n=2.38, Ea=(18990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H6OH + H2 <=> IC4H7OH + H""",
)

entry(
    index = 1250,
    label = "IC4H6OH + HO2 <=> IC4H7OH + O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.57e+13, 'cm^3/(mol*s)'),
        n = -0.315,
        Ea = (862, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H6OH + HO2 <=> IC4H7OH + O2""",
)

entry(
    index = 1251,
    label = "IC4H6OH + CH2O <=> IC4H7OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (6.3e+08, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (18190, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC4H6OH + CH2O <=> IC4H7OH + HCO""",
)

entry(
    index = 1252,
    label = "IC4H6OH + IC4H8 <=> IC4H7OH + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(470, 'cm^3/(mol*s)'), n=3.3, Ea=(19840, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H6OH + IC4H8 <=> IC4H7OH + IC4H7""",
)

entry(
    index = 1253,
    label = "IC4H7OH <=> IC4H6OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.902e+16, 's^-1'), n=-0.4, Ea=(89850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7OH <=> IC4H6OH + H""",
)

entry(
    index = 1254,
    label = "IC4H7OH + HO2 <=> IC4H6OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7644, 'cm^3/(mol*s)'), n=2.712, Ea=(13930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7OH + HO2 <=> IC4H6OH + H2O2""",
)

entry(
    index = 1255,
    label = "IC4H6OH <=> C3H4-A + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.244e+19, 's^-1'), n=-1.859, Ea=(57050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H6OH <=> C3H4-A + CH2OH""",
)

entry(
    index = 1256,
    label = "IC4H7O + O2 <=> IC3H5CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+10, 'cm^3/(mol*s)'), n=0, Ea=(1649, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + O2 <=> IC3H5CHO + HO2""",
)

entry(
    index = 1257,
    label = "IC4H7O + HO2 <=> IC3H5CHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + HO2 <=> IC3H5CHO + H2O2""",
)

entry(
    index = 1258,
    label = "IC4H7O + CH3 <=> IC3H5CHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + CH3 <=> IC3H5CHO + CH4""",
)

entry(
    index = 1259,
    label = "IC4H7O + O <=> IC3H5CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + O <=> IC3H5CHO + OH""",
)

entry(
    index = 1260,
    label = "IC4H7O + OH <=> IC3H5CHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.81e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + OH <=> IC3H5CHO + H2O""",
)

entry(
    index = 1261,
    label = "IC4H7O + H <=> IC3H5CHO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + H <=> IC3H5CHO + H2""",
)

entry(
    index = 1262,
    label = "IC3H5CHO + OH <=> IC3H5CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is IC3H5CHO + OH <=> IC3H5CO + H2O""",
)

entry(
    index = 1263,
    label = "IC3H5CHO + HO2 <=> IC3H5CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5CHO + HO2 <=> IC3H5CO + H2O2""",
)

entry(
    index = 1264,
    label = "IC3H5CHO + CH3 <=> IC3H5CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+12, 'cm^3/(mol*s)'), n=0, Ea=(8700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5CHO + CH3 <=> IC3H5CO + CH4""",
)

entry(
    index = 1265,
    label = "IC3H5CHO + O <=> IC3H5CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.18e+12, 'cm^3/(mol*s)'), n=0, Ea=(1389, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5CHO + O <=> IC3H5CO + OH""",
)

entry(
    index = 1266,
    label = "IC3H5CHO + O2 <=> IC3H5CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(40700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5CHO + O2 <=> IC3H5CO + HO2""",
)

entry(
    index = 1267,
    label = "IC3H5CHO + H <=> IC3H5CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(2600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5CHO + H <=> IC3H5CO + H2""",
)

entry(
    index = 1268,
    label = "IC3H5CO <=> C3H5-T + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.278e+20, 's^-1'), n=-1.89, Ea=(34460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5CO <=> C3H5-T + CO""",
)

entry(
    index = 1269,
    label = "TC3H6CHO + HO2 <=> TC3H6OCHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + HO2 <=> TC3H6OCHO + OH""",
)

entry(
    index = 1270,
    label = "TC3H6OCHO <=> CH3COCH3 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+13, 's^-1'), n=0, Ea=(9700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6OCHO <=> CH3COCH3 + HCO""",
)

entry(
    index = 1271,
    label = "TC3H6CHO <=> IC3H5CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.325e+14, 's^-1'), n=0.01, Ea=(39340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO <=> IC3H5CHO + H""",
)

entry(
    index = 1272,
    label = "TC3H6CHO <=> IC3H6CO + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.086e+14, 's^-1'), n=-0.072, Ea=(42410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO <=> IC3H6CO + H""",
)

entry(
    index = 1273,
    label = "TC3H6CHO + H2 <=> IC3H7CHO + H",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (216000, 'cm^3/(mol*s)'),
        n = 2.38,
        Ea = (18990, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + H2 <=> IC3H7CHO + H""",
)

entry(
    index = 1274,
    label = "IC4H7OOH <=> IC4H7O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.4e+15, 's^-1'), n=0, Ea=(45550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7OOH <=> IC4H7O + OH""",
)

entry(
    index = 1275,
    label = "IC4H7OH <=> IC4H7O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.969e+16, 's^-1'), n=-0.56, Ea=(105900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7OH <=> IC4H7O + H""",
)

entry(
    index = 1276,
    label = "IC4H8OH <=> IC4H7OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.708e+12, 's^-1'), n=0.277, Ea=(38850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8OH <=> IC4H7OH + H""",
)

entry(
    index = 1277,
    label = "IC4H7O + H2 <=> IC4H7OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.05e+06, 'cm^3/(mol*s)'), n=2, Ea=(17830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + H2 <=> IC4H7OH + H""",
)

entry(
    index = 1278,
    label = "IC4H7OH <=> IC4H7 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.31e+16, 's^-1'), n=-0.41, Ea=(79700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7OH <=> IC4H7 + OH""",
)

entry(
    index = 1279,
    label = "IC4H7O + CH2O <=> IC4H7OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+11, 'cm^3/(mol*s)'), n=0, Ea=(1280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + CH2O <=> IC4H7OH + HCO""",
)

entry(
    index = 1280,
    label = "TC3H6CHO + CH2O <=> IC3H7CHO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.52e+08, 'cm^3/(mol*s)'),
        n = 1.9,
        Ea = (18190, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + CH2O <=> IC3H7CHO + HCO""",
)

entry(
    index = 1281,
    label = "TC3H6CHO + IC4H8 <=> IC3H7CHO + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(470, 'cm^3/(mol*s)'), n=3.3, Ea=(19840, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + IC4H8 <=> IC3H7CHO + IC4H7""",
)

entry(
    index = 1282,
    label = "IC3H6CO + OH <=> IC3H7 + CO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.73e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H6CO + OH <=> IC3H7 + CO2""",
)

entry(
    index = 1283,
    label = "TC3H6OHCHO <=> TC3H6CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.99e+20, 's^-1'), n=-1.46, Ea=(87480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6OHCHO <=> TC3H6CHO + OH""",
)

entry(
    index = 1284,
    label = "TC3H6OHCHO <=> TC3H6OH + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.164e+23, 's^-1'), n=-1.9, Ea=(76850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6OHCHO <=> TC3H6OH + HCO""",
)

entry(
    index = 1285,
    label = "TC3H6OH <=> CH3COCH3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+13, 's^-1'), n=0, Ea=(21860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6OH <=> CH3COCH3 + H""",
)

entry(
    index = 1286,
    label = "TC3H6OH <=> IC3H5OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.196e+15, 's^-1'), n=-0.66, Ea=(40340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6OH <=> IC3H5OH + H""",
)

entry(
    index = 1287,
    label = "IC3H5OH <=> C3H5-T + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.369e+19, 's^-1'), n=-0.94, Ea=(109100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5OH <=> C3H5-T + OH""",
)

entry(
    index = 1288,
    label = "TC3H6O2CHO <=> TC3H6CHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.458e+25, 's^-1'), n=-4.065, Ea=(27080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6O2CHO <=> TC3H6CHO + O2""",
)

entry(
    index = 1289,
    label = "TC3H6O2CHO <=> IC3H5O2HCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(29880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6O2CHO <=> IC3H5O2HCHO""",
)

entry(
    index = 1290,
    label = "TC3H6O2CHO <=> TC3H6O2HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(25750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6O2CHO <=> TC3H6O2HCO""",
)

entry(
    index = 1291,
    label = "IC3H5O2HCHO <=> IC3H5CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.944e+20, 's^-1'), n=-2.44, Ea=(15030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H5O2HCHO <=> IC3H5CHO + HO2""",
)

entry(
    index = 1292,
    label = "TC3H6O2HCO => CH3COCH3 + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4.244e+18, 's^-1'), n=-1.43, Ea=(4800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6O2HCO => CH3COCH3 + CO + OH""",
)

entry(
    index = 1293,
    label = "TC3H6OH + O2 <=> CH3COCH3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.23e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6OH + O2 <=> CH3COCH3 + HO2""",
)

entry(
    index = 1294,
    label = "IC3H6CO + OH <=> TC3H6OH + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC3H6CO + OH <=> TC3H6OH + CO""",
)

entry(
    index = 1295,
    label = "TC3H6CHO + O2 <=> IC3H5CHO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.725e-19, 'cm^3/(mol*s)'), n=0, Ea=(7240, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + O2 <=> IC3H5CHO + HO2""",
)

entry(
    index = 1296,
    label = "TC3H6CHO + O2 => CH3COCH3 + CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.62e-20, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + O2 => CH3COCH3 + CO + OH""",
)

entry(
    index = 1297,
    label = "TC3H6CHO + HO2 <=> IC3H7CHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.675e+12, 'cm^3/(mol*s)'), n=0, Ea=(1310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + HO2 <=> IC3H7CHO + O2""",
)

entry(
    index = 1298,
    label = "TC3H6CHO + CH3 <=> IC3H5CHO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.01e+12, 'cm^3/(mol*s)'),
        n = -0.32,
        Ea = (-131, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is TC3H6CHO + CH3 <=> IC3H5CHO + CH4""",
)

entry(
    index = 1299,
    label = "TC4H8CHO <=> IC3H5CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(26290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H8CHO <=> IC3H5CHO + CH3""",
)

entry(
    index = 1300,
    label = "TC4H8CHO <=> IC4H8 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.52e+12, 's^-1'), n=0, Ea=(20090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is TC4H8CHO <=> IC4H8 + HCO""",
)

entry(
    index = 1301,
    label = "O2C4H8CHO <=> TC4H8CHO + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.519e+19, 's^-1'), n=-1.44, Ea=(34510, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C4H8CHO <=> TC4H8CHO + O2""",
)

entry(
    index = 1302,
    label = "O2C4H8CHO <=> O2HC4H8CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.16e+11, 's^-1'), n=0, Ea=(15360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C4H8CHO <=> O2HC4H8CO""",
)

entry(
    index = 1303,
    label = "O2HC4H8CO <=> IC4H8O2H-T + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.299e+22, 's^-1'), n=-2.72, Ea=(11760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2HC4H8CO <=> IC4H8O2H-T + CO""",
)

entry(
    index = 1304,
    label = "IC4H7O + IC4H8 <=> IC4H7OH + IC4H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.7e+11, 'cm^3/(mol*s)'), n=0, Ea=(4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7O + IC4H8 <=> IC4H7OH + IC4H7""",
)

entry(
    index = 1305,
    label = "IC4H6OH + HO2 => CH2CCH2OH + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.446e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H6OH + HO2 => CH2CCH2OH + CH2O + OH""",
)

entry(
    index = 1306,
    label = "IC4H8 + CH2CCH2OH <=> IC4H7 + C3H5OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.94e+11, 'cm^3/(mol*s)'), n=0, Ea=(20500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H8 + CH2CCH2OH <=> IC4H7 + C3H5OH""",
)

entry(
    index = 1307,
    label = "C3H5OH + HO2 <=> CH2CCH2OH + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.761e+09, 'cm^3/(mol*s)'),
        n = 0.28,
        Ea = (22590, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C3H5OH + HO2 <=> CH2CCH2OH + H2O2""",
)

entry(
    index = 1308,
    label = "C3H5OH + OH <=> CH2CCH2OH + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.06e+12, 'cm^3/(mol*s)'), n=0, Ea=(5960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5OH + OH <=> CH2CCH2OH + H2O""",
)

entry(
    index = 1309,
    label = "C3H5OH + H <=> CH2CCH2OH + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(390000, 'cm^3/(mol*s)'), n=2.5, Ea=(5821, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5OH + H <=> CH2CCH2OH + H2""",
)

entry(
    index = 1310,
    label = "C3H5OH + O2 <=> CH2CCH2OH + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(60690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5OH + O2 <=> CH2CCH2OH + HO2""",
)

entry(
    index = 1311,
    label = "C3H5OH + CH3 <=> CH2CCH2OH + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(8030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5OH + CH3 <=> CH2CCH2OH + CH4""",
)

entry(
    index = 1312,
    label = "IC4H7OH <=> CH2CCH2OH + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.247e+20, 's^-1'), n=-0.98, Ea=(98570, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is IC4H7OH <=> CH2CCH2OH + CH3""",
)

entry(
    index = 1313,
    label = "C3H5OH <=> CH2CCH2OH + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.835e+19, 's^-1'), n=-1.05, Ea=(111100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5OH <=> CH2CCH2OH + H""",
)

entry(
    index = 1314,
    label = "CH2CCH2OH + O2 => CH2OH + CO + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4.335e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CCH2OH + O2 => CH2OH + CO + CH2O""",
)

entry(
    index = 1315,
    label = "CH2CCH2OH <=> C2H2 + CH2OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.163e+40, 's^-1'), n=-8.31, Ea=(45110, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CCH2OH <=> C2H2 + CH2OH""",
)

entry(
    index = 1316,
    label = "CH2CCH2OH <=> C3H4-A + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.697e+16, 's^-1'), n=-1.11, Ea=(42580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH2CCH2OH <=> C3H4-A + OH""",
)

entry(
    index = 1317,
    label = "NC5H12 <=> C5H11-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.375e+17, 's^-1'), n=-0.36, Ea=(101200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 <=> C5H11-1 + H""",
)

entry(
    index = 1318,
    label = "NC5H12 <=> C5H11-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.325e+18, 's^-1'), n=-0.763, Ea=(98800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 <=> C5H11-2 + H""",
)

entry(
    index = 1319,
    label = "NC5H12 <=> C5H11-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.602e+18, 's^-1'), n=-0.758, Ea=(98790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 <=> C5H11-3 + H""",
)

entry(
    index = 1320,
    label = "NC5H12 <=> CH3 + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.101e+22, 's^-1'), n=-1.862, Ea=(89430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 <=> CH3 + PC4H9""",
)

entry(
    index = 1321,
    label = "NC5H12 <=> NC3H7 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.082e+24, 's^-1'), n=-2.269, Ea=(88440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 <=> NC3H7 + C2H5""",
)

entry(
    index = 1322,
    label = "NC5H12 + H <=> C5H11-1 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(188000, 'cm^3/(mol*s)'), n=2.75, Ea=(6280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + H <=> C5H11-1 + H2""",
)

entry(
    index = 1323,
    label = "NC5H12 + H <=> C5H11-2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + H <=> C5H11-2 + H2""",
)

entry(
    index = 1324,
    label = "NC5H12 + H <=> C5H11-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + H <=> C5H11-3 + H2""",
)

entry(
    index = 1325,
    label = "NC5H12 + OH <=> C5H11-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1590, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + OH <=> C5H11-1 + H2O""",
)

entry(
    index = 1326,
    label = "NC5H12 + OH <=> C5H11-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (9.34e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + OH <=> C5H11-2 + H2O""",
)

entry(
    index = 1327,
    label = "NC5H12 + OH <=> C5H11-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + OH <=> C5H11-3 + H2O""",
)

entry(
    index = 1328,
    label = "NC5H12 + O <=> C5H11-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.092e+06, 'cm^3/(mol*s)'),
        n = 2.424,
        Ea = (4766, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O <=> C5H11-1 + OH""",
)

entry(
    index = 1329,
    label = "NC5H12 + O <=> C5H11-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.189e+06, 'cm^3/(mol*s)'),
        n = 2.439,
        Ea = (2846, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O <=> C5H11-2 + OH""",
)

entry(
    index = 1330,
    label = "NC5H12 + O <=> C5H11-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (594600, 'cm^3/(mol*s)'),
        n = 2.439,
        Ea = (2846, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O <=> C5H11-3 + OH""",
)

entry(
    index = 1331,
    label = "NC5H12 + CH3 <=> C5H11-1 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.904, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3 <=> C5H11-1 + CH4""",
)

entry(
    index = 1332,
    label = "NC5H12 + CH3 <=> C5H11-2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (168000, 'cm^3/(mol*s)'),
        n = 2.133,
        Ea = (7574, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3 <=> C5H11-2 + CH4""",
)

entry(
    index = 1333,
    label = "NC5H12 + CH3 <=> C5H11-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(84000, 'cm^3/(mol*s)'), n=2.133, Ea=(7574, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3 <=> C5H11-3 + CH4""",
)

entry(
    index = 1334,
    label = "NC5H12 + HO2 <=> C5H11-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40.8, 'cm^3/(mol*s)'), n=3.59, Ea=(17160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + HO2 <=> C5H11-1 + H2O2""",
)

entry(
    index = 1335,
    label = "NC5H12 + HO2 <=> C5H11-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(126.4, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + HO2 <=> C5H11-2 + H2O2""",
)

entry(
    index = 1336,
    label = "NC5H12 + HO2 <=> C5H11-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(63.2, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + HO2 <=> C5H11-3 + H2O2""",
)

entry(
    index = 1337,
    label = "NC5H12 + CH3O2 <=> C5H11-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40.8, 'cm^3/(mol*s)'), n=3.59, Ea=(17160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3O2 <=> C5H11-1 + CH3O2H""",
)

entry(
    index = 1338,
    label = "NC5H12 + CH3O2 <=> C5H11-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(126.4, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3O2 <=> C5H11-2 + CH3O2H""",
)

entry(
    index = 1339,
    label = "NC5H12 + CH3O2 <=> C5H11-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(63.2, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3O2 <=> C5H11-3 + CH3O2H""",
)

entry(
    index = 1340,
    label = "NC5H12 + C2H5 <=> C5H11-1 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(13400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C2H5 <=> C5H11-1 + C2H6""",
)

entry(
    index = 1341,
    label = "NC5H12 + C2H5 <=> C5H11-2 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C2H5 <=> C5H11-2 + C2H6""",
)

entry(
    index = 1342,
    label = "NC5H12 + C2H5 <=> C5H11-3 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C2H5 <=> C5H11-3 + C2H6""",
)

entry(
    index = 1343,
    label = "NC5H12 + C2H3 <=> C5H11-1 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(18000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C2H3 <=> C5H11-1 + C2H4""",
)

entry(
    index = 1344,
    label = "NC5H12 + C2H3 <=> C5H11-2 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C2H3 <=> C5H11-2 + C2H4""",
)

entry(
    index = 1345,
    label = "NC5H12 + C2H3 <=> C5H11-3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C2H3 <=> C5H11-3 + C2H4""",
)

entry(
    index = 1346,
    label = "NC5H12 + C5H11-1 <=> C5H11-2 + NC5H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C5H11-1 <=> C5H11-2 + NC5H12""",
)

entry(
    index = 1347,
    label = "NC5H12 + C5H11-1 <=> C5H11-3 + NC5H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C5H11-1 <=> C5H11-3 + NC5H12""",
)

entry(
    index = 1348,
    label = "NC5H12 + C5H11-2 <=> C5H11-3 + NC5H12",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(12300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + C5H11-2 <=> C5H11-3 + NC5H12""",
)

entry(
    index = 1349,
    label = "NC5H12 + O2CHO <=> C5H11-1 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.68e+13, 'cm^3/(mol*s)'), n=0, Ea=(20440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2CHO <=> C5H11-1 + HO2CHO""",
)

entry(
    index = 1350,
    label = "NC5H12 + O2CHO <=> C5H11-2 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2CHO <=> C5H11-2 + HO2CHO""",
)

entry(
    index = 1351,
    label = "NC5H12 + O2CHO <=> C5H11-3 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2CHO <=> C5H11-3 + HO2CHO""",
)

entry(
    index = 1352,
    label = "NC5H12 + CH3O <=> C5H11-1 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3O <=> C5H11-1 + CH3OH""",
)

entry(
    index = 1353,
    label = "NC5H12 + CH3O <=> C5H11-2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3O <=> C5H11-2 + CH3OH""",
)

entry(
    index = 1354,
    label = "NC5H12 + CH3O <=> C5H11-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + CH3O <=> C5H11-3 + CH3OH""",
)

entry(
    index = 1355,
    label = "NC5H12 + O2 <=> C5H11-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(52800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2 <=> C5H11-1 + HO2""",
)

entry(
    index = 1356,
    label = "NC5H12 + O2 <=> C5H11-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(50160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2 <=> C5H11-2 + HO2""",
)

entry(
    index = 1357,
    label = "NC5H12 + O2 <=> C5H11-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H12 + O2 <=> C5H11-3 + HO2""",
)

entry(
    index = 1358,
    label = "C5H11-1 <=> C2H4 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.205e+12, 's^-1'), n=0.451, Ea=(29430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 <=> C2H4 + NC3H7""",
)

entry(
    index = 1359,
    label = "C5H11-1 <=> H + C5H10-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.354e+11, 's^-1'), n=0.608, Ea=(35640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 <=> H + C5H10-1""",
)

entry(
    index = 1360,
    label = "C5H11-1 <=> C5H11-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.875e+09, 's^-1'), n=0.353, Ea=(19760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 <=> C5H11-2""",
)

entry(
    index = 1361,
    label = "C5H11-2 <=> C3H6 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.223e+12, 's^-1'), n=0.635, Ea=(29360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 <=> C3H6 + C2H5""",
)

entry(
    index = 1362,
    label = "C5H11-2 <=> C5H10-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.353e+10, 's^-1'), n=1.011, Ea=(36680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 <=> C5H10-1 + H""",
)

entry(
    index = 1363,
    label = "C5H11-2 <=> C5H10-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.988e+11, 's^-1'), n=0.41, Ea=(35220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 <=> C5H10-2 + H""",
)

entry(
    index = 1364,
    label = "C5H11-3 <=> C4H8-1 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.343e+10, 's^-1'), n=1.119, Ea=(30460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 <=> C4H8-1 + CH3""",
)

entry(
    index = 1365,
    label = "C5H11-3 <=> C5H10-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.277e+11, 's^-1'), n=0.405, Ea=(35230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 <=> C5H10-2 + H""",
)

entry(
    index = 1366,
    label = "C5H10-1 <=> C2H5 + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.864e+21, 's^-1'), n=-2.086, Ea=(75060, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 <=> C2H5 + C3H5-A""",
)

entry(
    index = 1367,
    label = "C5H10-1 + H <=> C5H91-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + H <=> C5H91-3 + H2""",
)

entry(
    index = 1368,
    label = "C5H10-1 + H <=> C5H91-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + H <=> C5H91-4 + H2""",
)

entry(
    index = 1369,
    label = "C5H10-1 + H <=> C5H91-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665000, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + H <=> C5H91-5 + H2""",
)

entry(
    index = 1370,
    label = "C5H10-1 + O <=> C5H91-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(660000, 'cm^3/(mol*s)'), n=2.43, Ea=(1210, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + O <=> C5H91-3 + OH""",
)

entry(
    index = 1371,
    label = "C5H10-1 + O <=> C5H91-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(551000, 'cm^3/(mol*s)'), n=2.45, Ea=(2830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + O <=> C5H91-4 + OH""",
)

entry(
    index = 1372,
    label = "C5H10-1 + O <=> C5H91-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(980000, 'cm^3/(mol*s)'), n=2.43, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + O <=> C5H91-5 + OH""",
)

entry(
    index = 1373,
    label = "C5H10-1 + OH <=> C5H91-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + OH <=> C5H91-3 + H2O""",
)

entry(
    index = 1374,
    label = "C5H10-1 + OH <=> C5H91-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + OH <=> C5H91-4 + H2O""",
)

entry(
    index = 1375,
    label = "C5H10-1 + OH <=> C5H91-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + OH <=> C5H91-5 + H2O""",
)

entry(
    index = 1376,
    label = "C5H10-1 + CH3 <=> C5H91-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3 <=> C5H91-3 + CH4""",
)

entry(
    index = 1377,
    label = "C5H10-1 + CH3 <=> C5H91-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3 <=> C5H91-4 + CH4""",
)

entry(
    index = 1378,
    label = "C5H10-1 + CH3 <=> C5H91-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4521, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3 <=> C5H91-5 + CH4""",
)

entry(
    index = 1379,
    label = "C5H10-1 + O2 <=> C5H91-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(37220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + O2 <=> C5H91-3 + HO2""",
)

entry(
    index = 1380,
    label = "C5H10-1 + O2 <=> C5H91-4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(49640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + O2 <=> C5H91-4 + HO2""",
)

entry(
    index = 1381,
    label = "C5H10-1 + O2 <=> C5H91-5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + O2 <=> C5H91-5 + HO2""",
)

entry(
    index = 1382,
    label = "C5H10-1 + HO2 <=> C5H91-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + HO2 <=> C5H91-3 + H2O2""",
)

entry(
    index = 1383,
    label = "C5H10-1 + HO2 <=> C5H91-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + HO2 <=> C5H91-4 + H2O2""",
)

entry(
    index = 1384,
    label = "C5H10-1 + HO2 <=> C5H91-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + HO2 <=> C5H91-5 + H2O2""",
)

entry(
    index = 1385,
    label = "C5H10-1 + CH3O2 <=> C5H91-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3O2 <=> C5H91-3 + CH3O2H""",
)

entry(
    index = 1386,
    label = "C5H10-1 + CH3O2 <=> C5H91-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3O2 <=> C5H91-4 + CH3O2H""",
)

entry(
    index = 1387,
    label = "C5H10-1 + CH3O2 <=> C5H91-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3O2 <=> C5H91-5 + CH3O2H""",
)

entry(
    index = 1388,
    label = "C5H10-1 + CH3O <=> C5H91-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3O <=> C5H91-3 + CH3OH""",
)

entry(
    index = 1389,
    label = "C5H10-1 + CH3O <=> C5H91-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3O <=> C5H91-4 + CH3OH""",
)

entry(
    index = 1390,
    label = "C5H10-1 + CH3O <=> C5H91-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-1 + CH3O <=> C5H91-5 + CH3OH""",
)

entry(
    index = 1391,
    label = "C5H10-2 + H <=> C5H91-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(173000, 'cm^3/(mol*s)'), n=2.5, Ea=(2492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + H <=> C5H91-3 + H2""",
)

entry(
    index = 1392,
    label = "C5H10-2 + H <=> C5H92-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + H <=> C5H92-4 + H2""",
)

entry(
    index = 1393,
    label = "C5H10-2 + H <=> C5H92-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665100, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + H <=> C5H92-5 + H2""",
)

entry(
    index = 1394,
    label = "C5H10-2 + O <=> C5H91-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(441000, 'cm^3/(mol*s)'), n=2.42, Ea=(3150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + O <=> C5H91-3 + OH""",
)

entry(
    index = 1395,
    label = "C5H10-2 + O <=> C5H92-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(990000, 'cm^3/(mol*s)'), n=2.43, Ea=(1210, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + O <=> C5H92-4 + OH""",
)

entry(
    index = 1396,
    label = "C5H10-2 + O <=> C5H92-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(980000, 'cm^3/(mol*s)'), n=2.43, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + O <=> C5H92-5 + OH""",
)

entry(
    index = 1397,
    label = "C5H10-2 + OH <=> C5H91-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + OH <=> C5H91-3 + H2O""",
)

entry(
    index = 1398,
    label = "C5H10-2 + OH <=> C5H92-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + OH <=> C5H92-4 + H2O""",
)

entry(
    index = 1399,
    label = "C5H10-2 + OH <=> C5H92-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + OH <=> C5H92-5 + H2O""",
)

entry(
    index = 1400,
    label = "C5H10-2 + CH3 <=> C5H91-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.21, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3 <=> C5H91-3 + CH4""",
)

entry(
    index = 1401,
    label = "C5H10-2 + CH3 <=> C5H92-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3 <=> C5H92-4 + CH4""",
)

entry(
    index = 1402,
    label = "C5H10-2 + CH3 <=> C5H92-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4521, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3 <=> C5H92-5 + CH4""",
)

entry(
    index = 1403,
    label = "C5H10-2 + O2 <=> C5H91-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.3e+12, 'cm^3/(mol*s)'), n=0, Ea=(39900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + O2 <=> C5H91-3 + HO2""",
)

entry(
    index = 1404,
    label = "C5H10-2 + O2 <=> C5H92-4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(37220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + O2 <=> C5H92-4 + HO2""",
)

entry(
    index = 1405,
    label = "C5H10-2 + O2 <=> C5H92-5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + O2 <=> C5H92-5 + HO2""",
)

entry(
    index = 1406,
    label = "C5H10-2 + HO2 <=> C5H91-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9639, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + HO2 <=> C5H91-3 + H2O2""",
)

entry(
    index = 1407,
    label = "C5H10-2 + HO2 <=> C5H92-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + HO2 <=> C5H92-4 + H2O2""",
)

entry(
    index = 1408,
    label = "C5H10-2 + HO2 <=> C5H92-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + HO2 <=> C5H92-5 + H2O2""",
)

entry(
    index = 1409,
    label = "C5H10-2 + CH3O2 <=> C5H91-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9639, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3O2 <=> C5H91-3 + CH3O2H""",
)

entry(
    index = 1410,
    label = "C5H10-2 + CH3O2 <=> C5H92-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3O2 <=> C5H92-4 + CH3O2H""",
)

entry(
    index = 1411,
    label = "C5H10-2 + CH3O2 <=> C5H92-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3O2 <=> C5H92-5 + CH3O2H""",
)

entry(
    index = 1412,
    label = "C5H10-2 + CH3O <=> C5H91-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(90, 'cm^3/(mol*s)'), n=2.95, Ea=(11990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3O <=> C5H91-3 + CH3OH""",
)

entry(
    index = 1413,
    label = "C5H10-2 + CH3O <=> C5H92-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3O <=> C5H92-4 + CH3OH""",
)

entry(
    index = 1414,
    label = "C5H10-2 + CH3O <=> C5H92-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 + CH3O <=> C5H92-5 + CH3OH""",
)

entry(
    index = 1415,
    label = "C5H91-3 + HO2 <=> C5H9O1-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-3 + HO2 <=> C5H9O1-3 + OH""",
)

entry(
    index = 1416,
    label = "C5H91-3 + CH3O2 <=> C5H9O1-3 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-3 + CH3O2 <=> C5H9O1-3 + CH3O""",
)

entry(
    index = 1417,
    label = "C5H91-3 + C2H5O2 <=> C5H9O1-3 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-3 + C2H5O2 <=> C5H9O1-3 + C2H5O""",
)

entry(
    index = 1418,
    label = "C5H92-4 + HO2 <=> C5H9O2-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H92-4 + HO2 <=> C5H9O2-4 + OH""",
)

entry(
    index = 1419,
    label = "C5H92-4 + CH3O2 <=> C5H9O2-4 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H92-4 + CH3O2 <=> C5H9O2-4 + CH3O""",
)

entry(
    index = 1420,
    label = "C5H92-4 + C2H5O2 <=> C5H9O2-4 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H92-4 + C2H5O2 <=> C5H9O2-4 + C2H5O""",
)

entry(
    index = 1421,
    label = "C5H91-3 <=> C4H6 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.554e+14, 's^-1'), n=-0.52, Ea=(38520, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-3 <=> C4H6 + CH3""",
)

entry(
    index = 1422,
    label = "C5H91-3 <=> C5H81-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.733e+11, 's^-1'), n=0.636, Ea=(42640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-3 <=> C5H81-3 + H""",
)

entry(
    index = 1423,
    label = "C5H91-4 <=> C3H6 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.814e+11, 's^-1'), n=0.17, Ea=(35850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-4 <=> C3H6 + C2H3""",
)

entry(
    index = 1424,
    label = "C5H91-5 <=> C2H4 + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.157e+16, 's^-1'), n=-1.42, Ea=(17750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-5 <=> C2H4 + C3H5-A""",
)

entry(
    index = 1425,
    label = "C5H92-4 <=> C5H81-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.015e+15, 's^-1'), n=-0.34, Ea=(46030, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H92-4 <=> C5H81-3 + H""",
)

entry(
    index = 1426,
    label = "C5H92-5 <=> C2H4 + C3H5-S",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.899e+16, 's^-1'), n=-1.18, Ea=(42180, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H92-5 <=> C2H4 + C3H5-S""",
)

entry(
    index = 1427,
    label = "C5H9O1-3 <=> C2H3CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.131e+19, 's^-1'), n=-1.85, Ea=(10670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O1-3 <=> C2H3CHO + C2H5""",
)

entry(
    index = 1428,
    label = "C5H9O1-3 <=> C2H5CHO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.417e+18, 's^-1'), n=-1.56, Ea=(23340, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O1-3 <=> C2H5CHO + C2H3""",
)

entry(
    index = 1429,
    label = "C5H81-3 + OH <=> CH2O + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H81-3 + OH <=> CH2O + C4H71-3""",
)

entry(
    index = 1430,
    label = "C5H81-3 + OH <=> C2H3CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H81-3 + OH <=> C2H3CHO + C2H5""",
)

entry(
    index = 1431,
    label = "C5H81-3 + OH <=> CH3CHO + C3H5-S",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H81-3 + OH <=> CH3CHO + C3H5-S""",
)

entry(
    index = 1432,
    label = "C5H9O2-4 <=> SC3H5CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.983e+15, 's^-1'), n=-1.13, Ea=(9941, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O2-4 <=> SC3H5CHO + CH3""",
)

entry(
    index = 1433,
    label = "C5H9O2-4 <=> CH3CHO + C3H5-S",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.073e+22, 's^-1'), n=-2.66, Ea=(29650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O2-4 <=> CH3CHO + C3H5-S""",
)

entry(
    index = 1434,
    label = "C5H10-2 <=> CH3 + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.486e+19, 's^-1'), n=-1.367, Ea=(76320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10-2 <=> CH3 + C4H71-3""",
)

entry(
    index = 1435,
    label = "C5H11-1 + O2 <=> C5H10-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.837, 'cm^3/(mol*s)'), n=3.59, Ea=(11960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 + O2 <=> C5H10-1 + HO2""",
)

entry(
    index = 1436,
    label = "C5H11-2 + O2 <=> C5H10-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.535, 'cm^3/(mol*s)'), n=3.71, Ea=(9322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + O2 <=> C5H10-1 + HO2""",
)

entry(
    index = 1437,
    label = "C5H11-2 + O2 <=> C5H10-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.07, 'cm^3/(mol*s)'), n=3.71, Ea=(9322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + O2 <=> C5H10-2 + HO2""",
)

entry(
    index = 1438,
    label = "C5H11-3 + O2 <=> C5H10-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.14, 'cm^3/(mol*s)'), n=3.71, Ea=(9322, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 + O2 <=> C5H10-2 + HO2""",
)

entry(
    index = 1439,
    label = "C5H11-1 + HO2 <=> C5H11O-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 + HO2 <=> C5H11O-1 + OH""",
)

entry(
    index = 1440,
    label = "C5H11-2 + HO2 <=> C5H11O-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + HO2 <=> C5H11O-2 + OH""",
)

entry(
    index = 1441,
    label = "C5H11-3 + HO2 <=> C5H11O-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 + HO2 <=> C5H11O-3 + OH""",
)

entry(
    index = 1442,
    label = "C5H11-1 + CH3O2 <=> C5H11O-1 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 + CH3O2 <=> C5H11O-1 + CH3O""",
)

entry(
    index = 1443,
    label = "C5H11-2 + CH3O2 <=> C5H11O-2 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + CH3O2 <=> C5H11O-2 + CH3O""",
)

entry(
    index = 1444,
    label = "C5H11-3 + CH3O2 <=> C5H11O-3 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 + CH3O2 <=> C5H11O-3 + CH3O""",
)

entry(
    index = 1445,
    label = "C5H11O2-1 + NC5H12 <=> C5H11O2H-1 + C5H11-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + NC5H12 <=> C5H11O2H-1 + C5H11-1""",
)

entry(
    index = 1446,
    label = "C5H11O2-2 + NC5H12 <=> C5H11O2H-2 + C5H11-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + NC5H12 <=> C5H11O2H-2 + C5H11-1""",
)

entry(
    index = 1447,
    label = "C5H11O2-3 + NC5H12 <=> C5H11O2H-3 + C5H11-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 + NC5H12 <=> C5H11O2H-3 + C5H11-1""",
)

entry(
    index = 1448,
    label = "C5H11O2-1 + NC5H12 <=> C5H11O2H-1 + C5H11-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + NC5H12 <=> C5H11O2H-1 + C5H11-2""",
)

entry(
    index = 1449,
    label = "C5H11O2-2 + NC5H12 <=> C5H11O2H-2 + C5H11-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + NC5H12 <=> C5H11O2H-2 + C5H11-2""",
)

entry(
    index = 1450,
    label = "C5H11O2-3 + NC5H12 <=> C5H11O2H-3 + C5H11-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 + NC5H12 <=> C5H11O2H-3 + C5H11-2""",
)

entry(
    index = 1451,
    label = "C5H11O2-1 + NC5H12 <=> C5H11O2H-1 + C5H11-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + NC5H12 <=> C5H11O2H-1 + C5H11-3""",
)

entry(
    index = 1452,
    label = "C5H11O2-2 + NC5H12 <=> C5H11O2H-2 + C5H11-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + NC5H12 <=> C5H11O2H-2 + C5H11-3""",
)

entry(
    index = 1453,
    label = "C5H11O2-3 + NC5H12 <=> C5H11O2H-3 + C5H11-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 + NC5H12 <=> C5H11O2H-3 + C5H11-3""",
)

entry(
    index = 1454,
    label = "C5H11-1 + C5H11O2-1 <=> C5H11O-1 + C5H11O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 + C5H11O2-1 <=> C5H11O-1 + C5H11O-1""",
)

entry(
    index = 1455,
    label = "C5H11-1 + C5H11O2-2 <=> C5H11O-1 + C5H11O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 + C5H11O2-2 <=> C5H11O-1 + C5H11O-2""",
)

entry(
    index = 1456,
    label = "C5H11-1 + C5H11O2-3 <=> C5H11O-1 + C5H11O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-1 + C5H11O2-3 <=> C5H11O-1 + C5H11O-3""",
)

entry(
    index = 1457,
    label = "C5H11-2 + C5H11O2-1 <=> C5H11O-2 + C5H11O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + C5H11O2-1 <=> C5H11O-2 + C5H11O-1""",
)

entry(
    index = 1458,
    label = "C5H11-2 + C5H11O2-2 <=> C5H11O-2 + C5H11O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + C5H11O2-2 <=> C5H11O-2 + C5H11O-2""",
)

entry(
    index = 1459,
    label = "C5H11-2 + C5H11O2-3 <=> C5H11O-2 + C5H11O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-2 + C5H11O2-3 <=> C5H11O-2 + C5H11O-3""",
)

entry(
    index = 1460,
    label = "C5H11-3 + C5H11O2-1 <=> C5H11O-3 + C5H11O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 + C5H11O2-1 <=> C5H11O-3 + C5H11O-1""",
)

entry(
    index = 1461,
    label = "C5H11-3 + C5H11O2-2 <=> C5H11O-3 + C5H11O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 + C5H11O2-2 <=> C5H11O-3 + C5H11O-2""",
)

entry(
    index = 1462,
    label = "C5H11-3 + C5H11O2-3 <=> C5H11O-3 + C5H11O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11-3 + C5H11O2-3 <=> C5H11O-3 + C5H11O-3""",
)

entry(
    index = 1463,
    label = "C5H11O2-1 + C5H11O2-2 => O2 + C5H11O-1 + C5H11O-2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + C5H11O2-2 => O2 + C5H11O-1 + C5H11O-2""",
)

entry(
    index = 1464,
    label = "C5H11O2-1 + C5H11O2-3 => O2 + C5H11O-1 + C5H11O-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + C5H11O2-3 => O2 + C5H11O-1 + C5H11O-3""",
)

entry(
    index = 1465,
    label = "C5H11O2-1 + CH3O2 => O2 + C5H11O-1 + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + CH3O2 => O2 + C5H11O-1 + CH3O""",
)

entry(
    index = 1466,
    label = "C5H11O2-1 + C5H11O2-1 => O2 + C5H11O-1 + C5H11O-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + C5H11O2-1 => O2 + C5H11O-1 + C5H11O-1""",
)

entry(
    index = 1467,
    label = "H2O2 + C5H11O2-1 <=> HO2 + C5H11O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C5H11O2-1 <=> HO2 + C5H11O2H-1""",
)

entry(
    index = 1468,
    label = "C5H11O2-1 + HO2 <=> C5H11O2H-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 + HO2 <=> C5H11O2H-1 + O2""",
)

entry(
    index = 1469,
    label = "C5H11O2-2 + CH3O2 => O2 + C5H11O-2 + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + CH3O2 => O2 + C5H11O-2 + CH3O""",
)

entry(
    index = 1470,
    label = "C5H11O2-2 + C5H11O2-3 => C5H11O-2 + C5H11O-3 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + C5H11O2-3 => C5H11O-2 + C5H11O-3 + O2""",
)

entry(
    index = 1471,
    label = "C5H11O2-2 + C5H11O2-2 => O2 + C5H11O-2 + C5H11O-2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + C5H11O2-2 => O2 + C5H11O-2 + C5H11O-2""",
)

entry(
    index = 1472,
    label = "H2O2 + C5H11O2-2 <=> HO2 + C5H11O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C5H11O2-2 <=> HO2 + C5H11O2H-2""",
)

entry(
    index = 1473,
    label = "C5H11O2-2 + HO2 <=> C5H11O2H-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 + HO2 <=> C5H11O2H-2 + O2""",
)

entry(
    index = 1474,
    label = "C5H11O2-3 + CH3O2 => O2 + C5H11O-3 + CH3O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 + CH3O2 => O2 + C5H11O-3 + CH3O""",
)

entry(
    index = 1475,
    label = "C5H11O2-3 + C5H11O2-3 => O2 + C5H11O-3 + C5H11O-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 + C5H11O2-3 => O2 + C5H11O-3 + C5H11O-3""",
)

entry(
    index = 1476,
    label = "H2O2 + C5H11O2-3 <=> HO2 + C5H11O2H-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C5H11O2-3 <=> HO2 + C5H11O2H-3""",
)

entry(
    index = 1477,
    label = "C5H11O2-3 + HO2 <=> C5H11O2H-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 + HO2 <=> C5H11O2H-3 + O2""",
)

entry(
    index = 1478,
    label = "C5H11O2H-1 <=> C5H11O-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2H-1 <=> C5H11O-1 + OH""",
)

entry(
    index = 1479,
    label = "C5H11O2H-2 <=> C5H11O-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.45e+15, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2H-2 <=> C5H11O-2 + OH""",
)

entry(
    index = 1480,
    label = "C5H11O2H-3 <=> C5H11O-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.45e+15, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2H-3 <=> C5H11O-3 + OH""",
)

entry(
    index = 1481,
    label = "C5H11O-1 <=> CH2O + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.661e+20, 's^-1'), n=-2.247, Ea=(24960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O-1 <=> CH2O + PC4H9""",
)

entry(
    index = 1482,
    label = "C5H11O-2 <=> CH3CHO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.689e+22, 's^-1'), n=-2.601, Ea=(19550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O-2 <=> CH3CHO + NC3H7""",
)

entry(
    index = 1483,
    label = "C5H11O-3 <=> C2H5 + C2H5CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.239e+18, 's^-1'), n=-1.199, Ea=(18590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O-3 <=> C2H5 + C2H5CHO""",
)

entry(
    index = 1484,
    label = "C5H11O2-1 <=> C5H11-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.338e+20, 's^-1'), n=-1.62, Ea=(35830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 <=> C5H11-1 + O2""",
)

entry(
    index = 1485,
    label = "C5H11O2-2 <=> C5H11-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.087e+22, 's^-1'), n=-2.287, Ea=(38150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H11-2 + O2""",
)

entry(
    index = 1486,
    label = "C5H11O2-3 <=> C5H11-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.415e+22, 's^-1'), n=-2.282, Ea=(38150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 <=> C5H11-3 + O2""",
)

entry(
    index = 1487,
    label = "C5H11O2-1 <=> C5H10OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 <=> C5H10OOH1-2""",
)

entry(
    index = 1488,
    label = "C5H11O2-1 <=> C5H10OOH1-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 <=> C5H10OOH1-3""",
)

entry(
    index = 1489,
    label = "C5H11O2-1 <=> C5H10OOH1-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(19050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 <=> C5H10OOH1-4""",
)

entry(
    index = 1490,
    label = "C5H11O2-1 <=> C5H10OOH1-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.86e+08, 's^-1'), n=0, Ea=(25550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 <=> C5H10OOH1-5""",
)

entry(
    index = 1491,
    label = "C5H11O2-2 <=> C5H10OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H10OOH2-1""",
)

entry(
    index = 1492,
    label = "C5H11O2-2 <=> C5H10OOH2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H10OOH2-3""",
)

entry(
    index = 1493,
    label = "C5H11O2-2 <=> C5H10OOH2-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H10OOH2-4""",
)

entry(
    index = 1494,
    label = "C5H11O2-2 <=> C5H10OOH2-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.688e+09, 's^-1'), n=0, Ea=(22350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H10OOH2-5""",
)

entry(
    index = 1495,
    label = "C5H11O2-3 <=> C5H10OOH3-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 <=> C5H10OOH3-1""",
)

entry(
    index = 1496,
    label = "C5H11O2-3 <=> C5H10OOH3-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 <=> C5H10OOH3-2""",
)

entry(
    index = 1497,
    label = "C5H11O2-1 <=> C5H10-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-1 <=> C5H10-1 + HO2""",
)

entry(
    index = 1498,
    label = "C5H11O2-2 <=> C5H10-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.075e+42, 's^-1'), n=-9.41, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H10-1 + HO2""",
)

entry(
    index = 1499,
    label = "C5H11O2-2 <=> C5H10-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-2 <=> C5H10-2 + HO2""",
)

entry(
    index = 1500,
    label = "C5H11O2-3 <=> C5H10-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.009e+39, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H11O2-3 <=> C5H10-2 + HO2""",
)

entry(
    index = 1501,
    label = "C5H10OOH1-2 => C5H10O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-2 => C5H10O1-2 + OH""",
)

entry(
    index = 1502,
    label = "C5H10OOH1-3 => C5H10O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-3 => C5H10O1-3 + OH""",
)

entry(
    index = 1503,
    label = "C5H10OOH1-4 => C5H10O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-4 => C5H10O1-4 + OH""",
)

entry(
    index = 1504,
    label = "C5H10OOH1-5 => C5H10O1-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.172e+09, 's^-1'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-5 => C5H10O1-5 + OH""",
)

entry(
    index = 1505,
    label = "C5H10OOH2-1 => C5H10O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-1 => C5H10O1-2 + OH""",
)

entry(
    index = 1506,
    label = "C5H10OOH2-3 => C5H10O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-3 => C5H10O2-3 + OH""",
)

entry(
    index = 1507,
    label = "C5H10OOH2-4 => C5H10O2-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-4 => C5H10O2-4 + OH""",
)

entry(
    index = 1508,
    label = "C5H10OOH2-5 => C5H10O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(6000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-5 => C5H10O1-4 + OH""",
)

entry(
    index = 1509,
    label = "C5H10OOH3-2 => C5H10O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-2 => C5H10O2-3 + OH""",
)

entry(
    index = 1510,
    label = "C5H10OOH3-1 => C5H10O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-1 => C5H10O1-3 + OH""",
)

entry(
    index = 1511,
    label = "C5H10O1-2 + OH => CH2CO + NC3H7 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-2 + OH => CH2CO + NC3H7 + H2O""",
)

entry(
    index = 1512,
    label = "C5H10O1-3 + OH => C2H4 + C2H5CO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-3 + OH => C2H4 + C2H5CO + H2O""",
)

entry(
    index = 1513,
    label = "C5H10O1-4 + OH => CH3COCH2 + C2H4 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-4 + OH => CH3COCH2 + C2H4 + H2O""",
)

entry(
    index = 1514,
    label = "C5H10O1-5 + OH => CH2CH2CHO + C2H4 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-5 + OH => CH2CH2CHO + C2H4 + H2O""",
)

entry(
    index = 1515,
    label = "C5H10O2-3 + OH => CH3CHCO + C2H5 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-3 + OH => CH3CHCO + C2H5 + H2O""",
)

entry(
    index = 1516,
    label = "C5H10O2-4 + OH => CH3CO + C3H6 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-4 + OH => CH3CO + C3H6 + H2O""",
)

entry(
    index = 1517,
    label = "C5H10O1-2 + OH => C2H3CHO + C2H5 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-2 + OH => C2H3CHO + C2H5 + H2O""",
)

entry(
    index = 1518,
    label = "C5H10O1-3 + OH => HCO + C4H8-1 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-3 + OH => HCO + C4H8-1 + H2O""",
)

entry(
    index = 1519,
    label = "C5H10O1-4 + OH => CH2CHO + C3H6 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-4 + OH => CH2CHO + C3H6 + H2O""",
)

entry(
    index = 1520,
    label = "C5H10O1-5 + OH => CH2O + C4H71-3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-5 + OH => CH2O + C4H71-3 + H2O""",
)

entry(
    index = 1521,
    label = "C5H10O2-3 + OH => C2H3COCH3 + CH3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-3 + OH => C2H3COCH3 + CH3 + H2O""",
)

entry(
    index = 1522,
    label = "C5H10O2-4 + OH => CH3CHO + C3H5-S + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-4 + OH => CH3CHO + C3H5-S + H2O""",
)

entry(
    index = 1523,
    label = "C5H10O1-2 + HO2 => CH2CO + NC3H7 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-2 + HO2 => CH2CO + NC3H7 + H2O2""",
)

entry(
    index = 1524,
    label = "C5H10O1-3 + HO2 => C2H4 + C2H5CO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-3 + HO2 => C2H4 + C2H5CO + H2O2""",
)

entry(
    index = 1525,
    label = "C5H10O1-4 + HO2 => CH3COCH2 + C2H4 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-4 + HO2 => CH3COCH2 + C2H4 + H2O2""",
)

entry(
    index = 1526,
    label = "C5H10O1-5 + HO2 => CH2CH2CHO + C2H4 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-5 + HO2 => CH2CH2CHO + C2H4 + H2O2""",
)

entry(
    index = 1527,
    label = "C5H10O2-3 + HO2 => CH3CHCO + C2H5 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-3 + HO2 => CH3CHCO + C2H5 + H2O2""",
)

entry(
    index = 1528,
    label = "C5H10O2-4 + HO2 => CH3CO + C3H6 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-4 + HO2 => CH3CO + C3H6 + H2O2""",
)

entry(
    index = 1529,
    label = "C5H10O1-2 + HO2 => C2H3CHO + C2H5 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-2 + HO2 => C2H3CHO + C2H5 + H2O2""",
)

entry(
    index = 1530,
    label = "C5H10O1-3 + HO2 => HCO + C4H8-1 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-3 + HO2 => HCO + C4H8-1 + H2O2""",
)

entry(
    index = 1531,
    label = "C5H10O1-4 + HO2 => CH2CHO + C3H6 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-4 + HO2 => CH2CHO + C3H6 + H2O2""",
)

entry(
    index = 1532,
    label = "C5H10O1-5 + HO2 => CH2O + C4H71-3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O1-5 + HO2 => CH2O + C4H71-3 + H2O2""",
)

entry(
    index = 1533,
    label = "C5H10O2-3 + HO2 => C2H3COCH3 + CH3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-3 + HO2 => C2H3COCH3 + CH3 + H2O2""",
)

entry(
    index = 1534,
    label = "C5H10O2-4 + HO2 => CH3CHO + C3H5-S + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10O2-4 + HO2 => CH3CHO + C3H5-S + H2O2""",
)

entry(
    index = 1535,
    label = "C5H10OOH1-2 <=> C5H10-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.449e+17, 's^-1'), n=-1.556, Ea=(17980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-2 <=> C5H10-1 + HO2""",
)

entry(
    index = 1536,
    label = "C5H10OOH2-1 <=> C5H10-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.462e+19, 's^-1'), n=-2.231, Ea=(21050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-1 <=> C5H10-1 + HO2""",
)

entry(
    index = 1537,
    label = "C5H10OOH2-3 <=> C5H10-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.651e+19, 's^-1'), n=-2.455, Ea=(20680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-3 <=> C5H10-2 + HO2""",
)

entry(
    index = 1538,
    label = "C5H10OOH3-2 <=> C5H10-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.651e+19, 's^-1'), n=-2.455, Ea=(20680, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-2 <=> C5H10-2 + HO2""",
)

entry(
    index = 1539,
    label = "C5H10OOH1-3 => OH + CH2O + C4H8-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.283e+13, 's^-1'), n=-0.17, Ea=(30090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-3 => OH + CH2O + C4H8-1""",
)

entry(
    index = 1540,
    label = "C5H10OOH2-4 => OH + CH3CHO + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.312e+17, 's^-1'), n=-1.4, Ea=(27170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-4 => OH + CH3CHO + C3H6""",
)

entry(
    index = 1541,
    label = "C5H10OOH3-1 => OH + C2H5CHO + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.096e+18, 's^-1'), n=-1.73, Ea=(26820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-1 => OH + C2H5CHO + C2H4""",
)

entry(
    index = 1542,
    label = "C5H10OOH1-2O2 <=> C5H10OOH1-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.039e+22, 's^-1'), n=-2.295, Ea=(37970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-2O2 <=> C5H10OOH1-2 + O2""",
)

entry(
    index = 1543,
    label = "C5H10OOH1-3O2 <=> C5H10OOH1-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.039e+22, 's^-1'), n=-2.295, Ea=(37970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-3O2 <=> C5H10OOH1-3 + O2""",
)

entry(
    index = 1544,
    label = "C5H10OOH1-4O2 <=> C5H10OOH1-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.039e+22, 's^-1'), n=-2.295, Ea=(37970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-4O2 <=> C5H10OOH1-4 + O2""",
)

entry(
    index = 1545,
    label = "C5H10OOH1-5O2 <=> C5H10OOH1-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.777e+20, 's^-1'), n=-1.623, Ea=(35690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-5O2 <=> C5H10OOH1-5 + O2""",
)

entry(
    index = 1546,
    label = "C5H10OOH2-1O2 <=> C5H10OOH2-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.572e+20, 's^-1'), n=-1.62, Ea=(35650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-1O2 <=> C5H10OOH2-1 + O2""",
)

entry(
    index = 1547,
    label = "C5H10OOH2-3O2 <=> C5H10OOH2-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.251e+22, 's^-1'), n=-2.29, Ea=(37910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-3O2 <=> C5H10OOH2-3 + O2""",
)

entry(
    index = 1548,
    label = "C5H10OOH2-4O2 <=> C5H10OOH2-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.251e+22, 's^-1'), n=-2.29, Ea=(37910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-4O2 <=> C5H10OOH2-4 + O2""",
)

entry(
    index = 1549,
    label = "C5H10OOH2-5O2 <=> C5H10OOH2-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.572e+20, 's^-1'), n=-1.62, Ea=(35650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-5O2 <=> C5H10OOH2-5 + O2""",
)

entry(
    index = 1550,
    label = "C5H10OOH3-1O2 <=> C5H10OOH3-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.572e+20, 's^-1'), n=-1.62, Ea=(35650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-1O2 <=> C5H10OOH3-1 + O2""",
)

entry(
    index = 1551,
    label = "C5H10OOH3-2O2 <=> C5H10OOH3-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.251e+22, 's^-1'), n=-2.29, Ea=(37910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-2O2 <=> C5H10OOH3-2 + O2""",
)

entry(
    index = 1552,
    label = "C5H10OOH1-2O2 <=> NC5KET12 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-2O2 <=> NC5KET12 + OH""",
)

entry(
    index = 1553,
    label = "C5H10OOH1-3O2 <=> NC5KET13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(21400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-3O2 <=> NC5KET13 + OH""",
)

entry(
    index = 1554,
    label = "C5H10OOH1-4O2 <=> NC5KET14 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(19350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-4O2 <=> NC5KET14 + OH""",
)

entry(
    index = 1555,
    label = "C5H10OOH1-5O2 <=> NC5KET15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.907e+08, 's^-1'), n=0, Ea=(22550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH1-5O2 <=> NC5KET15 + OH""",
)

entry(
    index = 1556,
    label = "C5H10OOH2-1O2 <=> NC5KET21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-1O2 <=> NC5KET21 + OH""",
)

entry(
    index = 1557,
    label = "C5H10OOH2-3O2 <=> NC5KET23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-3O2 <=> NC5KET23 + OH""",
)

entry(
    index = 1558,
    label = "C5H10OOH2-4O2 <=> NC5KET24 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-4O2 <=> NC5KET24 + OH""",
)

entry(
    index = 1559,
    label = "C5H10OOH2-5O2 <=> NC5KET25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(16050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH2-5O2 <=> NC5KET25 + OH""",
)

entry(
    index = 1560,
    label = "C5H10OOH3-1O2 <=> NC5KET31 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-1O2 <=> NC5KET31 + OH""",
)

entry(
    index = 1561,
    label = "C5H10OOH3-2O2 <=> NC5KET32 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OOH3-2O2 <=> NC5KET32 + OH""",
)

entry(
    index = 1562,
    label = "NC5KET12 => NC3H7CHO + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET12 => NC3H7CHO + HCO + OH""",
)

entry(
    index = 1563,
    label = "NC5KET13 => C2H5CHO + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET13 => C2H5CHO + CH2CHO + OH""",
)

entry(
    index = 1564,
    label = "NC5KET14 => CH3CHO + CH2CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET14 => CH3CHO + CH2CH2CHO + OH""",
)

entry(
    index = 1565,
    label = "NC5KET15 => CH2O + C3H6CHO-1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET15 => CH2O + C3H6CHO-1 + OH""",
)

entry(
    index = 1566,
    label = "NC5KET21 => CH2O + NC3H7CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET21 => CH2O + NC3H7CO + OH""",
)

entry(
    index = 1567,
    label = "NC5KET23 => C2H5CHO + CH3CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET23 => C2H5CHO + CH3CO + OH""",
)

entry(
    index = 1568,
    label = "NC5KET24 => CH3CHO + CH3COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET24 => CH3CHO + CH3COCH2 + OH""",
)

entry(
    index = 1569,
    label = "NC5KET25 => CH2O + CH2CH2COCH3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET25 => CH2O + CH2CH2COCH3 + OH""",
)

entry(
    index = 1570,
    label = "NC5KET31 => CH2O + C2H5COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+16, 's^-1'), n=0, Ea=(42000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET31 => CH2O + C2H5COCH2 + OH""",
)

entry(
    index = 1571,
    label = "NC5KET32 => CH3CHO + C2H5CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.05e+16, 's^-1'), n=0, Ea=(41600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5KET32 => CH3CHO + C2H5CO + OH""",
)

entry(
    index = 1572,
    label = "C5H10OH-1 <=> C5H10-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(25830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OH-1 <=> C5H10-1 + OH""",
)

entry(
    index = 1573,
    label = "C5H10OH-2 <=> C5H10-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(25830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10OH-2 <=> C5H10-2 + OH""",
)

entry(
    index = 1574,
    label = "O2C5H10OH-1 <=> C5H10OH-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.058e+21, 's^-1'), n=-1.84, Ea=(37750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C5H10OH-1 <=> C5H10OH-1 + O2""",
)

entry(
    index = 1575,
    label = "O2C5H10OH-1 => NC3H7CHO + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C5H10OH-1 => NC3H7CHO + CH2O + OH""",
)

entry(
    index = 1576,
    label = "O2C5H10OH-2 <=> C5H10OH-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.385e+21, 's^-1'), n=-2.01, Ea=(37870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C5H10OH-2 <=> C5H10OH-2 + O2""",
)

entry(
    index = 1577,
    label = "O2C5H10OH-2 => C2H5CHO + CH3CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C5H10OH-2 => C2H5CHO + CH3CHO + OH""",
)

entry(
    index = 1578,
    label = "NC6H14 <=> C5H11-1 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.608e+22, 's^-1'), n=-1.61, Ea=(89350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 <=> C5H11-1 + CH3""",
)

entry(
    index = 1579,
    label = "NC6H14 <=> NC3H7 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.169e+24, 's^-1'), n=-2.19, Ea=(87840, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 <=> NC3H7 + NC3H7""",
)

entry(
    index = 1580,
    label = "NC6H14 <=> PC4H9 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.617e+24, 's^-1'), n=-2.21, Ea=(88580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 <=> PC4H9 + C2H5""",
)

entry(
    index = 1581,
    label = "NC6H14 <=> C6H13-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.356e+17, 's^-1'), n=-0.36, Ea=(101200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 <=> C6H13-1 + H""",
)

entry(
    index = 1582,
    label = "NC6H14 <=> C6H13-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.088e+18, 's^-1'), n=-0.7, Ea=(98710, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 <=> C6H13-2 + H""",
)

entry(
    index = 1583,
    label = "NC6H14 <=> C6H13-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.088e+18, 's^-1'), n=-0.7, Ea=(98710, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 <=> C6H13-3 + H""",
)

entry(
    index = 1584,
    label = "NC6H14 + H <=> C6H13-1 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(188000, 'cm^3/(mol*s)'), n=2.75, Ea=(6280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + H <=> C6H13-1 + H2""",
)

entry(
    index = 1585,
    label = "NC6H14 + H <=> C6H13-2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + H <=> C6H13-2 + H2""",
)

entry(
    index = 1586,
    label = "NC6H14 + H <=> C6H13-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + H <=> C6H13-3 + H2""",
)

entry(
    index = 1587,
    label = "NC6H14 + O <=> C6H13-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.092e+06, 'cm^3/(mol*s)'),
        n = 2.42,
        Ea = (4766, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O <=> C6H13-1 + OH""",
)

entry(
    index = 1588,
    label = "NC6H14 + O <=> C6H13-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.189e+06, 'cm^3/(mol*s)'),
        n = 2.44,
        Ea = (2846, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O <=> C6H13-2 + OH""",
)

entry(
    index = 1589,
    label = "NC6H14 + O <=> C6H13-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.189e+06, 'cm^3/(mol*s)'),
        n = 2.44,
        Ea = (2846, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O <=> C6H13-3 + OH""",
)

entry(
    index = 1590,
    label = "NC6H14 + OH <=> C6H13-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.57e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(954, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + OH <=> C6H13-1 + H2O""",
)

entry(
    index = 1591,
    label = "NC6H14 + OH <=> C6H13-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+06, 'cm^3/(mol*s)'), n=2, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + OH <=> C6H13-2 + H2O""",
)

entry(
    index = 1592,
    label = "NC6H14 + OH <=> C6H13-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+06, 'cm^3/(mol*s)'), n=2, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + OH <=> C6H13-3 + H2O""",
)

entry(
    index = 1593,
    label = "NC6H14 + CH3 <=> C6H13-1 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.904, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3 <=> C6H13-1 + CH4""",
)

entry(
    index = 1594,
    label = "NC6H14 + CH3 <=> C6H13-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(168000, 'cm^3/(mol*s)'), n=2.13, Ea=(7574, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3 <=> C6H13-3 + CH4""",
)

entry(
    index = 1595,
    label = "NC6H14 + CH3 <=> C6H13-2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(168000, 'cm^3/(mol*s)'), n=2.13, Ea=(7574, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3 <=> C6H13-2 + CH4""",
)

entry(
    index = 1596,
    label = "NC6H14 + HO2 <=> C6H13-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(81000, 'cm^3/(mol*s)'), n=2.5, Ea=(16690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + HO2 <=> C6H13-1 + H2O2""",
)

entry(
    index = 1597,
    label = "NC6H14 + HO2 <=> C6H13-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(117600, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + HO2 <=> C6H13-2 + H2O2""",
)

entry(
    index = 1598,
    label = "NC6H14 + HO2 <=> C6H13-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(117600, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + HO2 <=> C6H13-3 + H2O2""",
)

entry(
    index = 1599,
    label = "NC6H14 + CH3O <=> C6H13-1 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.16e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3O <=> C6H13-1 + CH3OH""",
)

entry(
    index = 1600,
    label = "NC6H14 + CH3O <=> C6H13-2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.19e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3O <=> C6H13-2 + CH3OH""",
)

entry(
    index = 1601,
    label = "NC6H14 + CH3O <=> C6H13-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.19e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3O <=> C6H13-3 + CH3OH""",
)

entry(
    index = 1602,
    label = "NC6H14 + O2 <=> C6H13-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(52800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2 <=> C6H13-1 + HO2""",
)

entry(
    index = 1603,
    label = "NC6H14 + O2 <=> C6H13-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(50160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2 <=> C6H13-2 + HO2""",
)

entry(
    index = 1604,
    label = "NC6H14 + O2 <=> C6H13-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(50160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2 <=> C6H13-3 + HO2""",
)

entry(
    index = 1605,
    label = "NC6H14 + C2H5 <=> C6H13-1 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(13400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C2H5 <=> C6H13-1 + C2H6""",
)

entry(
    index = 1606,
    label = "NC6H14 + C2H5 <=> C6H13-2 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C2H5 <=> C6H13-2 + C2H6""",
)

entry(
    index = 1607,
    label = "NC6H14 + C2H5 <=> C6H13-3 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C2H5 <=> C6H13-3 + C2H6""",
)

entry(
    index = 1608,
    label = "NC6H14 + C2H3 <=> C6H13-1 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(18000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C2H3 <=> C6H13-1 + C2H4""",
)

entry(
    index = 1609,
    label = "NC6H14 + C2H3 <=> C6H13-2 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C2H3 <=> C6H13-2 + C2H4""",
)

entry(
    index = 1610,
    label = "NC6H14 + C2H3 <=> C6H13-3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C2H3 <=> C6H13-3 + C2H4""",
)

entry(
    index = 1611,
    label = "NC6H14 + CH3O2 <=> C6H13-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(81000, 'cm^3/(mol*s)'), n=2.5, Ea=(16690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3O2 <=> C6H13-1 + CH3O2H""",
)

entry(
    index = 1612,
    label = "NC6H14 + CH3O2 <=> C6H13-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(117600, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3O2 <=> C6H13-2 + CH3O2H""",
)

entry(
    index = 1613,
    label = "NC6H14 + CH3O2 <=> C6H13-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(117600, 'cm^3/(mol*s)'), n=2.5, Ea=(14860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + CH3O2 <=> C6H13-3 + CH3O2H""",
)

entry(
    index = 1614,
    label = "NC6H14 + O2CHO <=> C6H13-1 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.68e+13, 'cm^3/(mol*s)'), n=0, Ea=(20440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2CHO <=> C6H13-1 + HO2CHO""",
)

entry(
    index = 1615,
    label = "NC6H14 + O2CHO <=> C6H13-2 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2CHO <=> C6H13-2 + HO2CHO""",
)

entry(
    index = 1616,
    label = "NC6H14 + O2CHO <=> C6H13-3 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + O2CHO <=> C6H13-3 + HO2CHO""",
)

entry(
    index = 1617,
    label = "C6H13O2-1 + NC6H14 <=> C6H13O2H-1 + C6H13-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + NC6H14 <=> C6H13O2H-1 + C6H13-1""",
)

entry(
    index = 1618,
    label = "C6H13O2-2 + NC6H14 <=> C6H13O2H-2 + C6H13-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + NC6H14 <=> C6H13O2H-2 + C6H13-1""",
)

entry(
    index = 1619,
    label = "C6H13O2-3 + NC6H14 <=> C6H13O2H-3 + C6H13-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + NC6H14 <=> C6H13O2H-3 + C6H13-1""",
)

entry(
    index = 1620,
    label = "C6H13O2-1 + NC6H14 <=> C6H13O2H-1 + C6H13-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + NC6H14 <=> C6H13O2H-1 + C6H13-2""",
)

entry(
    index = 1621,
    label = "C6H13O2-2 + NC6H14 <=> C6H13O2H-2 + C6H13-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + NC6H14 <=> C6H13O2H-2 + C6H13-2""",
)

entry(
    index = 1622,
    label = "C6H13O2-3 + NC6H14 <=> C6H13O2H-3 + C6H13-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + NC6H14 <=> C6H13O2H-3 + C6H13-2""",
)

entry(
    index = 1623,
    label = "C6H13O2-1 + NC6H14 <=> C6H13O2H-1 + C6H13-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + NC6H14 <=> C6H13O2H-1 + C6H13-3""",
)

entry(
    index = 1624,
    label = "C6H13O2-2 + NC6H14 <=> C6H13O2H-2 + C6H13-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + NC6H14 <=> C6H13O2H-2 + C6H13-3""",
)

entry(
    index = 1625,
    label = "C6H13O2-3 + NC6H14 <=> C6H13O2H-3 + C6H13-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + NC6H14 <=> C6H13O2H-3 + C6H13-3""",
)

entry(
    index = 1626,
    label = "NC6H14 + C6H13-1 <=> NC6H14 + C6H13-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C6H13-1 <=> NC6H14 + C6H13-2""",
)

entry(
    index = 1627,
    label = "NC6H14 + C6H13-1 <=> NC6H14 + C6H13-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C6H13-1 <=> NC6H14 + C6H13-3""",
)

entry(
    index = 1628,
    label = "NC6H14 + C6H13-2 <=> NC6H14 + C6H13-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6H14 + C6H13-2 <=> NC6H14 + C6H13-3""",
)

entry(
    index = 1629,
    label = "C6H13-1 + HO2 <=> C6H13O-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 + HO2 <=> C6H13O-1 + OH""",
)

entry(
    index = 1630,
    label = "C6H13-2 + HO2 <=> C6H13O-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + HO2 <=> C6H13O-2 + OH""",
)

entry(
    index = 1631,
    label = "C6H13-3 + HO2 <=> C6H13O-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + HO2 <=> C6H13O-3 + OH""",
)

entry(
    index = 1632,
    label = "C6H13-1 + CH3O2 <=> C6H13O-1 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 + CH3O2 <=> C6H13O-1 + CH3O""",
)

entry(
    index = 1633,
    label = "C6H13-2 + CH3O2 <=> C6H13O-2 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + CH3O2 <=> C6H13O-2 + CH3O""",
)

entry(
    index = 1634,
    label = "C6H13-3 + CH3O2 <=> C6H13O-3 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + CH3O2 <=> C6H13O-3 + CH3O""",
)

entry(
    index = 1635,
    label = "C6H13-1 + O2 <=> C6H12-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-19, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 + O2 <=> C6H12-1 + HO2""",
)

entry(
    index = 1636,
    label = "C6H13-2 + O2 <=> C6H12-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e-19, 'cm^3/(mol*s)'), n=0, Ea=(5020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + O2 <=> C6H12-1 + HO2""",
)

entry(
    index = 1637,
    label = "C6H13-2 + O2 <=> C6H12-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-19, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + O2 <=> C6H12-2 + HO2""",
)

entry(
    index = 1638,
    label = "C6H13-3 + O2 <=> C6H12-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-19, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + O2 <=> C6H12-2 + HO2""",
)

entry(
    index = 1639,
    label = "C6H13-3 + O2 <=> C6H12-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-19, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + O2 <=> C6H12-3 + HO2""",
)

entry(
    index = 1640,
    label = "C6H13-1 <=> C2H4 + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.391e+19, 's^-1'), n=-1.97, Ea=(30640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 <=> C2H4 + PC4H9""",
)

entry(
    index = 1641,
    label = "C6H13-1 <=> C6H12-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.621e+13, 's^-1'), n=-0.26, Ea=(36000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 <=> C6H12-1 + H""",
)

entry(
    index = 1642,
    label = "C6H13-2 <=> C3H6 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.827e+19, 's^-1'), n=-1.8, Ea=(30170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 <=> C3H6 + NC3H7""",
)

entry(
    index = 1643,
    label = "C6H13-2 <=> C6H12-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.248e+12, 's^-1'), n=0.09, Ea=(36820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 <=> C6H12-1 + H""",
)

entry(
    index = 1644,
    label = "C6H13-2 <=> C6H12-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.271e+13, 's^-1'), n=-0.09, Ea=(35650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 <=> C6H12-2 + H""",
)

entry(
    index = 1645,
    label = "C6H13-3 <=> C4H8-1 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.202e+19, 's^-1'), n=-1.76, Ea=(30450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 <=> C4H8-1 + C2H5""",
)

entry(
    index = 1646,
    label = "C6H13-3 <=> C5H10-1 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.671e+16, 's^-1'), n=-0.93, Ea=(31480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 <=> C5H10-1 + CH3""",
)

entry(
    index = 1647,
    label = "C6H13-3 <=> C6H12-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.271e+13, 's^-1'), n=-0.09, Ea=(35650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 <=> C6H12-2 + H""",
)

entry(
    index = 1648,
    label = "C6H13-3 <=> C6H12-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.659e+12, 's^-1'), n=-0.02, Ea=(35740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 <=> C6H12-3 + H""",
)

entry(
    index = 1649,
    label = "C6H13-1 <=> C6H13-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.478e+08, 's^-1'), n=1.62, Ea=(38760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 <=> C6H13-2""",
)

entry(
    index = 1650,
    label = "C6H13-1 <=> C6H13-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.495e+09, 's^-1'), n=0.97, Ea=(33760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 <=> C6H13-3""",
)

entry(
    index = 1651,
    label = "C6H12-1 <=> NC3H7 + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(71000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 <=> NC3H7 + C3H5-A""",
)

entry(
    index = 1652,
    label = "C6H12-2 <=> C2H5 + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(71000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 <=> C2H5 + C4H71-3""",
)

entry(
    index = 1653,
    label = "C6H12-3 <=> CH3 + C5H91-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(71000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 <=> CH3 + C5H91-3""",
)

entry(
    index = 1654,
    label = "C6H12-1 + OH => C5H11-1 + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + OH => C5H11-1 + CH2O""",
)

entry(
    index = 1655,
    label = "C6H12-2 + OH => PC4H9 + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + OH => PC4H9 + CH3CHO""",
)

entry(
    index = 1656,
    label = "C6H12-3 + OH => PC4H9 + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + OH => PC4H9 + CH3CHO""",
)

entry(
    index = 1657,
    label = "C6H12-1 + O => C5H11-1 + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O => C5H11-1 + HCO""",
)

entry(
    index = 1658,
    label = "C6H12-2 + O => PC4H9 + CH3CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O => PC4H9 + CH3CO""",
)

entry(
    index = 1659,
    label = "C6H12-1 + H <=> C6H111-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + H <=> C6H111-3 + H2""",
)

entry(
    index = 1660,
    label = "C6H12-1 + H <=> C6H111-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + H <=> C6H111-4 + H2""",
)

entry(
    index = 1661,
    label = "C6H12-1 + H <=> C6H111-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + H <=> C6H111-5 + H2""",
)

entry(
    index = 1662,
    label = "C6H12-1 + H <=> C6H111-6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665000, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + H <=> C6H111-6 + H2""",
)

entry(
    index = 1663,
    label = "C6H12-1 + OH <=> C6H111-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + OH <=> C6H111-3 + H2O""",
)

entry(
    index = 1664,
    label = "C6H12-1 + OH <=> C6H111-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + OH <=> C6H111-4 + H2O""",
)

entry(
    index = 1665,
    label = "C6H12-1 + OH <=> C6H111-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + OH <=> C6H111-5 + H2O""",
)

entry(
    index = 1666,
    label = "C6H12-1 + OH <=> C6H111-6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + OH <=> C6H111-6 + H2O""",
)

entry(
    index = 1667,
    label = "C6H12-1 + CH3 <=> C6H111-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3 <=> C6H111-3 + CH4""",
)

entry(
    index = 1668,
    label = "C6H12-1 + CH3 <=> C6H111-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3 <=> C6H111-4 + CH4""",
)

entry(
    index = 1669,
    label = "C6H12-1 + CH3 <=> C6H111-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3 <=> C6H111-5 + CH4""",
)

entry(
    index = 1670,
    label = "C6H12-1 + CH3 <=> C6H111-6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4521, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3 <=> C6H111-6 + CH4""",
)

entry(
    index = 1671,
    label = "C6H12-1 + HO2 <=> C6H111-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + HO2 <=> C6H111-3 + H2O2""",
)

entry(
    index = 1672,
    label = "C6H12-1 + HO2 <=> C6H111-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + HO2 <=> C6H111-4 + H2O2""",
)

entry(
    index = 1673,
    label = "C6H12-1 + HO2 <=> C6H111-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + HO2 <=> C6H111-5 + H2O2""",
)

entry(
    index = 1674,
    label = "C6H12-1 + HO2 <=> C6H111-6 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + HO2 <=> C6H111-6 + H2O2""",
)

entry(
    index = 1675,
    label = "C6H12-1 + CH3O2 <=> C6H111-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O2 <=> C6H111-3 + CH3O2H""",
)

entry(
    index = 1676,
    label = "C6H12-1 + CH3O2 <=> C6H111-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O2 <=> C6H111-4 + CH3O2H""",
)

entry(
    index = 1677,
    label = "C6H12-1 + CH3O2 <=> C6H111-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O2 <=> C6H111-5 + CH3O2H""",
)

entry(
    index = 1678,
    label = "C6H12-1 + CH3O2 <=> C6H111-6 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O2 <=> C6H111-6 + CH3O2H""",
)

entry(
    index = 1679,
    label = "C6H12-1 + CH3O <=> C6H111-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O <=> C6H111-3 + CH3OH""",
)

entry(
    index = 1680,
    label = "C6H12-1 + CH3O <=> C6H111-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O <=> C6H111-4 + CH3OH""",
)

entry(
    index = 1681,
    label = "C6H12-1 + CH3O <=> C6H111-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O <=> C6H111-5 + CH3OH""",
)

entry(
    index = 1682,
    label = "C6H12-1 + CH3O <=> C6H111-6 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + CH3O <=> C6H111-6 + CH3OH""",
)

entry(
    index = 1683,
    label = "C6H12-2 + H <=> C6H111-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(173000, 'cm^3/(mol*s)'), n=2.5, Ea=(2492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + H <=> C6H111-3 + H2""",
)

entry(
    index = 1684,
    label = "C6H12-2 + H <=> C6H112-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + H <=> C6H112-4 + H2""",
)

entry(
    index = 1685,
    label = "C6H12-2 + H <=> C6H112-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + H <=> C6H112-5 + H2""",
)

entry(
    index = 1686,
    label = "C6H12-2 + H <=> C6H112-6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665100, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + H <=> C6H112-6 + H2""",
)

entry(
    index = 1687,
    label = "C6H12-2 + OH <=> C6H111-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + OH <=> C6H111-3 + H2O""",
)

entry(
    index = 1688,
    label = "C6H12-2 + OH <=> C6H112-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + OH <=> C6H112-4 + H2O""",
)

entry(
    index = 1689,
    label = "C6H12-2 + OH <=> C6H112-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(467000, 'cm^3/(mol*s)'), n=1.61, Ea=(-35, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + OH <=> C6H112-5 + H2O""",
)

entry(
    index = 1690,
    label = "C6H12-2 + OH <=> C6H112-6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + OH <=> C6H112-6 + H2O""",
)

entry(
    index = 1691,
    label = "C6H12-2 + CH3 <=> C6H111-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.21, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3 <=> C6H111-3 + CH4""",
)

entry(
    index = 1692,
    label = "C6H12-2 + CH3 <=> C6H112-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3 <=> C6H112-4 + CH4""",
)

entry(
    index = 1693,
    label = "C6H12-2 + CH3 <=> C6H112-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3 <=> C6H112-5 + CH4""",
)

entry(
    index = 1694,
    label = "C6H12-2 + CH3 <=> C6H112-6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4521, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3 <=> C6H112-6 + CH4""",
)

entry(
    index = 1695,
    label = "C6H12-2 + HO2 <=> C6H111-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9639, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + HO2 <=> C6H111-3 + H2O2""",
)

entry(
    index = 1696,
    label = "C6H12-2 + HO2 <=> C6H112-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + HO2 <=> C6H112-4 + H2O2""",
)

entry(
    index = 1697,
    label = "C6H12-2 + HO2 <=> C6H112-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + HO2 <=> C6H112-5 + H2O2""",
)

entry(
    index = 1698,
    label = "C6H12-2 + HO2 <=> C6H112-6 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + HO2 <=> C6H112-6 + H2O2""",
)

entry(
    index = 1699,
    label = "C6H12-2 + CH3O2 <=> C6H111-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9639, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O2 <=> C6H111-3 + CH3O2H""",
)

entry(
    index = 1700,
    label = "C6H12-2 + CH3O2 <=> C6H112-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O2 <=> C6H112-4 + CH3O2H""",
)

entry(
    index = 1701,
    label = "C6H12-2 + CH3O2 <=> C6H112-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O2 <=> C6H112-5 + CH3O2H""",
)

entry(
    index = 1702,
    label = "C6H12-2 + CH3O2 <=> C6H112-6 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O2 <=> C6H112-6 + CH3O2H""",
)

entry(
    index = 1703,
    label = "C6H12-2 + CH3O <=> C6H111-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(90, 'cm^3/(mol*s)'), n=2.95, Ea=(11990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O <=> C6H111-3 + CH3OH""",
)

entry(
    index = 1704,
    label = "C6H12-2 + CH3O <=> C6H112-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O <=> C6H112-4 + CH3OH""",
)

entry(
    index = 1705,
    label = "C6H12-2 + CH3O <=> C6H112-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O <=> C6H112-5 + CH3OH""",
)

entry(
    index = 1706,
    label = "C6H12-2 + CH3O <=> C6H112-6 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + CH3O <=> C6H112-6 + CH3OH""",
)

entry(
    index = 1707,
    label = "C6H12-3 + H <=> C6H113-1 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.33e+06, 'cm^3/(mol*s)'),
        n = 2.54,
        Ea = (6756, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + H <=> C6H113-1 + H2""",
)

entry(
    index = 1708,
    label = "C6H12-3 + H <=> C6H112-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(675200, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + H <=> C6H112-4 + H2""",
)

entry(
    index = 1709,
    label = "C6H12-3 + OH <=> C6H113-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + OH <=> C6H113-1 + H2O""",
)

entry(
    index = 1710,
    label = "C6H12-3 + OH <=> C6H112-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(55280, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + OH <=> C6H112-4 + H2O""",
)

entry(
    index = 1711,
    label = "C6H12-3 + CH3 <=> C6H113-1 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.9042, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + CH3 <=> C6H113-1 + CH4""",
)

entry(
    index = 1712,
    label = "C6H12-3 + CH3 <=> C6H112-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.38, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + CH3 <=> C6H112-4 + CH4""",
)

entry(
    index = 1713,
    label = "C6H12-3 + HO2 <=> C6H113-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + HO2 <=> C6H113-1 + H2O2""",
)

entry(
    index = 1714,
    label = "C6H12-3 + HO2 <=> C6H112-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + HO2 <=> C6H112-4 + H2O2""",
)

entry(
    index = 1715,
    label = "C6H12-3 + CH3O2 <=> C6H113-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + CH3O2 <=> C6H113-1 + CH3O2H""",
)

entry(
    index = 1716,
    label = "C6H12-3 + CH3O2 <=> C6H112-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + CH3O2 <=> C6H112-4 + CH3O2H""",
)

entry(
    index = 1717,
    label = "C6H12-3 + CH3O <=> C6H113-1 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.34e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + CH3O <=> C6H113-1 + CH3OH""",
)

entry(
    index = 1718,
    label = "C6H12-3 + CH3O <=> C6H112-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(80, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + CH3O <=> C6H112-4 + CH3OH""",
)

entry(
    index = 1719,
    label = "C6H111-3 + HO2 <=> C6H11O1-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + HO2 <=> C6H11O1-3 + OH""",
)

entry(
    index = 1720,
    label = "C6H111-3 + CH3O2 <=> C6H11O1-3 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + CH3O2 <=> C6H11O1-3 + CH3O""",
)

entry(
    index = 1721,
    label = "C6H111-3 + C2H5O2 <=> C6H11O1-3 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + C2H5O2 <=> C6H11O1-3 + C2H5O""",
)

entry(
    index = 1722,
    label = "C6H111-6 <=> C6H111-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.113e+12, 's^-1'), n=0, Ea=(31700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-6 <=> C6H111-3""",
)

entry(
    index = 1723,
    label = "C6H112-4 + HO2 <=> C6H11O2-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + HO2 <=> C6H11O2-4 + OH""",
)

entry(
    index = 1724,
    label = "C6H112-4 + CH3O2 <=> C6H11O2-4 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + CH3O2 <=> C6H11O2-4 + CH3O""",
)

entry(
    index = 1725,
    label = "C6H112-4 + C2H5O2 <=> C6H11O2-4 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + C2H5O2 <=> C6H11O2-4 + C2H5O""",
)

entry(
    index = 1726,
    label = "C6H11O1-3 <=> C2H3CHO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.427e+20, 's^-1'), n=-2.04, Ea=(11230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O1-3 <=> C2H3CHO + NC3H7""",
)

entry(
    index = 1727,
    label = "C6H11O1-3 <=> NC3H7CHO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+18, 's^-1'), n=-1.63, Ea=(23410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O1-3 <=> NC3H7CHO + C2H3""",
)

entry(
    index = 1728,
    label = "C6H11O2-4 <=> SC3H5CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.449e+19, 's^-1'), n=-1.92, Ea=(10760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O2-4 <=> SC3H5CHO + C2H5""",
)

entry(
    index = 1729,
    label = "C6H11O2-4 <=> C2H5CHO + C3H5-S",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.346e+22, 's^-1'), n=-2.58, Ea=(29310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O2-4 <=> C2H5CHO + C3H5-S""",
)

entry(
    index = 1730,
    label = "C6H13O2-1 <=> C6H13-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.15e+20, 's^-1'), n=-1.71, Ea=(35790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 <=> C6H13-1 + O2""",
)

entry(
    index = 1731,
    label = "C6H13O2-2 <=> C6H13-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.18e+23, 's^-1'), n=-2.33, Ea=(38040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H13-2 + O2""",
)

entry(
    index = 1732,
    label = "C6H13O2-3 <=> C6H13-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.18e+23, 's^-1'), n=-2.33, Ea=(38040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H13-3 + O2""",
)

entry(
    index = 1733,
    label = "C6H13-1 + C6H13O2-1 <=> C6H13O-1 + C6H13O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 + C6H13O2-1 <=> C6H13O-1 + C6H13O-1""",
)

entry(
    index = 1734,
    label = "C6H13-1 + C6H13O2-2 <=> C6H13O-1 + C6H13O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 + C6H13O2-2 <=> C6H13O-1 + C6H13O-2""",
)

entry(
    index = 1735,
    label = "C6H13-1 + C6H13O2-3 <=> C6H13O-1 + C6H13O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-1 + C6H13O2-3 <=> C6H13O-1 + C6H13O-3""",
)

entry(
    index = 1736,
    label = "C6H13-2 + C6H13O2-1 <=> C6H13O-2 + C6H13O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + C6H13O2-1 <=> C6H13O-2 + C6H13O-1""",
)

entry(
    index = 1737,
    label = "C6H13-2 + C6H13O2-2 <=> C6H13O-2 + C6H13O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + C6H13O2-2 <=> C6H13O-2 + C6H13O-2""",
)

entry(
    index = 1738,
    label = "C6H13-2 + C6H13O2-3 <=> C6H13O-2 + C6H13O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-2 + C6H13O2-3 <=> C6H13O-2 + C6H13O-3""",
)

entry(
    index = 1739,
    label = "C6H13-3 + C6H13O2-1 <=> C6H13O-3 + C6H13O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + C6H13O2-1 <=> C6H13O-3 + C6H13O-1""",
)

entry(
    index = 1740,
    label = "C6H13-3 + C6H13O2-2 <=> C6H13O-3 + C6H13O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + C6H13O2-2 <=> C6H13O-3 + C6H13O-2""",
)

entry(
    index = 1741,
    label = "C6H13-3 + C6H13O2-3 <=> C6H13O-3 + C6H13O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13-3 + C6H13O2-3 <=> C6H13O-3 + C6H13O-3""",
)

entry(
    index = 1742,
    label = "C6H13O2-1 <=> C6H12-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 <=> C6H12-1 + HO2""",
)

entry(
    index = 1743,
    label = "C6H13O2-2 <=> C6H12-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.075e+42, 's^-1'), n=-9.41, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12-1 + HO2""",
)

entry(
    index = 1744,
    label = "C6H13O2-2 <=> C6H12-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.044e+38, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12-2 + HO2""",
)

entry(
    index = 1745,
    label = "C6H13O2-3 <=> C6H12-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.009e+39, 's^-1'), n=-8.11, Ea=(40490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H12-3 + HO2""",
)

entry(
    index = 1746,
    label = "C6H13O2-1 <=> C6H12OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 <=> C6H12OOH1-2""",
)

entry(
    index = 1747,
    label = "C6H13O2-1 <=> C6H12OOH1-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 <=> C6H12OOH1-3""",
)

entry(
    index = 1748,
    label = "C6H13O2-1 <=> C6H12OOH1-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(19050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 <=> C6H12OOH1-4""",
)

entry(
    index = 1749,
    label = "C6H13O2-1 <=> C6H12OOH1-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.906e+08, 's^-1'), n=0, Ea=(22050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 <=> C6H12OOH1-5""",
)

entry(
    index = 1750,
    label = "C6H13O2-2 <=> C6H12OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12OOH2-1""",
)

entry(
    index = 1751,
    label = "C6H13O2-2 <=> C6H12OOH2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12OOH2-3""",
)

entry(
    index = 1752,
    label = "C6H13O2-2 <=> C6H12OOH2-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12OOH2-4""",
)

entry(
    index = 1753,
    label = "C6H13O2-2 <=> C6H12OOH2-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(19050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12OOH2-5""",
)

entry(
    index = 1754,
    label = "C6H13O2-2 <=> C6H12OOH2-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.86e+08, 's^-1'), n=0, Ea=(25550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 <=> C6H12OOH2-6""",
)

entry(
    index = 1755,
    label = "C6H13O2-3 <=> C6H12OOH3-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+10, 's^-1'), n=0, Ea=(24400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H12OOH3-1""",
)

entry(
    index = 1756,
    label = "C6H13O2-3 <=> C6H12OOH3-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H12OOH3-2""",
)

entry(
    index = 1757,
    label = "C6H13O2-3 <=> C6H12OOH3-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H12OOH3-4""",
)

entry(
    index = 1758,
    label = "C6H13O2-3 <=> C6H12OOH3-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H12OOH3-5""",
)

entry(
    index = 1759,
    label = "C6H13O2-3 <=> C6H12OOH3-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.688e+09, 's^-1'), n=0, Ea=(22350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 <=> C6H12OOH3-6""",
)

entry(
    index = 1760,
    label = "C6H13O2-1 + HO2 <=> C6H13O2H-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + HO2 <=> C6H13O2H-1 + O2""",
)

entry(
    index = 1761,
    label = "C6H13O2-2 + HO2 <=> C6H13O2H-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + HO2 <=> C6H13O2H-2 + O2""",
)

entry(
    index = 1762,
    label = "C6H13O2-3 + HO2 <=> C6H13O2H-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + HO2 <=> C6H13O2H-3 + O2""",
)

entry(
    index = 1763,
    label = "C6H13O2-1 + H2O2 <=> C6H13O2H-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + H2O2 <=> C6H13O2H-1 + HO2""",
)

entry(
    index = 1764,
    label = "C6H13O2-2 + H2O2 <=> C6H13O2H-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + H2O2 <=> C6H13O2H-2 + HO2""",
)

entry(
    index = 1765,
    label = "C6H13O2-3 + H2O2 <=> C6H13O2H-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + H2O2 <=> C6H13O2H-3 + HO2""",
)

entry(
    index = 1766,
    label = "C6H13O2-1 + CH3O2 => C6H13O-1 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + CH3O2 => C6H13O-1 + CH3O + O2""",
)

entry(
    index = 1767,
    label = "C6H13O2-2 + CH3O2 => C6H13O-2 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + CH3O2 => C6H13O-2 + CH3O + O2""",
)

entry(
    index = 1768,
    label = "C6H13O2-3 + CH3O2 => C6H13O-3 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + CH3O2 => C6H13O-3 + CH3O + O2""",
)

entry(
    index = 1769,
    label = "C6H13O2-1 + C6H13O2-1 => O2 + C6H13O-1 + C6H13O-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + C6H13O2-1 => O2 + C6H13O-1 + C6H13O-1""",
)

entry(
    index = 1770,
    label = "C6H13O2-1 + C6H13O2-2 => O2 + C6H13O-1 + C6H13O-2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + C6H13O2-2 => O2 + C6H13O-1 + C6H13O-2""",
)

entry(
    index = 1771,
    label = "C6H13O2-1 + C6H13O2-3 => O2 + C6H13O-1 + C6H13O-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-1 + C6H13O2-3 => O2 + C6H13O-1 + C6H13O-3""",
)

entry(
    index = 1772,
    label = "C6H13O2-2 + C6H13O2-2 => O2 + C6H13O-2 + C6H13O-2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + C6H13O2-2 => O2 + C6H13O-2 + C6H13O-2""",
)

entry(
    index = 1773,
    label = "C6H13O2-2 + C6H13O2-3 => O2 + C6H13O-2 + C6H13O-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-2 + C6H13O2-3 => O2 + C6H13O-2 + C6H13O-3""",
)

entry(
    index = 1774,
    label = "C6H13O2-3 + C6H13O2-3 => O2 + C6H13O-3 + C6H13O-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H13O2-3 + C6H13O2-3 => O2 + C6H13O-3 + C6H13O-3""",
)

entry(
    index = 1775,
    label = "C6H13O2H-1 <=> C6H13O-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2H-1 <=> C6H13O-1 + OH""",
)

entry(
    index = 1776,
    label = "C6H13O2H-2 <=> C6H13O-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2H-2 <=> C6H13O-2 + OH""",
)

entry(
    index = 1777,
    label = "C6H13O2H-3 <=> C6H13O-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O2H-3 <=> C6H13O-3 + OH""",
)

entry(
    index = 1778,
    label = "C6H13O-1 <=> C5H11-1 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.226e+20, 's^-1'), n=-2.08, Ea=(24830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O-1 <=> C5H11-1 + CH2O""",
)

entry(
    index = 1779,
    label = "C6H13O-2 <=> PC4H9 + CH3CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.365e+22, 's^-1'), n=-2.61, Ea=(19620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O-2 <=> PC4H9 + CH3CHO""",
)

entry(
    index = 1780,
    label = "C6H13O-3 <=> C2H5CHO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.045e+17, 's^-1'), n=-1.17, Ea=(18170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H13O-3 <=> C2H5CHO + NC3H7""",
)

entry(
    index = 1781,
    label = "C6H12OOH1-2 => C6H12O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-2 => C6H12O1-2 + OH""",
)

entry(
    index = 1782,
    label = "C6H12OOH1-3 => C6H12O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-3 => C6H12O1-3 + OH""",
)

entry(
    index = 1783,
    label = "C6H12OOH1-4 => C6H12O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-4 => C6H12O1-4 + OH""",
)

entry(
    index = 1784,
    label = "C6H12OOH1-5 => C6H12O1-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.172e+09, 's^-1'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-5 => C6H12O1-5 + OH""",
)

entry(
    index = 1785,
    label = "C6H12OOH2-1 => C6H12O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-1 => C6H12O1-2 + OH""",
)

entry(
    index = 1786,
    label = "C6H12OOH2-3 => C6H12O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-3 => C6H12O2-3 + OH""",
)

entry(
    index = 1787,
    label = "C6H12OOH2-4 => C6H12O2-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-4 => C6H12O2-4 + OH""",
)

entry(
    index = 1788,
    label = "C6H12OOH2-5 => C6H12O2-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-5 => C6H12O2-5 + OH""",
)

entry(
    index = 1789,
    label = "C6H12OOH2-6 => C6H12O1-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.172e+09, 's^-1'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-6 => C6H12O1-5 + OH""",
)

entry(
    index = 1790,
    label = "C6H12OOH3-2 => C6H12O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-2 => C6H12O2-3 + OH""",
)

entry(
    index = 1791,
    label = "C6H12OOH3-4 => C6H12O3-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-4 => C6H12O3-4 + OH""",
)

entry(
    index = 1792,
    label = "C6H12OOH3-1 => C6H12O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-1 => C6H12O1-3 + OH""",
)

entry(
    index = 1793,
    label = "C6H12OOH3-5 => C6H12O2-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-5 => C6H12O2-4 + OH""",
)

entry(
    index = 1794,
    label = "C6H12OOH3-6 => C6H12O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-6 => C6H12O1-4 + OH""",
)

entry(
    index = 1795,
    label = "C6H12OOH1-3 => OH + CH2O + C5H10-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.703e+13, 's^-1'), n=-0.16, Ea=(30090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-3 => OH + CH2O + C5H10-1""",
)

entry(
    index = 1796,
    label = "C6H12OOH2-4 => OH + CH3CHO + C4H8-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.933e+18, 's^-1'), n=-1.7, Ea=(24080, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-4 => OH + CH3CHO + C4H8-1""",
)

entry(
    index = 1797,
    label = "C6H12OOH3-5 => OH + C2H5CHO + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.171e+17, 's^-1'), n=-1.31, Ea=(28880, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-5 => OH + C2H5CHO + C3H6""",
)

entry(
    index = 1798,
    label = "C6H12OOH3-1 => OH + CH3CHO + C2H4 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.326e+18, 's^-1'), n=-1.74, Ea=(27420, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-1 => OH + CH3CHO + C2H4 + C2H4""",
)

entry(
    index = 1799,
    label = "C6H12OOH1-2O2 <=> C6H12OOH1-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.868e+22, 's^-1'), n=-2.31, Ea=(37980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-2O2 <=> C6H12OOH1-2 + O2""",
)

entry(
    index = 1800,
    label = "C6H12OOH1-3O2 <=> C6H12OOH1-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.868e+22, 's^-1'), n=-2.31, Ea=(37980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-3O2 <=> C6H12OOH1-3 + O2""",
)

entry(
    index = 1801,
    label = "C6H12OOH1-4O2 <=> C6H12OOH1-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.868e+22, 's^-1'), n=-2.31, Ea=(37980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-4O2 <=> C6H12OOH1-4 + O2""",
)

entry(
    index = 1802,
    label = "C6H12OOH1-5O2 <=> C6H12OOH1-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.868e+22, 's^-1'), n=-2.31, Ea=(37980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-5O2 <=> C6H12OOH1-5 + O2""",
)

entry(
    index = 1803,
    label = "C6H12OOH2-1O2 <=> C6H12OOH2-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.463e+20, 's^-1'), n=-1.64, Ea=(35670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-1O2 <=> C6H12OOH2-1 + O2""",
)

entry(
    index = 1804,
    label = "C6H12OOH2-3O2 <=> C6H12OOH2-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.792e+22, 's^-1'), n=-2.33, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-3O2 <=> C6H12OOH2-3 + O2""",
)

entry(
    index = 1805,
    label = "C6H12OOH2-4O2 <=> C6H12OOH2-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.792e+22, 's^-1'), n=-2.33, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-4O2 <=> C6H12OOH2-4 + O2""",
)

entry(
    index = 1806,
    label = "C6H12OOH2-5O2 <=> C6H12OOH2-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.792e+22, 's^-1'), n=-2.33, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-5O2 <=> C6H12OOH2-5 + O2""",
)

entry(
    index = 1807,
    label = "C6H12OOH2-6O2 <=> C6H12OOH2-6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.463e+20, 's^-1'), n=-1.64, Ea=(35670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-6O2 <=> C6H12OOH2-6 + O2""",
)

entry(
    index = 1808,
    label = "C6H12OOH3-1O2 <=> C6H12OOH3-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.463e+20, 's^-1'), n=-1.64, Ea=(35670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-1O2 <=> C6H12OOH3-1 + O2""",
)

entry(
    index = 1809,
    label = "C6H12OOH3-2O2 <=> C6H12OOH3-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.792e+22, 's^-1'), n=-2.33, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-2O2 <=> C6H12OOH3-2 + O2""",
)

entry(
    index = 1810,
    label = "C6H12OOH3-4O2 <=> C6H12OOH3-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.792e+22, 's^-1'), n=-2.33, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-4O2 <=> C6H12OOH3-4 + O2""",
)

entry(
    index = 1811,
    label = "C6H12OOH3-5O2 <=> C6H12OOH3-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.792e+22, 's^-1'), n=-2.33, Ea=(37960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-5O2 <=> C6H12OOH3-5 + O2""",
)

entry(
    index = 1812,
    label = "C6H12OOH3-6O2 <=> C6H12OOH3-6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.463e+20, 's^-1'), n=-1.64, Ea=(35670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-6O2 <=> C6H12OOH3-6 + O2""",
)

entry(
    index = 1813,
    label = "C6H12OOH1-2O2 <=> NC6KET12 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-2O2 <=> NC6KET12 + OH""",
)

entry(
    index = 1814,
    label = "C6H12OOH1-3O2 <=> NC6KET13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(21400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-3O2 <=> NC6KET13 + OH""",
)

entry(
    index = 1815,
    label = "C6H12OOH1-4O2 <=> NC6KET14 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(19350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-4O2 <=> NC6KET14 + OH""",
)

entry(
    index = 1816,
    label = "C6H12OOH1-5O2 <=> NC6KET15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.906e+08, 's^-1'), n=0, Ea=(22550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH1-5O2 <=> NC6KET15 + OH""",
)

entry(
    index = 1817,
    label = "C6H12OOH2-1O2 <=> NC6KET21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-1O2 <=> NC6KET21 + OH""",
)

entry(
    index = 1818,
    label = "C6H12OOH2-3O2 <=> NC6KET23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-3O2 <=> NC6KET23 + OH""",
)

entry(
    index = 1819,
    label = "C6H12OOH2-4O2 <=> NC6KET24 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-4O2 <=> NC6KET24 + OH""",
)

entry(
    index = 1820,
    label = "C6H12OOH2-5O2 <=> NC6KET25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(16050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-5O2 <=> NC6KET25 + OH""",
)

entry(
    index = 1821,
    label = "C6H12OOH2-6O2 <=> NC6KET26 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.906e+08, 's^-1'), n=0, Ea=(22550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH2-6O2 <=> NC6KET26 + OH""",
)

entry(
    index = 1822,
    label = "C6H12OOH3-1O2 <=> NC6KET31 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-1O2 <=> NC6KET31 + OH""",
)

entry(
    index = 1823,
    label = "C6H12OOH3-2O2 <=> NC6KET32 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-2O2 <=> NC6KET32 + OH""",
)

entry(
    index = 1824,
    label = "C6H12OOH3-4O2 <=> NC6KET34 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-4O2 <=> NC6KET34 + OH""",
)

entry(
    index = 1825,
    label = "C6H12OOH3-5O2 <=> NC6KET35 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-5O2 <=> NC6KET35 + OH""",
)

entry(
    index = 1826,
    label = "C6H12OOH3-6O2 <=> NC6KET36 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(16050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OOH3-6O2 <=> NC6KET36 + OH""",
)

entry(
    index = 1827,
    label = "NC6KET12 => NC4H9CHO + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET12 => NC4H9CHO + HCO + OH""",
)

entry(
    index = 1828,
    label = "NC6KET13 => NC3H7CHO + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET13 => NC3H7CHO + CH2CHO + OH""",
)

entry(
    index = 1829,
    label = "NC6KET14 => C2H5CHO + CH2CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET14 => C2H5CHO + CH2CH2CHO + OH""",
)

entry(
    index = 1830,
    label = "NC6KET15 => CH3CHO + C3H6CHO-1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET15 => CH3CHO + C3H6CHO-1 + OH""",
)

entry(
    index = 1831,
    label = "NC6KET21 => CH2O + NC4H9CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET21 => CH2O + NC4H9CO + OH""",
)

entry(
    index = 1832,
    label = "NC6KET23 => NC3H7CHO + CH3CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET23 => NC3H7CHO + CH3CO + OH""",
)

entry(
    index = 1833,
    label = "NC6KET24 => C2H5CHO + CH3COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET24 => C2H5CHO + CH3COCH2 + OH""",
)

entry(
    index = 1834,
    label = "NC6KET25 => CH3CHO + CH2CH2COCH3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET25 => CH3CHO + CH2CH2COCH3 + OH""",
)

entry(
    index = 1835,
    label = "NC6KET26 => CH2O + C3H6COCH3-1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET26 => CH2O + C3H6COCH3-1 + OH""",
)

entry(
    index = 1836,
    label = "NC6KET31 => CH2O + NC3H7COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET31 => CH2O + NC3H7COCH2 + OH""",
)

entry(
    index = 1837,
    label = "NC6KET32 => CH3CHO + NC3H7CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET32 => CH3CHO + NC3H7CO + OH""",
)

entry(
    index = 1838,
    label = "NC6KET34 => C2H5CHO + C2H5CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET34 => C2H5CHO + C2H5CO + OH""",
)

entry(
    index = 1839,
    label = "NC6KET35 => CH3CHO + C2H5COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET35 => CH3CHO + C2H5COCH2 + OH""",
)

entry(
    index = 1840,
    label = "NC6KET36 => CH2O + C2H5COC2H4P + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6KET36 => CH2O + C2H5COC2H4P + OH""",
)

entry(
    index = 1841,
    label = "C6H12O1-2 + OH => C2H3CHO + NC3H7 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-2 + OH => C2H3CHO + NC3H7 + H2O""",
)

entry(
    index = 1842,
    label = "C6H12O1-3 + OH => C5H10-1 + HCO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-3 + OH => C5H10-1 + HCO + H2O""",
)

entry(
    index = 1843,
    label = "C6H12O1-4 + OH => C4H8-1 + CH2CHO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-4 + OH => C4H8-1 + CH2CHO + H2O""",
)

entry(
    index = 1844,
    label = "C6H12O1-5 + OH => C3H6 + CH2CH2CHO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-5 + OH => C3H6 + CH2CH2CHO + H2O""",
)

entry(
    index = 1845,
    label = "C6H12O2-3 + OH => C2H3COCH3 + C2H5 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-3 + OH => C2H3COCH3 + C2H5 + H2O""",
)

entry(
    index = 1846,
    label = "C6H12O2-4 + OH => C4H8-1 + CH3CO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-4 + OH => C4H8-1 + CH3CO + H2O""",
)

entry(
    index = 1847,
    label = "C6H12O2-5 + OH => C3H6 + CH3COCH2 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-5 + OH => C3H6 + CH3COCH2 + H2O""",
)

entry(
    index = 1848,
    label = "C6H12O3-4 + OH => C2H5COC2H3 + CH3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O3-4 + OH => C2H5COC2H3 + CH3 + H2O""",
)

entry(
    index = 1849,
    label = "C6H12O1-2 + OH => CH2CO + PC4H9 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-2 + OH => CH2CO + PC4H9 + H2O""",
)

entry(
    index = 1850,
    label = "C6H12O1-3 + OH => C2H4 + NC3H7CO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-3 + OH => C2H4 + NC3H7CO + H2O""",
)

entry(
    index = 1851,
    label = "C6H12O1-4 + OH => C2H4 + C2H5COCH2 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-4 + OH => C2H4 + C2H5COCH2 + H2O""",
)

entry(
    index = 1852,
    label = "C6H12O1-5 + OH => C2H4 + CH2CH2COCH3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-5 + OH => C2H4 + CH2CH2COCH3 + H2O""",
)

entry(
    index = 1853,
    label = "C6H12O2-3 + OH => CH3CHCO + NC3H7 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-3 + OH => CH3CHCO + NC3H7 + H2O""",
)

entry(
    index = 1854,
    label = "C6H12O2-4 + OH => C3H6 + C2H5CO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-4 + OH => C3H6 + C2H5CO + H2O""",
)

entry(
    index = 1855,
    label = "C6H12O2-5 + OH => CH3CHO + C4H71-3 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-5 + OH => CH3CHO + C4H71-3 + H2O""",
)

entry(
    index = 1856,
    label = "C6H12O3-4 + OH => C2H5CHO + C3H5-S + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O3-4 + OH => C2H5CHO + C3H5-S + H2O""",
)

entry(
    index = 1857,
    label = "C6H12O1-2 + HO2 => C2H3CHO + NC3H7 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-2 + HO2 => C2H3CHO + NC3H7 + H2O2""",
)

entry(
    index = 1858,
    label = "C6H12O1-3 + HO2 => C5H10-1 + HCO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-3 + HO2 => C5H10-1 + HCO + H2O2""",
)

entry(
    index = 1859,
    label = "C6H12O1-4 + HO2 => C4H8-1 + CH2CHO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-4 + HO2 => C4H8-1 + CH2CHO + H2O2""",
)

entry(
    index = 1860,
    label = "C6H12O1-5 + HO2 => C3H6 + CH2CH2CHO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-5 + HO2 => C3H6 + CH2CH2CHO + H2O2""",
)

entry(
    index = 1861,
    label = "C6H12O2-3 + HO2 => C2H3COCH3 + C2H5 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-3 + HO2 => C2H3COCH3 + C2H5 + H2O2""",
)

entry(
    index = 1862,
    label = "C6H12O2-4 + HO2 => C4H8-1 + CH3CO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-4 + HO2 => C4H8-1 + CH3CO + H2O2""",
)

entry(
    index = 1863,
    label = "C6H12O2-5 + HO2 => C3H6 + CH3COCH2 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-5 + HO2 => C3H6 + CH3COCH2 + H2O2""",
)

entry(
    index = 1864,
    label = "C6H12O3-4 + HO2 => C2H5COC2H3 + CH3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O3-4 + HO2 => C2H5COC2H3 + CH3 + H2O2""",
)

entry(
    index = 1865,
    label = "C6H12O1-2 + HO2 => CH2CO + PC4H9 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-2 + HO2 => CH2CO + PC4H9 + H2O2""",
)

entry(
    index = 1866,
    label = "C6H12O1-3 + HO2 => C2H4 + NC3H7CO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-3 + HO2 => C2H4 + NC3H7CO + H2O2""",
)

entry(
    index = 1867,
    label = "C6H12O1-4 + HO2 => C2H4 + C2H5COCH2 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-4 + HO2 => C2H4 + C2H5COCH2 + H2O2""",
)

entry(
    index = 1868,
    label = "C6H12O1-5 + HO2 => C2H4 + CH2CH2COCH3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O1-5 + HO2 => C2H4 + CH2CH2COCH3 + H2O2""",
)

entry(
    index = 1869,
    label = "C6H12O2-3 + HO2 => CH3CHCO + NC3H7 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-3 + HO2 => CH3CHCO + NC3H7 + H2O2""",
)

entry(
    index = 1870,
    label = "C6H12O2-4 + HO2 => C3H6 + C2H5CO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-4 + HO2 => C3H6 + C2H5CO + H2O2""",
)

entry(
    index = 1871,
    label = "C6H12O2-5 + HO2 => CH3CHO + C4H71-3 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O2-5 + HO2 => CH3CHO + C4H71-3 + H2O2""",
)

entry(
    index = 1872,
    label = "C6H12O3-4 + HO2 => C2H5CHO + C3H5-S + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12O3-4 + HO2 => C2H5CHO + C3H5-S + H2O2""",
)

entry(
    index = 1873,
    label = "C6H12-1 + OH <=> C6H12OH-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + OH <=> C6H12OH-1""",
)

entry(
    index = 1874,
    label = "C6H12OH-1 + O2 <=> O2C6H12OH-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-1 + O2 <=> O2C6H12OH-1""",
)

entry(
    index = 1875,
    label = "O2C6H12OH-1 <=> NC4H9CHO + CH2O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 's^-1'), n=0, Ea=(27800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C6H12OH-1 <=> NC4H9CHO + CH2O + OH""",
)

entry(
    index = 1876,
    label = "C6H12OH-2 <=> C6H12-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.065e+16, 's^-1'), n=-1, Ea=(29330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-2 <=> C6H12-2 + OH""",
)

entry(
    index = 1877,
    label = "O2C6H12OH-2 <=> C6H12OH-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.816e+21, 's^-1'), n=-2.02, Ea=(37830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C6H12OH-2 <=> C6H12OH-2 + O2""",
)

entry(
    index = 1878,
    label = "O2C6H12OH-2 => NC3H7CHO + CH3CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C6H12OH-2 => NC3H7CHO + CH3CHO + OH""",
)

entry(
    index = 1879,
    label = "C6H12OH-3 <=> C6H12-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.903e+15, 's^-1'), n=-0.93, Ea=(29420, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-3 <=> C6H12-3 + OH""",
)

entry(
    index = 1880,
    label = "O2C6H12OH-3 <=> C6H12OH-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.816e+21, 's^-1'), n=-2.02, Ea=(37830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C6H12OH-3 <=> C6H12OH-3 + O2""",
)

entry(
    index = 1881,
    label = "O2C6H12OH-3 => OH + C2H5CHO + C2H5CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C6H12OH-3 => OH + C2H5CHO + C2H5CHO""",
)

entry(
    index = 1882,
    label = "NC4H9CHO + O2 <=> NC4H9CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0.5, Ea=(42200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + O2 <=> NC4H9CO + HO2""",
)

entry(
    index = 1883,
    label = "NC4H9CHO + OH <=> NC4H9CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + OH <=> NC4H9CO + H2O""",
)

entry(
    index = 1884,
    label = "NC4H9CHO + H <=> NC4H9CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(4200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + H <=> NC4H9CO + H2""",
)

entry(
    index = 1885,
    label = "NC4H9CHO + O <=> NC4H9CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(1790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + O <=> NC4H9CO + OH""",
)

entry(
    index = 1886,
    label = "NC4H9CHO + HO2 <=> NC4H9CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + HO2 <=> NC4H9CO + H2O2""",
)

entry(
    index = 1887,
    label = "NC4H9CHO + CH3 <=> NC4H9CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + CH3 <=> NC4H9CO + CH4""",
)

entry(
    index = 1888,
    label = "NC4H9CHO + CH3O <=> NC4H9CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+11, 'cm^3/(mol*s)'), n=0, Ea=(1280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + CH3O <=> NC4H9CO + CH3OH""",
)

entry(
    index = 1889,
    label = "NC4H9CHO + CH3O2 <=> NC4H9CO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + CH3O2 <=> NC4H9CO + CH3O2H""",
)

entry(
    index = 1890,
    label = "NC4H9CHO + OH <=> C4H8CHO-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + OH <=> C4H8CHO-1 + H2O""",
)

entry(
    index = 1891,
    label = "NC4H9CHO + OH <=> C4H8CHO-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + OH <=> C4H8CHO-2 + H2O""",
)

entry(
    index = 1892,
    label = "NC4H9CHO + OH <=> C4H8CHO-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + OH <=> C4H8CHO-3 + H2O""",
)

entry(
    index = 1893,
    label = "NC4H9CHO + OH <=> C4H8CHO-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9CHO + OH <=> C4H8CHO-4 + H2O""",
)

entry(
    index = 1894,
    label = "NC4H9CO <=> PC4H9 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(9600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9CO <=> PC4H9 + CO""",
)

entry(
    index = 1895,
    label = "C4H8CHO-1 <=> C2H4 + CH2CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.984e+18, 's^-1'), n=-1.6, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8CHO-1 <=> C2H4 + CH2CH2CHO""",
)

entry(
    index = 1896,
    label = "C4H8CHO-2 <=> C3H6 + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.984e+14, 's^-1'), n=-0.76, Ea=(23320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8CHO-2 <=> C3H6 + CH2CHO""",
)

entry(
    index = 1897,
    label = "C4H8CHO-3 <=> C4H8-1 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.797e+14, 's^-1'), n=-0.72, Ea=(24350, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8CHO-3 <=> C4H8-1 + HCO""",
)

entry(
    index = 1898,
    label = "C4H8CHO-3 <=> AC3H5CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.643e+13, 's^-1'), n=-0.36, Ea=(30330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8CHO-3 <=> AC3H5CHO + CH3""",
)

entry(
    index = 1899,
    label = "C4H8CHO-4 <=> C2H3CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.863e+18, 's^-1'), n=-1.3, Ea=(30830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8CHO-4 <=> C2H3CHO + C2H5""",
)

entry(
    index = 1900,
    label = "AC3H5CHO <=> C3H5-A + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.813e+19, 's^-1'), n=-1.08, Ea=(68480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO <=> C3H5-A + HCO""",
)

entry(
    index = 1901,
    label = "AC3H5CHO + OH <=> AC3H5CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO + OH <=> AC3H5CO + H2O""",
)

entry(
    index = 1902,
    label = "AC3H5CHO + OH <=> C2H3CHCHO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO + OH <=> C2H3CHCHO + H2O""",
)

entry(
    index = 1903,
    label = "AC3H5CHO + HO2 <=> AC3H5CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO + HO2 <=> AC3H5CO + H2O2""",
)

entry(
    index = 1904,
    label = "AC3H5CHO + HO2 <=> C2H3CHCHO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9630, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO + HO2 <=> C2H3CHCHO + H2O2""",
)

entry(
    index = 1905,
    label = "AC3H5CHO + CH3O2 <=> AC3H5CO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO + CH3O2 <=> AC3H5CO + CH3O2H""",
)

entry(
    index = 1906,
    label = "AC3H5CHO + CH3O2 <=> C2H3CHCHO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CHO + CH3O2 <=> C2H3CHCHO + CH3O2H""",
)

entry(
    index = 1907,
    label = "AC3H5CO <=> C3H5-A + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.199e+15, 's^-1'), n=-1.09, Ea=(-330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is AC3H5CO <=> C3H5-A + CO""",
)

entry(
    index = 1908,
    label = "C2H3CHCHO + HO2 => C2H3CHO + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.91e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3CHCHO + HO2 => C2H3CHO + HCO + OH""",
)

entry(
    index = 1909,
    label = "NC3H7COCH3 + OH <=> C3H6COCH3-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.065e+07, 'cm^3/(mol*s)'),
        n = 1.73,
        Ea = (753, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + OH <=> C3H6COCH3-1 + H2O""",
)

entry(
    index = 1910,
    label = "NC3H7COCH3 + OH <=> C3H6COCH3-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.615e+07, 'cm^3/(mol*s)'),
        n = 1.64,
        Ea = (-247, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + OH <=> C3H6COCH3-2 + H2O""",
)

entry(
    index = 1911,
    label = "NC3H7COCH3 + OH <=> C3H6COCH3-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + OH <=> C3H6COCH3-3 + H2O""",
)

entry(
    index = 1912,
    label = "NC3H7COCH3 + OH <=> NC3H7COCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + OH <=> NC3H7COCH2 + H2O""",
)

entry(
    index = 1913,
    label = "NC3H7COCH3 + HO2 <=> C3H6COCH3-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + HO2 <=> C3H6COCH3-1 + H2O2""",
)

entry(
    index = 1914,
    label = "NC3H7COCH3 + HO2 <=> C3H6COCH3-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + HO2 <=> C3H6COCH3-2 + H2O2""",
)

entry(
    index = 1915,
    label = "NC3H7COCH3 + HO2 <=> C3H6COCH3-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + HO2 <=> C3H6COCH3-3 + H2O2""",
)

entry(
    index = 1916,
    label = "NC3H7COCH3 + HO2 <=> NC3H7COCH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(14690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + HO2 <=> NC3H7COCH2 + H2O2""",
)

entry(
    index = 1917,
    label = "NC3H7COCH3 + CH3O2 <=> C3H6COCH3-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + CH3O2 <=> C3H6COCH3-1 + CH3O2H""",
)

entry(
    index = 1918,
    label = "NC3H7COCH3 + CH3O2 <=> C3H6COCH3-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + CH3O2 <=> C3H6COCH3-2 + CH3O2H""",
)

entry(
    index = 1919,
    label = "NC3H7COCH3 + CH3O2 <=> C3H6COCH3-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + CH3O2 <=> C3H6COCH3-3 + CH3O2H""",
)

entry(
    index = 1920,
    label = "NC3H7COCH3 + CH3O2 <=> NC3H7COCH2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(17580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH3 + CH3O2 <=> NC3H7COCH2 + CH3O2H""",
)

entry(
    index = 1921,
    label = "C3H6COCH3-1 <=> C2H4 + CH3COCH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.904e+16, 's^-1'), n=-1.21, Ea=(27000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6COCH3-1 <=> C2H4 + CH3COCH2""",
)

entry(
    index = 1922,
    label = "C3H6COCH3-2 <=> C3H6 + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.719e+16, 's^-1'), n=-1.05, Ea=(25590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6COCH3-2 <=> C3H6 + CH3CO""",
)

entry(
    index = 1923,
    label = "C3H6COCH3-3 <=> C2H3COCH3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.615e+15, 's^-1'), n=-0.75, Ea=(32390, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6COCH3-3 <=> C2H3COCH3 + CH3""",
)

entry(
    index = 1924,
    label = "NC3H7COCH2 <=> NC3H7 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.226e+18, 's^-1'), n=-1.4, Ea=(43450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COCH2 <=> NC3H7 + CH2CO""",
)

entry(
    index = 1925,
    label = "C2H5COC2H5 + OH <=> C2H5COC2H4P + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + OH <=> C2H5COC2H4P + H2O""",
)

entry(
    index = 1926,
    label = "C2H5COC2H5 + OH <=> C2H5COC2H4S + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.69e+12, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + OH <=> C2H5COC2H4S + H2O""",
)

entry(
    index = 1927,
    label = "C2H5COC2H5 + HO2 <=> C2H5COC2H4P + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + HO2 <=> C2H5COC2H4P + H2O2""",
)

entry(
    index = 1928,
    label = "C2H5COC2H5 + HO2 <=> C2H5COC2H4S + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + HO2 <=> C2H5COC2H4S + H2O2""",
)

entry(
    index = 1929,
    label = "C2H5COC2H5 + O2 <=> C2H5COC2H4P + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(51310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + O2 <=> C2H5COC2H4P + HO2""",
)

entry(
    index = 1930,
    label = "C2H5COC2H5 + O2 <=> C2H5COC2H4S + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.1e+13, 'cm^3/(mol*s)'), n=0, Ea=(41970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + O2 <=> C2H5COC2H4S + HO2""",
)

entry(
    index = 1931,
    label = "C2H5COC2H5 + H <=> C2H5COC2H4P + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.832e+07, 'cm^3/(mol*s)'), n=2, Ea=(7700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + H <=> C2H5COC2H4P + H2""",
)

entry(
    index = 1932,
    label = "C2H5COC2H5 + H <=> C2H5COC2H4S + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.96e+06, 'cm^3/(mol*s)'), n=2, Ea=(3200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + H <=> C2H5COC2H4S + H2""",
)

entry(
    index = 1933,
    label = "C2H5COC2H5 + C2H3 <=> C2H5COC2H4P + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + C2H3 <=> C2H5COC2H4P + C2H4""",
)

entry(
    index = 1934,
    label = "C2H5COC2H5 + C2H3 <=> C2H5COC2H4S + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+11, 'cm^3/(mol*s)'), n=0, Ea=(5600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + C2H3 <=> C2H5COC2H4S + C2H4""",
)

entry(
    index = 1935,
    label = "C2H5COC2H5 + C2H5 <=> C2H5COC2H4P + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(13400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + C2H5 <=> C2H5COC2H4P + C2H6""",
)

entry(
    index = 1936,
    label = "C2H5COC2H5 + C2H5 <=> C2H5COC2H4S + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+10, 'cm^3/(mol*s)'), n=0, Ea=(8600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + C2H5 <=> C2H5COC2H4S + C2H6""",
)

entry(
    index = 1937,
    label = "C2H5COC2H5 + CH3O <=> C2H5COC2H4P + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.34e+11, 'cm^3/(mol*s)'), n=0, Ea=(6460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + CH3O <=> C2H5COC2H4P + CH3OH""",
)

entry(
    index = 1938,
    label = "C2H5COC2H5 + CH3O <=> C2H5COC2H4S + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.9e+11, 'cm^3/(mol*s)'), n=0, Ea=(2771, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + CH3O <=> C2H5COC2H4S + CH3OH""",
)

entry(
    index = 1939,
    label = "C2H5COC2H5 + CH3O2 <=> C2H5COC2H4P + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.02e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + CH3O2 <=> C2H5COC2H4P + CH3O2H""",
)

entry(
    index = 1940,
    label = "C2H5COC2H5 + CH3O2 <=> C2H5COC2H4S + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H5 + CH3O2 <=> C2H5COC2H4S + CH3O2H""",
)

entry(
    index = 1941,
    label = "C2H5COC2H4P <=> C2H5CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.769e+17, 's^-1'), n=-1.46, Ea=(29540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H4P <=> C2H5CO + C2H4""",
)

entry(
    index = 1942,
    label = "C2H5COC2H4S <=> C2H5COC2H3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.704e+16, 's^-1'), n=-0.82, Ea=(42130, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H4S <=> C2H5COC2H3 + H""",
)

entry(
    index = 1943,
    label = "C2H5COC2H3 + OH <=> C2H5COCH2 + CH2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+10, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + OH <=> C2H5COCH2 + CH2O""",
)

entry(
    index = 1944,
    label = "C2H5COC2H3 + OH <=> PC2H4COC2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (7.55e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + OH <=> PC2H4COC2H3 + H2O""",
)

entry(
    index = 1945,
    label = "C2H5COC2H3 + OH <=> SC2H4COC2H3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + OH <=> SC2H4COC2H3 + H2O""",
)

entry(
    index = 1946,
    label = "C2H5COC2H3 + HO2 <=> C2H5CO + CH2CHO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+09, 'cm^3/(mol*s)'), n=0, Ea=(7949, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + HO2 <=> C2H5CO + CH2CHO + OH""",
)

entry(
    index = 1947,
    label = "C2H5COC2H3 + HO2 <=> PC2H4COC2H3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + HO2 <=> PC2H4COC2H3 + H2O2""",
)

entry(
    index = 1948,
    label = "C2H5COC2H3 + HO2 <=> SC2H4COC2H3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + HO2 <=> SC2H4COC2H3 + H2O2""",
)

entry(
    index = 1949,
    label = "C2H5COC2H3 + CH3O2 <=> C2H5CO + CH2CHO + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.97e+11, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + CH3O2 <=> C2H5CO + CH2CHO + CH3O""",
)

entry(
    index = 1950,
    label = "C2H5COC2H3 + CH3O2 <=> PC2H4COC2H3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + CH3O2 <=> PC2H4COC2H3 + CH3O2H""",
)

entry(
    index = 1951,
    label = "C2H5COC2H3 + CH3O2 <=> SC2H4COC2H3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5COC2H3 + CH3O2 <=> SC2H4COC2H3 + CH3O2H""",
)

entry(
    index = 1952,
    label = "PC2H4COC2H3 <=> C2H3CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.263e+14, 's^-1'), n=0.38, Ea=(21460, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC2H4COC2H3 <=> C2H3CO + C2H4""",
)

entry(
    index = 1953,
    label = "SC2H4COC2H3 <=> CH3CHCO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.637e+16, 's^-1'), n=-0.74, Ea=(54590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is SC2H4COC2H3 <=> CH3CHCO + C2H3""",
)

entry(
    index = 1954,
    label = "NC7H16 <=> H + C7H15-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.68e+88, 's^-1'), n=-21.17, Ea=(142800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> H + C7H15-1""",
)

entry(
    index = 1955,
    label = "NC7H16 <=> H + C7H15-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+88, 's^-1'), n=-21.01, Ea=(139500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> H + C7H15-2""",
)

entry(
    index = 1956,
    label = "NC7H16 <=> H + C7H15-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+88, 's^-1'), n=-21.01, Ea=(139500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> H + C7H15-3""",
)

entry(
    index = 1957,
    label = "NC7H16 <=> H + C7H15-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.5e+87, 's^-1'), n=-21.01, Ea=(139500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> H + C7H15-4""",
)

entry(
    index = 1958,
    label = "NC7H16 <=> C6H13-1 + CH3",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(4.325e+24, 's^-1'), n=-2.12, Ea=(89900, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (4.963e+42, 'cm^3/(mol*s)'),
            n = -7.78,
            Ea = (42800, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.892,
        T3 = (1e+10, 'K'),
        T1 = (2.228, 'K'),
        T2 = (1.798e+09, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> C6H13-1 + CH3""",
)

entry(
    index = 1959,
    label = "NC7H16 <=> C5H11-1 + C2H5",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(6.818e+26, 's^-1'), n=-2.7, Ea=(88910, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (3.753e+48, 'cm^3/(mol*s)'),
            n = -9.46,
            Ea = (41310, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.09,
        T3 = (3.6556, 'K'),
        T1 = (1e+10, 'K'),
        T2 = (9.33e+09, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> C5H11-1 + C2H5""",
)

entry(
    index = 1960,
    label = "NC7H16 <=> PC4H9 + NC3H7",
    degeneracy = 1,
    kinetics = Troe(
        arrheniusHigh = Arrhenius(A=(1.362e+26, 's^-1'), n=-2.53, Ea=(88760, 'cal/mol'), T0=(1, 'K')),
        arrheniusLow = Arrhenius(
            A = (6.509e+48, 'cm^3/(mol*s)'),
            n = -9.57,
            Ea = (41290, 'cal/mol'),
            T0 = (1, 'K'),
        ),
        alpha = 0.911,
        T3 = (1e+10, 'K'),
        T1 = (22.382, 'K'),
        T2 = (5e+09, 'K'),
        efficiencies = {},
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 <=> PC4H9 + NC3H7""",
)

entry(
    index = 1961,
    label = "NC7H16 + H <=> C7H15-1 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(188000, 'cm^3/(mol*s)'), n=2.75, Ea=(6280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> C7H15-1 + H2""",
)

entry(
    index = 1962,
    label = "NC7H16 + H <=> C7H15-2 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> C7H15-2 + H2""",
)

entry(
    index = 1963,
    label = "NC7H16 + H <=> C7H15-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.6e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> C7H15-3 + H2""",
)

entry(
    index = 1964,
    label = "NC7H16 + H <=> C7H15-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + H <=> C7H15-4 + H2""",
)

entry(
    index = 1965,
    label = "NC7H16 + O <=> C7H15-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(193000, 'cm^3/(mol*s)'), n=2.68, Ea=(3716, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> C7H15-1 + OH""",
)

entry(
    index = 1966,
    label = "NC7H16 + O <=> C7H15-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95400, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> C7H15-2 + OH""",
)

entry(
    index = 1967,
    label = "NC7H16 + O <=> C7H15-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(95400, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> C7H15-3 + OH""",
)

entry(
    index = 1968,
    label = "NC7H16 + O <=> C7H15-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47700, 'cm^3/(mol*s)'), n=2.71, Ea=(2106, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O <=> C7H15-4 + OH""",
)

entry(
    index = 1969,
    label = "NC7H16 + OH <=> C7H15-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.57e+07, 'cm^3/(mol*s)'), n=1.8, Ea=(954, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> C7H15-1 + H2O""",
)

entry(
    index = 1970,
    label = "NC7H16 + OH <=> C7H15-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+06, 'cm^3/(mol*s)'), n=2, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> C7H15-2 + H2O""",
)

entry(
    index = 1971,
    label = "NC7H16 + OH <=> C7H15-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.9e+06, 'cm^3/(mol*s)'), n=2, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> C7H15-3 + H2O""",
)

entry(
    index = 1972,
    label = "NC7H16 + OH <=> C7H15-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.45e+06, 'cm^3/(mol*s)'), n=2, Ea=(-596, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + OH <=> C7H15-4 + H2O""",
)

entry(
    index = 1973,
    label = "NC7H16 + HO2 <=> C7H15-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40.8, 'cm^3/(mol*s)'), n=3.59, Ea=(17160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> C7H15-1 + H2O2""",
)

entry(
    index = 1974,
    label = "NC7H16 + HO2 <=> C7H15-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(126.4, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> C7H15-2 + H2O2""",
)

entry(
    index = 1975,
    label = "NC7H16 + HO2 <=> C7H15-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(126.4, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> C7H15-3 + H2O2""",
)

entry(
    index = 1976,
    label = "NC7H16 + HO2 <=> C7H15-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(63.2, 'cm^3/(mol*s)'), n=3.37, Ea=(13720, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + HO2 <=> C7H15-4 + H2O2""",
)

entry(
    index = 1977,
    label = "NC7H16 + CH3 <=> C7H15-1 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.904, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> C7H15-1 + CH4""",
)

entry(
    index = 1978,
    label = "NC7H16 + CH3 <=> C7H15-2 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(54100, 'cm^3/(mol*s)'), n=2.26, Ea=(7287, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> C7H15-2 + CH4""",
)

entry(
    index = 1979,
    label = "NC7H16 + CH3 <=> C7H15-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(54100, 'cm^3/(mol*s)'), n=2.26, Ea=(7287, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> C7H15-3 + CH4""",
)

entry(
    index = 1980,
    label = "NC7H16 + CH3 <=> C7H15-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27050, 'cm^3/(mol*s)'), n=2.26, Ea=(7287, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3 <=> C7H15-4 + CH4""",
)

entry(
    index = 1981,
    label = "NC7H16 + O2 <=> C7H15-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.2e+13, 'cm^3/(mol*s)'), n=0, Ea=(52800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> C7H15-1 + HO2""",
)

entry(
    index = 1982,
    label = "NC7H16 + O2 <=> C7H15-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(50150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> C7H15-2 + HO2""",
)

entry(
    index = 1983,
    label = "NC7H16 + O2 <=> C7H15-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+13, 'cm^3/(mol*s)'), n=0, Ea=(50150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> C7H15-3 + HO2""",
)

entry(
    index = 1984,
    label = "NC7H16 + O2 <=> C7H15-4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.4e+13, 'cm^3/(mol*s)'), n=0, Ea=(50150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2 <=> C7H15-4 + HO2""",
)

entry(
    index = 1985,
    label = "NC7H16 + C2H5 <=> C7H15-1 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(13400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H5 <=> C7H15-1 + C2H6""",
)

entry(
    index = 1986,
    label = "NC7H16 + C2H5 <=> C7H15-2 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H5 <=> C7H15-2 + C2H6""",
)

entry(
    index = 1987,
    label = "NC7H16 + C2H5 <=> C7H15-3 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H5 <=> C7H15-3 + C2H6""",
)

entry(
    index = 1988,
    label = "NC7H16 + C2H5 <=> C7H15-4 + C2H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H5 <=> C7H15-4 + C2H6""",
)

entry(
    index = 1989,
    label = "NC7H16 + CH3O <=> C7H15-1 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.16e+11, 'cm^3/(mol*s)'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O <=> C7H15-1 + CH3OH""",
)

entry(
    index = 1990,
    label = "NC7H16 + CH3O <=> C7H15-2 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.19e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O <=> C7H15-2 + CH3OH""",
)

entry(
    index = 1991,
    label = "NC7H16 + CH3O <=> C7H15-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.19e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O <=> C7H15-3 + CH3OH""",
)

entry(
    index = 1992,
    label = "NC7H16 + CH3O <=> C7H15-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.095e+11, 'cm^3/(mol*s)'), n=0, Ea=(5000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O <=> C7H15-4 + CH3OH""",
)

entry(
    index = 1993,
    label = "NC7H16 + C2H3 <=> C7H15-1 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(18000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H3 <=> C7H15-1 + C2H4""",
)

entry(
    index = 1994,
    label = "NC7H16 + C2H3 <=> C7H15-2 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H3 <=> C7H15-2 + C2H4""",
)

entry(
    index = 1995,
    label = "NC7H16 + C2H3 <=> C7H15-3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H3 <=> C7H15-3 + C2H4""",
)

entry(
    index = 1996,
    label = "NC7H16 + C2H3 <=> C7H15-4 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(16800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C2H3 <=> C7H15-4 + C2H4""",
)

entry(
    index = 1997,
    label = "NC7H16 + CH3O2 <=> C7H15-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.386, 'cm^3/(mol*s)'), n=3.97, Ea=(18280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O2 <=> C7H15-1 + CH3O2H""",
)

entry(
    index = 1998,
    label = "NC7H16 + CH3O2 <=> C7H15-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20.37, 'cm^3/(mol*s)'), n=3.58, Ea=(14810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O2 <=> C7H15-2 + CH3O2H""",
)

entry(
    index = 1999,
    label = "NC7H16 + CH3O2 <=> C7H15-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(20.37, 'cm^3/(mol*s)'), n=3.58, Ea=(14810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O2 <=> C7H15-3 + CH3O2H""",
)

entry(
    index = 2000,
    label = "NC7H16 + CH3O2 <=> C7H15-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(10.19, 'cm^3/(mol*s)'), n=3.58, Ea=(14810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + CH3O2 <=> C7H15-4 + CH3O2H""",
)

entry(
    index = 2001,
    label = "NC7H16 + O2CHO <=> C7H15-1 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.68e+13, 'cm^3/(mol*s)'), n=0, Ea=(20440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2CHO <=> C7H15-1 + HO2CHO""",
)

entry(
    index = 2002,
    label = "NC7H16 + O2CHO <=> C7H15-2 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2CHO <=> C7H15-2 + HO2CHO""",
)

entry(
    index = 2003,
    label = "NC7H16 + O2CHO <=> C7H15-3 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.12e+13, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2CHO <=> C7H15-3 + HO2CHO""",
)

entry(
    index = 2004,
    label = "NC7H16 + O2CHO <=> C7H15-4 + HO2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(17690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + O2CHO <=> C7H15-4 + HO2CHO""",
)

entry(
    index = 2005,
    label = "NC7H16 + C7H15O2-1 <=> C7H15-1 + C7H15O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-1 <=> C7H15-1 + C7H15O2H-1""",
)

entry(
    index = 2006,
    label = "NC7H16 + C7H15O2-2 <=> C7H15-1 + C7H15O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-2 <=> C7H15-1 + C7H15O2H-2""",
)

entry(
    index = 2007,
    label = "NC7H16 + C7H15O2-3 <=> C7H15-1 + C7H15O2H-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-3 <=> C7H15-1 + C7H15O2H-3""",
)

entry(
    index = 2008,
    label = "NC7H16 + C7H15O2-4 <=> C7H15-1 + C7H15O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.21e+13, 'cm^3/(mol*s)'), n=0, Ea=(20430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-4 <=> C7H15-1 + C7H15O2H-4""",
)

entry(
    index = 2009,
    label = "NC7H16 + C7H15O2-1 <=> C7H15-2 + C7H15O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-1 <=> C7H15-2 + C7H15O2H-1""",
)

entry(
    index = 2010,
    label = "NC7H16 + C7H15O2-2 <=> C7H15-2 + C7H15O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-2 <=> C7H15-2 + C7H15O2H-2""",
)

entry(
    index = 2011,
    label = "NC7H16 + C7H15O2-3 <=> C7H15-2 + C7H15O2H-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-3 <=> C7H15-2 + C7H15O2H-3""",
)

entry(
    index = 2012,
    label = "NC7H16 + C7H15O2-4 <=> C7H15-2 + C7H15O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-4 <=> C7H15-2 + C7H15O2H-4""",
)

entry(
    index = 2013,
    label = "NC7H16 + C7H15O2-1 <=> C7H15-3 + C7H15O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-1 <=> C7H15-3 + C7H15O2H-1""",
)

entry(
    index = 2014,
    label = "NC7H16 + C7H15O2-2 <=> C7H15-3 + C7H15O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-2 <=> C7H15-3 + C7H15O2H-2""",
)

entry(
    index = 2015,
    label = "NC7H16 + C7H15O2-3 <=> C7H15-3 + C7H15O2H-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-3 <=> C7H15-3 + C7H15O2H-3""",
)

entry(
    index = 2016,
    label = "NC7H16 + C7H15O2-4 <=> C7H15-3 + C7H15O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (8.064e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-4 <=> C7H15-3 + C7H15O2H-4""",
)

entry(
    index = 2017,
    label = "NC7H16 + C7H15O2-1 <=> C7H15-4 + C7H15O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-1 <=> C7H15-4 + C7H15O2H-1""",
)

entry(
    index = 2018,
    label = "NC7H16 + C7H15O2-2 <=> C7H15-4 + C7H15O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-2 <=> C7H15-4 + C7H15O2H-2""",
)

entry(
    index = 2019,
    label = "NC7H16 + C7H15O2-3 <=> C7H15-4 + C7H15O2H-3",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-3 <=> C7H15-4 + C7H15O2H-3""",
)

entry(
    index = 2020,
    label = "NC7H16 + C7H15O2-4 <=> C7H15-4 + C7H15O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.032e+12, 'cm^3/(mol*s)'),
        n = 0,
        Ea = (17700, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15O2-4 <=> C7H15-4 + C7H15O2H-4""",
)

entry(
    index = 2021,
    label = "NC7H16 + C7H15-1 <=> C7H15-2 + NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15-1 <=> C7H15-2 + NC7H16""",
)

entry(
    index = 2022,
    label = "NC7H16 + C7H15-1 <=> C7H15-3 + NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15-1 <=> C7H15-3 + NC7H16""",
)

entry(
    index = 2023,
    label = "NC7H16 + C7H15-1 <=> C7H15-4 + NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15-1 <=> C7H15-4 + NC7H16""",
)

entry(
    index = 2024,
    label = "NC7H16 + C7H15-2 <=> C7H15-3 + NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15-2 <=> C7H15-3 + NC7H16""",
)

entry(
    index = 2025,
    label = "NC7H16 + C7H15-2 <=> C7H15-4 + NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15-2 <=> C7H15-4 + NC7H16""",
)

entry(
    index = 2026,
    label = "NC7H16 + C7H15-3 <=> C7H15-4 + NC7H16",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 'cm^3/(mol*s)'), n=0, Ea=(10400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7H16 + C7H15-3 <=> C7H15-4 + NC7H16""",
)

entry(
    index = 2027,
    label = "C7H15-1 <=> C5H11-1 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.227e+19, 's^-1'), n=-1.91, Ea=(31400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 <=> C5H11-1 + C2H4""",
)

entry(
    index = 2028,
    label = "C7H15-1 <=> C7H14-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.648e+13, 's^-1'), n=-0.26, Ea=(36010, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 <=> C7H14-1 + H""",
)

entry(
    index = 2029,
    label = "C7H15-2 <=> PC4H9 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.764e+18, 's^-1'), n=-1.79, Ea=(31360, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 <=> PC4H9 + C3H6""",
)

entry(
    index = 2030,
    label = "C7H15-2 <=> C7H14-1 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.067e+12, 's^-1'), n=0.09, Ea=(36810, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 <=> C7H14-1 + H""",
)

entry(
    index = 2031,
    label = "C7H15-2 <=> C7H14-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.229e+13, 's^-1'), n=-0.08, Ea=(35640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 <=> C7H14-2 + H""",
)

entry(
    index = 2032,
    label = "C7H15-3 <=> C4H8-1 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.165e+18, 's^-1'), n=-1.71, Ea=(30960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 <=> C4H8-1 + NC3H7""",
)

entry(
    index = 2033,
    label = "C7H15-3 <=> C6H12-1 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.698e+17, 's^-1'), n=-1.35, Ea=(31480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 <=> C6H12-1 + CH3""",
)

entry(
    index = 2034,
    label = "C7H15-3 <=> C7H14-2 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.229e+13, 's^-1'), n=-0.08, Ea=(35640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 <=> C7H14-2 + H""",
)

entry(
    index = 2035,
    label = "C7H15-3 <=> C7H14-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.152e+12, 's^-1'), n=-0.02, Ea=(35730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 <=> C7H14-3 + H""",
)

entry(
    index = 2036,
    label = "C7H15-4 <=> C2H5 + C5H10-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.143e+18, 's^-1'), n=-1.34, Ea=(31430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 <=> C2H5 + C5H10-1""",
)

entry(
    index = 2037,
    label = "C7H15-4 <=> C7H14-3 + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.819e+13, 's^-1'), n=-0.02, Ea=(35730, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 <=> C7H14-3 + H""",
)

entry(
    index = 2038,
    label = "C7H15-1 + O2 <=> C7H14-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-09, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + O2 <=> C7H14-1 + HO2""",
)

entry(
    index = 2039,
    label = "C7H15-2 + O2 <=> C7H14-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.5e-09, 'cm^3/(mol*s)'), n=0, Ea=(5020, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + O2 <=> C7H14-1 + HO2""",
)

entry(
    index = 2040,
    label = "C7H15-2 + O2 <=> C7H14-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-09, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + O2 <=> C7H14-2 + HO2""",
)

entry(
    index = 2041,
    label = "C7H15-3 + O2 <=> C7H14-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-09, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + O2 <=> C7H14-2 + HO2""",
)

entry(
    index = 2042,
    label = "C7H15-3 + O2 <=> C7H14-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e-09, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + O2 <=> C7H14-3 + HO2""",
)

entry(
    index = 2043,
    label = "C7H15-4 + O2 <=> C7H14-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e-09, 'cm^3/(mol*s)'), n=0, Ea=(3000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + O2 <=> C7H14-3 + HO2""",
)

entry(
    index = 2044,
    label = "C7H15-1 <=> C7H15-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.386e+09, 's^-1'), n=0.98, Ea=(33760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 <=> C7H15-3""",
)

entry(
    index = 2045,
    label = "C7H15-1 <=> C7H15-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.541e+09, 's^-1'), n=0.35, Ea=(19760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 <=> C7H15-4""",
)

entry(
    index = 2046,
    label = "C7H15-2 <=> C7H15-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.587e+08, 's^-1'), n=1.39, Ea=(39700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 <=> C7H15-3""",
)

entry(
    index = 2047,
    label = "C7H15-1 <=> C7H15-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.478e+08, 's^-1'), n=1.62, Ea=(38760, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 <=> C7H15-2""",
)

entry(
    index = 2048,
    label = "C7H14-1 + H <=> C7H131-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + H <=> C7H131-3 + H2""",
)

entry(
    index = 2049,
    label = "C7H14-1 + H <=> C7H131-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + H <=> C7H131-4 + H2""",
)

entry(
    index = 2050,
    label = "C7H14-1 + H <=> C7H131-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + H <=> C7H131-5 + H2""",
)

entry(
    index = 2051,
    label = "C7H14-1 + H <=> C7H131-6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + H <=> C7H131-6 + H2""",
)

entry(
    index = 2052,
    label = "C7H14-1 + H <=> C7H131-7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665000, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + H <=> C7H131-7 + H2""",
)

entry(
    index = 2053,
    label = "C7H14-1 + OH <=> C7H131-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH <=> C7H131-3 + H2O""",
)

entry(
    index = 2054,
    label = "C7H14-1 + OH <=> C7H131-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH <=> C7H131-4 + H2O""",
)

entry(
    index = 2055,
    label = "C7H14-1 + OH <=> C7H131-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH <=> C7H131-5 + H2O""",
)

entry(
    index = 2056,
    label = "C7H14-1 + OH <=> C7H131-6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH <=> C7H131-6 + H2O""",
)

entry(
    index = 2057,
    label = "C7H14-1 + OH <=> C7H131-7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH <=> C7H131-7 + H2O""",
)

entry(
    index = 2058,
    label = "C7H14-1 + CH3 <=> C7H131-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3 <=> C7H131-3 + CH4""",
)

entry(
    index = 2059,
    label = "C7H14-1 + CH3 <=> C7H131-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3 <=> C7H131-4 + CH4""",
)

entry(
    index = 2060,
    label = "C7H14-1 + CH3 <=> C7H131-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3 <=> C7H131-5 + CH4""",
)

entry(
    index = 2061,
    label = "C7H14-1 + CH3 <=> C7H131-6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3 <=> C7H131-6 + CH4""",
)

entry(
    index = 2062,
    label = "C7H14-1 + CH3 <=> C7H131-7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4521, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3 <=> C7H131-7 + CH4""",
)

entry(
    index = 2063,
    label = "C7H14-1 + HO2 <=> C7H131-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H131-3 + H2O2""",
)

entry(
    index = 2064,
    label = "C7H14-1 + HO2 <=> C7H131-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H131-4 + H2O2""",
)

entry(
    index = 2065,
    label = "C7H14-1 + HO2 <=> C7H131-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H131-5 + H2O2""",
)

entry(
    index = 2066,
    label = "C7H14-1 + HO2 <=> C7H131-6 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H131-6 + H2O2""",
)

entry(
    index = 2067,
    label = "C7H14-1 + HO2 <=> C7H131-7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H131-7 + H2O2""",
)

entry(
    index = 2068,
    label = "C7H14-1 + CH3O2 <=> C7H131-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O2 <=> C7H131-3 + CH3O2H""",
)

entry(
    index = 2069,
    label = "C7H14-1 + CH3O2 <=> C7H131-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O2 <=> C7H131-4 + CH3O2H""",
)

entry(
    index = 2070,
    label = "C7H14-1 + CH3O2 <=> C7H131-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O2 <=> C7H131-5 + CH3O2H""",
)

entry(
    index = 2071,
    label = "C7H14-1 + CH3O2 <=> C7H131-6 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O2 <=> C7H131-6 + CH3O2H""",
)

entry(
    index = 2072,
    label = "C7H14-1 + CH3O2 <=> C7H131-7 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O2 <=> C7H131-7 + CH3O2H""",
)

entry(
    index = 2073,
    label = "C7H14-1 + CH3O <=> C7H131-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O <=> C7H131-3 + CH3OH""",
)

entry(
    index = 2074,
    label = "C7H14-1 + CH3O <=> C7H131-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O <=> C7H131-4 + CH3OH""",
)

entry(
    index = 2075,
    label = "C7H14-1 + CH3O <=> C7H131-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O <=> C7H131-5 + CH3OH""",
)

entry(
    index = 2076,
    label = "C7H14-1 + CH3O <=> C7H131-6 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O <=> C7H131-6 + CH3OH""",
)

entry(
    index = 2077,
    label = "C7H14-1 + CH3O <=> C7H131-7 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + CH3O <=> C7H131-7 + CH3OH""",
)

entry(
    index = 2078,
    label = "C7H14-2 + H <=> C7H131-3 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(173000, 'cm^3/(mol*s)'), n=2.5, Ea=(2492, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C7H131-3 + H2""",
)

entry(
    index = 2079,
    label = "C7H14-2 + H <=> C7H132-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C7H132-4 + H2""",
)

entry(
    index = 2080,
    label = "C7H14-2 + H <=> C7H132-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C7H132-5 + H2""",
)

entry(
    index = 2081,
    label = "C7H14-2 + H <=> C7H132-6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C7H132-6 + H2""",
)

entry(
    index = 2082,
    label = "C7H14-2 + H <=> C7H132-7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665100, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + H <=> C7H132-7 + H2""",
)

entry(
    index = 2083,
    label = "C7H14-2 + OH <=> C7H131-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH <=> C7H131-3 + H2O""",
)

entry(
    index = 2084,
    label = "C7H14-2 + OH <=> C7H132-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH <=> C7H132-4 + H2O""",
)

entry(
    index = 2085,
    label = "C7H14-2 + OH <=> C7H132-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(467000, 'cm^3/(mol*s)'), n=1.61, Ea=(-35, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH <=> C7H132-5 + H2O""",
)

entry(
    index = 2086,
    label = "C7H14-2 + OH <=> C7H132-6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(467000, 'cm^3/(mol*s)'), n=1.61, Ea=(-35, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH <=> C7H132-6 + H2O""",
)

entry(
    index = 2087,
    label = "C7H14-2 + OH <=> C7H132-7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH <=> C7H132-7 + H2O""",
)

entry(
    index = 2088,
    label = "C7H14-2 + CH3 <=> C7H131-3 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.21, 'cm^3/(mol*s)'), n=3.5, Ea=(5675, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3 <=> C7H131-3 + CH4""",
)

entry(
    index = 2089,
    label = "C7H14-2 + CH3 <=> C7H132-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3 <=> C7H132-4 + CH4""",
)

entry(
    index = 2090,
    label = "C7H14-2 + CH3 <=> C7H132-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3 <=> C7H132-5 + CH4""",
)

entry(
    index = 2091,
    label = "C7H14-2 + CH3 <=> C7H132-6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3 <=> C7H132-6 + CH4""",
)

entry(
    index = 2092,
    label = "C7H14-2 + CH3 <=> C7H132-7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.4521, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3 <=> C7H132-7 + CH4""",
)

entry(
    index = 2093,
    label = "C7H14-2 + HO2 <=> C7H131-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9639, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H131-3 + H2O2""",
)

entry(
    index = 2094,
    label = "C7H14-2 + HO2 <=> C7H132-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H132-4 + H2O2""",
)

entry(
    index = 2095,
    label = "C7H14-2 + HO2 <=> C7H132-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H132-5 + H2O2""",
)

entry(
    index = 2096,
    label = "C7H14-2 + HO2 <=> C7H132-6 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H132-6 + H2O2""",
)

entry(
    index = 2097,
    label = "C7H14-2 + HO2 <=> C7H132-7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H132-7 + H2O2""",
)

entry(
    index = 2098,
    label = "C7H14-2 + CH3O2 <=> C7H131-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9639, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O2 <=> C7H131-3 + CH3O2H""",
)

entry(
    index = 2099,
    label = "C7H14-2 + CH3O2 <=> C7H132-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O2 <=> C7H132-4 + CH3O2H""",
)

entry(
    index = 2100,
    label = "C7H14-2 + CH3O2 <=> C7H132-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O2 <=> C7H132-5 + CH3O2H""",
)

entry(
    index = 2101,
    label = "C7H14-2 + CH3O2 <=> C7H132-6 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O2 <=> C7H132-6 + CH3O2H""",
)

entry(
    index = 2102,
    label = "C7H14-2 + CH3O2 <=> C7H132-7 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O2 <=> C7H132-7 + CH3O2H""",
)

entry(
    index = 2103,
    label = "C7H14-2 + CH3O <=> C7H131-3 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(90, 'cm^3/(mol*s)'), n=2.95, Ea=(11990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O <=> C7H131-3 + CH3OH""",
)

entry(
    index = 2104,
    label = "C7H14-2 + CH3O <=> C7H132-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O <=> C7H132-4 + CH3OH""",
)

entry(
    index = 2105,
    label = "C7H14-2 + CH3O <=> C7H132-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O <=> C7H132-5 + CH3OH""",
)

entry(
    index = 2106,
    label = "C7H14-2 + CH3O <=> C7H132-6 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O <=> C7H132-6 + CH3OH""",
)

entry(
    index = 2107,
    label = "C7H14-2 + CH3O <=> C7H132-7 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.17e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + CH3O <=> C7H132-7 + CH3OH""",
)

entry(
    index = 2108,
    label = "C7H14-3 + H <=> C7H133-1 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665100, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + H <=> C7H133-1 + H2""",
)

entry(
    index = 2109,
    label = "C7H14-3 + H <=> C7H132-4 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + H <=> C7H132-4 + H2""",
)

entry(
    index = 2110,
    label = "C7H14-3 + H <=> C7H133-5 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + H <=> C7H133-5 + H2""",
)

entry(
    index = 2111,
    label = "C7H14-3 + H <=> C7H133-6 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.3e+06, 'cm^3/(mol*s)'), n=2.4, Ea=(4471, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + H <=> C7H133-6 + H2""",
)

entry(
    index = 2112,
    label = "C7H14-3 + H <=> C7H133-7 + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(665100, 'cm^3/(mol*s)'), n=2.54, Ea=(6756, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + H <=> C7H133-7 + H2""",
)

entry(
    index = 2113,
    label = "C7H14-3 + OH <=> C7H133-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + OH <=> C7H133-1 + H2O""",
)

entry(
    index = 2114,
    label = "C7H14-3 + OH <=> C7H132-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + OH <=> C7H132-4 + H2O""",
)

entry(
    index = 2115,
    label = "C7H14-3 + OH <=> C7H133-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + OH <=> C7H133-5 + H2O""",
)

entry(
    index = 2116,
    label = "C7H14-3 + OH <=> C7H133-6 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + OH <=> C7H133-6 + H2O""",
)

entry(
    index = 2117,
    label = "C7H14-3 + OH <=> C7H133-7 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.054e+10, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + OH <=> C7H133-7 + H2O""",
)

entry(
    index = 2118,
    label = "C7H14-3 + CH3 <=> C7H133-1 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.9042, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3 <=> C7H133-1 + CH4""",
)

entry(
    index = 2119,
    label = "C7H14-3 + CH3 <=> C7H132-4 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3 <=> C7H132-4 + CH4""",
)

entry(
    index = 2120,
    label = "C7H14-3 + CH3 <=> C7H133-5 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.69, 'cm^3/(mol*s)'), n=3.31, Ea=(4002, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3 <=> C7H133-5 + CH4""",
)

entry(
    index = 2121,
    label = "C7H14-3 + CH3 <=> C7H133-6 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.51, 'cm^3/(mol*s)'), n=3.46, Ea=(5481, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3 <=> C7H133-6 + CH4""",
)

entry(
    index = 2122,
    label = "C7H14-3 + CH3 <=> C7H133-7 + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(0.9042, 'cm^3/(mol*s)'), n=3.65, Ea=(7154, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3 <=> C7H133-7 + CH4""",
)

entry(
    index = 2123,
    label = "C7H14-3 + HO2 <=> C7H133-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H133-1 + H2O2""",
)

entry(
    index = 2124,
    label = "C7H14-3 + HO2 <=> C7H132-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H132-4 + H2O2""",
)

entry(
    index = 2125,
    label = "C7H14-3 + HO2 <=> C7H133-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H133-5 + H2O2""",
)

entry(
    index = 2126,
    label = "C7H14-3 + HO2 <=> C7H133-6 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H133-6 + H2O2""",
)

entry(
    index = 2127,
    label = "C7H14-3 + HO2 <=> C7H133-7 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H133-7 + H2O2""",
)

entry(
    index = 2128,
    label = "C7H14-3 + CH3O2 <=> C7H133-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O2 <=> C7H133-1 + CH3O2H""",
)

entry(
    index = 2129,
    label = "C7H14-3 + CH3O2 <=> C7H132-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O2 <=> C7H132-4 + CH3O2H""",
)

entry(
    index = 2130,
    label = "C7H14-3 + CH3O2 <=> C7H133-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O2 <=> C7H133-5 + CH3O2H""",
)

entry(
    index = 2131,
    label = "C7H14-3 + CH3O2 <=> C7H133-6 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O2 <=> C7H133-6 + CH3O2H""",
)

entry(
    index = 2132,
    label = "C7H14-3 + CH3O2 <=> C7H133-7 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(47600, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O2 <=> C7H133-7 + CH3O2H""",
)

entry(
    index = 2133,
    label = "C7H14-3 + CH3O <=> C7H133-1 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.34e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O <=> C7H133-1 + CH3OH""",
)

entry(
    index = 2134,
    label = "C7H14-3 + CH3O <=> C7H132-4 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O <=> C7H132-4 + CH3OH""",
)

entry(
    index = 2135,
    label = "C7H14-3 + CH3O <=> C7H133-5 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(40, 'cm^3/(mol*s)'), n=2.9, Ea=(8609, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O <=> C7H133-5 + CH3OH""",
)

entry(
    index = 2136,
    label = "C7H14-3 + CH3O <=> C7H133-6 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(4571, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O <=> C7H133-6 + CH3OH""",
)

entry(
    index = 2137,
    label = "C7H14-3 + CH3O <=> C7H133-7 + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.34e+11, 'cm^3/(mol*s)'), n=0, Ea=(6458, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + CH3O <=> C7H133-7 + CH3OH""",
)

entry(
    index = 2138,
    label = "C7H131-3 + HO2 <=> C7H13O1-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-3 + HO2 <=> C7H13O1-3 + OH""",
)

entry(
    index = 2139,
    label = "C7H131-3 + CH3O2 <=> C7H13O1-3 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-3 + CH3O2 <=> C7H13O1-3 + CH3O""",
)

entry(
    index = 2140,
    label = "C7H131-3 + C2H5O2 <=> C7H13O1-3 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-3 + C2H5O2 <=> C7H13O1-3 + C2H5O""",
)

entry(
    index = 2141,
    label = "C7H131-3 <=> C4H6 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.105e+19, 's^-1'), n=-1.53, Ea=(40700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-3 <=> C4H6 + NC3H7""",
)

entry(
    index = 2142,
    label = "C7H131-4 <=> C5H10-1 + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.492e+12, 's^-1'), n=0.03, Ea=(37300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-4 <=> C5H10-1 + C2H3""",
)

entry(
    index = 2143,
    label = "C7H131-4 <=> C5H81-3 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.289e+14, 's^-1'), n=-0.68, Ea=(22050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-4 <=> C5H81-3 + C2H5""",
)

entry(
    index = 2144,
    label = "C7H131-5 <=> C4H8-1 + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.177e+15, 's^-1'), n=-1.18, Ea=(17980, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-5 <=> C4H8-1 + C3H5-A""",
)

entry(
    index = 2145,
    label = "C7H131-6 <=> C3H6 + C4H71-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.238e+17, 's^-1'), n=-1.17, Ea=(30740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-6 <=> C3H6 + C4H71-4""",
)

entry(
    index = 2146,
    label = "C7H131-7 <=> C2H4 + C5H91-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.416e+17, 's^-1'), n=-1.56, Ea=(31180, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H131-7 <=> C2H4 + C5H91-5""",
)

entry(
    index = 2147,
    label = "C7H132-4 + HO2 <=> C7H13O2-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H132-4 + HO2 <=> C7H13O2-4 + OH""",
)

entry(
    index = 2148,
    label = "C7H132-4 + CH3O2 <=> C7H13O2-4 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H132-4 + CH3O2 <=> C7H13O2-4 + CH3O""",
)

entry(
    index = 2149,
    label = "C7H132-4 + C2H5O2 <=> C7H13O2-4 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H132-4 + C2H5O2 <=> C7H13O2-4 + C2H5O""",
)

entry(
    index = 2150,
    label = "C7H132-5 <=> C4H8-1 + C3H5-S",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.419e+17, 's^-1'), n=-1.49, Ea=(44260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H132-5 <=> C4H8-1 + C3H5-S""",
)

entry(
    index = 2151,
    label = "C7H132-6 <=> C3H6 + C4H71-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.189e+14, 's^-1'), n=-0.74, Ea=(17740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H132-6 <=> C3H6 + C4H71-3""",
)

entry(
    index = 2152,
    label = "C7H132-7 <=> C2H4 + C5H92-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.786e+17, 's^-1'), n=-1.57, Ea=(31160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H132-7 <=> C2H4 + C5H92-5""",
)

entry(
    index = 2153,
    label = "C7H133-1 <=> C2H4 + C5H91-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.007e+17, 's^-1'), n=-1.32, Ea=(43990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H133-1 <=> C2H4 + C5H91-1""",
)

entry(
    index = 2154,
    label = "C7H133-5 + HO2 <=> C7H13O3-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H133-5 + HO2 <=> C7H13O3-5 + OH""",
)

entry(
    index = 2155,
    label = "C7H133-5 + CH3O2 <=> C7H13O3-5 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H133-5 + CH3O2 <=> C7H13O3-5 + CH3O""",
)

entry(
    index = 2156,
    label = "C7H133-5 + C2H5O2 <=> C7H13O3-5 + C2H5O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.64e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H133-5 + C2H5O2 <=> C7H13O3-5 + C2H5O""",
)

entry(
    index = 2157,
    label = "C7H133-6 <=> C3H6 + C4H71-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.354e+16, 's^-1'), n=-0.92, Ea=(43560, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H133-6 <=> C3H6 + C4H71-1""",
)

entry(
    index = 2158,
    label = "C7H133-7 <=> C2H4 + C5H91-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.274e+15, 's^-1'), n=-1.2, Ea=(18050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H133-7 <=> C2H4 + C5H91-3""",
)

entry(
    index = 2159,
    label = "C7H13O1-3 <=> C2H3CHO + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.625e+19, 's^-1'), n=-1.96, Ea=(10850, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H13O1-3 <=> C2H3CHO + PC4H9""",
)

entry(
    index = 2160,
    label = "C7H13O1-3 <=> NC4H9CHO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.026e+18, 's^-1'), n=-1.51, Ea=(23300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H13O1-3 <=> NC4H9CHO + C2H3""",
)

entry(
    index = 2161,
    label = "C7H13O2-4 <=> SC3H5CHO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.319e+19, 's^-1'), n=-1.93, Ea=(11130, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H13O2-4 <=> SC3H5CHO + NC3H7""",
)

entry(
    index = 2162,
    label = "C7H13O2-4 <=> NC3H7CHO + C3H5-S",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.012e+22, 's^-1'), n=-2.46, Ea=(29190, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H13O2-4 <=> NC3H7CHO + C3H5-S""",
)

entry(
    index = 2163,
    label = "C7H13O3-5 <=> C4H7CHO1-1 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.033e+18, 's^-1'), n=-1.62, Ea=(10450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H13O3-5 <=> C4H7CHO1-1 + C2H5""",
)

entry(
    index = 2164,
    label = "C7H13O3-5 <=> C2H5CHO + C4H71-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.425e+21, 's^-1'), n=-2.43, Ea=(30090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H13O3-5 <=> C2H5CHO + C4H71-1""",
)

entry(
    index = 2165,
    label = "C7H14-1 + OH => CH2O + C6H13-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH => CH2O + C6H13-1""",
)

entry(
    index = 2166,
    label = "C7H14-1 + OH => CH3CHO + C5H11-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + OH => CH3CHO + C5H11-1""",
)

entry(
    index = 2167,
    label = "C7H14-2 + OH => CH3CHO + C5H11-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH => CH3CHO + C5H11-1""",
)

entry(
    index = 2168,
    label = "C7H14-2 + OH => C2H5CHO + PC4H9",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + OH => C2H5CHO + PC4H9""",
)

entry(
    index = 2169,
    label = "C7H14-3 + OH => C2H5CHO + PC4H9",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-4000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + OH => C2H5CHO + PC4H9""",
)

entry(
    index = 2170,
    label = "C7H14-1 + O => CH2CHO + C5H11-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + O => CH2CHO + C5H11-1""",
)

entry(
    index = 2171,
    label = "C7H14-2 + O => CH3CHO + C5H10-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + O => CH3CHO + C5H10-1""",
)

entry(
    index = 2172,
    label = "C7H14-3 + O => CH3CHO + C5H10-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + O => CH3CHO + C5H10-1""",
)

entry(
    index = 2173,
    label = "C7H14-1 <=> PC4H9 + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.17e+21, 's^-1'), n=-1.62, Ea=(75330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 <=> PC4H9 + C3H5-A""",
)

entry(
    index = 2174,
    label = "C7H14-2 <=> C4H71-3 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.735e+21, 's^-1'), n=-1.74, Ea=(75710, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 <=> C4H71-3 + NC3H7""",
)

entry(
    index = 2175,
    label = "C7H14-3 <=> C5H91-3 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.947e+21, 's^-1'), n=-1.85, Ea=(75790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 <=> C5H91-3 + C2H5""",
)

entry(
    index = 2176,
    label = "C7H15O2-1 <=> C7H15-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.657e+20, 's^-1'), n=-1.67, Ea=(35400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 <=> C7H15-1 + O2""",
)

entry(
    index = 2177,
    label = "C7H15O2-2 <=> C7H15-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.357e+23, 's^-1'), n=-2.36, Ea=(37670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H15-2 + O2""",
)

entry(
    index = 2178,
    label = "C7H15O2-3 <=> C7H15-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.357e+23, 's^-1'), n=-2.36, Ea=(37670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H15-3 + O2""",
)

entry(
    index = 2179,
    label = "C7H15O2-4 <=> C7H15-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.357e+23, 's^-1'), n=-2.36, Ea=(37670, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 <=> C7H15-4 + O2""",
)

entry(
    index = 2180,
    label = "C7H15-1 + C7H15O2-1 <=> C7H15O-1 + C7H15O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + C7H15O2-1 <=> C7H15O-1 + C7H15O-1""",
)

entry(
    index = 2181,
    label = "C7H15-1 + C7H15O2-2 <=> C7H15O-1 + C7H15O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + C7H15O2-2 <=> C7H15O-1 + C7H15O-2""",
)

entry(
    index = 2182,
    label = "C7H15-1 + C7H15O2-3 <=> C7H15O-1 + C7H15O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + C7H15O2-3 <=> C7H15O-1 + C7H15O-3""",
)

entry(
    index = 2183,
    label = "C7H15-1 + C7H15O2-4 <=> C7H15O-1 + C7H15O-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + C7H15O2-4 <=> C7H15O-1 + C7H15O-4""",
)

entry(
    index = 2184,
    label = "C7H15-2 + C7H15O2-1 <=> C7H15O-2 + C7H15O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + C7H15O2-1 <=> C7H15O-2 + C7H15O-1""",
)

entry(
    index = 2185,
    label = "C7H15-2 + C7H15O2-2 <=> C7H15O-2 + C7H15O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + C7H15O2-2 <=> C7H15O-2 + C7H15O-2""",
)

entry(
    index = 2186,
    label = "C7H15-2 + C7H15O2-3 <=> C7H15O-2 + C7H15O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + C7H15O2-3 <=> C7H15O-2 + C7H15O-3""",
)

entry(
    index = 2187,
    label = "C7H15-2 + C7H15O2-4 <=> C7H15O-2 + C7H15O-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + C7H15O2-4 <=> C7H15O-2 + C7H15O-4""",
)

entry(
    index = 2188,
    label = "C7H15-3 + C7H15O2-1 <=> C7H15O-3 + C7H15O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + C7H15O2-1 <=> C7H15O-3 + C7H15O-1""",
)

entry(
    index = 2189,
    label = "C7H15-3 + C7H15O2-2 <=> C7H15O-3 + C7H15O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + C7H15O2-2 <=> C7H15O-3 + C7H15O-2""",
)

entry(
    index = 2190,
    label = "C7H15-3 + C7H15O2-3 <=> C7H15O-3 + C7H15O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + C7H15O2-3 <=> C7H15O-3 + C7H15O-3""",
)

entry(
    index = 2191,
    label = "C7H15-3 + C7H15O2-4 <=> C7H15O-3 + C7H15O-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + C7H15O2-4 <=> C7H15O-3 + C7H15O-4""",
)

entry(
    index = 2192,
    label = "C7H15-4 + C7H15O2-1 <=> C7H15O-4 + C7H15O-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + C7H15O2-1 <=> C7H15O-4 + C7H15O-1""",
)

entry(
    index = 2193,
    label = "C7H15-4 + C7H15O2-2 <=> C7H15O-4 + C7H15O-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + C7H15O2-2 <=> C7H15O-4 + C7H15O-2""",
)

entry(
    index = 2194,
    label = "C7H15-4 + C7H15O2-3 <=> C7H15O-4 + C7H15O-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + C7H15O2-3 <=> C7H15O-4 + C7H15O-3""",
)

entry(
    index = 2195,
    label = "C7H15-4 + C7H15O2-4 <=> C7H15O-4 + C7H15O-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + C7H15O2-4 <=> C7H15O-4 + C7H15O-4""",
)

entry(
    index = 2196,
    label = "C7H15-1 + HO2 <=> C7H15O-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + HO2 <=> C7H15O-1 + OH""",
)

entry(
    index = 2197,
    label = "C7H15-2 + HO2 <=> C7H15O-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + HO2 <=> C7H15O-2 + OH""",
)

entry(
    index = 2198,
    label = "C7H15-3 + HO2 <=> C7H15O-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + HO2 <=> C7H15O-3 + OH""",
)

entry(
    index = 2199,
    label = "C7H15-4 + HO2 <=> C7H15O-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + HO2 <=> C7H15O-4 + OH""",
)

entry(
    index = 2200,
    label = "C7H15-1 + CH3O2 <=> C7H15O-1 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-1 + CH3O2 <=> C7H15O-1 + CH3O""",
)

entry(
    index = 2201,
    label = "C7H15-2 + CH3O2 <=> C7H15O-2 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-2 + CH3O2 <=> C7H15O-2 + CH3O""",
)

entry(
    index = 2202,
    label = "C7H15-3 + CH3O2 <=> C7H15O-3 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-3 + CH3O2 <=> C7H15O-3 + CH3O""",
)

entry(
    index = 2203,
    label = "C7H15-4 + CH3O2 <=> C7H15O-4 + CH3O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7e+12, 'cm^3/(mol*s)'), n=0, Ea=(-1000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15-4 + CH3O2 <=> C7H15O-4 + CH3O""",
)

entry(
    index = 2204,
    label = "C7H15O2-1 <=> C7H14-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.0044e+39, 's^-1'), n=-8.11, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 <=> C7H14-1 + HO2""",
)

entry(
    index = 2205,
    label = "C7H15O2-2 <=> C7H14-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.0075e+43, 's^-1'), n=-9.41, Ea=(42490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14-1 + HO2""",
)

entry(
    index = 2206,
    label = "C7H15O2-2 <=> C7H14-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.0044e+39, 's^-1'), n=-8.11, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14-2 + HO2""",
)

entry(
    index = 2207,
    label = "C7H15O2-3 <=> C7H14-2 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.0044e+39, 's^-1'), n=-8.11, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14-2 + HO2""",
)

entry(
    index = 2208,
    label = "C7H15O2-3 <=> C7H14-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.0044e+39, 's^-1'), n=-8.11, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14-3 + HO2""",
)

entry(
    index = 2209,
    label = "C7H15O2-4 <=> C7H14-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.009e+39, 's^-1'), n=-8.11, Ea=(41490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 <=> C7H14-3 + HO2""",
)

entry(
    index = 2210,
    label = "C7H15O2-1 <=> C7H14OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 <=> C7H14OOH1-2""",
)

entry(
    index = 2211,
    label = "C7H15O2-1 <=> C7H14OOH1-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 <=> C7H14OOH1-3""",
)

entry(
    index = 2212,
    label = "C7H15O2-1 <=> C7H14OOH1-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(18650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 <=> C7H14OOH1-4""",
)

entry(
    index = 2213,
    label = "C7H15O2-1 <=> C7H14OOH1-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.906e+08, 's^-1'), n=0, Ea=(21650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 <=> C7H14OOH1-5""",
)

entry(
    index = 2214,
    label = "C7H15O2-2 <=> C7H14OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14OOH2-1""",
)

entry(
    index = 2215,
    label = "C7H15O2-2 <=> C7H14OOH2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14OOH2-3""",
)

entry(
    index = 2216,
    label = "C7H15O2-2 <=> C7H14OOH2-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14OOH2-4""",
)

entry(
    index = 2217,
    label = "C7H15O2-2 <=> C7H14OOH2-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(18650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14OOH2-5""",
)

entry(
    index = 2218,
    label = "C7H15O2-2 <=> C7H14OOH2-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.906e+08, 's^-1'), n=0, Ea=(21650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 <=> C7H14OOH2-6""",
)

entry(
    index = 2219,
    label = "C7H15O2-3 <=> C7H14OOH3-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+10, 's^-1'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14OOH3-1""",
)

entry(
    index = 2220,
    label = "C7H15O2-3 <=> C7H14OOH3-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14OOH3-2""",
)

entry(
    index = 2221,
    label = "C7H15O2-3 <=> C7H14OOH3-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14OOH3-4""",
)

entry(
    index = 2222,
    label = "C7H15O2-3 <=> C7H14OOH3-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14OOH3-5""",
)

entry(
    index = 2223,
    label = "C7H15O2-3 <=> C7H14OOH3-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(18650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14OOH3-6""",
)

entry(
    index = 2224,
    label = "C7H15O2-3 <=> C7H14OOH3-7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.86e+08, 's^-1'), n=0, Ea=(25150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 <=> C7H14OOH3-7""",
)

entry(
    index = 2225,
    label = "C7H15O2-4 <=> C7H14OOH4-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.376e+09, 's^-1'), n=0, Ea=(21950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 <=> C7H14OOH4-1""",
)

entry(
    index = 2226,
    label = "C7H15O2-4 <=> C7H14OOH4-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+10, 's^-1'), n=0, Ea=(20450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 <=> C7H14OOH4-2""",
)

entry(
    index = 2227,
    label = "C7H15O2-4 <=> C7H14OOH4-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 <=> C7H14OOH4-3""",
)

entry(
    index = 2228,
    label = "C7H15O2-1 + HO2 <=> C7H15O2H-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 + HO2 <=> C7H15O2H-1 + O2""",
)

entry(
    index = 2229,
    label = "C7H15O2-2 + HO2 <=> C7H15O2H-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 + HO2 <=> C7H15O2H-2 + O2""",
)

entry(
    index = 2230,
    label = "C7H15O2-3 + HO2 <=> C7H15O2H-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 + HO2 <=> C7H15O2H-3 + O2""",
)

entry(
    index = 2231,
    label = "C7H15O2-4 + HO2 <=> C7H15O2H-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.75e+10, 'cm^3/(mol*s)'), n=0, Ea=(-3275, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 + HO2 <=> C7H15O2H-4 + O2""",
)

entry(
    index = 2232,
    label = "H2O2 + C7H15O2-1 <=> HO2 + C7H15O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C7H15O2-1 <=> HO2 + C7H15O2H-1""",
)

entry(
    index = 2233,
    label = "H2O2 + C7H15O2-2 <=> HO2 + C7H15O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C7H15O2-2 <=> HO2 + C7H15O2H-2""",
)

entry(
    index = 2234,
    label = "H2O2 + C7H15O2-3 <=> HO2 + C7H15O2H-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C7H15O2-3 <=> HO2 + C7H15O2H-3""",
)

entry(
    index = 2235,
    label = "H2O2 + C7H15O2-4 <=> HO2 + C7H15O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is H2O2 + C7H15O2-4 <=> HO2 + C7H15O2H-4""",
)

entry(
    index = 2236,
    label = "C7H15O2-1 + CH3O2 => C7H15O-1 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 + CH3O2 => C7H15O-1 + CH3O + O2""",
)

entry(
    index = 2237,
    label = "C7H15O2-2 + CH3O2 => C7H15O-2 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 + CH3O2 => C7H15O-2 + CH3O + O2""",
)

entry(
    index = 2238,
    label = "C7H15O2-3 + CH3O2 => C7H15O-3 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 + CH3O2 => C7H15O-3 + CH3O + O2""",
)

entry(
    index = 2239,
    label = "C7H15O2-4 + CH3O2 => C7H15O-4 + CH3O + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 + CH3O2 => C7H15O-4 + CH3O + O2""",
)

entry(
    index = 2240,
    label = "C7H15O2-1 + C7H15O2-1 => O2 + C7H15O-1 + C7H15O-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 + C7H15O2-1 => O2 + C7H15O-1 + C7H15O-1""",
)

entry(
    index = 2241,
    label = "C7H15O2-1 + C7H15O2-2 => C7H15O-1 + C7H15O-2 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 + C7H15O2-2 => C7H15O-1 + C7H15O-2 + O2""",
)

entry(
    index = 2242,
    label = "C7H15O2-1 + C7H15O2-3 => C7H15O-1 + C7H15O-3 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 + C7H15O2-3 => C7H15O-1 + C7H15O-3 + O2""",
)

entry(
    index = 2243,
    label = "C7H15O2-1 + C7H15O2-4 => C7H15O-1 + C7H15O-4 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-1 + C7H15O2-4 => C7H15O-1 + C7H15O-4 + O2""",
)

entry(
    index = 2244,
    label = "C7H15O2-2 + C7H15O2-2 => O2 + C7H15O-2 + C7H15O-2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 + C7H15O2-2 => O2 + C7H15O-2 + C7H15O-2""",
)

entry(
    index = 2245,
    label = "C7H15O2-2 + C7H15O2-3 => C7H15O-2 + C7H15O-3 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 + C7H15O2-3 => C7H15O-2 + C7H15O-3 + O2""",
)

entry(
    index = 2246,
    label = "C7H15O2-2 + C7H15O2-4 => C7H15O-2 + C7H15O-4 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-2 + C7H15O2-4 => C7H15O-2 + C7H15O-4 + O2""",
)

entry(
    index = 2247,
    label = "C7H15O2-3 + C7H15O2-3 => O2 + C7H15O-3 + C7H15O-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 + C7H15O2-3 => O2 + C7H15O-3 + C7H15O-3""",
)

entry(
    index = 2248,
    label = "C7H15O2-3 + C7H15O2-4 => C7H15O-3 + C7H15O-4 + O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-3 + C7H15O2-4 => C7H15O-3 + C7H15O-4 + O2""",
)

entry(
    index = 2249,
    label = "C7H15O2-4 + C7H15O2-4 => O2 + C7H15O-4 + C7H15O-4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.4e+16, 'cm^3/(mol*s)'),
        n = -1.61,
        Ea = (1860, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C7H15O2-4 + C7H15O2-4 => O2 + C7H15O-4 + C7H15O-4""",
)

entry(
    index = 2250,
    label = "C7H15O2H-1 <=> C7H15O-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2H-1 <=> C7H15O-1 + OH""",
)

entry(
    index = 2251,
    label = "C7H15O2H-2 <=> C7H15O-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2H-2 <=> C7H15O-2 + OH""",
)

entry(
    index = 2252,
    label = "C7H15O2H-3 <=> C7H15O-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2H-3 <=> C7H15O-3 + OH""",
)

entry(
    index = 2253,
    label = "C7H15O2H-4 <=> C7H15O-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O2H-4 <=> C7H15O-4 + OH""",
)

entry(
    index = 2254,
    label = "C7H15O-1 <=> CH2O + C6H13-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.023e+20, 's^-1'), n=-2.18, Ea=(24830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O-1 <=> CH2O + C6H13-1""",
)

entry(
    index = 2255,
    label = "C7H15O-2 <=> CH3CHO + C5H11-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.475e+22, 's^-1'), n=-2.55, Ea=(20370, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O-2 <=> CH3CHO + C5H11-1""",
)

entry(
    index = 2256,
    label = "C7H15O-3 <=> C2H5CHO + PC4H9",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.573e+17, 's^-1'), n=-1.16, Ea=(19370, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O-3 <=> C2H5CHO + PC4H9""",
)

entry(
    index = 2257,
    label = "C7H15O-4 <=> NC3H7CHO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.769e+22, 's^-1'), n=-2.58, Ea=(21820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H15O-4 <=> NC3H7CHO + NC3H7""",
)

entry(
    index = 2258,
    label = "C7H14-1 + HO2 <=> C7H14OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2500, 'cm^3/(mol*s)'), n=2.5, Ea=(11200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H14OOH1-2""",
)

entry(
    index = 2259,
    label = "C7H14-1 + HO2 <=> C7H14OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=2.5, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-1 + HO2 <=> C7H14OOH2-1""",
)

entry(
    index = 2260,
    label = "C7H14-2 + HO2 <=> C7H14OOH2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=2.5, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H14OOH2-3""",
)

entry(
    index = 2261,
    label = "C7H14-2 + HO2 <=> C7H14OOH3-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=2.5, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-2 + HO2 <=> C7H14OOH3-2""",
)

entry(
    index = 2262,
    label = "C7H14-3 + HO2 <=> C7H14OOH3-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=2.5, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H14OOH3-4""",
)

entry(
    index = 2263,
    label = "C7H14-3 + HO2 <=> C7H14OOH4-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2700, 'cm^3/(mol*s)'), n=2.5, Ea=(10500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14-3 + HO2 <=> C7H14OOH4-3""",
)

entry(
    index = 2264,
    label = "C7H14OOH1-2 => C7H14O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-2 => C7H14O1-2 + OH""",
)

entry(
    index = 2265,
    label = "C7H14OOH1-3 => C7H14O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-3 => C7H14O1-3 + OH""",
)

entry(
    index = 2266,
    label = "C7H14OOH1-4 => C7H14O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-4 => C7H14O1-4 + OH""",
)

entry(
    index = 2267,
    label = "C7H14OOH1-5 => C7H14O1-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.172e+09, 's^-1'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-5 => C7H14O1-5 + OH""",
)

entry(
    index = 2268,
    label = "C7H14OOH2-1 => C7H14O1-2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-1 => C7H14O1-2 + OH""",
)

entry(
    index = 2269,
    label = "C7H14OOH2-3 => C7H14O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-3 => C7H14O2-3 + OH""",
)

entry(
    index = 2270,
    label = "C7H14OOH2-4 => C7H14O2-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-4 => C7H14O2-4 + OH""",
)

entry(
    index = 2271,
    label = "C7H14OOH2-5 => C7H14O2-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-5 => C7H14O2-5 + OH""",
)

entry(
    index = 2272,
    label = "C7H14OOH2-6 => C7H14O2-6 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.172e+09, 's^-1'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-6 => C7H14O2-6 + OH""",
)

entry(
    index = 2273,
    label = "C7H14OOH3-1 => C7H14O1-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-1 => C7H14O1-3 + OH""",
)

entry(
    index = 2274,
    label = "C7H14OOH3-2 => C7H14O2-3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-2 => C7H14O2-3 + OH""",
)

entry(
    index = 2275,
    label = "C7H14OOH3-4 => C7H14O3-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-4 => C7H14O3-4 + OH""",
)

entry(
    index = 2276,
    label = "C7H14OOH3-5 => C7H14O3-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-5 => C7H14O3-5 + OH""",
)

entry(
    index = 2277,
    label = "C7H14OOH3-6 => C7H14O2-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-6 => C7H14O2-5 + OH""",
)

entry(
    index = 2278,
    label = "C7H14OOH3-7 => C7H14O1-5 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.172e+09, 's^-1'), n=0, Ea=(1800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-7 => C7H14O1-5 + OH""",
)

entry(
    index = 2279,
    label = "C7H14OOH4-1 => C7H14O1-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-1 => C7H14O1-4 + OH""",
)

entry(
    index = 2280,
    label = "C7H14OOH4-2 => C7H14O2-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-2 => C7H14O2-4 + OH""",
)

entry(
    index = 2281,
    label = "C7H14OOH4-3 => C7H14O3-4 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-3 => C7H14O3-4 + OH""",
)

entry(
    index = 2282,
    label = "C7H14OOH1-3 => OH + CH2O + C6H12-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.147e+09, 's^-1'), n=1.23, Ea=(30370, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-3 => OH + CH2O + C6H12-1""",
)

entry(
    index = 2283,
    label = "C7H14OOH2-4 => OH + CH3CHO + C5H10-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.548e+12, 's^-1'), n=0.59, Ea=(30090, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-4 => OH + CH3CHO + C5H10-1""",
)

entry(
    index = 2284,
    label = "C7H14OOH3-1 => OH + NC4H9CHO + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(8.182e+13, 's^-1'), n=-0.13, Ea=(31330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-1 => OH + NC4H9CHO + C2H4""",
)

entry(
    index = 2285,
    label = "C7H14OOH3-5 => OH + C2H5CHO + C4H8-1",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.661e+13, 's^-1'), n=0.13, Ea=(30430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-5 => OH + C2H5CHO + C4H8-1""",
)

entry(
    index = 2286,
    label = "C7H14OOH4-2 => OH + NC3H7CHO + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6.186e+13, 's^-1'), n=0.09, Ea=(30840, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-2 => OH + NC3H7CHO + C3H6""",
)

entry(
    index = 2287,
    label = "C7H14OOH1-3 <=> C4H7OOH1-4 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.609e+12, 's^-1'), n=0.54, Ea=(27740, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-3 <=> C4H7OOH1-4 + NC3H7""",
)

entry(
    index = 2288,
    label = "C7H14OOH1-4 => C5H10-1 + C2H4 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.446e+11, 's^-1'), n=0.69, Ea=(30820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-4 => C5H10-1 + C2H4 + HO2""",
)

entry(
    index = 2289,
    label = "C7H14OOH1-4 <=> C5H9OOH1-5 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.721e+12, 's^-1'), n=0.51, Ea=(27900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-4 <=> C5H9OOH1-5 + C2H5""",
)

entry(
    index = 2290,
    label = "C7H14OOH2-4 <=> C5H9OOH1-4 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.08e+12, 's^-1'), n=0.32, Ea=(29230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-4 <=> C5H9OOH1-4 + C2H5""",
)

entry(
    index = 2291,
    label = "C7H14OOH2-5 <=> C4H8-1 + C3H6OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.856e+13, 's^-1'), n=0.03, Ea=(31380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-5 <=> C4H8-1 + C3H6OOH2-1""",
)

entry(
    index = 2292,
    label = "C7H14OOH2-5 <=> C6H11OOH1-5 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.417e+10, 's^-1'), n=0.75, Ea=(30260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-5 <=> C6H11OOH1-5 + CH3""",
)

entry(
    index = 2293,
    label = "C7H14OOH3-5 <=> C6H11OOH1-4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.417e+10, 's^-1'), n=0.75, Ea=(30260, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-5 <=> C6H11OOH1-4 + CH3""",
)

entry(
    index = 2294,
    label = "C7H14OOH3-6 <=> C4H8OOH2-1 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.578e+11, 's^-1'), n=0.6, Ea=(29170, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-6 <=> C4H8OOH2-1 + C3H6""",
)

entry(
    index = 2295,
    label = "C7H14OOH4-1 <=> C5H10OOH2-1 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.643e+12, 's^-1'), n=0.44, Ea=(29320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-1 <=> C5H10OOH2-1 + C2H4""",
)

entry(
    index = 2296,
    label = "C7H14OOH1-2O2 <=> C7H14OOH1-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.367e+23, 's^-1'), n=-2.37, Ea=(37640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-2O2 <=> C7H14OOH1-2 + O2""",
)

entry(
    index = 2297,
    label = "C7H14OOH1-3O2 <=> C7H14OOH1-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.367e+23, 's^-1'), n=-2.37, Ea=(37640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-3O2 <=> C7H14OOH1-3 + O2""",
)

entry(
    index = 2298,
    label = "C7H14OOH1-4O2 <=> C7H14OOH1-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.367e+23, 's^-1'), n=-2.37, Ea=(37640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-4O2 <=> C7H14OOH1-4 + O2""",
)

entry(
    index = 2299,
    label = "C7H14OOH1-5O2 <=> C7H14OOH1-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.367e+23, 's^-1'), n=-2.37, Ea=(37640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-5O2 <=> C7H14OOH1-5 + O2""",
)

entry(
    index = 2300,
    label = "C7H14OOH2-1O2 <=> C7H14OOH2-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.323e+20, 's^-1'), n=-1.65, Ea=(35280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-1O2 <=> C7H14OOH2-1 + O2""",
)

entry(
    index = 2301,
    label = "C7H14OOH2-3O2 <=> C7H14OOH2-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-3O2 <=> C7H14OOH2-3 + O2""",
)

entry(
    index = 2302,
    label = "C7H14OOH2-4O2 <=> C7H14OOH2-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-4O2 <=> C7H14OOH2-4 + O2""",
)

entry(
    index = 2303,
    label = "C7H14OOH2-5O2 <=> C7H14OOH2-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-5O2 <=> C7H14OOH2-5 + O2""",
)

entry(
    index = 2304,
    label = "C7H14OOH2-6O2 <=> C7H14OOH2-6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-6O2 <=> C7H14OOH2-6 + O2""",
)

entry(
    index = 2305,
    label = "C7H14OOH3-1O2 <=> C7H14OOH3-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.323e+20, 's^-1'), n=-1.65, Ea=(35280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-1O2 <=> C7H14OOH3-1 + O2""",
)

entry(
    index = 2306,
    label = "C7H14OOH3-2O2 <=> C7H14OOH3-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-2O2 <=> C7H14OOH3-2 + O2""",
)

entry(
    index = 2307,
    label = "C7H14OOH3-4O2 <=> C7H14OOH3-4 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-4O2 <=> C7H14OOH3-4 + O2""",
)

entry(
    index = 2308,
    label = "C7H14OOH3-5O2 <=> C7H14OOH3-5 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-5O2 <=> C7H14OOH3-5 + O2""",
)

entry(
    index = 2309,
    label = "C7H14OOH3-6O2 <=> C7H14OOH3-6 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.389e+23, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-6O2 <=> C7H14OOH3-6 + O2""",
)

entry(
    index = 2310,
    label = "C7H14OOH3-7O2 <=> C7H14OOH3-7 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.323e+20, 's^-1'), n=-1.65, Ea=(35280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-7O2 <=> C7H14OOH3-7 + O2""",
)

entry(
    index = 2311,
    label = "C7H14OOH4-1O2 <=> C7H14OOH4-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.668e+20, 's^-1'), n=-1.65, Ea=(35280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-1O2 <=> C7H14OOH4-1 + O2""",
)

entry(
    index = 2312,
    label = "C7H14OOH4-2O2 <=> C7H14OOH4-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.969e+22, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-2O2 <=> C7H14OOH4-2 + O2""",
)

entry(
    index = 2313,
    label = "C7H14OOH4-3O2 <=> C7H14OOH4-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.969e+22, 's^-1'), n=-2.38, Ea=(37600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-3O2 <=> C7H14OOH4-3 + O2""",
)

entry(
    index = 2314,
    label = "C7H14OOH1-2O2 <=> NC7KET12 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-2O2 <=> NC7KET12 + OH""",
)

entry(
    index = 2315,
    label = "C7H14OOH1-3O2 <=> NC7KET13 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-3O2 <=> NC7KET13 + OH""",
)

entry(
    index = 2316,
    label = "C7H14OOH1-4O2 <=> NC7KET14 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(18950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-4O2 <=> NC7KET14 + OH""",
)

entry(
    index = 2317,
    label = "C7H14OOH1-5O2 <=> NC7KET15 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.906e+08, 's^-1'), n=0, Ea=(22150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH1-5O2 <=> NC7KET15 + OH""",
)

entry(
    index = 2318,
    label = "C7H14OOH2-1O2 <=> NC7KET21 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-1O2 <=> NC7KET21 + OH""",
)

entry(
    index = 2319,
    label = "C7H14OOH2-3O2 <=> NC7KET23 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-3O2 <=> NC7KET23 + OH""",
)

entry(
    index = 2320,
    label = "C7H14OOH2-4O2 <=> NC7KET24 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-4O2 <=> NC7KET24 + OH""",
)

entry(
    index = 2321,
    label = "C7H14OOH2-5O2 <=> NC7KET25 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(15650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-5O2 <=> NC7KET25 + OH""",
)

entry(
    index = 2322,
    label = "C7H14OOH2-6O2 <=> NC7KET26 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.953e+08, 's^-1'), n=0, Ea=(18650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH2-6O2 <=> NC7KET26 + OH""",
)

entry(
    index = 2323,
    label = "C7H14OOH3-1O2 <=> NC7KET31 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-1O2 <=> NC7KET31 + OH""",
)

entry(
    index = 2324,
    label = "C7H14OOH3-2O2 <=> NC7KET32 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-2O2 <=> NC7KET32 + OH""",
)

entry(
    index = 2325,
    label = "C7H14OOH3-4O2 <=> NC7KET34 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-4O2 <=> NC7KET34 + OH""",
)

entry(
    index = 2326,
    label = "C7H14OOH3-5O2 <=> NC7KET35 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-5O2 <=> NC7KET35 + OH""",
)

entry(
    index = 2327,
    label = "C7H14OOH3-6O2 <=> NC7KET36 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(15650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-6O2 <=> NC7KET36 + OH""",
)

entry(
    index = 2328,
    label = "C7H14OOH3-7O2 <=> NC7KET37 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.953e+08, 's^-1'), n=0, Ea=(18650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH3-7O2 <=> NC7KET37 + OH""",
)

entry(
    index = 2329,
    label = "C7H14OOH4-1O2 <=> NC7KET41 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(15650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-1O2 <=> NC7KET41 + OH""",
)

entry(
    index = 2330,
    label = "C7H14OOH4-2O2 <=> NC7KET42 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-2O2 <=> NC7KET42 + OH""",
)

entry(
    index = 2331,
    label = "C7H14OOH4-3O2 <=> NC7KET43 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OOH4-3O2 <=> NC7KET43 + OH""",
)

entry(
    index = 2332,
    label = "NC7KET12 => NC5H11CHO + HCO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET12 => NC5H11CHO + HCO + OH""",
)

entry(
    index = 2333,
    label = "NC7KET13 => NC4H9CHO + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET13 => NC4H9CHO + CH2CHO + OH""",
)

entry(
    index = 2334,
    label = "NC7KET14 => NC3H7CHO + CH2CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET14 => NC3H7CHO + CH2CH2CHO + OH""",
)

entry(
    index = 2335,
    label = "NC7KET15 => C2H5CHO + C3H6CHO-1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET15 => C2H5CHO + C3H6CHO-1 + OH""",
)

entry(
    index = 2336,
    label = "NC7KET21 => CH2O + NC5H11CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET21 => CH2O + NC5H11CO + OH""",
)

entry(
    index = 2337,
    label = "NC7KET23 => NC4H9CHO + CH3CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET23 => NC4H9CHO + CH3CO + OH""",
)

entry(
    index = 2338,
    label = "NC7KET24 => NC3H7CHO + CH3COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET24 => NC3H7CHO + CH3COCH2 + OH""",
)

entry(
    index = 2339,
    label = "NC7KET25 => C2H5CHO + CH2CH2COCH3 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET25 => C2H5CHO + CH2CH2COCH3 + OH""",
)

entry(
    index = 2340,
    label = "NC7KET26 => CH3CHO + C3H6COCH3-1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET26 => CH3CHO + C3H6COCH3-1 + OH""",
)

entry(
    index = 2341,
    label = "NC7KET31 => CH2O + NC4H9COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET31 => CH2O + NC4H9COCH2 + OH""",
)

entry(
    index = 2342,
    label = "NC7KET32 => CH3CHO + NC4H9CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET32 => CH3CHO + NC4H9CO + OH""",
)

entry(
    index = 2343,
    label = "NC7KET34 => NC3H7CHO + C2H5CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET34 => NC3H7CHO + C2H5CO + OH""",
)

entry(
    index = 2344,
    label = "NC7KET35 => C2H5CHO + C2H5COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET35 => C2H5CHO + C2H5COCH2 + OH""",
)

entry(
    index = 2345,
    label = "NC7KET36 => CH3CHO + C2H5COC2H4P + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET36 => CH3CHO + C2H5COC2H4P + OH""",
)

entry(
    index = 2346,
    label = "NC7KET37 => CH2O + C3H6COC2H5-1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET37 => CH2O + C3H6COC2H5-1 + OH""",
)

entry(
    index = 2347,
    label = "NC7KET41 => CH2O + NC3H7COC2H4P + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET41 => CH2O + NC3H7COC2H4P + OH""",
)

entry(
    index = 2348,
    label = "NC7KET42 => CH3CHO + NC3H7COCH2 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET42 => CH3CHO + NC3H7COCH2 + OH""",
)

entry(
    index = 2349,
    label = "NC7KET43 => C2H5CHO + NC3H7CO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC7KET43 => C2H5CHO + NC3H7CO + OH""",
)

entry(
    index = 2350,
    label = "C7H14O1-2 + OH => PC4H9 + C2H3CHO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-2 + OH => PC4H9 + C2H3CHO + H2O""",
)

entry(
    index = 2351,
    label = "C7H14O1-3 + OH => C6H12-1 + HCO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-3 + OH => C6H12-1 + HCO + H2O""",
)

entry(
    index = 2352,
    label = "C7H14O1-4 + OH => C5H10-1 + CH2CHO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-4 + OH => C5H10-1 + CH2CHO + H2O""",
)

entry(
    index = 2353,
    label = "C7H14O1-5 + OH => C4H8-1 + CH2CH2CHO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-5 + OH => C4H8-1 + CH2CH2CHO + H2O""",
)

entry(
    index = 2354,
    label = "C7H14O2-3 + OH => C2H3COCH3 + NC3H7 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-3 + OH => C2H3COCH3 + NC3H7 + H2O""",
)

entry(
    index = 2355,
    label = "C7H14O2-4 + OH => CH3CO + C5H10-1 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-4 + OH => CH3CO + C5H10-1 + H2O""",
)

entry(
    index = 2356,
    label = "C7H14O2-5 + OH => CH3COCH2 + C4H8-1 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-5 + OH => CH3COCH2 + C4H8-1 + H2O""",
)

entry(
    index = 2357,
    label = "C7H14O2-6 + OH => CH2CH2COCH3 + C3H6 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-6 + OH => CH2CH2COCH3 + C3H6 + H2O""",
)

entry(
    index = 2358,
    label = "C7H14O3-4 + OH => C2H5COC2H3 + C2H5 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-4 + OH => C2H5COC2H3 + C2H5 + H2O""",
)

entry(
    index = 2359,
    label = "C7H14O3-5 + OH => C2H5CO + C4H8-1 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-5 + OH => C2H5CO + C4H8-1 + H2O""",
)

entry(
    index = 2360,
    label = "C7H14O1-2 + OH => CH2CO + C5H11-1 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-2 + OH => CH2CO + C5H11-1 + H2O""",
)

entry(
    index = 2361,
    label = "C7H14O1-3 + OH => C2H4 + NC4H9CO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-3 + OH => C2H4 + NC4H9CO + H2O""",
)

entry(
    index = 2362,
    label = "C7H14O1-4 + OH => C2H4 + NC3H7COCH2 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-4 + OH => C2H4 + NC3H7COCH2 + H2O""",
)

entry(
    index = 2363,
    label = "C7H14O1-5 + OH => C2H4 + C2H5COC2H4P + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-5 + OH => C2H4 + C2H5COC2H4P + H2O""",
)

entry(
    index = 2364,
    label = "C7H14O2-3 + OH => CH3CHCO + PC4H9 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-3 + OH => CH3CHCO + PC4H9 + H2O""",
)

entry(
    index = 2365,
    label = "C7H14O2-4 + OH => C3H6 + NC3H7CO + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-4 + OH => C3H6 + NC3H7CO + H2O""",
)

entry(
    index = 2366,
    label = "C7H14O2-5 + OH => C3H6 + C2H5COCH2 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-5 + OH => C3H6 + C2H5COCH2 + H2O""",
)

entry(
    index = 2367,
    label = "C7H14O2-6 + OH => CH3CHO + C5H91-4 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-6 + OH => CH3CHO + C5H91-4 + H2O""",
)

entry(
    index = 2368,
    label = "C7H14O3-4 + OH => C2H5CHCO + NC3H7 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-4 + OH => C2H5CHCO + NC3H7 + H2O""",
)

entry(
    index = 2369,
    label = "C7H14O3-5 + OH => C2H5CHO + C4H71-2 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-5 + OH => C2H5CHO + C4H71-2 + H2O""",
)

entry(
    index = 2370,
    label = "C7H14O1-2 + HO2 => PC4H9 + C2H3CHO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-2 + HO2 => PC4H9 + C2H3CHO + H2O2""",
)

entry(
    index = 2371,
    label = "C7H14O1-3 + HO2 => C6H12-1 + HCO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-3 + HO2 => C6H12-1 + HCO + H2O2""",
)

entry(
    index = 2372,
    label = "C7H14O1-4 + HO2 => C5H10-1 + CH2CHO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-4 + HO2 => C5H10-1 + CH2CHO + H2O2""",
)

entry(
    index = 2373,
    label = "C7H14O1-5 + HO2 => C4H8-1 + CH2CH2CHO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-5 + HO2 => C4H8-1 + CH2CH2CHO + H2O2""",
)

entry(
    index = 2374,
    label = "C7H14O2-3 + HO2 => C2H3COCH3 + NC3H7 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-3 + HO2 => C2H3COCH3 + NC3H7 + H2O2""",
)

entry(
    index = 2375,
    label = "C7H14O2-4 + HO2 => CH3CO + C5H10-1 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-4 + HO2 => CH3CO + C5H10-1 + H2O2""",
)

entry(
    index = 2376,
    label = "C7H14O2-5 + HO2 => CH3COCH2 + C4H8-1 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-5 + HO2 => CH3COCH2 + C4H8-1 + H2O2""",
)

entry(
    index = 2377,
    label = "C7H14O2-6 + HO2 => CH2CH2COCH3 + C3H6 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-6 + HO2 => CH2CH2COCH3 + C3H6 + H2O2""",
)

entry(
    index = 2378,
    label = "C7H14O3-4 + HO2 => C2H5COC2H3 + C2H5 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-4 + HO2 => C2H5COC2H3 + C2H5 + H2O2""",
)

entry(
    index = 2379,
    label = "C7H14O3-5 + HO2 => C2H5CO + C4H8-1 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-5 + HO2 => C2H5CO + C4H8-1 + H2O2""",
)

entry(
    index = 2380,
    label = "C7H14O1-2 + HO2 => CH2CO + C5H11-1 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-2 + HO2 => CH2CO + C5H11-1 + H2O2""",
)

entry(
    index = 2381,
    label = "C7H14O1-3 + HO2 => C2H4 + NC4H9CO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-3 + HO2 => C2H4 + NC4H9CO + H2O2""",
)

entry(
    index = 2382,
    label = "C7H14O1-4 + HO2 => C2H4 + NC3H7COCH2 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-4 + HO2 => C2H4 + NC3H7COCH2 + H2O2""",
)

entry(
    index = 2383,
    label = "C7H14O1-5 + HO2 => C2H4 + C2H5COC2H4P + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O1-5 + HO2 => C2H4 + C2H5COC2H4P + H2O2""",
)

entry(
    index = 2384,
    label = "C7H14O2-3 + HO2 => CH3CHCO + PC4H9 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-3 + HO2 => CH3CHCO + PC4H9 + H2O2""",
)

entry(
    index = 2385,
    label = "C7H14O2-4 + HO2 => C3H6 + NC3H7CO + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-4 + HO2 => C3H6 + NC3H7CO + H2O2""",
)

entry(
    index = 2386,
    label = "C7H14O2-5 + HO2 => C3H6 + C2H5COCH2 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-5 + HO2 => C3H6 + C2H5COCH2 + H2O2""",
)

entry(
    index = 2387,
    label = "C7H14O2-6 + HO2 => CH3CHO + C5H91-4 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O2-6 + HO2 => CH3CHO + C5H91-4 + H2O2""",
)

entry(
    index = 2388,
    label = "C7H14O3-4 + HO2 => C2H5CHCO + NC3H7 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-4 + HO2 => C2H5CHCO + NC3H7 + H2O2""",
)

entry(
    index = 2389,
    label = "C7H14O3-5 + HO2 => C2H5CHO + C4H71-2 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14O3-5 + HO2 => C2H5CHO + C4H71-2 + H2O2""",
)

entry(
    index = 2390,
    label = "C7H14OH-1 <=> C7H14-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.43e+14, 's^-1'), n=-0.53, Ea=(27830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OH-1 <=> C7H14-1 + OH""",
)

entry(
    index = 2391,
    label = "O2C7H14OH-1 <=> C7H14OH-1 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.983e+21, 's^-1'), n=-1.98, Ea=(37820, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C7H14OH-1 <=> C7H14OH-1 + O2""",
)

entry(
    index = 2392,
    label = "O2C7H14OH-1 => NC5H11CHO + CH2O + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C7H14OH-1 => NC5H11CHO + CH2O + OH""",
)

entry(
    index = 2393,
    label = "C7H14OH-2 <=> C7H14-2 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.044e+16, 's^-1'), n=-1, Ea=(29400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OH-2 <=> C7H14-2 + OH""",
)

entry(
    index = 2394,
    label = "O2C7H14OH-2 <=> C7H14OH-2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.775e+21, 's^-1'), n=-2.06, Ea=(37860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C7H14OH-2 <=> C7H14OH-2 + O2""",
)

entry(
    index = 2395,
    label = "O2C7H14OH-1 => NC4H9CHO + CH3CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C7H14OH-1 => NC4H9CHO + CH3CHO + OH""",
)

entry(
    index = 2396,
    label = "C7H14OH-3 <=> C7H14-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.773e+15, 's^-1'), n=-0.93, Ea=(29490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C7H14OH-3 <=> C7H14-3 + OH""",
)

entry(
    index = 2397,
    label = "O2C7H14OH-3 <=> C7H14OH-3 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.775e+21, 's^-1'), n=-2.06, Ea=(37860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C7H14OH-3 <=> C7H14OH-3 + O2""",
)

entry(
    index = 2398,
    label = "O2C7H14OH-3 => NC3H7CHO + C2H5CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2C7H14OH-3 => NC3H7CHO + C2H5CHO + OH""",
)

entry(
    index = 2399,
    label = "NC5H11CHO + O2 <=> NC5H11CO + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0.5, Ea=(42200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + O2 <=> NC5H11CO + HO2""",
)

entry(
    index = 2400,
    label = "NC5H11CHO + OH <=> NC5H11CO + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.69e+10, 'cm^3/(mol*s)'),
        n = 0.76,
        Ea = (-340, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + OH <=> NC5H11CO + H2O""",
)

entry(
    index = 2401,
    label = "NC5H11CHO + H <=> NC5H11CO + H2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+13, 'cm^3/(mol*s)'), n=0, Ea=(4200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + H <=> NC5H11CO + H2""",
)

entry(
    index = 2402,
    label = "NC5H11CHO + O <=> NC5H11CO + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5e+12, 'cm^3/(mol*s)'), n=0, Ea=(1790, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + O <=> NC5H11CO + OH""",
)

entry(
    index = 2403,
    label = "NC5H11CHO + HO2 <=> NC5H11CO + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + HO2 <=> NC5H11CO + H2O2""",
)

entry(
    index = 2404,
    label = "NC5H11CHO + CH3 <=> NC5H11CO + CH4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.7e+12, 'cm^3/(mol*s)'), n=0, Ea=(8440, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3 <=> NC5H11CO + CH4""",
)

entry(
    index = 2405,
    label = "NC5H11CHO + CH3O <=> NC5H11CO + CH3OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.15e+11, 'cm^3/(mol*s)'), n=0, Ea=(1280, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O <=> NC5H11CO + CH3OH""",
)

entry(
    index = 2406,
    label = "NC5H11CHO + CH3O2 <=> NC5H11CO + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(9500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O2 <=> NC5H11CO + CH3O2H""",
)

entry(
    index = 2407,
    label = "NC5H11CHO + OH <=> C5H10CHO-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + OH <=> C5H10CHO-1 + H2O""",
)

entry(
    index = 2408,
    label = "NC5H11CHO + OH <=> C5H10CHO-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + OH <=> C5H10CHO-2 + H2O""",
)

entry(
    index = 2409,
    label = "NC5H11CHO + OH <=> C5H10CHO-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + OH <=> C5H10CHO-3 + H2O""",
)

entry(
    index = 2410,
    label = "NC5H11CHO + OH <=> C5H10CHO-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + OH <=> C5H10CHO-4 + H2O""",
)

entry(
    index = 2411,
    label = "NC5H11CHO + OH <=> C5H10CHO-5 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + OH <=> C5H10CHO-5 + H2O""",
)

entry(
    index = 2412,
    label = "NC5H11CO <=> C5H11-1 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(9600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CO <=> C5H11-1 + CO""",
)

entry(
    index = 2413,
    label = "NC5H11CHO + HO2 <=> C5H10CHO-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27600, 'cm^3/(mol*s)'), n=2.55, Ea=(16480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + HO2 <=> C5H10CHO-1 + H2O2""",
)

entry(
    index = 2414,
    label = "NC5H11CHO + HO2 <=> C5H10CHO-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14750, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + HO2 <=> C5H10CHO-2 + H2O2""",
)

entry(
    index = 2415,
    label = "NC5H11CHO + HO2 <=> C5H10CHO-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14750, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + HO2 <=> C5H10CHO-3 + H2O2""",
)

entry(
    index = 2416,
    label = "NC5H11CHO + HO2 <=> C5H10CHO-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14750, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + HO2 <=> C5H10CHO-4 + H2O2""",
)

entry(
    index = 2417,
    label = "NC5H11CHO + HO2 <=> C5H10CHO-5 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(29500, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + HO2 <=> C5H10CHO-5 + H2O2""",
)

entry(
    index = 2418,
    label = "NC5H11CHO + CH3O2 <=> C5H10CHO-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O2 <=> C5H10CHO-1 + CH3O2H""",
)

entry(
    index = 2419,
    label = "NC5H11CHO + CH3O2 <=> C5H10CHO-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O2 <=> C5H10CHO-2 + CH3O2H""",
)

entry(
    index = 2420,
    label = "NC5H11CHO + CH3O2 <=> C5H10CHO-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O2 <=> C5H10CHO-3 + CH3O2H""",
)

entry(
    index = 2421,
    label = "NC5H11CHO + CH3O2 <=> C5H10CHO-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O2 <=> C5H10CHO-4 + CH3O2H""",
)

entry(
    index = 2422,
    label = "NC5H11CHO + CH3O2 <=> C5H10CHO-5 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.98e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC5H11CHO + CH3O2 <=> C5H10CHO-5 + CH3O2H""",
)

entry(
    index = 2423,
    label = "C5H10CHO-1 <=> C2H4 + C3H6CHO-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.68e+18, 's^-1'), n=-1.58, Ea=(30410, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-1 <=> C2H4 + C3H6CHO-1""",
)

entry(
    index = 2424,
    label = "C5H10CHO-2 <=> C3H6 + CH2CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.38e+17, 's^-1'), n=-1.31, Ea=(31970, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-2 <=> C3H6 + CH2CH2CHO""",
)

entry(
    index = 2425,
    label = "C5H10CHO-3 <=> C4H8-1 + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.27e+16, 's^-1'), n=-1.43, Ea=(25990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-3 <=> C4H8-1 + CH2CHO""",
)

entry(
    index = 2426,
    label = "C5H10CHO-3 <=> C4H7CHO1-4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.373e+14, 's^-1'), n=-0.56, Ea=(31320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-3 <=> C4H7CHO1-4 + CH3""",
)

entry(
    index = 2427,
    label = "C5H10CHO-4 <=> AC3H5CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.188e+17, 's^-1'), n=-1.37, Ea=(33230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-4 <=> AC3H5CHO + C2H5""",
)

entry(
    index = 2428,
    label = "C5H10CHO-4 <=> C5H10-1 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.057e+14, 's^-1'), n=-0.41, Ea=(26330, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-4 <=> C5H10-1 + HCO""",
)

entry(
    index = 2429,
    label = "C5H10CHO-5 <=> C2H3CHO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.564e+19, 's^-1'), n=-1.53, Ea=(33310, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H10CHO-5 <=> C2H3CHO + NC3H7""",
)

entry(
    index = 2430,
    label = "C4H7CHO1-4 + OH <=> C4H7CO1-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.37e+12, 'cm^3/(mol*s)'), n=0, Ea=(-616, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + OH <=> C4H7CO1-4 + H2O""",
)

entry(
    index = 2431,
    label = "C4H7CHO1-4 + OH <=> C4H6CHO1-43 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.08e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + OH <=> C4H6CHO1-43 + H2O""",
)

entry(
    index = 2432,
    label = "C4H7CHO1-4 + OH <=> C4H6CHO1-44 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.67e+07, 'cm^3/(mol*s)'),
        n = 1.61,
        Ea = (-35, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + OH <=> C4H6CHO1-44 + H2O""",
)

entry(
    index = 2433,
    label = "C4H7CHO1-4 + OH <=> CH3CHO + CH2CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + OH <=> CH3CHO + CH2CH2CHO""",
)

entry(
    index = 2434,
    label = "C4H7CHO1-4 + HO2 <=> C4H7CO1-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + HO2 <=> C4H7CO1-4 + H2O2""",
)

entry(
    index = 2435,
    label = "C4H7CHO1-4 + HO2 <=> C4H6CHO1-43 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + HO2 <=> C4H6CHO1-43 + H2O2""",
)

entry(
    index = 2436,
    label = "C4H7CHO1-4 + HO2 <=> C4H6CHO1-44 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(14750, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + HO2 <=> C4H6CHO1-44 + H2O2""",
)

entry(
    index = 2437,
    label = "C4H7CHO1-4 + CH3O2 <=> C4H7CO1-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + CH3O2 <=> C4H7CO1-4 + CH3O2H""",
)

entry(
    index = 2438,
    label = "C4H7CHO1-4 + CH3O2 <=> C4H6CHO1-43 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + CH3O2 <=> C4H6CHO1-43 + CH3O2H""",
)

entry(
    index = 2439,
    label = "C4H7CHO1-4 + CH3O2 <=> C4H6CHO1-44 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-4 + CH3O2 <=> C4H6CHO1-44 + CH3O2H""",
)

entry(
    index = 2440,
    label = "C4H7CO1-4 <=> C4H71-4 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.372e+18, 's^-1'), n=-1.76, Ea=(15230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CO1-4 <=> C4H71-4 + CO""",
)

entry(
    index = 2441,
    label = "C4H6CHO1-43 <=> C4H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.834e+15, 's^-1'), n=-0.79, Ea=(33540, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6CHO1-43 <=> C4H6 + HCO""",
)

entry(
    index = 2442,
    label = "C4H6CHO1-44 <=> C2H3CHO + C2H3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.118e+14, 's^-1'), n=-0.39, Ea=(37160, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6CHO1-44 <=> C2H3CHO + C2H3""",
)

entry(
    index = 2443,
    label = "NC4H9COCH3 + OH <=> C4H8COCH3-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.065e+07, 'cm^3/(mol*s)'),
        n = 1.73,
        Ea = (753, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + OH <=> C4H8COCH3-1 + H2O""",
)

entry(
    index = 2444,
    label = "NC4H9COCH3 + OH <=> C4H8COCH3-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.615e+07, 'cm^3/(mol*s)'),
        n = 1.64,
        Ea = (-247, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + OH <=> C4H8COCH3-2 + H2O""",
)

entry(
    index = 2445,
    label = "NC4H9COCH3 + OH <=> C4H8COCH3-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.615e+07, 'cm^3/(mol*s)'),
        n = 1.64,
        Ea = (-247, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + OH <=> C4H8COCH3-3 + H2O""",
)

entry(
    index = 2446,
    label = "NC4H9COCH3 + OH <=> C4H8COCH3-4 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + OH <=> C4H8COCH3-4 + H2O""",
)

entry(
    index = 2447,
    label = "NC4H9COCH3 + OH <=> NC4H9COCH2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + OH <=> NC4H9COCH2 + H2O""",
)

entry(
    index = 2448,
    label = "NC4H9COCH3 + HO2 <=> C4H8COCH3-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + HO2 <=> C4H8COCH3-1 + H2O2""",
)

entry(
    index = 2449,
    label = "NC4H9COCH3 + HO2 <=> C4H8COCH3-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + HO2 <=> C4H8COCH3-2 + H2O2""",
)

entry(
    index = 2450,
    label = "NC4H9COCH3 + HO2 <=> C4H8COCH3-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + HO2 <=> C4H8COCH3-3 + H2O2""",
)

entry(
    index = 2451,
    label = "NC4H9COCH3 + HO2 <=> C4H8COCH3-4 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + HO2 <=> C4H8COCH3-4 + H2O2""",
)

entry(
    index = 2452,
    label = "NC4H9COCH3 + HO2 <=> NC4H9COCH2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(14690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + HO2 <=> NC4H9COCH2 + H2O2""",
)

entry(
    index = 2453,
    label = "NC4H9COCH3 + CH3O2 <=> C4H8COCH3-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + CH3O2 <=> C4H8COCH3-1 + CH3O2H""",
)

entry(
    index = 2454,
    label = "NC4H9COCH3 + CH3O2 <=> C4H8COCH3-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + CH3O2 <=> C4H8COCH3-2 + CH3O2H""",
)

entry(
    index = 2455,
    label = "NC4H9COCH3 + CH3O2 <=> C4H8COCH3-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + CH3O2 <=> C4H8COCH3-3 + CH3O2H""",
)

entry(
    index = 2456,
    label = "NC4H9COCH3 + CH3O2 <=> C4H8COCH3-4 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + CH3O2 <=> C4H8COCH3-4 + CH3O2H""",
)

entry(
    index = 2457,
    label = "NC4H9COCH3 + CH3O2 <=> NC4H9COCH2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(17580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH3 + CH3O2 <=> NC4H9COCH2 + CH3O2H""",
)

entry(
    index = 2458,
    label = "C4H8COCH3-1 <=> CH2CH2COCH3 + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.133e+18, 's^-1'), n=-1.59, Ea=(30910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8COCH3-1 <=> CH2CH2COCH3 + C2H4""",
)

entry(
    index = 2459,
    label = "C4H8COCH3-2 <=> C3H6 + CH3COCH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.175e+15, 's^-1'), n=-0.79, Ea=(26220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8COCH3-2 <=> C3H6 + CH3COCH2""",
)

entry(
    index = 2460,
    label = "C4H8COCH3-3 <=> C4H8-1 + CH3CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.714e+18, 's^-1'), n=-1.61, Ea=(28250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8COCH3-3 <=> C4H8-1 + CH3CO""",
)

entry(
    index = 2461,
    label = "C4H8COCH3-4 <=> C2H3COCH3 + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.691e+19, 's^-1'), n=-1.61, Ea=(33750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H8COCH3-4 <=> C2H3COCH3 + C2H5""",
)

entry(
    index = 2462,
    label = "NC4H9COCH2 <=> PC4H9 + CH2CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.554e+18, 's^-1'), n=-1.41, Ea=(43140, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC4H9COCH2 <=> PC4H9 + CH2CO""",
)

entry(
    index = 2463,
    label = "C4H7OOH1-4 <=> C4H7O1-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.021e+20, 's^-1'), n=-1.53, Ea=(47040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7OOH1-4 <=> C4H7O1-4 + OH""",
)

entry(
    index = 2464,
    label = "C5H9OOH1-4 <=> C5H9O1-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.178e+20, 's^-1'), n=-1.38, Ea=(46050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9OOH1-4 <=> C5H9O1-4 + OH""",
)

entry(
    index = 2465,
    label = "C5H9OOH1-5 <=> C5H9O1-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.594e+20, 's^-1'), n=-1.5, Ea=(46990, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9OOH1-5 <=> C5H9O1-5 + OH""",
)

entry(
    index = 2466,
    label = "C6H11OOH1-4 <=> C6H11O1-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.234e+20, 's^-1'), n=-1.39, Ea=(46050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11OOH1-4 <=> C6H11O1-4 + OH""",
)

entry(
    index = 2467,
    label = "C6H11OOH1-5 <=> C6H11O1-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.234e+20, 's^-1'), n=-1.39, Ea=(46050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11OOH1-5 <=> C6H11O1-5 + OH""",
)

entry(
    index = 2468,
    label = "C4H7O1-4 <=> CH2O + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.412e+16, 's^-1'), n=-1.14, Ea=(7550, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7O1-4 <=> CH2O + C3H5-A""",
)

entry(
    index = 2469,
    label = "C5H9O1-4 <=> CH3CHO + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.715e+20, 's^-1'), n=-2.43, Ea=(5890, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O1-4 <=> CH3CHO + C3H5-A""",
)

entry(
    index = 2470,
    label = "C5H9O1-4 <=> AC3H5CHO + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.112e+17, 's^-1'), n=-1.21, Ea=(17960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O1-4 <=> AC3H5CHO + CH3""",
)

entry(
    index = 2471,
    label = "C5H9O1-5 <=> CH2O + C4H71-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.31e+17, 's^-1'), n=-1.33, Ea=(17940, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H9O1-5 <=> CH2O + C4H71-4""",
)

entry(
    index = 2472,
    label = "C6H11O1-4 <=> AC3H5CHO + C2H5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.268e+20, 's^-1'), n=-2.1, Ea=(18870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O1-4 <=> AC3H5CHO + C2H5""",
)

entry(
    index = 2473,
    label = "C6H11O1-4 <=> C2H5CHO + C3H5-A",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.211e+21, 's^-1'), n=-2.46, Ea=(5641, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O1-4 <=> C2H5CHO + C3H5-A""",
)

entry(
    index = 2474,
    label = "C6H11O1-5 <=> C4H7CHO1-4 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.399e+17, 's^-1'), n=-1.3, Ea=(16960, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O1-5 <=> C4H7CHO1-4 + CH3""",
)

entry(
    index = 2475,
    label = "C6H11O1-5 <=> CH3CHO + C4H71-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.999e+22, 's^-1'), n=-2.58, Ea=(18530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H11O1-5 <=> CH3CHO + C4H71-4""",
)

entry(
    index = 2476,
    label = "C5H91-1 <=> C2H2 + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.755e+15, 's^-1'), n=-0.67, Ea=(30800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-1 <=> C2H2 + NC3H7""",
)

entry(
    index = 2477,
    label = "C4H7CHO1-1 + OH <=> C4H7CO1-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.37e+12, 'cm^3/(mol*s)'), n=0, Ea=(-616, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + OH <=> C4H7CO1-1 + H2O""",
)

entry(
    index = 2478,
    label = "C4H7CHO1-1 + OH <=> C4H6CHO1-14 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (5.27e+09, 'cm^3/(mol*s)'),
        n = 0.97,
        Ea = (1586, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + OH <=> C4H6CHO1-14 + H2O""",
)

entry(
    index = 2479,
    label = "C4H7CHO1-1 + OH <=> C4H6CHO1-13 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.12e+06, 'cm^3/(mol*s)'), n=2, Ea=(-298, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + OH <=> C4H6CHO1-13 + H2O""",
)

entry(
    index = 2480,
    label = "C4H7CHO1-1 + OH <=> C2H5CHO + CH2CHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + OH <=> C2H5CHO + CH2CHO""",
)

entry(
    index = 2481,
    label = "C4H7CHO1-1 + HO2 <=> C4H7CO1-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(11920, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + HO2 <=> C4H7CO1-1 + H2O2""",
)

entry(
    index = 2482,
    label = "C4H7CHO1-1 + HO2 <=> C4H6CHO1-14 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(27600, 'cm^3/(mol*s)'), n=2.55, Ea=(16480, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + HO2 <=> C4H6CHO1-14 + H2O2""",
)

entry(
    index = 2483,
    label = "C4H7CHO1-1 + HO2 <=> C4H6CHO1-13 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9640, 'cm^3/(mol*s)'), n=2.6, Ea=(13910, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + HO2 <=> C4H6CHO1-13 + H2O2""",
)

entry(
    index = 2484,
    label = "C4H7CHO1-1 + CH3O2 <=> C4H7CO1-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.8e+12, 'cm^3/(mol*s)'), n=0, Ea=(13600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + CH3O2 <=> C4H7CO1-1 + CH3O2H""",
)

entry(
    index = 2485,
    label = "C4H7CHO1-1 + CH3O2 <=> C4H6CHO1-14 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.03e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + CH3O2 <=> C4H6CHO1-14 + CH3O2H""",
)

entry(
    index = 2486,
    label = "C4H7CHO1-1 + CH3O2 <=> C4H6CHO1-13 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CHO1-1 + CH3O2 <=> C4H6CHO1-13 + CH3O2H""",
)

entry(
    index = 2487,
    label = "C4H7CO1-1 <=> C4H71-1 + CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.51e+20, 's^-1'), n=-2.12, Ea=(40320, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H7CO1-1 <=> C4H71-1 + CO""",
)

entry(
    index = 2488,
    label = "C4H6CHO1-14 <=> C2H4 + CHCHCHO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.05e+16, 's^-1'), n=-1.33, Ea=(46870, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6CHO1-14 <=> C2H4 + CHCHCHO""",
)

entry(
    index = 2489,
    label = "C4H6CHO1-13 <=> C4H6 + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.954e+17, 's^-1'), n=-1.28, Ea=(46230, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H6CHO1-13 <=> C4H6 + HCO""",
)

entry(
    index = 2490,
    label = "NC3H7COC2H5 + OH <=> C3H6COC2H5-1 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (2.065e+07, 'cm^3/(mol*s)'),
        n = 1.73,
        Ea = (753, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + OH <=> C3H6COC2H5-1 + H2O""",
)

entry(
    index = 2491,
    label = "NC3H7COC2H5 + OH <=> C3H6COC2H5-2 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (3.615e+07, 'cm^3/(mol*s)'),
        n = 1.64,
        Ea = (-247, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + OH <=> C3H6COC2H5-2 + H2O""",
)

entry(
    index = 2492,
    label = "NC3H7COC2H5 + OH <=> C3H6COC2H5-3 + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + OH <=> C3H6COC2H5-3 + H2O""",
)

entry(
    index = 2493,
    label = "NC3H7COC2H5 + OH <=> NC3H7COC2H4P + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.1e+11, 'cm^3/(mol*s)'), n=0, Ea=(1192, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + OH <=> NC3H7COC2H4P + H2O""",
)

entry(
    index = 2494,
    label = "NC3H7COC2H5 + OH <=> NC3H7COC2H4S + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.45e+11, 'cm^3/(mol*s)'), n=0, Ea=(-228, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + OH <=> NC3H7COC2H4S + H2O""",
)

entry(
    index = 2495,
    label = "NC3H7COC2H5 + HO2 <=> C3H6COC2H5-1 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(16490, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + HO2 <=> C3H6COC2H5-1 + H2O2""",
)

entry(
    index = 2496,
    label = "NC3H7COC2H5 + HO2 <=> C3H6COC2H5-2 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.6e+12, 'cm^3/(mol*s)'), n=0, Ea=(17700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + HO2 <=> C3H6COC2H5-2 + H2O2""",
)

entry(
    index = 2497,
    label = "NC3H7COC2H5 + HO2 <=> C3H6COC2H5-3 + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + HO2 <=> C3H6COC2H5-3 + H2O2""",
)

entry(
    index = 2498,
    label = "NC3H7COC2H5 + HO2 <=> NC3H7COC2H4P + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(23800, 'cm^3/(mol*s)'), n=2.55, Ea=(14690, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + HO2 <=> NC3H7COC2H4P + H2O2""",
)

entry(
    index = 2499,
    label = "NC3H7COC2H5 + HO2 <=> NC3H7COC2H4S + H2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(8698, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + HO2 <=> NC3H7COC2H4S + H2O2""",
)

entry(
    index = 2500,
    label = "NC3H7COC2H5 + CH3O2 <=> C3H6COC2H5-1 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(19380, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + CH3O2 <=> C3H6COC2H5-1 + CH3O2H""",
)

entry(
    index = 2501,
    label = "NC3H7COC2H5 + CH3O2 <=> C3H6COC2H5-2 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.99e+12, 'cm^3/(mol*s)'), n=0, Ea=(17050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + CH3O2 <=> C3H6COC2H5-2 + CH3O2H""",
)

entry(
    index = 2502,
    label = "NC3H7COC2H5 + CH3O2 <=> C3H6COC2H5-3 + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + CH3O2 <=> C3H6COC2H5-3 + CH3O2H""",
)

entry(
    index = 2503,
    label = "NC3H7COC2H5 + CH3O2 <=> NC3H7COC2H4P + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.01e+12, 'cm^3/(mol*s)'), n=0, Ea=(17580, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + CH3O2 <=> NC3H7COC2H4P + CH3O2H""",
)

entry(
    index = 2504,
    label = "NC3H7COC2H5 + CH3O2 <=> NC3H7COC2H4S + CH3O2H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H5 + CH3O2 <=> NC3H7COC2H4S + CH3O2H""",
)

entry(
    index = 2505,
    label = "C3H6COC2H5-1 <=> C2H4 + C2H5COCH2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.218e+15, 's^-1'), n=-0.84, Ea=(23590, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6COC2H5-1 <=> C2H4 + C2H5COCH2""",
)

entry(
    index = 2506,
    label = "C3H6COC2H5-2 <=> C3H6 + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.051e+16, 's^-1'), n=-1.11, Ea=(26150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6COC2H5-2 <=> C3H6 + C2H5CO""",
)

entry(
    index = 2507,
    label = "C3H6COC2H5-3 <=> C2H5COC2H3 + CH3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.915e+15, 's^-1'), n=-0.68, Ea=(32300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H6COC2H5-3 <=> C2H5COC2H3 + CH3""",
)

entry(
    index = 2508,
    label = "NC3H7COC2H4P <=> NC3H7CO + C2H4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.398e+17, 's^-1'), n=-1.45, Ea=(26040, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H4P <=> NC3H7CO + C2H4""",
)

entry(
    index = 2509,
    label = "NC3H7COC2H4S <=> CH3CHCO + NC3H7",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.973e+19, 's^-1'), n=-1.49, Ea=(42860, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7COC2H4S <=> CH3CHCO + NC3H7""",
)

entry(
    index = 2510,
    label = "CHCHCHO + OH <=> CH2CHO + HCO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CHCHCHO + OH <=> CH2CHO + HCO""",
)

entry(
    index = 2511,
    label = "O2 + C6H12-1 => CH2O + NC4H9CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(37000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is O2 + C6H12-1 => CH2O + NC4H9CHO""",
)

entry(
    index = 2512,
    label = "C6H12-1 <=> C3H6 + C3H6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 's^-1'), n=0, Ea=(58000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 <=> C3H6 + C3H6""",
)

entry(
    index = 2513,
    label = "C6H111-3 + H <=> C6H12-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + H <=> C6H12-1""",
)

entry(
    index = 2514,
    label = "C5H91-5 + CH3 <=> C6H12-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H91-5 + CH3 <=> C6H12-1""",
)

entry(
    index = 2515,
    label = "PC4H9 + C2H3 <=> C6H12-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is PC4H9 + C2H3 <=> C6H12-1""",
)

entry(
    index = 2516,
    label = "C4H71-4 + C2H5 <=> C6H12-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-4 + C2H5 <=> C6H12-1""",
)

entry(
    index = 2517,
    label = "C6H111-3 + H <=> C6H12-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + H <=> C6H12-2""",
)

entry(
    index = 2518,
    label = "C6H112-4 + H <=> C6H12-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + H <=> C6H12-2""",
)

entry(
    index = 2519,
    label = "C5H92-5 + CH3 <=> C6H12-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C5H92-5 + CH3 <=> C6H12-2""",
)

entry(
    index = 2520,
    label = "NC3H7 + C3H5-T <=> C6H12-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC3H7 + C3H5-T <=> C6H12-2""",
)

entry(
    index = 2521,
    label = "C6H112-4 + H <=> C6H12-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+14, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + H <=> C6H12-3""",
)

entry(
    index = 2522,
    label = "C4H71-1 + C2H5 <=> C6H12-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+13, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-1 + C2H5 <=> C6H12-3""",
)

entry(
    index = 2523,
    label = "C6H12-1 + O2 <=> C6H111-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(37220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O2 <=> C6H111-3 + HO2""",
)

entry(
    index = 2524,
    label = "C6H12-1 + O2 <=> C6H111-4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(49640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O2 <=> C6H111-4 + HO2""",
)

entry(
    index = 2525,
    label = "C6H12-1 + O2 <=> C6H111-5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(49640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O2 <=> C6H111-5 + HO2""",
)

entry(
    index = 2526,
    label = "C6H12-1 + O2 <=> C6H111-6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O2 <=> C6H111-6 + HO2""",
)

entry(
    index = 2527,
    label = "C6H12-1 + O <=> C6H111-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(660000, 'cm^3/(mol*s)'), n=2.43, Ea=(1210, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O <=> C6H111-3 + OH""",
)

entry(
    index = 2528,
    label = "C6H12-1 + O <=> C6H111-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(551000, 'cm^3/(mol*s)'), n=2.45, Ea=(2830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O <=> C6H111-4 + OH""",
)

entry(
    index = 2529,
    label = "C6H12-1 + O <=> C6H111-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(551000, 'cm^3/(mol*s)'), n=2.45, Ea=(2830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O <=> C6H111-5 + OH""",
)

entry(
    index = 2530,
    label = "C6H12-1 + O <=> C6H111-6 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(980000, 'cm^3/(mol*s)'), n=2.43, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + O <=> C6H111-6 + OH""",
)

entry(
    index = 2531,
    label = "C6H12-2 + O2 <=> C6H111-3 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(5.5e+12, 'cm^3/(mol*s)'), n=0, Ea=(39900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O2 <=> C6H111-3 + HO2""",
)

entry(
    index = 2532,
    label = "C6H12-2 + O2 <=> C6H112-4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.2e+12, 'cm^3/(mol*s)'), n=0, Ea=(37220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O2 <=> C6H112-4 + HO2""",
)

entry(
    index = 2533,
    label = "C6H12-2 + O2 <=> C6H112-5 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+13, 'cm^3/(mol*s)'), n=0, Ea=(49640, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O2 <=> C6H112-5 + HO2""",
)

entry(
    index = 2534,
    label = "C6H12-2 + O2 <=> C6H112-6 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O2 <=> C6H112-6 + HO2""",
)

entry(
    index = 2535,
    label = "C6H12-2 + O <=> C6H111-3 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(441000, 'cm^3/(mol*s)'), n=2.42, Ea=(3150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O <=> C6H111-3 + OH""",
)

entry(
    index = 2536,
    label = "C6H12-2 + O <=> C6H112-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(660000, 'cm^3/(mol*s)'), n=2.43, Ea=(1210, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O <=> C6H112-4 + OH""",
)

entry(
    index = 2537,
    label = "C6H12-2 + O <=> C6H112-5 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(551000, 'cm^3/(mol*s)'), n=2.45, Ea=(2830, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O <=> C6H112-5 + OH""",
)

entry(
    index = 2538,
    label = "C6H12-2 + O <=> C6H112-6 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(980000, 'cm^3/(mol*s)'), n=2.43, Ea=(4750, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + O <=> C6H112-6 + OH""",
)

entry(
    index = 2539,
    label = "C6H12-3 + O2 <=> C6H112-4 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.4e+12, 'cm^3/(mol*s)'), n=0, Ea=(37220, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + O2 <=> C6H112-4 + HO2""",
)

entry(
    index = 2540,
    label = "C6H12-3 + O2 <=> C6H113-1 + HO2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+13, 'cm^3/(mol*s)'), n=0, Ea=(52290, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + O2 <=> C6H113-1 + HO2""",
)

entry(
    index = 2541,
    label = "C6H12-3 + O <=> C6H112-4 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.32e+06, 'cm^3/(mol*s)'),
        n = 2.43,
        Ea = (1210, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + O <=> C6H112-4 + OH""",
)

entry(
    index = 2542,
    label = "C6H12-3 + O <=> C6H113-1 + OH",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (1.96e+06, 'cm^3/(mol*s)'),
        n = 2.43,
        Ea = (4750, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + O <=> C6H113-1 + OH""",
)

entry(
    index = 2543,
    label = "C2H5 + C4H6 <=> C6H111-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(8.5e+10, 'cm^3/(mol*s)'), n=0, Ea=(8300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C4H6 <=> C6H111-3""",
)

entry(
    index = 2544,
    label = "CH3 + C5H81-4 <=> C6H111-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(7800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C5H81-4 <=> C6H111-4""",
)

entry(
    index = 2545,
    label = "C3H5-A + C3H6 <=> C6H111-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(16900, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-A + C3H6 <=> C6H111-5""",
)

entry(
    index = 2546,
    label = "C4H71-4 + C2H4 <=> C6H111-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(8200, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-4 + C2H4 <=> C6H111-6""",
)

entry(
    index = 2547,
    label = "CH3 + C5H81-3 <=> C6H112-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6e+10, 'cm^3/(mol*s)'), n=0, Ea=(7500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C5H81-3 <=> C6H112-4""",
)

entry(
    index = 2548,
    label = "C3H5-T + C3H6 <=> C6H112-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(3100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C3H5-T + C3H6 <=> C6H112-5""",
)

entry(
    index = 2549,
    label = "C4H71-3 + C2H4 <=> C6H112-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+11, 'cm^3/(mol*s)'), n=0, Ea=(13050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-3 + C2H4 <=> C6H112-6""",
)

entry(
    index = 2550,
    label = "C4H71-1 + C2H4 <=> C6H113-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.3e+11, 'cm^3/(mol*s)'), n=0, Ea=(3100, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C4H71-1 + C2H4 <=> C6H113-1""",
)

entry(
    index = 2551,
    label = "C6H111-6 <=> C6H111-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+10, 's^-1'), n=0.67, Ea=(36000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-6 <=> C6H111-4""",
)

entry(
    index = 2552,
    label = "C6H111-6 <=> C6H111-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.56e+10, 's^-1'), n=0.88, Ea=(37300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-6 <=> C6H111-5""",
)

entry(
    index = 2553,
    label = "C6H112-6 <=> C6H112-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+10, 's^-1'), n=0.67, Ea=(28400, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-6 <=> C6H112-4""",
)

entry(
    index = 2554,
    label = "C6H112-6 <=> C6H112-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.56e+10, 's^-1'), n=0.88, Ea=(37300, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-6 <=> C6H112-5""",
)

entry(
    index = 2555,
    label = "C6H113-1 <=> C6H112-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.56e+10, 's^-1'), n=0.88, Ea=(29600, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113-1 <=> C6H112-4""",
)

entry(
    index = 2556,
    label = "C6H111-3 + O2 <=> C6H112O2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + O2 <=> C6H112O2-1""",
)

entry(
    index = 2557,
    label = "C6H111-3 + O2 <=> C6H111O2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-3 + O2 <=> C6H111O2-3""",
)

entry(
    index = 2558,
    label = "C6H111-4 + O2 <=> C6H111O2-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-4 + O2 <=> C6H111O2-4""",
)

entry(
    index = 2559,
    label = "C6H111-5 + O2 <=> C6H111O2-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-5 + O2 <=> C6H111O2-5""",
)

entry(
    index = 2560,
    label = "C6H111-6 + O2 <=> C6H111O2-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-6 + O2 <=> C6H111O2-6""",
)

entry(
    index = 2561,
    label = "C6H112-4 + O2 <=> C6H113O2-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + O2 <=> C6H113O2-2""",
)

entry(
    index = 2562,
    label = "C6H112-4 + O2 <=> C6H112O2-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + O2 <=> C6H112O2-4""",
)

entry(
    index = 2563,
    label = "C6H112-5 + O2 <=> C6H112O2-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-5 + O2 <=> C6H112O2-5""",
)

entry(
    index = 2564,
    label = "C6H112-6 + O2 <=> C6H112O2-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-6 + O2 <=> C6H112O2-6""",
)

entry(
    index = 2565,
    label = "C6H113-1 + O2 <=> C6H113O2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113-1 + O2 <=> C6H113O2-1""",
)

entry(
    index = 2566,
    label = "C6H111O2-3 => C6H101-3 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.004e+39, 's^-1'), n=-8.11, Ea=(40500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-3 => C6H101-3 + HO2""",
)

entry(
    index = 2567,
    label = "C6H111O2-4 => C6H101-3 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.04e+38, 's^-1'), n=-8.11, Ea=(38500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-4 => C6H101-3 + HO2""",
)

entry(
    index = 2568,
    label = "C6H111O2-5 => C6H101-4 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.004e+39, 's^-1'), n=-8.11, Ea=(40500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-5 => C6H101-4 + HO2""",
)

entry(
    index = 2569,
    label = "C6H111O2-6 => C6H101-5 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.004e+39, 's^-1'), n=-8.11, Ea=(40500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-6 => C6H101-5 + HO2""",
)

entry(
    index = 2570,
    label = "C6H111O2-5 => C6H101-5 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.008e+43, 's^-1'), n=-9.41, Ea=(41500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-5 => C6H101-5 + HO2""",
)

entry(
    index = 2571,
    label = "C6H112O2-4 => C6H101-3 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.04e+38, 's^-1'), n=-8.11, Ea=(40500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-4 => C6H101-3 + HO2""",
)

entry(
    index = 2572,
    label = "C6H112O2-5 => C6H102-4 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.04e+38, 's^-1'), n=-8.11, Ea=(37500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-5 => C6H102-4 + HO2""",
)

entry(
    index = 2573,
    label = "C6H112O2-6 => C6H101-4 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.04e+38, 's^-1'), n=-8.11, Ea=(40500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-6 => C6H101-4 + HO2""",
)

entry(
    index = 2574,
    label = "C6H113O2-2 => C6H101-3 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5.08e+42, 's^-1'), n=-9.41, Ea=(41500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113O2-2 => C6H101-3 + HO2""",
)

entry(
    index = 2575,
    label = "C6H113O2-1 => C6H101-3 + HO2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.04e+38, 's^-1'), n=-8.11, Ea=(37500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113O2-1 => C6H101-3 + HO2""",
)

entry(
    index = 2576,
    label = "C6H111O2-3 <=> C6H101OOH3-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-3 <=> C6H101OOH3-4""",
)

entry(
    index = 2577,
    label = "C6H111O2-3 <=> C6H101OOH3-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-3 <=> C6H101OOH3-5""",
)

entry(
    index = 2578,
    label = "C6H111O2-3 <=> C6H101OOH3-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.376e+09, 's^-1'), n=0, Ea=(21950, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-3 <=> C6H101OOH3-6""",
)

entry(
    index = 2579,
    label = "C6H111O2-4 <=> C6H101OOH4-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(24450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-4 <=> C6H101OOH4-3""",
)

entry(
    index = 2580,
    label = "C6H111O2-4 <=> C6H101OOH4-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-4 <=> C6H101OOH4-5""",
)

entry(
    index = 2581,
    label = "C6H111O2-4 <=> C6H101OOH4-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+10, 's^-1'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-4 <=> C6H101OOH4-6""",
)

entry(
    index = 2582,
    label = "C6H111O2-5 <=> C6H101OOH5-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-5 <=> C6H101OOH5-3""",
)

entry(
    index = 2583,
    label = "C6H111O2-5 <=> C6H101OOH5-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-5 <=> C6H101OOH5-4""",
)

entry(
    index = 2584,
    label = "C6H111O2-5 <=> C6H101OOH5-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-5 <=> C6H101OOH5-6""",
)

entry(
    index = 2585,
    label = "C6H111O2-6 <=> C6H101OOH6-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.125e+09, 's^-1'), n=0, Ea=(16650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-6 <=> C6H101OOH6-3""",
)

entry(
    index = 2586,
    label = "C6H111O2-6 <=> C6H101OOH6-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(20450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-6 <=> C6H101OOH6-4""",
)

entry(
    index = 2587,
    label = "C6H111O2-6 <=> C6H101OOH6-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-6 <=> C6H101OOH6-5""",
)

entry(
    index = 2588,
    label = "C6H112O2-4 <=> C6H102OOH4-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-4 <=> C6H102OOH4-5""",
)

entry(
    index = 2589,
    label = "C6H112O2-4 <=> C6H102OOH4-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.75e+10, 's^-1'), n=0, Ea=(24000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-4 <=> C6H102OOH4-6""",
)

entry(
    index = 2590,
    label = "C6H112O2-5 <=> C6H102OOH5-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(24450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-5 <=> C6H102OOH5-4""",
)

entry(
    index = 2591,
    label = "C6H112O2-5 <=> C6H102OOH5-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-5 <=> C6H102OOH5-6""",
)

entry(
    index = 2592,
    label = "C6H112O2-6 <=> C6H102OOH6-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(18450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-6 <=> C6H102OOH6-4""",
)

entry(
    index = 2593,
    label = "C6H112O2-6 <=> C6H102OOH6-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-6 <=> C6H102OOH6-5""",
)

entry(
    index = 2594,
    label = "C6H113O2-1 <=> C6H103OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(24450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113O2-1 <=> C6H103OOH1-2""",
)

entry(
    index = 2595,
    label = "C6H113O2-2 <=> C6H103OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3e+11, 's^-1'), n=0, Ea=(29000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113O2-2 <=> C6H103OOH2-1""",
)

entry(
    index = 2596,
    label = "C6H101OOH3-4 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-4 => ETES1 + OH""",
)

entry(
    index = 2597,
    label = "C6H101OOH3-5 => MVOX + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-5 => MVOX + OH""",
)

entry(
    index = 2598,
    label = "C6H101OOH3-4 => VTHF + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-4 => VTHF + OH""",
)

entry(
    index = 2599,
    label = "C6H101OOH4-3 => EDHF + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-3 => EDHF + OH""",
)

entry(
    index = 2600,
    label = "C6H101OOH4-3 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-3 => ETES1 + OH""",
)

entry(
    index = 2601,
    label = "C6H101OOH4-5 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-5 => ETES1 + OH""",
)

entry(
    index = 2602,
    label = "C6H101OOH4-6 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-6 => ETES1 + OH""",
)

entry(
    index = 2603,
    label = "C6H101OOH5-3 => MVOX + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-3 => MVOX + OH""",
)

entry(
    index = 2604,
    label = "C6H101OOH5-4 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-4 => ETES1 + OH""",
)

entry(
    index = 2605,
    label = "C6H101OOH5-6 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-6 => ETES1 + OH""",
)

entry(
    index = 2606,
    label = "C6H101OOH6-3 => VTHF + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(9.375e+09, 's^-1'), n=0, Ea=(7000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-3 => VTHF + OH""",
)

entry(
    index = 2607,
    label = "C6H101OOH6-4 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-4 => ETES1 + OH""",
)

entry(
    index = 2608,
    label = "C6H101OOH6-5 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-5 => ETES1 + OH""",
)

entry(
    index = 2609,
    label = "C6H102OOH4-5 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-5 => ETES1 + OH""",
)

entry(
    index = 2610,
    label = "C6H102OOH4-6 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-6 => ETES1 + OH""",
)

entry(
    index = 2611,
    label = "C6H102OOH5-4 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-4 => ETES1 + OH""",
)

entry(
    index = 2612,
    label = "C6H102OOH5-6 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-6 => ETES1 + OH""",
)

entry(
    index = 2613,
    label = "C6H102OOH6-4 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(7.5e+10, 's^-1'), n=0, Ea=(15250, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-4 => ETES1 + OH""",
)

entry(
    index = 2614,
    label = "C6H102OOH6-5 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-5 => ETES1 + OH""",
)

entry(
    index = 2615,
    label = "C6H103OOH2-1 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH2-1 => ETES1 + OH""",
)

entry(
    index = 2616,
    label = "C6H103OOH1-2 => ETES1 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(6e+11, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH1-2 => ETES1 + OH""",
)

entry(
    index = 2617,
    label = "HO2 + C6H101-3 <=> C6H101OOH3-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.85e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-3 <=> C6H101OOH3-4""",
)

entry(
    index = 2618,
    label = "HO2 + C6H101-3 <=> C6H101OOH4-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.85e+11, 'cm^3/(mol*s)'), n=0, Ea=(9000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-3 <=> C6H101OOH4-3""",
)

entry(
    index = 2619,
    label = "HO2 + C6H101-4 <=> C6H101OOH4-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(11800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-4 <=> C6H101OOH4-5""",
)

entry(
    index = 2620,
    label = "HO2 + C6H101-4 <=> C6H101OOH5-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.4e+11, 'cm^3/(mol*s)'), n=0, Ea=(11800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-4 <=> C6H101OOH5-4""",
)

entry(
    index = 2621,
    label = "HO2 + C6H101-5 <=> C6H101OOH5-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(13700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-5 <=> C6H101OOH5-6""",
)

entry(
    index = 2622,
    label = "HO2 + C6H101-5 <=> C6H101OOH6-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(13700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-5 <=> C6H101OOH6-5""",
)

entry(
    index = 2623,
    label = "HO2 + C6H102-4 <=> C6H102OOH4-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(8800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H102-4 <=> C6H102OOH4-5""",
)

entry(
    index = 2624,
    label = "HO2 + C6H102-4 <=> C6H102OOH5-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(3.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(7800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H102-4 <=> C6H102OOH5-4""",
)

entry(
    index = 2625,
    label = "HO2 + C6H101-4 <=> C6H102OOH5-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(13700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-4 <=> C6H102OOH5-6""",
)

entry(
    index = 2626,
    label = "HO2 + C6H101-4 <=> C6H102OOH6-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(13700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-4 <=> C6H102OOH6-5""",
)

entry(
    index = 2627,
    label = "HO2 + C6H101-3 <=> C6H103OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.85e+11, 'cm^3/(mol*s)'), n=0, Ea=(10000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-3 <=> C6H103OOH2-1""",
)

entry(
    index = 2628,
    label = "HO2 + C6H101-3 <=> C6H103OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.85e+11, 'cm^3/(mol*s)'), n=0, Ea=(9000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C6H101-3 <=> C6H103OOH1-2""",
)

entry(
    index = 2629,
    label = "C6H12-1 + HO2 <=> C6H12OOH1-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(13700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + HO2 <=> C6H12OOH1-2""",
)

entry(
    index = 2630,
    label = "C6H12-1 + HO2 <=> C6H12OOH2-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(13700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-1 + HO2 <=> C6H12OOH2-1""",
)

entry(
    index = 2631,
    label = "C6H12-2 + HO2 <=> C6H12OOH2-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(11800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + HO2 <=> C6H12OOH2-3""",
)

entry(
    index = 2632,
    label = "C6H12-2 + HO2 <=> C6H12OOH3-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(11800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-2 + HO2 <=> C6H12OOH3-2""",
)

entry(
    index = 2633,
    label = "C6H12-3 + HO2 <=> C6H12OOH3-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(6.8e+11, 'cm^3/(mol*s)'), n=0, Ea=(11800, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + HO2 <=> C6H12OOH3-4""",
)

entry(
    index = 2634,
    label = "C6H101OOH3-4 => CH3 + AC3H5OOH + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(33500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-4 => CH3 + AC3H5OOH + C2H2""",
)

entry(
    index = 2635,
    label = "C6H101OOH3-5 => OH + C3H6 + C2H3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(23000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-5 => OH + C3H6 + C2H3CHO""",
)

entry(
    index = 2636,
    label = "C6H101OOH3-6 => HO2 + C4H6 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-6 => HO2 + C4H6 + C2H4""",
)

entry(
    index = 2637,
    label = "C6H101OOH4-3 => C2H5 + CH3O2H + C3H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(35000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-3 => C2H5 + CH3O2H + C3H2""",
)

entry(
    index = 2638,
    label = "C6H101OOH4-5 => C3H5-A + AC3H5OOH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(23000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-5 => C3H5-A + AC3H5OOH""",
)

entry(
    index = 2639,
    label = "C6H101OOH4-6 => OH + C2H3COCH3 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-6 => OH + C2H3COCH3 + C2H4""",
)

entry(
    index = 2640,
    label = "C6H101OOH5-3 => OH + C4H6 + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(35000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-3 => OH + C4H6 + CH3CHO""",
)

entry(
    index = 2641,
    label = "C6H101OOH5-4 => CH3 + CH3O2H + C2H2 + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(33500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-4 => CH3 + CH3O2H + C2H2 + C2H2""",
)

entry(
    index = 2642,
    label = "C6H101OOH5-6 => C4H71-4 + C2H3OOH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-6 => C4H71-4 + C2H3OOH""",
)

entry(
    index = 2643,
    label = "C6H101OOH6-3 => HO2 + C4H6 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(35000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-3 => HO2 + C4H6 + C2H4""",
)

entry(
    index = 2644,
    label = "C6H101OOH6-4 => OH + C5H81-4 + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(29500, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-4 => OH + C5H81-4 + CH2O""",
)

entry(
    index = 2645,
    label = "C6H101OOH6-5 => C3H5-A + AC3H5OOH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(23000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-5 => C3H5-A + AC3H5OOH""",
)

entry(
    index = 2646,
    label = "C6H102OOH4-6 => OH + C2H3COCH3 + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-6 => OH + C2H3COCH3 + C2H4""",
)

entry(
    index = 2647,
    label = "C6H102OOH5-6 => C4H71-3 + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(22000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-6 => C4H71-3 + CH2CHO + OH""",
)

entry(
    index = 2648,
    label = "C6H102OOH6-4 => OH + C5H81-3 + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+13, 's^-1'), n=0, Ea=(35000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-4 => OH + C5H81-3 + CH2O""",
)

entry(
    index = 2649,
    label = "C6H102OOH6-5 => C3H5-T + AC3H5OOH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-5 => C3H5-T + AC3H5OOH""",
)

entry(
    index = 2650,
    label = "C6H103OOH2-1 => C4H71-4 + CH2CHO + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH2-1 => C4H71-4 + CH2CHO + OH""",
)

entry(
    index = 2651,
    label = "C6H101OOH3-4 + O2 <=> C6H101OOH3-4O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-4 + O2 <=> C6H101OOH3-4O2""",
)

entry(
    index = 2652,
    label = "C6H101OOH3-5 + O2 <=> C6H101OOH3-5O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-5 + O2 <=> C6H101OOH3-5O2""",
)

entry(
    index = 2653,
    label = "C6H101OOH3-6 + O2 <=> C6H101OOH3-6O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-6 + O2 <=> C6H101OOH3-6O2""",
)

entry(
    index = 2654,
    label = "C6H101OOH4-3 + O2 <=> C6H101OOH4-3O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-3 + O2 <=> C6H101OOH4-3O2""",
)

entry(
    index = 2655,
    label = "C6H101OOH4-5 + O2 <=> C6H101OOH4-5O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-5 + O2 <=> C6H101OOH4-5O2""",
)

entry(
    index = 2656,
    label = "C6H101OOH4-6 + O2 <=> C6H101OOH4-6O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-6 + O2 <=> C6H101OOH4-6O2""",
)

entry(
    index = 2657,
    label = "C6H101OOH5-3 + O2 <=> C6H101OOH5-3O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-3 + O2 <=> C6H101OOH5-3O2""",
)

entry(
    index = 2658,
    label = "C6H101OOH5-4 + O2 <=> C6H101OOH5-4O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-4 + O2 <=> C6H101OOH5-4O2""",
)

entry(
    index = 2659,
    label = "C6H101OOH5-6 + O2 <=> C6H101OOH5-6O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-6 + O2 <=> C6H101OOH5-6O2""",
)

entry(
    index = 2660,
    label = "C6H101OOH6-3 + O2 <=> C6H101OOH6-3O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-3 + O2 <=> C6H101OOH6-3O2""",
)

entry(
    index = 2661,
    label = "C6H101OOH6-4 + O2 <=> C6H101OOH6-4O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-4 + O2 <=> C6H101OOH6-4O2""",
)

entry(
    index = 2662,
    label = "C6H101OOH6-5 + O2 <=> C6H101OOH6-5O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-5 + O2 <=> C6H101OOH6-5O2""",
)

entry(
    index = 2663,
    label = "C6H102OOH4-5 + O2 <=> C6H102OOH4-5O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-5 + O2 <=> C6H102OOH4-5O2""",
)

entry(
    index = 2664,
    label = "C6H102OOH4-6 + O2 <=> C6H102OOH4-6O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-6 + O2 <=> C6H102OOH4-6O2""",
)

entry(
    index = 2665,
    label = "C6H102OOH5-4 + O2 <=> C6H102OOH5-4O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-4 + O2 <=> C6H102OOH5-4O2""",
)

entry(
    index = 2666,
    label = "C6H102OOH5-6 + O2 <=> C6H102OOH5-6O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-6 + O2 <=> C6H102OOH5-6O2""",
)

entry(
    index = 2667,
    label = "C6H102OOH6-4 + O2 <=> C6H102OOH6-4O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-4 + O2 <=> C6H102OOH6-4O2""",
)

entry(
    index = 2668,
    label = "C6H102OOH6-5 + O2 <=> C6H102OOH6-5O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(7.54e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-5 + O2 <=> C6H102OOH6-5O2""",
)

entry(
    index = 2669,
    label = "C6H103OOH2-1 + O2 <=> C6H103OOH2-1O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4.52e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH2-1 + O2 <=> C6H103OOH2-1O2""",
)

entry(
    index = 2670,
    label = "C6H103OOH1-2 + O2 <=> C6H103OOH1-2O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH1-2 + O2 <=> C6H103OOH1-2O2""",
)

entry(
    index = 2671,
    label = "C6H101OOH3-4O2 => NC6D1KET34 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(21450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-4O2 => NC6D1KET34 + OH""",
)

entry(
    index = 2672,
    label = "C6H101OOH3-5O2 => NC6D1KET35 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(15450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-5O2 => NC6D1KET35 + OH""",
)

entry(
    index = 2673,
    label = "C6H101OOH3-6O2 => NC6D1KET36 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.563e+09, 's^-1'), n=0, Ea=(13650, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-6O2 => NC6D1KET36 + OH""",
)

entry(
    index = 2674,
    label = "C6H101OOH4-3O2 => NC6D1KET43 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-3O2 => NC6D1KET43 + OH""",
)

entry(
    index = 2675,
    label = "C6H101OOH4-5O2 => NC6D1KET45 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-5O2 => NC6D1KET45 + OH""",
)

entry(
    index = 2676,
    label = "C6H101OOH4-6O2 => NC6D1KET46 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH4-6O2 => NC6D1KET46 + OH""",
)

entry(
    index = 2677,
    label = "C6H101OOH5-3O2 => NC6D1KET53 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(17450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-3O2 => NC6D1KET53 + OH""",
)

entry(
    index = 2678,
    label = "C6H101OOH5-4O2 => NC6D1KET54 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-4O2 => NC6D1KET54 + OH""",
)

entry(
    index = 2679,
    label = "C6H101OOH5-6O2 => NC6D1KET56 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH5-6O2 => NC6D1KET56 + OH""",
)

entry(
    index = 2680,
    label = "C6H101OOH6-3O2 => NC6D1KET63 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3.12e+09, 's^-1'), n=0, Ea=(22150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-3O2 => NC6D1KET63 + OH""",
)

entry(
    index = 2681,
    label = "C6H101OOH6-4O2 => NC6D1KET64 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-4O2 => NC6D1KET64 + OH""",
)

entry(
    index = 2682,
    label = "C6H101OOH6-5O2 => NC6D1KET65 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26150, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH6-5O2 => NC6D1KET65 + OH""",
)

entry(
    index = 2683,
    label = "C6H102OOH4-5O2 => NC6D2KET45 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(21450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-5O2 => NC6D2KET45 + OH""",
)

entry(
    index = 2684,
    label = "C6H102OOH4-6O2 => NC6D2KET46 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.25e+10, 's^-1'), n=0, Ea=(15450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH4-6O2 => NC6D2KET46 + OH""",
)

entry(
    index = 2685,
    label = "C6H102OOH5-4O2 => NC6D2KET54 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-4O2 => NC6D2KET54 + OH""",
)

entry(
    index = 2686,
    label = "C6H102OOH5-6O2 => NC6D2KET56 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(23450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH5-6O2 => NC6D2KET56 + OH""",
)

entry(
    index = 2687,
    label = "C6H102OOH6-4O2 => NC6D2KET64 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2.5e+10, 's^-1'), n=0, Ea=(21000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-4O2 => NC6D2KET64 + OH""",
)

entry(
    index = 2688,
    label = "C6H102OOH6-5O2 => NC6D2KET65 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H102OOH6-5O2 => NC6D2KET65 + OH""",
)

entry(
    index = 2689,
    label = "C6H103OOH2-1O2 => NC6D3KET21 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(21450, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH2-1O2 => NC6D3KET21 + OH""",
)

entry(
    index = 2690,
    label = "C6H103OOH1-2O2 => NC6D3KET12 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 's^-1'), n=0, Ea=(26000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H103OOH1-2O2 => NC6D3KET12 + OH""",
)

entry(
    index = 2691,
    label = "NC6D1KET34 => OH + C2H3 + CO + C2H5CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET34 => OH + C2H3 + CO + C2H5CHO""",
)

entry(
    index = 2692,
    label = "NC6D1KET35 => OH + CH2CHO + C2H3COCH3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+15, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET35 => OH + CH2CHO + C2H3COCH3""",
)

entry(
    index = 2693,
    label = "NC6D1KET35 => OH + CH3 + CO + C2H3COCH3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+15, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET35 => OH + CH3 + CO + C2H3COCH3""",
)

entry(
    index = 2694,
    label = "NC6D1KET36 => OH + C2H3 + CO + CH2O + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET36 => OH + C2H3 + CO + CH2O + C2H4""",
)

entry(
    index = 2695,
    label = "NC6D1KET43 => OH + C3H5O + C2H3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET43 => OH + C3H5O + C2H3CHO""",
)

entry(
    index = 2696,
    label = "NC6D1KET45 => OH + C2H3 + CH3CHO + CH2CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET45 => OH + C2H3 + CH3CHO + CH2CO""",
)

entry(
    index = 2697,
    label = "NC6D1KET46 => OH + C3H5-A + CO + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET46 => OH + C3H5-A + CO + CH3CHO""",
)

entry(
    index = 2698,
    label = "NC6D1KET53 => OH + C2H5CO + C2H3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET53 => OH + C2H5CO + C2H3CHO""",
)

entry(
    index = 2699,
    label = "NC6D1KET54 => OH + C3H5-A + CO + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET54 => OH + C3H5-A + CO + CH3CHO""",
)

entry(
    index = 2700,
    label = "NC6D1KET56 => OH + C3H5-A + CH2CO + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET56 => OH + C3H5-A + CH2CO + CH2O""",
)

entry(
    index = 2701,
    label = "NC6D1KET63 => OH + C2H3CHO + HCO + C2H4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET63 => OH + C2H3CHO + HCO + C2H4""",
)

entry(
    index = 2702,
    label = "NC6D1KET64 => OH + C3H5-A + CH3CHO + CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+15, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET64 => OH + C3H5-A + CH3CHO + CO""",
)

entry(
    index = 2703,
    label = "NC6D1KET64 => OH + CH2CHO + AC3H5CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(5e+15, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET64 => OH + CH2CHO + AC3H5CHO""",
)

entry(
    index = 2704,
    label = "NC6D1KET65 => OH + HCO + CH2CO + C3H6",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D1KET65 => OH + HCO + CH2CO + C3H6""",
)

entry(
    index = 2705,
    label = "NC6D2KET45 => OH + C3H5-T + CO + CH3CHO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D2KET45 => OH + C3H5-T + CO + CH3CHO""",
)

entry(
    index = 2706,
    label = "NC6D2KET46 => OH + C3H5-T + CH2CO + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D2KET46 => OH + C3H5-T + CH2CO + CH2O""",
)

entry(
    index = 2707,
    label = "NC6D2KET54 => OH + CH2CHO + C2H3COCH3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D2KET54 => OH + CH2CHO + C2H3COCH3""",
)

entry(
    index = 2708,
    label = "NC6D2KET56 => OH + C3H5-T + CH2O + CH2CO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D2KET56 => OH + C3H5-T + CH2O + CH2CO""",
)

entry(
    index = 2709,
    label = "NC6D2KET64 => OH + CH2CHO + C2H3COCH3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D2KET64 => OH + CH2CHO + C2H3COCH3""",
)

entry(
    index = 2710,
    label = "NC6D2KET65 => OH + C4H71-3 + CO + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D2KET65 => OH + C4H71-3 + CO + CH2O""",
)

entry(
    index = 2711,
    label = "NC6D3KET21 => OH + C4H71-4 + CO + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D3KET21 => OH + C4H71-4 + CO + CH2O""",
)

entry(
    index = 2712,
    label = "NC6D3KET12 => OH + HCO + C4H7CHO1-4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is NC6D3KET12 => OH + HCO + C4H7CHO1-4""",
)

entry(
    index = 2713,
    label = "C6H12OH-1 => CH3CHO + PC4H9",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-1 => CH3CHO + PC4H9""",
)

entry(
    index = 2714,
    label = "C6H12OH-1 => C2H5CHO + NC3H7",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-1 => C2H5CHO + NC3H7""",
)

entry(
    index = 2715,
    label = "C6H12OH-2 => NC3H7CHO + C2H5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-2 => NC3H7CHO + C2H5""",
)

entry(
    index = 2716,
    label = "C6H12OH-2 => NC4H9CHO + CH3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1.5e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-2 => NC4H9CHO + CH3""",
)

entry(
    index = 2717,
    label = "C6H12OH-3 => NC3H7CHO + C2H5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(3e+13, 's^-1'), n=0, Ea=(30000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12OH-3 => NC3H7CHO + C2H5""",
)

entry(
    index = 2718,
    label = "C6H12-3 + O <=> NC3H7 + C2H5CO",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H12-3 + O <=> NC3H7 + C2H5CO""",
)

entry(
    index = 2719,
    label = "C6H111-4 + HO2 <=> C6H111O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-4 + HO2 <=> C6H111O2H-4""",
)

entry(
    index = 2720,
    label = "C6H111-5 + HO2 <=> C6H111O2H-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-5 + HO2 <=> C6H111O2H-5""",
)

entry(
    index = 2721,
    label = "C6H111-6 + HO2 <=> C6H111O2H-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111-6 + HO2 <=> C6H111O2H-6""",
)

entry(
    index = 2722,
    label = "C6H112-5 + HO2 <=> C6H112O2H-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-5 + HO2 <=> C6H112O2H-5""",
)

entry(
    index = 2723,
    label = "C6H112-6 + HO2 <=> C6H112O2H-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-6 + HO2 <=> C6H112O2H-6""",
)

entry(
    index = 2724,
    label = "C6H112-4 + HO2 <=> C6H113O2H-2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + HO2 <=> C6H113O2H-2""",
)

entry(
    index = 2725,
    label = "C6H112-4 + HO2 <=> C6H112O2H-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(4e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112-4 + HO2 <=> C6H112O2H-4""",
)

entry(
    index = 2726,
    label = "C6H113-1 + HO2 <=> C6H113O2H-1",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113-1 + HO2 <=> C6H113O2H-1""",
)

entry(
    index = 2727,
    label = "C6H112O2-1 => C2H3COC3H7 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 's^-1'), n=0, Ea=(38000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2-1 => C2H3COC3H7 + OH""",
)

entry(
    index = 2728,
    label = "C6H111O2-3 => C2H3COC3H7 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4e+12, 's^-1'), n=0, Ea=(38000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2-3 => C2H3COC3H7 + OH""",
)

entry(
    index = 2729,
    label = "C6H101OOH3-6 => C2H3COC3H7 + OH",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+11, 's^-1'), n=0, Ea=(16700, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101OOH3-6 => C2H3COC3H7 + OH""",
)

entry(
    index = 2730,
    label = "C6H111O2H-4 => OH + C2H5CHO + C3H5-A",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2H-4 => OH + C2H5CHO + C3H5-A""",
)

entry(
    index = 2731,
    label = "C6H111O2H-5 => OH + CH3CHO + C4H71-4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2H-5 => OH + CH3CHO + C4H71-4""",
)

entry(
    index = 2732,
    label = "C6H111O2H-6 => OH + CH2O + C5H91-5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H111O2H-6 => OH + CH2O + C5H91-5""",
)

entry(
    index = 2733,
    label = "C6H112O2H-4 => OH + C2H3COCH3 + C2H5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2H-4 => OH + C2H3COCH3 + C2H5""",
)

entry(
    index = 2734,
    label = "C6H112O2H-5 => OH + CH3CHO + C4H71-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2H-5 => OH + CH3CHO + C4H71-3""",
)

entry(
    index = 2735,
    label = "C6H112O2H-6 => OH + CH2O + C5H91-5",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H112O2H-6 => OH + CH2O + C5H91-5""",
)

entry(
    index = 2736,
    label = "C6H113O2H-2 => OH + CH3CHO + C4H71-4",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113O2H-2 => OH + CH3CHO + C4H71-4""",
)

entry(
    index = 2737,
    label = "C6H113O2H-1 => OH + CH2O + C5H91-3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+16, 's^-1'), n=0, Ea=(39000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H113O2H-1 => OH + CH2O + C5H91-3""",
)

entry(
    index = 2738,
    label = "O2 + ETES1 => HO2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.045e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (40722.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O2 + ETES1 => HO2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2739,
    label = "H + ETES1 => H2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.574e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (3950.57, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + ETES1 => H2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2740,
    label = "OH + ETES1 => H2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (4.793e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-2259.83, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + ETES1 => H2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2741,
    label = "O + ETES1 => OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.624e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (2579.54, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + ETES1 => OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2742,
    label = "HO2 + ETES1 => H2O2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(646500, 'cm^3/(mol*s)'), n=2, Ea=(11887.7, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + ETES1 => H2O2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2743,
    label = "HCO + ETES1 => CH2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.516e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (12360.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + ETES1 => CH2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2744,
    label = "CH3 + ETES1 => CH4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(468400, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + ETES1 => CH4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2745,
    label = "C2H5 + ETES1 => C2H6 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(275800, 'cm^3/(mol*s)'), n=2, Ea=(7658.07, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + ETES1 => C2H6 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2746,
    label = "CH3O + ETES1 => CH3OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(513600, 'cm^3/(mol*s)'), n=2, Ea=(1583.56, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + ETES1 => CH3OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2747,
    label = "CH3O2 + ETES1 => CH3O2H + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(913300, 'cm^3/(mol*s)'), n=2, Ea=(12360.4, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + ETES1 => CH3O2H + C2H3COCH3 + C2H3""",
)

entry(
    index = 2748,
    label = "O2 + MVOX => HO2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.045e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (40722.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O2 + MVOX => HO2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2749,
    label = "H + MVOX => H2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.574e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (3950.57, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + MVOX => H2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2750,
    label = "OH + MVOX => H2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (4.793e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-2259.83, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + MVOX => H2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2751,
    label = "O + MVOX => OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.624e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (2579.54, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + MVOX => OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2752,
    label = "HO2 + MVOX => H2O2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(646500, 'cm^3/(mol*s)'), n=2, Ea=(11887.7, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + MVOX => H2O2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2753,
    label = "HCO + MVOX => CH2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.516e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (12360.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + MVOX => CH2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2754,
    label = "CH3 + MVOX => CH4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(468400, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + MVOX => CH4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2755,
    label = "C2H5 + MVOX => C2H6 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(275800, 'cm^3/(mol*s)'), n=2, Ea=(7658.07, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + MVOX => C2H6 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2756,
    label = "CH3O2 + MVOX => CH3O2H + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(913300, 'cm^3/(mol*s)'), n=2, Ea=(12360.4, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + MVOX => CH3O2H + C2H3COCH3 + C2H3""",
)

entry(
    index = 2757,
    label = "C2H3 + MVOX => C2H4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(813900, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + MVOX => C2H4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2758,
    label = "CH3O + MVOX => CH3OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(513600, 'cm^3/(mol*s)'), n=2, Ea=(1583.56, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + MVOX => CH3OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2759,
    label = "O2 + VTHF => HO2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.045e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (40722.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O2 + VTHF => HO2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2760,
    label = "H + VTHF => H2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.574e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (3950.57, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + VTHF => H2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2761,
    label = "OH + VTHF => H2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (4.793e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-2259.83, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + VTHF => H2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2762,
    label = "O + VTHF => OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.624e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (2579.54, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + VTHF => OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2763,
    label = "HO2 + VTHF => H2O2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(646500, 'cm^3/(mol*s)'), n=2, Ea=(11887.7, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + VTHF => H2O2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2764,
    label = "HCO + VTHF => CH2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.516e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (12360.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + VTHF => CH2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2765,
    label = "CH3 + VTHF => CH4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(468400, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + VTHF => CH4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2766,
    label = "C2H5 + VTHF => C2H6 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(275800, 'cm^3/(mol*s)'), n=2, Ea=(7658.07, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + VTHF => C2H6 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2767,
    label = "CH3O + VTHF => CH3OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(513600, 'cm^3/(mol*s)'), n=2, Ea=(1583.56, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + VTHF => CH3OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2768,
    label = "CH3O2 + VTHF => CH3O2H + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(913300, 'cm^3/(mol*s)'), n=2, Ea=(12360.4, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + VTHF => CH3O2H + C2H3COCH3 + C2H3""",
)

entry(
    index = 2769,
    label = "C2H3 + VTHF => C2H4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(813900, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + VTHF => C2H4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2770,
    label = "O2 + EDHF => HO2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.045e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (40722.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O2 + EDHF => HO2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2771,
    label = "H + EDHF => H2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (2.574e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (3950.57, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + EDHF => H2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2772,
    label = "OH + EDHF => H2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (4.793e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-2259.83, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + EDHF => H2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2773,
    label = "O + EDHF => OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.624e+07, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (2579.54, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + EDHF => OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2774,
    label = "HO2 + EDHF => H2O2 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(646500, 'cm^3/(mol*s)'), n=2, Ea=(11887.7, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + EDHF => H2O2 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2775,
    label = "HCO + EDHF => CH2O + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.516e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (12360.4, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is HCO + EDHF => CH2O + C2H3COCH3 + C2H3""",
)

entry(
    index = 2776,
    label = "CH3 + EDHF => CH4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(468400, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + EDHF => CH4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2777,
    label = "C2H5 + EDHF => C2H6 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(275800, 'cm^3/(mol*s)'), n=2, Ea=(7658.07, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + EDHF => C2H6 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2778,
    label = "C2H3 + EDHF => C2H4 + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(813900, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + EDHF => C2H4 + C2H3COCH3 + C2H3""",
)

entry(
    index = 2779,
    label = "CH3O2 + EDHF => CH3O2H + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(913300, 'cm^3/(mol*s)'), n=2, Ea=(12360.4, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + EDHF => CH3O2H + C2H3COCH3 + C2H3""",
)

entry(
    index = 2780,
    label = "CH3O + EDHF => CH3OH + C2H3COCH3 + C2H3",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(513600, 'cm^3/(mol*s)'), n=2, Ea=(1583.56, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + EDHF => CH3OH + C2H3COCH3 + C2H3""",
)

entry(
    index = 2781,
    label = "O2 + C5H81-4 => HO2 + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (5.111e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (40722.5, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O2 + C5H81-4 => HO2 + C3H5-A + C2H2""",
)

entry(
    index = 2782,
    label = "H + C5H81-4 => H2 + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (6.435e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (3950.57, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is H + C5H81-4 => H2 + C3H5-A + C2H2""",
)

entry(
    index = 2783,
    label = "OH + C5H81-4 => H2O + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (1.198e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (-2259.83, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is OH + C5H81-4 => H2O + C3H5-A + C2H2""",
)

entry(
    index = 2784,
    label = "O + C5H81-4 => OH + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(
        A = (4.06e+06, 'cm^3/(mol*s)'),
        n = 2,
        Ea = (2579.54, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is O + C5H81-4 => OH + C3H5-A + C2H2""",
)

entry(
    index = 2785,
    label = "HO2 + C5H81-4 => H2O2 + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(161600, 'cm^3/(mol*s)'), n=2, Ea=(11887.7, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HO2 + C5H81-4 => H2O2 + C3H5-A + C2H2""",
)

entry(
    index = 2786,
    label = "HCO + C5H81-4 => CH2O + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(378900, 'cm^3/(mol*s)'), n=2, Ea=(12360.4, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is HCO + C5H81-4 => CH2O + C3H5-A + C2H2""",
)

entry(
    index = 2787,
    label = "CH3 + C5H81-4 => CH4 + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(117100, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3 + C5H81-4 => CH4 + C3H5-A + C2H2""",
)

entry(
    index = 2788,
    label = "C2H5 + C5H81-4 => C2H6 + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(68950, 'cm^3/(mol*s)'), n=2, Ea=(7658.07, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H5 + C5H81-4 => C2H6 + C3H5-A + C2H2""",
)

entry(
    index = 2789,
    label = "CH3O + C5H81-4 => CH3OH + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(128400, 'cm^3/(mol*s)'), n=2, Ea=(1583.56, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O + C5H81-4 => CH3OH + C3H5-A + C2H2""",
)

entry(
    index = 2790,
    label = "CH3O2 + C5H81-4 => CH3O2H + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(228300, 'cm^3/(mol*s)'), n=2, Ea=(12360.4, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is CH3O2 + C5H81-4 => CH3O2H + C3H5-A + C2H2""",
)

entry(
    index = 2791,
    label = "C2H3 + C5H81-4 => C2H4 + C3H5-A + C2H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(203500, 'cm^3/(mol*s)'), n=2, Ea=(4871.29, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C2H3 + C5H81-4 => C2H4 + C3H5-A + C2H2""",
)

entry(
    index = 2792,
    label = "C6H101-3 + H <=> C6H111-3",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + H <=> C6H111-3""",
)

entry(
    index = 2793,
    label = "C6H101-3 + H <=> C6H111-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + H <=> C6H111-4""",
)

entry(
    index = 2794,
    label = "C6H101-3 + H <=> C6H113-1",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.25e+11, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (1230, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + H <=> C6H113-1""",
)

entry(
    index = 2795,
    label = "C6H101-3 + H <=> C6H112-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + H <=> C6H112-4""",
)

entry(
    index = 2796,
    label = "C6H101-4 + H <=> C6H111-4",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + H <=> C6H111-4""",
)

entry(
    index = 2797,
    label = "C6H101-4 + H <=> C6H111-5",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + H <=> C6H111-5""",
)

entry(
    index = 2798,
    label = "C6H101-4 + H <=> C6H112-5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.25e+11, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (1230, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + H <=> C6H112-5""",
)

entry(
    index = 2799,
    label = "C6H101-4 + H <=> C6H112-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + H <=> C6H112-6""",
)

entry(
    index = 2800,
    label = "C6H101-5 + H <=> C6H111-5",
    degeneracy = 1,
    kinetics = Arrhenius(
        A = (4.25e+11, 'cm^3/(mol*s)'),
        n = 0.5,
        Ea = (1230, 'cal/mol'),
        T0 = (1, 'K'),
    ),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + H <=> C6H111-5""",
)

entry(
    index = 2801,
    label = "C6H101-5 + H <=> C6H111-6",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.5e+11, 'cm^3/(mol*s)'), n=0.5, Ea=(2620, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + H <=> C6H111-6""",
)

entry(
    index = 2802,
    label = "C6H101-3 + O => C5H91-3 + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + O => C5H91-3 + HCO""",
)

entry(
    index = 2803,
    label = "C6H101-4 + O => C5H92-5 + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + O => C5H92-5 + HCO""",
)

entry(
    index = 2804,
    label = "C6H101-5 + O => C5H91-5 + HCO",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(2e+11, 'cm^3/(mol*s)'), n=0, Ea=(-1050, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + O => C5H91-5 + HCO""",
)

entry(
    index = 2805,
    label = "C6H101-3 + OH => C5H91-3 + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + OH => C5H91-3 + CH2O""",
)

entry(
    index = 2806,
    label = "C6H101-4 + OH => C5H92-5 + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + OH => C5H92-5 + CH2O""",
)

entry(
    index = 2807,
    label = "C6H101-5 + OH => C5H91-5 + CH2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(1e+12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + OH => C5H91-5 + CH2O""",
)

entry(
    index = 2808,
    label = "C6H101-3 + OH => C2H3 + C4H6 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + OH => C2H3 + C4H6 + H2O""",
)

entry(
    index = 2809,
    label = "C6H101-4 + OH => C2H3 + C4H6 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + OH => C2H3 + C4H6 + H2O""",
)

entry(
    index = 2810,
    label = "C6H101-5 + OH => C2H3 + C4H6 + H2O",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(27640, 'cm^3/(mol*s)'), n=2.64, Ea=(-1919, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + OH => C2H3 + C4H6 + H2O""",
)

entry(
    index = 2811,
    label = "C6H101-3 + HO2 => C2H3 + C4H6 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + HO2 => C2H3 + C4H6 + H2O2""",
)

entry(
    index = 2812,
    label = "C6H101-4 + HO2 => C2H3 + C4H6 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + HO2 => C2H3 + C4H6 + H2O2""",
)

entry(
    index = 2813,
    label = "C6H101-5 + HO2 => C2H3 + C4H6 + H2O2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(4820, 'cm^3/(mol*s)'), n=2.55, Ea=(10530, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + HO2 => C2H3 + C4H6 + H2O2""",
)

entry(
    index = 2814,
    label = "C6H101-3 + H => C2H3 + C4H6 + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-3 + H => C2H3 + C4H6 + H2""",
)

entry(
    index = 2815,
    label = "C6H101-4 + H => C2H3 + C4H6 + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-4 + H => C2H3 + C4H6 + H2""",
)

entry(
    index = 2816,
    label = "C6H101-5 + H => C2H3 + C4H6 + H2",
    degeneracy = 1,
    reversible = False,
    kinetics = Arrhenius(A=(337600, 'cm^3/(mol*s)'), n=2.36, Ea=(207, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""The chemkin file reaction is C6H101-5 + H => C2H3 + C4H6 + H2""",
)
