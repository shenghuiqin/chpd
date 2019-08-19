species(
    label = '[CH2]C(C)C(C)[C]=C(19128)',
    structure = SMILES('[CH2]C(C)C(C)[C]=C'),
    E0 = (341.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,366.246,366.56],'cm^-1')),
        HinderedRotor(inertia=(0.0609728,'amu*angstrom^2'), symmetry=1, barrier=(5.82104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00125472,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610627,'amu*angstrom^2'), symmetry=1, barrier=(5.82248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.061059,'amu*angstrom^2'), symmetry=1, barrier=(5.82062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256176,'amu*angstrom^2'), symmetry=1, barrier=(24.4238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479762,0.068443,-3.63188e-05,3.14928e-09,2.78028e-12,41194.8,31.1779], Tmin=(100,'K'), Tmax=(1065.73,'K')), NASAPolynomial(coeffs=[12.6904,0.0376393,-1.4113e-05,2.48885e-09,-1.68411e-13,37738.7,-32.5125], Tmin=(1065.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C3H6(72)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,-487.138,-4.54468], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'H(3)',
    structure = SMILES('[H]'),
    E0 = (211.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25472.7,-0.459566], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=[C]C(C)C(=C)C(24257)',
    structure = SMILES('C=[C]C(C)C(=C)C'),
    E0 = (253.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,283.157,283.339],'cm^-1')),
        HinderedRotor(inertia=(0.180377,'amu*angstrom^2'), symmetry=1, barrier=(10.2792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180653,'amu*angstrom^2'), symmetry=1, barrier=(10.2792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180637,'amu*angstrom^2'), symmetry=1, barrier=(10.2799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180655,'amu*angstrom^2'), symmetry=1, barrier=(10.279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.945107,0.0708261,-5.52303e-05,2.52471e-08,-5.13744e-12,30652.1,24.6123], Tmin=(100,'K'), Tmax=(1105.31,'K')), NASAPolynomial(coeffs=[8.2329,0.0444522,-1.94386e-05,3.65931e-09,-2.54663e-13,29041,-11.2768], Tmin=(1105.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)C(C)=C=C(24258)',
    structure = SMILES('[CH2]C(C)C(C)=C=C'),
    E0 = (258.332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,315.045],'cm^-1')),
        HinderedRotor(inertia=(0.116865,'amu*angstrom^2'), symmetry=1, barrier=(8.23005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116912,'amu*angstrom^2'), symmetry=1, barrier=(8.23275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116827,'amu*angstrom^2'), symmetry=1, barrier=(8.23016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.422467,'amu*angstrom^2'), symmetry=1, barrier=(29.7604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.360798,0.0716729,-5.22077e-05,2.05829e-08,-3.34596e-12,31208.1,28.5446], Tmin=(100,'K'), Tmax=(1440.07,'K')), NASAPolynomial(coeffs=[14.1578,0.0333498,-1.22898e-05,2.10331e-09,-1.37861e-13,27234.3,-43.0497], Tmin=(1440.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(C)C([CH2])C(24259)',
    structure = SMILES('C#CC(C)C([CH2])C'),
    E0 = (268.504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,369.695],'cm^-1')),
        HinderedRotor(inertia=(0.0853198,'amu*angstrom^2'), symmetry=1, barrier=(8.2503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818101,'amu*angstrom^2'), symmetry=1, barrier=(78.9858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00124238,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853722,'amu*angstrom^2'), symmetry=1, barrier=(8.24652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818485,'amu*angstrom^2'), symmetry=1, barrier=(78.9851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.624277,0.0660428,-3.15507e-05,-6.31353e-09,8.07105e-12,32422.5,27.9924], Tmin=(100,'K'), Tmax=(939.752,'K')), NASAPolynomial(coeffs=[13.0084,0.0333095,-1.1193e-05,1.8676e-09,-1.23705e-13,29212.7,-35.6784], Tmin=(939.752,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=[C]C(C)C=C(24201)',
    structure = SMILES('C=[C]C(C)C=C'),
    E0 = (293.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3010,987.5,1337.5,450,1655,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,272.335,275.352],'cm^-1')),
        HinderedRotor(inertia=(0.2164,'amu*angstrom^2'), symmetry=1, barrier=(11.4749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214712,'amu*angstrom^2'), symmetry=1, barrier=(11.4778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213609,'amu*angstrom^2'), symmetry=1, barrier=(11.4757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55688,0.0537055,-3.52508e-05,1.22891e-08,-1.83166e-12,35329.9,21.6469], Tmin=(100,'K'), Tmax=(1475.3,'K')), NASAPolynomial(coeffs=[9.43814,0.032337,-1.35245e-05,2.47131e-09,-1.67972e-13,33004.4,-19.4404], Tmin=(1475.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C3H6(T)(143)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238389,'amu*angstrom^2'), symmetry=1, barrier=(5.48103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26842e-06,-1.44873e-08,5.74941e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.81,'K')), NASAPolynomial(coeffs=[5.93909,0.0171892,-6.69152e-06,1.21546e-09,-8.39795e-14,33151.2,-4.14888], Tmin=(1046.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C)C=C=C(24260)',
    structure = SMILES('[CH2]C(C)C=C=C'),
    E0 = (297.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,338.333],'cm^-1')),
        HinderedRotor(inertia=(0.112667,'amu*angstrom^2'), symmetry=1, barrier=(9.15834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112727,'amu*angstrom^2'), symmetry=1, barrier=(9.15333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361205,'amu*angstrom^2'), symmetry=1, barrier=(29.3265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35132,0.0500663,-1.67046e-05,-1.18626e-08,7.86829e-12,35870,24.2252], Tmin=(100,'K'), Tmax=(992.295,'K')), NASAPolynomial(coeffs=[11.2642,0.0279033,-1.01035e-05,1.77617e-09,-1.21375e-13,33026.5,-27.9366], Tmin=(992.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C(C)[C](C)C(24261)',
    structure = SMILES('C=[C]C(C)[C](C)C'),
    E0 = (321.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25184,0.0658973,-4.17849e-05,1.51318e-08,-2.57135e-12,38790,27.0662], Tmin=(100,'K'), Tmax=(1206.15,'K')), NASAPolynomial(coeffs=[6.09294,0.0498426,-2.18188e-05,4.09605e-09,-2.83955e-13,37622.2,2.80334], Tmin=(1206.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)[C](C)C=C(19125)',
    structure = SMILES('[CH2]C(C)[C](C)C=C'),
    E0 = (235.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,360,370,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,3590.93,3591.31],'cm^-1')),
        HinderedRotor(inertia=(0.0847374,'amu*angstrom^2'), symmetry=1, barrier=(84.2868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866637,'amu*angstrom^2'), symmetry=1, barrier=(17.6046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765841,'amu*angstrom^2'), symmetry=1, barrier=(17.6082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.6666,'amu*angstrom^2'), symmetry=1, barrier=(84.3023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.66569,'amu*angstrom^2'), symmetry=1, barrier=(84.2815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852679,0.0577386,-3.72909e-06,-3.16845e-08,1.54818e-11,28451.7,27.6176], Tmin=(100,'K'), Tmax=(981.615,'K')), NASAPolynomial(coeffs=[12.1951,0.037873,-1.36438e-05,2.39936e-09,-1.64405e-13,24955.3,-33.3599], Tmin=(981.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_T)"""),
)

species(
    label = '[CH]=CC(C)C([CH2])C(19130)',
    structure = SMILES('[CH]=CC(C)C([CH2])C'),
    E0 = (350.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,395.955],'cm^-1')),
        HinderedRotor(inertia=(0.0737522,'amu*angstrom^2'), symmetry=1, barrier=(8.20663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737676,'amu*angstrom^2'), symmetry=1, barrier=(8.20739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0010752,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737811,'amu*angstrom^2'), symmetry=1, barrier=(8.20712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243539,'amu*angstrom^2'), symmetry=1, barrier=(27.0943,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.48393,0.0668613,-2.65614e-05,-1.06134e-08,8.51368e-12,42309.2,31.091], Tmin=(100,'K'), Tmax=(1002.84,'K')), NASAPolynomial(coeffs=[13.7533,0.0359125,-1.31437e-05,2.32069e-09,-1.58743e-13,38542.6,-38.4735], Tmin=(1002.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][C](C)C(C)C(24262)',
    structure = SMILES('C=[C][C](C)C(C)C'),
    E0 = (268.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.750608,0.0612011,-1.76528e-05,-1.32307e-08,7.61486e-12,32394,25.9861], Tmin=(100,'K'), Tmax=(1069.6,'K')), NASAPolynomial(coeffs=[11.783,0.0398118,-1.55206e-05,2.80756e-09,-1.93081e-13,28897.4,-33.2941], Tmin=(1069.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C](C)C(C)C=C(19127)',
    structure = SMILES('[CH2][C](C)C(C)C=C'),
    E0 = (288.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00661,0.0665493,-4.22222e-05,1.50731e-08,-2.39043e-12,34862.6,29.9408], Tmin=(100,'K'), Tmax=(1360.66,'K')), NASAPolynomial(coeffs=[8.44539,0.0446813,-1.81148e-05,3.26161e-09,-2.20258e-13,32838.3,-8.238], Tmin=(1360.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C([CH2])C=C(18971)',
    structure = SMILES('[CH2]C(C)C([CH2])C=C'),
    E0 = (308.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2654.03],'cm^-1')),
        HinderedRotor(inertia=(0.269031,'amu*angstrom^2'), symmetry=1, barrier=(77.7226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00281092,'amu*angstrom^2'), symmetry=1, barrier=(14.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611302,'amu*angstrom^2'), symmetry=1, barrier=(14.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611227,'amu*angstrom^2'), symmetry=1, barrier=(14.0533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137233,'amu*angstrom^2'), symmetry=1, barrier=(77.7265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3397.52,'J/mol'), sigma=(6.23225,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.69 K, Pc=31.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611931,0.064599,-2.09315e-05,-1.73992e-08,1.16317e-11,37251.2,32.01], Tmin=(100,'K'), Tmax=(956.077,'K')), NASAPolynomial(coeffs=[13.091,0.0357298,-1.2257e-05,2.08624e-09,-1.40241e-13,33798.3,-33.2126], Tmin=(956.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C(C)C=C(18362)',
    structure = SMILES('[CH2]C([CH2])C(C)C=C'),
    E0 = (308.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2654.03],'cm^-1')),
        HinderedRotor(inertia=(0.269031,'amu*angstrom^2'), symmetry=1, barrier=(77.7226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00281092,'amu*angstrom^2'), symmetry=1, barrier=(14.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611302,'amu*angstrom^2'), symmetry=1, barrier=(14.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611227,'amu*angstrom^2'), symmetry=1, barrier=(14.0533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137233,'amu*angstrom^2'), symmetry=1, barrier=(77.7265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.611931,0.064599,-2.09315e-05,-1.73992e-08,1.16317e-11,37251.2,31.3169], Tmin=(100,'K'), Tmax=(956.077,'K')), NASAPolynomial(coeffs=[13.091,0.0357298,-1.2257e-05,2.08624e-09,-1.40241e-13,33798.3,-33.9057], Tmin=(956.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([C]=C)C(C)C(19129)',
    structure = SMILES('[CH2]C([C]=C)C(C)C'),
    E0 = (341.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,366.246,366.56],'cm^-1')),
        HinderedRotor(inertia=(0.0609728,'amu*angstrom^2'), symmetry=1, barrier=(5.82104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00125472,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0610627,'amu*angstrom^2'), symmetry=1, barrier=(5.82248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.061059,'amu*angstrom^2'), symmetry=1, barrier=(5.82062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256176,'amu*angstrom^2'), symmetry=1, barrier=(24.4238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479762,0.068443,-3.63188e-05,3.14928e-09,2.78028e-12,41194.8,30.4848], Tmin=(100,'K'), Tmax=(1065.73,'K')), NASAPolynomial(coeffs=[12.6904,0.0376393,-1.4113e-05,2.48885e-09,-1.68411e-13,37738.7,-33.2056], Tmin=(1065.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=[C]C(C)C(C)C(24263)',
    structure = SMILES('[CH]=[C]C(C)C(C)C'),
    E0 = (383.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0180511,'amu*angstrom^2'), symmetry=1, barrier=(4.09088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.52733,'amu*angstrom^2'), symmetry=1, barrier=(12.1244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.527333,'amu*angstrom^2'), symmetry=1, barrier=(12.1244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0534996,'amu*angstrom^2'), symmetry=1, barrier=(12.1244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0535041,'amu*angstrom^2'), symmetry=1, barrier=(12.1242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.326143,0.0708978,-4.20859e-05,9.30529e-09,2.93826e-13,46253.9,29.6648], Tmin=(100,'K'), Tmax=(1179.89,'K')), NASAPolynomial(coeffs=[14.1064,0.0366302,-1.43484e-05,2.57556e-09,-1.75034e-13,42135.5,-42.7688], Tmin=(1179.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)C(=C)C(19134)',
    structure = SMILES('C=CC(C)C(=C)C'),
    E0 = (16.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60745,0.070627,-4.57496e-05,1.53761e-08,-2.1561e-12,2064.29,25.3447], Tmin=(100,'K'), Tmax=(1597.77,'K')), NASAPolynomial(coeffs=[13.4412,0.0384975,-1.55858e-05,2.79019e-09,-1.86787e-13,-2036.74,-42.5846], Tmin=(1597.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(C)C(C)C(24264)',
    structure = SMILES('C=C=C(C)C(C)C'),
    E0 = (53.2498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.385266,0.0694695,-3.78221e-05,5.54959e-09,1.3869e-12,6542.91,25.9834], Tmin=(100,'K'), Tmax=(1169.84,'K')), NASAPolynomial(coeffs=[13.7833,0.0376259,-1.49012e-05,2.69373e-09,-1.83911e-13,2452.44,-44.841], Tmin=(1169.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CH2(S)(23)',
    structure = SMILES('[CH2]'),
    E0 = (419.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.36,2789.41,2993.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(419.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]CC(C)[C]=C(19002)',
    structure = SMILES('[CH2]CC(C)[C]=C'),
    E0 = (371.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2140.01,4000],'cm^-1')),
        HinderedRotor(inertia=(0.226736,'amu*angstrom^2'), symmetry=1, barrier=(9.03751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.586158,'amu*angstrom^2'), symmetry=1, barrier=(23.3629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.586135,'amu*angstrom^2'), symmetry=1, barrier=(23.3629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.586207,'amu*angstrom^2'), symmetry=1, barrier=(23.3629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.55,'J/mol'), sigma=(5.8838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.98 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14098,0.0550014,-2.90653e-05,4.40999e-09,8.41902e-13,44734.5,26.8529], Tmin=(100,'K'), Tmax=(1211.95,'K')), NASAPolynomial(coeffs=[11.4531,0.0312433,-1.23797e-05,2.22819e-09,-1.51275e-13,41480.2,-27.9932], Tmin=(1211.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C)C[C]=C(24254)',
    structure = SMILES('[CH2]C(C)C[C]=C'),
    E0 = (369.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,731.424,2713.7],'cm^-1')),
        HinderedRotor(inertia=(0.553745,'amu*angstrom^2'), symmetry=1, barrier=(12.7317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.55533,'amu*angstrom^2'), symmetry=1, barrier=(12.7681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847991,'amu*angstrom^2'), symmetry=1, barrier=(1.9497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.49886,'amu*angstrom^2'), symmetry=1, barrier=(80.4458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980639,0.0578249,-3.78769e-05,1.33248e-08,-1.94065e-12,44592.6,27.3941], Tmin=(100,'K'), Tmax=(1585.57,'K')), NASAPolynomial(coeffs=[12.5497,0.0286387,-1.02656e-05,1.71532e-09,-1.10141e-13,40923.9,-33.7527], Tmin=(1585.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=[C]C(C)C[CH]C(19169)',
    structure = SMILES('C=[C]C(C)C[CH]C'),
    E0 = (336.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0597,0.0647849,-3.73125e-05,1.09539e-08,-1.35798e-12,40570.8,30.3163], Tmin=(100,'K'), Tmax=(1715.74,'K')), NASAPolynomial(coeffs=[11.094,0.0413914,-1.68606e-05,3.00709e-09,-2.00057e-13,37127.5,-23.5102], Tmin=(1715.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)[CH]CC(24265)',
    structure = SMILES('C=[C]C(C)[CH]CC'),
    E0 = (336.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751194,0.0665364,-3.88394e-05,1.11087e-08,-1.29102e-12,40598.1,31.271], Tmin=(100,'K'), Tmax=(1913.45,'K')), NASAPolynomial(coeffs=[15.9081,0.0348515,-1.40008e-05,2.45473e-09,-1.60346e-13,34797.7,-51.6879], Tmin=(1913.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C)C(=C)[CH]C(24172)',
    structure = SMILES('[CH2]C(=CC)C([CH2])C'),
    E0 = (237.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0358237,'amu*angstrom^2'), symmetry=1, barrier=(17.0825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.92092,'amu*angstrom^2'), symmetry=1, barrier=(90.1497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.742218,'amu*angstrom^2'), symmetry=1, barrier=(17.065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172372,'amu*angstrom^2'), symmetry=1, barrier=(3.96316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.93149,'amu*angstrom^2'), symmetry=1, barrier=(90.3926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3479.64,'J/mol'), sigma=(6.29859,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=543.51 K, Pc=31.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441935,0.0682285,-3.06307e-05,-5.44961e-09,6.38237e-12,28690.7,29.4733], Tmin=(100,'K'), Tmax=(1022.38,'K')), NASAPolynomial(coeffs=[13.4828,0.0369601,-1.37359e-05,2.43177e-09,-1.6594e-13,24991.8,-38.7786], Tmin=(1022.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)[CH]C(=C)C(24266)',
    structure = SMILES('[CH2]C(C)=CC([CH2])C'),
    E0 = (233.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441935,0.0682285,-3.06307e-05,-5.44961e-09,6.38237e-12,28187.5,29.4733], Tmin=(100,'K'), Tmax=(1022.38,'K')), NASAPolynomial(coeffs=[13.4828,0.0369601,-1.37359e-05,2.43177e-09,-1.6594e-13,24488.6,-38.7786], Tmin=(1022.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1CC(C)C1C(24267)',
    structure = SMILES('C=C1CC(C)C1C'),
    E0 = (35.7798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.99122,-0.0100132,0.000144281,-1.38728e-07,3.23306e-11,3980.95,-16.4526], Tmin=(100,'K'), Tmax=(1679.8,'K')), NASAPolynomial(coeffs=[69.3651,0.0437352,-7.79586e-05,1.86251e-08,-1.37959e-12,-43496.9,-415.639], Tmin=(1679.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.7798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane)"""),
)

species(
    label = 'CH2(19)',
    structure = SMILES('[CH2]'),
    E0 = (381.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1032.72,2936.3,3459],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8328,0.000224446,4.68033e-06,-6.04743e-09,2.59009e-12,45920.8,1.40666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16229,0.00281798,-7.56235e-07,5.05446e-11,5.65236e-15,46099.1,4.77656], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(381.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=[C]C(C)[CH]C(24204)',
    structure = SMILES('C=[C]C(C)[CH]C'),
    E0 = (360.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1922.94,4000],'cm^-1')),
        HinderedRotor(inertia=(0.132455,'amu*angstrom^2'), symmetry=1, barrier=(9.24278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364554,'amu*angstrom^2'), symmetry=1, barrier=(25.4103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368751,'amu*angstrom^2'), symmetry=1, barrier=(25.4216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134024,'amu*angstrom^2'), symmetry=1, barrier=(9.25903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27264,0.0530151,-2.86303e-05,6.96782e-09,-5.78367e-13,43441.7,27.1711], Tmin=(100,'K'), Tmax=(1635.89,'K')), NASAPolynomial(coeffs=[13.4028,0.0288102,-1.14381e-05,2.00001e-09,-1.30696e-13,38743,-39.5512], Tmin=(1635.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C(C)[CH]C(433)',
    structure = SMILES('[CH2]C(C)[CH]C'),
    E0 = (225.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1442.63],'cm^-1')),
        HinderedRotor(inertia=(0.759922,'amu*angstrom^2'), symmetry=1, barrier=(17.4721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00105696,'amu*angstrom^2'), symmetry=1, barrier=(7.72029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00522788,'amu*angstrom^2'), symmetry=1, barrier=(7.72135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.69593,'amu*angstrom^2'), symmetry=1, barrier=(84.9767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80045,0.0421185,-1.07299e-05,-7.7653e-09,4.08528e-12,27247.4,23.8048], Tmin=(100,'K'), Tmax=(1119.01,'K')), NASAPolynomial(coeffs=[7.75182,0.0320962,-1.23773e-05,2.20163e-09,-1.48913e-13,25211,-8.72391], Tmin=(1119.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Cs_S)"""),
)

species(
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.69489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49898e-06,-1.43376e-09,2.58636e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.04,'K')), NASAPolynomial(coeffs=[2.9759,0.00164141,-7.19722e-07,1.25378e-10,-7.91526e-15,-1025.84,5.53757], Tmin=(1817.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (341.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (465.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (481.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (493.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (386.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (458.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (450.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (465.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (444.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (523.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (456.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (459.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (483.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (493.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (385.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (424.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (416.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (645.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (430.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (419.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (790.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (789.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (501.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (536.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (435.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (486.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (349.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (741.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (627.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C3H6(72)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=[C]C(C)C(=C)C(24257)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]C(C)C(C)=C=C(24258)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.51e+07,'cm^3/(mol*s)'), n=1.64, Ea=(11.8407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2579 used for Cds-CsCs_Ca;HJ
Exact match found for rate rule [Cds-CsCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC(C)C([CH2])C(24259)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C3H6(72)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', 'C=[C]C(C)C=C(24201)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C3H6(T)(143)', 'CH3CHCCH2(18175)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(20.1294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH2]C(C)C=C=C(24260)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(10800,'cm^3/(mol*s)'), n=2.41, Ea=(32.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 597 used for Cds-CsH_Ca;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=[C]C(C)[C](C)C(24261)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['[CH2]C(C)[C](C)C=C(19125)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_Cs2]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(C)C([CH2])C(19130)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=[C][C](C)C(C)C(24262)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['[CH2][C](C)C(C)C=C(19127)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['[CH2]C(C)C([CH2])C=C(18971)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['[CH2]C([CH2])C(C)C=C(18362)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([C]=C)C(C)C(19129)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(228000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]C(C)C(C)C(24263)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(816000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C3H6(T)(143)', 'C=[C][CH]C(18176)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=CC(C)C(=C)C(19134)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=C=C(C)C(C)C(24264)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', '[CH2]CC(C)[C]=C(19002)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(S)(23)', '[CH2]C(C)C[C]=C(24254)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=[C]C(C)C[CH]C(19169)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=[C]C(C)[CH]CC(24265)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['[CH2]C(C)C(=C)[CH]C(24172)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['[CH2]C(C)[CH]C(=C)C(24266)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C)C(C)[C]=C(19128)'],
    products = ['C=C1CC(C)C1C(24267)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH2(19)', 'C=[C]C(C)[CH]C(24204)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C)[CH]C(433)', 'H2CC(41)'],
    products = ['[CH2]C(C)C(C)[C]=C(19128)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4248',
    isomers = [
        '[CH2]C(C)C(C)[C]=C(19128)',
    ],
    reactants = [
        ('C3H6(72)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4248',
    Tmin = (300,'K'),
    Tmax = (2000,'K'),
    Tcount = 8,
    Tlist = ([302.47,323.145,369.86,455.987,609.649,885.262,1353.64,1896.74],'K'),
    Pmin = (0.01,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.0125282,0.0667467,1,14.982,79.8202],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

