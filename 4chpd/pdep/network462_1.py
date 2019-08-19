species(
    label = '[CH2]OC([O])C(OO)C1CCCC1(2209)',
    structure = SMILES('[CH2]OC([O])C(OO)C1CCCC1'),
    E0 = (-174.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.651,0.113218,-7.92988e-05,2.76429e-08,-3.90559e-12,-20832.9,48.093], Tmin=(100,'K'), Tmax=(1632.5,'K')), NASAPolynomial(coeffs=[24.2519,0.049751,-2.09834e-05,3.82887e-09,-2.58767e-13,-29290.3,-89.5693], Tmin=(1632.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CCOJ) + radical(CsJOCH3)"""),
)

species(
    label = 'CH2O(35)',
    structure = SMILES('C=O'),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=CC(OO)C1CCCC1(1993)',
    structure = SMILES('O=CC(OO)C1CCCC1'),
    E0 = (-329.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4681.46,'J/mol'), sigma=(7.80337,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=731.23 K, Pc=22.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0321018,0.063013,3.67252e-05,-9.13792e-08,3.86252e-11,-39459.5,37.6681], Tmin=(100,'K'), Tmax=(989.421,'K')), NASAPolynomial(coeffs=[20.4293,0.0387808,-1.48147e-05,2.82869e-09,-2.07746e-13,-46345.9,-74.9226], Tmin=(989.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopentane)"""),
)

species(
    label = 'H(19)',
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
    label = '[CH2]OC(=O)C(OO)C1CCCC1(3871)',
    structure = SMILES('[CH2]OC(=O)C(OO)C1CCCC1'),
    E0 = (-375.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (173.186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23717,0.0912959,-6.58144e-06,-5.53241e-08,2.65775e-11,-44986.6,45.847], Tmin=(100,'K'), Tmax=(1023.19,'K')), NASAPolynomial(coeffs=[23.8484,0.0464986,-1.90031e-05,3.65227e-09,-2.64792e-13,-52908.6,-89.3782], Tmin=(1023.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-375.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(CsJOC(O)C)"""),
)

species(
    label = '[O][CH]C(OO)C1CCCC1(2032)',
    structure = SMILES('[O][CH]C(OO)C1CCCC1'),
    E0 = (-3.06575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.720628,0.0925877,-6.32804e-05,2.17495e-08,-3.03318e-12,-189.696,39.0194], Tmin=(100,'K'), Tmax=(1656.29,'K')), NASAPolynomial(coeffs=[20.6219,0.0410454,-1.66022e-05,2.96149e-09,-1.97353e-13,-7259.65,-74.715], Tmin=(1656.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.06575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'OO[CH]C1CCCC1(1997)',
    structure = SMILES('OO[CH]C1CCCC1'),
    E0 = (-24.2356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.743437,0.0552412,1.5949e-05,-5.51478e-08,2.35154e-11,-2783.25,28.8778], Tmin=(100,'K'), Tmax=(1007,'K')), NASAPolynomial(coeffs=[14.777,0.037646,-1.46669e-05,2.74074e-09,-1.956e-13,-7543.87,-48.5281], Tmin=(1007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.2356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + ring(Cyclopentane) + radical(CCsJOOH)"""),
)

species(
    label = 'CH2OCHO(3872)',
    structure = SMILES('[CH2]OC=O'),
    E0 = (-170.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1010.94],'cm^-1')),
        HinderedRotor(inertia=(1.46401,'amu*angstrom^2'), symmetry=1, barrier=(33.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395799,'amu*angstrom^2'), symmetry=1, barrier=(33.6572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79262,0.0222941,-9.08435e-06,-2.40312e-09,1.8225e-12,-20487.5,13.1527], Tmin=(100,'K'), Tmax=(1184.74,'K')), NASAPolynomial(coeffs=[8.25878,0.0116085,-5.39244e-06,1.05479e-09,-7.52244e-14,-22328,-16.4462], Tmin=(1184.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""CH2OCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[C](O)C(OO)C1CCCC1(3873)',
    structure = SMILES('[CH2]O[C](O)C(OO)C1CCCC1'),
    E0 = (-195.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42232,0.121809,-9.4856e-05,3.71099e-08,-5.78563e-12,-23258.2,50.7015], Tmin=(100,'K'), Tmax=(1528.89,'K')), NASAPolynomial(coeffs=[28.1201,0.0419018,-1.64595e-05,2.9256e-09,-1.95927e-13,-32597.4,-109.615], Tmin=(1528.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CsJOCH3) + radical(Cs_P)"""),
)

species(
    label = 'CO[C]([O])C(OO)C1CCCC1(3874)',
    structure = SMILES('CO[C]([O])C(OO)C1CCCC1'),
    E0 = (-157.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,360,370,350,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17889,0.108819,-7.28242e-05,2.4163e-08,-3.28392e-12,-18765.6,46.2785], Tmin=(100,'K'), Tmax=(1658.5,'K')), NASAPolynomial(coeffs=[21.4451,0.0542541,-2.34733e-05,4.32534e-09,-2.93595e-13,-26269.9,-74.3146], Tmin=(1658.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2]OC(O)[C](OO)C1CCCC1(3875)',
    structure = SMILES('[CH2]OC(O)[C](OO)C1CCCC1'),
    E0 = (-213.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75336,0.117233,-8.86073e-05,3.4184e-08,-5.37463e-12,-25497.9,47.8808], Tmin=(100,'K'), Tmax=(1481.18,'K')), NASAPolynomial(coeffs=[22.8247,0.0508588,-2.13905e-05,3.93048e-09,-2.6834e-13,-32778.9,-80.3497], Tmin=(1481.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(C2CsJOOH) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]OC(O)C(OO)[C]1CCCC1(3876)',
    structure = SMILES('[CH2]OC(O)C(OO)[C]1CCCC1'),
    E0 = (-209.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55932,0.104444,-5.16534e-05,-3.17503e-09,7.13986e-12,-25014.6,51.7364], Tmin=(100,'K'), Tmax=(1104.14,'K')), NASAPolynomial(coeffs=[21.9476,0.0505568,-2.09302e-05,3.92589e-09,-2.75902e-13,-32111.8,-72.6318], Tmin=(1104.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(C2CJCOOH) + radical(CsJOCH3)"""),
)

species(
    label = 'COC([O])[C](OO)C1CCCC1(3877)',
    structure = SMILES('COC([O])[C](OO)C1CCCC1'),
    E0 = (-175.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.47832,0.103959,-6.58994e-05,2.06689e-08,-2.71879e-12,-21007.1,43.337], Tmin=(100,'K'), Tmax=(1630.84,'K')), NASAPolynomial(coeffs=[16.1787,0.063104,-2.83214e-05,5.30729e-09,-3.63911e-13,-26440,-45.1701], Tmin=(1630.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CCOJ) + radical(C2CsJOOH)"""),
)

species(
    label = '[CH2]OC(O)C(O[O])C1CCCC1(3878)',
    structure = SMILES('[CH2]OC(O)C(O[O])C1CCCC1'),
    E0 = (-248.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18625,0.117404,-8.90012e-05,3.42994e-08,-5.2961e-12,-29670.8,49.9234], Tmin=(100,'K'), Tmax=(1537.7,'K')), NASAPolynomial(coeffs=[26.3411,0.0431963,-1.66123e-05,2.91527e-09,-1.93635e-13,-38444.1,-99.9795], Tmin=(1537.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(ROOJ) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]OC(O)C(OO)C1[CH]CCC1(3879)',
    structure = SMILES('[CH2]OC(O)C(OO)C1[CH]CCC1'),
    E0 = (-206.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.06484,0.11751,-8.77897e-05,3.28329e-08,-4.92832e-12,-24561.6,51.2823], Tmin=(100,'K'), Tmax=(1568.59,'K')), NASAPolynomial(coeffs=[26.422,0.0448687,-1.83264e-05,3.31105e-09,-2.23271e-13,-33498.7,-98.9749], Tmin=(1568.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-206.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(CsJOCH3)"""),
)

species(
    label = 'COC([O])C(OO)[C]1CCCC1(3880)',
    structure = SMILES('COC([O])C(OO)[C]1CCCC1'),
    E0 = (-171.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20577,0.100817,-5.76224e-05,1.43163e-08,-1.1484e-12,-20480.5,50.5891], Tmin=(100,'K'), Tmax=(1587.84,'K')), NASAPolynomial(coeffs=[23.7739,0.0500438,-2.114e-05,3.81974e-09,-2.54767e-13,-29945.5,-86.2982], Tmin=(1587.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CCOJ) + radical(C2CJCOOH)"""),
)

species(
    label = '[CH2]OC(O)C(OO)C1C[CH]CC1(3881)',
    structure = SMILES('[CH2]OC(O)C(OO)C1C[CH]CC1'),
    E0 = (-214.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56464,0.114379,-8.48414e-05,3.23383e-08,-5.05277e-12,-25628.3,49.4707], Tmin=(100,'K'), Tmax=(1481.67,'K')), NASAPolynomial(coeffs=[21.4684,0.0521977,-2.18906e-05,4.01392e-09,-2.73632e-13,-32453.7,-70.7063], Tmin=(1481.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(602.799,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CsJOCH3) + radical(cyclopentane)"""),
)

species(
    label = 'COC([O])C(O[O])C1CCCC1(3882)',
    structure = SMILES('COC([O])C(O[O])C1CCCC1'),
    E0 = (-210.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,492.5,1135,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.941538,0.104411,-6.69848e-05,2.13872e-08,-2.80919e-12,-25178.4,45.4947], Tmin=(100,'K'), Tmax=(1697.84,'K')), NASAPolynomial(coeffs=[20.0198,0.0550275,-2.33557e-05,4.25598e-09,-2.86693e-13,-32296.2,-66.7276], Tmin=(1697.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(611.114,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = 'COC([O])C(OO)C1[CH]CCC1(3883)',
    structure = SMILES('COC([O])C(OO)C1[CH]CCC1'),
    E0 = (-168.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3615,1310,387.5,850,1000,2750,2807.14,2864.29,2921.43,2978.57,3035.71,3092.86,3150,900,928.571,957.143,985.714,1014.29,1042.86,1071.43,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.813931,0.104485,-6.57735e-05,1.99856e-08,-2.47463e-12,-20069.8,46.8273], Tmin=(100,'K'), Tmax=(1786.52,'K')), NASAPolynomial(coeffs=[21.3314,0.0549015,-2.41423e-05,4.45025e-09,-3.00654e-13,-27982.4,-72.8612], Tmin=(1786.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(Cs_S) + radical(CCOJ)"""),
)

species(
    label = 'COC([O])C(OO)C1C[CH]CC1(3884)',
    structure = SMILES('COC([O])C(OO)C1C[CH]CC1'),
    E0 = (-176.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.276949,0.101003,-6.19235e-05,1.86777e-08,-2.36549e-12,-21138.2,44.8774], Tmin=(100,'K'), Tmax=(1668.54,'K')), NASAPolynomial(coeffs=[15.2433,0.0637958,-2.84744e-05,5.31301e-09,-3.63012e-13,-26317.4,-37.9443], Tmin=(1668.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-176.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CCOJ) + radical(cyclopentane)"""),
)

species(
    label = '[CH2][O](354)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = 'COC(=O)C(OO)C1CCCC1(3885)',
    structure = SMILES('COC(=O)C(OO)C1CCCC1'),
    E0 = (-577.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.931038,0.084927,9.10588e-06,-6.61361e-08,2.90375e-11,-69257.5,45.0749], Tmin=(100,'K'), Tmax=(1028.09,'K')), NASAPolynomial(coeffs=[21.0563,0.0526411,-2.14962e-05,4.09787e-09,-2.94567e-13,-76593.2,-75.2994], Tmin=(1028.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-577.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(606.956,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-OsHHH) + group(Cds-OdCsOs) + ring(Cyclopentane)"""),
)

species(
    label = 'OH(21)',
    structure = SMILES('[OH]'),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]OC1OOC1C1CCCC1(3886)',
    structure = SMILES('[CH2]OC1OOC1C1CCCC1'),
    E0 = (-71.5006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.506572,0.0791709,4.55458e-07,-5.05291e-08,2.2765e-11,-8420.05,36.2107], Tmin=(100,'K'), Tmax=(1036.37,'K')), NASAPolynomial(coeffs=[18.7671,0.0490196,-1.99325e-05,3.77269e-09,-2.69317e-13,-14790.7,-68.9238], Tmin=(1036.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.5006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(569.541,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsHHH) + ring(Cyclopentane) + ring(12dioxetane) + radical(CsJOCH3)"""),
)

species(
    label = '[O]C1OCOC1C1CCCC1(3887)',
    structure = SMILES('[O]C1OCOC1C1CCCC1'),
    E0 = (-364.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (157.187,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0155438,0.0470342,0.000130138,-2.29967e-07,1.01538e-10,-43679.9,33.649], Tmin=(100,'K'), Tmax=(898.146,'K')), NASAPolynomial(coeffs=[32.5009,0.0119097,5.60304e-06,-1.54689e-09,1.06267e-13,-53945,-144.361], Tmin=(898.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-364.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(577.856,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + ring(1,3-Dioxolane) + ring(Cyclopentane) + radical(CCOJ)"""),
)

species(
    label = 'OOC(C1CCCC1)C1OCO1(2212)',
    structure = SMILES('OOC(C1CCCC1)C1OCO1'),
    E0 = (-431.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (174.194,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24276,0.0738776,8.52509e-05,-1.76485e-07,7.57455e-11,-51700.9,42.4083], Tmin=(100,'K'), Tmax=(957.744,'K')), NASAPolynomial(coeffs=[33.6741,0.0311787,-9.39595e-06,1.8286e-09,-1.47683e-13,-63119.1,-149.231], Tmin=(957.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-431.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(615.271,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-OsOsHH) + ring(Cyclobutane) + ring(Cyclopentane)"""),
)

species(
    label = 'CH2(T)(33)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C([O])C(OO)C1CCCC1(3888)',
    structure = SMILES('[O]C([O])C(OO)C1CCCC1'),
    E0 = (-157.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (160.168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.378104,0.0915676,-5.59849e-05,1.62732e-08,-1.9005e-12,-18829.2,41.5586], Tmin=(100,'K'), Tmax=(1916.7,'K')), NASAPolynomial(coeffs=[22.1577,0.0445372,-1.91792e-05,3.47147e-09,-2.3074e-13,-27468.1,-81.8258], Tmin=(1916.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + ring(Cyclopentane) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'O(20)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]O[CH]C(OO)C1CCCC1(2039)',
    structure = SMILES('[CH2]O[CH]C(OO)C1CCCC1'),
    E0 = (-20.0954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,2750,2800,2850,2900,2950,3000,3050,3100,3150,900,925,950,975,1000,1025,1050,1075,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.195,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57447,0.10489,-5.31414e-05,-7.71822e-09,1.08919e-11,-2200.76,44.1225], Tmin=(100,'K'), Tmax=(1019.54,'K')), NASAPolynomial(coeffs=[22.9803,0.0444107,-1.69147e-05,3.0888e-09,-2.1655e-13,-9071.26,-83.9544], Tmin=(1019.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.0954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(582.013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-OsHHH) + ring(Cyclopentane) + radical(CCsJOCs) + radical(CsJOCC)"""),
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
    E0 = (-174.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-117.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-94.0875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-163.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-16.9415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (19.2556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-91.3187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-74.1438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-66.8213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-116.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-149.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-110.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-95.2726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-117.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-107.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-109.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (189.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-111.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-39.5957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-113.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-166.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (223.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (222.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['CH2O(35)', 'O=CC(OO)C1CCCC1(1993)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(19)', '[CH2]OC(=O)C(OO)C1CCCC1(3871)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2O(35)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2330,'cm^3/(mol*s)'), n=3.17, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OO[CH]C1CCCC1(1997)', 'CH2OCHO(3872)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.09579,'m^3/(mol*s)'), n=1.29734, Ea=(31.9242,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ] for rate rule [CO-NdH_O;CsJ-OsCsH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['[CH2]O[C](O)C(OO)C1CCCC1(3873)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CO[C]([O])C(OO)C1CCCC1(3874)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.71e+08,'s^-1'), n=1.45, Ea=(176.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['[CH2]OC(O)[C](OO)C1CCCC1(3875)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0250998,'s^-1'), n=4.2935, Ea=(83.6591,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['[CH2]OC(O)C(OO)[C]1CCCC1(3876)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0, Ea=(100.834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;Cs_H_out_Cs2] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_Cs2_cy5]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['COC([O])[C](OO)C1CCCC1(3877)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.19773e+08,'s^-1'), n=1.45, Ea=(108.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS_OCs;Y_rad_out;Cs_H_out_OOH/Cs] for rate rule [R4H_SSS_OCs;C_rad_out_2H;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['[CH2]OC(O)C(O[O])C1CCCC1(3878)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['[CH2]OC(O)C(OO)C1[CH]CCC1(3879)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.6e+11,'s^-1'), n=0, Ea=(25.7316,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_CCC;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['COC([O])C(OO)[C]1CCCC1(3880)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.174e+08,'s^-1'), n=1.15, Ea=(64.3081,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_SSSS_OCC;Y_rad_out;Cs_H_out_Cs2] for rate rule [R5H_SSSS_OCC;C_rad_out_2H;Cs_H_out_Cs2_cy5]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['[CH2]OC(O)C(OO)C1C[CH]CC1(3881)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.24e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['COC([O])C(O[O])C1CCCC1(3882)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 140 used for R6H_SSSSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['COC([O])C(OO)C1[CH]CCC1(3883)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(138.3,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6H_SSSSS;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['COC([O])C(OO)C1C[CH]CC1(3884)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(107887,'s^-1'), n=1.4, Ea=(65.6261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7H;C_rad_out_2H;Cs_H_out_H/NonDeC] + [R7H_OOCCCC(Cs/Cs);Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R7H_OOCCCC(Cs/Cs);C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][O](354)', '[O][CH]C(OO)C1CCCC1(2032)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.35017e+06,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['COC(=O)C(OO)C1CCCC1(3885)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['OH(21)', '[CH2]OC1OOC1C1CCCC1(3886)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(135.382,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO_SS;Y_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['OH(21)', '[O]C1OCOC1C1CCCC1(3887)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.13e+10,'s^-1','*|/',1.41), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4OO;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    products = ['OOC(C1CCCC1)C1OCO1(2212)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(T)(33)', '[O]C([O])C(OO)C1CCCC1(3888)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(109477,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O(20)', '[CH2]O[CH]C(OO)C1CCCC1(2039)'],
    products = ['[CH2]OC([O])C(OO)C1CCCC1(2209)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '462',
    isomers = [
        '[CH2]OC([O])C(OO)C1CCCC1(2209)',
    ],
    reactants = [
        ('CH2O(35)', 'O=CC(OO)C1CCCC1(1993)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '462',
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

