species(
    label = 'CCCC(C)O[O](591)',
    structure = SMILES('CCCC(C)O[O]'),
    E0 = (-117.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,492.5,1135,1000,1380,1390,370,380,2900,435,287.073,287.346],'cm^-1')),
        HinderedRotor(inertia=(0.00203368,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154178,'amu*angstrom^2'), symmetry=1, barrier=(9.0585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155102,'amu*angstrom^2'), symmetry=1, barrier=(9.06012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155066,'amu*angstrom^2'), symmetry=1, barrier=(9.05818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155803,'amu*angstrom^2'), symmetry=1, barrier=(9.06601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.819521,0.0706875,-5.22391e-05,2.12644e-08,-3.71969e-12,-13963.6,27.1405], Tmin=(100,'K'), Tmax=(1290.5,'K')), NASAPolynomial(coeffs=[10.5251,0.0406044,-1.72721e-05,3.20061e-09,-2.20297e-13,-16468.6,-22.1584], Tmin=(1290.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ)"""),
)

species(
    label = 'C=CCCC(134)',
    structure = SMILES('C=CCCC'),
    E0 = (-40.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,434.463,434.474],'cm^-1')),
        HinderedRotor(inertia=(0.0736945,'amu*angstrom^2'), symmetry=1, barrier=(9.87096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.073698,'amu*angstrom^2'), symmetry=1, barrier=(9.87114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0736934,'amu*angstrom^2'), symmetry=1, barrier=(9.87114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3052.11,'J/mol'), sigma=(5.53315,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.73 K, Pc=40.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9015,0.0362536,1.29132e-05,-3.55767e-08,1.43068e-11,-4763.1,19.8563], Tmin=(100,'K'), Tmax=(1027.61,'K')), NASAPolynomial(coeffs=[9.28067,0.0304042,-1.19376e-05,2.20664e-09,-1.55067e-13,-7487.42,-21.8213], Tmin=(1027.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'HO2(10)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CCCCO[O](164)',
    structure = SMILES('CCCCO[O]'),
    E0 = (-80.6218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,232.129,232.181],'cm^-1')),
        HinderedRotor(inertia=(0.157672,'amu*angstrom^2'), symmetry=1, barrier=(6.03016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157561,'amu*angstrom^2'), symmetry=1, barrier=(6.02839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157484,'amu*angstrom^2'), symmetry=1, barrier=(6.02954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157681,'amu*angstrom^2'), symmetry=1, barrier=(6.02962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3638.29,'J/mol'), sigma=(6.34008,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.29 K, Pc=32.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56797,0.0567642,-4.41695e-05,2.06267e-08,-4.35472e-12,-9611.76,22.5652], Tmin=(100,'K'), Tmax=(1060.64,'K')), NASAPolynomial(coeffs=[6.8104,0.0369938,-1.62098e-05,3.05291e-09,-2.12531e-13,-10723.8,-3.03526], Tmin=(1060.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.6218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ)"""),
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
    label = 'CCC[C](C)OO(5835)',
    structure = SMILES('CCC[C](C)OO'),
    E0 = (-82.1506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895111,0.0741148,-6.19761e-05,3.15338e-08,-7.31058e-12,-9773.52,26.4261], Tmin=(100,'K'), Tmax=(971.947,'K')), NASAPolynomial(coeffs=[7.14942,0.0483755,-2.22529e-05,4.28742e-09,-3.024e-13,-10989.3,-3.56935], Tmin=(971.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.1506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOOH)"""),
)

species(
    label = 'CC[CH]C(C)OO(5836)',
    structure = SMILES('CC[CH]C(C)OO'),
    E0 = (-68.6347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.673471,0.075176,-5.93309e-05,2.55748e-08,-4.71043e-12,-8136.97,29.2135], Tmin=(100,'K'), Tmax=(1237.03,'K')), NASAPolynomial(coeffs=[11.139,0.0413353,-1.82964e-05,3.46035e-09,-2.41178e-13,-10726.2,-23.5029], Tmin=(1237.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.6347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(CCC)OO(331)',
    structure = SMILES('[CH2]C(CCC)OO'),
    E0 = (-55.0871,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424678,0.0777384,-6.31189e-05,2.74429e-08,-4.97666e-12,-6496.19,29.796], Tmin=(100,'K'), Tmax=(1280.05,'K')), NASAPolynomial(coeffs=[13.2212,0.0377512,-1.6261e-05,3.03877e-09,-2.10445e-13,-9772.24,-35.0991], Tmin=(1280.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.0871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCOOH)"""),
)

species(
    label = 'C[CH]CC(C)OO(1959)',
    structure = SMILES('C[CH]CC(C)OO'),
    E0 = (-74.6034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944516,0.0720429,-5.76665e-05,2.76736e-08,-6.01122e-12,-8866.76,28.6842], Tmin=(100,'K'), Tmax=(1032.91,'K')), NASAPolynomial(coeffs=[7.46856,0.0467781,-2.09767e-05,3.99298e-09,-2.79667e-13,-10214.5,-3.00175], Tmin=(1032.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.6034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCC(C)OO(1958)',
    structure = SMILES('[CH2]CCC(C)OO'),
    E0 = (-63.8033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.14,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.591036,0.0750202,-5.79103e-05,2.39063e-08,-4.15849e-12,-7551.39,28.9886], Tmin=(100,'K'), Tmax=(1314.99,'K')), NASAPolynomial(coeffs=[12.364,0.0392092,-1.70615e-05,3.19732e-09,-2.21461e-13,-10647.7,-31.0332], Tmin=(1314.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.8033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C[CH]CCC(139)',
    structure = SMILES('C[CH]CCC'),
    E0 = (26.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1584.16,1584.16],'cm^-1')),
        HinderedRotor(inertia=(0.289393,'amu*angstrom^2'), symmetry=1, barrier=(8.10552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0042711,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289397,'amu*angstrom^2'), symmetry=1, barrier=(8.10552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289396,'amu*angstrom^2'), symmetry=1, barrier=(8.10551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.1408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77405,0.0441717,-1.45099e-05,-1.38039e-09,1.20138e-12,3254.59,22.595], Tmin=(100,'K'), Tmax=(1461.29,'K')), NASAPolynomial(coeffs=[8.55193,0.0345019,-1.37025e-05,2.41137e-09,-1.59036e-13,325.233,-15.9206], Tmin=(1461.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]O[O](97)',
    structure = SMILES('C[CH]O[O]'),
    E0 = (162.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.659845,'amu*angstrom^2'), symmetry=1, barrier=(15.1711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000844769,'amu*angstrom^2'), symmetry=1, barrier=(5.4037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88063,0.0212973,-8.38042e-06,-1.15683e-09,1.18565e-12,19629.9,16.0899], Tmin=(100,'K'), Tmax=(1195.85,'K')), NASAPolynomial(coeffs=[6.93453,0.0132692,-5.24916e-06,9.65724e-10,-6.67381e-14,18264.7,-5.84667], Tmin=(1195.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'npropyl(83)',
    structure = SMILES('[CH2]CC'),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928812,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092914,'amu*angstrom^2'), symmetry=1, barrier=(2.13628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02815,0.0147023,2.4051e-05,-3.66738e-08,1.38611e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.464,'K')), NASAPolynomial(coeffs=[6.16543,0.0184495,-6.79029e-06,1.23049e-09,-8.63866e-14,9095.06,-6.67607], Tmin=(984.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CCC[CH]O[O](257)',
    structure = SMILES('CCC[CH]O[O]'),
    E0 = (107.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.119013,'amu*angstrom^2'), symmetry=1, barrier=(2.73635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120007,'amu*angstrom^2'), symmetry=1, barrier=(2.75919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121685,'amu*angstrom^2'), symmetry=1, barrier=(2.79779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481289,'amu*angstrom^2'), symmetry=1, barrier=(11.0658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80219,0.0553016,-3.67985e-05,-1.19632e-08,2.59165e-11,13057,21.9547], Tmin=(100,'K'), Tmax=(531.462,'K')), NASAPolynomial(coeffs=[5.78891,0.0367206,-1.66008e-05,3.14923e-09,-2.19362e-13,12471.9,3.72319], Tmin=(531.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH)"""),
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
    label = 'CCC[C](C)O[O](5693)',
    structure = SMILES('CCC[C](C)O[O]'),
    E0 = (69.8541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,492.5,1135,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,247.255,1440.44],'cm^-1')),
        HinderedRotor(inertia=(0.153406,'amu*angstrom^2'), symmetry=1, barrier=(7.05468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155857,'amu*angstrom^2'), symmetry=1, barrier=(7.06985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158265,'amu*angstrom^2'), symmetry=1, barrier=(7.04549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15756,'amu*angstrom^2'), symmetry=1, barrier=(7.03921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159956,'amu*angstrom^2'), symmetry=1, barrier=(7.03974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46528,0.064758,-3.2121e-05,-3.76886e-08,4.8021e-11,8483.83,25.2528], Tmin=(100,'K'), Tmax=(505.746,'K')), NASAPolynomial(coeffs=[5.54392,0.046813,-2.13501e-05,4.07342e-09,-2.84875e-13,7888.22,6.54639], Tmin=(505.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.8541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(C)O[O](218)',
    structure = SMILES('[CH2]C(C)O[O]'),
    E0 = (145.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.220885,'amu*angstrom^2'), symmetry=1, barrier=(5.07858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221154,'amu*angstrom^2'), symmetry=1, barrier=(5.08476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221416,'amu*angstrom^2'), symmetry=1, barrier=(5.09078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.31,'J/mol'), sigma=(6.05676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.08 K, Pc=35.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89966,0.048151,-5.47371e-05,3.64967e-08,-1.00252e-11,17622.9,21.4278], Tmin=(100,'K'), Tmax=(880.217,'K')), NASAPolynomial(coeffs=[7.85206,0.0211006,-8.63848e-06,1.58105e-09,-1.08162e-13,16575,-6.52949], Tmin=(880.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'C2H5(29)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CC[CH]C(C)O[O](5694)',
    structure = SMILES('CC[CH]C(C)O[O]'),
    E0 = (83.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,492.5,1135,1000,1380,1390,370,380,2900,435,325.846,325.868],'cm^-1')),
        HinderedRotor(inertia=(0.00158639,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100096,'amu*angstrom^2'), symmetry=1, barrier=(7.53428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.099607,'amu*angstrom^2'), symmetry=1, barrier=(7.5206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0999526,'amu*angstrom^2'), symmetry=1, barrier=(7.52373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09997,'amu*angstrom^2'), symmetry=1, barrier=(7.52812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925438,0.0715654,-6.33827e-05,3.35041e-08,-7.76446e-12,10134.5,29.0807], Tmin=(100,'K'), Tmax=(997.362,'K')), NASAPolynomial(coeffs=[8.36553,0.041726,-1.85047e-05,3.50595e-09,-2.44995e-13,8650.42,-6.79377], Tmin=(997.362,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CC(C)O[O](220)',
    structure = SMILES('[CH2]CC(C)O[O]'),
    E0 = (111.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,373.188],'cm^-1')),
        HinderedRotor(inertia=(0.0624435,'amu*angstrom^2'), symmetry=1, barrier=(6.1951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985009,'amu*angstrom^2'), symmetry=1, barrier=(9.37411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0596551,'amu*angstrom^2'), symmetry=1, barrier=(6.17246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620122,'amu*angstrom^2'), symmetry=1, barrier=(6.24353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53632,0.0561472,-4.67691e-05,2.24435e-08,-4.6632e-12,13555.3,24.1169], Tmin=(100,'K'), Tmax=(1107.31,'K')), NASAPolynomial(coeffs=[8.27729,0.0317965,-1.37829e-05,2.58384e-09,-1.79462e-13,12062.5,-9.09151], Tmin=(1107.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'C[CH]CC(C)O[O](1941)',
    structure = SMILES('C[CH]CC(C)O[O]'),
    E0 = (77.4013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,492.5,1135,1000,1380,1390,370,380,2900,435,212.94,1691.28],'cm^-1')),
        HinderedRotor(inertia=(0.00371782,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207373,'amu*angstrom^2'), symmetry=1, barrier=(6.67246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20737,'amu*angstrom^2'), symmetry=1, barrier=(6.67241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207352,'amu*angstrom^2'), symmetry=1, barrier=(6.67241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20737,'amu*angstrom^2'), symmetry=1, barrier=(6.6724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40703,0.0646267,-3.95767e-05,-1.32593e-08,2.63781e-11,9395.3,27.864], Tmin=(100,'K'), Tmax=(542.574,'K')), NASAPolynomial(coeffs=[5.60228,0.0456864,-2.03571e-05,3.84796e-09,-2.67996e-13,8763.59,8.56327], Tmin=(542.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.4013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(CCC)O[O](259)',
    structure = SMILES('[CH2]C(CCC)O[O]'),
    E0 = (96.9176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,341.753,342.286],'cm^-1')),
        HinderedRotor(inertia=(0.100562,'amu*angstrom^2'), symmetry=1, barrier=(8.34038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100627,'amu*angstrom^2'), symmetry=1, barrier=(8.34316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100614,'amu*angstrom^2'), symmetry=1, barrier=(8.34225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100471,'amu*angstrom^2'), symmetry=1, barrier=(8.34355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3800.62,'J/mol'), sigma=(6.68442,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.65 K, Pc=28.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.761822,0.0732136,-6.43675e-05,3.22252e-08,-6.86808e-12,11771.4,29.3513], Tmin=(100,'K'), Tmax=(1094.47,'K')), NASAPolynomial(coeffs=[10.3082,0.0383245,-1.65516e-05,3.09981e-09,-2.15291e-13,9681.76,-17.5662], Tmin=(1094.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.9176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]CCC(C)O[O](221)',
    structure = SMILES('[CH2]CCC(C)O[O]'),
    E0 = (88.2014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,289.047,289.076],'cm^-1')),
        HinderedRotor(inertia=(0.00201822,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00201781,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136593,'amu*angstrom^2'), symmetry=1, barrier=(8.09619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136588,'amu*angstrom^2'), symmetry=1, barrier=(8.09605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136575,'amu*angstrom^2'), symmetry=1, barrier=(8.09618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92117,0.0705479,-5.92358e-05,2.86917e-08,-6.02768e-12,10716.6,28.5715], Tmin=(100,'K'), Tmax=(1094.77,'K')), NASAPolynomial(coeffs=[9.2491,0.0401197,-1.75444e-05,3.30343e-09,-2.30025e-13,8893.22,-12.3599], Tmin=(1094.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.2014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'CCC(C)O[O](716)',
    structure = SMILES('CCC(C)O[O]'),
    E0 = (-93.2647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,492.5,1135,1000,1380,1390,370,380,2900,435,314.511],'cm^-1')),
        HinderedRotor(inertia=(0.126123,'amu*angstrom^2'), symmetry=1, barrier=(8.85324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126139,'amu*angstrom^2'), symmetry=1, barrier=(8.85306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126136,'amu*angstrom^2'), symmetry=1, barrier=(8.85339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126152,'amu*angstrom^2'), symmetry=1, barrier=(8.85338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41175,0.0565301,-4.05021e-05,1.58077e-08,-2.63472e-12,-11123.9,22.7699], Tmin=(100,'K'), Tmax=(1352.72,'K')), NASAPolynomial(coeffs=[9.67023,0.0321097,-1.34229e-05,2.4622e-09,-1.68292e-13,-13358.2,-19.5675], Tmin=(1352.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.2647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ)"""),
)

species(
    label = 'CC=CCC(460)',
    structure = SMILES('CC=CCC'),
    E0 = (-52.5474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94607,0.0364834,8.35486e-06,-2.83273e-08,1.11812e-11,-6238.66,19.2593], Tmin=(100,'K'), Tmax=(1059.79,'K')), NASAPolynomial(coeffs=[8.36587,0.0318147,-1.27245e-05,2.34965e-09,-1.63922e-13,-8697.93,-17.2682], Tmin=(1059.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.5474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'O(4)',
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
    label = 'CCCC(C)[O](5837)',
    structure = SMILES('CCCC(C)[O]'),
    E0 = (-110.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,1036.77],'cm^-1')),
        HinderedRotor(inertia=(0.198631,'amu*angstrom^2'), symmetry=1, barrier=(4.56693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800707,'amu*angstrom^2'), symmetry=1, barrier=(18.4098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800555,'amu*angstrom^2'), symmetry=1, barrier=(18.4063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00162151,'amu*angstrom^2'), symmetry=1, barrier=(18.4107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.1402,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.935002,0.0586257,-2.83009e-05,9.62055e-10,2.31732e-12,-13135.6,25.0232], Tmin=(100,'K'), Tmax=(1158.21,'K')), NASAPolynomial(coeffs=[12.0463,0.0337058,-1.34516e-05,2.44427e-09,-1.67469e-13,-16611.9,-34.1103], Tmin=(1158.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ)"""),
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
    E0 = (8.89301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (43.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (15.4637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (23.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-39.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-7.29841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (17.7402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (249.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (244.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (281.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (253.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (295.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (248.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (289.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (308.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (299.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (339.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (326.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (6.80101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (132.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CCCC(C)O[O](591)'],
    products = ['C=CCCC(134)', 'HO2(10)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.937e+09,'s^-1'), n=1.17, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for R2OO_2H_HNd
Exact match found for rate rule [R2OO_2H_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCCC(C)O[O](591)'],
    products = ['CCC[C](C)OO(5835)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.97e+09,'s^-1'), n=1.01, Ea=(160.958,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 245 used for R3H_SS_O;O_rad_out;Cs_H_out_Cs2
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC[CH]C(C)OO(5836)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 323 used for R4H_SSS;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CCCC(C)O[O](591)'],
    products = ['[CH2]C(CCC)OO(331)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.18e+08,'s^-1'), n=1.06, Ea=(140.206,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 247 used for R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CCCC(C)O[O](591)'],
    products = ['C[CH]CC(C)OO(1959)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.76e+10,'s^-1'), n=0.21, Ea=(77.404,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 261 C5H11O2-11 <=> C5H11O2-12 in intra_H_migration/training
This reaction matched rate rule [R5H_SSSS_OCC_C;O_rad_out;Cs_H_out_H/NonDeC]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CCC(C)OO(1958)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(135.127,'s^-1'), n=2.18479, Ea=(56.5049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;C_rad_out_2H;XH_out] for rate rule [R6H_SSSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O2(2)', 'C[CH]CCC(139)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.54e+12,'cm^3/(mol*s)','+|-',1e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 112 used for O2_birad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]O[O](97)', 'npropyl(83)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.15e+14,'cm^3/(mol*s)','*|/',2), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_rad/H/NonDe;C_rad/H2/Cs] for rate rule [C_rad/H/CsO;C_rad/H2/Cs]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CCC[CH]O[O](257)', 'CH3(17)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.26832e+07,'m^3/(mol*s)'), n=-0.3275, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;C_rad/H/CsO] + [C_methyl;C_rad/H/NonDe] for rate rule [C_methyl;C_rad/H/CsO]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'CCC[C](C)O[O](5693)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.68468e+06,'m^3/(mol*s)'), n=-0.0135, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/NonDe;Y_rad] for rate rule [C_rad/NonDeCO;H_rad]
Euclidian distance = 2.2360679775
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)O[O](218)', 'C2H5(29)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.73e+14,'cm^3/(mol*s)'), n=-0.699, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 11 used for C_rad/H2/Cs;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', 'CC[CH]C(C)O[O](5694)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC(C)O[O](220)', 'CH3(17)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'C[CH]CC(C)O[O](1941)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2]C(CCC)O[O](259)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]CCC(C)O[O](221)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CCCCO[O](164)', 'CH2(S)(23)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/CsO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(S)(23)', 'CCC(C)O[O](716)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CCCC(C)O[O](591)'],
    products = ['CC=CCC(460)', 'HO2(10)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.3e+09,'s^-1'), n=1.01, Ea=(123.846,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 6 used for R2OO_HNd_HNd
Exact match found for rate rule [R2OO_HNd_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O(4)', 'CCCC(C)[O](5837)'],
    products = ['CCCC(C)O[O](591)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '243',
    isomers = [
        'CCCC(C)O[O](591)',
    ],
    reactants = [
        ('C=CCCC(134)', 'HO2(10)'),
        ('CCCCO[O](164)', 'CH2(S)(23)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '243',
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

