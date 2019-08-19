species(
    label = 'S(226)(225)',
    structure = SMILES('C=CC(C=CO)O[O]'),
    E0 = (-74.3664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.843722,'amu*angstrom^2'), symmetry=1, barrier=(19.3988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843369,'amu*angstrom^2'), symmetry=1, barrier=(19.3907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842685,'amu*angstrom^2'), symmetry=1, barrier=(19.375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.842892,'amu*angstrom^2'), symmetry=1, barrier=(19.3797,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4426.56,'J/mol'), sigma=(7.01777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.42 K, Pc=29.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0486304,0.0746373,-5.27035e-05,2.94134e-09,7.56643e-12,-8790.83,31.4628], Tmin=(100,'K'), Tmax=(957.415,'K')), NASAPolynomial(coeffs=[21.1448,0.0152251,-4.62649e-06,8.02403e-10,-5.80316e-14,-14146.9,-76.2715], Tmin=(957.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.3664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'S(233)(232)',
    structure = SMILES('C=CC(C=C[O])OO'),
    E0 = (-84.9085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06174,'amu*angstrom^2'), symmetry=1, barrier=(24.4114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06067,'amu*angstrom^2'), symmetry=1, barrier=(24.3869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06152,'amu*angstrom^2'), symmetry=1, barrier=(24.4065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06091,'amu*angstrom^2'), symmetry=1, barrier=(24.3925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4426.56,'J/mol'), sigma=(7.01777,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.42 K, Pc=29.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.155241,0.0717014,-4.59171e-05,1.17414e-09,6.2146e-12,-10062.4,31.7334], Tmin=(100,'K'), Tmax=(1015.29,'K')), NASAPolynomial(coeffs=[19.988,0.0197165,-7.7504e-06,1.48253e-09,-1.0822e-13,-15437.5,-70.8872], Tmin=(1015.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.9085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'O2(2)(2)',
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
    label = 'C5H7O(224)(223)',
    structure = SMILES('[CH2]C=CC=CO'),
    E0 = (-32.8413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.31979,'amu*angstrom^2'), symmetry=1, barrier=(30.3447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32385,'amu*angstrom^2'), symmetry=1, barrier=(30.4378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31983,'amu*angstrom^2'), symmetry=1, barrier=(30.3455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25554,0.0396093,3.48643e-05,-9.07575e-08,4.31029e-11,-3831.66,21.0706], Tmin=(100,'K'), Tmax=(908.805,'K')), NASAPolynomial(coeffs=[21.7542,0.00426328,2.6289e-06,-6.6847e-10,4.32897e-14,-9823.72,-88.3316], Tmin=(908.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.8413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'HO2(8)(9)',
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
    label = 'C5H6O(219)(218)',
    structure = SMILES('[CH2]C=CC=C[O]'),
    E0 = (108.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.59878,'amu*angstrom^2'), symmetry=1, barrier=(36.759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59525,'amu*angstrom^2'), symmetry=1, barrier=(36.6779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.673,0.0325848,3.83615e-05,-8.44753e-08,3.82236e-11,13165.2,20.9814], Tmin=(100,'K'), Tmax=(922.022,'K')), NASAPolynomial(coeffs=[18.3994,0.00812125,-9.23824e-08,-9.0725e-11,1.79375e-15,8036.17,-69.4433], Tmin=(922.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C5H6O(217)(216)',
    structure = SMILES('C=CC=CC=O'),
    E0 = (-29.5668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.9508,'amu*angstrom^2'), symmetry=1, barrier=(21.8608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.95303,'amu*angstrom^2'), symmetry=1, barrier=(21.912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3519.98,'J/mol'), sigma=(5.63814,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.81 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58677,0.0449815,-2.33922e-05,-1.17435e-10,2.4605e-12,-3462.46,19.7432], Tmin=(100,'K'), Tmax=(1159.43,'K')), NASAPolynomial(coeffs=[12.6828,0.020301,-9.05764e-06,1.75766e-09,-1.25367e-13,-6949.62,-39.3724], Tmin=(1159.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.5668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=CC1OOC1[CH]O(7310)',
    structure = SMILES('C=CC1OOC1[CH]O'),
    E0 = (35.5786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688995,0.0722601,-6.63359e-05,3.21359e-08,-6.37009e-12,4398.76,25.9518], Tmin=(100,'K'), Tmax=(1196.02,'K')), NASAPolynomial(coeffs=[13.3659,0.0298633,-1.31637e-05,2.49749e-09,-1.74897e-13,1366.39,-37.4761], Tmin=(1196.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.5786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1OOC1C=CO(7311)',
    structure = SMILES('[CH2]C1OOC1C=CO'),
    E0 = (20.2782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.21556,0.0666164,-2.2954e-05,-3.13783e-08,2.04949e-11,2590.33,27.2503], Tmin=(100,'K'), Tmax=(946.263,'K')), NASAPolynomial(coeffs=[22.5542,0.0131975,-3.28361e-06,5.63434e-10,-4.42212e-14,-3473.36,-88.9885], Tmin=(946.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxetane) + radical(CJCOOH)"""),
)

species(
    label = 'C=C[C](C=CO)OO(7312)',
    structure = SMILES('[CH2]C=C(C=CO)OO'),
    E0 = (-122.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02633,0.0886849,-9.03721e-05,4.42705e-08,-8.11954e-12,-14502.3,32.557], Tmin=(100,'K'), Tmax=(1526.08,'K')), NASAPolynomial(coeffs=[24.0308,0.010558,-1.34358e-06,3.27916e-11,3.11741e-15,-20700.4,-94.1712], Tmin=(1526.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC([C]=CO)OO(7313)',
    structure = SMILES('C=CC([C]=CO)OO'),
    E0 = (11.4707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973199,'amu*angstrom^2'), symmetry=1, barrier=(22.3758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975122,'amu*angstrom^2'), symmetry=1, barrier=(22.42,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97585,'amu*angstrom^2'), symmetry=1, barrier=(22.4367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973399,'amu*angstrom^2'), symmetry=1, barrier=(22.3804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973652,'amu*angstrom^2'), symmetry=1, barrier=(22.3862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288843,0.0832995,-7.43853e-05,2.48884e-08,-4.42352e-13,1543.83,32.3891], Tmin=(100,'K'), Tmax=(991.52,'K')), NASAPolynomial(coeffs=[21.9196,0.0157435,-5.52363e-06,1.00413e-09,-7.21722e-14,-3943.48,-80.0275], Tmin=(991.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.4707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C=CO)OO(7314)',
    structure = SMILES('C=[C]C(C=CO)OO'),
    E0 = (11.4707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.973199,'amu*angstrom^2'), symmetry=1, barrier=(22.3758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975122,'amu*angstrom^2'), symmetry=1, barrier=(22.42,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97585,'amu*angstrom^2'), symmetry=1, barrier=(22.4367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973399,'amu*angstrom^2'), symmetry=1, barrier=(22.3804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973652,'amu*angstrom^2'), symmetry=1, barrier=(22.3862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.288843,0.0832995,-7.43853e-05,2.48884e-08,-4.42352e-13,1543.83,32.3891], Tmin=(100,'K'), Tmax=(991.52,'K')), NASAPolynomial(coeffs=[21.9196,0.0157435,-5.52363e-06,1.00413e-09,-7.21722e-14,-3943.48,-80.0275], Tmin=(991.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.4707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(C=[C]O)OO(7315)',
    structure = SMILES('C=CC(C=[C]O)OO'),
    E0 = (13.3731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,196.644,196.695],'cm^-1')),
        HinderedRotor(inertia=(0.717383,'amu*angstrom^2'), symmetry=1, barrier=(19.7085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717032,'amu*angstrom^2'), symmetry=1, barrier=(19.709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716877,'amu*angstrom^2'), symmetry=1, barrier=(19.7084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717213,'amu*angstrom^2'), symmetry=1, barrier=(19.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.715781,'amu*angstrom^2'), symmetry=1, barrier=(19.7096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.188399,0.0820328,-8.14369e-05,3.99251e-08,-7.6156e-12,1767.72,34.6894], Tmin=(100,'K'), Tmax=(1284.54,'K')), NASAPolynomial(coeffs=[20.3188,0.0181744,-6.86746e-06,1.22421e-09,-8.3546e-14,-3500.75,-69.381], Tmin=(1284.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.3731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC(C=CO)OO(7316)',
    structure = SMILES('[CH]=CC(C=CO)OO'),
    E0 = (20.7251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.01974,'amu*angstrom^2'), symmetry=1, barrier=(23.4457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0167,'amu*angstrom^2'), symmetry=1, barrier=(23.376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01619,'amu*angstrom^2'), symmetry=1, barrier=(23.3641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01594,'amu*angstrom^2'), symmetry=1, barrier=(23.3584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01794,'amu*angstrom^2'), symmetry=1, barrier=(23.4043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.316025,0.0820629,-6.57016e-05,1.22947e-08,4.90282e-12,2659.63,32.4161], Tmin=(100,'K'), Tmax=(968.721,'K')), NASAPolynomial(coeffs=[23.2597,0.0135615,-4.29836e-06,7.76659e-10,-5.76541e-14,-3261.51,-87.5597], Tmin=(968.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.7251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = '[CH]=CC(C=C)O[O](2622)',
    structure = SMILES('[CH]=CC(C=C)O[O]'),
    E0 = (381.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.603227,'amu*angstrom^2'), symmetry=1, barrier=(13.8694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602694,'amu*angstrom^2'), symmetry=1, barrier=(13.8571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603105,'amu*angstrom^2'), symmetry=1, barrier=(13.8666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31348,0.0588054,-5.20174e-05,2.42616e-08,-4.65637e-12,45983.4,27.4687], Tmin=(100,'K'), Tmax=(1226.44,'K')), NASAPolynomial(coeffs=[11.5055,0.0255643,-1.13616e-05,2.16191e-09,-1.51506e-13,43483.5,-23.7822], Tmin=(1226.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'H(3)(3)',
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
    label = 'C=CC([CH]C=O)O[O](7127)',
    structure = SMILES('C=CC(C=C[O])O[O]'),
    E0 = (67.0962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.835242,'amu*angstrom^2'), symmetry=1, barrier=(19.2039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.835286,'amu*angstrom^2'), symmetry=1, barrier=(19.2049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834793,'amu*angstrom^2'), symmetry=1, barrier=(19.1935,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438301,0.0679004,-4.99817e-05,9.80629e-09,2.66748e-12,8207.2,31.4755], Tmin=(100,'K'), Tmax=(1014.94,'K')), NASAPolynomial(coeffs=[18.1121,0.0185466,-7.04297e-06,1.30894e-09,-9.36704e-14,3574.06,-59.2032], Tmin=(1014.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.0962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=COJ)"""),
)

species(
    label = 'CHCHOH(49)(49)',
    structure = SMILES('[CH]=CO'),
    E0 = (120.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,458.887,459.577,459.893,460.324],'cm^-1')),
        HinderedRotor(inertia=(0.000141718,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1868.27,'J/mol'), sigma=(4.162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99181,0.0473233,-6.66059e-05,4.68997e-08,-1.30686e-11,15200.1,31.4259], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.78246,0.00524797,-1.71857e-06,2.59722e-10,-1.48227e-14,12883.6,-21.0851], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(298,'K'), Tmax=(6000,'K'), E0=(120.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C[CH]O[O](2619)',
    structure = SMILES('C=C[CH]O[O]'),
    E0 = (199.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,316.03],'cm^-1')),
        HinderedRotor(inertia=(0.156402,'amu*angstrom^2'), symmetry=1, barrier=(11.1962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158062,'amu*angstrom^2'), symmetry=1, barrier=(11.1979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21608,0.0376502,-3.45386e-05,1.68249e-08,-3.30513e-12,24056.5,18.0455], Tmin=(100,'K'), Tmax=(1220.77,'K')), NASAPolynomial(coeffs=[9.39745,0.0141196,-5.62587e-06,1.03554e-09,-7.16583e-14,22303.1,-18.0331], Tmin=(1220.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C2H3(28)(29)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]O[CH]C=CO(7317)',
    structure = SMILES('[O]O[CH]C=CO'),
    E0 = (-14.9999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.12219,'amu*angstrom^2'), symmetry=1, barrier=(25.8013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12188,'amu*angstrom^2'), symmetry=1, barrier=(25.7941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12189,'amu*angstrom^2'), symmetry=1, barrier=(25.7944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41255,0.0421239,-1.65469e-06,-4.44282e-08,2.4907e-11,-1697.1,22.1858], Tmin=(100,'K'), Tmax=(913.122,'K')), NASAPolynomial(coeffs=[19.9116,0.000280065,2.70015e-06,-6.02097e-10,3.85777e-14,-6709.41,-74.3272], Tmin=(913.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-14.9999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C=C[C](C=CO)O[O](7318)',
    structure = SMILES('[CH2]C=C(C=CO)O[O]'),
    E0 = (29.7654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.927284,'amu*angstrom^2'), symmetry=1, barrier=(21.3201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927702,'amu*angstrom^2'), symmetry=1, barrier=(21.3297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922899,'amu*angstrom^2'), symmetry=1, barrier=(21.2193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925227,'amu*angstrom^2'), symmetry=1, barrier=(21.2728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.556183,0.0828208,-8.78495e-05,4.5187e-08,-8.70714e-12,3758.88,31.6178], Tmin=(100,'K'), Tmax=(1449.19,'K')), NASAPolynomial(coeffs=[21.7621,0.0100903,-1.05035e-06,-4.19022e-11,9.42282e-15,-1541.25,-80.3035], Tmin=(1449.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.7654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC([C]=CO)O[O](7319)',
    structure = SMILES('C=CC([C]=CO)O[O]'),
    E0 = (163.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.816308,'amu*angstrom^2'), symmetry=1, barrier=(18.7685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816587,'amu*angstrom^2'), symmetry=1, barrier=(18.7749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816317,'amu*angstrom^2'), symmetry=1, barrier=(18.7687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816037,'amu*angstrom^2'), symmetry=1, barrier=(18.7623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27508,0.082612,-8.90201e-05,4.67281e-08,-9.39164e-12,19825.2,33.1011], Tmin=(100,'K'), Tmax=(1267.51,'K')), NASAPolynomial(coeffs=[21.29,0.0125184,-3.65738e-06,5.61344e-10,-3.5577e-14,14522.2,-75.4042], Tmin=(1267.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C(C=CO)O[O](7320)',
    structure = SMILES('C=[C]C(C=CO)O[O]'),
    E0 = (163.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.816308,'amu*angstrom^2'), symmetry=1, barrier=(18.7685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816587,'amu*angstrom^2'), symmetry=1, barrier=(18.7749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816317,'amu*angstrom^2'), symmetry=1, barrier=(18.7687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.816037,'amu*angstrom^2'), symmetry=1, barrier=(18.7623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.27508,0.082612,-8.90201e-05,4.67281e-08,-9.39164e-12,19825.2,33.1011], Tmin=(100,'K'), Tmax=(1267.51,'K')), NASAPolynomial(coeffs=[21.29,0.0125184,-3.65738e-06,5.61344e-10,-3.5577e-14,14522.2,-75.4042], Tmin=(1267.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = 'C=CC(C=[C]O)O[O](7321)',
    structure = SMILES('C=CC(C=[C]O)O[O]'),
    E0 = (165.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,492.5,1135,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.659057,'amu*angstrom^2'), symmetry=1, barrier=(15.153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65986,'amu*angstrom^2'), symmetry=1, barrier=(15.1715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658688,'amu*angstrom^2'), symmetry=1, barrier=(15.1445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658118,'amu*angstrom^2'), symmetry=1, barrier=(15.1314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.232554,0.0766692,-8.03485e-05,4.23417e-08,-8.71809e-12,20031.3,33.9328], Tmin=(100,'K'), Tmax=(1190.28,'K')), NASAPolynomial(coeffs=[17.564,0.0184262,-6.95067e-06,1.23242e-09,-8.37701e-14,15905.4,-52.7006], Tmin=(1190.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC(C=CO)O[O](7322)',
    structure = SMILES('[CH]=CC(C=CO)O[O]'),
    E0 = (172.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.873818,'amu*angstrom^2'), symmetry=1, barrier=(20.0908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87399,'amu*angstrom^2'), symmetry=1, barrier=(20.0948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.873884,'amu*angstrom^2'), symmetry=1, barrier=(20.0923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874142,'amu*angstrom^2'), symmetry=1, barrier=(20.0982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.529529,0.0839975,-8.9257e-05,4.53846e-08,-8.72092e-12,20951,33.9475], Tmin=(100,'K'), Tmax=(1397.96,'K')), NASAPolynomial(coeffs=[23.0484,0.00960847,-2.00756e-06,2.33017e-10,-1.26692e-14,15035.5,-85.2809], Tmin=(1397.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1[CH]C(O)OO1(7323)',
    structure = SMILES('C=CC1[CH]C(O)OO1'),
    E0 = (-56.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03075,0.047502,2.30603e-05,-7.57959e-08,3.67907e-11,-6731.23,29.2158], Tmin=(100,'K'), Tmax=(906.734,'K')), NASAPolynomial(coeffs=[19.0555,0.0152985,-1.93271e-06,1.24922e-10,-7.80083e-15,-11944.9,-66.7032], Tmin=(906.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'OC=CC1[CH]COO1(7304)',
    structure = SMILES('OC=CC1[CH]COO1'),
    E0 = (-65.5556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632245,0.048697,4.47556e-05,-1.18543e-07,5.78327e-11,-7739.21,28.0429], Tmin=(100,'K'), Tmax=(888.409,'K')), NASAPolynomial(coeffs=[26.4697,0.00192401,6.28422e-06,-1.54405e-09,1.09004e-13,-15075.1,-109], Tmin=(888.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.5556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'C=CC=C=CO(6203)',
    structure = SMILES('C=CC=C=CO'),
    E0 = (25.7073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.24386,'amu*angstrom^2'), symmetry=1, barrier=(28.5988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24416,'amu*angstrom^2'), symmetry=1, barrier=(28.6057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02364,0.048058,2.85443e-06,-5.66469e-08,3.08127e-11,3215.32,19.3983], Tmin=(100,'K'), Tmax=(909.916,'K')), NASAPolynomial(coeffs=[22.142,0.00160986,2.95307e-06,-6.91045e-10,4.50216e-14,-2548.21,-91.0442], Tmin=(909.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.7073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CC=CO(6189)',
    structure = SMILES('C=C=CC=CO'),
    E0 = (25.7073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.24386,'amu*angstrom^2'), symmetry=1, barrier=(28.5988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24416,'amu*angstrom^2'), symmetry=1, barrier=(28.6057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02364,0.048058,2.85443e-06,-5.66469e-08,3.08127e-11,3215.32,19.3983], Tmin=(100,'K'), Tmax=(909.916,'K')), NASAPolynomial(coeffs=[22.142,0.00160986,2.95307e-06,-6.91045e-10,4.50216e-14,-2548.21,-91.0442], Tmin=(909.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.7073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(CC=O)O[O](7324)',
    structure = SMILES('C=CC(CC=O)O[O]'),
    E0 = (-72.3984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,2782.5,750,1395,475,1775,1000,353.04,353.04],'cm^-1')),
        HinderedRotor(inertia=(0.109645,'amu*angstrom^2'), symmetry=1, barrier=(9.69757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109645,'amu*angstrom^2'), symmetry=1, barrier=(9.69757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109645,'amu*angstrom^2'), symmetry=1, barrier=(9.69757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109645,'amu*angstrom^2'), symmetry=1, barrier=(9.69757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04483,0.0695462,-6.7204e-05,3.72372e-08,-8.84566e-12,-8604.9,29.7646], Tmin=(100,'K'), Tmax=(984.825,'K')), NASAPolynomial(coeffs=[9.14489,0.0366469,-1.70949e-05,3.31663e-09,-2.34895e-13,-10200.3,-9.18967], Tmin=(984.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.3984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'O(4)(4)',
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
    label = 'C=CC([O])C=CO(7325)',
    structure = SMILES('C=CC([O])C=CO'),
    E0 = (-67.5153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.921306,'amu*angstrom^2'), symmetry=1, barrier=(21.1826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920284,'amu*angstrom^2'), symmetry=1, barrier=(21.1592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.9211,'amu*angstrom^2'), symmetry=1, barrier=(21.1779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656471,0.0569123,-9.55029e-06,-4.15863e-08,2.36966e-11,-7984.58,27.5687], Tmin=(100,'K'), Tmax=(939.754,'K')), NASAPolynomial(coeffs=[21.389,0.0104166,-1.97682e-06,3.16563e-10,-2.7257e-14,-13724.9,-80.9743], Tmin=(939.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.5153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1OC=CC1OO(7326)',
    structure = SMILES('[CH2]C1OC=CC1OO'),
    E0 = (-55.7385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.358509,0.0602682,2.23335e-06,-6.47907e-08,3.5083e-11,-6554.17,23.7053], Tmin=(100,'K'), Tmax=(911.417,'K')), NASAPolynomial(coeffs=[24.134,0.00855494,7.2148e-07,-3.24956e-10,2.06631e-14,-13074.1,-100.785], Tmin=(911.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.7385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CC(C=C=O)OO(7327)',
    structure = SMILES('C=CC(C=C=O)OO'),
    E0 = (-108.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.767141,'amu*angstrom^2'), symmetry=1, barrier=(17.6381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766964,'amu*angstrom^2'), symmetry=1, barrier=(17.634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767061,'amu*angstrom^2'), symmetry=1, barrier=(17.6362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7213,'amu*angstrom^2'), symmetry=1, barrier=(39.576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.62386,0.0840019,-0.000100833,4.65803e-08,5.07864e-12,-12896.3,26.1513], Tmin=(100,'K'), Tmax=(570.48,'K')), NASAPolynomial(coeffs=[10.3855,0.0336501,-1.60139e-05,3.05624e-09,-2.11505e-13,-14304.5,-18.0445], Tmin=(570.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC([O])C=C[O](7328)',
    structure = SMILES('C=CC([O])C=C[O]'),
    E0 = (73.9473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,252.813,254.996,255.279,258.365],'cm^-1')),
        HinderedRotor(inertia=(0.478068,'amu*angstrom^2'), symmetry=1, barrier=(22.2341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.482104,'amu*angstrom^2'), symmetry=1, barrier=(22.2014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06895,0.049925,-6.06005e-06,-3.55219e-08,1.90346e-11,9012.46,27.4986], Tmin=(100,'K'), Tmax=(965.867,'K')), NASAPolynomial(coeffs=[18.1508,0.0140799,-4.58732e-06,8.68398e-10,-6.66198e-14,4084.94,-62.7445], Tmin=(965.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.9473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'CHCHO(45)(45)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C[CH]OO(7329)',
    structure = SMILES('C=C[CH]OO'),
    E0 = (47.4703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.577285,'amu*angstrom^2'), symmetry=1, barrier=(13.2729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0696899,'amu*angstrom^2'), symmetry=1, barrier=(27.4887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19648,'amu*angstrom^2'), symmetry=1, barrier=(27.5094,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75473,0.0434201,-3.68207e-05,1.56822e-08,-2.64238e-12,5795.15,18.9538], Tmin=(100,'K'), Tmax=(1428.93,'K')), NASAPolynomial(coeffs=[12.5007,0.0133387,-5.24287e-06,9.49426e-10,-6.47747e-14,2724.13,-36.7246], Tmin=(1428.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.4703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CCJO)"""),
)

species(
    label = '[O]C=C[CH]OO(7330)',
    structure = SMILES('O=C[CH][CH]OO'),
    E0 = (107.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3050,390,425,1340,1360,335,370,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.29733,'amu*angstrom^2'), symmetry=1, barrier=(6.8362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296448,'amu*angstrom^2'), symmetry=1, barrier=(6.81593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00276163,'amu*angstrom^2'), symmetry=1, barrier=(6.81912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27018,'amu*angstrom^2'), symmetry=1, barrier=(52.1958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29753,0.0702206,-0.000130518,1.25515e-07,-4.52445e-11,12972.3,23.6343], Tmin=(100,'K'), Tmax=(863.388,'K')), NASAPolynomial(coeffs=[5.25405,0.0274692,-1.38164e-05,2.64251e-09,-1.80005e-13,13199.3,10.3988], Tmin=(863.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJCO) + radical(CCsJOOH)"""),
)

species(
    label = 'C=C[C](C=C[O])OO(7331)',
    structure = SMILES('[CH2]C=C(C=C[O])OO'),
    E0 = (19.2233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,350,440,435,1725,195.675,195.68],'cm^-1')),
        HinderedRotor(inertia=(1.05578,'amu*angstrom^2'), symmetry=1, barrier=(28.6881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05588,'amu*angstrom^2'), symmetry=1, barrier=(28.688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05576,'amu*angstrom^2'), symmetry=1, barrier=(28.688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05579,'amu*angstrom^2'), symmetry=1, barrier=(28.6879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24957,0.0717151,-5.28315e-05,7.27808e-09,5.15717e-12,2456.77,29.3752], Tmin=(100,'K'), Tmax=(963.047,'K')), NASAPolynomial(coeffs=[19.3787,0.0167924,-5.49285e-06,9.56321e-10,-6.75213e-14,-2365.19,-68.0972], Tmin=(963.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.2233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC([C]=C[O])OO(7332)',
    structure = SMILES('C=CC([C]=C[O])OO'),
    E0 = (152.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,182.399,182.476,182.637],'cm^-1')),
        HinderedRotor(inertia=(1.01092,'amu*angstrom^2'), symmetry=1, barrier=(23.8341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00828,'amu*angstrom^2'), symmetry=1, barrier=(23.8334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00605,'amu*angstrom^2'), symmetry=1, barrier=(23.8352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0079,'amu*angstrom^2'), symmetry=1, barrier=(23.8336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0165855,0.077819,-7.54149e-05,3.56287e-08,-6.54345e-12,18547.1,32.8311], Tmin=(100,'K'), Tmax=(1331.64,'K')), NASAPolynomial(coeffs=[20.2913,0.0168171,-6.69975e-06,1.2271e-09,-8.48686e-14,13138.6,-70.9591], Tmin=(1331.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C=C[O])OO(7333)',
    structure = SMILES('C=[C]C(C=C[O])OO'),
    E0 = (152.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,182.399,182.476,182.637],'cm^-1')),
        HinderedRotor(inertia=(1.01092,'amu*angstrom^2'), symmetry=1, barrier=(23.8341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00828,'amu*angstrom^2'), symmetry=1, barrier=(23.8334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00605,'amu*angstrom^2'), symmetry=1, barrier=(23.8352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0079,'amu*angstrom^2'), symmetry=1, barrier=(23.8336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0165855,0.077819,-7.54149e-05,3.56287e-08,-6.54345e-12,18547.1,32.8311], Tmin=(100,'K'), Tmax=(1331.64,'K')), NASAPolynomial(coeffs=[20.2913,0.0168171,-6.69975e-06,1.2271e-09,-8.48686e-14,13138.6,-70.9591], Tmin=(1331.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC(C=C[O])OO(7334)',
    structure = SMILES('[CH]=CC(C=C[O])OO'),
    E0 = (162.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08412,'amu*angstrom^2'), symmetry=1, barrier=(24.926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08365,'amu*angstrom^2'), symmetry=1, barrier=(24.9153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0836,'amu*angstrom^2'), symmetry=1, barrier=(24.914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08405,'amu*angstrom^2'), symmetry=1, barrier=(24.9245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0690303,0.0753674,-6.30554e-05,1.91425e-08,6.83373e-14,19657.9,32.4462], Tmin=(100,'K'), Tmax=(1032.85,'K')), NASAPolynomial(coeffs=[20.3137,0.0167415,-6.63574e-06,1.26493e-09,-9.18044e-14,14421,-70.9837], Tmin=(1032.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C=[C][O])OO(7335)',
    structure = SMILES('C=CC([CH][C]=O)OO'),
    E0 = (135.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,373.197,373.197],'cm^-1')),
        HinderedRotor(inertia=(0.000944745,'amu*angstrom^2'), symmetry=1, barrier=(10.7267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246002,'amu*angstrom^2'), symmetry=1, barrier=(24.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282501,'amu*angstrom^2'), symmetry=1, barrier=(27.9204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108533,'amu*angstrom^2'), symmetry=1, barrier=(10.7267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282499,'amu*angstrom^2'), symmetry=1, barrier=(27.9204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787601,0.0723492,-7.44551e-05,4.03924e-08,-8.92467e-12,16406.3,33.2364], Tmin=(100,'K'), Tmax=(1082.93,'K')), NASAPolynomial(coeffs=[12.818,0.0279128,-1.29049e-05,2.5013e-09,-1.77338e-13,13800.7,-25.7618], Tmin=(1082.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC(OO)C1[CH]O1(7336)',
    structure = SMILES('C=CC(OO)C1[CH]O1'),
    E0 = (59.2119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.770044,0.0839297,-8.25157e-05,4.06494e-08,-7.52722e-12,7311.28,33.0797], Tmin=(100,'K'), Tmax=(1535.9,'K')), NASAPolynomial(coeffs=[20.4884,0.0151199,-2.18302e-06,8.08931e-11,3.90375e-15,2367.03,-73.4395], Tmin=(1535.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCsJO)"""),
)

species(
    label = 'OOC1[CH]COC=C1(7337)',
    structure = SMILES('OOC1[CH]COC=C1'),
    E0 = (-55.2732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901038,0.0473732,3.00093e-05,-8.39279e-08,3.88774e-11,-6517.12,24.9568], Tmin=(100,'K'), Tmax=(932.93,'K')), NASAPolynomial(coeffs=[20.8236,0.0149064,-2.92822e-06,4.48808e-10,-3.66854e-14,-12538.8,-82.1249], Tmin=(932.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.2732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CCJCOOH)"""),
)

species(
    label = 'C=CC1C=COO1(7338)',
    structure = SMILES('C=CC1C=COO1'),
    E0 = (63.5798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76332,0.0305859,4.36057e-05,-8.15802e-08,3.37742e-11,7744.15,23.8196], Tmin=(100,'K'), Tmax=(971.171,'K')), NASAPolynomial(coeffs=[16.3464,0.016559,-5.83444e-06,1.16891e-09,-9.19988e-14,2740.59,-57.2859], Tmin=(971.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.5798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolene)"""),
)

species(
    label = '[CH]=CC(C=C)OO(2618)',
    structure = SMILES('[CH]=CC(C=C)OO'),
    E0 = (229.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.99082,'amu*angstrom^2'), symmetry=1, barrier=(22.7809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139979,'amu*angstrom^2'), symmetry=1, barrier=(14.7244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989926,'amu*angstrom^2'), symmetry=1, barrier=(22.7603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99012,'amu*angstrom^2'), symmetry=1, barrier=(22.7648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885317,0.0642782,-5.35503e-05,2.24348e-08,-3.78936e-12,27720,28.2489], Tmin=(100,'K'), Tmax=(1397.88,'K')), NASAPolynomial(coeffs=[14.8117,0.0244282,-1.07892e-05,2.04154e-09,-1.42186e-13,23826.5,-43.6028], Tmin=(1397.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1C(C=O)C1OO(7339)',
    structure = SMILES('[CH2]C1C(C=O)C1OO'),
    E0 = (1.02761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.555589,0.066533,-4.58079e-05,9.66992e-09,1.69707e-12,255.597,28.8955], Tmin=(100,'K'), Tmax=(1034.03,'K')), NASAPolynomial(coeffs=[15.5876,0.024928,-9.45396e-06,1.70465e-09,-1.18083e-13,-3737.61,-48.4054], Tmin=(1034.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.02761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC(=CC=O)OO(7340)',
    structure = SMILES('C=CC(=CC=O)OO'),
    E0 = (-118.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.773573,'amu*angstrom^2'), symmetry=1, barrier=(17.786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.774129,'amu*angstrom^2'), symmetry=1, barrier=(17.7987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773446,'amu*angstrom^2'), symmetry=1, barrier=(17.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773597,'amu*angstrom^2'), symmetry=1, barrier=(17.7865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851437,0.0758053,-8.47487e-05,5.20451e-08,-1.33985e-11,-14200.3,25.6827], Tmin=(100,'K'), Tmax=(922.369,'K')), NASAPolynomial(coeffs=[10.362,0.0345603,-1.76725e-05,3.56275e-09,-2.57486e-13,-15954.7,-19.4313], Tmin=(922.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'O=CC=COO(7341)',
    structure = SMILES('O=CC=COO'),
    E0 = (-165.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.893552,'amu*angstrom^2'), symmetry=1, barrier=(20.5445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.893859,'amu*angstrom^2'), symmetry=1, barrier=(20.5516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.893539,'amu*angstrom^2'), symmetry=1, barrier=(20.5442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9884,0.0465945,-4.23753e-05,1.92808e-08,-3.60921e-12,-19851.2,19.5603], Tmin=(100,'K'), Tmax=(1242.75,'K')), NASAPolynomial(coeffs=[10.237,0.0200451,-1.03304e-05,2.09062e-09,-1.51153e-13,-21901.4,-22.0273], Tmin=(1242.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C[C](CC=O)OO(7342)',
    structure = SMILES('[CH2]C=C(CC=O)OO'),
    E0 = (-67.7289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.18945,'amu*angstrom^2'), symmetry=1, barrier=(4.35582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94352,'amu*angstrom^2'), symmetry=1, barrier=(21.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0212669,'amu*angstrom^2'), symmetry=1, barrier=(21.6924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.943382,'amu*angstrom^2'), symmetry=1, barrier=(21.6902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0212957,'amu*angstrom^2'), symmetry=1, barrier=(21.6865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847733,0.0656244,-4.74457e-05,1.58907e-08,-2.11217e-12,-8030.04,29.1779], Tmin=(100,'K'), Tmax=(1733.6,'K')), NASAPolynomial(coeffs=[18.4157,0.025089,-1.23722e-05,2.40287e-09,-1.67112e-13,-14121.2,-65.2435], Tmin=(1733.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.7289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC(C[C]=O)OO(7343)',
    structure = SMILES('C=CC(C[C]=O)OO'),
    E0 = (-64.4425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1855,455,950,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,352.454,352.455],'cm^-1')),
        HinderedRotor(inertia=(0.248513,'amu*angstrom^2'), symmetry=1, barrier=(21.9069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109568,'amu*angstrom^2'), symmetry=1, barrier=(9.6587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109568,'amu*angstrom^2'), symmetry=1, barrier=(9.65873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115719,'amu*angstrom^2'), symmetry=1, barrier=(10.2009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248513,'amu*angstrom^2'), symmetry=1, barrier=(21.9069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.67681,0.0765099,-7.91173e-05,4.43925e-08,-1.02951e-11,-7633.78,31.1423], Tmin=(100,'K'), Tmax=(1026.69,'K')), NASAPolynomial(coeffs=[11.8814,0.0328573,-1.53412e-05,2.98089e-09,-2.11438e-13,-9934.54,-23.2086], Tmin=(1026.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.4425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C(CC=O)OO(7344)',
    structure = SMILES('C=[C]C(CC=O)OO'),
    E0 = (13.4387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,357.37,357.37],'cm^-1')),
        HinderedRotor(inertia=(0.393225,'amu*angstrom^2'), symmetry=1, barrier=(35.6373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962822,'amu*angstrom^2'), symmetry=1, barrier=(8.72588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962821,'amu*angstrom^2'), symmetry=1, barrier=(8.72588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962824,'amu*angstrom^2'), symmetry=1, barrier=(8.72588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222971,'amu*angstrom^2'), symmetry=1, barrier=(20.2075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75045,0.0776018,-8.6341e-05,5.53913e-08,-1.50708e-11,1728.04,30.5431], Tmin=(100,'K'), Tmax=(873.223,'K')), NASAPolynomial(coeffs=[9.36257,0.0381534,-1.85798e-05,3.66043e-09,-2.61028e-13,223.928,-9.83819], Tmin=(873.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.4387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(CC=O)OO(7345)',
    structure = SMILES('[CH]=CC(CC=O)OO'),
    E0 = (22.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,315.906],'cm^-1')),
        HinderedRotor(inertia=(0.356917,'amu*angstrom^2'), symmetry=1, barrier=(27.4272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133243,'amu*angstrom^2'), symmetry=1, barrier=(9.40498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160487,'amu*angstrom^2'), symmetry=1, barrier=(11.3069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128084,'amu*angstrom^2'), symmetry=1, barrier=(9.3954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373366,'amu*angstrom^2'), symmetry=1, barrier=(27.4244,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74245,0.0762422,-7.76669e-05,4.33073e-08,-1.01014e-11,2842.85,30.4942], Tmin=(100,'K'), Tmax=(1013.54,'K')), NASAPolynomial(coeffs=[11.1518,0.03516,-1.68656e-05,3.31366e-09,-2.36369e-13,732.828,-19.8649], Tmin=(1013.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'O=CC1C[CH]C1OO(7346)',
    structure = SMILES('O=CC1C[CH]C1OO'),
    E0 = (0.982394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767617,0.0610525,-3.46427e-05,3.90184e-09,1.91137e-12,242.863,29.6548], Tmin=(100,'K'), Tmax=(1166.23,'K')), NASAPolynomial(coeffs=[14.8303,0.0281638,-1.20774e-05,2.28448e-09,-1.60386e-13,-4080.69,-44.826], Tmin=(1166.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.982394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCOOH)"""),
)

species(
    label = 'C=CC1OC1C=O(7274)',
    structure = SMILES('C=CC1OC1C=O'),
    E0 = (-116.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80385,0.0294689,5.26188e-05,-1.00687e-07,4.49928e-11,-13902.8,23.9649], Tmin=(100,'K'), Tmax=(900.611,'K')), NASAPolynomial(coeffs=[17.2226,0.0109598,2.16104e-07,-2.86516e-10,2.04043e-14,-18706.6,-60.0588], Tmin=(900.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide)"""),
)

species(
    label = 'C=CC([CH]OO)C=O(7347)',
    structure = SMILES('C=CC([CH]OO)C=O'),
    E0 = (-29.7367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,318.985,318.985],'cm^-1')),
        HinderedRotor(inertia=(0.648152,'amu*angstrom^2'), symmetry=1, barrier=(46.8005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989153,'amu*angstrom^2'), symmetry=1, barrier=(7.14265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989256,'amu*angstrom^2'), symmetry=1, barrier=(7.14405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989236,'amu*angstrom^2'), symmetry=1, barrier=(7.1427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648133,'amu*angstrom^2'), symmetry=1, barrier=(46.7997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.371224,0.0868126,-0.000121597,1.02351e-07,-3.55773e-11,-3452.52,30.1896], Tmin=(100,'K'), Tmax=(773.633,'K')), NASAPolynomial(coeffs=[8.0553,0.0406649,-1.96779e-05,3.80041e-09,-2.65344e-13,-4449.4,-3.66819], Tmin=(773.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.7367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOOH)"""),
)

species(
    label = 'C=CC([O])C(O)C=O(7348)',
    structure = SMILES('C=CC([O])C(O)C=O'),
    E0 = (-235.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.628821,0.0617216,-2.9539e-05,-9.58627e-09,8.73739e-12,-28183.3,34.3448], Tmin=(100,'K'), Tmax=(1020.65,'K')), NASAPolynomial(coeffs=[17.5293,0.0214623,-8.54586e-06,1.63598e-09,-1.18922e-13,-32986.2,-54.1636], Tmin=(1020.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-235.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
)

species(
    label = 'HCO(14)(15)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH]C(C=C)OO(7349)',
    structure = SMILES('[CH]C(C=C)OO'),
    E0 = (330.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,370.383,370.545,1574.99,1575.07],'cm^-1')),
        HinderedRotor(inertia=(0.10756,'amu*angstrom^2'), symmetry=1, barrier=(10.4728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252793,'amu*angstrom^2'), symmetry=1, barrier=(24.6289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252883,'amu*angstrom^2'), symmetry=1, barrier=(24.6321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252711,'amu*angstrom^2'), symmetry=1, barrier=(24.6317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03126,0.0569342,-4.91055e-05,2.09313e-08,-3.51543e-12,39850.8,25.778], Tmin=(100,'K'), Tmax=(1436.62,'K')), NASAPolynomial(coeffs=[15.656,0.0162146,-6.58991e-06,1.20202e-09,-8.21755e-14,35648.7,-50.0765], Tmin=(1436.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (-41.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (35.5786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (47.6766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (83.2026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (71.2366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (71.2366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (46.4133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (37.0786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (4.32405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (409.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (278.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (323.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (274.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (241.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (375.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (375.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (377.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (384.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-16.2088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-18.3298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (104.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (104.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (69.5004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (175.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-4.16563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (136.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (176.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (161.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (28.8125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (44.5109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (138.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (102.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (111.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (278.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (296.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (396.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (231.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (365.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (364.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (373.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (347.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (128.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-0.199202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (94.5038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (472.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (23.2475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (102.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (129.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (-18.8344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (90.4148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (92.4308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (154.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (29.7469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (67.0016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (40.2629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (-39.2652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (130.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (16.7625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (362.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'C5H7O(224)(223)'],
    products = ['S(226)(225)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.9277e+07,'m^3/(mol*s)'), n=-0.00388889, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CdCd;Y_rad] for rate rule [C_rad/H/CdCd;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(226)(225)'],
    products = ['C=CC1OOC1[CH]O(7310)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(109.945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 109.8 to 109.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(226)(225)'],
    products = ['[CH2]C1OOC1C=CO(7311)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['S(226)(225)'],
    products = ['C=C[C](C=CO)OO(7312)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['S(226)(225)'],
    products = ['C=CC([C]=CO)OO(7313)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['S(226)(225)'],
    products = ['C=[C]C(C=CO)OO(7314)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=CC(C=[C]O)OO(7315)'],
    products = ['S(226)(225)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 3.31662479036
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['S(226)(225)'],
    products = ['[CH]=CC(C=CO)OO(7316)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['S(226)(225)'],
    products = ['S(233)(232)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(146928,'s^-1'), n=1.47286, Ea=(78.6904,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6H;O_rad_out;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)(5)', '[CH]=CC(C=C)O[O](2622)'],
    products = ['S(226)(225)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.44289e+07,'m^3/(mol*s)'), n=0.0225, Ea=(0.0523,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_pri_rad] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)(3)', 'C=CC([CH]C=O)O[O](7127)'],
    products = ['S(226)(225)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHCHOH(49)(49)', 'C=C[CH]O[O](2619)'],
    products = ['S(226)(225)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/OneDe] for rate rule [Cd_pri_rad;C_rad/H/OneDeO]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C2H3(28)(29)', '[O]O[CH]C=CO(7317)'],
    products = ['S(226)(225)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)(3)', 'C=C[C](C=CO)O[O](7318)'],
    products = ['S(226)(225)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', 'C=CC([C]=CO)O[O](7319)'],
    products = ['S(226)(225)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', 'C=[C]C(C=CO)O[O](7320)'],
    products = ['S(226)(225)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)(3)', 'C=CC(C=[C]O)O[O](7321)'],
    products = ['S(226)(225)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)(3)', '[CH]=CC(C=CO)O[O](7322)'],
    products = ['S(226)(225)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(226)(225)'],
    products = ['C=CC1[CH]C(O)OO1(7323)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.94e+07,'s^-1'), n=0.93, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_HNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(226)(225)'],
    products = ['OC=CC1[CH]COO1(7304)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['S(226)(225)'],
    products = ['HO2(8)(9)', 'C=CC=C=CO(6203)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(226)(225)'],
    products = ['HO2(8)(9)', 'C=C=CC=CO(6189)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(226)(225)'],
    products = ['C=CC(CC=O)O[O](7324)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)(4)', 'C=CC([O])C=CO(7325)'],
    products = ['S(226)(225)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['S(233)(232)'],
    products = ['[CH2]C1OC=CC1OO(7326)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.249e+08,'s^-1'), n=0.846, Ea=(80.7428,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)(3)', 'C=CC(C=C=O)OO(7327)'],
    products = ['S(233)(232)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=CC(C=[C]O)OO(7315)'],
    products = ['S(233)(232)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=CC([C]=CO)OO(7313)'],
    products = ['S(233)(232)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['S(233)(232)'],
    products = ['C=C[C](C=CO)OO(7312)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4H_SDS;O_rad_out;Cs_H_out] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(C=CO)OO(7314)'],
    products = ['S(233)(232)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Cd_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=CC(C=CO)OO(7316)'],
    products = ['S(233)(232)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(173703,'s^-1'), n=1.89007, Ea=(118.15,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['OH(5)(5)', 'C=CC([O])C=C[O](7328)'],
    products = ['S(233)(232)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 103 used for O_pri_rad;O_rad/NonDe
Exact match found for rate rule [O_pri_rad;O_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['HO2(8)(9)', 'C5H6O(219)(218)'],
    products = ['S(233)(232)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.03e+12,'cm^3/(mol*s)','*|/',5), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [O_rad/NonDe;C_sec_rad] for rate rule [O_rad/NonDe;C_rad/H/CdCd]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)(3)', 'C=CC([CH]C=O)O[O](7127)'],
    products = ['S(233)(232)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['CHCHO(45)(45)', 'C=C[CH]OO(7329)'],
    products = ['S(233)(232)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;Cd_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C2H3(28)(29)', '[O]C=C[CH]OO(7330)'],
    products = ['S(233)(232)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)(3)', 'C=C[C](C=C[O])OO(7331)'],
    products = ['S(233)(232)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)(3)', 'C=CC([C]=C[O])OO(7332)'],
    products = ['S(233)(232)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(3)(3)', 'C=[C]C(C=C[O])OO(7333)'],
    products = ['S(233)(232)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)(3)', '[CH]=CC(C=C[O])OO(7334)'],
    products = ['S(233)(232)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(3)(3)', 'C=CC(C=[C][O])OO(7335)'],
    products = ['S(233)(232)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['S(233)(232)'],
    products = ['C=CC(OO)C1[CH]O1(7336)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(9.85157e+11,'s^-1'), n=0.224969, Ea=(213.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;doublebond_intra_pri_HNd_Cs;radadd_intra] + [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['S(233)(232)'],
    products = ['OOC1[CH]COC=C1(7337)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(84.7093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['S(233)(232)'],
    products = ['OH(5)(5)', 'C=CC1C=COO1(7338)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(179.412,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OOH] for rate rule [R4OO_SDS;Y_rad_intra;OOH]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O(4)(4)', '[CH]=CC(C=C)OO(2618)'],
    products = ['S(233)(232)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction46',
    reactants = ['S(233)(232)'],
    products = ['[CH2]C1C(C=O)C1OO(7339)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(108.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_csHDe]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['H(3)(3)', 'C=CC(=CC=O)OO(7340)'],
    products = ['S(233)(232)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CdOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C2H3(28)(29)', 'O=CC=COO(7341)'],
    products = ['S(233)(232)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(0.00977588,'m^3/(mol*s)'), n=2.40996, Ea=(8.29121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ-H] for rate rule [Cds-OsH_Cds;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction49',
    reactants = ['HO2(8)(9)', 'C5H6O(217)(216)'],
    products = ['S(233)(232)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(0.910519,'m^3/(mol*s)'), n=2.0547, Ea=(8.05595,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-OneDeH;YJ] for rate rule [Cds-CdH_Cds-COH;OJ-O2s]
Euclidian distance = 5.09901951359
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C[C](CC=O)OO(7342)'],
    products = ['S(233)(232)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.4044e+07,'s^-1'), n=1.72645, Ea=(158.144,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_noH;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_OneDe;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_OneDe/O;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=CC(C[C]=O)OO(7343)'],
    products = ['S(233)(232)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=[C]C(CC=O)OO(7344)'],
    products = ['S(233)(232)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=CC(CC=O)O[O](7324)'],
    products = ['S(233)(232)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]=CC(CC=O)OO(7345)'],
    products = ['S(233)(232)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['S(233)(232)'],
    products = ['O=CC1C[CH]C1OO(7346)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3.72602e+07,'s^-1'), n=1.21458, Ea=(125.171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_csHCO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['S(233)(232)'],
    products = ['OH(5)(5)', 'C=CC1OC1C=O(7274)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.38e+12,'s^-1','*|/',1.2), n=0, Ea=(45.6433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOH] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOH]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=CC([CH]OO)C=O(7347)'],
    products = ['S(233)(232)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction58',
    reactants = ['S(233)(232)'],
    products = ['C=CC([O])C(O)C=O(7348)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(2.69e+11,'s^-1'), n=0, Ea=(101.671,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OOH_S;C_rad_out_1H] for rate rule [R2OOH_S;C_rad_out_H/OneDe]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['HCO(14)(15)', '[CH]C(C=C)OO(7349)'],
    products = ['S(233)(232)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '601',
    isomers = [
        'S(226)(225)',
        'S(233)(232)',
    ],
    reactants = [
        ('O2(2)(2)', 'C5H7O(224)(223)'),
        ('HO2(8)(9)', 'C5H6O(219)(218)'),
        ('HO2(8)(9)', 'C5H6O(217)(216)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '601',
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

