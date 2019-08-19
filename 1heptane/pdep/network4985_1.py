species(
    label = 'O=CCC[CH]OC1[CH]CC1(22017)',
    structure = SMILES('O=CCC[CH]OC1[CH]CC1'),
    E0 = (74.5989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05483,0.095959,-6.98693e-05,2.55588e-08,-3.74931e-12,9166.85,41.337], Tmin=(100,'K'), Tmax=(1606.83,'K')), NASAPolynomial(coeffs=[22.8074,0.0365569,-1.44166e-05,2.55163e-09,-1.69721e-13,1498.34,-85.1016], Tmin=(1606.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.5989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'C1=CCC1(2460)',
    structure = SMILES('C1=CCC1'),
    E0 = (144.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,204.189,1057.03,1057.56,1057.77,1058.14,1058.14,1058.17,1058.23,1058.24,1058.29,1058.3,3484.72],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3170.39,'J/mol'), sigma=(5.50295,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.21 K, Pc=43.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39621,-0.00343725,9.09465e-05,-1.13275e-07,4.24111e-11,17412.4,10.7744], Tmin=(100,'K'), Tmax=(950.238,'K')), NASAPolynomial(coeffs=[9.42742,0.0139644,-4.06892e-06,7.75034e-10,-6.20421e-14,14334.3,-28.18], Tmin=(950.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'O=CCCC=O(5767)',
    structure = SMILES('O=CCCC=O'),
    E0 = (-309.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.18601,'amu*angstrom^2'), symmetry=1, barrier=(4.27673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186498,'amu*angstrom^2'), symmetry=1, barrier=(4.28796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186203,'amu*angstrom^2'), symmetry=1, barrier=(4.28117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3653.08,'J/mol'), sigma=(5.8998,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.60 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98042,0.0494518,-6.06096e-05,5.52155e-08,-2.1639e-11,-37204.9,20.2227], Tmin=(100,'K'), Tmax=(774.657,'K')), NASAPolynomial(coeffs=[2.99129,0.0355848,-1.7014e-05,3.28724e-09,-2.30123e-13,-37102,17.2786], Tmin=(774.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH)"""),
)

species(
    label = '[O]C1CCC1OC1[CH]CC1(31316)',
    structure = SMILES('[O]C1CCC1OC1[CH]CC1'),
    E0 = (180.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.198954,0.0549421,6.63492e-05,-1.25203e-07,5.11721e-11,21855.6,36.0224], Tmin=(100,'K'), Tmax=(979.253,'K')), NASAPolynomial(coeffs=[21.8571,0.0368308,-1.36789e-05,2.6483e-09,-1.98849e-13,14240.5,-85.2356], Tmin=(979.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[O]C1CC[CH]OC2CCC12(31307)',
    structure = SMILES('[O]C1CC[CH]OC2CCC12'),
    E0 = (-62.0786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.70964,0.118701,-0.000109566,5.18786e-08,-9.84299e-12,-7255.1,2.75366], Tmin=(100,'K'), Tmax=(1264.57,'K')), NASAPolynomial(coeffs=[22.8895,0.0408909,-1.7271e-05,3.22176e-09,-2.23804e-13,-13476.6,-121.697], Tmin=(1264.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.0786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(PolycyclicRing) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
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
    label = 'O=CCC[CH]OC1=CCC1(31317)',
    structure = SMILES('O=CCC[CH]OC1=CCC1'),
    E0 = (-25.3872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.13745,0.0978728,-6.69086e-05,1.16322e-08,3.69267e-12,-2855.32,37.2593], Tmin=(100,'K'), Tmax=(1042.17,'K')), NASAPolynomial(coeffs=[22.9941,0.032779,-1.28376e-05,2.38711e-09,-1.69084e-13,-9380.04,-87.3303], Tmin=(1042.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.3872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=CCC[CH]OC1C=CC1(31318)',
    structure = SMILES('O=CCC[CH]OC1C=CC1'),
    E0 = (4.95509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.961082,0.0956405,-7.40823e-05,2.88174e-08,-4.48283e-12,785.691,38.8382], Tmin=(100,'K'), Tmax=(1525,'K')), NASAPolynomial(coeffs=[22.3709,0.0344419,-1.3887e-05,2.50253e-09,-1.68921e-13,-6330.58,-83.5713], Tmin=(1525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.95509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CCC=COC1[CH]CC1(31319)',
    structure = SMILES('O=CCC=COC1[CH]CC1'),
    E0 = (-22.9193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.172,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.189073,0.0601379,5.91867e-05,-1.27018e-07,5.37688e-11,-2577.05,36.234], Tmin=(100,'K'), Tmax=(978.039,'K')), NASAPolynomial(coeffs=[27.0988,0.0268925,-9.99992e-06,2.05699e-09,-1.63229e-13,-11662.4,-113.967], Tmin=(978.039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.9193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCO)"""),
)

species(
    label = '[CH]1[CH]CC1(20477)',
    structure = SMILES('[CH]1[CH]CC1'),
    E0 = (389.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,180,388.387,1395.58,1395.59,1395.59,1395.59,1395.6,1395.6,1395.6,1395.61,1395.61,1802.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61381,-0.00133339,6.67733e-05,-7.44496e-08,2.51138e-11,46874.5,15.2546], Tmin=(100,'K'), Tmax=(1012.32,'K')), NASAPolynomial(coeffs=[3.9691,0.0238094,-9.81743e-06,1.89417e-09,-1.37335e-13,45442.3,6.81781], Tmin=(1012.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'CH2CHO(40)',
    structure = SMILES('[CH2]C=O'),
    E0 = (1.22925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,526.75,532.597,975.94,1639.13,1641.45],'cm^-1')),
        HinderedRotor(inertia=(0.00114821,'amu*angstrom^2'), symmetry=1, barrier=(2.1986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(1.22925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=COC1[CH]CC1(23711)',
    structure = SMILES('C=COC1[CH]CC1'),
    E0 = (112.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21728,0.0313461,8.90159e-05,-1.50615e-07,6.33681e-11,13642.5,24.8798], Tmin=(100,'K'), Tmax=(936.671,'K')), NASAPolynomial(coeffs=[23.7084,0.011405,-9.27447e-07,1.4578e-10,-2.29605e-14,6090.52,-99.9771], Tmin=(936.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCO)"""),
)

species(
    label = 'O=CCC[CH]O[C]1CCC1(31320)',
    structure = SMILES('O=CCC[CH]O[C]1CCC1'),
    E0 = (55.4007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.249829,0.0992375,-8.64567e-05,4.32714e-08,-9.41344e-12,6811.37,35.8793], Tmin=(100,'K'), Tmax=(1057.61,'K')), NASAPolynomial(coeffs=[11.0958,0.0563274,-2.55982e-05,4.90943e-09,-3.45451e-13,4411.5,-19.4924], Tmin=(1057.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.4007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCsJOCs) + radical(C2CsJOCs)"""),
)

species(
    label = 'O=CCC[CH]OC1C[CH]C1(31321)',
    structure = SMILES('O=CCC[CH]OC1C[CH]C1'),
    E0 = (62.5353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.244071,0.092394,-6.62533e-05,2.45393e-08,-3.79924e-12,7674.1,38.2386], Tmin=(100,'K'), Tmax=(1458.21,'K')), NASAPolynomial(coeffs=[15.6945,0.0486728,-2.12788e-05,3.97759e-09,-2.74049e-13,3025.79,-44.6677], Tmin=(1458.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.5353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCsJOCs) + radical(cyclobutane)"""),
)

species(
    label = 'O=CC[CH]COC1[CH]CC1(31322)',
    structure = SMILES('O=CC[CH]COC1[CH]CC1'),
    E0 = (94.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.326338,0.0812728,-3.47646e-05,-5.45256e-09,5.66695e-12,11478.3,41.7531], Tmin=(100,'K'), Tmax=(1148.26,'K')), NASAPolynomial(coeffs=[16.8761,0.0452853,-1.90237e-05,3.56384e-09,-2.4892e-13,5949.6,-50.4886], Tmin=(1148.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = 'O=CCCCO[C]1[CH]CC1(31323)',
    structure = SMILES('O=CCCCO[C]1[CH]CC1'),
    E0 = (74.8467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32667,0.0864229,-5.79456e-05,2.03577e-08,-3.10523e-12,9128.64,37.5242], Tmin=(100,'K'), Tmax=(1407.61,'K')), NASAPolynomial(coeffs=[11.2324,0.0554322,-2.49209e-05,4.71673e-09,-3.27306e-13,6058.43,-18.8183], Tmin=(1407.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(C2CsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'O=C[CH]CCOC1[CH]CC1(31324)',
    structure = SMILES('[O]C=CCCOC1[CH]CC1'),
    E0 = (37.3945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.477827,0.0749183,1.12514e-05,-7.20651e-08,3.36476e-11,4679.79,39.4988], Tmin=(100,'K'), Tmax=(980.579,'K')), NASAPolynomial(coeffs=[23.1455,0.0345346,-1.26079e-05,2.37663e-09,-1.74759e-13,-2644.49,-87.7301], Tmin=(980.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.3945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ) + radical(CCJCO)"""),
)

species(
    label = 'O=CCCCOC1[CH]C[CH]1(31325)',
    structure = SMILES('O=CCCCOC1[CH]C[CH]1'),
    E0 = (94.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.326338,0.0812728,-3.47646e-05,-5.45256e-09,5.66695e-12,11478.3,41.7531], Tmin=(100,'K'), Tmax=(1148.26,'K')), NASAPolynomial(coeffs=[16.8761,0.0452853,-1.90237e-05,3.56384e-09,-2.4892e-13,5949.6,-50.4886], Tmin=(1148.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = 'O=[C]CCCOC1[CH]CC1(31326)',
    structure = SMILES('O=[C]CCCOC1[CH]CC1'),
    E0 = (54.1035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.743815,0.090905,-6.21063e-05,2.11997e-08,-2.9159e-12,6689.08,41.5955], Tmin=(100,'K'), Tmax=(1691.18,'K')), NASAPolynomial(coeffs=[21.6621,0.03791,-1.5102e-05,2.67048e-09,-1.7679e-13,-889.402,-78.2729], Tmin=(1691.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.1035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=CC[CH][CH]OC1CCC1(31327)',
    structure = SMILES('O=CC[CH][CH]OC1CCC1'),
    E0 = (74.5989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,3000,3050,390,425,1340,1360,335,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05483,0.095959,-6.98693e-05,2.55588e-08,-3.74931e-12,9166.85,41.337], Tmin=(100,'K'), Tmax=(1606.83,'K')), NASAPolynomial(coeffs=[22.8074,0.0365569,-1.44166e-05,2.55163e-09,-1.69721e-13,1498.34,-85.1016], Tmin=(1606.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.5989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'O=CCCCOC1[CH][CH]C1(31328)',
    structure = SMILES('O=CCCCOC1[CH][CH]C1'),
    E0 = (81.9813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2782.5,750,1395,475,1775,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.140878,0.0843,-5.07751e-05,1.45593e-08,-1.67096e-12,10014.2,40.953], Tmin=(100,'K'), Tmax=(1964.37,'K')), NASAPolynomial(coeffs=[21.8459,0.0395284,-1.65871e-05,2.95645e-09,-1.94288e-13,1376.27,-79.9653], Tmin=(1964.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.9813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(CCJCO)"""),
)

species(
    label = 'O=C[CH]C[CH]OC1CCC1(31329)',
    structure = SMILES('[O]C=CC[CH]OC1CCC1'),
    E0 = (17.9485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.966302,0.087073,-1.63244e-05,-4.90788e-08,2.69821e-11,2357.37,38.2007], Tmin=(100,'K'), Tmax=(975.061,'K')), NASAPolynomial(coeffs=[25.2,0.0317198,-1.11489e-05,2.06487e-09,-1.51059e-13,-5216.79,-100.049], Tmin=(975.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.9485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ) + radical(CCsJOCs)"""),
)

species(
    label = 'O=[C]CC[CH]OC1CCC1(31330)',
    structure = SMILES('O=[C]CC[CH]OC1CCC1'),
    E0 = (34.6575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.893533,0.0994747,-7.89098e-05,3.24598e-08,-5.43208e-12,4351.02,39.0527], Tmin=(100,'K'), Tmax=(1406.19,'K')), NASAPolynomial(coeffs=[19.5217,0.0414019,-1.69625e-05,3.09079e-09,-2.10685e-13,-1390.5,-66.3983], Tmin=(1406.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.6575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCCJ=O) + radical(CCsJOCs)"""),
)

species(
    label = '[O][CH]CCC=O(5764)',
    structure = SMILES('[O][CH]CCC=O'),
    E0 = (19.0721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,180,180,2092.49],'cm^-1')),
        HinderedRotor(inertia=(0.166661,'amu*angstrom^2'), symmetry=1, barrier=(3.83187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164989,'amu*angstrom^2'), symmetry=1, barrier=(3.79342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166144,'amu*angstrom^2'), symmetry=1, barrier=(3.81999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56305,0.0641173,-0.000107371,1.07453e-07,-4.10375e-11,2371.43,23.494], Tmin=(100,'K'), Tmax=(845.992,'K')), NASAPolynomial(coeffs=[1.74033,0.0387928,-1.90537e-05,3.64333e-09,-2.50242e-13,3217.68,27.8472], Tmin=(845.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.0721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH]1CCC(O1)OC1[CH]CC1(31331)',
    structure = SMILES('[CH]1CCC(O1)OC1[CH]CC1'),
    E0 = (62.7367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249545,0.0543555,7.16903e-05,-1.38894e-07,5.98864e-11,7706.43,32.9586], Tmin=(100,'K'), Tmax=(930.237,'K')), NASAPolynomial(coeffs=[22.1413,0.0320377,-8.12566e-06,1.29905e-09,-9.42641e-14,526.231,-87.775], Tmin=(930.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.7367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Cyclobutane) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH]1CC[CH]OC2CCC2O1(31332)',
    structure = SMILES('[CH]1CC[CH]OC2CCC2O1'),
    E0 = (-72.1838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89375,0.129747,-0.000137459,7.58435e-08,-1.6793e-11,-8469.45,2.01572], Tmin=(100,'K'), Tmax=(1091.01,'K')), NASAPolynomial(coeffs=[21.5319,0.0438612,-1.93782e-05,3.69005e-09,-2.59479e-13,-13581,-113.04], Tmin=(1091.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.1838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + polycyclic(PolycyclicRing) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(C=C)O[CH]CCC=O(21400)',
    structure = SMILES('[CH2]C(C=C)O[CH]CCC=O'),
    E0 = (54.4759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4019.24,'J/mol'), sigma=(6.94615,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=627.80 K, Pc=27.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25782,0.12272,-0.000145179,9.67452e-08,-2.65659e-11,6735.11,41.1267], Tmin=(100,'K'), Tmax=(878.177,'K')), NASAPolynomial(coeffs=[14.3527,0.0516159,-2.3728e-05,4.54658e-09,-3.18864e-13,3993.34,-32.1569], Tmin=(878.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.4759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]CC=CO[CH]CCC=O(31333)',
    structure = SMILES('[CH2]CC=CO[CH]CCC=O'),
    E0 = (24.6358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8088,0.118486,-0.00011632,5.82789e-08,-1.15392e-11,3179.96,42.035], Tmin=(100,'K'), Tmax=(1227.61,'K')), NASAPolynomial(coeffs=[24.0206,0.0343252,-1.34852e-05,2.43351e-09,-1.6648e-13,-3161.77,-87.874], Tmin=(1227.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.6358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(RCCJ)"""),
)

species(
    label = 'O=CCCCOC1=CCC1(31334)',
    structure = SMILES('O=CCCCOC1=CCC1'),
    E0 = (-219.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.761927,0.0883621,-3.9545e-05,-1.17668e-08,1.05031e-11,-26191.9,36.439], Tmin=(100,'K'), Tmax=(1042.83,'K')), NASAPolynomial(coeffs=[20.577,0.039121,-1.56216e-05,2.92486e-09,-2.07477e-13,-32415.6,-75.9055], Tmin=(1042.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-219.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'O=CCC=COC1CCC1(31335)',
    structure = SMILES('O=CCC=COC1CCC1'),
    E0 = (-222.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383454,0.0653392,5.06632e-05,-1.17762e-07,5.00739e-11,-26613.7,34.4358], Tmin=(100,'K'), Tmax=(984.702,'K')), NASAPolynomial(coeffs=[26.6302,0.0310294,-1.19654e-05,2.42444e-09,-1.87969e-13,-35590.4,-114.04], Tmin=(984.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-222.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
)

species(
    label = 'O=CCCCOC1C=CC1(31336)',
    structure = SMILES('O=CCCCOC1C=CC1'),
    E0 = (-175.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612821,0.0880275,-5.26251e-05,1.18131e-08,-3.49332e-14,-20930.8,38.1452], Tmin=(100,'K'), Tmax=(1285.51,'K')), NASAPolynomial(coeffs=[18.9949,0.0426818,-1.79931e-05,3.33276e-09,-2.29296e-13,-27266.4,-66.4101], Tmin=(1285.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'CO(12)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(762.44,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35308e-07,1.51269e-10,-9.88872e-15,-14292.7,6.51157], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CC[CH]OC1[CH]CC1(20828)',
    structure = SMILES('CC[CH]OC1[CH]CC1'),
    E0 = (181.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.17,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3708.65,'J/mol'), sigma=(6.76005,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=579.28 K, Pc=27.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251153,0.0653452,4.10522e-07,-4.60496e-08,2.18024e-11,21940.8,32.3404], Tmin=(100,'K'), Tmax=(999.287,'K')), NASAPolynomial(coeffs=[17.5581,0.035322,-1.34457e-05,2.50458e-09,-1.79378e-13,16522,-60.9498], Tmin=(999.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = 'O=CCC[CH]O[CH]C1CC1(31337)',
    structure = SMILES('O=CCC[CH]O[CH]C1CC1'),
    E0 = (63.1736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17282,0.109672,-0.000101785,4.99196e-08,-9.9149e-12,7787.68,38.729], Tmin=(100,'K'), Tmax=(1204.91,'K')), NASAPolynomial(coeffs=[19.26,0.04184,-1.73398e-05,3.19638e-09,-2.20491e-13,2863.78,-63.656], Tmin=(1204.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.1736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C1CC1O[CH]CCC=O(21999)',
    structure = SMILES('[CH2]C1CC1O[CH]CCC=O'),
    E0 = (75.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01029,0.103973,-9.04519e-05,4.20246e-08,-7.94695e-12,9224.81,39.4432], Tmin=(100,'K'), Tmax=(1259.72,'K')), NASAPolynomial(coeffs=[18.4856,0.0420674,-1.6737e-05,3.01291e-09,-2.04736e-13,4313,-59.1143], Tmin=(1259.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C(CC=O)OC1[CH]CC1(31338)',
    structure = SMILES('[CH2]C(CC=O)OC1[CH]CC1'),
    E0 = (94.1021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729584,0.0953678,-6.99315e-05,2.60443e-08,-3.94781e-12,11494.9,40.109], Tmin=(100,'K'), Tmax=(1532.16,'K')), NASAPolynomial(coeffs=[19.8405,0.0416652,-1.73558e-05,3.16755e-09,-2.15032e-13,5191.62,-67.9065], Tmin=(1532.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.1021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = 'O=CCCC1OC2CCC12(31287)',
    structure = SMILES('O=CCCC1OC2CCC12'),
    E0 = (-148.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0295529,0.0744545,-1.81058e-05,-1.82675e-08,8.87828e-12,-17682.7,31.1938], Tmin=(100,'K'), Tmax=(1167.19,'K')), NASAPolynomial(coeffs=[15.4716,0.0491845,-2.14255e-05,4.07396e-09,-2.86248e-13,-23198.5,-54.1141], Tmin=(1167.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(519.654,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + polycyclic(s2_4_4_ane)"""),
)

species(
    label = 'C=COC[CH]OC1[CH]CC1(31339)',
    structure = SMILES('C=COC[CH]OC1[CH]CC1'),
    E0 = (116.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2418,0.0903708,-1.44532e-05,-6.01827e-08,3.33214e-11,14186.4,39.022], Tmin=(100,'K'), Tmax=(956.205,'K')), NASAPolynomial(coeffs=[28.7545,0.0254851,-7.72172e-06,1.3964e-09,-1.05419e-13,5679.72,-118.834], Tmin=(956.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'OC=CC[CH]OC1[CH]CC1(31340)',
    structure = SMILES('OC=CC[CH]OC1[CH]CC1'),
    E0 = (76.3879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.18,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17861,0.0887882,-1.10327e-05,-6.47485e-08,3.54928e-11,9396.69,40.0484], Tmin=(100,'K'), Tmax=(948.006,'K')), NASAPolynomial(coeffs=[28.8898,0.0239494,-6.5904e-06,1.14978e-09,-8.73077e-14,908.272,-118.11], Tmin=(948.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.3879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCsJOCs)"""),
)

species(
    label = '[CH]CCC=O(18484)',
    structure = SMILES('[CH]CCC=O'),
    E0 = (221.112,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,338.309,338.317,1520,1520.03],'cm^-1')),
        HinderedRotor(inertia=(0.00147232,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09193,'amu*angstrom^2'), symmetry=1, barrier=(7.47087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919401,'amu*angstrom^2'), symmetry=1, barrier=(7.47082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.39867,0.0383489,-2.66828e-05,1.05264e-08,-1.89984e-12,26648.5,18.3896], Tmin=(100,'K'), Tmax=(1182.31,'K')), NASAPolynomial(coeffs=[5.79639,0.0268536,-1.20987e-05,2.30283e-09,-1.60963e-13,25845.1,1.42861], Tmin=(1182.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C1[CH]CC1(23382)',
    structure = SMILES('[O]C1[CH]CC1'),
    E0 = (271.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88605,0.00621875,8.18905e-05,-1.10737e-07,4.26794e-11,32730.5,17.8324], Tmin=(100,'K'), Tmax=(952.303,'K')), NASAPolynomial(coeffs=[12.0158,0.0143731,-4.20109e-06,8.09184e-10,-6.53698e-14,28883.1,-36.8383], Tmin=(952.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CCJCO)"""),
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
    E0 = (74.5989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (180.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (120.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (197.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (222.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (195.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (115.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (149.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (228.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (233.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (202.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (232.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (202.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (140.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (171.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (99.3263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (126.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (172.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (202.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (408.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (132.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (119.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (178.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (147.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (137.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (82.9669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (82.9669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (419.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (234.534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (235.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (252.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (82.8832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (429.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (220.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (492.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['C1=CCC1(2460)', 'O=CCCC=O(5767)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['[O]C1CCC1OC1[CH]CC1(31316)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.13771e+06,'s^-1'), n=1.58803, Ea=(105.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 101.2 to 105.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['[O]C1CC[CH]OC2CCC12(31307)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.65621e+07,'s^-1'), n=0.823355, Ea=(45.7375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_csHNd] + [R8;multiplebond_intra;radadd_intra_cs] for rate rule [R8;carbonylbond_intra_H;radadd_intra_csHNd]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'O=CCC[CH]OC1=CCC1(31317)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'O=CCC[CH]OC1C=CC1(31318)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'O=CCC=COC1[CH]CC1(31319)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=CCCC=O(5767)', '[CH]1[CH]CC1(20477)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.6e+10,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2CHO(40)', 'C=COC1[CH]CC1(23711)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0114756,'m^3/(mol*s)'), n=2.44484, Ea=(36.1431,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CsJ-OneDeHH] for rate rule [Cds-HH_Cds-OsH;CsJ-COHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCC[CH]O[C]1CCC1(31320)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCC[CH]OC1C[CH]C1(31321)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=CC[CH]COC1[CH]CC1(31322)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCCCO[C]1[CH]CC1(31323)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.2544e+06,'s^-1'), n=1.86276, Ea=(157.734,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_O;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=C[CH]CCOC1[CH]CC1(31324)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=CCCCOC1[CH]C[CH]1(31325)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00456,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=[C]CCCOC1[CH]CC1(31326)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(27900,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeO;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeO;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=CC[CH][CH]OC1CCC1(31327)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.1128,'s^-1'), n=3.28, Ea=(24.7274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=CCCCOC1[CH][CH]C1(31328)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(264929,'s^-1'), n=1.581, Ea=(44.7981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_SSSS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=C[CH]C[CH]OC1CCC1(31329)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.86164e+06,'s^-1'), n=1.64667, Ea=(98.1149,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R6HJ_3;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=[C]CC[CH]OC1CCC1(31330)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(60440.4,'s^-1'), n=2.15184, Ea=(128.168,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeC;XH_out] for rate rule [R7HJ_3;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][CH]CCC=O(5764)', '[CH]1[CH]CC1(20477)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['[CH]1CCC(O1)OC1[CH]CC1(31331)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['[CH]1CC[CH]OC2CCC2O1(31332)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.95444e+07,'s^-1'), n=0.733383, Ea=(44.7109,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;multiplebond_intra;radadd_intra_csHCs] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)O[CH]CCC=O(21400)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CC=CO[CH]CCC=O(31333)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R4_Cs_HH_D;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCCCOC1=CCC1(31334)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCC=COC1CCC1(31335)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCCCOC1C=CC1(31336)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CO(12)', 'CC[CH]OC1[CH]CC1(20828)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCC[CH]O[CH]C1CC1(31337)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C1CC1O[CH]CCC=O(21999)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(CC=O)OC1[CH]CC1(31338)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    products = ['O=CCCC1OC2CCC12(31287)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=COC[CH]OC1[CH]CC1(31339)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['OC=CC[CH]OC1[CH]CC1(31340)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]CCC=O(18484)', '[O]C1[CH]CC1(23382)'],
    products = ['O=CCC[CH]OC1[CH]CC1(22017)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '4985',
    isomers = [
        'O=CCC[CH]OC1[CH]CC1(22017)',
    ],
    reactants = [
        ('C1=CCC1(2460)', 'O=CCCC=O(5767)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4985',
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

