species(
    label = 'C5H6O(219)(218)',
    structure = SMILES('[CH2]C=CC=C[O]'),
    E0 = (108.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.59991,'amu*angstrom^2'), symmetry=1, barrier=(36.785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59845,'amu*angstrom^2'), symmetry=1, barrier=(36.7515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3794.25,'J/mol'), sigma=(6.17562,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.65 K, Pc=36.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67303,0.0325845,3.83629e-05,-8.44772e-08,3.82244e-11,13165.2,20.9813], Tmin=(100,'K'), Tmax=(922.018,'K')), NASAPolynomial(coeffs=[18.3993,0.00812139,-9.24666e-08,-9.07048e-11,1.79206e-15,8036.21,-69.4428], Tmin=(922.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=COJ)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = 'C4H6(123)(122)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.210311,'amu*angstrom^2'), symmetry=1, barrier=(25.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779031,'amu*angstrom^2'), symmetry=1, barrier=(93.4717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56318,0.0223429,1.87067e-05,-3.93099e-08,1.63982e-11,33100.5,13.4097], Tmin=(100,'K'), Tmax=(974.264,'K')), NASAPolynomial(coeffs=[9.82995,0.0151966,-5.22272e-06,9.67656e-10,-7.0786e-14,30607.7,-26.985], Tmin=(974.264,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = 'C=[C]C=CC=O(4715)',
    structure = SMILES('C=C=CC=C[O]'),
    E0 = (167.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40235,'amu*angstrom^2'), symmetry=1, barrier=(32.2429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43985,0.0410476,6.30691e-06,-5.03167e-08,2.59189e-11,20212.2,19.3136], Tmin=(100,'K'), Tmax=(926.458,'K')), NASAPolynomial(coeffs=[18.7967,0.00545165,2.41126e-07,-1.15506e-10,3.7087e-15,15307.7,-72.2094], Tmin=(926.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=CC=C[C]=O(4549)',
    structure = SMILES('C=C[CH]C=C=O'),
    E0 = (103.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,292.831,292.831],'cm^-1')),
        HinderedRotor(inertia=(0.484811,'amu*angstrom^2'), symmetry=1, barrier=(29.5007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48481,'amu*angstrom^2'), symmetry=1, barrier=(29.5007,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78546,0.0413554,-1.91973e-05,-5.00058e-09,4.9111e-12,12521.9,18.5543], Tmin=(100,'K'), Tmax=(1025.01,'K')), NASAPolynomial(coeffs=[11.5734,0.0180947,-7.01485e-06,1.29194e-09,-9.09454e-14,9730.71,-32.7362], Tmin=(1025.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(C=CCJC=C=O)"""),
)

species(
    label = 'C=C[C]=CC=O(4714)',
    structure = SMILES('C=C[C]=CC=O'),
    E0 = (169.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.18163,'amu*angstrom^2'), symmetry=1, barrier=(27.168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17982,'amu*angstrom^2'), symmetry=1, barrier=(27.1265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83685,0.0478006,-4.01103e-05,1.75891e-08,-3.19807e-12,20455.1,18.4329], Tmin=(100,'K'), Tmax=(1277.85,'K')), NASAPolynomial(coeffs=[10.001,0.0222446,-1.01115e-05,1.93838e-09,-1.36142e-13,18368.6,-22.9561], Tmin=(1277.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=[C]C=O(4716)',
    structure = SMILES('C=CC=C=C[O]'),
    E0 = (167.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40235,'amu*angstrom^2'), symmetry=1, barrier=(32.2429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43985,0.0410476,6.30691e-06,-5.03167e-08,2.59189e-11,20212.2,19.3136], Tmin=(100,'K'), Tmax=(926.458,'K')), NASAPolynomial(coeffs=[18.7967,0.00545165,2.41126e-07,-1.15506e-10,3.7087e-15,15307.7,-72.2094], Tmin=(926.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CHCHO(45)(45)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (150.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3215,0.0099565,8.23489e-06,-1.92986e-08,8.37464e-12,18150.6,9.47663], Tmin=(100,'K'), Tmax=(969.366,'K')), NASAPolynomial(coeffs=[7.73867,0.0039342,-1.33182e-06,2.69015e-10,-2.15546e-14,16720.8,-14.6539], Tmin=(969.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C[CH2](924)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,229.711,230.18,230.787],'cm^-1')),
        HinderedRotor(inertia=(1.33306,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,44095,-3.44606], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C=CC=CO(4279)',
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
    label = 'C5H6O(222)(221)',
    structure = SMILES('CC=CC=C=O'),
    E0 = (-7.61157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.540838,'amu*angstrom^2'), symmetry=1, barrier=(12.4349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.542177,'amu*angstrom^2'), symmetry=1, barrier=(12.4657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3093,'J/mol'), sigma=(5.67324,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=483.12 K, Pc=38.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25024,0.0422182,-2.96346e-05,1.21411e-08,-2.32488e-12,-855.642,18.4018], Tmin=(100,'K'), Tmax=(1109.93,'K')), NASAPolynomial(coeffs=[5.49731,0.0305163,-1.382e-05,2.64226e-09,-1.85352e-13,-1576.44,2.39787], Tmin=(1109.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.61157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'CC=C=CC=O(4730)',
    structure = SMILES('CC=C=CC=O'),
    E0 = (23.2144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4813,0.0418731,-2.50262e-05,6.67668e-09,-7.13521e-13,2837.61,16.5114], Tmin=(100,'K'), Tmax=(1950.41,'K')), NASAPolynomial(coeffs=[10.4372,0.0255568,-1.24779e-05,2.38761e-09,-1.63759e-13,-265.864,-27.186], Tmin=(1950.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.2144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CCC=O(4731)',
    structure = SMILES('C=C=CCC=O'),
    E0 = (47.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9532,0.0363655,-2.7601e-06,-1.67698e-08,6.90597e-12,5786.48,19.8864], Tmin=(100,'K'), Tmax=(1171.94,'K')), NASAPolynomial(coeffs=[11.261,0.0239241,-1.15738e-05,2.31612e-09,-1.67335e-13,2277.56,-32.158], Tmin=(1171.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC=C=CO(4293)',
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
    label = 'C=CCC=C=O(4732)',
    structure = SMILES('C=CCC=C=O'),
    E0 = (2.26599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91283,0.0462809,-3.56684e-05,1.505e-08,-2.70323e-12,347.216,19.5763], Tmin=(100,'K'), Tmax=(1270.75,'K')), NASAPolynomial(coeffs=[8.598,0.0252377,-1.0829e-05,2.01866e-09,-1.3953e-13,-1351.82,-14.2776], Tmin=(1270.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.26599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C5H6O(231)(230)',
    structure = SMILES('C1C=COCC=1'),
    E0 = (-51.2755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3681.99,'J/mol'), sigma=(6.05729,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=575.12 K, Pc=37.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19667,0.0134732,9.65046e-05,-1.46039e-07,6.00172e-11,-6077.91,13.1812], Tmin=(100,'K'), Tmax=(930.415,'K')), NASAPolynomial(coeffs=[19.7128,0.0053276,1.36432e-06,-2.88683e-10,8.84395e-15,-12244.2,-85.6821], Tmin=(930.415,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.2755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(1,3-Cyclohexadiene)"""),
)

species(
    label = 'C5H6O(218)(217)',
    structure = SMILES('O=CC1C=CC1'),
    E0 = (30.9926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.88,'J/mol'), sigma=(5.9564,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.82 K, Pc=39.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34792,0.0256198,2.5475e-05,-4.66541e-08,1.81982e-11,3796.43,18.009], Tmin=(100,'K'), Tmax=(1008.25,'K')), NASAPolynomial(coeffs=[9.5562,0.0229525,-9.13368e-06,1.73733e-09,-1.25454e-13,1024.91,-23.3619], Tmin=(1008.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.9926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'C5H6O(220)(219)',
    structure = SMILES('C=CC1C=CO1'),
    E0 = (62.7482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3476.3,'J/mol'), sigma=(5.7538,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=542.99 K, Pc=41.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90702,0.0189559,8.69326e-05,-1.40405e-07,5.88222e-11,7647.19,17.8452], Tmin=(100,'K'), Tmax=(933.337,'K')), NASAPolynomial(coeffs=[21.9812,0.00224299,2.38776e-06,-4.41596e-10,1.7573e-14,880.759,-93.7905], Tmin=(933.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.7482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[O][CH]C1C=CC1(4733)',
    structure = SMILES('[O][CH]C1C=CC1'),
    E0 = (358.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06048,0.0350304,-2.44794e-06,-1.67858e-08,7.45533e-12,43240.9,20.6853], Tmin=(100,'K'), Tmax=(1075.63,'K')), NASAPolynomial(coeffs=[9.30928,0.024255,-9.98643e-06,1.87205e-09,-1.31702e-13,40745.4,-19.1657], Tmin=(1075.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C1C=CO1(4734)',
    structure = SMILES('[CH2][CH]C1C=CO1'),
    E0 = (343.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83042,0.0199859,8.64565e-05,-1.43068e-07,6.07842e-11,41461.9,22.0919], Tmin=(100,'K'), Tmax=(927.24,'K')), NASAPolynomial(coeffs=[23.0959,-0.000391489,3.98272e-06,-7.73085e-10,4.1412e-14,34450.6,-95.4373], Tmin=(927.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = 'C[C]=CC=C[O](4735)',
    structure = SMILES('C[C]=CC=C[O]'),
    E0 = (228.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36036,0.0461037,-9.87951e-06,-2.87128e-08,1.68452e-11,27577.1,21.0184], Tmin=(100,'K'), Tmax=(937.937,'K')), NASAPolynomial(coeffs=[16.4745,0.0116893,-2.88753e-06,4.67116e-10,-3.47805e-14,23420.4,-57.9745], Tmin=(937.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CC=[C]O(4287)',
    structure = SMILES('[CH2]C=CC=[C]O'),
    E0 = (206.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.17639,'amu*angstrom^2'), symmetry=1, barrier=(27.0474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17538,'amu*angstrom^2'), symmetry=1, barrier=(27.0242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17862,'amu*angstrom^2'), symmetry=1, barrier=(27.0988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47241,0.04129,8.21339e-06,-5.21908e-08,2.69244e-11,24989,23.4205], Tmin=(100,'K'), Tmax=(909.844,'K')), NASAPolynomial(coeffs=[17.5834,0.00844141,-2.48295e-07,-1.09576e-10,6.9686e-15,20485.3,-61.4227], Tmin=(909.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CC=[C]C=C[O](4736)',
    structure = SMILES('CC=[C]C=C[O]'),
    E0 = (189.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30837,'amu*angstrom^2'), symmetry=1, barrier=(30.082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30876,'amu*angstrom^2'), symmetry=1, barrier=(30.0909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31692,0.0480145,-1.39882e-05,-2.55089e-08,1.62953e-11,22905.7,20.2981], Tmin=(100,'K'), Tmax=(918.987,'K')), NASAPolynomial(coeffs=[16.1654,0.0122081,-2.59012e-06,3.51538e-10,-2.4112e-14,18959.5,-56.7046], Tmin=(918.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=C[C]=CO(4285)',
    structure = SMILES('[CH2]C=C[C]=CO'),
    E0 = (166.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.59361,'amu*angstrom^2'), symmetry=1, barrier=(36.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58873,'amu*angstrom^2'), symmetry=1, barrier=(36.5281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58986,'amu*angstrom^2'), symmetry=1, barrier=(36.5539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18563,0.0460273,6.37531e-06,-5.89599e-08,3.19557e-11,20100.1,20.92], Tmin=(100,'K'), Tmax=(890.779,'K')), NASAPolynomial(coeffs=[20.2965,0.00421406,2.68736e-06,-7.44084e-10,5.3374e-14,14949.6,-78.8679], Tmin=(890.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C[CH]C=[C]C=O(4722)',
    structure = SMILES('CC=C[C]=C[O]'),
    E0 = (189.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30837,'amu*angstrom^2'), symmetry=1, barrier=(30.082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30876,'amu*angstrom^2'), symmetry=1, barrier=(30.0909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31692,0.0480145,-1.39882e-05,-2.55089e-08,1.62953e-11,22905.7,20.2981], Tmin=(100,'K'), Tmax=(918.987,'K')), NASAPolynomial(coeffs=[16.1654,0.0122081,-2.59012e-06,3.51538e-10,-2.4112e-14,18959.5,-56.7046], Tmin=(918.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=[C]C=CO(4284)',
    structure = SMILES('[CH2]C=[C]C=CO'),
    E0 = (166.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.59361,'amu*angstrom^2'), symmetry=1, barrier=(36.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58873,'amu*angstrom^2'), symmetry=1, barrier=(36.5281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58986,'amu*angstrom^2'), symmetry=1, barrier=(36.5539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18563,0.0460273,6.37531e-06,-5.89599e-08,3.19557e-11,20100.1,20.92], Tmin=(100,'K'), Tmax=(890.779,'K')), NASAPolynomial(coeffs=[20.2965,0.00421406,2.68736e-06,-7.44084e-10,5.3374e-14,14949.6,-78.8679], Tmin=(890.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C[CH]C=C[C]=O(4725)',
    structure = SMILES('C[CH]C=C[C]=O'),
    E0 = (161.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,206.708],'cm^-1')),
        HinderedRotor(inertia=(0.561588,'amu*angstrom^2'), symmetry=1, barrier=(16.7763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566559,'amu*angstrom^2'), symmetry=1, barrier=(16.7671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581782,'amu*angstrom^2'), symmetry=1, barrier=(16.753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99911,0.0478572,-3.55846e-05,1.33834e-08,-2.13732e-12,19519.3,18.7486], Tmin=(100,'K'), Tmax=(1380,'K')), NASAPolynomial(coeffs=[9.08188,0.0273272,-1.32691e-05,2.60279e-09,-1.8429e-13,17564.5,-17.7028], Tmin=(1380,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_S) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2][C]=CC=CO(4286)',
    structure = SMILES('[CH2][C]=CC=CO'),
    E0 = (205.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.29765,'amu*angstrom^2'), symmetry=1, barrier=(29.8356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30302,'amu*angstrom^2'), symmetry=1, barrier=(29.9589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29819,'amu*angstrom^2'), symmetry=1, barrier=(29.848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23068,0.0441116,1.042e-05,-6.19274e-08,3.23181e-11,24771.3,21.6337], Tmin=(100,'K'), Tmax=(903.577,'K')), NASAPolynomial(coeffs=[20.5447,0.0037976,2.33142e-06,-6.1477e-10,4.15718e-14,19436.4,-79.7943], Tmin=(903.577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=CC1[CH]C1(4737)',
    structure = SMILES('[O]C=CC1[CH]C1'),
    E0 = (272.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15684,0.0207567,6.23783e-05,-1.02283e-07,4.24349e-11,32907.4,21.863], Tmin=(100,'K'), Tmax=(939.562,'K')), NASAPolynomial(coeffs=[16.3185,0.0107621,-1.96197e-06,3.44391e-10,-3.21619e-14,28026.3,-57.3899], Tmin=(939.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=COJ) + radical(cyclopropane)"""),
)

species(
    label = '[CH2]C=CC1[CH]O1(4738)',
    structure = SMILES('C=C[CH]C1[CH]O1'),
    E0 = (277.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38484,0.0445863,5.20922e-06,-5.30939e-08,2.91953e-11,33473.7,18.5154], Tmin=(100,'K'), Tmax=(873.225,'K')), NASAPolynomial(coeffs=[17.1888,0.00970641,6.84957e-07,-4.42987e-10,3.67884e-14,29283.4,-63.7767], Tmin=(873.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CCJCO) + radical(CCsJO)"""),
)

species(
    label = '[O]C1[CH]C=CC1(4739)',
    structure = SMILES('[O]C1[CH]C=CC1'),
    E0 = (196.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16702,0.0208235,6.33905e-05,-9.95311e-08,3.98738e-11,23753.7,15.9164], Tmin=(100,'K'), Tmax=(960.854,'K')), NASAPolynomial(coeffs=[15.1483,0.0157363,-5.0897e-06,1.00586e-09,-8.04007e-14,18999.3,-57.9518], Tmin=(960.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1[CH]C=CO1(4740)',
    structure = SMILES('[CH2]C1[CH]C=CO1'),
    E0 = (162.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55136,0.0322472,5.03143e-05,-1.06749e-07,4.9413e-11,19650.3,13.2885], Tmin=(100,'K'), Tmax=(898.403,'K')), NASAPolynomial(coeffs=[21.3454,0.00192799,4.41388e-06,-1.06308e-09,7.24434e-14,13760.7,-93.0698], Tmin=(898.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CJC(C)OC) + radical(C=CCJCO)"""),
)

species(
    label = 'CH2(17)(18)',
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
    label = '[CH]C=CC=O(4637)',
    structure = SMILES('[CH]C=CC=O'),
    E0 = (253.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13955,'amu*angstrom^2'), symmetry=1, barrier=(49.1924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13866,'amu*angstrom^2'), symmetry=1, barrier=(49.172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45042,0.0332273,-1.49834e-05,1.57226e-09,2.41248e-13,30523.8,16.9403], Tmin=(100,'K'), Tmax=(1750.2,'K')), NASAPolynomial(coeffs=[11.9107,0.0188928,-8.94295e-06,1.65009e-09,-1.09645e-13,26096.4,-37.1832], Tmin=(1750.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=CC=C[CH2](922)',
    structure = SMILES('[CH]=CC=C[CH2]'),
    E0 = (423.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.52606,'amu*angstrom^2'), symmetry=1, barrier=(35.0872,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53677,'amu*angstrom^2'), symmetry=1, barrier=(35.3333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22719,0.0272982,2.2821e-05,-5.20811e-08,2.2993e-11,50955.4,18.1254], Tmin=(100,'K'), Tmax=(937.462,'K')), NASAPolynomial(coeffs=[12.1491,0.0145308,-4.06041e-06,6.79572e-10,-4.92004e-14,47795.9,-36.031], Tmin=(937.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=[C]CC=O(4741)',
    structure = SMILES('[CH2]C=[C]CC=O'),
    E0 = (260.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1226.02],'cm^-1')),
        HinderedRotor(inertia=(0.0379246,'amu*angstrom^2'), symmetry=1, barrier=(20.3083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8822,'amu*angstrom^2'), symmetry=1, barrier=(20.2835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.883578,'amu*angstrom^2'), symmetry=1, barrier=(20.3152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87594,0.0388417,-9.80423e-06,-9.75618e-09,4.58851e-12,31374.9,21.867], Tmin=(100,'K'), Tmax=(1218.79,'K')), NASAPolynomial(coeffs=[11.5565,0.0236147,-1.14252e-05,2.26795e-09,-1.62452e-13,27786.4,-31.7922], Tmin=(1218.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C[CH]C[C]=O(4719)',
    structure = SMILES('[CH2]C=CC[C]=O'),
    E0 = (182.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,315.059],'cm^-1')),
        HinderedRotor(inertia=(0.00169838,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312441,'amu*angstrom^2'), symmetry=1, barrier=(22.0073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312443,'amu*angstrom^2'), symmetry=1, barrier=(22.0073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84068,0.0373102,-1.18216e-06,-2.21776e-08,9.72894e-12,22011.4,22.3283], Tmin=(100,'K'), Tmax=(1088.78,'K')), NASAPolynomial(coeffs=[12.2928,0.0212541,-9.84398e-06,1.97424e-09,-1.44531e-13,18411.1,-35.0678], Tmin=(1088.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=CCC=O(4742)',
    structure = SMILES('[CH2][C]=CCC=O'),
    E0 = (260.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1226.02],'cm^-1')),
        HinderedRotor(inertia=(0.0379246,'amu*angstrom^2'), symmetry=1, barrier=(20.3083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8822,'amu*angstrom^2'), symmetry=1, barrier=(20.2835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.883578,'amu*angstrom^2'), symmetry=1, barrier=(20.3152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87594,0.0388417,-9.80423e-06,-9.75618e-09,4.58851e-12,31374.9,21.867], Tmin=(100,'K'), Tmax=(1218.79,'K')), NASAPolynomial(coeffs=[11.5565,0.0236147,-1.14252e-05,2.26795e-09,-1.62452e-13,27786.4,-31.7922], Tmin=(1218.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]C1C=O(4743)',
    structure = SMILES('[CH2]C1[CH]C1C=O'),
    E0 = (301.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26061,0.0255862,3.26211e-05,-6.21187e-08,2.63556e-11,36383.8,22.8444], Tmin=(100,'K'), Tmax=(939.291,'K')), NASAPolynomial(coeffs=[11.7735,0.017094,-4.9499e-06,8.39146e-10,-6.04672e-14,33184.3,-29.9729], Tmin=(939.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCJCC=O) + radical(Isobutyl)"""),
)

species(
    label = '[CH]1[CH]OCC=C1(4744)',
    structure = SMILES('[CH]1C=C[CH]OC1'),
    E0 = (107.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40456,0.0204729,5.43974e-05,-8.33941e-08,3.32093e-11,13013.2,13.3209], Tmin=(100,'K'), Tmax=(945.232,'K')), NASAPolynomial(coeffs=[10.6722,0.0223543,-7.09444e-06,1.23988e-09,-8.89261e-14,9803.27,-34.8122], Tmin=(945.232,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(C=CCJCO) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=[C]C[CH]C=O(4718)',
    structure = SMILES('C=[C]CC=C[O]'),
    E0 = (259.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,193.632,193.635,193.659],'cm^-1')),
        HinderedRotor(inertia=(0.750811,'amu*angstrom^2'), symmetry=1, barrier=(19.9639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750884,'amu*angstrom^2'), symmetry=1, barrier=(19.9639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58152,0.0414199,-1.96358e-06,-3.23802e-08,1.68456e-11,31279.7,23.247], Tmin=(100,'K'), Tmax=(955.649,'K')), NASAPolynomial(coeffs=[14.9255,0.0142834,-4.44361e-06,7.93615e-10,-5.85182e-14,27418,-47.3852], Tmin=(955.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C=CC[C]=C[O](4745)',
    structure = SMILES('C=CC[C]=C[O]'),
    E0 = (259.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,193.632,193.635,193.659],'cm^-1')),
        HinderedRotor(inertia=(0.750811,'amu*angstrom^2'), symmetry=1, barrier=(19.9639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750884,'amu*angstrom^2'), symmetry=1, barrier=(19.9639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58152,0.0414199,-1.96358e-06,-3.23802e-08,1.68456e-11,31279.7,23.247], Tmin=(100,'K'), Tmax=(955.649,'K')), NASAPolynomial(coeffs=[14.9255,0.0142834,-4.44361e-06,7.93615e-10,-5.85182e-14,27418,-47.3852], Tmin=(955.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC[CH]C=O(4721)',
    structure = SMILES('[CH]=CCC=C[O]'),
    E0 = (268.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.96306,'amu*angstrom^2'), symmetry=1, barrier=(22.1426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966461,'amu*angstrom^2'), symmetry=1, barrier=(22.2208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54104,0.0403358,6.2145e-06,-4.43749e-08,2.19672e-11,32396.1,23.322], Tmin=(100,'K'), Tmax=(945.635,'K')), NASAPolynomial(coeffs=[16.3545,0.0119522,-3.13312e-06,5.46156e-10,-4.2351e-14,28061.9,-55.4194], Tmin=(945.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CCC=[C][O](4746)',
    structure = SMILES('C=CC[CH][C]=O'),
    E0 = (204.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,212.24,2238.04],'cm^-1')),
        HinderedRotor(inertia=(0.00370789,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506147,'amu*angstrom^2'), symmetry=1, barrier=(16.4674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810513,'amu*angstrom^2'), symmetry=1, barrier=(26.4327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70712,0.0447707,-2.9728e-05,8.35892e-09,-2.88686e-13,24669.6,23.1654], Tmin=(100,'K'), Tmax=(1111.25,'K')), NASAPolynomial(coeffs=[10.3514,0.0211694,-8.01303e-06,1.41651e-09,-9.57856e-14,22284.4,-21.5375], Tmin=(1111.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C=CC=CO(4292)',
    structure = SMILES('[CH]C=CC=CO'),
    E0 = (219.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.96466,'amu*angstrom^2'), symmetry=1, barrier=(45.1714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96743,'amu*angstrom^2'), symmetry=1, barrier=(45.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95958,'amu*angstrom^2'), symmetry=1, barrier=(45.0546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.987439,0.0483026,1.24103e-05,-6.58071e-08,3.34542e-11,26559.5,21.5833], Tmin=(100,'K'), Tmax=(913.7,'K')), NASAPolynomial(coeffs=[20.8413,0.00862734,-1.03498e-08,-1.58151e-10,9.01953e-15,20959.5,-83.1994], Tmin=(913.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC1[CH]C1[O](4747)',
    structure = SMILES('C=CC1[CH]C1[O]'),
    E0 = (372.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23175,0.019923,6.11621e-05,-9.86406e-08,4.05252e-11,44836.1,21.9469], Tmin=(100,'K'), Tmax=(943.727,'K')), NASAPolynomial(coeffs=[15.4391,0.0120362,-2.74289e-06,5.0206e-10,-4.31025e-14,40201.7,-52.3524], Tmin=(943.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1C=CC1[O](4748)',
    structure = SMILES('[CH2]C1C=CC1[O]'),
    E0 = (375.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20752,0.0212675,5.7579e-05,-9.53896e-08,3.97479e-11,45225.6,21.204], Tmin=(100,'K'), Tmax=(934.84,'K')), NASAPolynomial(coeffs=[15.0981,0.0124013,-2.46938e-06,4.00886e-10,-3.38393e-14,40792.8,-50.9353], Tmin=(934.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Isobutyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C[C]=CC=O(4717)',
    structure = SMILES('[CH2]C[C]=CC=O'),
    E0 = (303.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,390.712],'cm^-1')),
        HinderedRotor(inertia=(0.0824497,'amu*angstrom^2'), symmetry=1, barrier=(8.79403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0817205,'amu*angstrom^2'), symmetry=1, barrier=(8.8011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0818596,'amu*angstrom^2'), symmetry=1, barrier=(8.78964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1176,0.0468879,-4.07166e-05,2.28379e-08,-6.27065e-12,36514.9,22.196], Tmin=(100,'K'), Tmax=(802.944,'K')), NASAPolynomial(coeffs=[4.61844,0.0344296,-1.74432e-05,3.51459e-09,-2.54306e-13,36113.3,10.6798], Tmin=(802.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC=[C]C=O(4720)',
    structure = SMILES('[CH2]CC=C=C[O]'),
    E0 = (261.973,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.980754,'amu*angstrom^2'), symmetry=1, barrier=(22.5495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978519,'amu*angstrom^2'), symmetry=1, barrier=(22.4981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44143,0.0451847,-1.09227e-05,-2.37769e-08,1.39006e-11,31610.3,23.1768], Tmin=(100,'K'), Tmax=(959.347,'K')), NASAPolynomial(coeffs=[15.1908,0.0145156,-4.6525e-06,8.31928e-10,-6.07334e-14,27745.4,-48.979], Tmin=(959.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]CC=C[C]=O(4724)',
    structure = SMILES('[CH2]CC=C[C]=O'),
    E0 = (225.858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.570437,'amu*angstrom^2'), symmetry=1, barrier=(13.1155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569851,'amu*angstrom^2'), symmetry=1, barrier=(13.102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569841,'amu*angstrom^2'), symmetry=1, barrier=(13.1018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7789,0.0541679,-5.48797e-05,3.25938e-08,-8.48714e-12,27239.9,22.4585], Tmin=(100,'K'), Tmax=(892.974,'K')), NASAPolynomial(coeffs=[7.03472,0.0306247,-1.53317e-05,3.06817e-09,-2.2096e-13,26301.3,-2.30275], Tmin=(892.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(RCCJ)"""),
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
    label = '[CH]=C[CH]CC=O(4723)',
    structure = SMILES('[CH]C=CCC=O'),
    E0 = (241.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,399.151,399.151,399.152,399.154],'cm^-1')),
        HinderedRotor(inertia=(0.433696,'amu*angstrom^2'), symmetry=1, barrier=(49.032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433689,'amu*angstrom^2'), symmetry=1, barrier=(49.0321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.433717,'amu*angstrom^2'), symmetry=1, barrier=(49.0322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90949,0.0366928,6.54401e-06,-2.56812e-08,9.48674e-12,29131,22.0986], Tmin=(100,'K'), Tmax=(1159.72,'K')), NASAPolynomial(coeffs=[10.0066,0.0308595,-1.44886e-05,2.83732e-09,-2.02355e-13,25767.1,-24.5711], Tmin=(1159.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=CC=[C]C[O](2834)',
    structure = SMILES('C=CC=[C]C[O]'),
    E0 = (368.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,283.498,283.498,283.498],'cm^-1')),
        HinderedRotor(inertia=(0.113439,'amu*angstrom^2'), symmetry=1, barrier=(6.46974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113438,'amu*angstrom^2'), symmetry=1, barrier=(6.46974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08753,0.0455227,-3.82959e-05,2.05796e-08,-5.13328e-12,44363,21.7679], Tmin=(100,'K'), Tmax=(904.956,'K')), NASAPolynomial(coeffs=[5.40354,0.0308652,-1.39998e-05,2.68053e-09,-1.88414e-13,43762.8,6.10133], Tmin=(904.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'C=C[C]=CC[O](2835)',
    structure = SMILES('C=C[C]=CC[O]'),
    E0 = (329.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.156902,'amu*angstrom^2'), symmetry=1, barrier=(3.6075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.08362,'amu*angstrom^2'), symmetry=1, barrier=(93.8904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15325,0.0453843,-3.03584e-05,-2.52903e-09,1.30214e-11,39686.9,20.6968], Tmin=(100,'K'), Tmax=(566.82,'K')), NASAPolynomial(coeffs=[5.31978,0.0310279,-1.35096e-05,2.52177e-09,-1.74319e-13,39199.6,6.08568], Tmin=(566.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CCOJ)"""),
)

species(
    label = 'C=[C]C=CC[O](2836)',
    structure = SMILES('C=C=C[CH]C[O]'),
    E0 = (306.851,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,540,610,2055,3010,987.5,1337.5,450,1655,297.103,297.104,297.105],'cm^-1')),
        HinderedRotor(inertia=(0.785507,'amu*angstrom^2'), symmetry=1, barrier=(49.2031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615863,'amu*angstrom^2'), symmetry=1, barrier=(38.5768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87266,0.0513611,-4.52783e-05,2.42771e-08,-5.92264e-12,36978.4,19.3534], Tmin=(100,'K'), Tmax=(929.364,'K')), NASAPolynomial(coeffs=[6.12819,0.0330454,-1.57167e-05,3.07157e-09,-2.18343e-13,36187.4,-0.86533], Tmin=(929.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=CC=CC[O](2837)',
    structure = SMILES('[CH]=CC=CC[O]'),
    E0 = (377.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,280.173,281.892],'cm^-1')),
        HinderedRotor(inertia=(0.175586,'amu*angstrom^2'), symmetry=1, barrier=(9.8612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17786,'amu*angstrom^2'), symmetry=1, barrier=(9.86969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00808,0.0450931,-3.31887e-05,1.34218e-08,-2.34876e-12,45480.6,21.9681], Tmin=(100,'K'), Tmax=(1279.23,'K')), NASAPolynomial(coeffs=[7.89491,0.0266858,-1.16047e-05,2.17342e-09,-1.50488e-13,43974.4,-7.88216], Tmin=(1279.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ)"""),
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
    E0 = (390.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (347.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (398.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (395.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (527.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (133.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (133.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (186.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (197.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (186.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (197.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (117.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (116.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (116.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (358.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (343.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (316.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (369.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (382.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (359.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (293.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (269.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (194.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (238.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (334.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (334.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (252.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (224.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (634.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (666.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (394.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (294.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (407.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (339.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (196.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (416.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (416.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (413.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (351.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (371.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (372.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (375.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (465.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (438.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (364.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (460.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (108.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (285.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (531.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (522.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (410.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (410.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['H(3)(3)', 'C=[C]C=CC=O(4715)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)(3)', 'C=CC=C[C]=O(4549)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', 'C=C[C]=CC=O(4714)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', 'C=CC=[C]C=O(4716)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCHO(45)(45)', '[CH]=C[CH2](924)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(556926,'m^3/(mol*s)'), n=0.4, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;Cd_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -2.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C5H6O(219)(218)'],
    products = ['C=C=CC=CO(4279)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C5H6O(219)(218)'],
    products = ['C5H6O(222)(221)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C5H6O(219)(218)'],
    products = ['CC=C=CC=O(4730)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C5H6O(219)(218)'],
    products = ['C=C=CCC=O(4731)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C5H6O(219)(218)'],
    products = ['C=CC=C=CO(4293)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C5H6O(219)(218)'],
    products = ['C=CCC=C=O(4732)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C5H6O(219)(218)'],
    products = ['C5H6O(231)(230)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SDSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C5H6O(219)(218)'],
    products = ['C5H6O(218)(217)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_H/OneDe;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C5H6O(219)(218)'],
    products = ['C5H6O(220)(219)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C5H6O(219)(218)'],
    products = ['[O][CH]C1C=CC1(4733)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(250.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 249.9 to 250.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C5H6O(219)(218)'],
    products = ['[CH2][CH]C1C=CO1(4734)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(235.25,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 233.7 to 235.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C5H6O(219)(218)'],
    products = ['C[C]=CC=C[O](4735)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C=CC=[C]O(4287)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC=[C]C=C[O](4736)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 288 used for R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=C[C]=CO(4285)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]C=[C]C=O(4722)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.33e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 284 used for R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C=[C]C=CO(4284)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SDS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]C=C[C]=O(4725)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R5H_DSMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C]=CC=CO(4286)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out;XH_out] for rate rule [R5H_DSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C5H6O(219)(218)'],
    products = ['[O]C=CC1[CH]C1(4737)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""From training reaction 142 used for R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C5H6O(219)(218)'],
    products = ['[CH2]C=CC1[CH]O1(4738)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri_HCd;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HCd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C5H6O(219)(218)'],
    products = ['[O]C1[CH]C=CC1(4739)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.0909e+10,'s^-1'), n=0.248798, Ea=(144.352,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SD_D;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C5H6O(219)(218)'],
    products = ['[CH2]C1[CH]C=CO1(4740)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.56048e+11,'s^-1'), n=0.0862073, Ea=(116.147,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(17)(18)', '[CH]C=CC=O(4637)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(4)(4)', '[CH]=CC=C[CH2](922)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C=[C]CC=O(4741)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.82494e+10,'s^-1'), n=0.9, Ea=(134.725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/OneDe] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C[CH]C[C]=O(4719)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C]=CCC=O(4742)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/CO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C5H6O(219)(218)'],
    products = ['[CH2]C1[CH]C1C=O(4743)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C5H6O(219)(218)'],
    products = ['[CH]1[CH]OCC=C1(4744)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]C[CH]C=O(4718)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=CC[C]=C[O](4745)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.448e+10,'s^-1'), n=0.82, Ea=(156.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 154 used for R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=CC[CH]C=O(4721)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=CCC=[C][O](4746)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]C=CC=CO(4292)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C5H6O(219)(218)'],
    products = ['C=CC1[CH]C1[O](4747)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(263.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_csHDe] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic
Ea raised from 261.5 to 263.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C5H6O(219)(218)'],
    products = ['[CH2]C1C=CC1[O](4748)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(266.73,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 264.9 to 266.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C[C]=CC=O(4717)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC=[C]C=O(4720)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.65652e+07,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]CC=C[C]=O(4724)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.61619e+08,'s^-1'), n=1.1965, Ea=(138.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;XH_out] for rate rule [R4H_SDS;CO_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CO(10)(11)', 'C4H6(123)(122)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.236794,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C5H6O(219)(218)'],
    products = ['C5H6O(217)(216)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.01e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for Y_12_01
Exact match found for rate rule [Y_12_01]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C[CH]CC=O(4723)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R4HJ_2;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=CC=[C]C[O](2834)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C[C]=CC[O](2835)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=[C]C=CC[O](2836)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.22e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=CC=CC[O](2837)'],
    products = ['C5H6O(219)(218)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R5H_DSMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = '331',
    isomers = [
        'C5H6O(219)(218)',
    ],
    reactants = [
        ('CO(10)(11)', 'C4H6(123)(122)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '331',
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

