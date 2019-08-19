species(
    label = 'O=C[CH]OO[CH]C=O(1312)',
    structure = SMILES('O=C[CH]OO[CH]C=O'),
    E0 = (-78.8979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,3000,3050,390,425,1340,1360,335,370,325.144],'cm^-1')),
        HinderedRotor(inertia=(0.00160997,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647473,'amu*angstrom^2'), symmetry=1, barrier=(48.0889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.650852,'amu*angstrom^2'), symmetry=1, barrier=(48.0893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646056,'amu*angstrom^2'), symmetry=1, barrier=(48.114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65351,'amu*angstrom^2'), symmetry=1, barrier=(48.0675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04095,0.071564,-8.42897e-05,5.59919e-08,-1.57016e-11,-9388.28,27.2427], Tmin=(100,'K'), Tmax=(849.401,'K')), NASAPolynomial(coeffs=[9.08222,0.0336952,-1.74139e-05,3.5021e-09,-2.52265e-13,-10754.3,-10.2391], Tmin=(849.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.8979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(OCJC=O)"""),
)

species(
    label = 'OCHCHO(48)',
    structure = SMILES('O=CC=O'),
    E0 = (-225.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,857.339,857.946,861.419,1717.87],'cm^-1')),
        HinderedRotor(inertia=(0.00964493,'amu*angstrom^2'), symmetry=1, barrier=(5.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-225.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]C1[CH]OOC1C=O(4837)',
    structure = SMILES('[O]C1[CH]OOC1C=O'),
    E0 = (3.2769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47472,0.0386033,2.76591e-05,-7.70304e-08,3.7127e-11,500.986,25.4803], Tmin=(100,'K'), Tmax=(897.943,'K')), NASAPolynomial(coeffs=[18.3154,0.00920124,5.72238e-07,-3.44542e-10,2.51187e-14,-4362.47,-64.1938], Tmin=(897.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.2769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CC(C)OJ) + radical(CCsJOOC)"""),
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
    label = 'O=C=COO[CH]C=O(4838)',
    structure = SMILES('O=C=COO[CH]C=O'),
    E0 = (-35.1931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,258.357],'cm^-1')),
        HinderedRotor(inertia=(0.341242,'amu*angstrom^2'), symmetry=1, barrier=(16.1599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.508631,'amu*angstrom^2'), symmetry=1, barrier=(24.1205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09935,'amu*angstrom^2'), symmetry=1, barrier=(52.0507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09975,'amu*angstrom^2'), symmetry=1, barrier=(52.0517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00557,0.069656,-8.47173e-05,5.37072e-08,-1.37302e-11,-4128.09,26.7388], Tmin=(100,'K'), Tmax=(947.032,'K')), NASAPolynomial(coeffs=[11.9153,0.0235759,-1.17303e-05,2.32713e-09,-1.6661e-13,-6194.43,-25.3003], Tmin=(947.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.1931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + radical(OCJC=O)"""),
)

species(
    label = 'O=[C]COO[CH]C=O(4839)',
    structure = SMILES('O=[C]COO[CH]C=O'),
    E0 = (-62.7118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,190.418],'cm^-1')),
        HinderedRotor(inertia=(0.00463783,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121965,'amu*angstrom^2'), symmetry=1, barrier=(44.4609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.967925,'amu*angstrom^2'), symmetry=1, barrier=(24.9275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72433,'amu*angstrom^2'), symmetry=1, barrier=(44.4618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72633,'amu*angstrom^2'), symmetry=1, barrier=(44.4608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16776,0.0683521,-7.54406e-05,4.56234e-08,-1.1616e-11,-7445.58,27.7074], Tmin=(100,'K'), Tmax=(929.608,'K')), NASAPolynomial(coeffs=[9.66362,0.0317959,-1.6455e-05,3.32267e-09,-2.40215e-13,-9025.16,-12.6602], Tmin=(929.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.7118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(OCJC=O)"""),
)

species(
    label = 'O=[C][CH]OOCC=O(4840)',
    structure = SMILES('O=[C][CH]OOCC=O'),
    E0 = (-62.7118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1855,455,950,2750,2850,1437.5,1250,1305,750,350,190.418],'cm^-1')),
        HinderedRotor(inertia=(0.00463783,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121965,'amu*angstrom^2'), symmetry=1, barrier=(44.4609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.967925,'amu*angstrom^2'), symmetry=1, barrier=(24.9275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72433,'amu*angstrom^2'), symmetry=1, barrier=(44.4618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.72633,'amu*angstrom^2'), symmetry=1, barrier=(44.4608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16776,0.0683521,-7.54406e-05,4.56234e-08,-1.1616e-11,-7445.58,27.7074], Tmin=(100,'K'), Tmax=(929.608,'K')), NASAPolynomial(coeffs=[9.66362,0.0317959,-1.6455e-05,3.32267e-09,-2.40215e-13,-9025.16,-12.6602], Tmin=(929.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.7118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(OCJC=O)"""),
)

species(
    label = '[O]C=C[O](2537)',
    structure = SMILES('[O]C=C[O]'),
    E0 = (-18.8461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,714.947,715.113],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7018,-0.0027881,5.18222e-05,-5.96735e-08,2.03775e-11,-2247.83,13.3237], Tmin=(100,'K'), Tmax=(1035.3,'K')), NASAPolynomial(coeffs=[6.65087,0.0113883,-5.76543e-06,1.26602e-09,-9.87778e-14,-4228.83,-7.62441], Tmin=(1035.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'O=C[CH]OOC1[CH]O1(4841)',
    structure = SMILES('O=C[CH]OOC1[CH]O1'),
    E0 = (22.2613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475488,0.073177,-7.53467e-05,3.59003e-08,-5.18931e-12,2808.38,27.6865], Tmin=(100,'K'), Tmax=(901.657,'K')), NASAPolynomial(coeffs=[15.9612,0.0178903,-5.68393e-06,8.90322e-10,-5.63652e-14,-529.358,-48.443], Tmin=(901.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.2613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(OCJC=O) + radical(CCsJO)"""),
)

species(
    label = 'O=CC1O[CH][CH]OO1(4842)',
    structure = SMILES('O=CC1O[CH][CH]OO1'),
    E0 = (-41.5739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.958087,0.0725943,-7.2052e-05,3.52809e-08,-6.19296e-12,-4789.9,30.9941], Tmin=(100,'K'), Tmax=(1768.25,'K')), NASAPolynomial(coeffs=[15.0212,0.0108842,1.98154e-06,-8.07013e-10,6.51576e-14,-6444.56,-43.9041], Tmin=(1768.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.5739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(124trioxane) + radical(CCsJOOC) + radical(CCsJOCs)"""),
)

species(
    label = 'O=C=COOCC=O(4843)',
    structure = SMILES('O=C=COOCC=O'),
    E0 = (-173.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.984847,0.0673052,-7.09585e-05,3.82469e-08,-8.29601e-12,-20780.8,27.0481], Tmin=(100,'K'), Tmax=(1108.25,'K')), NASAPolynomial(coeffs=[13.2673,0.0229741,-1.09567e-05,2.15274e-09,-1.53807e-13,-23503.2,-33.4699], Tmin=(1108.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'O=CC1OOC1C=O(1313)',
    structure = SMILES('O=CC1OOC1C=O'),
    E0 = (-248.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68928,0.0376839,9.73746e-06,-3.91228e-08,1.68674e-11,-29850.6,24.7644], Tmin=(100,'K'), Tmax=(1026.29,'K')), NASAPolynomial(coeffs=[14.0811,0.0198348,-8.67732e-06,1.74773e-09,-1.3048e-13,-33997.6,-43.1526], Tmin=(1026.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-OdCsH) + ring(12dioxetane)"""),
)

species(
    label = 'CHCHO(47)',
    structure = SMILES('[CH]=C[O]'),
    E0 = (245.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(245.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]O[CH]C=O(3604)',
    structure = SMILES('[O]O[CH]C=O'),
    E0 = (38.8345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.304071,'amu*angstrom^2'), symmetry=1, barrier=(6.9912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17371,'amu*angstrom^2'), symmetry=1, barrier=(56.7408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48924,0.0359331,-5.02567e-05,4.15061e-08,-1.40856e-11,4722.53,16.7287], Tmin=(100,'K'), Tmax=(782.809,'K')), NASAPolynomial(coeffs=[5.94539,0.0159797,-7.6283e-06,1.46013e-09,-1.01349e-13,4251.69,1.34978], Tmin=(782.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.8345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]C1OOC1C=O(4844)',
    structure = SMILES('[O][CH]C1OOC1C=O'),
    E0 = (77.4075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48752,0.0609228,-6.91061e-05,4.84894e-08,-1.49278e-11,9395.37,24.8239], Tmin=(100,'K'), Tmax=(767.113,'K')), NASAPolynomial(coeffs=[6.57267,0.0344085,-1.72635e-05,3.43757e-09,-2.46311e-13,8615.15,1.63891], Tmin=(767.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.4075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C1[CH]OOC=CO1(4845)',
    structure = SMILES('[O]C1[CH]OOC=CO1'),
    E0 = (32.2129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0738,0.00764031,0.00019839,-3.15174e-07,1.36856e-10,4033.02,25.569], Tmin=(100,'K'), Tmax=(903.019,'K')), NASAPolynomial(coeffs=[45.1698,-0.0383337,2.66675e-05,-5.24166e-09,3.44985e-13,-10020.3,-216.388], Tmin=(903.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.2129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCOJ) + radical(CCsJOOC)"""),
)

species(
    label = 'O=C[CH]OOC=[C]O(4846)',
    structure = SMILES('O=C[CH]OOC=[C]O'),
    E0 = (23.0753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,350,500,795,815,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.79245,'amu*angstrom^2'), symmetry=1, barrier=(41.2119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01591,'amu*angstrom^2'), symmetry=1, barrier=(23.3578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01582,'amu*angstrom^2'), symmetry=1, barrier=(23.3557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01589,'amu*angstrom^2'), symmetry=1, barrier=(23.3573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01586,'amu*angstrom^2'), symmetry=1, barrier=(23.3567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146751,0.079208,-9.15107e-05,5.11105e-08,-1.10006e-11,2918.84,32.1325], Tmin=(100,'K'), Tmax=(1146.77,'K')), NASAPolynomial(coeffs=[19.0571,0.0132475,-5.2328e-06,9.53389e-10,-6.61339e-14,-1418.32,-61.6889], Tmin=(1146.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.0753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(OCJC=O)"""),
)

species(
    label = 'O=C[CH]OO[C]=CO(4847)',
    structure = SMILES('O=C[CH]OO[C]=CO'),
    E0 = (23.0753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,350,500,795,815,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.79245,'amu*angstrom^2'), symmetry=1, barrier=(41.2119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01591,'amu*angstrom^2'), symmetry=1, barrier=(23.3578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01582,'amu*angstrom^2'), symmetry=1, barrier=(23.3557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01589,'amu*angstrom^2'), symmetry=1, barrier=(23.3573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01586,'amu*angstrom^2'), symmetry=1, barrier=(23.3567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146751,0.079208,-9.15107e-05,5.11105e-08,-1.10006e-11,2918.84,32.1325], Tmin=(100,'K'), Tmax=(1146.77,'K')), NASAPolynomial(coeffs=[19.0571,0.0132475,-5.2328e-06,9.53389e-10,-6.61339e-14,-1418.32,-61.6889], Tmin=(1146.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.0753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OOCC=O(4848)',
    structure = SMILES('[O]C=[C]OOCC=O'),
    E0 = (26.0531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1685,370,2782.5,750,1395,475,1775,1000,2750,2850,1437.5,1250,1305,750,350,486.504,624.043],'cm^-1')),
        HinderedRotor(inertia=(0.854158,'amu*angstrom^2'), symmetry=1, barrier=(19.6388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30824,'amu*angstrom^2'), symmetry=1, barrier=(30.079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.852047,'amu*angstrom^2'), symmetry=1, barrier=(19.5902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853003,'amu*angstrom^2'), symmetry=1, barrier=(19.6122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633345,0.0686012,-6.91325e-05,3.41194e-08,-6.59641e-12,3259.24,32.0415], Tmin=(100,'K'), Tmax=(1258.76,'K')), NASAPolynomial(coeffs=[16.9451,0.016767,-7.36449e-06,1.40579e-09,-9.92418e-14,-847.284,-50.4069], Tmin=(1258.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.0531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'O=[C][CH]OOC=CO(4849)',
    structure = SMILES('O=[C][CH]OOC=CO'),
    E0 = (-61.9979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3615,1277.5,1000,1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.28105,'amu*angstrom^2'), symmetry=1, barrier=(29.4538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27784,'amu*angstrom^2'), symmetry=1, barrier=(29.3802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28214,'amu*angstrom^2'), symmetry=1, barrier=(29.479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27927,'amu*angstrom^2'), symmetry=1, barrier=(29.4129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28176,'amu*angstrom^2'), symmetry=1, barrier=(29.4702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.498778,0.083475,-9.36432e-05,4.87024e-08,-9.49192e-12,-7281.43,32.0038], Tmin=(100,'K'), Tmax=(1382.19,'K')), NASAPolynomial(coeffs=[24.1534,0.0048998,-5.21436e-07,1.4575e-12,1.50893e-15,-13405.3,-92.4088], Tmin=(1382.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-61.9979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(OCJC=O)"""),
)

species(
    label = '[CH]1[CH]OOC=COO1(4850)',
    structure = SMILES('[CH]1[CH]OOC=COO1'),
    E0 = (314.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67845,0.033686,3.56675e-05,-7.004e-08,2.83354e-11,37929.9,21.8742], Tmin=(100,'K'), Tmax=(1000.62,'K')), NASAPolynomial(coeffs=[15.3134,0.0213337,-9.00674e-06,1.82555e-09,-1.38601e-13,33091,-54.4599], Tmin=(1000.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(CCsJOOC) + radical(CCsJOOC)"""),
)

species(
    label = 'O=C=COOC=CO(4851)',
    structure = SMILES('O=C=COOC=CO'),
    E0 = (-172.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.905527,0.0848116,-9.64592e-05,4.96322e-08,-9.3264e-12,-20606.3,32.1664], Tmin=(100,'K'), Tmax=(1527.12,'K')), NASAPolynomial(coeffs=[25.4064,-0.000277783,3.00226e-06,-7.21765e-10,5.20007e-14,-26757.1,-99.7397], Tmin=(1527.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'C1=COO1(3609)',
    structure = SMILES('C1=COO1'),
    E0 = (143.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31094,-0.00195743,7.27275e-05,-1.01582e-07,4.10226e-11,17349.4,10.0685], Tmin=(100,'K'), Tmax=(927.669,'K')), NASAPolynomial(coeffs=[13.7806,-0.00309387,3.40676e-06,-6.26613e-10,3.47414e-14,13513.3,-49.8617], Tmin=(927.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene)"""),
)

species(
    label = 'O=CC1OC=COO1(4852)',
    structure = SMILES('O=CC1OC=COO1'),
    E0 = (-302.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3226,0.0341215,5.443e-05,-1.15135e-07,5.26334e-11,-36318.5,22.446], Tmin=(100,'K'), Tmax=(912.946,'K')), NASAPolynomial(coeffs=[24.0606,-0.000597569,4.83265e-06,-1.04303e-09,6.5833e-14,-43175,-99.9945], Tmin=(912.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-302.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(124trioxene)"""),
)

species(
    label = '[O]C(C=O)O[CH]C=O(1311)',
    structure = SMILES('[O]C=COC([O])C=O'),
    E0 = (-290.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386483,0.0624346,-2.59843e-05,-3.11712e-08,2.19752e-11,-34774.4,29.7735], Tmin=(100,'K'), Tmax=(928.8,'K')), NASAPolynomial(coeffs=[24.7494,0.00133247,1.92598e-06,-4.08465e-10,2.24894e-14,-41190.2,-96.1391], Tmin=(928.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=COJ)"""),
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
    label = '[CH]=COO[CH]C=O(4726)',
    structure = SMILES('[CH]=COO[CH]C=O'),
    E0 = (239.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(1.35389,'amu*angstrom^2'), symmetry=1, barrier=(31.1286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.55003,'amu*angstrom^2'), symmetry=1, barrier=(58.6301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07673,'amu*angstrom^2'), symmetry=1, barrier=(24.7561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0777,'amu*angstrom^2'), symmetry=1, barrier=(24.7785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21264,0.0614842,-6.35255e-05,3.33307e-08,-7.0227e-12,28871.8,25.7253], Tmin=(100,'K'), Tmax=(1140.52,'K')), NASAPolynomial(coeffs=[12.864,0.0206207,-9.78214e-06,1.91607e-09,-1.36656e-13,26214,-32.0179], Tmin=(1140.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(Cds_P)"""),
)

species(
    label = '[O]C=[C]OOC=CO(4853)',
    structure = SMILES('[O]C=[C]OOC=CO'),
    E0 = (26.7669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.11764,'amu*angstrom^2'), symmetry=1, barrier=(25.6968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11757,'amu*angstrom^2'), symmetry=1, barrier=(25.6951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11713,'amu*angstrom^2'), symmetry=1, barrier=(25.6849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11791,'amu*angstrom^2'), symmetry=1, barrier=(25.7031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38017,0.0873461,-9.81342e-05,4.91358e-08,-8.87255e-12,3439.75,37.6183], Tmin=(100,'K'), Tmax=(1629.96,'K')), NASAPolynomial(coeffs=[26.4302,-0.0025714,4.55563e-06,-1.02142e-09,7.14565e-14,-2747.67,-101.308], Tmin=(1629.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.7669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C1=COOC=COO1(4854)',
    structure = SMILES('C1=COOC=COO1'),
    E0 = (88.6062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94233,0.0153194,0.000106484,-1.60398e-07,6.50548e-11,10758.4,20.7601], Tmin=(100,'K'), Tmax=(944.394,'K')), NASAPolynomial(coeffs=[22.0859,0.0061429,1.2085e-07,5.82072e-11,-2.10368e-14,3558.19,-93.2457], Tmin=(944.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.6062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclooctane)"""),
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
    E0 = (-78.8979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (4.03547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (201.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (89.7819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-1.53589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-37.6922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (104.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-3.19106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-39.6729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-70.9901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (284.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (77.4075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (32.2129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-78.7966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (185.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (173.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (70.3617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (17.2377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (314.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-53.9246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (125.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-71.3667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (234.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (482.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (144.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (88.6062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['OCHCHO(48)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[O]C1[CH]OOC1C=O(4837)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.27145e+08,'s^-1'), n=0.688067, Ea=(82.9333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_csHDe] + [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra_H;radadd_intra_csHDe]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'O=C=COO[CH]C=O(4838)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.3e-15,'cm^3/(molecule*s)'), n=1.43, Ea=(25.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ck_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C]COO[CH]C=O(4839)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.67823e+08,'s^-1'), n=1.48018, Ea=(152.494,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out] for rate rule [R2H_S;CO_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C][CH]OOCC=O(4840)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(60863,'s^-1'), n=1.86284, Ea=(61.1759,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;CO_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C=C[O](2537)', '[O]C=C[O](2537)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['O=C[CH]OOC1[CH]O1(4841)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.85443e+11,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_csHO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['O=CC1O[CH][CH]OO1(4842)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.01802e+09,'s^-1'), n=0.463766, Ea=(75.7068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHDe] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_csHCO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['O=C=COOCC=O(4843)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['O=CC1OOC1C=O(1313)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CHCHO(47)', '[O]O[CH]C=O(3604)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[O][CH]C1OOC1C=O(4844)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(156.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHDe]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 154.3 to 156.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[O]C1[CH]OOC=CO1(4845)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.6502e+11,'s^-1'), n=0.157, Ea=(111.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7plus;multiplebond_intra;radadd_intra_O] + [R8;multiplebond_intra;radadd_intra] for rate rule [R8;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 103.3 to 111.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['OCHCHO(48)', '[O]C=C[O](2537)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(165.35,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_R;O_rad/OneDe] for rate rule [Od_CO-DeH;O_rad/OneDe]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C[CH]OOC=[C]O(4846)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C[CH]OO[C]=CO(4847)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C=[C]OOCC=O(4848)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C][CH]OOC=CO(4849)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(364667,'s^-1'), n=1.22214, Ea=(79.2357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;Y_rad_out;XH_out] for rate rule [R7HJ_1;CO_rad_out;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[CH]1[CH]OOC=COO1(4850)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(393.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 387.9 to 393.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['O=C=COOC=CO(4851)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[O]C=C[O](2537)', 'C1=COO1(3609)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.22709e+11,'s^-1'), n=0, Ea=(203.965,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO_SD;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation
Ea raised from 202.3 to 204.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['O=CC1OC=COO1(4852)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['[O]C(C=O)O[CH]C=O(1311)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', '[CH]=COO[CH]C=O(4726)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C=[C]OOC=CO(4853)'],
    products = ['O=C[CH]OO[CH]C=O(1312)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.78681e+06,'s^-1'), n=1.58912, Ea=(118.22,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C[CH]OO[CH]C=O(1312)'],
    products = ['C1=COOC=COO1(4854)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.35773e+13,'s^-1'), n=0.0154583, Ea=(167.504,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Ypri_rad_out] for rate rule [R8;O_rad;Opri_rad]
Euclidian distance = 1.73205080757
family: Birad_recombination
Ea raised from 159.7 to 167.5 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = '898',
    isomers = [
        'O=C[CH]OO[CH]C=O(1312)',
    ],
    reactants = [
        ('OCHCHO(48)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '898',
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

