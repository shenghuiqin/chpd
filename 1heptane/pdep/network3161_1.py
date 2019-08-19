species(
    label = 'O=C[CH]OC=[C]O(17417)',
    structure = SMILES('[O]C=COC=[C]O'),
    E0 = (-66.5647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.167,'amu*angstrom^2'), symmetry=1, barrier=(26.8316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1677,'amu*angstrom^2'), symmetry=1, barrier=(26.8477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16732,'amu*angstrom^2'), symmetry=1, barrier=(26.8389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.553089,0.056016,-8.78773e-06,-5.70682e-08,3.45145e-11,-7863.07,27.6718], Tmin=(100,'K'), Tmax=(896.458,'K')), NASAPolynomial(coeffs=[27.3263,-0.00879962,8.22804e-06,-1.72352e-09,1.17099e-13,-14859.1,-110.814], Tmin=(896.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.5647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'HCCOH(50)',
    structure = SMILES('C#CO'),
    E0 = (80.0402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,287.694,287.7,1588.42],'cm^-1')),
        HinderedRotor(inertia=(0.002037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(80.0402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[O][CH]C1OC=C1O(17579)',
    structure = SMILES('[O][CH]C1OC=C1O'),
    E0 = (-5.87992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702831,0.0530335,-5.18465e-06,-5.41443e-08,3.09696e-11,-570.334,22.7833], Tmin=(100,'K'), Tmax=(919.178,'K')), NASAPolynomial(coeffs=[25.6058,-0.00424099,4.89789e-06,-9.80275e-10,6.12933e-14,-7306.9,-107.002], Tmin=(919.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.87992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = '[O]C=COC=C=O(4755)',
    structure = SMILES('[O]C=COC=C=O'),
    E0 = (-124.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27169,'amu*angstrom^2'), symmetry=1, barrier=(29.2387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27825,'amu*angstrom^2'), symmetry=1, barrier=(29.3895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01813,0.0513186,-1.99376e-05,-2.98551e-08,2.06906e-11,-14893.2,23.6763], Tmin=(100,'K'), Tmax=(911.12,'K')), NASAPolynomial(coeffs=[21.6464,-0.00098806,3.19512e-06,-6.98058e-10,4.5656e-14,-20240.1,-82.6368], Tmin=(911.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + radical(C=COJ)"""),
)

species(
    label = '[O]C=COC#CO(17580)',
    structure = SMILES('[O]C=COC#CO'),
    E0 = (-29.6554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29855,'amu*angstrom^2'), symmetry=1, barrier=(29.8563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29464,'amu*angstrom^2'), symmetry=1, barrier=(29.7664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29682,'amu*angstrom^2'), symmetry=1, barrier=(29.8166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998246,0.0535163,-3.12077e-05,-1.26849e-08,1.2594e-11,-3447.18,23.877], Tmin=(100,'K'), Tmax=(942.978,'K')), NASAPolynomial(coeffs=[20.2843,0.00272106,2.59151e-07,-5.34875e-11,-1.49188e-15,-8463.31,-75.3458], Tmin=(942.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.6554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtOs) + radical(C=COJ)"""),
)

species(
    label = 'O=C=COC=[C]O(17581)',
    structure = SMILES('O=C=COC=[C]O'),
    E0 = (-26.5515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2120,512.5,787.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00349,'amu*angstrom^2'), symmetry=1, barrier=(23.0723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00407,'amu*angstrom^2'), symmetry=1, barrier=(23.0854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00248,'amu*angstrom^2'), symmetry=1, barrier=(23.049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0707664,0.0688726,-8.13176e-05,4.35154e-08,-8.46187e-12,-3036.84,28.7925], Tmin=(100,'K'), Tmax=(1489.35,'K')), NASAPolynomial(coeffs=[20.6547,-0.000733753,3.21231e-06,-7.79348e-10,5.72734e-14,-7579.56,-73.3791], Tmin=(1489.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.5515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + radical(C=CJO)"""),
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
    label = '[O]C=CO[CH]C=O(4757)',
    structure = SMILES('[O]C=COC=C[O]'),
    E0 = (-164.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39708,'amu*angstrom^2'), symmetry=1, barrier=(32.1216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39822,'amu*angstrom^2'), symmetry=1, barrier=(32.1478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75379,0.0473164,2.13002e-05,-8.91986e-08,4.5705e-11,-19687,24.5388], Tmin=(100,'K'), Tmax=(907.428,'K')), NASAPolynomial(coeffs=[28.1154,-0.00907459,8.35811e-06,-1.69861e-09,1.11423e-13,-27296.7,-119.376], Tmin=(907.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=CO[C]=CO(17582)',
    structure = SMILES('[O]C=CO[C]=CO'),
    E0 = (-66.5647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.167,'amu*angstrom^2'), symmetry=1, barrier=(26.8316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1677,'amu*angstrom^2'), symmetry=1, barrier=(26.8477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16732,'amu*angstrom^2'), symmetry=1, barrier=(26.8389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.553089,0.056016,-8.78773e-06,-5.70682e-08,3.45145e-11,-7863.07,27.6718], Tmin=(100,'K'), Tmax=(896.458,'K')), NASAPolynomial(coeffs=[27.3263,-0.00879962,8.22804e-06,-1.72352e-09,1.17099e-13,-14859.1,-110.814], Tmin=(896.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.5647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'O[C]=COC=[C]O(17583)',
    structure = SMILES('O[C]=COC=[C]O'),
    E0 = (31.7169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.98824,'amu*angstrom^2'), symmetry=1, barrier=(22.7216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989693,'amu*angstrom^2'), symmetry=1, barrier=(22.755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990268,'amu*angstrom^2'), symmetry=1, barrier=(22.7682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989689,'amu*angstrom^2'), symmetry=1, barrier=(22.7549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05264,0.0812503,-9.68049e-05,5.08412e-08,-9.49275e-12,4022.37,34.4643], Tmin=(100,'K'), Tmax=(1608.97,'K')), NASAPolynomial(coeffs=[23.9626,-0.00537383,6.73212e-06,-1.49736e-09,1.06128e-13,-864.593,-88.2884], Tmin=(1608.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.7169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=CO[C]=CO(17584)',
    structure = SMILES('O[C]=CO[C]=CO'),
    E0 = (31.7169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.98824,'amu*angstrom^2'), symmetry=1, barrier=(22.7216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989693,'amu*angstrom^2'), symmetry=1, barrier=(22.755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990268,'amu*angstrom^2'), symmetry=1, barrier=(22.7682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989689,'amu*angstrom^2'), symmetry=1, barrier=(22.7549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05264,0.0812503,-9.68049e-05,5.08412e-08,-9.49275e-12,4022.37,35.1574], Tmin=(100,'K'), Tmax=(1608.97,'K')), NASAPolynomial(coeffs=[23.9626,-0.00537383,6.73212e-06,-1.49736e-09,1.06128e-13,-864.593,-87.5952], Tmin=(1608.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.7169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC=CO(17585)',
    structure = SMILES('[O]C=[C]OC=CO'),
    E0 = (-66.5647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.167,'amu*angstrom^2'), symmetry=1, barrier=(26.8316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1677,'amu*angstrom^2'), symmetry=1, barrier=(26.8477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16732,'amu*angstrom^2'), symmetry=1, barrier=(26.8389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.553089,0.056016,-8.78773e-06,-5.70682e-08,3.45145e-11,-7863.07,27.6718], Tmin=(100,'K'), Tmax=(896.458,'K')), NASAPolynomial(coeffs=[27.3263,-0.00879962,8.22804e-06,-1.72352e-09,1.17099e-13,-14859.1,-110.814], Tmin=(896.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.5647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]OC=CO(4762)',
    structure = SMILES('O=[C][CH]OC=CO'),
    E0 = (-112.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13064,'amu*angstrom^2'), symmetry=1, barrier=(25.9956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12849,'amu*angstrom^2'), symmetry=1, barrier=(25.9463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13079,'amu*angstrom^2'), symmetry=1, barrier=(25.9992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13151,'amu*angstrom^2'), symmetry=1, barrier=(26.0157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95756,0.0896085,-0.000105421,5.31741e-08,-9.46183e-12,-13255,32.3237], Tmin=(100,'K'), Tmax=(1704.65,'K')), NASAPolynomial(coeffs=[27.3513,-0.00942452,8.34875e-06,-1.7333e-09,1.17969e-13,-18850.9,-111.812], Tmin=(1704.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(CCsJOC(O))"""),
)

species(
    label = 'O[C]=[C]OC=CO(17586)',
    structure = SMILES('O[C]=[C]OC=CO'),
    E0 = (31.7169,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3580,3650,1210,1345,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.98824,'amu*angstrom^2'), symmetry=1, barrier=(22.7216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989693,'amu*angstrom^2'), symmetry=1, barrier=(22.755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990268,'amu*angstrom^2'), symmetry=1, barrier=(22.7682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.989689,'amu*angstrom^2'), symmetry=1, barrier=(22.7549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05264,0.0812503,-9.68049e-05,5.08412e-08,-9.49275e-12,4022.37,35.1574], Tmin=(100,'K'), Tmax=(1608.97,'K')), NASAPolynomial(coeffs=[23.9626,-0.00537383,6.73212e-06,-1.49736e-09,1.06128e-13,-864.593,-87.5952], Tmin=(1608.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.7169,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CJO)"""),
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
    label = '[O]C=[C]O(7867)',
    structure = SMILES('[O]C=[C]O'),
    E0 = (79.4349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,3010,987.5,1337.5,450,1655,472.648],'cm^-1')),
        HinderedRotor(inertia=(0.000754551,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.49143,0.00601759,2.14158e-05,-2.72001e-08,9.0609e-12,9576.42,16.4916], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[6.08617,0.0113044,-5.70157e-06,1.19907e-09,-8.97995e-14,8107.39,-0.339445], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.4349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]O(172)',
    structure = SMILES('[CH]=[C]O'),
    E0 = (348.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1058.91,1059.8],'cm^-1')),
        HinderedRotor(inertia=(0.315045,'amu*angstrom^2'), symmetry=1, barrier=(7.24351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22339,0.0198737,-3.18466e-05,3.10436e-08,-1.17427e-11,41917.3,13.7606], Tmin=(100,'K'), Tmax=(829.31,'K')), NASAPolynomial(coeffs=[3.78101,0.0111997,-5.33323e-06,1.02849e-09,-7.13893e-14,42030.6,12.4155], Tmin=(829.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=COC1[CH]O1(17587)',
    structure = SMILES('O[C]=COC1[CH]O1'),
    E0 = (18.4996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.95699,0.0852782,-9.94837e-05,5.15574e-08,-9.30974e-12,2479.25,33.4834], Tmin=(100,'K'), Tmax=(1732.09,'K')), NASAPolynomial(coeffs=[21.1968,-0.00395136,8.75745e-06,-2.02293e-09,1.44095e-13,-177.49,-75.4547], Tmin=(1732.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.4996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(C=CJO)"""),
)

species(
    label = '[O]C1[CH]OC=C1O(17588)',
    structure = SMILES('[O]C1[CH]OC=C1O'),
    E0 = (-91.6008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.942405,0.0406902,4.16745e-05,-1.16841e-07,5.84365e-11,-10881.7,19.232], Tmin=(100,'K'), Tmax=(887.645,'K')), NASAPolynomial(coeffs=[29.6464,-0.0155853,1.3288e-05,-2.77845e-09,1.90837e-13,-18856.3,-132.043], Tmin=(887.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.6008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=C=COC=CO(4765)',
    structure = SMILES('O=C=COC=CO'),
    E0 = (-266.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.949793,0.0762712,-8.48029e-05,4.21361e-08,-7.47565e-12,-31821.5,29.3556], Tmin=(100,'K'), Tmax=(1696.54,'K')), NASAPolynomial(coeffs=[22.4404,-0.001969,4.7904e-06,-1.09314e-09,7.68104e-14,-36434.7,-86.0582], Tmin=(1696.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'OC1=COC=CO1(17589)',
    structure = SMILES('OC1=COC=CO1'),
    E0 = (-253.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59712,0.0424579,-1.26946e-05,-2.23102e-08,1.4114e-11,-30410.7,18.6884], Tmin=(100,'K'), Tmax=(934.026,'K')), NASAPolynomial(coeffs=[15.4359,0.00952959,-2.10892e-06,3.22909e-10,-2.42315e-14,-34144.7,-53.2811], Tmin=(934.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(14dioxin)"""),
)

species(
    label = '[O]C(C=O)C=[C]O(17416)',
    structure = SMILES('[O]C(C=O)C=[C]O'),
    E0 = (4.34714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888694,0.0577859,-5.88272e-05,2.92357e-08,-5.5371e-12,643.923,31.1961], Tmin=(100,'K'), Tmax=(1404.02,'K')), NASAPolynomial(coeffs=[16.5952,0.00895309,-2.29127e-06,3.18393e-10,-1.90382e-14,-3363.85,-48.4746], Tmin=(1404.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.34714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=OCOJ) + radical(C=CJO)"""),
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
    label = '[CH]=COC=[C]O(6252)',
    structure = SMILES('[CH]=COC=[C]O'),
    E0 = (247.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06326,'amu*angstrom^2'), symmetry=1, barrier=(24.4465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06135,'amu*angstrom^2'), symmetry=1, barrier=(24.4025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05996,'amu*angstrom^2'), symmetry=1, barrier=(24.3706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0663482,0.0643695,-7.13374e-05,3.57433e-08,-6.42812e-12,29978.9,29.0424], Tmin=(100,'K'), Tmax=(1661.83,'K')), NASAPolynomial(coeffs=[19.3834,-0.000519482,3.54642e-06,-8.41927e-10,6.01868e-14,26010.2,-67.1612], Tmin=(1661.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'O=[C]COC=[C]O(4758)',
    structure = SMILES('O=[C]COC=[C]O'),
    E0 = (-66.4735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,1855,455,950,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01195,'amu*angstrom^2'), symmetry=1, barrier=(23.2668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01202,'amu*angstrom^2'), symmetry=1, barrier=(23.2684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01146,'amu*angstrom^2'), symmetry=1, barrier=(23.2554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0111,'amu*angstrom^2'), symmetry=1, barrier=(23.2471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.826517,0.0755393,-8.32937e-05,4.0935e-08,-7.22246e-12,-7794.63,31.9118], Tmin=(100,'K'), Tmax=(1687.69,'K')), NASAPolynomial(coeffs=[23.0497,-0.00193049,4.11914e-06,-9.25482e-10,6.43462e-14,-12880.1,-86.963], Tmin=(1687.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-66.4735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = 'O=CCO[C]=[C]O(17590)',
    structure = SMILES('O=CCO[C]=[C]O'),
    E0 = (18.5998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1670,1700,300,440,3615,1277.5,1000,2782.5,750,1395,475,1775,1000,238.085,239.05,239.44],'cm^-1')),
        HinderedRotor(inertia=(0.462686,'amu*angstrom^2'), symmetry=1, barrier=(18.738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.464712,'amu*angstrom^2'), symmetry=1, barrier=(18.7493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.462808,'amu*angstrom^2'), symmetry=1, barrier=(18.7489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.464571,'amu*angstrom^2'), symmetry=1, barrier=(18.7479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0472902,0.0699094,-7.72357e-05,3.92064e-08,-7.30037e-12,2399.17,31.5443], Tmin=(100,'K'), Tmax=(1538.14,'K')), NASAPolynomial(coeffs=[21.1545,0.00186483,1.71005e-06,-4.66685e-10,3.4667e-14,-2596.11,-74.9071], Tmin=(1538.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.5998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=COCC=O(4761)',
    structure = SMILES('[O][C]=COCC=O'),
    E0 = (-79.6818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,205.789,205.951,207.632,208.569],'cm^-1')),
        HinderedRotor(inertia=(0.767041,'amu*angstrom^2'), symmetry=1, barrier=(23.4449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78147,'amu*angstrom^2'), symmetry=1, barrier=(23.4432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779098,'amu*angstrom^2'), symmetry=1, barrier=(23.4295,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940589,0.0520301,-1.54337e-05,-3.36728e-08,2.11811e-11,-9459.25,26.2725], Tmin=(100,'K'), Tmax=(929.61,'K')), NASAPolynomial(coeffs=[21.5233,0.0026118,1.14034e-06,-2.59288e-10,1.28692e-14,-14977.5,-80.6228], Tmin=(929.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.6818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'OC1=CO[CH][CH]O1(17591)',
    structure = SMILES('O[C]1[CH]OC=CO1'),
    E0 = (-84.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990373,0.0424865,3.25261e-05,-1.07324e-07,5.61102e-11,-10009.4,16.9229], Tmin=(100,'K'), Tmax=(869.318,'K')), NASAPolynomial(coeffs=[28.2496,-0.0148646,1.40188e-05,-3.04827e-09,2.16078e-13,-17321.1,-125.565], Tmin=(869.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(23dihydro14dioxin) + radical(Cs_P) + radical(CCsJOC(O))"""),
)

species(
    label = 'O=C=COCC=O(4766)',
    structure = SMILES('O=C=COCC=O'),
    E0 = (-279.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982001,0.0544209,-3.02474e-05,-1.2637e-08,1.22567e-11,-33486,22.3888], Tmin=(100,'K'), Tmax=(942.02,'K')), NASAPolynomial(coeffs=[19.2668,0.0064494,-1.10464e-06,1.7272e-10,-1.57766e-14,-38247.3,-71.72], Tmin=(942.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'O=CC1OC=C1O(17592)',
    structure = SMILES('O=CC1OC=C1O'),
    E0 = (-333.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21467,0.0315067,6.61471e-05,-1.35134e-07,6.14745e-11,-39930.4,20.7433], Tmin=(100,'K'), Tmax=(918.098,'K')), NASAPolynomial(coeffs=[28.8021,-0.0111508,9.16152e-06,-1.76728e-09,1.10214e-13,-48263.8,-127.789], Tmin=(918.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutene)"""),
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
    E0 = (-66.5647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-5.87992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (117.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (197.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (217.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (107.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (96.2145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (167.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (194.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (182.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (327.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-33.5245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (95.5433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (325.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (329.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (115.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-28.4903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-58.1967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-58.3641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (247.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (490.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (144.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (67.1005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (178.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-39.6825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-56.4619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-27.3397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-59.0335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['HCCOH(50)', 'OCHCHO(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['[O][CH]C1OC=C1O(17579)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(60.6848,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 60.2 to 60.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[O]C=COC=C=O(4755)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[O]C=COC#CO(17580)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'O=C=COC=[C]O(17581)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HCCOH(50)', '[O]C=C[O](2537)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(46.5868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Ct_Ct;OJ_sec] for rate rule [Ct_Ct;O_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['[O]C=CO[CH]C=O(4757)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=CO[C]=CO(17582)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O[C]=COC=[C]O(17583)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O[C]=CO[C]=CO(17584)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=[C]OC=CO(17585)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;Cd_rad_out_Cd;Cd_H_out_singleNd]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['O=[C][CH]OC=CO(4762)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O[C]=[C]OC=CO(17586)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(313415,'s^-1'), n=1.7968, Ea=(63.8264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;XH_out] + [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHCHO(47)', '[O]C=[C]O(7867)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.63841e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/OneDe;Y_rad] for rate rule [O_rad/OneDe;Cd_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]O(172)', '[O]C=C[O](2537)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O_rad/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['O[C]=COC1[CH]O1(17587)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['[O]C1[CH]OC=C1O(17588)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.53e+07,'s^-1'), n=1.05, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleNd]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['O=C=COC=CO(4765)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['OC1=COC=CO1(17589)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6;CdsingleND_rad_out;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['[O]C(C=O)C=[C]O(17416)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O(4)', '[CH]=COC=[C]O(6252)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]O(172)', 'OCHCHO(48)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.2635,'m^3/(mol*s)'), n=2.48, Ea=(21.7568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-DeH;YJ] for rate rule [Od_CO-DeH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=[C]COC=[C]O(4758)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=CCO[C]=[C]O(17590)'],
    products = ['O=C[CH]OC=[C]O(17417)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Cd_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['[O][C]=COCC=O(4761)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_1H;XH_out] for rate rule [R5HJ_3;C_rad_out_H/OneDe;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['OC1=CO[CH][CH]O1(17591)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.43857e+20,'s^-1'), n=-1.34645, Ea=(10.1028,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cdsingleNd] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cdsingleNd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['O=C=COCC=O(4766)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O=C[CH]OC=[C]O(17417)'],
    products = ['O=CC1OC=C1O(17592)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriND_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

network(
    label = '3161',
    isomers = [
        'O=C[CH]OC=[C]O(17417)',
    ],
    reactants = [
        ('HCCOH(50)', 'OCHCHO(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '3161',
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

