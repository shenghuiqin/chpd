species(
    label = '[CH2]C(O)C=[C]O(1257)',
    structure = SMILES('[CH2]C(O)C=[C]O'),
    E0 = (53.2121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,340.613],'cm^-1')),
        HinderedRotor(inertia=(0.167235,'amu*angstrom^2'), symmetry=1, barrier=(13.7673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167187,'amu*angstrom^2'), symmetry=1, barrier=(13.7671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167209,'amu*angstrom^2'), symmetry=1, barrier=(13.7671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167247,'amu*angstrom^2'), symmetry=1, barrier=(13.7672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.372397,0.0645561,-6.7258e-05,3.36046e-08,-6.27118e-12,6543.59,29.9932], Tmin=(100,'K'), Tmax=(1500.67,'K')), NASAPolynomial(coeffs=[18.5905,0.0069173,-5.70202e-07,-5.26759e-11,7.46619e-15,2097.98,-61.8877], Tmin=(1500.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CJCO)"""),
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
    label = 'CH2CHOH(42)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C(O)C=C=O(5748)',
    structure = SMILES('[CH2]C(O)C=C=O'),
    E0 = (-70.9418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.513805,'amu*angstrom^2'), symmetry=1, barrier=(11.8134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513547,'amu*angstrom^2'), symmetry=1, barrier=(11.8075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.513802,'amu*angstrom^2'), symmetry=1, barrier=(11.8133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09397,0.0684221,-0.000104183,8.46073e-08,-2.68003e-11,-8431.98,22.2013], Tmin=(100,'K'), Tmax=(877.307,'K')), NASAPolynomial(coeffs=[9.6909,0.0210003,-9.03908e-06,1.62111e-09,-1.07068e-13,-9623.88,-16.3446], Tmin=(877.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-70.9418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C=[C]O(8492)',
    structure = SMILES('C=C(O)C=[C]O'),
    E0 = (-89.0201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.999644,'amu*angstrom^2'), symmetry=1, barrier=(22.9838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998892,'amu*angstrom^2'), symmetry=1, barrier=(22.9665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999021,'amu*angstrom^2'), symmetry=1, barrier=(22.9695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383059,0.0711318,-8.04555e-05,4.12386e-08,-7.59281e-12,-10527,25.9678], Tmin=(100,'K'), Tmax=(1617.64,'K')), NASAPolynomial(coeffs=[20.9015,-0.000781229,4.10752e-06,-9.80442e-10,7.0779e-14,-14890.4,-79.1573], Tmin=(1617.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.0201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)C#CO(8493)',
    structure = SMILES('[CH2]C(O)C#CO'),
    E0 = (46.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,878.846],'cm^-1')),
        HinderedRotor(inertia=(0.235002,'amu*angstrom^2'), symmetry=1, barrier=(14.2152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618138,'amu*angstrom^2'), symmetry=1, barrier=(14.2122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618625,'amu*angstrom^2'), symmetry=1, barrier=(14.2234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0599738,'amu*angstrom^2'), symmetry=1, barrier=(32.7983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.385,0.0523673,-4.84832e-05,1.99743e-08,-2.20397e-12,5720.87,24.3182], Tmin=(100,'K'), Tmax=(975.911,'K')), NASAPolynomial(coeffs=[13.3741,0.0135123,-4.5707e-06,7.75857e-10,-5.23304e-14,2891.03,-35.7393], Tmin=(975.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(46.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CtH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CJCO)"""),
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
    label = '[CH2][CH]O(284)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (143.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.217851,'amu*angstrom^2'), symmetry=1, barrier=(5.00882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197382,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84769,0.0229973,-1.85068e-05,7.77211e-09,-1.30725e-12,17300.5,13.0245], Tmin=(100,'K'), Tmax=(1418.34,'K')), NASAPolynomial(coeffs=[7.97636,0.00853337,-3.21015e-06,5.82141e-10,-3.99268e-14,15845.7,-13.5107], Tmin=(1418.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = 'OH(5)',
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
    label = 'C=CC=[C]O(6299)',
    structure = SMILES('C=CC=[C]O'),
    E0 = (124.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.01736,'amu*angstrom^2'), symmetry=1, barrier=(23.3911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02157,'amu*angstrom^2'), symmetry=1, barrier=(23.4879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85281,0.034453,6.0896e-06,-4.53851e-08,2.39475e-11,15108,18.6912], Tmin=(100,'K'), Tmax=(905.911,'K')), NASAPolynomial(coeffs=[16.7736,0.00282296,1.74834e-06,-4.54027e-10,3.0294e-14,10999.1,-59.5759], Tmin=(905.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)[CH]C=O(5749)',
    structure = SMILES('[CH2]C(O)C=C[O]'),
    E0 = (-45.0694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,342.207,343.883],'cm^-1')),
        HinderedRotor(inertia=(0.202523,'amu*angstrom^2'), symmetry=1, barrier=(17.2094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203817,'amu*angstrom^2'), symmetry=1, barrier=(17.2066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206057,'amu*angstrom^2'), symmetry=1, barrier=(17.2149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18004,0.0487819,-1.27556e-05,-2.98389e-08,1.8178e-11,-5306.9,25.3697], Tmin=(100,'K'), Tmax=(936.532,'K')), NASAPolynomial(coeffs=[18.491,0.00827652,-1.42509e-06,2.11294e-10,-1.83928e-14,-10015.5,-64.838], Tmin=(936.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.0694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(C=COJ)"""),
)

species(
    label = 'C[C](O)C=[C]O(8494)',
    structure = SMILES('C[C](O)C=[C]O'),
    E0 = (18.2509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972592,0.061107,-6.05567e-05,2.61298e-08,-2.92074e-12,2309.07,26.6231], Tmin=(100,'K'), Tmax=(941.955,'K')), NASAPolynomial(coeffs=[15.5035,0.0124202,-3.75765e-06,6.0293e-10,-3.99375e-14,-1005.99,-45.677], Tmin=(941.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.2509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C2CsJOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(O)[C]=CO(8495)',
    structure = SMILES('[CH2]C(O)[C]=CO'),
    E0 = (51.3097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.740968,'amu*angstrom^2'), symmetry=1, barrier=(17.0363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740852,'amu*angstrom^2'), symmetry=1, barrier=(17.0337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740913,'amu*angstrom^2'), symmetry=1, barrier=(17.0351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740964,'amu*angstrom^2'), symmetry=1, barrier=(17.0362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.14328,0.0705709,-7.61139e-05,3.81902e-08,-7.02559e-12,6337.98,29.1924], Tmin=(100,'K'), Tmax=(1568.99,'K')), NASAPolynomial(coeffs=[20.7818,0.00328938,1.53133e-06,-4.61961e-10,3.51086e-14,1486.95,-75.7183], Tmin=(1568.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.3097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C=[C]O(8496)',
    structure = SMILES('CC([O])C=[C]O'),
    E0 = (71.9838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,292.99,292.999],'cm^-1')),
        HinderedRotor(inertia=(0.249523,'amu*angstrom^2'), symmetry=1, barrier=(15.2041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249513,'amu*angstrom^2'), symmetry=1, barrier=(15.2037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249573,'amu*angstrom^2'), symmetry=1, barrier=(15.2042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28424,0.0484372,-1.77753e-05,-2.09587e-08,1.42468e-11,8765.8,25.8885], Tmin=(100,'K'), Tmax=(938.001,'K')), NASAPolynomial(coeffs=[16.8127,0.0105554,-2.51251e-06,4.00854e-10,-3.00345e-14,4606.03,-54.6784], Tmin=(938.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.9838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC(O)[C]=[C]O(8497)',
    structure = SMILES('CC(O)[C]=[C]O'),
    E0 = (79.4648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,279.906],'cm^-1')),
        HinderedRotor(inertia=(0.230552,'amu*angstrom^2'), symmetry=1, barrier=(12.818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230553,'amu*angstrom^2'), symmetry=1, barrier=(12.8181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230553,'amu*angstrom^2'), symmetry=1, barrier=(12.8181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230552,'amu*angstrom^2'), symmetry=1, barrier=(12.8181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722331,0.0613,-6.31822e-05,3.2134e-08,-6.22339e-12,9684.62,28.5735], Tmin=(100,'K'), Tmax=(1385.35,'K')), NASAPolynomial(coeffs=[16.8048,0.00996794,-2.30053e-06,2.85017e-10,-1.55335e-14,5698.49,-52.5616], Tmin=(1385.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.4648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=C(O)[CH][CH]O(5780)',
    structure = SMILES('[CH2]C(O)=C[CH]O'),
    E0 = (-121.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18723,0.0477655,-7.2595e-06,-3.90019e-08,2.27369e-11,-14484.8,24.1675], Tmin=(100,'K'), Tmax=(915.501,'K')), NASAPolynomial(coeffs=[19.3094,0.00581796,4.67114e-07,-2.06695e-10,1.24237e-14,-19363.2,-70.1831], Tmin=(915.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])C=CO(8498)',
    structure = SMILES('[CH2]C([O])C=CO'),
    E0 = (43.8287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.87872,'amu*angstrom^2'), symmetry=1, barrier=(20.2035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880118,'amu*angstrom^2'), symmetry=1, barrier=(20.2356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880096,'amu*angstrom^2'), symmetry=1, barrier=(20.2351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924302,0.0518325,-1.04575e-05,-4.09345e-08,2.45243e-11,5396.84,24.6862], Tmin=(100,'K'), Tmax=(917.65,'K')), NASAPolynomial(coeffs=[21.7504,0.00267698,1.85252e-06,-4.46979e-10,2.75119e-14,-177.936,-83.5471], Tmin=(917.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.8287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC(O)[CH][C]=O(3213)',
    structure = SMILES('CC(O)[CH][C]=O'),
    E0 = (-39.7112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69544,0.0523151,-5.3797e-05,3.1682e-08,-7.76293e-12,-4694.53,25.5875], Tmin=(100,'K'), Tmax=(975.767,'K')), NASAPolynomial(coeffs=[8.61451,0.0239517,-1.01956e-05,1.89276e-09,-1.30732e-13,-6044.82,-7.62335], Tmin=(975.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.7112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C=CO(5784)',
    structure = SMILES('C=C(O)C=CO'),
    E0 = (-328.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579041,0.0556429,-5.74162e-06,-5.97426e-08,3.54133e-11,-39399.5,19.3796], Tmin=(100,'K'), Tmax=(893.729,'K')), NASAPolynomial(coeffs=[26.6725,-0.00674585,7.6733e-06,-1.64806e-09,1.12979e-13,-46236,-115.728], Tmin=(893.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC(O)C=C=O(3214)',
    structure = SMILES('CC(O)C=C=O'),
    E0 = (-282.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30006,0.0625659,-8.19601e-05,6.20721e-08,-1.91467e-11,-33886.4,20.8302], Tmin=(100,'K'), Tmax=(814.265,'K')), NASAPolynomial(coeffs=[8.6719,0.025141,-1.07859e-05,1.97211e-09,-1.33496e-13,-35046.8,-12.9734], Tmin=(814.265,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'O[C]=CC[CH]O(1255)',
    structure = SMILES('O[C]=CC[CH]O'),
    E0 = (34.9002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1685,370,3010,987.5,1337.5,450,1655,311.231,311.254],'cm^-1')),
        HinderedRotor(inertia=(0.181118,'amu*angstrom^2'), symmetry=1, barrier=(12.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181105,'amu*angstrom^2'), symmetry=1, barrier=(12.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181108,'amu*angstrom^2'), symmetry=1, barrier=(12.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181104,'amu*angstrom^2'), symmetry=1, barrier=(12.4496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4129.77,'J/mol'), sigma=(6.69615,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=645.06 K, Pc=31.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86678,0.0658295,-7.55763e-05,4.25469e-08,-8.80929e-12,4313.16,26.5646], Tmin=(100,'K'), Tmax=(926.301,'K')), NASAPolynomial(coeffs=[14.768,0.0138416,-4.41148e-06,7.00613e-10,-4.46589e-14,1392.83,-41.2987], Tmin=(926.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.9002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCsJOH)"""),
)

species(
    label = 'OC1=CC(O)C1(8499)',
    structure = SMILES('OC1=CC(O)C1'),
    E0 = (-242.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70002,0.0298393,4.56077e-05,-9.38361e-08,4.18318e-11,-29111.9,20.1871], Tmin=(100,'K'), Tmax=(928.648,'K')), NASAPolynomial(coeffs=[20.1069,0.00386597,1.45045e-06,-3.18036e-10,1.40431e-14,-34829.3,-79.6294], Tmin=(928.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
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
    label = 'O[C]=C[CH]O(7929)',
    structure = SMILES('O[C]=C[CH]O'),
    E0 = (1.23461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,451.152],'cm^-1')),
        HinderedRotor(inertia=(0.136799,'amu*angstrom^2'), symmetry=1, barrier=(19.541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134977,'amu*angstrom^2'), symmetry=1, barrier=(19.5306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134838,'amu*angstrom^2'), symmetry=1, barrier=(19.5524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06242,0.0302119,8.48388e-06,-4.48245e-08,2.30673e-11,229.844,21.4611], Tmin=(100,'K'), Tmax=(911.13,'K')), NASAPolynomial(coeffs=[16.1372,0.00129417,1.97256e-06,-4.61847e-10,2.97312e-14,-3699.44,-52.6196], Tmin=(911.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.23461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=CCJO)"""),
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
    E0 = (53.2121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (171.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (132.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (273.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (218.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (246.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (154.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (215.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (151.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (285.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (210.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (231.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (203.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (463.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (126.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (491.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (116.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (78.1854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (211.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (61.4965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (382.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['HCCOH(50)', 'CH2CHOH(42)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[CH2]C(O)C=C=O(5748)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.185e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C(O)C=[C]O(8492)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-OneDeOs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]C(O)C#CO(8493)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C]O(172)', 'CH2CHOH(42)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(8.71744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HCCOH(50)', '[CH2][CH]O(284)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(5)', 'C=CC=[C]O(6299)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(47.2566,'m^3/(mol*s)'), n=1.48028, Ea=(1.31417,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;OJ_pri] + [Cds-OneDeH_Cds;YJ] for rate rule [Cds-OneDeH_Cds;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['[CH2]C(O)[CH]C=O(5749)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['C[C](O)C=[C]O(8494)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)[C]=CO(8495)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC([O])C=[C]O(8496)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC(O)[C]=[C]O(8497)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['C=C(O)[CH][CH]O(5780)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])C=CO(8498)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleNd]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['CC(O)[CH][C]=O(3213)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]O(284)', '[CH]=[C]O(172)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['C=C(O)C=CO(5784)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['CC(O)C=C=O(3214)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['O[C]=CC[CH]O(1255)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)C=[C]O(1257)'],
    products = ['OC1=CC(O)C1(8499)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(19)', 'O[C]=C[CH]O(7929)'],
    products = ['[CH2]C(O)C=[C]O(1257)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '842',
    isomers = [
        '[CH2]C(O)C=[C]O(1257)',
    ],
    reactants = [
        ('HCCOH(50)', 'CH2CHOH(42)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '842',
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

