species(
    label = '[CH2]C(O)C[CH]OO(11218)',
    structure = SMILES('[CH2]C(O)C[CH]OO'),
    E0 = (-4.92534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.251556,0.0879027,-0.00011764,8.8609e-08,-2.71829e-11,-462.437,29.6743], Tmin=(100,'K'), Tmax=(793.541,'K')), NASAPolynomial(coeffs=[10.9451,0.0339958,-1.57348e-05,2.99067e-09,-2.07504e-13,-2159.47,-19.4423], Tmin=(793.541,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.92534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(CJCO)"""),
)

species(
    label = 'CH2CHOOH(64)',
    structure = SMILES('C=COO'),
    E0 = (-53.0705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.754187,'amu*angstrom^2'), symmetry=1, barrier=(17.3402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754176,'amu*angstrom^2'), symmetry=1, barrier=(17.34,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79821,0.0274377,-2.03468e-05,7.62127e-09,-1.12671e-12,-6315.53,9.11829], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.82519,0.0274167,-2.04456e-05,7.77399e-09,-1.18661e-12,-6325.27,8.96641], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), E0=(-53.0705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CH2CHOOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C=C(O)C[CH]OO(16114)',
    structure = SMILES('C=C(O)C[CH]OO'),
    E0 = (-130.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0455811,'amu*angstrom^2'), symmetry=1, barrier=(16.6987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726491,'amu*angstrom^2'), symmetry=1, barrier=(16.7035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308966,'amu*angstrom^2'), symmetry=1, barrier=(16.6609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11197,'amu*angstrom^2'), symmetry=1, barrier=(41.3241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724899,'amu*angstrom^2'), symmetry=1, barrier=(16.6669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283316,0.0780214,-8.35263e-05,4.53785e-08,-9.70102e-12,-15597,27.0916], Tmin=(100,'K'), Tmax=(1142.34,'K')), NASAPolynomial(coeffs=[16.6227,0.0208076,-8.39905e-06,1.53443e-09,-1.05762e-13,-19330,-53.9109], Tmin=(1142.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C(O)C=COO(16115)',
    structure = SMILES('[CH2]C(O)C=COO'),
    E0 = (-62.0372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,3615,1277.5,1000,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.670381,'amu*angstrom^2'), symmetry=1, barrier=(15.4134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672623,'amu*angstrom^2'), symmetry=1, barrier=(15.4649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0229231,'amu*angstrom^2'), symmetry=1, barrier=(15.4655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.669864,'amu*angstrom^2'), symmetry=1, barrier=(15.4015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962915,'amu*angstrom^2'), symmetry=1, barrier=(15.4423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259057,0.0735875,-7.3771e-05,3.67776e-08,-7.13607e-12,-7319.32,31.3114], Tmin=(100,'K'), Tmax=(1263.95,'K')), NASAPolynomial(coeffs=[18.2887,0.0165295,-6.05737e-06,1.0623e-09,-7.18733e-14,-11877,-59.8945], Tmin=(1263.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.0372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]OO(104)',
    structure = SMILES('[CH2][CH]OO'),
    E0 = (224.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00920734,'amu*angstrom^2'), symmetry=1, barrier=(3.53679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00921023,'amu*angstrom^2'), symmetry=1, barrier=(3.53685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43223,'amu*angstrom^2'), symmetry=1, barrier=(32.9297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52392,0.0279389,-1.79903e-05,3.58397e-09,4.58838e-13,27095.5,18.6054], Tmin=(100,'K'), Tmax=(1150.47,'K')), NASAPolynomial(coeffs=[9.41961,0.0107482,-4.42273e-06,8.47943e-10,-6.05179e-14,25059.9,-17.5802], Tmin=(1150.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(CJCOOH)"""),
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
    label = '[CH2]C(O)CC=O(2311)',
    structure = SMILES('[CH2]C(O)CC=O'),
    E0 = (-187.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,309.222],'cm^-1')),
        HinderedRotor(inertia=(0.00176301,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114062,'amu*angstrom^2'), symmetry=1, barrier=(7.73966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114063,'amu*angstrom^2'), symmetry=1, barrier=(7.73967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114065,'amu*angstrom^2'), symmetry=1, barrier=(7.73964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4805,0.0589791,-6.40814e-05,4.16954e-08,-1.15216e-11,-22521.7,23.683], Tmin=(100,'K'), Tmax=(862.75,'K')), NASAPolynomial(coeffs=[7.76492,0.0298419,-1.34219e-05,2.54901e-09,-1.77961e-13,-23606.1,-5.7078], Tmin=(862.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO)"""),
)

species(
    label = 'C=CC[CH]OO(1931)',
    structure = SMILES('C=CC[CH]OO'),
    E0 = (83.7382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,180,1544.02],'cm^-1')),
        HinderedRotor(inertia=(0.68936,'amu*angstrom^2'), symmetry=1, barrier=(15.8497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.685354,'amu*angstrom^2'), symmetry=1, barrier=(15.7576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818608,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60568,0.0554494,-4.7278e-05,2.24922e-08,-4.60428e-12,10155.2,22.2921], Tmin=(100,'K'), Tmax=(1121.97,'K')), NASAPolynomial(coeffs=[8.60356,0.0305008,-1.39231e-05,2.67291e-09,-1.88084e-13,8584.96,-12.2739], Tmin=(1121.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.7382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOH)"""),
)

species(
    label = 'C[C](O)C[CH]OO(11214)',
    structure = SMILES('C[C](O)C[CH]OO'),
    E0 = (-39.8866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.136054,0.0929362,-0.000140904,1.20619e-07,-4.10499e-11,-4665.75,28.87], Tmin=(100,'K'), Tmax=(830.32,'K')), NASAPolynomial(coeffs=[9.0714,0.0372835,-1.75897e-05,3.32291e-09,-2.2759e-13,-5715,-9.95938], Tmin=(830.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.8866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(O)[CH]COO(16116)',
    structure = SMILES('[CH2]C(O)[CH]COO'),
    E0 = (6.90712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.289869,0.0869842,-0.00011596,8.67542e-08,-2.64181e-11,959.389,31.6929], Tmin=(100,'K'), Tmax=(799.421,'K')), NASAPolynomial(coeffs=[10.9889,0.0334415,-1.54787e-05,2.94513e-09,-2.04563e-13,-750.943,-17.5269], Tmin=(799.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.90712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCO) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(O)CCO[O](16117)',
    structure = SMILES('[CH2]C(O)CCO[O]'),
    E0 = (-41.5031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3615,1277.5,1000,3000,3100,440,815,1455,1000,492.5,1135,1000,282.891,285.789],'cm^-1')),
        HinderedRotor(inertia=(0.136798,'amu*angstrom^2'), symmetry=1, barrier=(7.90344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138946,'amu*angstrom^2'), symmetry=1, barrier=(7.89744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138959,'amu*angstrom^2'), symmetry=1, barrier=(7.89298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210228,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137005,'amu*angstrom^2'), symmetry=1, barrier=(7.89846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476647,0.0820401,-0.000107405,8.07371e-08,-2.4781e-11,-4869.04,29.4724], Tmin=(100,'K'), Tmax=(793.067,'K')), NASAPolynomial(coeffs=[10.206,0.0329678,-1.45901e-05,2.71525e-09,-1.86019e-13,-6412.25,-15.2105], Tmin=(793.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.5031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C[CH]OO(1990)',
    structure = SMILES('CC([O])C[CH]OO'),
    E0 = (13.8464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4219.77,'J/mol'), sigma=(7.06962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=659.12 K, Pc=27.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.743165,0.0765238,-8.37386e-05,5.28236e-08,-1.40297e-11,1778.42,27.0919], Tmin=(100,'K'), Tmax=(897.413,'K')), NASAPolynomial(coeffs=[9.69334,0.0366321,-1.70631e-05,3.29385e-09,-2.32244e-13,171.958,-15.1191], Tmin=(897.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.8464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCsJOOH)"""),
)

species(
    label = 'CC(O)[CH][CH]OO(11217)',
    structure = SMILES('CC(O)[CH][CH]OO'),
    E0 = (-16.0995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386205,0.0853251,-0.000116283,9.26132e-08,-3.03695e-11,-1811.7,30.8573], Tmin=(100,'K'), Tmax=(766.548,'K')), NASAPolynomial(coeffs=[9.41267,0.0364213,-1.70607e-05,3.2532e-09,-2.25757e-13,-3142.61,-9.94512], Tmin=(766.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.0995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2][C](O)CCOO(16118)',
    structure = SMILES('[CH2][C](O)CCOO'),
    E0 = (-16.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0113054,0.0952729,-0.000143356,1.19057e-07,-3.93105e-11,-1892.5,29.884], Tmin=(100,'K'), Tmax=(828.049,'K')), NASAPolynomial(coeffs=[10.7306,0.034152,-1.59158e-05,2.99232e-09,-2.04478e-13,-3355.02,-18.002], Tmin=(828.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])CCOO(11220)',
    structure = SMILES('[CH2]C([O])CCOO'),
    E0 = (36.8529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.570797,0.0792216,-8.77716e-05,5.3773e-08,-1.35704e-11,4552.69,28.1914], Tmin=(100,'K'), Tmax=(951.205,'K')), NASAPolynomial(coeffs=[11.5139,0.0332059,-1.52102e-05,2.9194e-09,-2.05391e-13,2470.78,-24.0556], Tmin=(951.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.8529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC(O)C[CH]O[O](11221)',
    structure = SMILES('CC(O)C[CH]O[O]'),
    E0 = (-64.5097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596576,0.0800672,-0.000106437,8.45825e-08,-2.76841e-11,-7641.13,28.5543], Tmin=(100,'K'), Tmax=(799.055,'K')), NASAPolynomial(coeffs=[8.59936,0.0360033,-1.62059e-05,3.03157e-09,-2.07917e-13,-8792.28,-7.45967], Tmin=(799.055,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.5097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(O)CCOO(16119)',
    structure = SMILES('C=C(O)CCOO'),
    E0 = (-319.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0368762,0.0781008,-7.63264e-05,3.75394e-08,-7.2251e-12,-38265,27.829], Tmin=(100,'K'), Tmax=(1270.33,'K')), NASAPolynomial(coeffs=[18.6329,0.0195464,-7.18643e-06,1.25522e-09,-8.44771e-14,-42989.6,-66.3358], Tmin=(1270.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC(O)C=COO(11223)',
    structure = SMILES('CC(O)C=COO'),
    E0 = (-273.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.508002,0.0672988,-5.03851e-05,1.32319e-08,7.40714e-13,-32775.7,29.7818], Tmin=(100,'K'), Tmax=(1041.69,'K')), NASAPolynomial(coeffs=[16.7455,0.0214935,-8.25139e-06,1.51424e-09,-1.06384e-13,-37056.3,-53.5266], Tmin=(1041.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-273.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'O[CH]CC[CH]OO(13278)',
    structure = SMILES('O[CH]CC[CH]OO'),
    E0 = (-23.5736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0573856,0.0947098,-0.000144506,1.23317e-07,-4.17287e-11,-2700.96,29.6226], Tmin=(100,'K'), Tmax=(833.251,'K')), NASAPolynomial(coeffs=[9.48692,0.0367478,-1.73097e-05,3.2642e-09,-2.23208e-13,-3831.66,-11.5047], Tmin=(833.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.5736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCsJOH) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C(O)C([CH2])OO(11192)',
    structure = SMILES('[CH2]C(O)C([CH2])OO'),
    E0 = (7.8118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3615,1310,387.5,850,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3615,1277.5,1000,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4236.76,'J/mol'), sigma=(7.09873,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.77 K, Pc=26.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242615,0.0849607,-0.00010035,6.27998e-08,-1.57545e-11,1072.94,31.2075], Tmin=(100,'K'), Tmax=(969.299,'K')), NASAPolynomial(coeffs=[14.0024,0.0281761,-1.24722e-05,2.35648e-09,-1.64398e-13,-1594.43,-34.7459], Tmin=(969.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.8118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCO) + radical(CJCOOH)"""),
)

species(
    label = 'OOC1CC(O)C1(13279)',
    structure = SMILES('OOC1CC(O)C1'),
    E0 = (-259.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05737,0.0506926,-7.93698e-07,-3.60721e-08,1.73415e-11,-31099.9,24.5197], Tmin=(100,'K'), Tmax=(1001.25,'K')), NASAPolynomial(coeffs=[15.7231,0.024289,-9.45629e-06,1.80099e-09,-1.31175e-13,-35650,-54.309], Tmin=(1001.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-259.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]C(O)CC([O])O(16120)',
    structure = SMILES('[CH2]C(O)CC([O])O'),
    E0 = (-239.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.105,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.601662,0.080641,-0.000110432,9.02046e-08,-3.04e-11,-28693.1,31.5195], Tmin=(100,'K'), Tmax=(777.619,'K')), NASAPolynomial(coeffs=[8.38508,0.0364249,-1.70797e-05,3.26107e-09,-2.26389e-13,-29777.3,-3.26097], Tmin=(777.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJCO)"""),
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
    label = 'O[CH]C[CH]OO(7926)',
    structure = SMILES('O[CH]C[CH]OO'),
    E0 = (0.206615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3615,1310,387.5,850,1000,3615,1277.5,1000,300.439,1898.37],'cm^-1')),
        HinderedRotor(inertia=(0.150217,'amu*angstrom^2'), symmetry=1, barrier=(9.57581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149534,'amu*angstrom^2'), symmetry=1, barrier=(9.58091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149712,'amu*angstrom^2'), symmetry=1, barrier=(9.58507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712464,'amu*angstrom^2'), symmetry=1, barrier=(45.5624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717609,'amu*angstrom^2'), symmetry=1, barrier=(45.5655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0779,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.733362,0.0795627,-0.000129255,1.13114e-07,-3.84742e-11,135.083,24.9514], Tmin=(100,'K'), Tmax=(853.56,'K')), NASAPolynomial(coeffs=[8.50021,0.028446,-1.35586e-05,2.54671e-09,-1.72795e-13,-654.612,-8.14842], Tmin=(853.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.206615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCsJOOH) + radical(CCsJOH)"""),
)

species(
    label = '[CH]OO(168)',
    structure = SMILES('[CH]OO'),
    E0 = (320.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00655552,'amu*angstrom^2'), symmetry=1, barrier=(3.39708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148635,'amu*angstrom^2'), symmetry=1, barrier=(16.5265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06576,0.0188176,-2.02137e-05,1.03271e-08,-2.01695e-12,38527.3,11.7484], Tmin=(100,'K'), Tmax=(1262.94,'K')), NASAPolynomial(coeffs=[8.15329,0.00270404,-1.0753e-06,2.24439e-10,-1.70827e-14,37242.2,-13.9836], Tmin=(1262.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])O(2453)',
    structure = SMILES('[CH2]C([CH2])O'),
    E0 = (132.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.00138042,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853629,'amu*angstrom^2'), symmetry=1, barrier=(7.40251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853972,'amu*angstrom^2'), symmetry=1, barrier=(7.39995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05987,0.0394237,-3.55746e-05,1.71593e-08,-3.30943e-12,15990.3,17.5693], Tmin=(100,'K'), Tmax=(1252.95,'K')), NASAPolynomial(coeffs=[10.0721,0.0138446,-4.95163e-06,8.65246e-10,-5.82583e-14,13982.5,-22.8917], Tmin=(1252.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJCO) + radical(CJCO)"""),
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
    E0 = (-4.92534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (92.2615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (156.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (94.8037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (121.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-4.92534,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (196.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (97.7542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (139.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (124.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (152.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (113.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (137.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (154.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (68.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (368.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (58.4748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (58.4748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (153.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (166.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (3.35898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (87.8198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (381.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (452.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['CH2CHOOH(64)', 'CH2CHOH(42)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=C(O)C[CH]OO(16114)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]C(O)C=COO(16115)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]OO(104)', 'CH2CHOH(42)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(8.71744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CHOOH(64)', '[CH2][CH]O(284)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00668046,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', '[CH2]C(O)CC=O(2311)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(154.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 150.6 to 154.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(5)', 'C=CC[CH]OO(1931)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00674586,'m^3/(mol*s)'), n=2.3625, Ea=(84.2031,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;OJ] for rate rule [Cds-CsH_Cds-HH;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['C[C](O)C[CH]OO(11214)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)[CH]COO(16116)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3114.97,'s^-1'), n=2.79033, Ea=(132.312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)CCO[O](16117)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 242 used for R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC([O])C[CH]OO(1990)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['CC(O)[CH][CH]OO(11217)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['[CH2][C](O)CCOO(16118)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])CCOO(11220)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.94e+12,'s^-1'), n=-0.24, Ea=(117.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['CC(O)C[CH]O[O](11221)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]OO(104)', '[CH2][CH]O(284)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['C=C(O)CCOO(16119)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['CC(O)C=COO(11223)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['O[CH]CC[CH]OO(13278)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)C([CH2])OO(11192)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['OOC1CC(O)C1(13279)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(O)C[CH]OO(11218)'],
    products = ['[CH2]C(O)CC([O])O(16120)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_H/NonDeC] for rate rule [ROOH;C_rad_out_H/NonDeC]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(19)', 'O[CH]C[CH]OO(7926)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]OO(168)', '[CH2]C([CH2])O(2453)'],
    products = ['[CH2]C(O)C[CH]OO(11218)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2631',
    isomers = [
        '[CH2]C(O)C[CH]OO(11218)',
    ],
    reactants = [
        ('CH2CHOOH(64)', 'CH2CHOH(42)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2631',
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

