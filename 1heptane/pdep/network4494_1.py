species(
    label = 'C#C[CH]CC([CH2])O(26526)',
    structure = SMILES('C#C[CH]CC([CH2])O'),
    E0 = (287.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,371.501,3606.23],'cm^-1')),
        HinderedRotor(inertia=(0.114857,'amu*angstrom^2'), symmetry=1, barrier=(11.2788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114719,'amu*angstrom^2'), symmetry=1, barrier=(11.2797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837763,'amu*angstrom^2'), symmetry=1, barrier=(81.7469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837239,'amu*angstrom^2'), symmetry=1, barrier=(81.7482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837332,'amu*angstrom^2'), symmetry=1, barrier=(81.7359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486979,0.0764537,-8.14569e-05,4.42898e-08,-8.33599e-12,34687.6,27.8559], Tmin=(100,'K'), Tmax=(832.808,'K')), NASAPolynomial(coeffs=[13.3624,0.0253093,-8.60445e-06,1.39334e-09,-8.85159e-14,32172.2,-34.1318], Tmin=(832.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO) + radical(Sec_Propargyl)"""),
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
    label = 'CH2CHCCH(26391)',
    structure = SMILES('C#CC=C'),
    E0 = (274.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.46338,'amu*angstrom^2'), symmetry=1, barrier=(33.6459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87083,0.0182042,1.06711e-05,-2.72492e-08,1.19478e-11,33023.8,11.2934], Tmin=(100,'K'), Tmax=(955.249,'K')), NASAPolynomial(coeffs=[8.52653,0.0108962,-3.56564e-06,6.31243e-10,-4.51891e-14,31196.2,-19.6435], Tmin=(955.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH2CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C1[CH]CC(O)C1(27847)',
    structure = SMILES('[CH]C1=CCC(O)C1'),
    E0 = (181.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.388,0.0406734,3.55189e-05,-7.25226e-08,2.98525e-11,21897.1,23.571], Tmin=(100,'K'), Tmax=(977.557,'K')), NASAPolynomial(coeffs=[13.8879,0.0303405,-1.12535e-05,2.08497e-09,-1.50238e-13,17503,-46.4247], Tmin=(977.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH]CC(=C)O(27848)',
    structure = SMILES('C#C[CH]CC(=C)O'),
    E0 = (161.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3615,1277.5,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,347.607,349.16],'cm^-1')),
        HinderedRotor(inertia=(0.188638,'amu*angstrom^2'), symmetry=1, barrier=(16.6837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85356,'amu*angstrom^2'), symmetry=1, barrier=(75.7529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194696,'amu*angstrom^2'), symmetry=1, barrier=(16.7365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267148,'amu*angstrom^2'), symmetry=1, barrier=(23.9723,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162327,0.0745926,-7.55548e-05,3.83664e-08,-7.31337e-12,19582.9,27.7193], Tmin=(100,'K'), Tmax=(1489.94,'K')), NASAPolynomial(coeffs=[18.482,0.0128342,-1.59595e-06,1.40068e-12,8.63281e-15,15326.3,-65.3034], Tmin=(1489.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC=CC([CH2])O(27849)',
    structure = SMILES('C#CC=CC([CH2])O'),
    E0 = (250.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,216.188],'cm^-1')),
        HinderedRotor(inertia=(0.617148,'amu*angstrom^2'), symmetry=1, barrier=(20.4714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618865,'amu*angstrom^2'), symmetry=1, barrier=(20.4677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615663,'amu*angstrom^2'), symmetry=1, barrier=(20.4693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616711,'amu*angstrom^2'), symmetry=1, barrier=(20.4685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.618576,0.0645046,-4.64046e-05,8.23683e-09,2.91874e-12,30296.9,27.6453], Tmin=(100,'K'), Tmax=(1009.95,'K')), NASAPolynomial(coeffs=[17.1916,0.0183745,-6.86668e-06,1.2648e-09,-9.00721e-14,25954.4,-57.4], Tmin=(1009.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CJCO)"""),
)

species(
    label = '[CH]=[C]C=C(4699)',
    structure = SMILES('[CH]=C=C[CH2]'),
    E0 = (451.584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,180,1024.85,1025.53,1026.61],'cm^-1')),
        HinderedRotor(inertia=(0.00938781,'amu*angstrom^2'), symmetry=1, barrier=(7.01846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76805,0.020302,8.75519e-06,-2.87666e-08,1.37354e-11,54363.7,13.5565], Tmin=(100,'K'), Tmax=(915.031,'K')), NASAPolynomial(coeffs=[9.46747,0.00887314,-1.78262e-06,2.38534e-10,-1.6263e-14,52390.1,-22.2544], Tmin=(915.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C=CJ) + radical(Allyl_P)"""),
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
    label = 'C#C[CH]CC=C(26789)',
    structure = SMILES('C#C[CH]CC=C'),
    E0 = (376.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,410.876,410.906],'cm^-1')),
        HinderedRotor(inertia=(0.0812548,'amu*angstrom^2'), symmetry=1, barrier=(9.73456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63521,'amu*angstrom^2'), symmetry=1, barrier=(76.1042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635141,'amu*angstrom^2'), symmetry=1, barrier=(76.1043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59773,0.0469098,-2.16011e-05,-7.36055e-09,7.51589e-12,45316,21.3454], Tmin=(100,'K'), Tmax=(911.508,'K')), NASAPolynomial(coeffs=[10.894,0.0220604,-6.94892e-06,1.11507e-09,-7.2543e-14,42958.8,-26.2758], Tmin=(911.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]C[C](C)O(27850)',
    structure = SMILES('C#C[CH]C[C](C)O'),
    E0 = (252.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396281,0.0811552,-0.000103371,7.42862e-08,-2.12312e-11,30483.3,26.9651], Tmin=(100,'K'), Tmax=(929.011,'K')), NASAPolynomial(coeffs=[11.3516,0.0288444,-1.06084e-05,1.7619e-09,-1.11687e-13,28669.6,-23.8872], Tmin=(929.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C2CsJOH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C([CH2])O(27795)',
    structure = SMILES('C#CC[CH]C([CH2])O'),
    E0 = (341.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,761.987,762.073],'cm^-1')),
        HinderedRotor(inertia=(0.00675102,'amu*angstrom^2'), symmetry=1, barrier=(2.79662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.57133,'amu*angstrom^2'), symmetry=1, barrier=(13.136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.571201,'amu*angstrom^2'), symmetry=1, barrier=(13.133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126311,'amu*angstrom^2'), symmetry=1, barrier=(2.90414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.23578,'amu*angstrom^2'), symmetry=1, barrier=(74.397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.477195,0.072943,-7.37923e-05,4.01619e-08,-8.70248e-12,41148.6,30.5408], Tmin=(100,'K'), Tmax=(1125.51,'K')), NASAPolynomial(coeffs=[14.2973,0.0238272,-8.33413e-06,1.38937e-09,-9.0253e-14,38037.6,-37.767], Tmin=(1125.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO) + radical(CCJCO)"""),
)

species(
    label = 'C#C[CH]CC(C)[O](27231)',
    structure = SMILES('C#C[CH]CC(C)[O]'),
    E0 = (306.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,374.452,374.487,374.553],'cm^-1')),
        HinderedRotor(inertia=(0.00120219,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721652,'amu*angstrom^2'), symmetry=1, barrier=(7.18241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768161,'amu*angstrom^2'), symmetry=1, barrier=(76.4321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768038,'amu*angstrom^2'), symmetry=1, barrier=(76.4331,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.570647,0.0701303,-6.64419e-05,3.48881e-08,-7.37471e-12,36945.9,26.721], Tmin=(100,'K'), Tmax=(1148.47,'K')), NASAPolynomial(coeffs=[13.2425,0.0259954,-8.79773e-06,1.42659e-09,-9.07558e-14,34035.3,-36.1676], Tmin=(1148.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#C[CH][CH]C(C)O(27851)',
    structure = SMILES('[CH]=C=C[CH]C(C)O'),
    E0 = (199.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258,0.0752143,-7.11155e-05,3.52159e-08,-6.93066e-12,24098.9,26.7979], Tmin=(100,'K'), Tmax=(1234.1,'K')), NASAPolynomial(coeffs=[16.1057,0.023849,-8.68381e-06,1.49045e-09,-9.87677e-14,20187.3,-52.9919], Tmin=(1234.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C#CCC[C]([CH2])O(27852)',
    structure = SMILES('C#CCC[C]([CH2])O'),
    E0 = (317.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,360,370,350,317.243,318.543],'cm^-1')),
        HinderedRotor(inertia=(0.125787,'amu*angstrom^2'), symmetry=1, barrier=(9.03171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131434,'amu*angstrom^2'), symmetry=1, barrier=(9.03141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126111,'amu*angstrom^2'), symmetry=1, barrier=(9.01821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00167432,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.987449,'amu*angstrom^2'), symmetry=1, barrier=(68.3906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263127,0.0852674,-0.000113667,8.42699e-08,-2.49113e-11,38350.2,28.2402], Tmin=(100,'K'), Tmax=(876.429,'K')), NASAPolynomial(coeffs=[11.7683,0.0292045,-1.16338e-05,2.03104e-09,-1.33165e-13,36470,-24.9694], Tmin=(876.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[C]#CCCC([CH2])O(27853)',
    structure = SMILES('[C]#CCCC([CH2])O'),
    E0 = (478.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665547,0.0769556,-8.56424e-05,4.64093e-08,-6.56984e-12,57641.2,28.1409], Tmin=(100,'K'), Tmax=(702.505,'K')), NASAPolynomial(coeffs=[11.4066,0.0288502,-1.07981e-05,1.83241e-09,-1.18752e-13,55809.9,-22.1781], Tmin=(702.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CJCO)"""),
)

species(
    label = 'C#CCCC([CH2])[O](27854)',
    structure = SMILES('C#CCCC([CH2])[O]'),
    E0 = (371.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,328.221,330.229,3959.27],'cm^-1')),
        HinderedRotor(inertia=(1.09716,'amu*angstrom^2'), symmetry=1, barrier=(85.264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231271,'amu*angstrom^2'), symmetry=1, barrier=(17.8829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231248,'amu*angstrom^2'), symmetry=1, barrier=(17.8804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09305,'amu*angstrom^2'), symmetry=1, barrier=(85.2543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.545453,0.0729651,-7.22327e-05,3.89307e-08,-8.45299e-12,44808.2,27.6089], Tmin=(100,'K'), Tmax=(1115.15,'K')), NASAPolynomial(coeffs=[13.4377,0.0267208,-1.00288e-05,1.74341e-09,-1.1613e-13,41932.8,-35.994], Tmin=(1115.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[C]#C[CH]CC(C)O(27855)',
    structure = SMILES('[C]#C[CH]CC(C)O'),
    E0 = (412.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86306,0.0719536,-7.14684e-05,2.98942e-08,8.19068e-13,49771.4,26.6425], Tmin=(100,'K'), Tmax=(694.749,'K')), NASAPolynomial(coeffs=[11.0047,0.0284683,-9.76158e-06,1.56089e-09,-9.70926e-14,48002.6,-21.18], Tmin=(694.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C(O)CC1[C]=C1(27856)',
    structure = SMILES('[CH2]C(O)CC1[C]=C1'),
    E0 = (460.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640982,0.0763381,-8.4932e-05,5.25045e-08,-1.32658e-11,55547.3,25.6682], Tmin=(100,'K'), Tmax=(955.216,'K')), NASAPolynomial(coeffs=[11.5426,0.0306868,-1.32443e-05,2.47188e-09,-1.7116e-13,53464.7,-26.4265], Tmin=(955.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(CJCO)"""),
)

species(
    label = 'OC1C[CH][C]=CC1(27831)',
    structure = SMILES('OC1C[CH][C]=CC1'),
    E0 = (182.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81546,0.030713,4.93896e-05,-8.29019e-08,3.27923e-11,22068.4,19.3898], Tmin=(100,'K'), Tmax=(981.202,'K')), NASAPolynomial(coeffs=[13.2749,0.0267867,-1.0022e-05,1.90936e-09,-1.40769e-13,17759.8,-46.174], Tmin=(981.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C#CCCC(=C)O(27857)',
    structure = SMILES('C#CCCC(=C)O'),
    E0 = (15.2447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.263879,0.0747624,-6.96489e-05,3.26367e-08,-5.80568e-12,2003.33,28.259], Tmin=(100,'K'), Tmax=(1570.01,'K')), NASAPolynomial(coeffs=[19.1933,0.015179,-3.15759e-06,3.41272e-10,-1.63849e-14,-2872.43,-70.458], Tmin=(1570.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.2447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC=CC(C)O(27858)',
    structure = SMILES('C#CC=CC(C)O'),
    E0 = (39.2295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777089,0.0592355,-2.63419e-05,-1.13829e-08,9.29804e-12,4844.51,26.443], Tmin=(100,'K'), Tmax=(1004.4,'K')), NASAPolynomial(coeffs=[16.3906,0.0221349,-8.39009e-06,1.5623e-09,-1.12017e-13,442.993,-55.2496], Tmin=(1004.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.2295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[CH]CC[CH]O(26522)',
    structure = SMILES('C#C[CH]CC[CH]O'),
    E0 = (268.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315033,0.082964,-0.000107124,7.72319e-08,-2.20461e-11,32448.2,27.7267], Tmin=(100,'K'), Tmax=(934.417,'K')), NASAPolynomial(coeffs=[11.7595,0.0283216,-1.03358e-05,1.70495e-09,-1.07451e-13,30556.1,-25.3898], Tmin=(934.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
)

species(
    label = 'C#CC([CH2])C([CH2])O(26527)',
    structure = SMILES('C#CC([CH2])C([CH2])O'),
    E0 = (337.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,283.786],'cm^-1')),
        HinderedRotor(inertia=(0.00209405,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00209922,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147656,'amu*angstrom^2'), symmetry=1, barrier=(8.42456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32971,'amu*angstrom^2'), symmetry=1, barrier=(75.9643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3299,'amu*angstrom^2'), symmetry=1, barrier=(75.9629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3829.04,'J/mol'), sigma=(6.57805,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=598.09 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.508899,0.0748716,-7.55395e-05,3.64921e-08,-5.00133e-12,40671.1,29.4412], Tmin=(100,'K'), Tmax=(837.637,'K')), NASAPolynomial(coeffs=[13.9899,0.0236927,-7.52394e-06,1.16864e-09,-7.25831e-14,37949.7,-35.972], Tmin=(837.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = 'C#CC1CC(O)C1(26530)',
    structure = SMILES('C#CC1CC(O)C1'),
    E0 = (78.6025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47129,0.0389168,3.27457e-05,-7.40327e-08,3.25117e-11,9560.02,22.4385], Tmin=(100,'K'), Tmax=(940.673,'K')), NASAPolynomial(coeffs=[15.5327,0.0211964,-6.08629e-06,1.03482e-09,-7.5117e-14,5053.14,-54.4339], Tmin=(940.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.6025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
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
    label = 'C#C[CH]C[CH]O(26949)',
    structure = SMILES('C#C[CH]C[CH]O'),
    E0 = (292.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,303.289,305.12],'cm^-1')),
        HinderedRotor(inertia=(0.00388476,'amu*angstrom^2'), symmetry=1, barrier=(10.6425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161691,'amu*angstrom^2'), symmetry=1, barrier=(10.6485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05523,'amu*angstrom^2'), symmetry=1, barrier=(68.8377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0553,'amu*angstrom^2'), symmetry=1, barrier=(68.8341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99162,0.0678312,-9.20601e-05,6.75543e-08,-1.91876e-11,35284.2,23.0522], Tmin=(100,'K'), Tmax=(981.921,'K')), NASAPolynomial(coeffs=[10.6377,0.020251,-6.71881e-06,1.01928e-09,-5.96874e-14,33789.3,-21.2748], Tmin=(981.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
)

species(
    label = 'C3H2(81)',
    structure = SMILES('[CH]C#C'),
    E0 = (530.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.621997,'amu*angstrom^2'), symmetry=1, barrier=(14.3009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6874,0.0353211,-7.31245e-05,7.03801e-08,-2.45313e-11,63833.3,8.93418], Tmin=(100,'K'), Tmax=(908.454,'K')), NASAPolynomial(coeffs=[4.78002,0.0100611,-4.92178e-06,8.868e-10,-5.66859e-14,64115.2,2.68365], Tmin=(908.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCCH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=[C]C1CC(O)C1(27859)',
    structure = SMILES('[CH]=[C]C1CC(O)C1'),
    E0 = (398.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33726,0.0438712,1.45536e-05,-5.10882e-08,2.28118e-11,48048.4,25.9885], Tmin=(100,'K'), Tmax=(977.127,'K')), NASAPolynomial(coeffs=[14.9018,0.0236144,-8.49547e-06,1.57952e-09,-1.15009e-13,43713.7,-47.7547], Tmin=(977.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(O)C[C]=C=C(27860)',
    structure = SMILES('[CH2]C#CCC([CH2])O'),
    E0 = (279.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.736995,0.0689778,-6.58087e-05,3.465e-08,-7.41136e-12,33693.4,29.3461], Tmin=(100,'K'), Tmax=(1126.54,'K')), NASAPolynomial(coeffs=[12.4809,0.0272782,-1.02848e-05,1.79158e-09,-1.19414e-13,31047.4,-28.7109], Tmin=(1126.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(CJCO)"""),
)

species(
    label = '[CH]=C=[C]CC(C)O(27861)',
    structure = SMILES('[CH]=C=[C]CC(C)O'),
    E0 = (320.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,232.725],'cm^-1')),
        HinderedRotor(inertia=(0.265273,'amu*angstrom^2'), symmetry=1, barrier=(10.2264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267413,'amu*angstrom^2'), symmetry=1, barrier=(10.2301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270879,'amu*angstrom^2'), symmetry=1, barrier=(10.2241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.503014,'amu*angstrom^2'), symmetry=1, barrier=(19.0081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729621,0.0724498,-7.61438e-05,4.50656e-08,-1.09018e-11,38619.5,28.3833], Tmin=(100,'K'), Tmax=(997.746,'K')), NASAPolynomial(coeffs=[11.3901,0.0297114,-1.18909e-05,2.13325e-09,-1.44378e-13,36492.2,-23.0231], Tmin=(997.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)C=C[C]=C(27803)',
    structure = SMILES('[CH2]C(O)[CH]C=C=C'),
    E0 = (256.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314026,0.0767939,-7.37798e-05,3.66188e-08,-7.26109e-12,30963.5,26.6045], Tmin=(100,'K'), Tmax=(1216.86,'K')), NASAPolynomial(coeffs=[15.9443,0.0254149,-1.04459e-05,1.92081e-09,-1.32498e-13,27159.5,-51.8703], Tmin=(1216.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2][C]=CCC(=C)O(20261)',
    structure = SMILES('[CH2][C]=CCC(=C)O'),
    E0 = (227.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,324.578,324.738],'cm^-1')),
        HinderedRotor(inertia=(0.197307,'amu*angstrom^2'), symmetry=1, barrier=(14.7404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60053,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197311,'amu*angstrom^2'), symmetry=1, barrier=(14.7429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197134,'amu*angstrom^2'), symmetry=1, barrier=(14.7402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538996,0.0646783,-3.62488e-05,-6.59434e-09,9.27722e-12,27497.6,27.8633], Tmin=(100,'K'), Tmax=(965.006,'K')), NASAPolynomial(coeffs=[17.9619,0.0189214,-6.25678e-06,1.10175e-09,-7.81432e-14,22902.9,-61.9555], Tmin=(965.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([O])CC=C=C(27862)',
    structure = SMILES('[CH2]C([O])CC=C=C'),
    E0 = (369.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,697.133,700.993],'cm^-1')),
        HinderedRotor(inertia=(0.6503,'amu*angstrom^2'), symmetry=1, barrier=(14.9517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.414086,'amu*angstrom^2'), symmetry=1, barrier=(14.8925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.649599,'amu*angstrom^2'), symmetry=1, barrier=(14.9355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70047,0.0687854,-5.99086e-05,2.74151e-08,-5.09341e-12,44593.3,27.826], Tmin=(100,'K'), Tmax=(1280,'K')), NASAPolynomial(coeffs=[13.9927,0.0272473,-1.12314e-05,2.06247e-09,-1.4176e-13,41190.4,-39.5828], Tmin=(1280,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'C=C=CCC(=C)O(20265)',
    structure = SMILES('C=C=CCC(=C)O'),
    E0 = (14.7699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613394,0.0621937,-2.89783e-05,-1.4181e-08,1.1956e-11,1909.43,25.8958], Tmin=(100,'K'), Tmax=(963.233,'K')), NASAPolynomial(coeffs=[18.053,0.0186319,-6.08307e-06,1.07763e-09,-7.7264e-14,-2789.08,-64.5367], Tmin=(963.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.7699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (287.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (356.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (384.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (465.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (321.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (420.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (488.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (390.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (466.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (444.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (405.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (442.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (625.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (473.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (487.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (595.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (513.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (346.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (376.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (365.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (445.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (497.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (295.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (674.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (662.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (398.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (419.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (364.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (331.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (320.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (470.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (312.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['CH2CHOH(42)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['[CH]=C1[CH]CC(O)C1(27847)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(69.2285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs2H] for rate rule [R6;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C#C[CH]CC(=C)O(27848)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC=CC([CH2])O(27849)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CHOH(42)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(8.71744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]O(284)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(5)', 'C#C[CH]CC=C(26789)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00674586,'m^3/(mol*s)'), n=2.3625, Ea=(84.2031,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;OJ] for rate rule [Cds-CsH_Cds-HH;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#C[CH]C[C](C)O(27850)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC[CH]C([CH2])O(27795)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]CC(C)[O](27231)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#C[CH][CH]C(C)O(27851)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CCC[C]([CH2])O(27852)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CCCC([CH2])O(27853)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CCCC([CH2])[O](27854)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#C[CH]CC(C)O(27855)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.31883e+06,'s^-1'), n=1.02765, Ea=(75.0925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]O(284)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['[CH2]C(O)CC1[C]=C1(27856)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['OC1C[CH][C]=CC1(27831)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.33254e+09,'s^-1'), n=0.487896, Ea=(59.5573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_linear;multiplebond_intra;radadd_intra_cs2H] + [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#CCCC(=C)O(27857)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#CC=CC(C)O(27858)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([CH2])C([CH2])O(26527)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#CC1CC(O)C1(26530)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(19)', 'C#C[CH]C[CH]O(26949)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H2(81)', '[CH2]C([CH2])O(2453)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.26928e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['[CH]=[C]C1CC(O)C1(27859)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(111.24,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 106.8 to 111.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['[CH2]C(O)C[C]=C=C(27860)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C=[C]CC(C)O(27861)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['[CH2]C(O)C=C[C]=C(27803)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/O)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['[CH2][C]=CCC(=C)O(20261)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([O])CC=C=C(27862)'],
    products = ['C#C[CH]CC([CH2])O(26526)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3027.78,'s^-1'), n=2.16706, Ea=(100.661,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;O_rad_out;Cd_H_out_singleH] + [R6H;O_rad_out;XH_out] for rate rule [R6H;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C=C=CCC(=C)O(20265)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

network(
    label = '4494',
    isomers = [
        'C#C[CH]CC([CH2])O(26526)',
    ],
    reactants = [
        ('CH2CHOH(42)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4494',
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

