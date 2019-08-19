species(
    label = 'C#C[CH]CC[CH]O(26522)',
    structure = SMILES('C#C[CH]CC[CH]O'),
    E0 = (268.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315033,0.082964,-0.000107124,7.72319e-08,-2.20461e-11,32448.2,27.7267], Tmin=(100,'K'), Tmax=(934.417,'K')), NASAPolynomial(coeffs=[11.7595,0.0283216,-1.03358e-05,1.70495e-09,-1.07451e-13,30556.1,-25.3898], Tmin=(934.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
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
    label = '[CH]=C1[CH]CCC1O(27779)',
    structure = SMILES('[CH]C1=CCCC1O'),
    E0 = (180.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43292,0.0375114,4.76904e-05,-8.72145e-08,3.54658e-11,21857,24.3187], Tmin=(100,'K'), Tmax=(972.343,'K')), NASAPolynomial(coeffs=[14.9541,0.0285292,-1.04043e-05,1.94868e-09,-1.42808e-13,17022.7,-51.8717], Tmin=(972.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#C[CH]CCC=O(27877)',
    structure = SMILES('C#C[CH]CCC=O'),
    E0 = (165.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,219.396,3437.59],'cm^-1')),
        HinderedRotor(inertia=(3.42184,'amu*angstrom^2'), symmetry=1, barrier=(79.6459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146839,'amu*angstrom^2'), symmetry=1, barrier=(4.62717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72261,'amu*angstrom^2'), symmetry=1, barrier=(79.4265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.98546,'amu*angstrom^2'), symmetry=1, barrier=(79.4829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28214,0.0629549,-6.13837e-05,3.13522e-08,-4.34165e-12,19992,24.0587], Tmin=(100,'K'), Tmax=(699.137,'K')), NASAPolynomial(coeffs=[8.22901,0.0317132,-1.25989e-05,2.23037e-09,-1.49112e-13,18812.8,-8.45592], Tmin=(699.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC=CC[CH]O(27878)',
    structure = SMILES('C#CC=CC[CH]O'),
    E0 = (232.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,232.731,232.732],'cm^-1')),
        HinderedRotor(inertia=(0.501975,'amu*angstrom^2'), symmetry=1, barrier=(19.2938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501965,'amu*angstrom^2'), symmetry=1, barrier=(19.2939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502001,'amu*angstrom^2'), symmetry=1, barrier=(19.2939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502024,'amu*angstrom^2'), symmetry=1, barrier=(19.2938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700799,0.0705183,-7.07593e-05,3.72297e-08,-7.85576e-12,28084.6,25.7036], Tmin=(100,'K'), Tmax=(1144.62,'K')), NASAPolynomial(coeffs=[13.9663,0.0241607,-1.00086e-05,1.84638e-09,-1.27588e-13,25047.8,-40.0864], Tmin=(1144.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CCsJOH)"""),
)

species(
    label = 'C#C[CH]CC=CO(27879)',
    structure = SMILES('C#C[CH]CC=CO'),
    E0 = (167.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16227,'amu*angstrom^2'), symmetry=1, barrier=(26.7229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16369,'amu*angstrom^2'), symmetry=1, barrier=(26.7555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.164,'amu*angstrom^2'), symmetry=1, barrier=(26.7627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16328,'amu*angstrom^2'), symmetry=1, barrier=(26.7462,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.72845,0.0776919,-7.80043e-05,3.82441e-08,-6.88619e-12,20305,29.5664], Tmin=(100,'K'), Tmax=(1640.16,'K')), NASAPolynomial(coeffs=[19.5509,0.0100547,4.79555e-07,-4.14924e-10,3.64806e-14,16098.1,-70.8491], Tmin=(1640.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
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
    label = 'C#C[CH]CCC[O](27880)',
    structure = SMILES('C#C[CH]CCC[O]'),
    E0 = (314.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2175,525,750,770,3400,2100,3025,407.5,1350,352.5,180,180,1822.37,1822.41],'cm^-1')),
        HinderedRotor(inertia=(0.0616152,'amu*angstrom^2'), symmetry=1, barrier=(1.41666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221805,'amu*angstrom^2'), symmetry=1, barrier=(5.09974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84913,'amu*angstrom^2'), symmetry=1, barrier=(65.507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84568,'amu*angstrom^2'), symmetry=1, barrier=(65.4279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.963477,0.0712485,-8.72816e-05,6.86596e-08,-2.23209e-11,37883.6,26.7128], Tmin=(100,'K'), Tmax=(863.44,'K')), NASAPolynomial(coeffs=[6.44804,0.0379443,-1.57067e-05,2.80484e-09,-1.86678e-13,37230.9,2.76277], Tmin=(863.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCOJ)"""),
)

species(
    label = 'C#CC[CH]C[CH]O(27876)',
    structure = SMILES('C#CC[CH]C[CH]O'),
    E0 = (316.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.609708,0.0821102,-0.000120755,1.0313e-07,-3.44452e-11,38235.1,29.3958], Tmin=(100,'K'), Tmax=(888.463,'K')), NASAPolynomial(coeffs=[6.95242,0.0366634,-1.55099e-05,2.75991e-09,-1.8149e-13,37774.8,3.29771], Tmin=(888.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOH) + radical(RCCJCC)"""),
)

species(
    label = 'C#C[CH]C[CH]CO(27881)',
    structure = SMILES('C#C[CH]C[CH]CO'),
    E0 = (288.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785696,0.0688944,-6.67094e-05,3.44198e-08,-6.43695e-12,34792.6,28.9872], Tmin=(100,'K'), Tmax=(879.195,'K')), NASAPolynomial(coeffs=[11.7647,0.0271249,-9.40354e-06,1.54994e-09,-9.98127e-14,32545.9,-24.3649], Tmin=(879.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCJCO)"""),
)

species(
    label = 'C#CCC[CH][CH]O(27882)',
    structure = SMILES('C#CCC[CH][CH]O'),
    E0 = (322.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568897,0.07617,-8.71368e-05,5.58555e-08,-1.44488e-11,38897.9,29.4773], Tmin=(100,'K'), Tmax=(941.658,'K')), NASAPolynomial(coeffs=[11.8297,0.0283357,-1.09392e-05,1.9093e-09,-1.26588e-13,36777.2,-24.1725], Tmin=(941.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = 'C#C[CH][CH]CCO(27883)',
    structure = SMILES('[CH]=C=C[CH]CCO'),
    E0 = (236.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8645,0.0660232,-5.77993e-05,2.79248e-08,-5.54905e-12,28504.2,26.6593], Tmin=(100,'K'), Tmax=(1197.32,'K')), NASAPolynomial(coeffs=[11.9932,0.0288446,-1.12223e-05,1.99092e-09,-1.34081e-13,25839.3,-29.0345], Tmin=(1197.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(C=C=CJ)"""),
)

species(
    label = '[C]#CCCC[CH]O(27884)',
    structure = SMILES('[C]#CCCC[CH]O'),
    E0 = (459.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2175,525,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.340336,0.0856411,-0.000121043,9.60467e-08,-2.98568e-11,55408.2,28.5405], Tmin=(100,'K'), Tmax=(913.515,'K')), NASAPolynomial(coeffs=[9.84588,0.0317676,-1.24645e-05,2.12685e-09,-1.36145e-13,54182.8,-13.6604], Tmin=(913.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CCsJOH)"""),
)

species(
    label = 'C#CCCC[CH][O](27885)',
    structure = SMILES('C#CCCC[CH][O]'),
    E0 = (348.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2175,525,750,770,3400,2100,3025,407.5,1350,352.5,247.658,247.665,2405.02,2405.02],'cm^-1')),
        HinderedRotor(inertia=(1.41192,'amu*angstrom^2'), symmetry=1, barrier=(61.4519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255937,'amu*angstrom^2'), symmetry=1, barrier=(11.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255936,'amu*angstrom^2'), symmetry=1, barrier=(11.1399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41184,'amu*angstrom^2'), symmetry=1, barrier=(61.4519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71339,0.0789426,-0.000109312,9.24154e-08,-3.14525e-11,41990.3,27.3208], Tmin=(100,'K'), Tmax=(851.612,'K')), NASAPolynomial(coeffs=[6.54824,0.0390935,-1.72062e-05,3.15555e-09,-2.12731e-13,41447.7,2.7574], Tmin=(851.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[C]#C[CH]CCCO(27886)',
    structure = SMILES('[C]#C[CH]CCCO'),
    E0 = (425.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2175,525,3615,1277.5,1000,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564297,0.0782895,-0.000100397,7.44147e-08,-2.18088e-11,51302.7,28.0241], Tmin=(100,'K'), Tmax=(953.061,'K')), NASAPolynomial(coeffs=[9.73843,0.0306273,-1.09688e-05,1.77678e-09,-1.10131e-13,49969.9,-13.6126], Tmin=(953.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'O[CH]CCC1[C]=C1(27887)',
    structure = SMILES('O[CH]CCC1[C]=C1'),
    E0 = (442.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.739066,0.0788424,-9.16962e-05,5.11072e-08,-6.01688e-12,53296.4,24.6166], Tmin=(100,'K'), Tmax=(621.025,'K')), NASAPolynomial(coeffs=[9.69799,0.0341703,-1.52739e-05,2.85872e-09,-1.96643e-13,51932.3,-16.3604], Tmin=(621.025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(cyclopropenyl-vinyl) + radical(CCsJOH)"""),
)

species(
    label = 'OC1C=[C][CH]CC1(27810)',
    structure = SMILES('OC1[CH][C]=CCC1'),
    E0 = (161.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39497,0.039938,3.25428e-05,-7.04147e-08,2.93764e-11,19475.4,19.8979], Tmin=(100,'K'), Tmax=(982.149,'K')), NASAPolynomial(coeffs=[15.348,0.0252758,-9.46042e-06,1.80757e-09,-1.33721e-13,14701,-57.5191], Tmin=(982.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C#CC=CCCO(27888)',
    structure = SMILES('C#CC=CCCO'),
    E0 = (52.2087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596713,0.0660543,-5.388e-05,2.25876e-08,-3.79044e-12,6409.08,26.5356], Tmin=(100,'K'), Tmax=(1425.15,'K')), NASAPolynomial(coeffs=[15.7931,0.0234012,-8.98573e-06,1.58614e-09,-1.06278e-13,2077.76,-52.1615], Tmin=(1425.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCC=CO(27889)',
    structure = SMILES('C#CCCC=CO'),
    E0 = (21.0107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.511666,0.0624063,-1.95078e-05,-3.40483e-08,2.24785e-11,2666,25.2653], Tmin=(100,'K'), Tmax=(905.939,'K')), NASAPolynomial(coeffs=[20.5665,0.0127115,-1.55754e-06,8.29573e-11,-4.25595e-15,-2562.1,-78.3065], Tmin=(905.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.0107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCCC=O(27890)',
    structure = SMILES('C#CCCCC=O'),
    E0 = (19.2217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25835,0.0626665,-5.63014e-05,3.08809e-08,-7.38124e-12,2408.55,24.2911], Tmin=(100,'K'), Tmax=(976.655,'K')), NASAPolynomial(coeffs=[7.7987,0.0358794,-1.51597e-05,2.7972e-09,-1.92382e-13,1131.04,-7.10774], Tmin=(976.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.2217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC([CH2])C[CH]O(26523)',
    structure = SMILES('C#CC([CH2])C[CH]O'),
    E0 = (318.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,1901.59,1901.64],'cm^-1')),
        HinderedRotor(inertia=(0.456642,'amu*angstrom^2'), symmetry=1, barrier=(12.9392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454656,'amu*angstrom^2'), symmetry=1, barrier=(12.941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42,'amu*angstrom^2'), symmetry=1, barrier=(69.1444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.455212,'amu*angstrom^2'), symmetry=1, barrier=(12.9401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43737,'amu*angstrom^2'), symmetry=1, barrier=(69.1489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3806.13,'J/mol'), sigma=(6.54554,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.51 K, Pc=30.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.609239,0.077667,-8.53754e-05,4.36253e-08,-4.61321e-12,38420,28.3647], Tmin=(100,'K'), Tmax=(708.798,'K')), NASAPolynomial(coeffs=[12.1751,0.027106,-9.50391e-06,1.54207e-09,-9.68423e-14,36411,-26.0596], Tmin=(708.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOH) + radical(Isobutyl)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3806.13,'J/mol'), sigma=(6.54554,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.51 K, Pc=30.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486979,0.0764537,-8.14569e-05,4.42898e-08,-8.33599e-12,34687.6,27.8559], Tmin=(100,'K'), Tmax=(832.808,'K')), NASAPolynomial(coeffs=[13.3624,0.0253093,-8.60445e-06,1.39334e-09,-8.85159e-14,32172.2,-34.1318], Tmin=(832.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1CCC1O(26531)',
    structure = SMILES('C#CC1CCC1O'),
    E0 = (78.6025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47129,0.0389168,3.27457e-05,-7.40327e-08,3.25117e-11,9560.02,22.4385], Tmin=(100,'K'), Tmax=(940.673,'K')), NASAPolynomial(coeffs=[15.5327,0.0211964,-6.08629e-06,1.03482e-09,-7.5117e-14,5053.14,-54.4339], Tmin=(940.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.6025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
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
    label = '[CH2]C[CH]O(562)',
    structure = SMILES('[CH2]C[CH]O'),
    E0 = (112.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,4000],'cm^-1')),
        HinderedRotor(inertia=(1.33788,'amu*angstrom^2'), symmetry=1, barrier=(30.7606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145425,'amu*angstrom^2'), symmetry=1, barrier=(3.51086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15145,'amu*angstrom^2'), symmetry=1, barrier=(3.48214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4465,0.0325713,-2.2601e-05,8.44361e-09,-1.31566e-12,13607.8,17.8684], Tmin=(100,'K'), Tmax=(1466.6,'K')), NASAPolynomial(coeffs=[8.1555,0.0170006,-6.67571e-06,1.20454e-09,-8.16704e-14,11933.3,-11.8605], Tmin=(1466.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C=CC[CH2](26618)',
    structure = SMILES('C#C[CH]C[CH2]'),
    E0 = (477.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0141456,'amu*angstrom^2'), symmetry=1, barrier=(6.90141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.94797,'amu*angstrom^2'), symmetry=1, barrier=(67.7796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0405762,'amu*angstrom^2'), symmetry=1, barrier=(67.6916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76795,0.0430915,-3.64272e-05,1.77221e-08,-3.47093e-12,57486.3,19.7996], Tmin=(100,'K'), Tmax=(1346.56,'K')), NASAPolynomial(coeffs=[9.46877,0.0175599,-5.02768e-06,7.11764e-10,-4.08662e-14,55653.2,-18.7496], Tmin=(1346.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=[C]C1CCC1O(27891)',
    structure = SMILES('[CH]=[C]C1CCC1O'),
    E0 = (398.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33726,0.0438712,1.45536e-05,-5.10882e-08,2.28118e-11,48048.4,25.9885], Tmin=(100,'K'), Tmax=(977.127,'K')), NASAPolynomial(coeffs=[14.9018,0.0236144,-8.49547e-06,1.57952e-09,-1.15009e-13,43713.7,-47.7547], Tmin=(977.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C=[C]CC[CH]O(27892)',
    structure = SMILES('[CH2]C#CCC[CH]O'),
    E0 = (260.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784352,0.0727062,-8.08166e-05,5.24009e-08,-1.40048e-11,31444.7,28.4432], Tmin=(100,'K'), Tmax=(904.261,'K')), NASAPolynomial(coeffs=[10.0436,0.0317472,-1.28723e-05,2.30831e-09,-1.55545e-13,29770.1,-15.2956], Tmin=(904.261,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CCCO(27893)',
    structure = SMILES('[CH]=C=[C]CCCO'),
    E0 = (332.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,540,610,2055,3120,650,792.5,1650,3615,1277.5,1000,1685,370,274.342,274.551],'cm^-1')),
        HinderedRotor(inertia=(0.152564,'amu*angstrom^2'), symmetry=1, barrier=(8.11102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152453,'amu*angstrom^2'), symmetry=1, barrier=(8.11414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151255,'amu*angstrom^2'), symmetry=1, barrier=(8.10567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298433,'amu*angstrom^2'), symmetry=1, barrier=(15.8617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.740406,0.0744948,-8.64657e-05,5.87748e-08,-1.64661e-11,40137.5,28.6921], Tmin=(100,'K'), Tmax=(864.068,'K')), NASAPolynomial(coeffs=[9.82508,0.032437,-1.34502e-05,2.43683e-09,-1.64978e-13,38567.7,-13.8084], Tmin=(864.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CC[CH]O(27825)',
    structure = SMILES('C=[C]C=CC[CH]O'),
    E0 = (254.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,193.955,193.989,193.999],'cm^-1')),
        HinderedRotor(inertia=(0.522666,'amu*angstrom^2'), symmetry=1, barrier=(13.9536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522574,'amu*angstrom^2'), symmetry=1, barrier=(13.9535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.47971,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522428,'amu*angstrom^2'), symmetry=1, barrier=(13.9534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430935,0.077993,-8.62319e-05,5.20268e-08,-1.25918e-11,30764.2,27.3759], Tmin=(100,'K'), Tmax=(1006.21,'K')), NASAPolynomial(coeffs=[13.2122,0.0271836,-1.04881e-05,1.84252e-09,-1.23158e-13,28192.1,-34.3653], Tmin=(1006.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCsJOH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CCC=CO(23619)',
    structure = SMILES('[CH2][C]=CCC=CO'),
    E0 = (233.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.894772,'amu*angstrom^2'), symmetry=1, barrier=(20.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895795,'amu*angstrom^2'), symmetry=1, barrier=(20.5961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894823,'amu*angstrom^2'), symmetry=1, barrier=(20.5737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895314,'amu*angstrom^2'), symmetry=1, barrier=(20.585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.50828,0.0616464,-1.79746e-05,-3.27952e-08,2.05646e-11,28195.8,27.7757], Tmin=(100,'K'), Tmax=(935.225,'K')), NASAPolynomial(coeffs=[20.5413,0.014238,-3.32393e-06,5.2059e-10,-3.88035e-14,22775,-76.4793], Tmin=(935.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CCCC=O(20129)',
    structure = SMILES('[CH2][C]=CCCC=O'),
    E0 = (230.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.643459,'amu*angstrom^2'), symmetry=1, barrier=(14.7944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00282082,'amu*angstrom^2'), symmetry=1, barrier=(14.7926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0238838,'amu*angstrom^2'), symmetry=1, barrier=(0.549137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.73376,'amu*angstrom^2'), symmetry=1, barrier=(39.8625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49273,0.0592251,-4.54673e-05,1.99228e-08,-3.89805e-12,27775.1,25.92], Tmin=(100,'K'), Tmax=(1134.32,'K')), NASAPolynomial(coeffs=[7.56511,0.0378119,-1.71511e-05,3.28073e-09,-2.30222e-13,26397.5,-4.14108], Tmin=(1134.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CCC=CO(23622)',
    structure = SMILES('C=C=CCC=CO'),
    E0 = (20.5358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580229,0.0591917,-1.08129e-05,-4.02349e-08,2.31783e-11,2607.74,25.8168], Tmin=(100,'K'), Tmax=(936.018,'K')), NASAPolynomial(coeffs=[20.6413,0.0139334,-3.14137e-06,4.94372e-10,-3.77492e-14,-2920.67,-79.1103], Tmin=(936.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.5358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=CCCC=O(20134)',
    structure = SMILES('C=C=CCCC=O'),
    E0 = (17.4773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55019,0.0570762,-3.98606e-05,1.49463e-08,-2.44505e-12,2187.17,24.0022], Tmin=(100,'K'), Tmax=(1332.5,'K')), NASAPolynomial(coeffs=[8.59192,0.0359379,-1.60654e-05,3.04139e-09,-2.11499e-13,310.535,-11.9915], Tmin=(1332.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.4773,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    E0 = (268.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (353.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (413.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (447.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (385.749,'kJ/mol'),
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
    E0 = (344.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (403.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (442.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (397.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (446.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (411.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (606.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (409.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (485.713,'kJ/mol'),
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
    E0 = (494.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (295.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (332.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (368.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (307.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (478.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (445.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (276.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (643.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (683.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (398.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (400.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (377.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (313.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (301.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (420.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (303.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (293.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['CH2CHOH(42)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['[CH]=C1[CH]CCC1O(27779)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs] for rate rule [R6;triplebond_intra_H;radadd_intra_csHNd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C#C[CH]CCC=O(27877)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC=CC[CH]O(27878)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#C[CH]CC=CO(27879)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]O(284)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CHOH(42)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00668046,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C#C[CH]CCC[O](27880)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(153000,'s^-1'), n=2.26, Ea=(88.9937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 325 used for R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC[CH]C[CH]O(27876)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]C[CH]CO(27881)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.4e-20,'s^-1'), n=9.13, Ea=(108.784,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 341 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CCC[CH][CH]O(27882)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C#C[CH][CH]CCO(27883)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CCCC[CH]O(27884)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CCCC[CH][O](27885)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(60863,'s^-1'), n=1.86284, Ea=(61.1759,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;O_rad_out;Cs_H_out_H/Ct]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#C[CH]CCCO(27886)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1098,'s^-1'), n=2.21, Ea=(60.1659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R6HJ_2;Ct_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]O(284)', '[CH]=[C]C=C(4699)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['O[CH]CCC1[C]=C1(27887)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['OC1C=[C][CH]CC1(27810)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHNd] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_csHO]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C#CC=CCCO(27888)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C#CCCC=CO(27889)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C#CCCCC=O(27890)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([CH2])C[CH]O(26523)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC([CH2])O(26526)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C#CC1CCC1O(26531)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Cpri_rad_out_single] + [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/OneDe;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H2(81)', '[CH2]C[CH]O(562)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['HCOH(T)(285)', '[CH]=C=CC[CH2](26618)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['[CH]=[C]C1CCC1O(27891)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.73e+06,'s^-1'), n=1.31, Ea=(129.889,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 124.8 to 129.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C=C=[C]CC[CH]O(27892)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C=[C]CCCO(27893)'],
    products = ['C#C[CH]CC[CH]O(26522)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS;Cd_rad_out_double;Cs_H_out_H/NonDeO]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C=[C]C=CC[CH]O(27825)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['[CH2][C]=CCC=CO(23619)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['[CH2][C]=CCCC=O(20129)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_5;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C=C=CCC=CO(23622)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC[CH]O(26522)'],
    products = ['C=C=CCCC=O(20134)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

network(
    label = '4490',
    isomers = [
        'C#C[CH]CC[CH]O(26522)',
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
    label = '4490',
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

