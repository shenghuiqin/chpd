species(
    label = '[CH]=C(C=C)C([O])CCC=O(29718)',
    structure = SMILES('[CH]=C(C=C)C([O])CCC=O'),
    E0 = (207.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.87783,0.102965,-9.68611e-05,4.72441e-08,-9.26957e-12,25192.9,39.0173], Tmin=(100,'K'), Tmax=(1222.51,'K')), NASAPolynomial(coeffs=[19.3753,0.0366964,-1.555e-05,2.90241e-09,-2.01711e-13,20241,-62.7613], Tmin=(1222.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
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
    label = '[CH2]C1C=C1C([O])CCC=O(30610)',
    structure = SMILES('[CH2]C1C=C1C([O])CCC=O'),
    E0 = (311.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.543741,0.102669,-0.000104912,5.96941e-08,-1.40354e-11,37645.1,37.0153], Tmin=(100,'K'), Tmax=(1016.69,'K')), NASAPolynomial(coeffs=[14.1817,0.0447348,-1.94366e-05,3.64653e-09,-2.53563e-13,34650.9,-34.2704], Tmin=(1016.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CC(C)OJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1C([CH2])OC1CCC=O(30402)',
    structure = SMILES('[CH]=C1C([CH2])OC1CCC=O'),
    E0 = (252.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.857902,0.0929358,-6.86392e-05,2.18941e-08,-1.84268e-12,30516,39.2577], Tmin=(100,'K'), Tmax=(1186.68,'K')), NASAPolynomial(coeffs=[21.4491,0.0337716,-1.41127e-05,2.64284e-09,-1.84692e-13,24093.3,-76.9335], Tmin=(1186.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CJC(C)OC) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C1CCC([O])O1(30611)',
    structure = SMILES('[CH]=C(C=C)C1CCC([O])O1'),
    E0 = (191.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00637309,0.0694304,4.86143e-06,-6.2106e-08,3.12796e-11,23187.7,32.5382], Tmin=(100,'K'), Tmax=(933.093,'K')), NASAPolynomial(coeffs=[19.6872,0.0316516,-9.30189e-06,1.52331e-09,-1.05121e-13,17486.7,-71.9153], Tmin=(933.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1=CC([O])CCC1[O](30545)',
    structure = SMILES('C=CC1=CC([O])CCC1[O]'),
    E0 = (144.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.155591,0.0635196,3.89888e-05,-1.0471e-07,4.65046e-11,17499.5,32.3189], Tmin=(100,'K'), Tmax=(964.591,'K')), NASAPolynomial(coeffs=[25.5204,0.025515,-8.38546e-06,1.62029e-09,-1.26234e-13,9360.87,-107.138], Tmin=(964.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclohexene) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = '[CH]=C(C=C)C(=O)CCC=O(30612)',
    structure = SMILES('[CH]=C(C=C)C(=O)CCC=O'),
    E0 = (49.9124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.156,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.33261,0.102404,-0.000114438,7.26081e-08,-1.93074e-11,6153.02,33.7926], Tmin=(100,'K'), Tmax=(897.474,'K')), NASAPolynomial(coeffs=[11.9887,0.0474875,-2.26515e-05,4.42575e-09,-3.14248e-13,3941.45,-24.3175], Tmin=(897.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.9124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'CH2CH2CHO(560)',
    structure = SMILES('[CH2]CC=O'),
    E0 = (11.2619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.221237,'amu*angstrom^2'), symmetry=1, barrier=(5.08666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221178,'amu*angstrom^2'), symmetry=1, barrier=(5.08532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76345,0.0293577,-2.47892e-05,1.49239e-08,-4.38497e-12,1397.08,14.4322], Tmin=(100,'K'), Tmax=(767.858,'K')), NASAPolynomial(coeffs=[4.24224,0.0216537,-9.73869e-06,1.85596e-09,-1.3004e-13,1169.99,7.6886], Tmin=(767.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.2619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH2CH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C(C=C)C=O(26939)',
    structure = SMILES('[CH]=C(C=C)C=O'),
    E0 = (227.064,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.06459,'amu*angstrom^2'), symmetry=1, barrier=(24.4769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06108,'amu*angstrom^2'), symmetry=1, barrier=(24.3963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30551,0.0499053,-3.22785e-05,9.60926e-10,4.21178e-12,27414.7,19.3092], Tmin=(100,'K'), Tmax=(1029.15,'K')), NASAPolynomial(coeffs=[15.6535,0.0127556,-5.26642e-06,1.03998e-09,-7.72081e-14,23475.5,-55.1139], Tmin=(1029.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
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
    label = 'C2H3(30)',
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
    label = 'C#CC([O])CCC=O(21967)',
    structure = SMILES('C#CC([O])CCC=O'),
    E0 = (76.6225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2175,525,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,949.418,949.547],'cm^-1')),
        HinderedRotor(inertia=(0.575444,'amu*angstrom^2'), symmetry=1, barrier=(13.2306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0876365,'amu*angstrom^2'), symmetry=1, barrier=(13.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575373,'amu*angstrom^2'), symmetry=1, barrier=(13.229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.25958,'amu*angstrom^2'), symmetry=1, barrier=(74.9442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64548,0.0736264,-7.22051e-05,3.86986e-08,-8.49089e-12,9336.52,29.3349], Tmin=(100,'K'), Tmax=(1091.24,'K')), NASAPolynomial(coeffs=[12.4568,0.0303312,-1.2692e-05,2.34033e-09,-1.61292e-13,6758.72,-28.6792], Tmin=(1091.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.6225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(C=C)[C](O)CCC=O(30613)',
    structure = SMILES('[CH]C(C=C)=C(O)CCC=O'),
    E0 = (47.2058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5843,0.11196,-0.000104668,5.05787e-08,-9.71385e-12,5887.72,38.8036], Tmin=(100,'K'), Tmax=(1261.71,'K')), NASAPolynomial(coeffs=[22.6335,0.0351824,-1.33906e-05,2.34947e-09,-1.57586e-13,-223.482,-83.6633], Tmin=(1261.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.2058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)C(O)[CH]CC=O(30614)',
    structure = SMILES('[CH]=C(C=C)C(O)[CH]CC=O'),
    E0 = (177.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.956108,0.103044,-9.87042e-05,4.87652e-08,-9.61692e-12,21533.8,41.9865], Tmin=(100,'K'), Tmax=(1223.71,'K')), NASAPolynomial(coeffs=[20.2142,0.0338438,-1.38814e-05,2.55498e-09,-1.7641e-13,16352.4,-64.4223], Tmin=(1223.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(=C)[C]([O])CCC=O(30615)',
    structure = SMILES('[CH2]C(C=C)=C([O])CCC=O'),
    E0 = (-34.1747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10022,0.104207,-9.94245e-05,4.91911e-08,-9.68448e-12,-3919.94,37.7373], Tmin=(100,'K'), Tmax=(1230.23,'K')), NASAPolynomial(coeffs=[20.7242,0.0332462,-1.29032e-05,2.30471e-09,-1.56473e-13,-9289.72,-72.0746], Tmin=(1230.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(=C)C([O])CCC=O(30616)',
    structure = SMILES('C=[C]C(=C)C([O])CCC=O'),
    E0 = (159.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.754229,0.104441,-0.000103806,5.56368e-08,-1.21526e-11,19399.9,37.771], Tmin=(100,'K'), Tmax=(1097.69,'K')), NASAPolynomial(coeffs=[16.5585,0.0413521,-1.75932e-05,3.27577e-09,-2.27125e-13,15599.1,-47.3664], Tmin=(1097.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(C=C)C(O)C[CH]C=O(30617)',
    structure = SMILES('[CH]=C(C=C)C(O)CC=C[O]'),
    E0 = (120.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.3484,0.0993633,-6.1257e-05,-7.87089e-09,1.46069e-11,14746,40.607], Tmin=(100,'K'), Tmax=(956.852,'K')), NASAPolynomial(coeffs=[27.7195,0.0210951,-6.35597e-06,1.113e-09,-8.15747e-14,7203.48,-108.692], Tmin=(956.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = 'C=CC(=C)C([O])[CH]CC=O(30618)',
    structure = SMILES('C=CC(=C)C([O])[CH]CC=O'),
    E0 = (160.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.975585,0.098433,-8.57465e-05,3.79779e-08,-6.69565e-12,19525.8,41.4681], Tmin=(100,'K'), Tmax=(1363.73,'K')), NASAPolynomial(coeffs=[21.637,0.0321074,-1.27934e-05,2.31433e-09,-1.57775e-13,13358.3,-74.6393], Tmin=(1363.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C([C]=C)C(O)CCC=O(30619)',
    structure = SMILES('[CH]C(=C=C)C(O)CCC=O'),
    E0 = (152.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.524635,0.104749,-0.00010003,5.44636e-08,-1.26344e-11,18502.5,38.8492], Tmin=(100,'K'), Tmax=(1011.21,'K')), NASAPolynomial(coeffs=[12.3431,0.0538507,-2.45301e-05,4.68967e-09,-3.29219e-13,15900,-23.3738], Tmin=(1011.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CC(=C)C([O])CCC=O(30620)',
    structure = SMILES('[CH]=CC(=C)C([O])CCC=O'),
    E0 = (207.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.87783,0.102965,-9.68611e-05,4.72441e-08,-9.26957e-12,25192.9,39.0173], Tmin=(100,'K'), Tmax=(1222.51,'K')), NASAPolynomial(coeffs=[19.3753,0.0366964,-1.555e-05,2.90241e-09,-2.01711e-13,20241,-62.7613], Tmin=(1222.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C=C)C(O)CC[C]=O(30621)',
    structure = SMILES('[CH]=C(C=C)C(O)CC[C]=O'),
    E0 = (137.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0636,0.109318,-0.000115688,6.39257e-08,-1.40934e-11,16730.6,40.6969], Tmin=(100,'K'), Tmax=(1100.67,'K')), NASAPolynomial(coeffs=[19.3465,0.0351441,-1.46032e-05,2.69958e-09,-1.8684e-13,12237.7,-59.7278], Tmin=(1100.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC(=C)C([O])C[CH]C=O(30622)',
    structure = SMILES('C=CC(=C)C([O])CC=C[O]'),
    E0 = (104.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.11561,0.0918909,-3.88442e-05,-3.00987e-08,2.20504e-11,12726.9,39.1764], Tmin=(100,'K'), Tmax=(958.876,'K')), NASAPolynomial(coeffs=[27.4251,0.0221174,-6.79435e-06,1.22203e-09,-9.12749e-14,4987.76,-109.132], Tmin=(958.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC(=[CH])C(O)CCC=O(30623)',
    structure = SMILES('[CH]=CC(=[CH])C(O)CCC=O'),
    E0 = (224.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.983579,0.108899,-0.000113786,6.23449e-08,-1.37192e-11,27206.6,39.9958], Tmin=(100,'K'), Tmax=(1097.87,'K')), NASAPolynomial(coeffs=[18.6409,0.0373985,-1.60963e-05,3.02437e-09,-2.11073e-13,22897.6,-56.5134], Tmin=(1097.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(=C)C([O])CC[C]=O(30624)',
    structure = SMILES('C=CC(=C)C([O])CC[C]=O'),
    E0 = (120.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.967159,0.103474,-9.89993e-05,4.90483e-08,-9.71328e-12,14717.4,39.7533], Tmin=(100,'K'), Tmax=(1218.11,'K')), NASAPolynomial(coeffs=[20.0233,0.0345447,-1.41181e-05,2.59245e-09,-1.78732e-13,9603.7,-65.6545], Tmin=(1218.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(CCC=O)C1[CH]CC=1(30625)',
    structure = SMILES('[O]C(CCC=O)C1[CH]CC=1'),
    E0 = (174.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.203758,0.0774602,-3.0018e-05,-1.35805e-08,9.55872e-12,21152.1,35.8605], Tmin=(100,'K'), Tmax=(1073.65,'K')), NASAPolynomial(coeffs=[18.2195,0.0380768,-1.58667e-05,3.01086e-09,-2.1397e-13,15510,-62.1823], Tmin=(1073.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(cyclobutene-allyl)"""),
)

species(
    label = '[CH]=C1[CH]COC1CCC=O(30596)',
    structure = SMILES('[CH]C1=CCOC1CCC=O'),
    E0 = (79.3679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.340019,0.0811186,-3.16495e-05,-1.00621e-08,7.57518e-12,9714.02,35.8957], Tmin=(100,'K'), Tmax=(1111.35,'K')), NASAPolynomial(coeffs=[16.7592,0.0457003,-1.91073e-05,3.56715e-09,-2.49175e-13,4300.02,-55.6618], Tmin=(1111.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.3679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(25dihydrofuran) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=C)C1CC[CH]OO1(30626)',
    structure = SMILES('[CH]=C(C=C)C1CC[CH]OO1'),
    E0 = (394.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.433329,0.0735877,1.5538e-05,-8.98752e-08,4.57131e-11,47633,32.1781], Tmin=(100,'K'), Tmax=(910.592,'K')), NASAPolynomial(coeffs=[25.7982,0.0216981,-3.32184e-06,3.19869e-10,-2.14972e-14,40229.8,-106.336], Tmin=(910.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(473.925,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(12dioxane) + radical(Cds_P) + radical(CCsJOOC)"""),
)

species(
    label = 'C=CC1=CO[CH]CCC1[O](30365)',
    structure = SMILES('C=CC1=CO[CH]CCC1[O]'),
    E0 = (137.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979768,0.0390274,0.000211396,-3.69639e-07,1.64919e-10,16804.3,33.8246], Tmin=(100,'K'), Tmax=(901.571,'K')), NASAPolynomial(coeffs=[58.7046,-0.0390148,3.05164e-05,-6.12255e-09,4.05499e-13,-1547.85,-290.028], Tmin=(901.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(CCsJOC(O)) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC(=C)C(=O)CCC=O(30627)',
    structure = SMILES('C=CC(=C)C(=O)CCC=O'),
    E0 = (-197.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.168441,0.0978991,-9.47053e-05,5.16418e-08,-1.19892e-11,-23570.8,32.7946], Tmin=(100,'K'), Tmax=(1008.8,'K')), NASAPolynomial(coeffs=[11.9124,0.0499975,-2.34797e-05,4.57236e-09,-3.2456e-13,-26008.2,-25.5942], Tmin=(1008.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-197.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = '[CH]=C(C=C)C([O])CC(26557)',
    structure = SMILES('[CH]=C(C=C)C([O])CC'),
    E0 = (314.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,330.797,330.835,330.996],'cm^-1')),
        HinderedRotor(inertia=(0.20041,'amu*angstrom^2'), symmetry=1, barrier=(15.5681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200262,'amu*angstrom^2'), symmetry=1, barrier=(15.5681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200222,'amu*angstrom^2'), symmetry=1, barrier=(15.5671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.200509,'amu*angstrom^2'), symmetry=1, barrier=(15.5678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3956.66,'J/mol'), sigma=(6.72118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=618.02 K, Pc=29.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183153,0.0790989,-4.8097e-05,7.17233e-10,6.72178e-12,37994.2,32.2454], Tmin=(100,'K'), Tmax=(1004.03,'K')), NASAPolynomial(coeffs=[19.5827,0.026973,-9.99208e-06,1.82298e-09,-1.28835e-13,32683.4,-69.8747], Tmin=(1004.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1=COC1CCC=O(30387)',
    structure = SMILES('C=CC1=COC1CCC=O'),
    E0 = (-117.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.842832,0.0802827,7.37923e-07,-7.4949e-08,3.88023e-11,-13951.7,33.1098], Tmin=(100,'K'), Tmax=(949.976,'K')), NASAPolynomial(coeffs=[29.1195,0.0189896,-4.90613e-06,8.90924e-10,-7.21211e-14,-22571.4,-125.309], Tmin=(949.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C(C=C)C([O])COC=C(30628)',
    structure = SMILES('[CH]=C(C=C)C([O])COC=C'),
    E0 = (249.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65917,0.103912,-6.21579e-05,-1.45516e-08,1.87704e-11,30239.1,38.8671], Tmin=(100,'K'), Tmax=(946.713,'K')), NASAPolynomial(coeffs=[30.7488,0.0171902,-4.30045e-06,7.22619e-10,-5.5536e-14,21853,-127.59], Tmin=(946.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C(C=C)C([O])CC=CO(30629)',
    structure = SMILES('[CH]=C(C=C)C([O])CC=CO'),
    E0 = (209.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59394,0.102301,-5.86095e-05,-1.93279e-08,2.10547e-11,25449.3,39.8865], Tmin=(100,'K'), Tmax=(937.021,'K')), NASAPolynomial(coeffs=[30.8941,0.0156386,-3.16034e-06,4.74003e-10,-3.72638e-14,17077,-126.923], Tmin=(937.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
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
    label = '[CH]=C([CH]CCC=O)C=C(30630)',
    structure = SMILES('[CH]C(C=C)=CCCC=O'),
    E0 = (261.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.239431,0.089217,-6.81566e-05,2.76496e-08,-4.66976e-12,31638.6,33.9157], Tmin=(100,'K'), Tmax=(1368.55,'K')), NASAPolynomial(coeffs=[15.3461,0.0436637,-1.82277e-05,3.32754e-09,-2.26712e-13,27372.7,-46.1653], Tmin=(1368.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(AllylJ2_triplet)"""),
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
    E0 = (207.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (311.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (330.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (234.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (230.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (293.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (266.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (207.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (320.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (383.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (333.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (283.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (400.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (430.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (310.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (252.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (353.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (611.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (240.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (246.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (319.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (270.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (470.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (346.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (264.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (394.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (298.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (271.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (553.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (216.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (563.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (353.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (504.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['O=CCCC=O(5767)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH2]C1C=C1C([O])CCC=O(30610)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(103.848,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 14 used for R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH
Exact match found for rate rule [R4_D_D;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C1C([CH2])OC1CCC=O(30402)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C(C=C)C1CCC([O])O1(30611)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(132522,'s^-1'), n=1.48406, Ea=(26.3859,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC1=CC([O])CCC1[O](30545)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(24134.6,'s^-1'), n=1.34484, Ea=(22.9508,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_RSSR;multiplebond_intra;radadd_intra] for rate rule [R7_DSSS_CO;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 4.89897948557
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]=C(C=C)C(=O)CCC=O(30612)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2826 used for CO-CdCs_O;HJ
Exact match found for rate rule [CO-CdCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CH2CHO(560)', '[CH]=C(C=C)C=O(26939)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdH_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C(4699)', 'O=CCCC=O(5767)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0403742,'m^3/(mol*s)'), n=2.2105, Ea=(66.2944,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 62.2 to 66.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][CH]CCC=O(5764)', 'CH2CHCCH(26391)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0345031,'m^3/(mol*s)'), n=2.33733, Ea=(26.7584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-Cd_Ct-H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C2H3(30)', 'C#CC([O])CCC=O(21967)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(94600,'cm^3/(mol*s)'), n=2.41, Ea=(20.7945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2340 used for Ct-Cs_Ct-H;CdsJ-H
Exact match found for rate rule [Ct-Cs_Ct-H;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C(C=C)[C](O)CCC=O(30613)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=C)C(O)[CH]CC=O(30614)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC(=C)[C]([O])CCC=O(30615)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=[C]C(=C)C([O])CCC=O(30616)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C(C=C)C(O)C[CH]C=O(30617)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(363473,'s^-1'), n=1.92229, Ea=(102.145,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC(=C)C([O])[CH]CC=O(30618)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C([C]=C)C(O)CCC=O(30619)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=CC(=C)C([O])CCC=O(30620)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.27529e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;Cd_H_out_singleH] for rate rule [R4H_DSD;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C(C=C)C(O)CC[C]=O(30621)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(214655,'s^-1'), n=1.70206, Ea=(32.737,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;O_rad_out;XH_out] for rate rule [R5H_CCC_O;O_rad_out;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC(=C)C([O])C[CH]C=O(30622)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(13813.5,'s^-1'), n=1.88327, Ea=(38.7799,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/CO]
Euclidian distance = 3.74165738677
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=CC(=[CH])C(O)CCC=O(30623)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC(=C)C([O])CC[C]=O(30624)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(13085.7,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]C=C(4699)', '[O][CH]CCC=O(5764)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[O]C(CCC=O)C1[CH]CC=1(30625)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D_D;doublebond_intra_pri;radadd_intra_cdsingleH] for rate rule [R4_D_D;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C1[CH]COC1CCC=O(30596)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['[CH]=C(C=C)C1CC[CH]OO1(30626)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(463580,'s^-1'), n=1.14062, Ea=(186.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 183.0 to 186.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC1=CO[CH]CCC1[O](30365)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.41957e+11,'s^-1'), n=0.2055, Ea=(90.797,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cdsingleH] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC(=C)C(=O)CCC=O(30627)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CO(12)', '[CH]=C(C=C)C([O])CC(26557)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(274200,'cm^3/(mol*s)'), n=2.53, Ea=(357.732,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 4 used for CO;C_pri/NonDeC
Exact match found for rate rule [CO;C_pri/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    products = ['C=CC1=COC1CCC=O(30387)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C(C=C)C([O])COC=C(30628)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 7 used for R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C
Exact match found for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(C=C)C([O])CC=CO(30629)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_H;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction33',
    reactants = ['O(4)', '[CH]=C([CH]CCC=O)C=C(30630)'],
    products = ['[CH]=C(C=C)C([O])CCC=O(29718)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '4737',
    isomers = [
        '[CH]=C(C=C)C([O])CCC=O(29718)',
    ],
    reactants = [
        ('O=CCCC=O(5767)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4737',
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

