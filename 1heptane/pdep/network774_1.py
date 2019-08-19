species(
    label = 'CC([O])O[CH]C=O(1190)',
    structure = SMILES('CC([O])OC=C[O]'),
    E0 = (-227.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,246.728,246.765,246.797,246.827,246.862],'cm^-1')),
        HinderedRotor(inertia=(0.537772,'amu*angstrom^2'), symmetry=1, barrier=(23.2474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537774,'amu*angstrom^2'), symmetry=1, barrier=(23.2476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537755,'amu*angstrom^2'), symmetry=1, barrier=(23.2485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.602458,0.0587073,-1.79627e-05,-3.27306e-08,2.04842e-11,-27208.5,26.1681], Tmin=(100,'K'), Tmax=(944.968,'K')), NASAPolynomial(coeffs=[21.8024,0.00871087,-1.68463e-06,2.90518e-10,-2.61427e-14,-32989.5,-84.2978], Tmin=(944.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-227.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCOJ)"""),
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
    label = 'CH3CHO(37)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CC1OC([CH][O])O1(4425)',
    structure = SMILES('CC1OC([CH][O])O1'),
    E0 = (-78.8963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73423,0.0498958,-3.07851e-05,8.53854e-09,-9.29534e-13,-9408.41,23.4051], Tmin=(100,'K'), Tmax=(2068.45,'K')), NASAPolynomial(coeffs=[16.4092,0.021517,-1.02052e-05,1.90555e-09,-1.27843e-13,-15479.3,-58.0588], Tmin=(2068.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.8963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = 'CC(=O)OC=C[O](4426)',
    structure = SMILES('CC(=O)OC=C[O]'),
    E0 = (-417.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0554,0.0494914,-6.04122e-06,-3.88539e-08,2.13698e-11,-50081.5,23.2762], Tmin=(100,'K'), Tmax=(948.622,'K')), NASAPolynomial(coeffs=[19.5462,0.00928542,-2.17874e-06,3.96353e-10,-3.36394e-14,-55288.8,-73.9119], Tmin=(948.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-417.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'CC([O])OC=C=O(4427)',
    structure = SMILES('CC([O])OC=C=O'),
    E0 = (-187.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,361.1,361.101,361.157,361.176],'cm^-1')),
        HinderedRotor(inertia=(0.188154,'amu*angstrom^2'), symmetry=1, barrier=(17.4042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188066,'amu*angstrom^2'), symmetry=1, barrier=(17.4051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188131,'amu*angstrom^2'), symmetry=1, barrier=(17.4046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.751499,0.0639771,-6.31218e-05,3.08052e-08,-5.86287e-12,-22409.7,25.0314], Tmin=(100,'K'), Tmax=(1285.02,'K')), NASAPolynomial(coeffs=[16.5451,0.0148153,-5.73592e-06,1.03376e-09,-7.0893e-14,-26468.8,-55.1243], Tmin=(1285.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + radical(CCOJ)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]C=COC=O(830)',
    structure = SMILES('[O]C=COC=O'),
    E0 = (-361.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40903,'amu*angstrom^2'), symmetry=1, barrier=(32.3963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41581,'amu*angstrom^2'), symmetry=1, barrier=(32.5523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64366,0.0361947,7.84074e-06,-4.97625e-08,2.51765e-11,-43382.7,18.8534], Tmin=(100,'K'), Tmax=(938.333,'K')), NASAPolynomial(coeffs=[19.2694,0.000634451,1.42126e-06,-2.52749e-10,9.92051e-15,-48432.8,-74.3425], Tmin=(938.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-361.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(C=COJ)"""),
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
    label = 'C[C](O)OC=C[O](4428)',
    structure = SMILES('C[C](O)OC=C[O]'),
    E0 = (-247.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311307,0.0617877,-1.47335e-05,-4.73121e-08,2.88608e-11,-29655.3,27.0408], Tmin=(100,'K'), Tmax=(921.236,'K')), NASAPolynomial(coeffs=[26.1138,0.000294792,3.09843e-06,-6.63288e-10,4.02709e-14,-36554,-106.965], Tmin=(921.236,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-247.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Cs_P) + radical(C=COJ)"""),
)

species(
    label = 'CC([O])OC=[C]O(4429)',
    structure = SMILES('CC([O])OC=[C]O'),
    E0 = (-129.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,273.092,273.092,273.093,273.093],'cm^-1')),
        HinderedRotor(inertia=(0.352588,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352588,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352589,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352589,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.386876,0.076503,-7.89971e-05,3.84613e-08,-6.97281e-12,-15349.8,31.4519], Tmin=(100,'K'), Tmax=(1548.03,'K')), NASAPolynomial(coeffs=[22.2956,0.00664902,-4.15201e-07,-7.24954e-11,7.97959e-15,-21025.2,-83.5382], Tmin=(1548.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(O)O[CH]C=O(1252)',
    structure = SMILES('[CH2]C(O)OC=C[O]'),
    E0 = (-241.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.119904,0.0662626,-2.40814e-05,-4.03739e-08,2.7147e-11,-28885.7,27.6738], Tmin=(100,'K'), Tmax=(917.327,'K')), NASAPolynomial(coeffs=[27.1486,-0.00121599,3.87825e-06,-8.23872e-10,5.20689e-14,-35964.3,-111.945], Tmin=(917.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])O[C]=CO(4430)',
    structure = SMILES('CC([O])O[C]=CO'),
    E0 = (-129.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,273.092,273.092,273.093,273.093],'cm^-1')),
        HinderedRotor(inertia=(0.352588,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352588,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352589,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352589,'amu*angstrom^2'), symmetry=1, barrier=(18.6601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.386876,0.076503,-7.89971e-05,3.84613e-08,-6.97281e-12,-15349.8,31.4519], Tmin=(100,'K'), Tmax=(1548.03,'K')), NASAPolynomial(coeffs=[22.2956,0.00664902,-4.15201e-07,-7.24954e-11,7.97959e-15,-21025.2,-83.5382], Tmin=(1548.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=CJO)"""),
)

species(
    label = 'CC(O)O[C]=C[O](4431)',
    structure = SMILES('CC(O)O[C]=C[O]'),
    E0 = (-213.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,299.173,299.173,299.173,299.174],'cm^-1')),
        HinderedRotor(inertia=(0.29285,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292851,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29285,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292849,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481425,0.0628491,-3.13386e-05,-2.04695e-08,1.68954e-11,-25516.8,28.8704], Tmin=(100,'K'), Tmax=(932.504,'K')), NASAPolynomial(coeffs=[22.2007,0.00667976,-4.96756e-07,2.6317e-11,-5.67288e-15,-31176,-83.0193], Tmin=(932.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'CC(O)OC=[C][O](4432)',
    structure = SMILES('CC(O)OC=[C][O]'),
    E0 = (-213.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,1380,1390,370,380,2900,435,299.173,299.173,299.173,299.174],'cm^-1')),
        HinderedRotor(inertia=(0.29285,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292851,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29285,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.292849,'amu*angstrom^2'), symmetry=1, barrier=(18.6003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481425,0.0628491,-3.13386e-05,-2.04695e-08,1.68954e-11,-25516.8,28.8704], Tmin=(100,'K'), Tmax=(932.504,'K')), NASAPolynomial(coeffs=[22.2007,0.00667976,-4.96756e-07,2.6317e-11,-5.67288e-15,-31176,-83.0193], Tmin=(932.504,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'C[C]([O])OC=CO(4433)',
    structure = SMILES('CC(=O)O[CH][CH]O'),
    E0 = (-249.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861956,0.0735279,-9.50323e-05,6.8783e-08,-2.03197e-11,-29849.6,26.8384], Tmin=(100,'K'), Tmax=(821.992,'K')), NASAPolynomial(coeffs=[10.0592,0.0287744,-1.33692e-05,2.55478e-09,-1.78211e-13,-31361.7,-15.7306], Tmin=(821.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + radical(CCsJOH) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH2]C([O])OC=CO(4434)',
    structure = SMILES('[CH2]C([O])OC=CO'),
    E0 = (-157.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.992118,'amu*angstrom^2'), symmetry=1, barrier=(22.8107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993702,'amu*angstrom^2'), symmetry=1, barrier=(22.8472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993736,'amu*angstrom^2'), symmetry=1, barrier=(22.8479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.992137,'amu*angstrom^2'), symmetry=1, barrier=(22.8112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28474,0.0861001,-9.28524e-05,4.54608e-08,-8.08676e-12,-18694.9,32.1907], Tmin=(100,'K'), Tmax=(1643.86,'K')), NASAPolynomial(coeffs=[25.6072,0.000869774,2.98096e-06,-7.29425e-10,5.18918e-14,-24861.8,-102.779], Tmin=(1643.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = 'C[CH][O](176)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOH) + radical(CCOJ)"""),
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
    label = 'CC([O])[O](2201)',
    structure = SMILES('CC([O])[O]'),
    E0 = (-5.39175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2002.46,2002.55],'cm^-1')),
        HinderedRotor(inertia=(0.236654,'amu*angstrom^2'), symmetry=1, barrier=(5.44115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27157,0.0126589,2.96378e-06,-5.02044e-09,1.04673e-12,-675.327,10.9923], Tmin=(100,'K'), Tmax=(2092.59,'K')), NASAPolynomial(coeffs=[11.6822,0.0123424,-6.73655e-06,1.23259e-09,-7.81597e-14,-6809.02,-37.4771], Tmin=(2092.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.39175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'CC([O])OC1[CH]O1(4435)',
    structure = SMILES('CC([O])OC1[CH]O1'),
    E0 = (-74.8437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.993813,0.063591,-6.45611e-05,3.70881e-08,-8.6119e-12,-8891.22,26.1972], Tmin=(100,'K'), Tmax=(1047.19,'K')), NASAPolynomial(coeffs=[11.2438,0.0244382,-8.47802e-06,1.38396e-09,-8.80294e-14,-11037.9,-23.7257], Tmin=(1047.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.8437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CCOJ)"""),
)

species(
    label = 'CC1O[CH]C([O])O1(4436)',
    structure = SMILES('CC1O[CH]C([O])O1'),
    E0 = (-165.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.33335,0.0596715,-6.00175e-05,3.1061e-08,-5.8363e-12,-19802,24.5558], Tmin=(100,'K'), Tmax=(1625.87,'K')), NASAPolynomial(coeffs=[12.9378,0.011108,9.81056e-07,-5.91096e-10,5.17018e-14,-21580.5,-35.2445], Tmin=(1625.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(1,3-Dioxolane) + radical(CCsJOCs) + radical(CCOJ)"""),
)

species(
    label = 'CC(=O)OC=CO(4437)',
    structure = SMILES('CC(=O)OC=CO'),
    E0 = (-558.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641227,0.056489,-9.51379e-06,-4.5038e-08,2.61347e-11,-67078.5,23.3529], Tmin=(100,'K'), Tmax=(927.414,'K')), NASAPolynomial(coeffs=[22.8349,0.00553841,4.79189e-07,-1.66541e-10,6.6316e-15,-73120.5,-92.4268], Tmin=(927.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-558.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CC(O)OC=C=O(4438)',
    structure = SMILES('CC(O)OC=C=O'),
    E0 = (-413.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.519584,0.0652785,-4.62877e-05,7.38979e-10,7.89987e-12,-49543.4,24.9984], Tmin=(100,'K'), Tmax=(947.69,'K')), NASAPolynomial(coeffs=[19.9604,0.0104896,-2.72571e-06,4.54533e-10,-3.40037e-14,-54452.6,-74.2079], Tmin=(947.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH)"""),
)

species(
    label = 'CH2(S)(23)',
    structure = SMILES('[CH2]'),
    E0 = (419.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.36,2789.41,2993.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(200,'K'), Tmax=(3000,'K'), E0=(419.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[O]CO[CH]C=O(971)',
    structure = SMILES('[O]C=COC[O]'),
    E0 = (-179.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.40766,'amu*angstrom^2'), symmetry=1, barrier=(32.365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40681,'amu*angstrom^2'), symmetry=1, barrier=(32.3454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4132.15,'J/mol'), sigma=(6.57465,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=645.43 K, Pc=32.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28928,0.0348407,4.20423e-05,-1.036e-07,4.89469e-11,-21510,20.2394], Tmin=(100,'K'), Tmax=(916.929,'K')), NASAPolynomial(coeffs=[26.1607,-0.00887106,7.56513e-06,-1.47529e-09,9.27251e-14,-28794.6,-112.445], Tmin=(916.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(OCOJ)"""),
)

species(
    label = 'CC1OC=COO1(4439)',
    structure = SMILES('CC1OC=COO1'),
    E0 = (-221.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42236,0.0301566,7.11353e-05,-1.33832e-07,5.97258e-11,-26577.2,19.1779], Tmin=(100,'K'), Tmax=(909.348,'K')), NASAPolynomial(coeffs=[24.0119,0.000488408,5.1058e-06,-1.13761e-09,7.31019e-14,-33567.3,-103.502], Tmin=(909.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene)"""),
)

species(
    label = 'CC([O])C([O])C=O(1189)',
    structure = SMILES('CC([O])C([O])C=O'),
    E0 = (-95.3615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13498,0.0552044,-4.16136e-05,1.44536e-08,-1.61151e-12,-11359.7,29.6857], Tmin=(100,'K'), Tmax=(1185.98,'K')), NASAPolynomial(coeffs=[13.766,0.0207402,-8.31549e-06,1.52102e-09,-1.04838e-13,-14927.9,-35.8185], Tmin=(1185.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.3615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(C=OCOJ)"""),
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
    label = 'C[CH]OC=C[O](2538)',
    structure = SMILES('C[CH]OC=C[O]'),
    E0 = (-58.9096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.08807,'amu*angstrom^2'), symmetry=1, barrier=(25.0168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0912,'amu*angstrom^2'), symmetry=1, barrier=(25.0888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09168,'amu*angstrom^2'), symmetry=1, barrier=(25.0999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.555272,0.0558089,-4.93052e-06,-5.81091e-08,3.35106e-11,-6942.47,21.7927], Tmin=(100,'K'), Tmax=(909.121,'K')), NASAPolynomial(coeffs=[26.2914,-0.00399611,5.58725e-06,-1.17534e-09,7.70395e-14,-13829.9,-112.061], Tmin=(909.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.9096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOC(O)) + radical(C=COJ)"""),
)

species(
    label = '[CH]=COC(C)[O](3218)',
    structure = SMILES('[CH]=COC(C)[O]'),
    E0 = (87.0633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,325.011,325.138,325.282],'cm^-1')),
        HinderedRotor(inertia=(0.251601,'amu*angstrom^2'), symmetry=1, barrier=(18.8467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251564,'amu*angstrom^2'), symmetry=1, barrier=(18.8472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251192,'amu*angstrom^2'), symmetry=1, barrier=(18.8471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13663,0.0536379,-3.41505e-05,2.91542e-10,5.10229e-12,10582.7,23.3849], Tmin=(100,'K'), Tmax=(996.725,'K')), NASAPolynomial(coeffs=[15.7423,0.0148021,-5.47126e-06,1.01833e-09,-7.36361e-14,6688.61,-51.9597], Tmin=(996.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'CC([O])OC[C]=O(4440)',
    structure = SMILES('CC([O])OC[C]=O'),
    E0 = (-159.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,337.49,337.516,337.553,1524.89],'cm^-1')),
        HinderedRotor(inertia=(0.0915187,'amu*angstrom^2'), symmetry=1, barrier=(7.39815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0915162,'amu*angstrom^2'), symmetry=1, barrier=(7.39876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0915264,'amu*angstrom^2'), symmetry=1, barrier=(7.3993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0166303,'amu*angstrom^2'), symmetry=1, barrier=(27.4404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25608,0.0481103,1.50348e-06,-1.15898e-07,1.21773e-10,-19169,24.3742], Tmin=(100,'K'), Tmax=(423.479,'K')), NASAPolynomial(coeffs=[4.9392,0.0383713,-1.92728e-05,3.82344e-09,-2.72574e-13,-19536.2,12.0831], Tmin=(423.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CsCJ=O)"""),
)

species(
    label = 'C[C]([O])OCC=O(4441)',
    structure = SMILES('C[C]([O])OCC=O'),
    E0 = (-109.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,261.544,262.405,1907.17,1907.44],'cm^-1')),
        HinderedRotor(inertia=(0.168465,'amu*angstrom^2'), symmetry=1, barrier=(8.32286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.479557,'amu*angstrom^2'), symmetry=1, barrier=(23.4617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1698,'amu*angstrom^2'), symmetry=1, barrier=(8.31715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909255,'amu*angstrom^2'), symmetry=1, barrier=(23.4532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23134,0.0486805,-6.61278e-07,-1.11855e-07,1.18456e-10,-13085.3,24.4688], Tmin=(100,'K'), Tmax=(426.096,'K')), NASAPolynomial(coeffs=[5.02504,0.0381957,-1.91663e-05,3.79898e-09,-2.70625e-13,-13466.3,11.6973], Tmin=(426.096,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])OCC=O(4442)',
    structure = SMILES('[CH2]C([O])OCC=O'),
    E0 = (-102.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,276.628,276.633,1432.33,1432.35],'cm^-1')),
        HinderedRotor(inertia=(0.163412,'amu*angstrom^2'), symmetry=1, barrier=(8.87488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163394,'amu*angstrom^2'), symmetry=1, barrier=(8.87517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00609609,'amu*angstrom^2'), symmetry=1, barrier=(8.87532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519756,'amu*angstrom^2'), symmetry=1, barrier=(28.2277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7915,0.0577736,-3.94012e-05,-2.86421e-08,4.81056e-11,-12305.7,25.9086], Tmin=(100,'K'), Tmax=(483.472,'K')), NASAPolynomial(coeffs=[6.07054,0.0366874,-1.83972e-05,3.64265e-09,-2.59292e-13,-12886.7,6.64414], Tmin=(483.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJCO) + radical(CCOJ)"""),
)

species(
    label = 'CC1O[CH][CH]OO1(4443)',
    structure = SMILES('CC1O[CH][CH]OO1'),
    E0 = (39.4351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14305,0.0452353,2.61745e-05,-8.99987e-08,4.72924e-11,4862.89,20.5244], Tmin=(100,'K'), Tmax=(853.028,'K')), NASAPolynomial(coeffs=[21.2185,0.00390968,5.97681e-06,-1.63571e-09,1.24922e-13,-483.529,-84.3992], Tmin=(853.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.4351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(124trioxane) + radical(CCsJOOC) + radical(CCsJOCs)"""),
)

species(
    label = 'CC(=O)OCC=O(4444)',
    structure = SMILES('CC(=O)OCC=O'),
    E0 = (-530.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67609,0.0436138,-1.58222e-05,-4.62299e-09,3.10658e-12,-63663.7,24.4019], Tmin=(100,'K'), Tmax=(1218.64,'K')), NASAPolynomial(coeffs=[10.8398,0.0269605,-1.18488e-05,2.24339e-09,-1.56575e-13,-66894,-25.7099], Tmin=(1218.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-530.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'CC1OC(C=O)O1(1194)',
    structure = SMILES('CC1OC(C=O)O1'),
    E0 = (-406.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89033,0.032606,2.52235e-05,-5.1419e-08,1.9981e-11,-48752.4,22.6429], Tmin=(100,'K'), Tmax=(1043.26,'K')), NASAPolynomial(coeffs=[12.5488,0.0248339,-1.11841e-05,2.25248e-09,-1.66843e-13,-52777.2,-37.8607], Tmin=(1043.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane)"""),
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
    E0 = (-227.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-78.8963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-159.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (56.9368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-193.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-120.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-69.3265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (33.6978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-88.6849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (21.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-89.6099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-169.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-139.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-76.4852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (138.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (240.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-45.5702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-165.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-202.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-202.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (240.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-219.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (86.4371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (184.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (330.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-45.9427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-26.2428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (50.1686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-21.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (39.4351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-138.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-219.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['OCHCHO(48)', 'CH3CHO(37)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC1OC([CH][O])O1(4425)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(148.467,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 147.2 to 148.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'CC(=O)OC=C[O](4426)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'CC([O])OC=C=O(4427)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.54607e+06,'m^3/(mol*s)'), n=0.466452, Ea=(32.4943,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(17)', '[O]C=COC=O(830)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-NdH_O;CsJ-HHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3CHO(37)', '[O]C=C[O](2537)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.148573,'m^3/(mol*s)'), n=1.743, Ea=(77.4229,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;O_rad/OneDe] + [CO-CsH_O;OJ_sec] for rate rule [CO-CsH_O;O_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['C[C](O)OC=C[O](4428)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC([O])OC=[C]O(4429)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['[CH2]C(O)O[CH]C=O(1252)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC([O])O[C]=CO(4430)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC(O)O[C]=C[O](4431)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.3925e+09,'s^-1'), n=0.84206, Ea=(123.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS_O(Cs)Cs;Y_rad_out;XH_out] for rate rule [R4H_SSS_O(Cs)Cs;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CC(O)OC=[C][O](4432)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['C[C]([O])OC=CO(4433)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;XH_out] for rate rule [R5H_SMSS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([O])OC=CO(4434)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(12584.6,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;C_rad_out_2H;XH_out] for rate rule [R6H_RSSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH][O](176)', '[O]C=C[O](2537)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['CHCHO(47)', 'CC([O])[O](2201)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.75733e+07,'m^3/(mol*s)'), n=-0.180109, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/NonDe] for rate rule [Cd_rad;O_rad/NonDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC([O])OC1[CH]O1(4435)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri_HNd_O;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC1O[CH]C([O])O1(4436)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(61.4715,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 59.1 to 61.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC(=O)OC=CO(4437)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC(O)OC=C=O(4438)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', '[O]CO[CH]C=O(971)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC1OC=COO1(4439)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC([O])C([O])C=O(1189)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', 'C[CH]OC=C[O](2538)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['O(4)', '[CH]=COC(C)[O](3218)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[CH][O](176)', 'OCHCHO(48)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.2635,'m^3/(mol*s)'), n=2.48, Ea=(21.7568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-DeH;YJ] for rate rule [Od_CO-DeH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC([O])OC[C]=O(4440)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[C]([O])OCC=O(4441)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([O])OCC=O(4442)'],
    products = ['CC([O])O[CH]C=O(1190)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(46663.9,'s^-1'), n=1.985, Ea=(81.1696,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC1O[CH][CH]OO1(4443)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(266.798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 265.4 to 266.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC(=O)OCC=O(4444)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CC([O])O[CH]C=O(1190)'],
    products = ['CC1OC(C=O)O1(1194)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '774',
    isomers = [
        'CC([O])O[CH]C=O(1190)',
    ],
    reactants = [
        ('OCHCHO(48)', 'CH3CHO(37)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '774',
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

