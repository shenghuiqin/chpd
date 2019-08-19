species(
    label = 'CCCC(=O)CO[O](2417)',
    structure = SMILES('CCCC(=O)CO[O]'),
    E0 = (-203.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.760279,0.0800334,-7.74372e-05,2.88842e-08,6.68474e-12,-24393.3,28.4196], Tmin=(100,'K'), Tmax=(572.595,'K')), NASAPolynomial(coeffs=[7.66965,0.0451413,-2.10701e-05,4.05085e-09,-2.84048e-13,-25403.8,-2.97623], Tmin=(572.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C(=O)CCC(2148)',
    structure = SMILES('C=C([O])CCC'),
    E0 = (-117.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,335.188,336.29,336.687],'cm^-1')),
        HinderedRotor(inertia=(0.00148026,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981052,'amu*angstrom^2'), symmetry=1, barrier=(7.83672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992182,'amu*angstrom^2'), symmetry=1, barrier=(7.84609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.1,'J/mol'), sigma=(6.27192,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.76 K, Pc=33.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26061,0.0511177,-1.77487e-05,-1.26003e-08,8.41823e-12,-13971.9,23.7131], Tmin=(100,'K'), Tmax=(999.148,'K')), NASAPolynomial(coeffs=[12.5223,0.0258178,-9.4698e-06,1.69475e-09,-1.17565e-13,-17209.9,-35.5505], Tmin=(999.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCCC1([O])COO1(2919)',
    structure = SMILES('CCCC1([O])COO1'),
    E0 = (-79.8716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.340004,0.0728112,-5.29432e-05,1.94924e-08,-2.91285e-12,-9468.44,25.7144], Tmin=(100,'K'), Tmax=(1557.17,'K')), NASAPolynomial(coeffs=[16.5416,0.0311931,-1.28531e-05,2.3287e-09,-1.57264e-13,-14514.2,-59.6242], Tmin=(1557.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.8716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'CCCC(=O)[CH]OO(2920)',
    structure = SMILES('CCCC(=O)[CH]OO'),
    E0 = (-217.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.160229,0.0914291,-0.000113204,8.57963e-08,-2.76466e-11,-25995.9,29.5093], Tmin=(100,'K'), Tmax=(745.276,'K')), NASAPolynomial(coeffs=[8.5943,0.0461594,-2.20844e-05,4.28277e-09,-3.01399e-13,-27252.9,-8.7], Tmin=(745.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-217.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(OCJC=O)"""),
)

species(
    label = 'CC[CH]C(=O)COO(2921)',
    structure = SMILES('CCC=C([O])COO'),
    E0 = (-210.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35415,0.0795831,-7.27017e-05,3.59412e-08,-7.31857e-12,-25237.5,31.0475], Tmin=(100,'K'), Tmax=(1164.28,'K')), NASAPolynomial(coeffs=[13.4516,0.0345848,-1.47271e-05,2.7444e-09,-1.90255e-13,-28287.3,-34.1321], Tmin=(1164.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]CC(=O)COO(2922)',
    structure = SMILES('C[CH]CC(=O)COO'),
    E0 = (-155.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.567171,0.0809397,-7.97257e-05,4.5057e-08,-1.08997e-11,-18621.8,30.7924], Tmin=(100,'K'), Tmax=(969.458,'K')), NASAPolynomial(coeffs=[9.95241,0.0422144,-1.98053e-05,3.8499e-09,-2.72908e-13,-20441.5,-14.1945], Tmin=(969.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-155.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]CCC(=O)COO(2923)',
    structure = SMILES('[CH2]CCC(=O)COO'),
    E0 = (-150.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.300972,0.088133,-0.000103024,7.22304e-08,-2.15609e-11,-17971,31.0402], Tmin=(100,'K'), Tmax=(800.349,'K')), NASAPolynomial(coeffs=[9.0265,0.0445232,-2.12891e-05,4.14575e-09,-2.93143e-13,-19367.7,-9.11193], Tmin=(800.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-150.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJ)"""),
)

species(
    label = 'C2H5(29)',
    structure = SMILES('C[CH2]'),
    E0 = (107.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1190.6,1642.82,1642.96,3622.23,3622.39],'cm^-1')),
        HinderedRotor(inertia=(0.866817,'amu*angstrom^2'), symmetry=1, barrier=(19.9298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(107.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(=O)CO[O](2924)',
    structure = SMILES('[CH2]C(=O)CO[O]'),
    E0 = (18.6224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,229.744],'cm^-1')),
        HinderedRotor(inertia=(0.559066,'amu*angstrom^2'), symmetry=1, barrier=(12.854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.55905,'amu*angstrom^2'), symmetry=1, barrier=(12.8537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0314479,'amu*angstrom^2'), symmetry=1, barrier=(12.8542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45265,0.0560523,-6.56697e-05,3.91582e-08,-9.18676e-12,2331.54,21.428], Tmin=(100,'K'), Tmax=(1042.64,'K')), NASAPolynomial(coeffs=[12.2101,0.0147824,-6.29693e-06,1.19535e-09,-8.42293e-14,88.2915,-30.9199], Tmin=(1042.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.6224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C2JC=O) + radical(ROOJ)"""),
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
    label = '[CH2]CC(=O)CO[O](2925)',
    structure = SMILES('[CH2]CC(=O)CO[O]'),
    E0 = (31.6511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,238.19,238.19],'cm^-1')),
        HinderedRotor(inertia=(0.135987,'amu*angstrom^2'), symmetry=1, barrier=(5.47484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135987,'amu*angstrom^2'), symmetry=1, barrier=(5.47484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135987,'amu*angstrom^2'), symmetry=1, barrier=(5.47484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135987,'amu*angstrom^2'), symmetry=1, barrier=(5.47484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.776331,0.077407,-0.000115774,9.87877e-08,-3.36983e-11,3916.58,26.6748], Tmin=(100,'K'), Tmax=(822.639,'K')), NASAPolynomial(coeffs=[8.15596,0.0318438,-1.50425e-05,2.85066e-09,-1.95857e-13,3029.98,-5.49597], Tmin=(822.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.6511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(ROOJ)"""),
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
    label = 'C[CH]CC(=O)CO[O](2926)',
    structure = SMILES('C[CH]CC(=O)CO[O]'),
    E0 = (-3.81614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.62653,0.0798094,-9.33123e-05,6.64134e-08,-2.01064e-11,-342.583,31.3348], Tmin=(100,'K'), Tmax=(790.835,'K')), NASAPolynomial(coeffs=[8.39188,0.0405324,-1.88138e-05,3.61131e-09,-2.53106e-13,-1570.79,-4.30623], Tmin=(790.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.81614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(ROOJ)"""),
)

species(
    label = 'npropyl(83)',
    structure = SMILES('[CH2]CC'),
    E0 = (87.0621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0928812,'amu*angstrom^2'), symmetry=1, barrier=(2.13552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.092914,'amu*angstrom^2'), symmetry=1, barrier=(2.13628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0877,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.02815,0.0147023,2.4051e-05,-3.66738e-08,1.38611e-11,10512.1,12.4699], Tmin=(100,'K'), Tmax=(984.464,'K')), NASAPolynomial(coeffs=[6.16543,0.0184495,-6.79029e-06,1.23049e-09,-8.63866e-14,9095.06,-6.67607], Tmin=(984.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.0621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""npropyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]OC[C]=O(407)',
    structure = SMILES('[O]OC[C]=O'),
    E0 = (55.0207,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,492.5,1135,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.654,'amu*angstrom^2'), symmetry=1, barrier=(15.0368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34052,'amu*angstrom^2'), symmetry=1, barrier=(7.82922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0355,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63093,0.0324223,-3.97506e-05,2.80677e-08,-8.25459e-12,6664.69,17.1477], Tmin=(100,'K'), Tmax=(818.661,'K')), NASAPolynomial(coeffs=[6.29653,0.0145127,-6.93676e-06,1.34708e-09,-9.50121e-14,6064.49,0.19668], Tmin=(818.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.0207,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(ROOJ)"""),
)

species(
    label = 'CC[CH]C(=O)CO[O](2927)',
    structure = SMILES('CCC=C([O])CO[O]'),
    E0 = (-58.9272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3010,987.5,1337.5,450,1655,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,294.116,294.122,3815.49],'cm^-1')),
        HinderedRotor(inertia=(0.177638,'amu*angstrom^2'), symmetry=1, barrier=(10.9092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177695,'amu*angstrom^2'), symmetry=1, barrier=(10.9093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177674,'amu*angstrom^2'), symmetry=1, barrier=(10.9092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177684,'amu*angstrom^2'), symmetry=1, barrier=(10.9091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640894,0.0756591,-7.60292e-05,4.33127e-08,-1.02458e-11,-6967.81,30.7824], Tmin=(100,'K'), Tmax=(1009.5,'K')), NASAPolynomial(coeffs=[11.0892,0.0342594,-1.45147e-05,2.68918e-09,-1.85617e-13,-9077.35,-19.7235], Tmin=(1009.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.9272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CCC(=O)CO[O](2928)',
    structure = SMILES('[CH2]CCC(=O)CO[O]'),
    E0 = (1.52816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,492.5,1135,1000,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323268,0.0876914,-0.000120267,1.00207e-07,-3.44813e-11,309.66,31.6998], Tmin=(100,'K'), Tmax=(786.804,'K')), NASAPolynomial(coeffs=[7.95413,0.041928,-1.97344e-05,3.7675e-09,-2.61332e-13,-675.419,-1.91415], Tmin=(786.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.52816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]O[O](92)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (206.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2219.2,2221.31],'cm^-1')),
        HinderedRotor(inertia=(0.156598,'amu*angstrom^2'), symmetry=1, barrier=(14.2064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45042,0.0115293,-8.64543e-06,3.98163e-09,-7.73155e-13,24893.2,10.7919], Tmin=(100,'K'), Tmax=(1212.82,'K')), NASAPolynomial(coeffs=[5.07483,0.00617178,-2.01911e-06,3.39169e-10,-2.23136e-14,24499.1,2.64172], Tmin=(1212.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'CCC[C]=O(2154)',
    structure = SMILES('CCC[C]=O'),
    E0 = (-67.1428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.225289,'amu*angstrom^2'), symmetry=1, barrier=(5.17985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225185,'amu*angstrom^2'), symmetry=1, barrier=(5.17744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225204,'amu*angstrom^2'), symmetry=1, barrier=(5.17788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32545,0.0389501,-2.48868e-05,8.84133e-09,-1.41386e-12,-8017.15,18.4614], Tmin=(100,'K'), Tmax=(1320.47,'K')), NASAPolynomial(coeffs=[6.17416,0.0272916,-1.16432e-05,2.15506e-09,-1.47972e-13,-9033.58,-1.17631], Tmin=(1320.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.1428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'CCCC(=O)[CH]O[O](2929)',
    structure = SMILES('CCCC(=O)[CH]O[O]'),
    E0 = (-65.2334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249699,0.0901115,-0.000126928,1.08403e-07,-3.78243e-11,-7718.06,29.933], Tmin=(100,'K'), Tmax=(804.629,'K')), NASAPolynomial(coeffs=[7.4553,0.0436829,-2.06004e-05,3.92159e-09,-2.71026e-13,-8534.23,-1.12967], Tmin=(804.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.2334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(ROOJ) + radical(OCJC=O)"""),
)

species(
    label = 'CCC[C]1COOO1(2930)',
    structure = SMILES('CCC[C]1COOO1'),
    E0 = (53.2981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.67154,0.0646153,-2.60944e-05,-1.4946e-08,1.23132e-11,6538.11,26.3135], Tmin=(100,'K'), Tmax=(906.505,'K')), NASAPolynomial(coeffs=[13.5297,0.0313096,-9.75579e-06,1.55252e-09,-1.00583e-13,3244.16,-39.7678], Tmin=(906.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(123trioxolane) + radical(C2CsJOO)"""),
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
    label = 'CCC(=O)CO[O](2931)',
    structure = SMILES('CCC(=O)CO[O]'),
    E0 = (-179.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,186.374,186.385],'cm^-1')),
        HinderedRotor(inertia=(0.229232,'amu*angstrom^2'), symmetry=1, barrier=(5.6488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229246,'amu*angstrom^2'), symmetry=1, barrier=(5.64877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229217,'amu*angstrom^2'), symmetry=1, barrier=(5.64867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229256,'amu*angstrom^2'), symmetry=1, barrier=(5.64879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.097,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.994422,0.071405,-9.3043e-05,7.56608e-08,-2.58628e-11,-21538.3,25.2611], Tmin=(100,'K'), Tmax=(768.046,'K')), NASAPolynomial(coeffs=[7.04794,0.0361439,-1.68847e-05,3.2248e-09,-2.24244e-13,-22358.1,-1.62904], Tmin=(768.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(ROOJ)"""),
)

species(
    label = 'C=C(CO[O])OCC(2932)',
    structure = SMILES('C=C(CO[O])OCC'),
    E0 = (-156.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0819025,0.0828197,-7.78057e-05,3.85945e-08,-7.73188e-12,-18729.8,30.0281], Tmin=(100,'K'), Tmax=(1197.65,'K')), NASAPolynomial(coeffs=[15.6613,0.0307865,-1.26367e-05,2.31852e-09,-1.59564e-13,-22461.5,-47.9434], Tmin=(1197.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-156.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'CCC=C(O)CO[O](2933)',
    structure = SMILES('CCC=C(O)CO[O]'),
    E0 = (-196.732,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,301.9,302.772],'cm^-1')),
        HinderedRotor(inertia=(0.167588,'amu*angstrom^2'), symmetry=1, barrier=(10.4152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160978,'amu*angstrom^2'), symmetry=1, barrier=(10.4622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162002,'amu*angstrom^2'), symmetry=1, barrier=(10.4519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160996,'amu*angstrom^2'), symmetry=1, barrier=(10.4147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159498,'amu*angstrom^2'), symmetry=1, barrier=(10.4196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0430865,0.0824486,-7.86753e-05,3.96972e-08,-8.03728e-12,-23515.1,31.4199], Tmin=(100,'K'), Tmax=(1192.2,'K')), NASAPolynomial(coeffs=[16.0021,0.0289037,-1.13063e-05,2.02511e-09,-1.37588e-13,-27320.4,-48.3789], Tmin=(1192.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ)"""),
)

species(
    label = 'C=C(CCC)OO[O](2416)',
    structure = SMILES('C=C(CCC)OO[O]'),
    E0 = (57.6159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390086,0.071619,-5.66437e-05,2.34223e-08,-3.92445e-12,7065.93,34.2599], Tmin=(100,'K'), Tmax=(1413.66,'K')), NASAPolynomial(coeffs=[15.5298,0.0287799,-1.11874e-05,1.98526e-09,-1.33334e-13,2785.52,-44.0213], Tmin=(1413.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.6159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'CCCC(O)=CO[O](2934)',
    structure = SMILES('CCCC(O)=CO[O]'),
    E0 = (-187.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,236.398,237.124],'cm^-1')),
        HinderedRotor(inertia=(0.286247,'amu*angstrom^2'), symmetry=1, barrier=(11.355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.286206,'amu*angstrom^2'), symmetry=1, barrier=(11.3442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280903,'amu*angstrom^2'), symmetry=1, barrier=(11.3398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282294,'amu*angstrom^2'), symmetry=1, barrier=(11.3513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283256,'amu*angstrom^2'), symmetry=1, barrier=(11.3488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.123,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.597996,0.0875453,-8.55518e-05,4.22925e-08,-8.10657e-12,-22332.3,33.5492], Tmin=(100,'K'), Tmax=(1336.22,'K')), NASAPolynomial(coeffs=[21.1472,0.0195194,-5.89739e-06,9.09602e-10,-5.68885e-14,-27881.9,-76.6825], Tmin=(1336.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
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
    label = 'CCCC(=O)C[O](2935)',
    structure = SMILES('CCCC(=O)C[O]'),
    E0 = (-183.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,342.322,342.537,1682.05,3816.06],'cm^-1')),
        HinderedRotor(inertia=(0.121961,'amu*angstrom^2'), symmetry=1, barrier=(10.1654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12232,'amu*angstrom^2'), symmetry=1, barrier=(10.1615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122253,'amu*angstrom^2'), symmetry=1, barrier=(10.1678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121933,'amu*angstrom^2'), symmetry=1, barrier=(10.1638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10929,0.0684512,-6.86881e-05,4.46276e-08,-1.29268e-11,-21969.6,26.3665], Tmin=(100,'K'), Tmax=(809.044,'K')), NASAPolynomial(coeffs=[6.53895,0.0416062,-1.89161e-05,3.61426e-09,-2.5329e-13,-22848.1,1.32211], Tmin=(809.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + radical(C=OCOJ)"""),
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
    E0 = (-125.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-79.8716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-44.3082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-152.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-124.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-120.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (126.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (167.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (207.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (142.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (155.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (213.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (139.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (149.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (53.2981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (239.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (156.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-68.8848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (261.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-44.0401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (59.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', '[CH2]C(=O)CCC(2148)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(311.435,'m^3/(mol*s)'), n=1.07778, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;O2_birad] for rate rule [C_rad/H2/CO;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -3.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CCCC(=O)CO[O](2417)'],
    products = ['CCCC1([O])COO1(2919)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(123.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_Nd;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 122.0 to 123.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CCCC(=O)CO[O](2417)'],
    products = ['CCCC(=O)[CH]OO(2920)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.32e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_O;O_rad_out;Cs_H_out_H/OneDe] for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CCCC(=O)CO[O](2417)'],
    products = ['CC[CH]C(=O)COO(2921)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_SSSS;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CCCC(=O)CO[O](2417)'],
    products = ['C[CH]CC(=O)COO(2922)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.12e+09,'s^-1'), n=0, Ea=(79.7052,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 141 used for R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R6H_SSSSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CCCC(=O)CO[O](2417)'],
    products = ['[CH2]CCC(=O)COO(2923)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(271800,'s^-1'), n=1.51, Ea=(83.4708,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 266 used for R7H_OOCs4;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R7H_OOCs4;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5(29)', '[CH2]C(=O)CO[O](2924)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.33778e+08,'m^3/(mol*s)'), n=-0.3495, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;C_rad/H2/Cs] for rate rule [C_rad/H2/CO;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH2]CC(=O)CO[O](2925)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.23e+15,'cm^3/(mol*s)'), n=-0.562, Ea=(0.085772,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 10 used for C_methyl;C_rad/H2/Cs
Exact match found for rate rule [C_rad/H2/Cs;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C[CH]CC(=O)CO[O](2926)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['npropyl(83)', '[O]OC[C]=O(407)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.12e+14,'cm^3/(mol*s)','*|/',3), n=-0.5, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 79 used for C_rad/H2/Cs;CO_rad/NonDe
Exact match found for rate rule [CO_rad/NonDe;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'CC[CH]C(=O)CO[O](2927)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/OneDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2]CCC(=O)CO[O](2928)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]O[O](92)', 'CCC[C]=O(2154)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.12e+14,'cm^3/(mol*s)','*|/',3), n=-0.5, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_pri_rad;CO_rad/NonDe] for rate rule [C_rad/H2/O;CO_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'CCCC(=O)[CH]O[O](2929)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CCCC(=O)CO[O](2417)'],
    products = ['CCC[C]1COOO1(2930)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(257.016,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 254.3 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2(S)(23)', 'CCC(=O)CO[O](2931)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C(CO[O])OCC(2932)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_Cs;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CCC=C(O)CO[O](2933)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(127.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_Cs;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(CCC)OO[O](2416)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CCCC(O)=CO[O](2934)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O(4)', 'CCCC(=O)C[O](2935)'],
    products = ['CCCC(=O)CO[O](2417)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '1046',
    isomers = [
        'CCCC(=O)CO[O](2417)',
    ],
    reactants = [
        ('O2(2)', '[CH2]C(=O)CCC(2148)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1046',
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

