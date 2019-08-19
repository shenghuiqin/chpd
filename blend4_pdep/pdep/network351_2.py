species(
    label = 'S(239)(238)',
    structure = SMILES('C=C(C=CCO)O[O]'),
    E0 = (-32.6335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.339356,'amu*angstrom^2'), symmetry=1, barrier=(7.80246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339181,'amu*angstrom^2'), symmetry=1, barrier=(7.79843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339524,'amu*angstrom^2'), symmetry=1, barrier=(7.80633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339606,'amu*angstrom^2'), symmetry=1, barrier=(7.80821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4425.45,'J/mol'), sigma=(7.00962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=691.25 K, Pc=29.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558009,0.0797378,-9.92164e-05,7.04331e-08,-2.04433e-11,-3804.62,29.9642], Tmin=(100,'K'), Tmax=(836.603,'K')), NASAPolynomial(coeffs=[10.4809,0.0322906,-1.41388e-05,2.63199e-09,-1.80924e-13,-5464.8,-16.137], Tmin=(836.603,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.6335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'S(243)(242)',
    structure = SMILES('[CH2]C1(C=CCO)OO1'),
    E0 = (16.7182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180,180,255.273,858.27,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.14781,'amu*angstrom^2'), symmetry=1, barrier=(3.39844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14781,'amu*angstrom^2'), symmetry=1, barrier=(3.39844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14781,'amu*angstrom^2'), symmetry=1, barrier=(3.39844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14781,'amu*angstrom^2'), symmetry=1, barrier=(3.39844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4244.91,'J/mol'), sigma=(6.9998,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=663.05 K, Pc=28.08 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0804752,0.0626991,7.45742e-06,-8.18699e-08,4.4549e-11,2174.01,30.773], Tmin=(100,'K'), Tmax=(898.816,'K')), NASAPolynomial(coeffs=[28.6529,-0.000386239,5.81314e-06,-1.34225e-09,9.16569e-14,-5550.31,-118.421], Tmin=(898.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.7182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(CJCOOH)"""),
)

species(
    label = 'O2(2)(2)',
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
    label = 'C5H7O(237)(236)',
    structure = SMILES('C=C=C[CH]CO'),
    E0 = (81.1461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1976,'amu*angstrom^2'), symmetry=1, barrier=(27.5353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19778,'amu*angstrom^2'), symmetry=1, barrier=(27.5393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19446,'amu*angstrom^2'), symmetry=1, barrier=(27.463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3721.65,'J/mol'), sigma=(6.20995,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=581.31 K, Pc=35.26 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14012,0.0586746,-4.97917e-05,2.20013e-08,-3.93766e-12,9866.02,21.1068], Tmin=(100,'K'), Tmax=(1327.36,'K')), NASAPolynomial(coeffs=[12.9976,0.0229418,-9.41102e-06,1.71991e-09,-1.17746e-13,6718.22,-39.4564], Tmin=(1327.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.1461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO)"""),
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
    label = 'S(244)(243)',
    structure = SMILES('C=C([O])C=CCO'),
    E0 = (-171.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,269.284,269.284,269.285],'cm^-1')),
        HinderedRotor(inertia=(0.204289,'amu*angstrom^2'), symmetry=1, barrier=(10.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20429,'amu*angstrom^2'), symmetry=1, barrier=(10.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20429,'amu*angstrom^2'), symmetry=1, barrier=(10.5123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4284.46,'J/mol'), sigma=(6.81655,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.22 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.633106,0.0679166,-6.69863e-05,3.47092e-08,-7.10702e-12,-20479.6,26.3839], Tmin=(100,'K'), Tmax=(1192.86,'K')), NASAPolynomial(coeffs=[14.8693,0.0201783,-6.95604e-06,1.15933e-09,-7.55888e-14,-23875.9,-44.8081], Tmin=(1192.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OOC1[CH]CO(5085)',
    structure = SMILES('C=C1OOC1[CH]CO'),
    E0 = (65.3564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849139,0.0564415,-1.67674e-05,-1.82082e-08,1.03402e-11,7985.05,31.7804], Tmin=(100,'K'), Tmax=(1056.34,'K')), NASAPolynomial(coeffs=[16.0597,0.0256156,-1.10102e-05,2.15052e-09,-1.56231e-13,3277.93,-49.5051], Tmin=(1056.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.3564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCOOH)"""),
)

species(
    label = 'C=C([C]=CCO)OO(5086)',
    structure = SMILES('C=C([C]=CCO)OO'),
    E0 = (14.3572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1685,370,208.32,208.449],'cm^-1')),
        HinderedRotor(inertia=(0.298314,'amu*angstrom^2'), symmetry=1, barrier=(9.18607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.621714,'amu*angstrom^2'), symmetry=1, barrier=(19.1538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298342,'amu*angstrom^2'), symmetry=1, barrier=(9.18599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298242,'amu*angstrom^2'), symmetry=1, barrier=(9.18603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10931,'amu*angstrom^2'), symmetry=1, barrier=(64.9333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0850797,0.091522,-0.000129883,1.02929e-07,-3.26419e-11,1862.59,30.4924], Tmin=(100,'K'), Tmax=(835.899,'K')), NASAPolynomial(coeffs=[11.0956,0.0330462,-1.45642e-05,2.67484e-09,-1.80673e-13,224.039,-19.4441], Tmin=(835.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.3572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(C=CCO)OO(5087)',
    structure = SMILES('[CH]=C(C=CCO)OO'),
    E0 = (62.4579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,228.737],'cm^-1')),
        HinderedRotor(inertia=(0.296487,'amu*angstrom^2'), symmetry=1, barrier=(11.0079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.532549,'amu*angstrom^2'), symmetry=1, barrier=(19.7724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296487,'amu*angstrom^2'), symmetry=1, barrier=(11.0079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296487,'amu*angstrom^2'), symmetry=1, barrier=(11.0079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.296487,'amu*angstrom^2'), symmetry=1, barrier=(11.0079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.280561,0.0860885,-0.000108238,7.43142e-08,-2.06285e-11,7642.12,30.6078], Tmin=(100,'K'), Tmax=(876.168,'K')), NASAPolynomial(coeffs=[12.3271,0.0310914,-1.40816e-05,2.67077e-09,-1.85933e-13,5531.2,-25.9171], Tmin=(876.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.4579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(C=[C]CO)OO(5088)',
    structure = SMILES('C=C(C=[C]CO)OO'),
    E0 = (53.2035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,1685,370,290.093,290.095],'cm^-1')),
        HinderedRotor(inertia=(0.164693,'amu*angstrom^2'), symmetry=1, barrier=(9.83513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164693,'amu*angstrom^2'), symmetry=1, barrier=(9.83514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164693,'amu*angstrom^2'), symmetry=1, barrier=(9.83513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164691,'amu*angstrom^2'), symmetry=1, barrier=(9.83512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310847,'amu*angstrom^2'), symmetry=1, barrier=(18.5632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.230682,0.0883233,-0.000120884,9.28541e-08,-2.8938e-11,6529.59,30.8517], Tmin=(100,'K'), Tmax=(782.692,'K')), NASAPolynomial(coeffs=[11.0119,0.0332203,-1.52718e-05,2.88921e-09,-1.99656e-13,4842.07,-18.5189], Tmin=(782.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.2035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=C(C=C[CH]O)OO(5089)',
    structure = SMILES('C=C([CH]C=CO)OO'),
    E0 = (-93.1899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.376842,0.0826182,-6.29312e-05,6.93009e-09,7.46735e-12,-11038.2,29.6914], Tmin=(100,'K'), Tmax=(956.183,'K')), NASAPolynomial(coeffs=[23.9706,0.0127988,-3.65428e-06,6.37344e-10,-4.77829e-14,-17158.7,-94.3371], Tmin=(956.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.1899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C(C=CC[O])OO(5090)',
    structure = SMILES('C=C(C=CC[O])OO'),
    E0 = (41.0669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.180763,'amu*angstrom^2'), symmetry=1, barrier=(4.1561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180995,'amu*angstrom^2'), symmetry=1, barrier=(4.16144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182589,'amu*angstrom^2'), symmetry=1, barrier=(4.19807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05983,'amu*angstrom^2'), symmetry=1, barrier=(24.3675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2826,0.0704986,-5.3295e-05,-2.67036e-08,5.30917e-11,5026.21,27.6086], Tmin=(100,'K'), Tmax=(498.162,'K')), NASAPolynomial(coeffs=[7.37732,0.0403911,-1.93383e-05,3.73158e-09,-2.60935e-13,4185.32,0.106927], Tmin=(498.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.0669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'OH(5)(5)',
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
    label = '[CH2]C=CC(=C)O[O](5091)',
    structure = SMILES('[CH2]C=CC(=C)O[O]'),
    E0 = (238.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,268.638],'cm^-1')),
        HinderedRotor(inertia=(0.182697,'amu*angstrom^2'), symmetry=1, barrier=(9.3465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18266,'amu*angstrom^2'), symmetry=1, barrier=(9.34623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39118,'amu*angstrom^2'), symmetry=1, barrier=(71.1869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24575,0.0580629,-5.18992e-05,2.54273e-08,-5.1033e-12,28793.2,25.2898], Tmin=(100,'K'), Tmax=(1188.71,'K')), NASAPolynomial(coeffs=[11.2048,0.0245509,-9.61146e-06,1.71104e-09,-1.15501e-13,26425.5,-24.4785], Tmin=(1188.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(ROOJ)"""),
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
    label = 'C=C(C=CC[O])O[O](5092)',
    structure = SMILES('C=C(C=CC[O])O[O]'),
    E0 = (193.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,492.5,1135,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,180,180,3264.18],'cm^-1')),
        HinderedRotor(inertia=(0.406239,'amu*angstrom^2'), symmetry=1, barrier=(9.34024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40579,'amu*angstrom^2'), symmetry=1, barrier=(9.3299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406476,'amu*angstrom^2'), symmetry=1, barrier=(9.34569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789085,0.0786864,-0.000118204,1.0545e-07,-3.74245e-11,23329,29.9857], Tmin=(100,'K'), Tmax=(832.75,'K')), NASAPolynomial(coeffs=[6.05905,0.0382242,-1.80342e-05,3.41311e-09,-2.34115e-13,22976.5,8.67944], Tmin=(832.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = 'CH2OH(22)(23)',
    structure = SMILES('[CH2]O'),
    E0 = (-28.7184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3003.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.057913,'amu*angstrom^2'), symmetry=1, barrier=(25.9304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.0339,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3467.15,'J/mol'), sigma=(3.69,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.47834,-0.0013507,2.78485e-05,-3.64869e-08,1.47907e-11,-3500.73,3.30913], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.09314,0.00594761,-2.06497e-06,3.23008e-10,-1.88126e-14,-4034.1,-1.84691], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-28.7184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""CH2OH""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]=CC(=C)O[O](5093)',
    structure = SMILES('[CH]=CC(=C)O[O]'),
    E0 = (403.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.613873,'amu*angstrom^2'), symmetry=1, barrier=(14.1141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618215,'amu*angstrom^2'), symmetry=1, barrier=(14.214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61681,0.0539303,-6.76015e-05,4.54309e-08,-1.21737e-11,48629.2,21.0018], Tmin=(100,'K'), Tmax=(912.619,'K')), NASAPolynomial(coeffs=[9.99994,0.0171879,-7.21242e-06,1.31777e-09,-8.97962e-14,47099.1,-18.6756], Tmin=(912.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(ROOJ)"""),
)

species(
    label = 'C=C(C=C[CH]O)O[O](5094)',
    structure = SMILES('C=C([CH]C=CO)O[O]'),
    E0 = (58.8148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.992527,'amu*angstrom^2'), symmetry=1, barrier=(22.8202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993107,'amu*angstrom^2'), symmetry=1, barrier=(22.8335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993093,'amu*angstrom^2'), symmetry=1, barrier=(22.8332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99239,'amu*angstrom^2'), symmetry=1, barrier=(22.817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.707558,0.0859234,-9.1243e-05,4.61643e-08,-8.77711e-12,7258.33,31.6443], Tmin=(100,'K'), Tmax=(1438.94,'K')), NASAPolynomial(coeffs=[23.7118,0.00886819,-1.35505e-06,8.83176e-11,-2.15336e-15,1180.44,-91.7514], Tmin=(1438.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(58.8148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C(C=[C]CO)O[O](5095)',
    structure = SMILES('C=C(C=[C]CO)O[O]'),
    E0 = (205.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.318583,'amu*angstrom^2'), symmetry=1, barrier=(7.32485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318296,'amu*angstrom^2'), symmetry=1, barrier=(7.31825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318589,'amu*angstrom^2'), symmetry=1, barrier=(7.325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317651,'amu*angstrom^2'), symmetry=1, barrier=(7.30341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.393687,0.0860213,-0.000130565,1.09251e-07,-3.59846e-11,24804.3,31.0188], Tmin=(100,'K'), Tmax=(854.669,'K')), NASAPolynomial(coeffs=[9.66284,0.0311289,-1.4022e-05,2.58556e-09,-1.74193e-13,23640.3,-9.78433], Tmin=(854.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]O[O](1206)',
    structure = SMILES('C=[C]O[O]'),
    E0 = (337.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1685,370,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0010132,'amu*angstrom^2'), symmetry=1, barrier=(11.504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.15933,0.0239696,-4.52783e-05,4.8193e-08,-1.85927e-11,40635.7,16.5561], Tmin=(100,'K'), Tmax=(873.243,'K')), NASAPolynomial(coeffs=[2.07186,0.015386,-7.23299e-06,1.35893e-09,-9.18705e-14,41342.9,24.6165], Tmin=(873.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO) + radical(ROOJ)"""),
)

species(
    label = '[CH]=CCO(4313)',
    structure = SMILES('[CH]=CCO'),
    E0 = (103.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350],'cm^-1')),
        HinderedRotor(inertia=(0.558245,'amu*angstrom^2'), symmetry=1, barrier=(12.8352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558882,'amu*angstrom^2'), symmetry=1, barrier=(12.8498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25457,0.0421352,-6.08019e-05,5.22479e-08,-1.79007e-11,12500.9,14.1261], Tmin=(100,'K'), Tmax=(845.742,'K')), NASAPolynomial(coeffs=[5.58243,0.0195677,-8.66594e-06,1.60488e-09,-1.08868e-13,12182.2,0.0724052], Tmin=(845.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P)"""),
)

species(
    label = 'C=C([C]=CCO)O[O](5096)',
    structure = SMILES('C=C([C]=CCO)O[O]'),
    E0 = (166.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.395992,'amu*angstrom^2'), symmetry=1, barrier=(9.10464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395707,'amu*angstrom^2'), symmetry=1, barrier=(9.09808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395505,'amu*angstrom^2'), symmetry=1, barrier=(9.09344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39582,'amu*angstrom^2'), symmetry=1, barrier=(9.10069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.312175,0.0883892,-0.000136288,1.14518e-07,-3.73732e-11,20134.6,30.4344], Tmin=(100,'K'), Tmax=(884.521,'K')), NASAPolynomial(coeffs=[9.54774,0.0313123,-1.35288e-05,2.42331e-09,-1.59628e-13,19099.7,-9.60274], Tmin=(884.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C(C=CCO)O[O](5097)',
    structure = SMILES('[CH]=C(C=CCO)O[O]'),
    E0 = (214.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.386895,'amu*angstrom^2'), symmetry=1, barrier=(8.89548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386694,'amu*angstrom^2'), symmetry=1, barrier=(8.89086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.38598,'amu*angstrom^2'), symmetry=1, barrier=(8.87444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386781,'amu*angstrom^2'), symmetry=1, barrier=(8.89286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370097,0.0847162,-0.000121504,9.5944e-08,-3.0223e-11,25920,31.0345], Tmin=(100,'K'), Tmax=(831.654,'K')), NASAPolynomial(coeffs=[11.0622,0.0288534,-1.27459e-05,2.34663e-09,-1.58759e-13,24295,-17.6554], Tmin=(831.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'OCC=C[C]1COO1(5098)',
    structure = SMILES('OC[CH]C=C1COO1'),
    E0 = (-10.7163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773377,0.0588682,-2.19811e-05,-1.1828e-08,7.76693e-12,-1162.43,28.3837], Tmin=(100,'K'), Tmax=(1083.85,'K')), NASAPolynomial(coeffs=[15.5741,0.0280049,-1.21498e-05,2.35053e-09,-1.68727e-13,-5766.33,-50.6509], Tmin=(1083.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.7163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C1[CH]C(CO)OO1(5081)',
    structure = SMILES('C=C1[CH]C(CO)OO1'),
    E0 = (-99.5721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.947957,0.0499112,1.28103e-05,-5.4795e-08,2.47629e-11,-11850.6,26.2925], Tmin=(100,'K'), Tmax=(987.955,'K')), NASAPolynomial(coeffs=[17.5818,0.0234045,-8.95125e-06,1.73122e-09,-1.28875e-13,-17130.3,-63.8415], Tmin=(987.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.5721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJCO)"""),
)

species(
    label = 'CC(=CC=CO)O[O](5099)',
    structure = SMILES('CC(=CC=CO)O[O]'),
    E0 = (-88.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.71528,0.0867709,-9.21173e-05,4.7338e-08,-9.15442e-12,-10514.7,31.4464], Tmin=(100,'K'), Tmax=(1420.87,'K')), NASAPolynomial(coeffs=[22.9017,0.0106467,-1.57882e-06,8.36699e-11,-1.19284e-16,-16253.1,-87.3643], Tmin=(1420.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(8)(9)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=C=C=CCO(4301)',
    structure = SMILES('C=C=C=CCO'),
    E0 = (134.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,540,563.333,586.667,610,1970,2140,2750,2850,1437.5,1250,1305,750,350,180],'cm^-1')),
        HinderedRotor(inertia=(0.742305,'amu*angstrom^2'), symmetry=1, barrier=(17.0671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744095,'amu*angstrom^2'), symmetry=1, barrier=(17.1082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79694,0.0502438,-4.52937e-05,2.33275e-08,-5.11722e-12,16211.8,20.5527], Tmin=(100,'K'), Tmax=(1065.41,'K')), NASAPolynomial(coeffs=[8.1957,0.02622,-1.14701e-05,2.16263e-09,-1.50798e-13,14848.3,-10.723], Tmin=(1065.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O]OC1=CC(CO)C1(5100)',
    structure = SMILES('[O]OC1=CC(CO)C1'),
    E0 = (14.8289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854536,0.0605431,-3.92932e-05,7.90972e-09,1.26999e-12,1904.11,29.4607], Tmin=(100,'K'), Tmax=(1075.96,'K')), NASAPolynomial(coeffs=[14.0764,0.0254148,-9.87352e-06,1.7961e-09,-1.2439e-13,-1753,-39.0679], Tmin=(1075.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.8289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = 'OC[CH]C1CC12OO2(5101)',
    structure = SMILES('OC[CH]C1CC12OO2'),
    E0 = (64.0663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396312,0.0567995,1.43427e-05,-7.5065e-08,3.70783e-11,7855.76,26.8686], Tmin=(100,'K'), Tmax=(939.667,'K')), NASAPolynomial(coeffs=[24.7185,0.00993154,-1.29984e-06,2.10508e-10,-2.3418e-14,783.029,-102.27], Tmin=(939.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.0663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s1_3_3_ane) + radical(CCJCO)"""),
)

species(
    label = 'C=C1OO1(1417)',
    structure = SMILES('C=C1OO1'),
    E0 = (186.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,180,731.306,912.669,914.038,920.398,2195.41],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46953,0.00675937,1.90271e-05,-2.5454e-08,8.80673e-12,22401.8,9.99112], Tmin=(100,'K'), Tmax=(1065.74,'K')), NASAPolynomial(coeffs=[5.78132,0.0110046,-5.13514e-06,1.03739e-09,-7.63733e-14,21175.2,-4.75203], Tmin=(1065.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropane)"""),
)

species(
    label = 'CC1([C]=CCO)OO1(5102)',
    structure = SMILES('CC1([C]=CCO)OO1'),
    E0 = (40.5974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180,180,319.427,802.847,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125188,0.0644426,-3.75107e-06,-6.5603e-08,3.75746e-11,5041.83,29.0906], Tmin=(100,'K'), Tmax=(899.258,'K')), NASAPolynomial(coeffs=[26.707,0.00303889,3.86922e-06,-9.69454e-10,6.68474e-14,-2036.99,-109.106], Tmin=(899.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.5974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(Cds_S)"""),
)

species(
    label = 'CC1(C=[C]CO)OO1(5103)',
    structure = SMILES('CC1(C=[C]CO)OO1'),
    E0 = (40.5974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180,180,319.427,802.847,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150608,'amu*angstrom^2'), symmetry=1, barrier=(3.46278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.125188,0.0644426,-3.75107e-06,-6.5603e-08,3.75746e-11,5041.83,29.0906], Tmin=(100,'K'), Tmax=(899.258,'K')), NASAPolynomial(coeffs=[26.707,0.00303889,3.86922e-06,-9.69454e-10,6.68474e-14,-2036.99,-109.106], Tmin=(899.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.5974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(Cds_S)"""),
)

species(
    label = 'CC1(C=C[CH]O)OO1(5104)',
    structure = SMILES('CC1([CH]C=CO)OO1'),
    E0 = (-122.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.560698,0.0613679,5.77509e-05,-1.65104e-07,8.19058e-11,-14544.9,25.4915], Tmin=(100,'K'), Tmax=(896.408,'K')), NASAPolynomial(coeffs=[41.2231,-0.0196416,1.68695e-05,-3.48153e-09,2.34999e-13,-26272.3,-195.151], Tmin=(896.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(dioxirane) + radical(C=CCJCO)"""),
)

species(
    label = 'CC1(C=CC[O])OO1(5105)',
    structure = SMILES('CC1(C=CC[O])OO1'),
    E0 = (28.4607,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,826.904,1370,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155009,'amu*angstrom^2'), symmetry=1, barrier=(3.56396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155009,'amu*angstrom^2'), symmetry=1, barrier=(3.56396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155009,'amu*angstrom^2'), symmetry=1, barrier=(3.56396,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.486746,0.0575082,7.23483e-06,-6.774e-08,3.5522e-11,3567.95,28.1786], Tmin=(100,'K'), Tmax=(914.261,'K')), NASAPolynomial(coeffs=[23.2971,0.00979745,5.41906e-08,-1.88987e-10,1.08646e-14,-2779.9,-91.7289], Tmin=(914.261,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.4607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C=CC1([CH2])OO1(5106)',
    structure = SMILES('[CH2]C1([CH]C=C)OO1'),
    E0 = (300.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420873,0.0482424,4.9203e-05,-1.31144e-07,6.37438e-11,36257.3,24.1024], Tmin=(100,'K'), Tmax=(901.615,'K')), NASAPolynomial(coeffs=[32.0105,-0.010257,1.0692e-05,-2.23015e-09,1.48775e-13,27242.4,-143.431], Tmin=(901.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(CJCOOH) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1(C=CC[O])OO1(5107)',
    structure = SMILES('[CH2]C1(C=CC[O])OO1'),
    E0 = (242.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419086,0.0602446,-5.92581e-06,-5.52662e-08,3.17443e-11,29303.1,30.4172], Tmin=(100,'K'), Tmax=(907.549,'K')), NASAPolynomial(coeffs=[24.0219,0.0059265,1.68911e-06,-5.05349e-10,3.37287e-14,22971.7,-92.4419], Tmin=(907.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(CJCOOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]1OO1(1419)',
    structure = SMILES('[CH2][C]1OO1'),
    E0 = (356.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,610.356,610.368,610.389,610.391,610.392],'cm^-1')),
        HinderedRotor(inertia=(0.000452438,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05171,0.0090744,3.70704e-05,-6.57102e-08,2.96296e-11,42936.5,14.6158], Tmin=(100,'K'), Tmax=(887.932,'K')), NASAPolynomial(coeffs=[12.8666,-0.00339148,4.4954e-06,-9.83912e-10,6.79446e-14,39942,-38.6165], Tmin=(887.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(dioxirane) + radical(CJCOOH) + radical(Cs_P)"""),
)

species(
    label = '[CH]=CC1([CH2])OO1(5108)',
    structure = SMILES('[CH]=CC1([CH2])OO1'),
    E0 = (452.976,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,345.676,345.676,345.716,345.72,345.73,345.754,345.824],'cm^-1')),
        HinderedRotor(inertia=(0.348345,'amu*angstrom^2'), symmetry=1, barrier=(29.5407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34846,'amu*angstrom^2'), symmetry=1, barrier=(29.5388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22299,0.0358145,4.32753e-05,-1.1297e-07,5.57192e-11,54604.3,21.5162], Tmin=(100,'K'), Tmax=(889.712,'K')), NASAPolynomial(coeffs=[27.8901,-0.014986,1.24394e-05,-2.58378e-09,1.76643e-13,47124.6,-119.389], Tmin=(889.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cds_P) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1(C=C[CH]O)OO1(5109)',
    structure = SMILES('[CH2]C1([CH]C=CO)OO1'),
    E0 = (91.3581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.63166,0.0641409,4.44785e-05,-1.52521e-07,7.81034e-11,11190.3,27.7421], Tmin=(100,'K'), Tmax=(893.001,'K')), NASAPolynomial(coeffs=[41.9749,-0.0235586,1.8531e-05,-3.80421e-09,2.58388e-13,-531.936,-196.016], Tmin=(893.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.3581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(dioxirane) + radical(CJCOOH) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1([C]=CCO)OO1(5110)',
    structure = SMILES('[CH2]C1([C]=CCO)OO1'),
    E0 = (254.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180,180,265.425,851.531,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0572488,0.0671779,-1.68816e-05,-5.32181e-08,3.38631e-11,30776.9,31.3305], Tmin=(100,'K'), Tmax=(892.636,'K')), NASAPolynomial(coeffs=[27.4495,-0.000861837,5.52121e-06,-1.28983e-09,9.00431e-14,23707.1,-109.919], Tmin=(892.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(Cds_S) + radical(CJCOOH)"""),
)

species(
    label = '[CH2]C1(C=[C]CO)OO1(5111)',
    structure = SMILES('[CH2]C1(C=[C]CO)OO1'),
    E0 = (254.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,1685,370,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,180,180,180,265.425,851.531,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148337,'amu*angstrom^2'), symmetry=1, barrier=(3.41056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0572488,0.0671779,-1.68816e-05,-5.32181e-08,3.38631e-11,30776.9,31.3305], Tmin=(100,'K'), Tmax=(892.636,'K')), NASAPolynomial(coeffs=[27.4495,-0.000861837,5.52121e-06,-1.28983e-09,9.00431e-14,23707.1,-109.919], Tmin=(892.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(CJCOOH) + radical(Cds_S)"""),
)

species(
    label = 'OCC1[CH]C2(C1)OO2(5112)',
    structure = SMILES('OCC1[CH]C2(C1)OO2'),
    E0 = (47.8804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507421,0.0566146,8.63827e-06,-6.28995e-08,3.0966e-11,5902.8,26.6269], Tmin=(100,'K'), Tmax=(953.958,'K')), NASAPolynomial(coeffs=[22.3126,0.0149795,-4.19251e-06,7.8494e-10,-6.3149e-14,-523.217,-89.4182], Tmin=(953.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.8804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + polycyclic(s1_3_4_ane) + radical(CCJCOOH)"""),
)

species(
    label = 'O(S)(1202)',
    structure = SMILES('O'),
    E0 = (432.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,51997.4,2.99252], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,51997.4,2.99252], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.331,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OCC=CC[C]1OO1(5113)',
    structure = SMILES('OCC=CC[C]1OO1'),
    E0 = (33.7157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759975,0.0604347,-2.72934e-05,-1.36813e-08,1.15684e-11,4181.56,31.367], Tmin=(100,'K'), Tmax=(947.113,'K')), NASAPolynomial(coeffs=[16.2897,0.0211385,-6.69679e-06,1.1305e-09,-7.79645e-14,60.6839,-48.9363], Tmin=(947.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.7157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(Cs_P)"""),
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
    label = 'OCC=C[C]1OO1(5114)',
    structure = SMILES('OCC=C[C]1OO1'),
    E0 = (56.3612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43077,0.0482972,-2.39038e-05,-9.40085e-09,9.1136e-12,6878.73,26.0474], Tmin=(100,'K'), Tmax=(931.188,'K')), NASAPolynomial(coeffs=[13.9207,0.015457,-4.52718e-06,7.27649e-10,-4.92289e-14,3650.34,-38.1637], Tmin=(931.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.3612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(Cs_P)"""),
)

species(
    label = '[CH]C1(C=CCO)OO1(5115)',
    structure = SMILES('[CH]C1(C=CCO)OO1'),
    E0 = (250.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.10294,0.0629213,1.50602e-06,-7.53051e-08,4.22852e-11,30346.6,29.7014], Tmin=(100,'K'), Tmax=(897.321,'K')), NASAPolynomial(coeffs=[28.8052,-0.00255994,6.54732e-06,-1.47163e-09,1.00686e-13,22680.8,-119.674], Tmin=(897.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(dioxirane) + radical(CCJ2_triplet)"""),
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
    E0 = (76.7209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (34.1631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (65.3564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (58.6658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (106.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (86.2437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (27.5324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (113.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (266.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (404.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (376.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (273.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (441.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (382.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (426.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (93.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (25.5241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (131.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (146.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (89.7795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (71.6765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (64.0663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (297.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (192.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (84.9059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (94.7498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (121.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (328.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (454.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (460.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (425.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (306.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (466.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (466.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (139.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (261.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (279.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (437.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (462.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'C5H7O(237)(236)'],
    products = ['S(239)(238)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.2234e+15,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/Cd;Y_rad] for rate rule [Cd_rad/Cd;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(239)(238)'],
    products = ['S(243)(242)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.1786e+12,'s^-1'), n=-0.0937843, Ea=(66.7967,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secDe;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secDe;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['S(239)(238)'],
    products = ['C=C1OOC1[CH]CO(5085)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.65e+06,'s^-1'), n=1.3, Ea=(97.99,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 94.9 to 98.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([C]=CCO)OO(5086)'],
    products = ['S(239)(238)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_S(Cd)SS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C(C=CCO)OO(5087)'],
    products = ['S(239)(238)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C(C=[C]CO)OO(5088)'],
    products = ['S(239)(238)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['S(239)(238)'],
    products = ['C=C(C=C[CH]O)OO(5089)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1098,'s^-1'), n=2.21, Ea=(60.1659,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H;O_rad_out;Cs_H_out_H/NonDeO] for rate rule [R6H_RSSMS;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C(C=CC[O])OO(5090)'],
    products = ['S(239)(238)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(627096,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;O_rad_out;XH_out] for rate rule [R7H;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OH(5)(5)', '[CH2]C=CC(=C)O[O](5091)'],
    products = ['S(239)(238)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.7e+13,'cm^3/(mol*s)','+|-',1e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(400,'K'), comment="""Estimated using template [O_pri_rad;C_pri_rad] for rate rule [O_pri_rad;C_rad/H2/Cd]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)(3)', 'C=C(C=CC[O])O[O](5092)'],
    products = ['S(239)(238)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2OH(22)(23)', '[CH]=CC(=C)O[O](5093)'],
    products = ['S(239)(238)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_rad] for rate rule [C_rad/H2/O;Cd_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)(3)', 'C=C(C=C[CH]O)O[O](5094)'],
    products = ['S(239)(238)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)(3)', 'C=C(C=[C]CO)O[O](5095)'],
    products = ['S(239)(238)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]O[O](1206)', '[CH]=CCO(4313)'],
    products = ['S(239)(238)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.68887e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_rad/NonDe] + [Cd_pri_rad;Cd_rad] for rate rule [Cd_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)(3)', 'C=C([C]=CCO)O[O](5096)'],
    products = ['S(239)(238)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)(3)', '[CH]=C(C=CCO)O[O](5097)'],
    products = ['S(239)(238)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['S(239)(238)'],
    products = ['OCC=C[C]1COO1(5098)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['S(239)(238)'],
    products = ['C=C1[CH]C(CO)OO1(5081)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.94e+07,'s^-1'), n=0.93, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(239)(238)'],
    products = ['CC(=CC=CO)O[O](5099)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.8582e+08,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH_end;CdH2_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(239)(238)'],
    products = ['HO2(8)(9)', 'C=C=C=CCO(4301)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction21',
    reactants = ['S(239)(238)'],
    products = ['[O]OC1=CC(CO)C1(5100)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH2_1;CdH(C)_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)(4)', 'S(244)(243)'],
    products = ['S(239)(238)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['S(243)(242)'],
    products = ['OC[CH]C1CC12OO2(5101)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.68e+09,'s^-1'), n=0.84, Ea=(47.3481,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 338 used for R4_S_D;doublebond_intra_HNd_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_HNd_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 46.4 to 47.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C1OO1(1417)', '[CH]=CCO(4313)'],
    products = ['S(243)(242)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00977588,'m^3/(mol*s)'), n=2.40996, Ea=(8.29121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ-H] for rate rule [Cds-OsOs_Cds;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC1([C]=CCO)OO1(5102)'],
    products = ['S(243)(242)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC1(C=[C]CO)OO1(5103)'],
    products = ['S(243)(242)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['S(243)(242)'],
    products = ['CC1(C=C[CH]O)OO1(5104)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(423178,'s^-1'), n=1.77, Ea=(78.0316,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;Cs_H_out] for rate rule [R5H_SSMS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC1(C=CC[O])OO1(5105)'],
    products = ['S(243)(242)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6H;O_rad_out;Cs_H_out_2H] for rate rule [R6H_RSMSR;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['OH(5)(5)', '[CH2]C=CC1([CH2])OO1(5106)'],
    products = ['S(243)(242)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.7e+13,'cm^3/(mol*s)','+|-',1e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(400,'K'), comment="""Estimated using template [O_pri_rad;C_pri_rad] for rate rule [O_pri_rad;C_rad/H2/Cd]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)(3)', '[CH2]C1(C=CC[O])OO1(5107)'],
    products = ['S(243)(242)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]1OO1(1419)', '[CH]=CCO(4313)'],
    products = ['S(243)(242)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.88428e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_pri_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2OH(22)(23)', '[CH]=CC1([CH2])OO1(5108)'],
    products = ['S(243)(242)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.41031e+09,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Cd_rad] for rate rule [C_rad/H2/O;Cd_pri_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)(3)', '[CH2]C1(C=C[CH]O)OO1(5109)'],
    products = ['S(243)(242)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;H_rad] for rate rule [C_rad/H/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)(3)', '[CH2]C1([C]=CCO)OO1(5110)'],
    products = ['S(243)(242)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)(3)', '[CH2]C1(C=[C]CO)OO1(5111)'],
    products = ['S(243)(242)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['S(243)(242)'],
    products = ['OCC1[CH]C2(C1)OO2(5112)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['S(244)(243)', 'O(S)(1202)'],
    products = ['S(243)(242)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.96556e+06,'m^3/(mol*s)'), n=0, Ea=(0.08368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [o_atom_singlet;multiplebond] for rate rule [o_atom_singlet;mb_carbonyl]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction38',
    reactants = ['OCC=CC[C]1OO1(5113)'],
    products = ['S(243)(242)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(17)(18)', 'OCC=C[C]1OO1(5114)'],
    products = ['S(243)(242)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)(3)', '[CH]C1(C=CCO)OO1(5115)'],
    products = ['S(243)(242)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '351',
    isomers = [
        'S(239)(238)',
        'S(243)(242)',
    ],
    reactants = [
        ('O2(2)(2)', 'C5H7O(237)(236)'),
        ('O(4)(4)', 'S(244)(243)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '351',
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

