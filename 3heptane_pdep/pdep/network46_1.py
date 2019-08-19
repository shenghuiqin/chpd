species(
    label = 'C[CH]CCOO(216)',
    structure = SMILES('C[CH]CCOO'),
    E0 = (-38.1802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,217.649,3236.35],'cm^-1')),
        HinderedRotor(inertia=(0.120436,'amu*angstrom^2'), symmetry=1, barrier=(4.12955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120396,'amu*angstrom^2'), symmetry=1, barrier=(4.13001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120996,'amu*angstrom^2'), symmetry=1, barrier=(4.13046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63323,'amu*angstrom^2'), symmetry=1, barrier=(55.1493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62971,'amu*angstrom^2'), symmetry=1, barrier=(55.16,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3638.29,'J/mol'), sigma=(6.34008,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.29 K, Pc=32.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95817,0.0528385,-1.6212e-05,-5.24154e-08,5.55562e-11,-4526.51,23.2677], Tmin=(100,'K'), Tmax=(485.277,'K')), NASAPolynomial(coeffs=[4.83968,0.0414029,-1.89326e-05,3.61988e-09,-2.53625e-13,-4951.19,9.95538], Tmin=(485.277,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-38.1802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-2.38914e-13,3.12709e-16,-1.33367e-19,1.7499e-23,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(4383.16,'K')), NASAPolynomial(coeffs=[2.50003,-3.04997e-08,1.01101e-11,-1.48797e-15,8.20356e-20,25472.7,-0.459785], Tmin=(4383.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CC=CCOO(237)',
    structure = SMILES('CC=CCOO'),
    E0 = (-111.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,247.061],'cm^-1')),
        HinderedRotor(inertia=(0.181143,'amu*angstrom^2'), symmetry=1, barrier=(7.82432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180553,'amu*angstrom^2'), symmetry=1, barrier=(7.82423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.840396,'amu*angstrom^2'), symmetry=1, barrier=(36.4233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31576,'amu*angstrom^2'), symmetry=1, barrier=(13.7024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80981,0.0510568,-3.39645e-05,1.18362e-08,-1.78618e-12,-13338.4,21.3408], Tmin=(100,'K'), Tmax=(1425.39,'K')), NASAPolynomial(coeffs=[8.40367,0.0325524,-1.44911e-05,2.72818e-09,-1.88681e-13,-15218.1,-12.8078], Tmin=(1425.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH)"""),
)

species(
    label = 'C=CCCOO(238)',
    structure = SMILES('C=CCCOO'),
    E0 = (-104.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,323.996,323.997],'cm^-1')),
        HinderedRotor(inertia=(0.111179,'amu*angstrom^2'), symmetry=1, barrier=(8.28184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331825,'amu*angstrom^2'), symmetry=1, barrier=(24.7181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111178,'amu*angstrom^2'), symmetry=1, barrier=(8.28184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419393,'amu*angstrom^2'), symmetry=1, barrier=(31.2412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37317,0.0554599,-4.0154e-05,1.50002e-08,-2.32081e-12,-12513.7,22.9714], Tmin=(100,'K'), Tmax=(1476.89,'K')), NASAPolynomial(coeffs=[11.7099,0.0274642,-1.17202e-05,2.16524e-09,-1.48175e-13,-15567,-30.9278], Tmin=(1476.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C3H6(59)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817957,3.34737e-05,-4.36194e-08,1.58214e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.27621e-14,-487.137,-4.54465], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]OO(239)',
    structure = SMILES('[CH2]OO'),
    E0 = (54.8878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.878705,'amu*angstrom^2'), symmetry=1, barrier=(20.2032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0333,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85643,0.0201434,-7.5271e-06,-9.6693e-09,6.63648e-12,6647.28,10.9504], Tmin=(100,'K'), Tmax=(928.338,'K')), NASAPolynomial(coeffs=[9.91701,0.00263656,-1.08425e-07,-1.03822e-11,-4.89454e-16,4779.82,-25.5853], Tmin=(928.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.8878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CsJOOH)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48579,0.001334,-4.70054e-06,5.64393e-09,-2.06324e-12,3411.96,1.99789], Tmin=(100,'K'), Tmax=(1005.24,'K')), NASAPolynomial(coeffs=[2.88226,0.00103867,-2.35641e-07,1.40204e-11,6.3479e-16,3669.56,5.59047], Tmin=(1005.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C[CH]CC[O](240)',
    structure = SMILES('C[CH]CC[O]'),
    E0 = (123.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,1820.69,1824.25],'cm^-1')),
        HinderedRotor(inertia=(0.228416,'amu*angstrom^2'), symmetry=1, barrier=(5.25173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227185,'amu*angstrom^2'), symmetry=1, barrier=(5.22342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00222396,'amu*angstrom^2'), symmetry=1, barrier=(5.25126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44242,0.0298368,-6.5991e-06,-1.96412e-09,6.28692e-13,14825.4,18.071], Tmin=(100,'K'), Tmax=(2180.07,'K')), NASAPolynomial(coeffs=[16.5628,0.0192345,-8.5728e-06,1.47379e-09,-9.05793e-14,5903.6,-62.7944], Tmin=(2180.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo library: CBS_QB3_1dHR + radical(RCCJC) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C[CH]C(86)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2031.22],'cm^-1')),
        HinderedRotor(inertia=(0.244949,'amu*angstrom^2'), symmetry=1, barrier=(5.63185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192353,'amu*angstrom^2'), symmetry=1, barrier=(5.63196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244945,'amu*angstrom^2'), symmetry=1, barrier=(5.63176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51473e-06,4.19257e-10,1.9052e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.58,'K')), NASAPolynomial(coeffs=[12.4235,0.0182237,-7.06296e-06,1.16765e-09,-7.11789e-14,25091.2,-39.6233], Tmin=(2154.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'HO2(10)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.8,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02957,-0.00263999,1.52235e-05,-1.71679e-08,6.26771e-12,322.677,4.84424], Tmin=(100,'K'), Tmax=(923.901,'K')), NASAPolynomial(coeffs=[4.1513,0.00191152,-4.11308e-07,6.35038e-11,-4.86452e-15,83.4341,3.09359], Tmin=(923.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C[CH]CCO[O](219)',
    structure = SMILES('C[CH]CCO[O]'),
    E0 = (113.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2117.31],'cm^-1')),
        HinderedRotor(inertia=(0.210814,'amu*angstrom^2'), symmetry=1, barrier=(4.84703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210864,'amu*angstrom^2'), symmetry=1, barrier=(4.84818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210865,'amu*angstrom^2'), symmetry=1, barrier=(4.8482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210852,'amu*angstrom^2'), symmetry=1, barrier=(4.8479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53015,0.0598849,-7.45192e-05,6.46902e-08,-2.34831e-11,13773.5,25.4288], Tmin=(100,'K'), Tmax=(824.028,'K')), NASAPolynomial(coeffs=[3.73803,0.0388311,-1.73791e-05,3.23958e-09,-2.21489e-13,13760.5,17.3338], Tmin=(824.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(ROOJ)"""),
)

species(
    label = 'C3H6(T)(82)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26859e-06,-1.44876e-08,5.7495e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93907,0.0171892,-6.69154e-06,1.21547e-09,-8.39798e-14,33151.2,-4.14876], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH]COO(241)',
    structure = SMILES('C[CH][CH]COO'),
    E0 = (162.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1946.26,1946.28],'cm^-1')),
        HinderedRotor(inertia=(0.340569,'amu*angstrom^2'), symmetry=1, barrier=(10.2209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57712,'amu*angstrom^2'), symmetry=1, barrier=(47.3369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57705,'amu*angstrom^2'), symmetry=1, barrier=(47.3367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00380237,'amu*angstrom^2'), symmetry=1, barrier=(10.221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34057,'amu*angstrom^2'), symmetry=1, barrier=(10.221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29793,0.0654316,-8.55343e-05,7.4494e-08,-2.70526e-11,19603.8,27.8082], Tmin=(100,'K'), Tmax=(804.143,'K')), NASAPolynomial(coeffs=[4.60467,0.0391518,-1.81748e-05,3.44673e-09,-2.38096e-13,19389.9,14.5524], Tmin=(804.143,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH]OO(242)',
    structure = SMILES('C[CH]C[CH]OO'),
    E0 = (150.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,215.112,1603.05],'cm^-1')),
        HinderedRotor(inertia=(0.168506,'amu*angstrom^2'), symmetry=1, barrier=(5.53198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303382,'amu*angstrom^2'), symmetry=1, barrier=(5.53226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74463,'amu*angstrom^2'), symmetry=1, barrier=(57.2691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168498,'amu*angstrom^2'), symmetry=1, barrier=(5.53193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74363,'amu*angstrom^2'), symmetry=1, barrier=(57.2691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.266,0.0662665,-8.68775e-05,7.58355e-08,-2.75565e-11,18181.7,25.7672], Tmin=(100,'K'), Tmax=(806.28,'K')), NASAPolynomial(coeffs=[4.55225,0.0397217,-1.84402e-05,3.49455e-09,-2.4123e-13,17984.7,12.6849], Tmin=(806.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2][CH]CCOO(243)',
    structure = SMILES('[CH2][CH]CCOO'),
    E0 = (167.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3025,407.5,1350,352.5,1842.9,1842.94],'cm^-1')),
        HinderedRotor(inertia=(0.206666,'amu*angstrom^2'), symmetry=1, barrier=(9.2458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206508,'amu*angstrom^2'), symmetry=1, barrier=(9.24398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206535,'amu*angstrom^2'), symmetry=1, barrier=(9.24525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10045,'amu*angstrom^2'), symmetry=1, barrier=(49.2346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0204282,'amu*angstrom^2'), symmetry=1, barrier=(49.2327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33189,0.0638614,-7.90214e-05,6.6096e-08,-2.36035e-11,20184.5,27.1692], Tmin=(100,'K'), Tmax=(784.688,'K')), NASAPolynomial(coeffs=[5.06478,0.0382979,-1.76624e-05,3.35252e-09,-2.32268e-13,19799.9,11.3473], Tmin=(784.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'CC[CH]COO(215)',
    structure = SMILES('CC[CH]COO'),
    E0 = (-32.2116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,255.525,3394.99],'cm^-1')),
        HinderedRotor(inertia=(0.262752,'amu*angstrom^2'), symmetry=1, barrier=(12.1744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00258184,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262746,'amu*angstrom^2'), symmetry=1, barrier=(12.1743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826362,'amu*angstrom^2'), symmetry=1, barrier=(38.2897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.826398,'amu*angstrom^2'), symmetry=1, barrier=(38.2897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39056,0.0615982,-5.23509e-05,2.61861e-08,-5.81249e-12,-3783.71,24.7522], Tmin=(100,'K'), Tmax=(1025.75,'K')), NASAPolynomial(coeffs=[7.56867,0.0375058,-1.71188e-05,3.28734e-09,-2.31432e-13,-5051.13,-5.21059], Tmin=(1025.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.2116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]CCCOO(217)',
    structure = SMILES('[CH2]CCCOO'),
    E0 = (-27.3802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,364.375,2004.58],'cm^-1')),
        HinderedRotor(inertia=(0.250731,'amu*angstrom^2'), symmetry=1, barrier=(23.6234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0760511,'amu*angstrom^2'), symmetry=1, barrier=(7.16522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0760487,'amu*angstrom^2'), symmetry=1, barrier=(7.1652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106701,'amu*angstrom^2'), symmetry=1, barrier=(10.0531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.457045,'amu*angstrom^2'), symmetry=1, barrier=(43.0627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38391,0.0606318,-4.84656e-05,2.17897e-08,-4.27182e-12,-3201.55,24.2498], Tmin=(100,'K'), Tmax=(1154.31,'K')), NASAPolynomial(coeffs=[8.58861,0.0356657,-1.60229e-05,3.05259e-09,-2.13777e-13,-4864.85,-11.5427], Tmin=(1154.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.3802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ)"""),
)

species(
    label = 'CCC[CH]OO(214)',
    structure = SMILES('CCC[CH]OO'),
    E0 = (-44.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,2546.25],'cm^-1')),
        HinderedRotor(inertia=(1.89684,'amu*angstrom^2'), symmetry=1, barrier=(43.6121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.542507,'amu*angstrom^2'), symmetry=1, barrier=(12.4733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167613,'amu*angstrom^2'), symmetry=1, barrier=(12.5115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0590487,'amu*angstrom^2'), symmetry=1, barrier=(43.6359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.544224,'amu*angstrom^2'), symmetry=1, barrier=(12.5128,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35904,0.0624115,-5.35472e-05,2.72539e-08,-6.17119e-12,-5205.8,22.7109], Tmin=(100,'K'), Tmax=(1006.82,'K')), NASAPolynomial(coeffs=[7.45804,0.0381808,-1.74476e-05,3.3506e-09,-2.35878e-13,-6433.92,-6.75477], Tmin=(1006.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH)"""),
)

species(
    label = 'CCCCO[O](107)',
    structure = SMILES('CCCCO[O]'),
    E0 = (-80.6218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,231.975,232.16],'cm^-1')),
        HinderedRotor(inertia=(0.157964,'amu*angstrom^2'), symmetry=1, barrier=(6.02847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156946,'amu*angstrom^2'), symmetry=1, barrier=(6.03016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157856,'amu*angstrom^2'), symmetry=1, barrier=(6.0287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157386,'amu*angstrom^2'), symmetry=1, barrier=(6.03044,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3638.29,'J/mol'), sigma=(6.34008,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=568.29 K, Pc=32.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56798,0.0567641,-4.41692e-05,2.06264e-08,-4.35462e-12,-9611.76,22.5651], Tmin=(100,'K'), Tmax=(1060.67,'K')), NASAPolynomial(coeffs=[6.81045,0.0369937,-1.62098e-05,3.05289e-09,-2.1253e-13,-10723.9,-3.03551], Tmin=(1060.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-80.6218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ)"""),
)

species(
    label = 'CC1CCO1(244)',
    structure = SMILES('CC1CCO1'),
    E0 = (-132.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60711,0.0163705,6.56447e-05,-9.99085e-08,4.19647e-11,-15901.6,14.3751], Tmin=(100,'K'), Tmax=(884.966,'K')), NASAPolynomial(coeffs=[10.4298,0.0189448,-3.01343e-06,2.48189e-10,-1.18933e-14,-18771.5,-30.8012], Tmin=(884.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]C(C)COO(245)',
    structure = SMILES('[CH2]C(C)COO'),
    E0 = (-36.5857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,296.018],'cm^-1')),
        HinderedRotor(inertia=(0.00192325,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100174,'amu*angstrom^2'), symmetry=1, barrier=(6.22901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100174,'amu*angstrom^2'), symmetry=1, barrier=(6.22891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39298,'amu*angstrom^2'), symmetry=1, barrier=(86.625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39264,'amu*angstrom^2'), symmetry=1, barrier=(86.6251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12376,0.0620874,-5.11548e-05,2.37099e-08,-4.62679e-12,-4295.76,24.8951], Tmin=(100,'K'), Tmax=(1197.69,'K')), NASAPolynomial(coeffs=[10.3239,0.031361,-1.26727e-05,2.28973e-09,-1.55655e-13,-6499.55,-21.1501], Tmin=(1197.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.5857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl)"""),
)

species(
    label = 'CC(O)CC[O](246)',
    structure = SMILES('CC(O)CC[O]'),
    E0 = (-250.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (89.1131,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07701,0.0502187,-9.62828e-06,-6.37655e-08,6.42015e-11,-30114.8,22.5732], Tmin=(100,'K'), Tmax=(472.853,'K')), NASAPolynomial(coeffs=[4.62326,0.0408783,-1.8697e-05,3.58058e-09,-2.51195e-13,-30491.9,10.754], Tmin=(472.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)"""),
)

species(
    label = 'CH2CH2OOH(55)',
    structure = SMILES('[CH2]COO'),
    E0 = (34.3631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.257945,'amu*angstrom^2'), symmetry=1, barrier=(5.93066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258128,'amu*angstrom^2'), symmetry=1, barrier=(5.93487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141801,'amu*angstrom^2'), symmetry=1, barrier=(16.6288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0599,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44846,0.0313149,-2.37472e-05,9.50342e-09,-1.55883e-12,4191.03,17.4641], Tmin=(100,'K'), Tmax=(1426.94,'K')), NASAPolynomial(coeffs=[8.63688,0.0139676,-5.51186e-06,9.83892e-10,-6.62178e-14,2424.92,-14.5917], Tmin=(1426.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.3631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""CH2CH2OOH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CHCH3(T)(80)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.415,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438701,'amu*angstrom^2'), symmetry=1, barrier=(26.7686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8236,-0.000909219,3.21369e-05,-3.73465e-08,1.33084e-11,41371.4,7.10957], Tmin=(100,'K'), Tmax=(960.825,'K')), NASAPolynomial(coeffs=[4.30495,0.00943054,-3.2755e-06,5.95101e-10,-4.2729e-14,40709.1,1.84155], Tmin=(960.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH3(18)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH]CCOO(247)',
    structure = SMILES('[CH]CCOO'),
    E0 = (247.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,320.321,320.324,1352.47,1352.47],'cm^-1')),
        HinderedRotor(inertia=(0.133672,'amu*angstrom^2'), symmetry=1, barrier=(9.73435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164293,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133695,'amu*angstrom^2'), symmetry=1, barrier=(9.73433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164264,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92275,0.0455979,-4.03983e-05,1.93855e-08,-3.83586e-12,29821.7,20.5537], Tmin=(100,'K'), Tmax=(1193.07,'K')), NASAPolynomial(coeffs=[9.4651,0.0203106,-8.60526e-06,1.62007e-09,-1.13197e-13,28022,-17.1652], Tmin=(1193.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]CCOO(248)',
    structure = SMILES('C[C]CCOO'),
    E0 = (215.589,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2489,0.063427,-5.89518e-05,3.08359e-08,-6.84941e-12,26026,22.858], Tmin=(100,'K'), Tmax=(1053.34,'K')), NASAPolynomial(coeffs=[9.4376,0.0323308,-1.46694e-05,2.80924e-09,-1.97541e-13,24300.9,-17.0733], Tmin=(1053.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49898e-06,-1.43375e-09,2.58635e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.04,'K')), NASAPolynomial(coeffs=[2.97591,0.0016414,-7.19719e-07,1.25377e-10,-7.91522e-15,-1025.85,5.53754], Tmin=(1817.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
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
    E0 = (105.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (109.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (97.9762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (151.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (258.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (325.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (339.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (374.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (362.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (378.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (126.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (120.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (128.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-4.59852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (37.6757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (123.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (73.1138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (378.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (383.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (427.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'CC=CCOO(237)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', 'C=CCCOO(238)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C3H6(59)', '[CH2]OO(239)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.381e+06,'cm^3/(mol*s)'), n=1.76, Ea=(37.1121,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2781 used for Cds-HH_Cds-Cs\H3/H;CsJ-OsHH
Exact match found for rate rule [Cds-HH_Cds-Cs\H3/H;CsJ-OsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OH(5)', 'C[CH]CC[O](240)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 103 used for O_pri_rad;O_rad/NonDe
Exact match found for rate rule [O_pri_rad;O_rad/NonDe]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C[CH]C(86)', 'HO2(10)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.58895e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [O_rad;C_rad/H2/Cs] + [O_rad/NonDe;Cs_rad] for rate rule [O_rad/NonDe;C_rad/H2/Cs]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C[CH]CCO[O](219)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C3H6(T)(82)', '[CH2]OO(239)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.40574e+06,'m^3/(mol*s)'), n=0.104005, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;Y_rad] for rate rule [C_rad/H2/O;Y_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'C[CH][CH]COO(241)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'C[CH]C[CH]OO(242)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CsO;Y_rad] for rate rule [C_rad/H/CsO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH]CCOO(243)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC[CH]COO(215)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CCCOO(217)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCC[CH]OO(214)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.85927e+09,'s^-1'), n=0.875, Ea=(172.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_1H;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CCCCO[O](107)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.788e+07,'s^-1'), n=1.26, Ea=(76.0233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 260 C4H9O2-13 <=> C4H9O2-14 in intra_H_migration/training
This reaction matched rate rule [R5H_SSSS_OCC;O_rad_out;Cs_H_out_H/NonDeC]
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]CCOO(216)'],
    products = ['OH(5)', 'CC1CCO1(244)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 9 used for R3OO_SS;C_rad/H/NonDeC_intra;OOH
Exact match found for rate rule [R3OO_SS;C_rad/H/NonDeC_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)COO(245)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH]CCOO(216)'],
    products = ['CC(O)CC[O](246)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.88e+10,'s^-1'), n=0, Ea=(111.294,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Matched reaction 4 C4H9O2 <=> C4H9O2-2 in intra_OH_migration/training
This reaction matched rate rule [R3OOH_SS;C_rad_out_H/NonDeC]
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH2CH2OOH(55)', 'CHCH3(T)(80)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH3(18)', '[CH]CCOO(247)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', 'C[C]CCOO(248)'],
    products = ['C[CH]CCOO(216)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '46',
    isomers = [
        'C[CH]CCOO(216)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '46',
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

