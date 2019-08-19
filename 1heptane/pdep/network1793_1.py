species(
    label = 'C=[C]OC([O])(CC)C(=C)O(8847)',
    structure = SMILES('C=[C]OC([O])(CC)C(=C)O'),
    E0 = (-96.7784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.3398,0.129292,-0.000139987,7.12127e-08,-1.33972e-11,-11347.9,45.1262], Tmin=(100,'K'), Tmax=(1507.99,'K')), NASAPolynomial(coeffs=[34.2002,0.00851664,1.23592e-06,-5.42778e-10,4.4204e-14,-20259.5,-143.411], Tmin=(1507.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.7784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(O)C(=O)CC(4626)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,332.629,332.633],'cm^-1')),
        HinderedRotor(inertia=(0.152178,'amu*angstrom^2'), symmetry=1, barrier=(11.9524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152152,'amu*angstrom^2'), symmetry=1, barrier=(11.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152157,'amu*angstrom^2'), symmetry=1, barrier=(11.9534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152209,'amu*angstrom^2'), symmetry=1, barrier=(11.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C(O)=C([O])CC(4557)',
    structure = SMILES('[CH2]C(O)=C([O])CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,402.631,402.631],'cm^-1')),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100034,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4303.36,'J/mol'), sigma=(7.05412,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.18 K, Pc=27.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1(O)OC1(CC)O[C]=C(9486)',
    structure = SMILES('[CH2]C1(O)OC1(CC)O[C]=C'),
    E0 = (-60.7774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35942,0.102623,-2.42071e-06,-1.24933e-07,7.29773e-11,-7045.37,36.9647], Tmin=(100,'K'), Tmax=(883.186,'K')), NASAPolynomial(coeffs=[46.8119,-0.0143023,1.65192e-05,-3.62437e-09,2.53886e-13,-19856.1,-217.505], Tmin=(883.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.7774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CJO) + radical(CJC(O)2C)"""),
)

species(
    label = '[CH2]C1(O)C(=C)OC1([O])CC(9487)',
    structure = SMILES('[CH2]C1(O)C(=C)OC1([O])CC'),
    E0 = (-103.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.24596,0.0929469,-3.06201e-05,-5.16118e-08,3.44474e-11,-12274.5,32.9656], Tmin=(100,'K'), Tmax=(907.619,'K')), NASAPolynomial(coeffs=[30.4149,0.0140412,-4.13096e-07,-2.01565e-10,1.44416e-14,-20518.9,-130.468], Tmin=(907.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=CC(C)(O)CJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'C#COC([O])(CC)C(=C)O(9488)',
    structure = SMILES('C#COC([O])(CC)C(=C)O'),
    E0 = (-127.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,947.941,1256.32,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159203,'amu*angstrom^2'), symmetry=1, barrier=(3.66039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.145,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.57124,0.133475,-0.00015048,7.80808e-08,-1.48836e-11,-14997.2,40.2816], Tmin=(100,'K'), Tmax=(1494.22,'K')), NASAPolynomial(coeffs=[36.3438,0.00345588,3.2996e-06,-9.0665e-10,6.80342e-14,-24339.3,-159.661], Tmin=(1494.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CC(C)(O)OJ)"""),
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
    label = 'C=[C]OC(=O)C(=C)O(9489)',
    structure = SMILES('C=[C]OC(=O)C(=C)O'),
    E0 = (-127.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,180,1600,1852.21,2658.74,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161871,'amu*angstrom^2'), symmetry=1, barrier=(3.72173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161871,'amu*angstrom^2'), symmetry=1, barrier=(3.72173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161871,'amu*angstrom^2'), symmetry=1, barrier=(3.72173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161871,'amu*angstrom^2'), symmetry=1, barrier=(3.72173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263466,0.0850513,-0.000116002,8.00385e-08,-2.16363e-11,-15158.6,29.2576], Tmin=(100,'K'), Tmax=(909.981,'K')), NASAPolynomial(coeffs=[15.0146,0.0202069,-9.10826e-06,1.7233e-09,-1.19776e-13,-17843.1,-40.516], Tmin=(909.981,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'CH2COH(99)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]OC(=O)CC(9490)',
    structure = SMILES('C=[C]OC(=O)CC'),
    E0 = (-135.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47841,0.0601081,-5.41933e-05,2.90853e-08,-6.92699e-12,-16239,24.9636], Tmin=(100,'K'), Tmax=(963.788,'K')), NASAPolynomial(coeffs=[7.25513,0.036133,-1.68793e-05,3.27466e-09,-2.31862e-13,-17352.6,-2.69266], Tmin=(963.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC([O])(CC)C(=C)O(9491)',
    structure = SMILES('[CH]=COC([O])(CC)C(=C)O'),
    E0 = (-89.4264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1021.56,1186.92,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.14344,0.137041,-0.000150111,7.55688e-08,-1.38736e-11,-10426.1,45.2957], Tmin=(100,'K'), Tmax=(1582.89,'K')), NASAPolynomial(coeffs=[37.2143,0.00330179,4.32209e-06,-1.13871e-09,8.39025e-14,-19857.6,-161.659], Tmin=(1582.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.4264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'C=[C]OC(O)([CH]C)C(=C)O(9492)',
    structure = SMILES('C=[C]OC(O)([CH]C)C(=C)O'),
    E0 = (-129.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,180,180,180,180,988.38,1208.86,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159505,'amu*angstrom^2'), symmetry=1, barrier=(3.66733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.123,0.107036,-4.83206e-05,-5.07378e-08,3.79261e-11,-15379.8,42.3682], Tmin=(100,'K'), Tmax=(913.021,'K')), NASAPolynomial(coeffs=[39.0627,0.000515528,5.24396e-06,-1.17849e-09,7.64879e-14,-25981.3,-169.453], Tmin=(913.021,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CJO)"""),
)

species(
    label = 'C=COC([O])([CH]C)C(=C)O(9493)',
    structure = SMILES('C=COC([O])([CH]C)C(=C)O'),
    E0 = (-136.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.34626,0.133532,-0.000141761,6.90295e-08,-1.21865e-11,-16087.9,48.1406], Tmin=(100,'K'), Tmax=(1672.33,'K')), NASAPolynomial(coeffs=[36.0791,0.00361666,4.56729e-06,-1.1833e-09,8.56466e-14,-24963,-153.786], Tmin=(1672.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = '[CH2]CC(O)(O[C]=C)C(=C)O(9494)',
    structure = SMILES('[CH2]CC(O)(O[C]=C)C(=C)O'),
    E0 = (-124.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,949.402,1251.29,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158678,'amu*angstrom^2'), symmetry=1, barrier=(3.64832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.31076,0.136445,-0.000147523,7.27521e-08,-1.3046e-11,-14644.1,49.5432], Tmin=(100,'K'), Tmax=(1629.33,'K')), NASAPolynomial(coeffs=[37.8002,0.00211787,4.62988e-06,-1.16043e-09,8.34736e-14,-24259.1,-161.57], Tmin=(1629.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC(O)(CC)C(=C)[O](9495)',
    structure = SMILES('C=[C]OC(O)(CC)C(=C)[O]'),
    E0 = (-192.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.75584,0.129117,-0.000135597,6.61275e-08,-1.18209e-11,-22779.7,47.2631], Tmin=(100,'K'), Tmax=(1620.05,'K')), NASAPolynomial(coeffs=[35.1152,0.00663838,2.34534e-06,-7.35725e-10,5.54852e-14,-31896.3,-148.286], Tmin=(1620.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)O[C]=C(9496)',
    structure = SMILES('[CH]=C(O)C(O)(CC)O[C]=C'),
    E0 = (-82.7285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,426.303,702.113,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.35585,0.138536,-0.000151619,7.55612e-08,-1.36939e-11,-9610.12,48.316], Tmin=(100,'K'), Tmax=(1608.77,'K')), NASAPolynomial(coeffs=[38.4119,0.00139735,4.96631e-06,-1.22818e-09,8.84073e-14,-19384.8,-165.96], Tmin=(1608.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.7285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(O)(CC)C(=C)O(9497)',
    structure = SMILES('[CH]=[C]OC(O)(CC)C(=C)O'),
    E0 = (-82.7285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,426.303,702.113,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155381,'amu*angstrom^2'), symmetry=1, barrier=(3.57252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.35585,0.138536,-0.000151619,7.55612e-08,-1.36939e-11,-9610.12,48.316], Tmin=(100,'K'), Tmax=(1608.77,'K')), NASAPolynomial(coeffs=[38.4119,0.00139735,4.96631e-06,-1.22818e-09,8.84073e-14,-19384.8,-165.96], Tmin=(1608.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.7285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CC([O])(OC=C)C(=C)O(9498)',
    structure = SMILES('[CH2]CC([O])(OC=C)C(=C)O'),
    E0 = (-131.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.09464,0.134914,-0.000145916,7.26646e-08,-1.31962e-11,-15460.3,46.509], Tmin=(100,'K'), Tmax=(1603.51,'K')), NASAPolynomial(coeffs=[36.6267,0.00399446,3.99724e-06,-1.07303e-09,7.91043e-14,-24747.9,-157.415], Tmin=(1603.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'C=COC([O])(CC)C(=C)[O](9499)',
    structure = SMILES('C=COC([O])(CC)C(=C)[O]'),
    E0 = (-198.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.54073,0.127596,-0.000134019,6.60686e-08,-1.19804e-11,-23595.8,44.2327], Tmin=(100,'K'), Tmax=(1590.52,'K')), NASAPolynomial(coeffs=[33.8943,0.0085787,1.68203e-06,-6.42029e-10,5.06473e-14,-32358,-143.852], Tmin=(1590.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)OC=C(9500)',
    structure = SMILES('[CH]=C(O)C([O])(CC)OC=C'),
    E0 = (-89.4264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1021.56,1186.92,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161409,'amu*angstrom^2'), symmetry=1, barrier=(3.71111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.14344,0.137041,-0.000150111,7.55688e-08,-1.38736e-11,-10426.1,45.2957], Tmin=(100,'K'), Tmax=(1582.89,'K')), NASAPolynomial(coeffs=[37.2143,0.00330179,4.32209e-06,-1.13871e-09,8.39025e-14,-19857.6,-161.659], Tmin=(1582.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.4264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)(O)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = 'C=[C]OC1(CC)OC[C]1O(9501)',
    structure = SMILES('C=[C]OC1(CC)OC[C]1O'),
    E0 = (-53.7243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.07949,0.125073,-0.000126737,6.11746e-08,-1.06966e-11,-6124.98,46.2679], Tmin=(100,'K'), Tmax=(1720.88,'K')), NASAPolynomial(coeffs=[29.4763,0.0116426,3.02027e-06,-1.05825e-09,8.2452e-14,-12427.4,-118.591], Tmin=(1720.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.7243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(C=CJO) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C1C[C](O)C([O])(CC)O1(9502)',
    structure = SMILES('C=C1C[C](O)C([O])(CC)O1'),
    E0 = (-185.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.67151,0.0780818,6.90679e-07,-7.11958e-08,3.67694e-11,-22177.1,34.2136], Tmin=(100,'K'), Tmax=(949.281,'K')), NASAPolynomial(coeffs=[27.3208,0.0211397,-5.73638e-06,1.02086e-09,-7.94451e-14,-30240.5,-113.855], Tmin=(949.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CC(C)(O)OJ)"""),
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
    label = 'C=[C]OC(C)([O])C(=C)O(9503)',
    structure = SMILES('C=[C]OC(C)([O])C(=C)O'),
    E0 = (-72.9981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1055.56,1155.43,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163712,'amu*angstrom^2'), symmetry=1, barrier=(3.76407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163712,'amu*angstrom^2'), symmetry=1, barrier=(3.76407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163712,'amu*angstrom^2'), symmetry=1, barrier=(3.76407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163712,'amu*angstrom^2'), symmetry=1, barrier=(3.76407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163712,'amu*angstrom^2'), symmetry=1, barrier=(3.76407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.74268,0.115066,-0.000127997,6.54651e-08,-1.22154e-11,-8508.3,40.7394], Tmin=(100,'K'), Tmax=(1553.63,'K')), NASAPolynomial(coeffs=[31.699,0.00247409,3.80097e-06,-9.98725e-10,7.40151e-14,-16323.6,-131.307], Tmin=(1553.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.9981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)(O)OJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OC(CC)(O1)C(=C)O(9504)',
    structure = SMILES('C=C1OC(CC)(O1)C(=C)O'),
    E0 = (-380.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7378,0.0844596,3.84134e-05,-1.51544e-07,7.57404e-11,-45458.5,30.2945], Tmin=(100,'K'), Tmax=(919.992,'K')), NASAPolynomial(coeffs=[44.2533,-0.0066259,9.40506e-06,-1.88557e-09,1.15812e-13,-58528.4,-212.792], Tmin=(919.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=[C]OC([C]=C)(CC)OO(9505)',
    structure = SMILES('C=[C]OC([C]=C)(CC)OO'),
    E0 = (194.081,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,812.91,1384.28,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155294,'amu*angstrom^2'), symmetry=1, barrier=(3.57053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.11138,0.127594,-0.000133045,6.55403e-08,-1.21114e-11,23623.2,44.6138], Tmin=(100,'K'), Tmax=(1485.07,'K')), NASAPolynomial(coeffs=[34.4333,0.0121004,-1.87916e-06,1.4355e-10,-5.65264e-15,14056.2,-146.031], Tmin=(1485.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.081,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OC([O])(CC)C(C)=O(9506)',
    structure = SMILES('C=[C]OC([O])(CC)C(C)=O'),
    E0 = (-110.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46269,0.0995604,-5.43691e-05,-2.1753e-08,2.14582e-11,-13048,40.674], Tmin=(100,'K'), Tmax=(938.561,'K')), NASAPolynomial(coeffs=[30.1532,0.0162802,-3.5169e-06,5.45886e-10,-4.23573e-14,-21249.3,-121.924], Tmin=(938.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJO)"""),
)

species(
    label = 'H2CC(41)',
    structure = SMILES('[C]=C'),
    E0 = (401.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2480.69,'J/mol'), sigma=(4.48499,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=387.48 K, Pc=62.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(401.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(O)C([O])([O])CC(5967)',
    structure = SMILES('C=C(O)C([O])([O])CC'),
    E0 = (-174.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,301.164,843.095,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.170865,'amu*angstrom^2'), symmetry=1, barrier=(3.92852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170865,'amu*angstrom^2'), symmetry=1, barrier=(3.92852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170865,'amu*angstrom^2'), symmetry=1, barrier=(3.92852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170865,'amu*angstrom^2'), symmetry=1, barrier=(3.92852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.418999,0.0924779,-0.000101782,5.1301e-08,-7.98404e-12,-20828.7,30.8276], Tmin=(100,'K'), Tmax=(882.458,'K')), NASAPolynomial(coeffs=[19.8992,0.0179469,-4.95538e-06,7.10878e-10,-4.28244e-14,-25098.7,-68.531], Tmin=(882.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)OJ) + radical(C=CC(C)(O)OJ)"""),
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
    label = 'C=[C]O[C](CC)C(=C)O(9507)',
    structure = SMILES('[CH2]C(O)=C(CC)O[C]=C'),
    E0 = (53.4535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.153,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.908638,0.098896,-9.82718e-05,5.04177e-08,-1.0187e-11,6613.38,39.257], Tmin=(100,'K'), Tmax=(1209.03,'K')), NASAPolynomial(coeffs=[20.6176,0.0276778,-9.91397e-06,1.69665e-09,-1.126e-13,1408.2,-68.6809], Tmin=(1209.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.4535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=C(O)CJ)"""),
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
    E0 = (-96.7784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-53.2648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-34.0184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (99.7058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-96.7784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (8.77306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (15.6855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (16.0104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-24.3606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-11.9843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-40.8983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (5.02924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-38.4199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-49.6883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-7.60337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-54.9923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-56.3862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-24.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (28.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-53.6832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (346.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-88.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (301.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (107.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (226.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (296.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['CH2CO(28)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['[CH2]C1(O)OC1(CC)O[C]=C(9486)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(43.5136,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['[CH2]C1(O)C(=C)OC1([O])CC(9487)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#COC([O])(CC)C(=C)O(9488)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(148.983,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from 145.1 to 149.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', 'C=[C]OC(=O)C(=C)O(9489)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-DeNd_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2COH(99)', 'C=[C]OC(=O)CC(9490)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(31600,'m^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-NdNd_O;CJ] for rate rule [CO-NdNd_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=COC([O])(CC)C(=C)O(9491)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]OC(O)([CH]C)C(=C)O(9492)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=COC([O])([CH]C)C(=C)O(9493)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.51492e+06,'s^-1'), n=1.58625, Ea=(84.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS_OCs;Y_rad_out;Cs_H_out_H/NonDeC] + [R4H_RSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_SSS_OCs;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC(O)(O[C]=C)C(=C)O(9494)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=[C]OC(O)(CC)C(=C)[O](9495)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(O)C(O)(CC)O[C]=C(9496)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]OC(O)(CC)C(=C)O(9497)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['[CH2]CC([O])(OC=C)C(=C)O(9498)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.13192e+07,'s^-1'), n=1.35333, Ea=(89.175,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS_OCC;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSSS_OCC;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=COC([O])(CC)C(=C)[O](9499)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C([O])(CC)OC=C(9500)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C][O](173)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=[C]OC1(CC)OC[C]1O(9501)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=C1C[C](O)C([O])(CC)O1(9502)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.6e+11,'s^-1'), n=0.27, Ea=(43.0952,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 8 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cddouble
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(S)(23)', 'C=[C]OC(C)([O])C(=C)O(9503)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=C1OC(CC)(O1)C(=C)O(9504)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]OC([C]=C)(CC)OO(9505)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    products = ['C=[C]OC([O])(CC)C(C)=O(9506)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H2CC(41)', 'C=C(O)C([O])([O])CC(5967)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2711.41,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(4)', 'C=[C]O[C](CC)C(=C)O(9507)'],
    products = ['C=[C]OC([O])(CC)C(=C)O(8847)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/ODMustO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '1793',
    isomers = [
        'C=[C]OC([O])(CC)C(=C)O(8847)',
    ],
    reactants = [
        ('CH2CO(28)', 'C=C(O)C(=O)CC(4626)'),
        ('CH2CO(28)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1793',
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

