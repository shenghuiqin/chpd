species(
    label = 'C=C[C]=CC([O])(CC)C(=C)O(29683)',
    structure = SMILES('C=C[C]=CC([O])(CC)C(=C)O'),
    E0 = (101.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,835.362,1371.51,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156093,'amu*angstrom^2'), symmetry=1, barrier=(3.58889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156093,'amu*angstrom^2'), symmetry=1, barrier=(3.58889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156093,'amu*angstrom^2'), symmetry=1, barrier=(3.58889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156093,'amu*angstrom^2'), symmetry=1, barrier=(3.58889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156093,'amu*angstrom^2'), symmetry=1, barrier=(3.58889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156093,'amu*angstrom^2'), symmetry=1, barrier=(3.58889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.2154,0.121361,-0.000111746,5.22779e-08,-9.63879e-12,12405.6,43.8919], Tmin=(100,'K'), Tmax=(1318.55,'K')), NASAPolynomial(coeffs=[26.4657,0.0343519,-1.2762e-05,2.23067e-09,-1.4963e-13,4842.18,-102.409], Tmin=(1318.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(101.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC(C)2OJ)"""),
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
    label = '[CH2]C1(O)OC1(C=[C]C=C)CC(30805)',
    structure = SMILES('[CH2]C1(O)OC1(C=[C]C=C)CC'),
    E0 = (146.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.59885,0.145167,-0.000151808,7.69576e-08,-1.43094e-11,17930.6,47.5649], Tmin=(100,'K'), Tmax=(1581.88,'K')), NASAPolynomial(coeffs=[33.5101,0.0159276,1.91551e-06,-9.65663e-10,8.18968e-14,9987.1,-140.765], Tmin=(1581.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1(O)C(C=C)=CC1([O])CC(29758)',
    structure = SMILES('[CH2]C1(O)C(C=C)=CC1([O])CC'),
    E0 = (170.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68917,0.106934,-6.20846e-05,-7.45733e-09,1.38812e-11,20689.9,41.687], Tmin=(100,'K'), Tmax=(965.876,'K')), NASAPolynomial(coeffs=[26.4887,0.0320415,-1.06927e-05,1.87745e-09,-1.32321e-13,13296.8,-103.37], Tmin=(965.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CC(C)2OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C1[C]=CC(CC)(O1)C(=C)O(30662)',
    structure = SMILES('[CH2]C1[C]=CC(CC)(O1)C(=C)O'),
    E0 = (94.7708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7809,0.105222,-4.99293e-05,-2.28429e-08,1.93095e-11,11626,39.4112], Tmin=(100,'K'), Tmax=(981.292,'K')), NASAPolynomial(coeffs=[28.8859,0.0298508,-1.05874e-05,1.97405e-09,-1.44928e-13,3217.66,-120.136], Tmin=(981.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.7708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(25dihydrofuran) + radical(Cds_S) + radical(CJC(C)OC)"""),
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
    label = 'C=CC#CC([O])(CC)C(=C)O(30806)',
    structure = SMILES('C=CC#CC([O])(CC)C(=C)O'),
    E0 = (72.8136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.182,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47896,0.103822,-6.64855e-05,4.70365e-09,7.45895e-12,8969.37,42.1519], Tmin=(100,'K'), Tmax=(1008.95,'K')), NASAPolynomial(coeffs=[24.6506,0.0339597,-1.27664e-05,2.34219e-09,-1.65882e-13,1979.94,-92.6485], Tmin=(1008.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.8136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C=C=CC([O])(CC)C(=C)O(30807)',
    structure = SMILES('C=C=C=CC([O])(CC)C(=C)O'),
    E0 = (131.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,374.152,762.639,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153945,'amu*angstrom^2'), symmetry=1, barrier=(3.5395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153945,'amu*angstrom^2'), symmetry=1, barrier=(3.5395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153945,'amu*angstrom^2'), symmetry=1, barrier=(3.5395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153945,'amu*angstrom^2'), symmetry=1, barrier=(3.5395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153945,'amu*angstrom^2'), symmetry=1, barrier=(3.5395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.182,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.61467,0.115648,-0.000106426,4.97604e-08,-9.28851e-12,16032.6,41.2175], Tmin=(100,'K'), Tmax=(1286.85,'K')), NASAPolynomial(coeffs=[23.2792,0.0382694,-1.62312e-05,3.03469e-09,-2.11058e-13,9625.64,-85.1592], Tmin=(1286.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C[C]=CC(=O)C(=C)O(30808)',
    structure = SMILES('C=C[C]=CC(=O)C(=C)O'),
    E0 = (21.0548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,375,552.5,462.5,1710,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.950941,'amu*angstrom^2'), symmetry=1, barrier=(21.864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951506,'amu*angstrom^2'), symmetry=1, barrier=(21.877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951277,'amu*angstrom^2'), symmetry=1, barrier=(21.8717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.951843,'amu*angstrom^2'), symmetry=1, barrier=(21.8848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.129,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.360548,0.0926055,-0.000100757,5.61544e-08,-1.2345e-11,2692.2,29.8418], Tmin=(100,'K'), Tmax=(1110.71,'K')), NASAPolynomial(coeffs=[18.2312,0.0256511,-1.03363e-05,1.88209e-09,-1.29334e-13,-1437.79,-61.8045], Tmin=(1110.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.0548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
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
    label = 'C=C[C]=CC(=O)CC(27172)',
    structure = SMILES('C=C[C]=CC(=O)CC'),
    E0 = (98.8568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,225.125,225.128,225.141],'cm^-1')),
        HinderedRotor(inertia=(0.261381,'amu*angstrom^2'), symmetry=1, barrier=(9.40576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261494,'amu*angstrom^2'), symmetry=1, barrier=(9.40617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261383,'amu*angstrom^2'), symmetry=1, barrier=(9.4053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957381,'amu*angstrom^2'), symmetry=1, barrier=(34.4381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660557,0.0792073,-8.33263e-05,5.56152e-08,-1.63174e-11,12004.8,26.9776], Tmin=(100,'K'), Tmax=(803.67,'K')), NASAPolynomial(coeffs=[7.35607,0.0458919,-2.11625e-05,4.06295e-09,-2.85329e-13,10928.3,-3.86284], Tmin=(803.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.8568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
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
    label = 'C=CC=[C]C([O])(CC)C(=C)O(30809)',
    structure = SMILES('C=CC=[C]C([O])(CC)C(=C)O'),
    E0 = (140.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,971.12,1240.32,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160067,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160067,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160067,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160067,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160067,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160067,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20138,0.1197,-0.000108043,4.89055e-08,-8.7046e-12,17078.3,44.724], Tmin=(100,'K'), Tmax=(1364.34,'K')), NASAPolynomial(coeffs=[27.442,0.032793,-1.2497e-05,2.21968e-09,-1.50183e-13,8989.38,-107.499], Tmin=(1364.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CC([O])(CC)C(=C)O(30810)',
    structure = SMILES('C=C=C[CH]C([O])(CC)C(=C)O'),
    E0 = (43.5591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78449,0.106569,-5.44444e-05,-1.68686e-08,1.69778e-11,5465.67,44.7717], Tmin=(100,'K'), Tmax=(980.743,'K')), NASAPolynomial(coeffs=[27.8788,0.0319896,-1.13501e-05,2.08129e-09,-1.50314e-13,-2584.43,-109.137], Tmin=(980.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.5591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC(O)C=C) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C[C]=CC(O)([CH]C)C(=C)O(30811)',
    structure = SMILES('C=C[C]=CC(O)([CH]C)C(=C)O'),
    E0 = (71.9779,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180,180,180,303.006,824.255,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152767,'amu*angstrom^2'), symmetry=1, barrier=(3.51242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.12561,0.119364,-9.65262e-05,2.61101e-08,2.54065e-12,8891.24,44.8415], Tmin=(100,'K'), Tmax=(974.042,'K')), NASAPolynomial(coeffs=[27.546,0.0304396,-1.02887e-05,1.7895e-09,-1.24184e-13,1549.08,-105.544], Tmin=(974.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.9779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CCJCO)"""),
)

species(
    label = 'C=C[C]=[C]C(O)(CC)C(=C)O(30812)',
    structure = SMILES('[CH2][CH]C#CC(O)(CC)C(=C)O'),
    E0 = (74.0418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.52541,0.124817,-0.000124196,6.36492e-08,-1.26208e-11,9156.05,48.0427], Tmin=(100,'K'), Tmax=(1332.12,'K')), NASAPolynomial(coeffs=[26.7274,0.0292512,-7.88479e-06,1.08618e-09,-6.22492e-14,2048.03,-98.9005], Tmin=(1332.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.0418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC=CC([O])(CC)C(=C)O(30813)',
    structure = SMILES('[CH]=CC=CC([O])(CC)C(=C)O'),
    E0 = (149.279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,418.581,715.153,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00856,0.116069,-9.19641e-05,2.81749e-08,-5.27663e-13,18184.4,43.9489], Tmin=(100,'K'), Tmax=(1049.85,'K')), NASAPolynomial(coeffs=[26.4593,0.0343106,-1.33058e-05,2.45598e-09,-1.7311e-13,10735.2,-101.786], Tmin=(1049.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=CC([O])([CH]C)C(=C)O(30814)',
    structure = SMILES('C=CC=CC([O])([CH]C)C(=C)O'),
    E0 = (102.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69634,0.106846,-6.51533e-05,-3.73911e-10,9.7943e-12,12499.3,44.9199], Tmin=(100,'K'), Tmax=(1007.94,'K')), NASAPolynomial(coeffs=[26.3558,0.0335556,-1.26848e-05,2.35956e-09,-1.69165e-13,4912.31,-100.221], Tmin=(1007.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)(C=[C]C=C)C(=C)O(30815)',
    structure = SMILES('[CH2]CC(O)(C=[C]C=C)C(=C)O'),
    E0 = (77.3222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,262.614,867.304,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15201,'amu*angstrom^2'), symmetry=1, barrier=(3.495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.65618,0.129601,-0.00013002,6.5607e-08,-1.2871e-11,9553.6,46.0419], Tmin=(100,'K'), Tmax=(1253.58,'K')), NASAPolynomial(coeffs=[28.8727,0.0289955,-9.63787e-06,1.58593e-09,-1.03282e-13,1648.86,-113.192], Tmin=(1253.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.3222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = 'C=C[C]=CC(O)(CC)C(=C)[O](30816)',
    structure = SMILES('C=C[C]=CC(O)(CC)C(=C)[O]'),
    E0 = (9.88064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09193,0.122166,-0.000117757,5.86268e-08,-1.1533e-11,1417.67,43.7285], Tmin=(100,'K'), Tmax=(1237.9,'K')), NASAPolynomial(coeffs=[24.6341,0.0358063,-1.31123e-05,2.27071e-09,-1.51605e-13,-5199.15,-90.913], Tmin=(1237.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.88064,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(O)C(O)(C=[C]C=C)CC(30817)',
    structure = SMILES('[CH]=C(O)C(O)(C=[C]C=C)CC'),
    E0 = (119.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,204.049,1532.77,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148857,'amu*angstrom^2'), symmetry=1, barrier=(3.42251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.73303,0.132013,-0.000135032,6.9361e-08,-1.38361e-11,14589.1,44.9327], Tmin=(100,'K'), Tmax=(1233.99,'K')), NASAPolynomial(coeffs=[29.1136,0.0287837,-9.55145e-06,1.57054e-09,-1.02309e-13,6729.33,-115.405], Tmin=(1233.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C=CC=C)C(=C)O(30818)',
    structure = SMILES('[CH2]CC([O])(C=CC=C)C(=C)O'),
    E0 = (107.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,934.901,1270.57,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158878,'amu*angstrom^2'), symmetry=1, barrier=(3.65291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158878,'amu*angstrom^2'), symmetry=1, barrier=(3.65291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158878,'amu*angstrom^2'), symmetry=1, barrier=(3.65291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158878,'amu*angstrom^2'), symmetry=1, barrier=(3.65291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158878,'amu*angstrom^2'), symmetry=1, barrier=(3.65291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158878,'amu*angstrom^2'), symmetry=1, barrier=(3.65291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88872,0.113174,-8.53814e-05,2.25488e-08,1.16606e-12,13146.9,44.9023], Tmin=(100,'K'), Tmax=(1047.26,'K')), NASAPolynomial(coeffs=[25.964,0.0349314,-1.36185e-05,2.5232e-09,-1.78282e-13,5769.93,-98.1249], Tmin=(1047.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=CC=CC([O])(CC)C(=C)[O](30819)',
    structure = SMILES('C=CC=CC([O])(CC)C(=C)[O]'),
    E0 = (39.9878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.67277,0.109687,-8.61448e-05,3.13078e-08,-3.69862e-12,5026.32,43.8487], Tmin=(100,'K'), Tmax=(1157.53,'K')), NASAPolynomial(coeffs=[23.6595,0.0385964,-1.5337e-05,2.80313e-09,-1.93641e-13,-1940.24,-86.8309], Tmin=(1157.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.9878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(C=CC=C)CC(30820)',
    structure = SMILES('[CH]=C(O)C([O])(C=CC=C)CC'),
    E0 = (149.279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,418.581,715.153,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155273,'amu*angstrom^2'), symmetry=1, barrier=(3.57002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00856,0.116069,-9.19641e-05,2.81749e-08,-5.27663e-13,18184.4,43.9489], Tmin=(100,'K'), Tmax=(1049.85,'K')), NASAPolynomial(coeffs=[26.4593,0.0343106,-1.33058e-05,2.45598e-09,-1.7311e-13,10735.2,-101.786], Tmin=(1049.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=[C][C]=CC(O)(CC)C(=C)O(30821)',
    structure = SMILES('C=C=[C][CH]C(O)(CC)C(=C)O'),
    E0 = (52.2983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17227,0.117165,-8.14847e-05,5.8042e-09,1.06964e-11,6529.02,45.4225], Tmin=(100,'K'), Tmax=(963.49,'K')), NASAPolynomial(coeffs=[29.5409,0.0280883,-9.10211e-06,1.59227e-09,-1.13183e-13,-1558.58,-116.653], Tmin=(963.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[CH]=C[C]=CC(O)(CC)C(=C)O(30822)',
    structure = SMILES('[CH]C=C=CC(O)(CC)C(=C)O'),
    E0 = (96.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.37355,0.126143,-0.000113627,5.23097e-08,-9.56798e-12,11853,44.7867], Tmin=(100,'K'), Tmax=(1320.61,'K')), NASAPolynomial(coeffs=[26.1651,0.0397035,-1.54467e-05,2.74701e-09,-1.8555e-13,4315.22,-100.832], Tmin=(1320.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(O)C([O])(C=C1[CH]C1)CC(30823)',
    structure = SMILES('C=C(O)C([O])([CH]C1=CC1)CC'),
    E0 = (132.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66676,0.0995916,-2.69896e-05,-5.10533e-08,3.05985e-11,16157.1,43.8348], Tmin=(100,'K'), Tmax=(959.929,'K')), NASAPolynomial(coeffs=[30.0128,0.0277562,-8.765e-06,1.59076e-09,-1.18672e-13,7302.75,-122.145], Tmin=(959.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.453,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=CCJC(O)C=C) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C[C]=CC1(CC)OC[C]1O(30824)',
    structure = SMILES('C=C[C]=CC1(CC)OC[C]1O'),
    E0 = (153.248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.42517,0.122425,-0.00012052,6.30198e-08,-1.27163e-11,18678.9,42.8395], Tmin=(100,'K'), Tmax=(1358.76,'K')), NASAPolynomial(coeffs=[23.8671,0.0325794,-7.59678e-06,8.74e-10,-4.17412e-14,12682.7,-87.8386], Tmin=(1358.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Oxetane) + radical(C=CJC=C) + radical(C2CsJOH)"""),
)

species(
    label = 'C=CC1=CC([O])(CC)[C](O)C1(29896)',
    structure = SMILES('C=CC1=CC([O])(CC)[C](O)C1'),
    E0 = (60.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36957,0.101408,-5.792e-05,-1.78975e-09,9.08849e-12,7431.14,39.7111], Tmin=(100,'K'), Tmax=(1012.03,'K')), NASAPolynomial(coeffs=[23.2485,0.0380425,-1.4302e-05,2.61161e-09,-1.83842e-13,710.371,-87.9376], Tmin=(1012.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentene) + radical(C2CsJOH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2][C](O)C1(C=C=CCO1)CC(30712)',
    structure = SMILES('[CH2][C](O)C1(C=C=CCO1)CC'),
    E0 = (151.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.80823,0.184401,-0.000281672,2.30882e-07,-7.47464e-11,18515.9,30.8173], Tmin=(100,'K'), Tmax=(835.404,'K')), NASAPolynomial(coeffs=[19.1038,0.0587464,-2.74156e-05,5.12689e-09,-3.48527e-13,15244.3,-72.2686], Tmin=(835.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C=C=CC(O)(CC)C(=C)O(30825)',
    structure = SMILES('C=C=C=CC(O)(CC)C(=C)O'),
    E0 = (-97.5342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.27039,0.125158,-0.000120838,5.88147e-08,-1.1255e-11,-11494.2,42.3585], Tmin=(100,'K'), Tmax=(1272.42,'K')), NASAPolynomial(coeffs=[26.8122,0.033733,-1.30614e-05,2.34639e-09,-1.60255e-13,-18895.2,-104.954], Tmin=(1272.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.5342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = 'C=C[C]=CC(C)([O])C(=C)O(30826)',
    structure = SMILES('C=C[C]=CC(C)([O])C(=C)O'),
    E0 = (124.959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180,180,180,247.318,882.939,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149974,'amu*angstrom^2'), symmetry=1, barrier=(3.44819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149974,'amu*angstrom^2'), symmetry=1, barrier=(3.44819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149974,'amu*angstrom^2'), symmetry=1, barrier=(3.44819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149974,'amu*angstrom^2'), symmetry=1, barrier=(3.44819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149974,'amu*angstrom^2'), symmetry=1, barrier=(3.44819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27528,0.10325,-8.69385e-05,3.10426e-08,-2.35307e-12,15230.1,38.2641], Tmin=(100,'K'), Tmax=(1032.96,'K')), NASAPolynomial(coeffs=[23.1755,0.0297977,-1.11046e-05,1.9964e-09,-1.38602e-13,9046.12,-85.972], Tmin=(1032.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CC1=CC(CC)(O1)C(=C)O(30505)',
    structure = SMILES('C=CC1=CC(CC)(O1)C(=C)O'),
    E0 = (-180.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.85436,0.0986119,-1.11189e-05,-7.79303e-08,4.27194e-11,-21433.3,35.5913], Tmin=(100,'K'), Tmax=(947.108,'K')), NASAPolynomial(coeffs=[34.5335,0.0197024,-4.56372e-06,8.10497e-10,-6.76884e-14,-31679.5,-155.686], Tmin=(947.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=[C]C(C=[C]C=C)(CC)OO(30827)',
    structure = SMILES('C=[C]C(C=[C]C=C)(CC)OO'),
    E0 = (395.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55444,0.121688,-0.000118068,6.09571e-08,-1.28254e-11,47826,41.4803], Tmin=(100,'K'), Tmax=(1135.73,'K')), NASAPolynomial(coeffs=[19.3251,0.0481498,-2.09415e-05,3.94364e-09,-2.75184e-13,43083.4,-61.9085], Tmin=(1135.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C]=CC([O])(CC)C(C)=O(30828)',
    structure = SMILES('C=C[C]=CC([O])(CC)C(C)=O'),
    E0 = (112.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35537,0.128054,-0.000149647,8.95715e-08,-1.61845e-11,13672.9,39.8284], Tmin=(100,'K'), Tmax=(646.112,'K')), NASAPolynomial(coeffs=[13.2668,0.0558807,-2.46947e-05,4.60222e-09,-3.16125e-13,11400.3,-27.2926], Tmin=(646.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=CJC=C)"""),
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
    label = 'C=C[C]=C[C](CC)C(=C)O(30829)',
    structure = SMILES('[CH2]C=[C]C=C(CC)C(=C)O'),
    E0 = (152.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (136.191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.51047,0.119883,-0.000114479,5.53911e-08,-1.03036e-11,18606.6,39.7253], Tmin=(100,'K'), Tmax=(1446.35,'K')), NASAPolynomial(coeffs=[27.7319,0.0261371,-6.77334e-06,9.14396e-10,-5.21775e-14,10915.6,-113.683], Tmin=(1446.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C1OC1(CC)C(=C)O(30709)',
    structure = SMILES('[CH2]C=[C]C1OC1(CC)C(=C)O'),
    E0 = (133.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.72666,0.129536,-0.000124057,5.90623e-08,-1.04787e-11,16385.5,48.37], Tmin=(100,'K'), Tmax=(1633.7,'K')), NASAPolynomial(coeffs=[29.823,0.022154,-2.29121e-06,-8.27128e-11,1.91105e-14,8791.35,-119.648], Tmin=(1633.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1(O)CC=C=CC1([O])CC(29941)',
    structure = SMILES('[CH2]C1(O)CC=C=CC1([O])CC'),
    E0 = (150.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.64751,0.178062,-0.000253016,1.92995e-07,-5.87755e-11,18345,29.2576], Tmin=(100,'K'), Tmax=(804.835,'K')), NASAPolynomial(coeffs=[20.8621,0.0562362,-2.59377e-05,4.87728e-09,-3.34997e-13,14400.2,-83.6625], Tmin=(804.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(C=CC(C)2OJ) + radical(CJC(C)2O)"""),
)

species(
    label = 'C=C(O)C([O])(C=C=[C]C)CC(30830)',
    structure = SMILES('C=C(O)C([O])([CH]C#CC)CC'),
    E0 = (156.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.22602,0.118324,-0.000111117,5.43603e-08,-1.03898e-11,19082,48.4635], Tmin=(100,'K'), Tmax=(1357.43,'K')), NASAPolynomial(coeffs=[25.3028,0.0317058,-9.32654e-06,1.38487e-09,-8.36827e-14,12114.8,-90.8942], Tmin=(1357.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C([O])([C]=C=CC)CC(30831)',
    structure = SMILES('C=C(O)C([O])(C#C[CH]C)CC'),
    E0 = (97.8981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.56772,0.110588,-8.55129e-05,2.46942e-08,1.15774e-12,11985.5,43.3387], Tmin=(100,'K'), Tmax=(965.407,'K')), NASAPolynomial(coeffs=[22.0563,0.0383517,-1.31243e-05,2.22302e-09,-1.4894e-13,6229.02,-75.9917], Tmin=(965.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.8981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C=C=CC(30832)',
    structure = SMILES('C=C(O)C([O])([CH]C)C=C=CC'),
    E0 = (154.866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,180,180,180,416.413,720.084,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155297,'amu*angstrom^2'), symmetry=1, barrier=(3.57059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155297,'amu*angstrom^2'), symmetry=1, barrier=(3.57059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155297,'amu*angstrom^2'), symmetry=1, barrier=(3.57059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155297,'amu*angstrom^2'), symmetry=1, barrier=(3.57059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155297,'amu*angstrom^2'), symmetry=1, barrier=(3.57059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155297,'amu*angstrom^2'), symmetry=1, barrier=(3.57059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68396,0.113346,-9.72649e-05,4.22912e-08,-7.34095e-12,18840.3,44.9149], Tmin=(100,'K'), Tmax=(1379.06,'K')), NASAPolynomial(coeffs=[24.15,0.0384145,-1.57622e-05,2.89125e-09,-1.98435e-13,11715,-88.0223], Tmin=(1379.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CC([O])(C=C=CC)C(=C)O(30833)',
    structure = SMILES('[CH2]CC([O])(C=C=CC)C(=C)O'),
    E0 = (160.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,396.298,741.368,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154853,'amu*angstrom^2'), symmetry=1, barrier=(3.56038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154853,'amu*angstrom^2'), symmetry=1, barrier=(3.56038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154853,'amu*angstrom^2'), symmetry=1, barrier=(3.56038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154853,'amu*angstrom^2'), symmetry=1, barrier=(3.56038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154853,'amu*angstrom^2'), symmetry=1, barrier=(3.56038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154853,'amu*angstrom^2'), symmetry=1, barrier=(3.56038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48971,0.115269,-0.000102823,4.72614e-08,-8.77392e-12,19470.9,43.5002], Tmin=(100,'K'), Tmax=(1282.52,'K')), NASAPolynomial(coeffs=[21.5901,0.0432867,-1.8636e-05,3.50023e-09,-2.43695e-13,13550.8,-73.5896], Tmin=(1282.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(C=C=CC)CC(30834)',
    structure = SMILES('C=C([O])C([O])(C=C=CC)CC'),
    E0 = (92.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1067.92,1152.96,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.16448,'amu*angstrom^2'), symmetry=1, barrier=(3.78172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16448,'amu*angstrom^2'), symmetry=1, barrier=(3.78172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16448,'amu*angstrom^2'), symmetry=1, barrier=(3.78172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16448,'amu*angstrom^2'), symmetry=1, barrier=(3.78172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16448,'amu*angstrom^2'), symmetry=1, barrier=(3.78172,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.925715,0.107833,-9.05414e-05,4.02466e-08,-7.41983e-12,11335,41.1881], Tmin=(100,'K'), Tmax=(1264.59,'K')), NASAPolynomial(coeffs=[17.3013,0.0501797,-2.21565e-05,4.1957e-09,-2.92893e-13,6725.06,-51.0254], Tmin=(1264.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(C=C=CC)CC(30835)',
    structure = SMILES('[CH]=C(O)C([O])(C=C=CC)CC'),
    E0 = (202.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1600,1605.09,2891.92,3200],'cm^-1')),
        HinderedRotor(inertia=(0.150377,'amu*angstrom^2'), symmetry=1, barrier=(3.45746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150377,'amu*angstrom^2'), symmetry=1, barrier=(3.45746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150377,'amu*angstrom^2'), symmetry=1, barrier=(3.45746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150377,'amu*angstrom^2'), symmetry=1, barrier=(3.45746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150377,'amu*angstrom^2'), symmetry=1, barrier=(3.45746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150377,'amu*angstrom^2'), symmetry=1, barrier=(3.45746,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57326,0.117739,-0.00010796,5.11021e-08,-9.75608e-12,24506.8,42.4168], Tmin=(100,'K'), Tmax=(1250.93,'K')), NASAPolynomial(coeffs=[21.7163,0.043268,-1.86602e-05,3.51089e-09,-2.44878e-13,18680.1,-75.1558], Tmin=(1250.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](O)C1(C=C(C=C)O1)CC(30475)',
    structure = SMILES('[CH2][C](O)C1(C=C(C=C)O1)CC'),
    E0 = (128.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.34997,0.121874,-9.27189e-05,1.48104e-08,8.34324e-12,15704.6,40.02], Tmin=(100,'K'), Tmax=(954.559,'K')), NASAPolynomial(coeffs=[30.3574,0.0264217,-8.1037e-06,1.37583e-09,-9.68609e-14,7564.9,-126.182], Tmin=(954.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C1[C]=CC1(30836)',
    structure = SMILES('C=C(O)C([O])(CC)C1[C]=CC1'),
    E0 = (214.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38442,0.0987682,-4.31571e-05,-2.3192e-08,1.82502e-11,26017,42.6795], Tmin=(100,'K'), Tmax=(978.297,'K')), NASAPolynomial(coeffs=[25.4609,0.0336771,-1.18496e-05,2.14991e-09,-1.53844e-13,18626.8,-97.1705], Tmin=(978.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CC(C)2OJ) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'CCC1([O])C=C=CCC[C]1O(30457)',
    structure = SMILES('CCC1([O])C=C=CCC[C]1O'),
    E0 = (184.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.762489,0.100076,-7.34208e-05,2.79755e-08,-4.41271e-12,22337.3,35.1644], Tmin=(100,'K'), Tmax=(1455.15,'K')), NASAPolynomial(coeffs=[17.8569,0.0488948,-2.06623e-05,3.80481e-09,-2.60135e-13,16918.5,-61.6477], Tmin=(1455.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C2CsJOH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C#CC([CH2])C[C](O)C(=O)CC(29690)',
    structure = SMILES('C#CC([CH2])CC(O)=C([O])CC'),
    E0 = (51.9968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,750,770,3400,2100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4871.92,'J/mol'), sigma=(8.06899,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=760.98 K, Pc=21.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.67145,0.129546,-0.000135374,7.23762e-08,-1.49057e-11,6508.66,46.9057], Tmin=(100,'K'), Tmax=(1284.85,'K')), NASAPolynomial(coeffs=[27.4131,0.0277953,-7.13789e-06,9.36913e-10,-5.15917e-14,-554.259,-103.176], Tmin=(1284.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.9968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
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
    label = '[CH]=C=CC([O])(CC)C(=C)O(30837)',
    structure = SMILES('[CH]=C=CC([O])(CC)C(=C)O'),
    E0 = (145.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,238.08,1495.9,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.144553,'amu*angstrom^2'), symmetry=1, barrier=(3.32356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144553,'amu*angstrom^2'), symmetry=1, barrier=(3.32356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144553,'amu*angstrom^2'), symmetry=1, barrier=(3.32356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144553,'amu*angstrom^2'), symmetry=1, barrier=(3.32356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144553,'amu*angstrom^2'), symmetry=1, barrier=(3.32356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34195,0.105501,-9.89517e-05,4.72164e-08,-8.88326e-12,17698,40.0072], Tmin=(100,'K'), Tmax=(1292.48,'K')), NASAPolynomial(coeffs=[23.0745,0.029936,-1.12536e-05,1.98111e-09,-1.33524e-13,11386.4,-84.0523], Tmin=(1292.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(C)2OJ) + radical(C=C=CJ)"""),
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
    E0 = (101.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (147.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (170.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (200.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (296.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (360.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (156.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (260.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (125.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (101.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (336.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (315.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (177.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (176.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (371.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (145.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (161.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (202.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (163.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (177.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (134.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (182.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (130.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (229.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (266.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (238.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (226.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (140.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (162.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (126.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (544.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (109.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (503.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (305.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (395.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (147.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (208.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (371.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (339.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (215.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (223.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (199.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (334.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (227.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (219.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (184.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (177.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (527.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2]C1(O)OC1(C=[C]C=C)CC(30805)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(46.2682,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2]C1(O)C(C=C)=CC1([O])CC(29758)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(69.0081,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_2H;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_2H_secNd;radadd_intra_cdsingleDe]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 67.6 to 69.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2]C1[C]=CC(CC)(O1)C(=C)O(30662)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.82166e+09,'s^-1'), n=0.527281, Ea=(99.1325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_O] + [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC#CC([O])(CC)C(=C)O(30806)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C=C=CC([O])(CC)C(=C)O(30807)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5(29)', 'C=C[C]=CC(=O)C(=C)O(30808)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-DeDe_O;CsJ-CsHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', 'C=C[C]=CC(=O)CC(27172)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CdCs_O;CdsJ-O2s]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]C=C(4699)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(11.9336,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 9.7 to 11.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=CC=[C]C([O])(CC)C(=C)O(30809)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=[C]C=CC([O])(CC)C(=C)O(30810)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=CC(O)([CH]C)C(=C)O(30811)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C[C]=[C]C(O)(CC)C(=C)O(30812)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC=CC([O])(CC)C(=C)O(30813)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=CC=CC([O])([CH]C)C(=C)O(30814)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC(O)(C=[C]C=C)C(=C)O(30815)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C[C]=CC(O)(CC)C(=C)[O](30816)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(O)C(O)(C=[C]C=C)CC(30817)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC([O])(C=CC=C)C(=C)O(30818)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.79808e+06,'s^-1'), n=1.3209, Ea=(69.9539,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_SSSR;C_rad_out_2H;XH_out] + [R5H_SSSD;Y_rad_out;Cd_H_out_single] for rate rule [R5H_SSSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=CC=CC([O])(CC)C(=C)[O](30819)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(O)C([O])(C=CC=C)CC(30820)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 3.60555127546
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=[C][C]=CC(O)(CC)C(=C)O(30821)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.658e+09,'s^-1'), n=0.699, Ea=(29.5516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_doubleC] for rate rule [R5HJ_3;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH]=C[C]=CC(O)(CC)C(=C)O(30822)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C(O)C([O])(C=C1[CH]C1)CC(30823)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C[C]=CC1(CC)OC[C]1O(30824)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=CC1=CC([O])(CC)[C](O)C1(29896)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.6e+10,'s^-1'), n=0.2, Ea=(38.9112,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS_D;doublebond_intra_secNd_2H;radadd_intra_cdsingle] for rate rule [R5_DS_D;doublebond_intra_secNd_2H;radadd_intra_cdsingleDe]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2][C](O)C1(C=C=CCO1)CC(30712)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C=C=CC(O)(CC)C(=C)O(30825)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(S)(23)', 'C=C[C]=CC(C)([O])C(=C)O(30826)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=CC1=CC(CC)(O1)C(=C)O(30505)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(C=[C]C=C)(CC)OO(30827)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C[C]=CC([O])(CC)C(C)=O(30828)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(4)', 'C=C[C]=C[C](CC)C(=C)O(30829)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2]C=[C]C1OC1(CC)C(=C)O(30709)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2]C1(O)CC=C=CC1([O])CC(29941)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.95984e+07,'s^-1'), n=0.765, Ea=(106.901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra_2H_secNd;radadd_intra_cs2H] + [R7_SMMS;doublebond_intra_2H;radadd_intra_cs2H] for rate rule [R7_SMMS;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C(O)C([O])(C=C=[C]C)CC(30830)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C(O)C([O])([C]=C=CC)CC(30831)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C(O)C([O])([CH]C)C=C=CC(30832)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(138.3,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R6H;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]CC([O])(C=C=CC)C(=C)O(30833)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 115 used for R7H;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R7H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])C([O])(C=C=CC)CC(30834)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5.85e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 143 used for R7H;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R7H;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C(O)C([O])(C=C=CC)CC(30835)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R7H;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['[CH2][C](O)C1(C=C(C=C)O1)CC(30475)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C=C(O)C([O])(CC)C1[C]=CC1(30836)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['CCC1([O])C=C=CCC[C]1O(30457)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(185000,'s^-1'), n=1.07, Ea=(83.0915,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 50 used for R7_linear;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R7_linear;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 80.8 to 83.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    products = ['C#CC([CH2])C[C](O)C(=O)CC(29690)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.213e+11,'s^-1'), n=0.07, Ea=(76.6885,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 3 used for 1_4_5_hexatriene
Exact match found for rate rule [1_4_5_hexatriene]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CH2(19)', '[CH]=C=CC([O])(CC)C(=C)O(30837)'],
    products = ['C=C[C]=CC([O])(CC)C(=C)O(29683)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4702',
    isomers = [
        'C=C[C]=CC([O])(CC)C(=C)O(29683)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'CH2CHCCH(26391)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4702',
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

