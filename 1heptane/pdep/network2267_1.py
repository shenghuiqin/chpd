species(
    label = 'C=C(O)C([O])(CC)C([O])C=O(12330)',
    structure = SMILES('C=C(O)C([O])(CC)C([O])C=O'),
    E0 = (-280.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,180,1093.77,1121.1,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163623,'amu*angstrom^2'), symmetry=1, barrier=(3.76202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163623,'amu*angstrom^2'), symmetry=1, barrier=(3.76202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163623,'amu*angstrom^2'), symmetry=1, barrier=(3.76202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163623,'amu*angstrom^2'), symmetry=1, barrier=(3.76202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163623,'amu*angstrom^2'), symmetry=1, barrier=(3.76202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163623,'amu*angstrom^2'), symmetry=1, barrier=(3.76202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.83581,0.113203,-0.000108403,5.21002e-08,-9.78493e-12,-33466.9,48.2139], Tmin=(100,'K'), Tmax=(1302.4,'K')), NASAPolynomial(coeffs=[25.9851,0.0277579,-9.99262e-06,1.72635e-09,-1.15445e-13,-40713.7,-93.3563], Tmin=(1302.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-280.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=OCOJ)"""),
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
    label = '[CH2]C1(O)OC1(CC)C([O])C=O(12672)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C([O])C=O'),
    E0 = (-233.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.74283,0.133864,-0.000142013,7.24358e-08,-1.36148e-11,-27822.9,47.9843], Tmin=(100,'K'), Tmax=(1542.23,'K')), NASAPolynomial(coeffs=[32.5039,0.0131861,1.2967e-06,-7.25184e-10,6.22077e-14,-35831.7,-132.307], Tmin=(1542.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(C=OCOJ) + radical(CJC(O)2C)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C1OC1[O](12673)',
    structure = SMILES('C=C(O)C([O])(CC)C1OC1[O]'),
    E0 = (-228.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42625,0.110614,-9.74902e-05,3.54949e-08,-1.12609e-12,-27290.4,41.8377], Tmin=(100,'K'), Tmax=(908.356,'K')), NASAPolynomial(coeffs=[21.8815,0.0328695,-1.02139e-05,1.60926e-09,-1.03139e-13,-32551.8,-74.0212], Tmin=(908.356,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CC(C)2OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1(O)OC(C=O)C1([O])CC(12446)',
    structure = SMILES('[CH2]C1(O)OC(C=O)C1([O])CC'),
    E0 = (-260.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-4.34394,0.13198,-0.000132897,6.37199e-08,-1.1137e-11,-30941.6,48.6237], Tmin=(100,'K'), Tmax=(1701.42,'K')), NASAPolynomial(coeffs=[32.0528,0.0126114,2.13995e-06,-8.68e-10,6.90275e-14,-38434.4,-131.936], Tmin=(1701.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CC(C)2OJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=C(O)C1(CC)OC([O])C1[O](12660)',
    structure = SMILES('C=C(O)C1(CC)OC([O])C1[O]'),
    E0 = (-232.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.47665,0.114312,-0.000107671,5.15194e-08,-9.35501e-12,-27750.7,43.8887], Tmin=(100,'K'), Tmax=(1547.09,'K')), NASAPolynomial(coeffs=[26.1608,0.0242062,-4.73358e-06,4.50971e-10,-1.81988e-14,-34689.2,-100.554], Tmin=(1547.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.871,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCOJ) + radical(CC(C)OJ)"""),
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
    label = 'C=C(O)C([O])(CC)C(=O)C=O(12674)',
    structure = SMILES('C=C(O)C([O])(CC)C(=O)C=O'),
    E0 = (-406.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,540.448,609.212,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15868,'amu*angstrom^2'), symmetry=1, barrier=(3.64837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15868,'amu*angstrom^2'), symmetry=1, barrier=(3.64837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15868,'amu*angstrom^2'), symmetry=1, barrier=(3.64837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15868,'amu*angstrom^2'), symmetry=1, barrier=(3.64837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15868,'amu*angstrom^2'), symmetry=1, barrier=(3.64837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15868,'amu*angstrom^2'), symmetry=1, barrier=(3.64837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.144,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.6445,0.131218,-0.000177461,1.28882e-07,-3.74861e-11,-48668.4,39.4871], Tmin=(100,'K'), Tmax=(840.444,'K')), NASAPolynomial(coeffs=[16.9244,0.0428364,-1.97112e-05,3.7428e-09,-2.59822e-13,-51789.5,-46.8683], Tmin=(840.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(C=CC(C)(C=O)OJ)"""),
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
    label = 'C=C(O)C(=O)C([O])C=O(12675)',
    structure = SMILES('C=C(O)C(=O)C([O])C=O'),
    E0 = (-366.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,742.458,742.46,4000],'cm^-1')),
        HinderedRotor(inertia=(0.720222,'amu*angstrom^2'), symmetry=1, barrier=(16.5593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720258,'amu*angstrom^2'), symmetry=1, barrier=(16.5602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409438,'amu*angstrom^2'), symmetry=1, barrier=(16.5592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.720254,'amu*angstrom^2'), symmetry=1, barrier=(16.5601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396368,0.0807286,-9.10626e-05,5.28766e-08,-1.22712e-11,-43942,31.2777], Tmin=(100,'K'), Tmax=(1044.37,'K')), NASAPolynomial(coeffs=[14.8007,0.025559,-1.18238e-05,2.29482e-09,-1.62954e-13,-46950.6,-38.8402], Tmin=(1044.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
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
    label = 'CCC(=O)C([O])C=O(11465)',
    structure = SMILES('CCC(=O)C([O])C=O'),
    E0 = (-291.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,375,552.5,462.5,1710,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,284.801,284.803,2650.74],'cm^-1')),
        HinderedRotor(inertia=(0.00207834,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102611,'amu*angstrom^2'), symmetry=1, barrier=(5.90595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102607,'amu*angstrom^2'), symmetry=1, barrier=(5.90596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102607,'amu*angstrom^2'), symmetry=1, barrier=(5.90594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53253,0.0632557,-4.2588e-05,-2.16731e-08,3.98272e-11,-34987.5,28.2196], Tmin=(100,'K'), Tmax=(507.35,'K')), NASAPolynomial(coeffs=[6.18267,0.0409947,-1.93496e-05,3.73841e-09,-2.62786e-13,-35644.7,7.1141], Tmin=(507.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-291.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
)

species(
    label = 'HCO(16)',
    structure = SMILES('[CH]=O'),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(O)C([O])(C=O)CC(6168)',
    structure = SMILES('C=C(O)C([O])(C=O)CC'),
    E0 = (-298.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1818.45,2694.67,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159255,'amu*angstrom^2'), symmetry=1, barrier=(3.66158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.608054,0.10466,-0.000128562,8.50037e-08,-2.25179e-11,-35713.8,33.3289], Tmin=(100,'K'), Tmax=(920.683,'K')), NASAPolynomial(coeffs=[15.433,0.0349675,-1.50146e-05,2.78299e-09,-1.91503e-13,-38667.5,-42.7339], Tmin=(920.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)[C](O)C=O(12676)',
    structure = SMILES('C=C(O)C([O])(CC)C(O)=C[O]'),
    E0 = (-431.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1717,0.114534,-7.78047e-05,-3.90975e-09,1.58711e-11,-51653.2,42.8198], Tmin=(100,'K'), Tmax=(952.416,'K')), NASAPolynomial(coeffs=[33.1022,0.0171892,-4.49938e-06,7.81357e-10,-6.0543e-14,-60676.4,-137.733], Tmin=(952.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-431.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C([O])C=O(12677)',
    structure = SMILES('C=C(O)C(O)([CH]C)C([O])C=O'),
    E0 = (-309.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180,180,526.44,609.615,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159393,'amu*angstrom^2'), symmetry=1, barrier=(3.66476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.44015,0.119172,-0.000119997,5.91932e-08,-1.11346e-11,-36950.8,51.6682], Tmin=(100,'K'), Tmax=(1417.39,'K')), NASAPolynomial(coeffs=[29.8434,0.0192512,-4.92489e-06,6.82021e-10,-4.05837e-14,-45217.1,-112.219], Tmin=(1417.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-309.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C(O)(CC)[C]([O])C=O(12678)',
    structure = SMILES('C=C(O)C(O)(CC)C([O])=C[O]'),
    E0 = (-522.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.37647,0.130569,-0.000135119,6.63021e-08,-1.21614e-11,-62582.7,47.4513], Tmin=(100,'K'), Tmax=(1516.59,'K')), NASAPolynomial(coeffs=[34.8051,0.0125997,-1.36419e-06,-3.4979e-13,5.46987e-15,-72178.3,-146.108], Tmin=(1516.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-522.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C(O)[C]=O(12679)',
    structure = SMILES('C=C(O)C([O])(CC)C(O)[C]=O'),
    E0 = (-363.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.84354,0.115638,-0.000108331,4.47434e-08,-5.23873e-12,-43543.6,48.2708], Tmin=(100,'K'), Tmax=(1006.52,'K')), NASAPolynomial(coeffs=[26.3166,0.02648,-9.36875e-06,1.65476e-09,-1.1472e-13,-50364.8,-93.4937], Tmin=(1006.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)C(O)C=O(12680)',
    structure = SMILES('C=C(O)C([O])([CH]C)C(O)C=O'),
    E0 = (-323.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57741,0.107472,-8.47062e-05,2.10604e-08,2.78943e-12,-38747.2,48.9936], Tmin=(100,'K'), Tmax=(995.845,'K')), NASAPolynomial(coeffs=[26.086,0.0270133,-9.69099e-06,1.75453e-09,-1.2442e-13,-45777,-91.9834], Tmin=(995.845,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)(C(=C)O)C([O])C=O(12681)',
    structure = SMILES('[CH2]CC(O)(C(=C)O)C([O])C=O'),
    E0 = (-303.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,559.058,577.855,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15876,'amu*angstrom^2'), symmetry=1, barrier=(3.65021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.27729,0.121454,-0.000126723,6.54899e-08,-1.30413e-11,-36318.9,50.3662], Tmin=(100,'K'), Tmax=(1266.91,'K')), NASAPolynomial(coeffs=[28.4548,0.0223008,-6.81272e-06,1.06881e-09,-6.8057e-14,-43935.5,-104.497], Tmin=(1266.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C([O])[C]=O(12682)',
    structure = SMILES('C=C(O)C(O)(CC)C([O])[C]=O'),
    E0 = (-349.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43375,0.124271,-0.000133554,7.06979e-08,-1.4326e-11,-41759.3,49.9574], Tmin=(100,'K'), Tmax=(1279.79,'K')), NASAPolynomial(coeffs=[29.2299,0.0201145,-5.39026e-06,7.65029e-10,-4.58295e-14,-49438.7,-108.952], Tmin=(1279.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-349.255,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C([O])C=O(12683)',
    structure = SMILES('C=C([O])C(O)(CC)C([O])C=O'),
    E0 = (-371.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71811,0.114073,-0.000114619,5.86826e-08,-1.17639e-11,-44454.6,48.0714], Tmin=(100,'K'), Tmax=(1222,'K')), NASAPolynomial(coeffs=[24.2567,0.0290495,-1.0254e-05,1.7462e-09,-1.15796e-13,-50802.9,-82.4498], Tmin=(1222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-371.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C([O])C=O(12684)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C([O])C=O'),
    E0 = (-262.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,180,1600,1746.05,2753.26,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155728,'amu*angstrom^2'), symmetry=1, barrier=(3.58049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98598,0.119593,-0.000117147,5.09111e-08,-6.46417e-12,-31299.4,47.9323], Tmin=(100,'K'), Tmax=(980.363,'K')), NASAPolynomial(coeffs=[27.0005,0.0248927,-8.31076e-06,1.42212e-09,-9.73176e-14,-38115.4,-97.1124], Tmin=(980.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-262.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([O])(C(=C)O)C(O)C=O(12685)',
    structure = SMILES('[CH2]CC([O])(C(=C)O)C(O)C=O'),
    E0 = (-318.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04201,0.11689,-0.000115149,5.63371e-08,-1.07109e-11,-38087.6,49.9602], Tmin=(100,'K'), Tmax=(1291.17,'K')), NASAPolynomial(coeffs=[27.4387,0.0255609,-9.04939e-06,1.55548e-09,-1.04035e-13,-45700.6,-99.8014], Tmin=(1291.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)C(O)C=O(12686)',
    structure = SMILES('C=C([O])C([O])(CC)C(O)C=O'),
    E0 = (-386.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.46956,0.109374,-0.000102664,4.91384e-08,-9.30223e-12,-46224,47.6161], Tmin=(100,'K'), Tmax=(1281.41,'K')), NASAPolynomial(coeffs=[23.2019,0.0323618,-1.25153e-05,2.23777e-09,-1.52123e-13,-52546.9,-77.5271], Tmin=(1281.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-386.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)C(O)C=O(12687)',
    structure = SMILES('[CH]=C(O)C([O])(CC)C(O)C=O'),
    E0 = (-276.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.11313,0.119242,-0.000119979,5.98934e-08,-1.16067e-11,-33052.4,48.8299], Tmin=(100,'K'), Tmax=(1268.01,'K')), NASAPolynomial(coeffs=[27.6247,0.0254339,-9.00848e-06,1.55027e-09,-1.03872e-13,-40594,-101.699], Tmin=(1268.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = 'CCC1(OC[C]1O)C([O])C=O(12688)',
    structure = SMILES('CCC1(OC[C]1O)C([O])C=O'),
    E0 = (-226.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63179,0.111785,-0.000112729,6.07132e-08,-1.28202e-11,-27071.6,43.4893], Tmin=(100,'K'), Tmax=(1224.92,'K')), NASAPolynomial(coeffs=[21.4534,0.0318268,-9.21455e-06,1.32706e-09,-7.76997e-14,-32384,-71.1662], Tmin=(1224.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(C=OCOJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C1[CH]OO1(12689)',
    structure = SMILES('C=C(O)C([O])(CC)C1[CH]OO1'),
    E0 = (-15.3284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.34917,0.113642,-0.00011131,5.68835e-08,-1.16689e-11,-1647.58,38.6782], Tmin=(100,'K'), Tmax=(1174.67,'K')), NASAPolynomial(coeffs=[20.4665,0.0393545,-1.64484e-05,3.04626e-09,-2.10897e-13,-6772.82,-70.0818], Tmin=(1174.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.3284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=CC(C)2OJ) + radical(CCsJOO)"""),
)

species(
    label = 'CCC1([O])[C](O)COC1C=O(12498)',
    structure = SMILES('CCC1([O])[C](O)COC1C=O'),
    E0 = (-320.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.855441,0.0935566,-5.2412e-05,-7.67446e-09,1.31187e-11,-38316.6,37.0963], Tmin=(100,'K'), Tmax=(925.936,'K')), NASAPolynomial(coeffs=[20.9109,0.0335358,-1.02728e-05,1.65248e-09,-1.09445e-13,-43805.3,-74.111], Tmin=(925.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-320.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(C2CsJOH) + radical(CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OO[CH]C1[O](12530)',
    structure = SMILES('C=C(O)C1(CC)OO[CH]C1[O]'),
    E0 = (-100.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.55036,0.122942,-0.000118647,5.62788e-08,-9.86424e-12,-11718.8,45.5187], Tmin=(100,'K'), Tmax=(1675.93,'K')), NASAPolynomial(coeffs=[28.6403,0.0187229,-8.55304e-07,-3.28447e-10,3.44797e-14,-18662.3,-114.929], Tmin=(1675.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCsJOOC) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C(O)(CC)C(=O)C=O(12690)',
    structure = SMILES('C=C(O)C(O)(CC)C(=O)C=O'),
    E0 = (-650.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94032,0.129496,-0.000150346,8.91927e-08,-2.08508e-11,-78031.2,39.5481], Tmin=(100,'K'), Tmax=(1046.39,'K')), NASAPolynomial(coeffs=[22.8302,0.0348045,-1.46029e-05,2.7075e-09,-1.87675e-13,-83215,-81.0782], Tmin=(1046.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H)"""),
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
    label = 'C=C(O)C(C)([O])C([O])C=O(12691)',
    structure = SMILES('C=C(O)C(C)([O])C([O])C=O'),
    E0 = (-256.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,180,180,180,554.265,583.221,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160032,'amu*angstrom^2'), symmetry=1, barrier=(3.67944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160032,'amu*angstrom^2'), symmetry=1, barrier=(3.67944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160032,'amu*angstrom^2'), symmetry=1, barrier=(3.67944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160032,'amu*angstrom^2'), symmetry=1, barrier=(3.67944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160032,'amu*angstrom^2'), symmetry=1, barrier=(3.67944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.831983,0.0943406,-8.09581e-05,2.74283e-08,-1.02586e-12,-30645.2,42.3576], Tmin=(100,'K'), Tmax=(1004.83,'K')), NASAPolynomial(coeffs=[22.6532,0.02329,-8.39084e-06,1.50614e-09,-1.05638e-13,-36497.7,-76.6949], Tmin=(1004.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=OCOJ)"""),
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
    label = 'C=C(O)C([O])(CC)C[O](5483)',
    structure = SMILES('C=C(O)C([O])(CC)C[O]'),
    E0 = (-181.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,350,440,435,1725,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1053.88,1168.22,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.164014,'amu*angstrom^2'), symmetry=1, barrier=(3.77101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164014,'amu*angstrom^2'), symmetry=1, barrier=(3.77101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164014,'amu*angstrom^2'), symmetry=1, barrier=(3.77101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164014,'amu*angstrom^2'), symmetry=1, barrier=(3.77101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164014,'amu*angstrom^2'), symmetry=1, barrier=(3.77101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4948.82,'J/mol'), sigma=(8.04162,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=772.99 K, Pc=21.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.300473,0.0946082,-9.2282e-05,4.87475e-08,-1.0534e-11,-21675.4,36.394], Tmin=(100,'K'), Tmax=(1106.33,'K')), NASAPolynomial(coeffs=[15.1617,0.0387046,-1.6487e-05,3.07458e-09,-2.1337e-13,-25096.7,-39.7647], Tmin=(1106.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)OOC1C=O(12340)',
    structure = SMILES('C=C(O)C1(CC)OOC1C=O'),
    E0 = (-352.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29143,0.0984332,-5.61118e-05,-7.923e-09,1.27268e-11,-42160.7,38.1222], Tmin=(100,'K'), Tmax=(984.517,'K')), NASAPolynomial(coeffs=[25.6038,0.0287527,-1.02697e-05,1.88239e-09,-1.35522e-13,-49375.3,-100.957], Tmin=(984.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-352.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxetane)"""),
)

species(
    label = 'C=[C]C(CC)(OO)C([O])C=O(12692)',
    structure = SMILES('C=[C]C(CC)(OO)C([O])C=O'),
    E0 = (14.6898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180,180,180,293.338,842.813,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152798,'amu*angstrom^2'), symmetry=1, barrier=(3.51312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19646,0.113767,-0.000115464,6.16157e-08,-1.3278e-11,1954.41,45.881], Tmin=(100,'K'), Tmax=(1116.07,'K')), NASAPolynomial(coeffs=[19.0343,0.0412586,-1.8011e-05,3.40297e-09,-2.38089e-13,-2561.32,-53.9421], Tmin=(1116.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.6898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C([O])(CC)C(C=O)OO(12693)',
    structure = SMILES('C=[C]C([O])(CC)C(C=O)OO'),
    E0 = (0.0593226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
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
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.896319,0.108518,-0.000101838,5.02285e-08,-1.01537e-11,182.698,45.2366], Tmin=(100,'K'), Tmax=(1171.79,'K')), NASAPolynomial(coeffs=[17.7227,0.0449604,-2.04777e-05,3.9398e-09,-2.7797e-13,-4180.78,-47.5408], Tmin=(1171.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.0593226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_S)"""),
)

species(
    label = 'CCC([O])(C(C)=O)C([O])C=O(12694)',
    structure = SMILES('CCC([O])(C(C)=O)C([O])C=O'),
    E0 = (-276.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.738628,0.110118,-0.000122572,7.76456e-08,-2.044e-11,-33090.8,44.5851], Tmin=(100,'K'), Tmax=(910.956,'K')), NASAPolynomial(coeffs=[13.1317,0.0492136,-2.22845e-05,4.25161e-09,-2.97869e-13,-35617.8,-21.0374], Tmin=(910.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(CC(C)(C=O)OJ) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)O[CH]C=O(12334)',
    structure = SMILES('C=C(O)C([O])(CC)OC=C[O]'),
    E0 = (-403.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(5298.89,'J/mol'), sigma=(8.39681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=827.67 K, Pc=20.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-5.51206,0.151582,-0.00016554,8.17514e-08,-1.45583e-11,-48180,51.0999], Tmin=(100,'K'), Tmax=(1664.79,'K')), NASAPolynomial(coeffs=[41.2097,-0.000334555,7.07124e-06,-1.68022e-09,1.19319e-13,-58240.7,-181.614], Tmin=(1664.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CC(C)(O)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C([O])=CO(12695)',
    structure = SMILES('C=C(O)C([O])(CC)C([O])=CO'),
    E0 = (-435.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,1011.34,1196.71,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161096,'amu*angstrom^2'), symmetry=1, barrier=(3.70392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161096,'amu*angstrom^2'), symmetry=1, barrier=(3.70392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161096,'amu*angstrom^2'), symmetry=1, barrier=(3.70392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161096,'amu*angstrom^2'), symmetry=1, barrier=(3.70392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161096,'amu*angstrom^2'), symmetry=1, barrier=(3.70392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161096,'amu*angstrom^2'), symmetry=1, barrier=(3.70392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.43189,0.131539,-0.00013633,6.68971e-08,-1.22681e-11,-52039.8,47.4542], Tmin=(100,'K'), Tmax=(1516.4,'K')), NASAPolynomial(coeffs=[35.1588,0.0123951,-1.31289e-06,-6.16114e-12,5.70543e-15,-61748.9,-148.213], Tmin=(1516.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-435.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C(O)[C](CC)C([O])C=O(12696)',
    structure = SMILES('[CH2]C(O)=C(CC)C([O])C=O'),
    E0 = (-179.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38908,0.106167,-0.000103889,5.16869e-08,-1.00719e-11,-21441.2,42.6044], Tmin=(100,'K'), Tmax=(1255.9,'K')), NASAPolynomial(coeffs=[23.3994,0.0272165,-9.59415e-06,1.63312e-09,-1.08233e-13,-27667.6,-82.6342], Tmin=(1255.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-179.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + radical(C=C(O)CJ) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C=O)CC(9323)',
    structure = SMILES('C=C(O)C([O])(C=C[O])CC'),
    E0 = (-216.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,180,942.414,1260.63,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.159601,'amu*angstrom^2'), symmetry=1, barrier=(3.66955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159601,'amu*angstrom^2'), symmetry=1, barrier=(3.66955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159601,'amu*angstrom^2'), symmetry=1, barrier=(3.66955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159601,'amu*angstrom^2'), symmetry=1, barrier=(3.66955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159601,'amu*angstrom^2'), symmetry=1, barrier=(3.66955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19184,0.0960429,-5.60068e-05,-7.81987e-09,1.27722e-11,-25886.2,39.2459], Tmin=(100,'K'), Tmax=(985.188,'K')), NASAPolynomial(coeffs=[26.1747,0.025039,-8.96606e-06,1.6711e-09,-1.223e-13,-33224.9,-102.252], Tmin=(985.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(C=COJ)"""),
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
    E0 = (-280.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-232.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-177.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-187.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-158.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-162.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-230.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-280.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-140.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-280.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-171.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-126.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-203.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-204.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-204.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-226.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-220.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-152.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-178.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-217.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-248.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-242.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-168.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-203.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-155.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-14.3265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-218.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-100.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-216.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (163.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (4.58029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-271.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (121.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (112.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-75.9339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-90.0526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-280.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (63.0323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (26.0858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['OCHCHO(48)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['[CH2]C1(O)OC1(CC)C([O])C=O(12672)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C([O])(CC)C1OC1[O](12673)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['[CH2]C1(O)OC(C=O)C1([O])CC(12446)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(92.4015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_O] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C1(CC)OC([O])C1[O](12660)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C(O)C([O])(CC)C(=O)C=O(12674)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5(29)', 'C=C(O)C(=O)C([O])C=O(12675)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C[O](2537)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.94e+07,'cm^3/(mol*s)'), n=1.88, Ea=(97.2008,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 93.4 to 97.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2COH(99)', 'CCC(=O)C([O])C=O(11465)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OCHCHO(48)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(130.13,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 124.5 to 130.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['HCO(16)', 'C=C(O)C([O])(C=O)CC(6168)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C([O])(CC)[C](O)C=O(12676)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C(O)([CH]C)C([O])C=O(12677)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C(O)(CC)[C]([O])C=O(12678)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C([O])(CC)C(O)[C]=O(12679)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C([O])([CH]C)C(O)C=O(12680)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC(O)(C(=C)O)C([O])C=O(12681)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C(O)(CC)C([O])[C]=O(12682)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C([O])C(O)(CC)C([O])C=O(12683)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C(O)(CC)C([O])C=O(12684)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['[CH2]CC([O])(C(=C)O)C(O)C=O(12685)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C([O])C([O])(CC)C(O)C=O(12686)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(22219.5,'s^-1'), n=1.84103, Ea=(37.7069,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_CCC;O_rad_out;XH_out] + [R5H_CCC(Cd);Y_rad_out;XH_out] for rate rule [R5H_CCC(Cd);O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['[CH]=C(O)C([O])(CC)C(O)C=O(12687)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C=C[O](2537)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['CCC1(OC[C]1O)C([O])C=O(12688)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C([O])(CC)C1[CH]OO1(12689)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['CCC1([O])[C](O)COC1C=O(12498)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.90996e+09,'s^-1'), n=0.623333, Ea=(61.7837,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C1(CC)OO[CH]C1[O](12530)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.81184e+09,'s^-1'), n=0.551229, Ea=(180.107,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 179.5 to 180.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C(O)(CC)C(=O)C=O(12690)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)([O])C([O])C=O(12691)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CO(12)', 'C=C(O)C([O])(CC)C[O](5483)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['C=C(O)C1(CC)OOC1C=O(12340)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(CC)(OO)C([O])C=O(12692)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C([O])(CC)C(C=O)OO(12693)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.95074e+10,'s^-1'), n=0, Ea=(112.549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OOH_SS;Y_rad_out] for rate rule [R3OOH_SS;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    products = ['CCC([O])(C(C)=O)C([O])C=O(12694)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)C([O])(CC)O[CH]C=O(12334)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C(O)C([O])(CC)C([O])=CO(12695)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(155.022,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 153.4 to 155.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['O(4)', 'C=C(O)[C](CC)C([O])C=O(12696)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['O(4)', 'C=C(O)C([O])([CH]C=O)CC(9323)'],
    products = ['C=C(O)C([O])(CC)C([O])C=O(12330)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '2267',
    isomers = [
        'C=C(O)C([O])(CC)C([O])C=O(12330)',
    ],
    reactants = [
        ('OCHCHO(48)', 'C=C(O)C(=O)CC(4626)'),
        ('OCHCHO(48)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2267',
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

