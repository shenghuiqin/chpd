species(
    label = 'C=C(O)C([O])(CC)C[CH]C(6359)',
    structure = SMILES('C=C(O)C([O])(CC)C[CH]C'),
    E0 = (-100.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1025.12,1193.47,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37847,0.111562,-9.67436e-05,4.49979e-08,-8.52632e-12,-11888.1,43.919], Tmin=(100,'K'), Tmax=(1257.13,'K')), NASAPolynomial(coeffs=[19.371,0.0455405,-1.79667e-05,3.2218e-09,-2.18469e-13,-17105,-60.9332], Tmin=(1257.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C3H6(72)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,-487.138,-4.54468], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C1(O)OC1(CC)C[CH]C(6971)',
    structure = SMILES('[CH2]C1(O)OC1(CC)C[CH]C'),
    E0 = (-54.2967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.31745,0.132569,-0.000131392,6.64025e-08,-1.26959e-11,-6242.63,43.8056], Tmin=(100,'K'), Tmax=(1473.54,'K')), NASAPolynomial(coeffs=[28.6801,0.0268144,-4.50359e-06,2.92598e-10,-3.31946e-15,-13621.1,-116.008], Tmin=(1473.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.2967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CJC(O)2C) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1(O)C(C)CC1([O])CC(6725)',
    structure = SMILES('[CH2]C1(O)C(C)CC1([O])CC'),
    E0 = (-30.2655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68847,0.109483,-7.41333e-05,1.37695e-08,3.72372e-12,-3421.71,37.0526], Tmin=(100,'K'), Tmax=(1029.75,'K')), NASAPolynomial(coeffs=[23.5548,0.0402161,-1.517e-05,2.74552e-09,-1.91097e-13,-10146.9,-92.8828], Tmin=(1029.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.2655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CC(C)2OJ)"""),
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
    label = 'C=C(O)C([O])(C=CC)CC(6972)',
    structure = SMILES('C=C(O)C([O])(C=CC)CC'),
    E0 = (-185.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,450.749,689.726,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156249,'amu*angstrom^2'), symmetry=1, barrier=(3.59247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156249,'amu*angstrom^2'), symmetry=1, barrier=(3.59247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156249,'amu*angstrom^2'), symmetry=1, barrier=(3.59247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156249,'amu*angstrom^2'), symmetry=1, barrier=(3.59247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156249,'amu*angstrom^2'), symmetry=1, barrier=(3.59247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156249,'amu*angstrom^2'), symmetry=1, barrier=(3.59247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48641,0.105689,-8.36706e-05,3.34415e-08,-5.33828e-12,-22114.3,40.8626], Tmin=(100,'K'), Tmax=(1492.24,'K')), NASAPolynomial(coeffs=[24.1007,0.0371022,-1.47266e-05,2.64029e-09,-1.78036e-13,-29750.7,-92.8221], Tmin=(1492.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=CCC([O])(CC)C(=C)O(6973)',
    structure = SMILES('C=CCC([O])(CC)C(=C)O'),
    E0 = (-167.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1035.82,1177.78,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161961,'amu*angstrom^2'), symmetry=1, barrier=(3.72379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161961,'amu*angstrom^2'), symmetry=1, barrier=(3.72379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161961,'amu*angstrom^2'), symmetry=1, barrier=(3.72379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161961,'amu*angstrom^2'), symmetry=1, barrier=(3.72379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161961,'amu*angstrom^2'), symmetry=1, barrier=(3.72379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161961,'amu*angstrom^2'), symmetry=1, barrier=(3.72379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.47803,0.10604,-7.64759e-05,1.85024e-08,1.89866e-12,-19895.3,42.0146], Tmin=(100,'K'), Tmax=(1029.45,'K')), NASAPolynomial(coeffs=[22.9602,0.0370836,-1.3884e-05,2.50202e-09,-1.73739e-13,-26304.6,-83.2865], Tmin=(1029.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CC(C)2OJ)"""),
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
    label = 'C=C(O)C(=O)C[CH]C(6974)',
    structure = SMILES('C=C(O)C(=O)C[CH]C'),
    E0 = (-182.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242852,0.0809486,-7.46424e-05,3.61903e-08,-7.13832e-12,-21794.5,29.6354], Tmin=(100,'K'), Tmax=(1206.94,'K')), NASAPolynomial(coeffs=[15.0055,0.0320234,-1.38385e-05,2.6051e-09,-1.81754e-13,-25358.1,-44.363], Tmin=(1206.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CCJCC=O)"""),
)

species(
    label = 'C3H6(T)(143)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238389,'amu*angstrom^2'), symmetry=1, barrier=(5.48103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26842e-06,-1.44873e-08,5.74941e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.81,'K')), NASAPolynomial(coeffs=[5.93909,0.0171892,-6.69152e-06,1.21546e-09,-8.39795e-14,33151.2,-4.14888], Tmin=(1046.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C[CH]CC(=O)CC(6975)',
    structure = SMILES('C[CH]CC(=O)CC'),
    E0 = (-107.489,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1685,0.0677878,-5.39461e-05,2.83459e-08,-7.19769e-12,-12830.7,27.2394], Tmin=(100,'K'), Tmax=(870.951,'K')), NASAPolynomial(coeffs=[5.13466,0.0495726,-2.25751e-05,4.33321e-09,-3.05084e-13,-13521.6,8.65301], Tmin=(870.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CCJCC=O)"""),
)

species(
    label = 'C=C(O)C([O])([CH]CC)CC(6800)',
    structure = SMILES('C=C(O)C([O])([CH]CC)CC'),
    E0 = (-95.0432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1064.56,1147.69,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60905,0.108824,-7.81789e-05,1.95543e-08,1.52878e-12,-11216.5,44.5071], Tmin=(100,'K'), Tmax=(1033.85,'K')), NASAPolynomial(coeffs=[22.9346,0.0394958,-1.47805e-05,2.65393e-09,-1.83543e-13,-17661.2,-81.3439], Tmin=(1033.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]CCC([O])(CC)C(=C)O(6976)',
    structure = SMILES('[CH2]CCC([O])(CC)C(=C)O'),
    E0 = (-89.6989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,180,180,180,180,1012.85,1201.66,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90868,0.116318,-0.000101977,4.63208e-08,-8.38509e-12,-10564.2,44.8803], Tmin=(100,'K'), Tmax=(1332.67,'K')), NASAPolynomial(coeffs=[24.0053,0.0385359,-1.44272e-05,2.52359e-09,-1.6889e-13,-17471,-87.5811], Tmin=(1332.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.6989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(RCCJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)C[CH]C(6977)',
    structure = SMILES('C=C(O)C(O)([CH]C)C[CH]C'),
    E0 = (-129.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,538.322,600.994,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99999,0.117727,-0.000108969,5.28002e-08,-1.01255e-11,-15371.2,47.435], Tmin=(100,'K'), Tmax=(1270.04,'K')), NASAPolynomial(coeffs=[23.9978,0.0358442,-1.2258e-05,2.03341e-09,-1.321e-13,-21974.7,-84.2032], Tmin=(1270.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = 'C=C(O)C(O)([CH][CH]C)CC(6978)',
    structure = SMILES('C=C(O)C(O)([CH][CH]C)CC'),
    E0 = (-129.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99999,0.117727,-0.000108969,5.28002e-08,-1.01255e-11,-15371.2,47.435], Tmin=(100,'K'), Tmax=(1270.04,'K')), NASAPolynomial(coeffs=[23.9978,0.0358442,-1.2258e-05,2.03341e-09,-1.321e-13,-21974.7,-84.2032], Tmin=(1270.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)CCC(6979)',
    structure = SMILES('C=C(O)C([O])([CH]C)CCC'),
    E0 = (-95.0432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1064.56,1147.69,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60905,0.108824,-7.81789e-05,1.95543e-08,1.52878e-12,-11216.5,44.5071], Tmin=(100,'K'), Tmax=(1033.85,'K')), NASAPolynomial(coeffs=[22.9346,0.0394958,-1.47805e-05,2.65393e-09,-1.83543e-13,-17661.2,-81.3439], Tmin=(1033.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)(C[CH]C)C(=C)O(6980)',
    structure = SMILES('[CH2]CC(O)(C[CH]C)C(=C)O'),
    E0 = (-124.355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,498.316,642.516,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8536,0.12016,-0.000116071,5.94474e-08,-1.21461e-11,-14738.5,46.1957], Tmin=(100,'K'), Tmax=(1187.86,'K')), NASAPolynomial(coeffs=[21.9915,0.0398639,-1.4675e-05,2.54033e-09,-1.69249e-13,-20403.4,-72.9481], Tmin=(1187.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C[CH]C(6981)',
    structure = SMILES('C=C([O])C(O)(CC)C[CH]C'),
    E0 = (-191.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31921,0.113029,-0.000104673,5.33633e-08,-1.11128e-11,-22873,43.993], Tmin=(100,'K'), Tmax=(1150.43,'K')), NASAPolynomial(coeffs=[17.7632,0.0466791,-1.8161e-05,3.2294e-09,-2.18026e-13,-27263.5,-50.7425], Tmin=(1150.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C[CH]C(6982)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C[CH]C'),
    E0 = (-82.5055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1600,1686.15,2816.21,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94568,0.122732,-0.000121557,6.3712e-08,-1.32909e-11,-9702.21,45.1425], Tmin=(100,'K'), Tmax=(1165.81,'K')), NASAPolynomial(coeffs=[22.2955,0.0395585,-1.45402e-05,2.51448e-09,-1.67467e-13,-15354.3,-75.526], Tmin=(1165.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.5055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CC([O])(CCC)C(=C)O(6983)',
    structure = SMILES('[CH2]CC([O])(CCC)C(=C)O'),
    E0 = (-89.6989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,180,180,180,180,1012.85,1201.66,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160573,'amu*angstrom^2'), symmetry=1, barrier=(3.69188,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.90868,0.116318,-0.000101977,4.63208e-08,-8.38509e-12,-10564.2,44.8803], Tmin=(100,'K'), Tmax=(1332.67,'K')), NASAPolynomial(coeffs=[24.0053,0.0385359,-1.44272e-05,2.52359e-09,-1.6889e-13,-17471,-87.5811], Tmin=(1332.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.6989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C([O])C([O])(CC)CCC(6984)',
    structure = SMILES('C=C([O])C([O])(CC)CCC'),
    E0 = (-157.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32932,0.108737,-8.9325e-05,3.8969e-08,-6.93056e-12,-18700.9,42.5102], Tmin=(100,'K'), Tmax=(1331.89,'K')), NASAPolynomial(coeffs=[19.8016,0.0452767,-1.78569e-05,3.19709e-09,-2.16233e-13,-24329.8,-65.4909], Tmin=(1331.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])(CC)CCC(6985)',
    structure = SMILES('[CH]=C(O)C([O])(CC)CCC'),
    E0 = (-47.8492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,350,440,435,1725,3120,650,792.5,1650,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,516.027,625.095,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157596,'amu*angstrom^2'), symmetry=1, barrier=(3.62343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97579,0.118625,-0.000106666,4.97196e-08,-9.22462e-12,-5529.07,43.7355], Tmin=(100,'K'), Tmax=(1302.8,'K')), NASAPolynomial(coeffs=[24.1084,0.0385385,-1.44567e-05,2.53429e-09,-1.70003e-13,-12325.5,-89.0052], Tmin=(1302.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.8492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]CC(O)(CC)C(=C)O(6986)',
    structure = SMILES('[CH2][CH]CC(O)(CC)C(=C)O'),
    E0 = (-124.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8536,0.12016,-0.000116071,5.94474e-08,-1.21461e-11,-14738.5,46.1957], Tmin=(100,'K'), Tmax=(1187.86,'K')), NASAPolynomial(coeffs=[21.9915,0.0398639,-1.4675e-05,2.54033e-09,-1.69249e-13,-20403.4,-72.9481], Tmin=(1187.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC1(CC)OC[C]1O(6987)',
    structure = SMILES('C[CH]CC1(CC)OC[C]1O'),
    E0 = (-47.2436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29167,0.111396,-0.00010489,5.7857e-08,-1.31037e-11,-5487.33,39.6242], Tmin=(100,'K'), Tmax=(1065.97,'K')), NASAPolynomial(coeffs=[15.4045,0.0487445,-1.67288e-05,2.7205e-09,-1.7267e-13,-9046.87,-41.9921], Tmin=(1065.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.2436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(RCCJC)"""),
)

species(
    label = 'CCC1([O])CC(C)C[C]1O(6702)',
    structure = SMILES('CCC1([O])CC(C)C[C]1O'),
    E0 = (-129.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.919339,0.0915417,-3.46749e-05,-1.95757e-08,1.39152e-11,-15376,36.9953], Tmin=(100,'K'), Tmax=(1010.74,'K')), NASAPolynomial(coeffs=[19.9882,0.0440714,-1.65705e-05,3.00813e-09,-2.10306e-13,-21404,-73.0075], Tmin=(1010.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(565.384,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=C(O)C(O)(C=CC)CC(6988)',
    structure = SMILES('C=C(O)C(O)(C=CC)CC'),
    E0 = (-414.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69317,0.110236,-8.22759e-05,2.41786e-08,-4.03963e-13,-49661.2,40.3709], Tmin=(100,'K'), Tmax=(1068.41,'K')), NASAPolynomial(coeffs=[23.6251,0.0388845,-1.50046e-05,2.73361e-09,-1.90067e-13,-56408.9,-89.7104], Tmin=(1068.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(O)(CC)C(=C)O(6989)',
    structure = SMILES('C=CCC(O)(CC)C(=C)O'),
    E0 = (-396.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.79539,0.111614,-7.73566e-05,1.03135e-08,7.16615e-12,-47436.9,41.938], Tmin=(100,'K'), Tmax=(973.27,'K')), NASAPolynomial(coeffs=[25.5103,0.0341835,-1.16427e-05,2.03036e-09,-1.40733e-13,-54399.9,-97.5218], Tmin=(973.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
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
    label = 'C=C(O)C(C)([O])C[CH]C(6990)',
    structure = SMILES('C=C(O)C(C)([O])C[CH]C'),
    E0 = (-76.7188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,465.963,677.395,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156672,'amu*angstrom^2'), symmetry=1, barrier=(3.6022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156672,'amu*angstrom^2'), symmetry=1, barrier=(3.6022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156672,'amu*angstrom^2'), symmetry=1, barrier=(3.6022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156672,'amu*angstrom^2'), symmetry=1, barrier=(3.6022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156672,'amu*angstrom^2'), symmetry=1, barrier=(3.6022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156672,'amu*angstrom^2'), symmetry=1, barrier=(3.6022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.789369,0.0974266,-8.50323e-05,3.95322e-08,-7.42844e-12,-9048.14,39.5609], Tmin=(100,'K'), Tmax=(1274.9,'K')), NASAPolynomial(coeffs=[18.3898,0.0372519,-1.4233e-05,2.51012e-09,-1.68647e-13,-13938.5,-57.6257], Tmin=(1274.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.7188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C(C)C([O])(CC)C(=C)O(6363)',
    structure = SMILES('[CH2]C(C)C([O])(CC)C(=C)O'),
    E0 = (-95.5573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,337.07,796.543,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153643,'amu*angstrom^2'), symmetry=1, barrier=(3.53256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4683.75,'J/mol'), sigma=(8.00266,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=731.59 K, Pc=20.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69604,0.112314,-8.62121e-05,2.53836e-08,4.95704e-13,-11276.4,43.8231], Tmin=(100,'K'), Tmax=(993.19,'K')), NASAPolynomial(coeffs=[22.7469,0.0391271,-1.38224e-05,2.396e-09,-1.62646e-13,-17377.3,-80.2039], Tmin=(993.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.5573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C1(CC)CC(C)O1(6369)',
    structure = SMILES('C=C(O)C1(CC)CC(C)O1'),
    E0 = (-352.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06185,0.0873362,5.10714e-06,-8.72703e-08,4.66795e-11,-42210.3,35.2076], Tmin=(100,'K'), Tmax=(900.24,'K')), NASAPolynomial(coeffs=[26.9604,0.0280298,-4.71912e-06,4.61965e-10,-2.68801e-14,-49897.8,-111.714], Tmin=(900.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-352.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = 'C=[C]C(CC)(C[CH]C)OO(6991)',
    structure = SMILES('C=[C]C(CC)(C[CH]C)OO'),
    E0 = (194.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,966.943,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151363,'amu*angstrom^2'), symmetry=1, barrier=(3.48014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.983594,0.114725,-0.000111695,6.32985e-08,-1.52453e-11,23544.5,42.4835], Tmin=(100,'K'), Tmax=(980.006,'K')), NASAPolynomial(coeffs=[12.7667,0.0586026,-2.57949e-05,4.86442e-09,-3.38978e-13,20849.4,-23.5761], Tmin=(980.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]CC([O])(CC)C(C)=O(6992)',
    structure = SMILES('C[CH]CC([O])(CC)C(C)=O'),
    E0 = (-96.8951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.935063,0.11629,-0.00013891,1.08275e-07,-3.61594e-11,-11483.4,42.6302], Tmin=(100,'K'), Tmax=(776.344,'K')), NASAPolynomial(coeffs=[8.20771,0.0642143,-2.86924e-05,5.38369e-09,-3.71166e-13,-12753.3,1.80061], Tmin=(776.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.8951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CC(C)(C=O)OJ)"""),
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
    label = 'C=C(O)[C](CC)C[CH]C(6993)',
    structure = SMILES('[CH2]C(O)=C(CC)C[CH]C'),
    E0 = (0.776236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    molecularWeight = (126.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16185,0.104189,-8.92945e-05,4.12023e-08,-7.702e-12,286.975,39.8092], Tmin=(100,'K'), Tmax=(1281.17,'K')), NASAPolynomial(coeffs=[19.1177,0.0408724,-1.51627e-05,2.62688e-09,-1.74539e-13,-4909.29,-63.0525], Tmin=(1281.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.776236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(RCCJC)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C([O])(CC)C(=C)O(6187)',
    structure = SMILES('[CH2]C([O])(CC)C(=C)O'),
    E0 = (-33.9413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1602.72,2893.45,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(3.42758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.676536,0.0908249,-8.54045e-05,4.10682e-08,-7.75587e-12,-3903.68,35.0552], Tmin=(100,'K'), Tmax=(1293.32,'K')), NASAPolynomial(coeffs=[20.7584,0.0245293,-8.51315e-06,1.43237e-09,-9.40947e-14,-9448.05,-73.8692], Tmin=(1293.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.9413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)2OJ)"""),
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
    E0 = (-100.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-52.8159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-30.2655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (31.9097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (47.1395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-46.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-41.3857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (43.9381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-100.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (63.9488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (57.9963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-24.1376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-25.2184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-48.3079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-40.6753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1.30858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-38.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-43.2565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-71.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-9.62641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-40.6051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (99.9221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (24.3934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-28.5342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-37.0989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-75.5258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (343.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (64.3779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-92.2147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (301.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (103.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (243.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (309.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C3H6(72)', 'C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['[CH2]C1(O)OC1(CC)C[CH]C(6971)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.52e+08,'s^-1'), n=0.89, Ea=(47.6831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['[CH2]C1(O)C(C)CC1([O])CC(6725)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.53628e+06,'s^-1'), n=1.28, Ea=(70.2335,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_csHNd] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 70.1 to 70.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C([O])(C=CC)CC(6972)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CCC([O])(CC)C(=C)O(6973)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H5(29)', 'C=C(O)C(=O)C[CH]C(6974)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.94e+10,'cm^3/(mol*s)'), n=0, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(333,'K'), Tmax=(363,'K'), comment="""Estimated using template [CO_O;CsJ-CsHH] for rate rule [CO-CdCs_O;CsJ-CsHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C3H6(T)(143)', 'C=C(O)C(=O)CC(4626)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.97e+07,'cm^3/(mol*s)'), n=1.88, Ea=(32.2168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdCs_O;YJ] for rate rule [CO-CdCs_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', 'C[CH]CC(=O)CC(6975)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.16e+10,'cm^3/(mol*s)'), n=0, Ea=(48.1578,'kJ/mol'), T0=(1,'K'), Tmin=(413,'K'), Tmax=(563,'K'), comment="""Estimated using template [CO-CsCs_O;CJ] for rate rule [CO-CsCs_O;CdsJ-O2s]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C3H6(72)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(78.4677,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 73.8 to 78.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C([O])([CH]CC)CC(6800)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CCC([O])(CC)C(=C)O(6976)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C=C(O)C(O)([CH][CH]C)CC(6978)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C(O)C([O])([CH]C)CCC(6979)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00228,'s^-1'), n=3.95, Ea=(46.7353,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC(O)(C[CH]C)C(=C)O(6980)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.6e-09,'s^-1'), n=5.55, Ea=(83.68,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 340 used for R4H_SSS;C_rad_out_2H;O_H_out
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C=C([O])C(O)(CC)C[CH]C(6981)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1070.11,'s^-1'), n=2.50856, Ea=(101.808,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SSS;O_rad_out;XH_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;O_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C(O)C(O)(CC)C[CH]C(6982)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC([O])(CCC)C(=C)O(6983)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 112 used for R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C=C([O])C([O])(CC)CCC(6984)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(O)C([O])(CC)CCC(6985)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(315594,'s^-1'), n=1.73223, Ea=(38.2227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['[CH2][CH]CC(O)(CC)C(=C)O(6986)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C3H6(T)(143)', '[CH2]C(O)=C([O])CC(4557)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C[CH]CC1(CC)OC[C]1O(6987)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.21748e+08,'s^-1'), n=0.95, Ea=(124.892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R4_S_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['CCC1([O])CC(C)C[C]1O(6702)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.51884e+07,'s^-1'), n=0.875, Ea=(71.9648,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra;radadd_intra_csHCs] + [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C=C(O)C(O)(C=CC)CC(6988)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C=CCC(O)(CC)C(=C)O(6989)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(S)(23)', 'C=C(O)C(C)([O])C[CH]C(6990)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C)C([O])(CC)C(=C)O(6363)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C=C(O)C1(CC)CC(C)O1(6369)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(CC)(C[CH]C)OO(6991)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.01978e+11,'s^-1'), n=0, Ea=(107.111,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OOH_S;Y_rad_out] for rate rule [R2OOH_S;Cd_rad_out_double]
Euclidian distance = 2.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    products = ['C[CH]CC([O])(CC)C(C)=O(6992)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', 'C=C(O)[C](CC)C[CH]C(6993)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/Cs2;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CHCH3(T)(95)', '[CH2]C([O])(CC)C(=C)O(6187)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1396',
    isomers = [
        'C=C(O)C([O])(CC)C[CH]C(6359)',
    ],
    reactants = [
        ('C3H6(72)', 'C=C(O)C(=O)CC(4626)'),
        ('C3H6(72)', '[CH2]C(O)=C([O])CC(4557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1396',
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

