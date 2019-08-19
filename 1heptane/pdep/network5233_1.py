species(
    label = '[CH2]OC(O)([CH]C)C(=C)O(5955)',
    structure = SMILES('[CH2]OC(O)([CH]C)C(=C)O'),
    E0 = (-232.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,425.199,700.913,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155275,'amu*angstrom^2'), symmetry=1, barrier=(3.57008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.12817,0.119085,-0.000124578,6.0517e-08,-1.07939e-11,-27627.8,44.1373], Tmin=(100,'K'), Tmax=(1617.67,'K')), NASAPolynomial(coeffs=[32.9242,0.00633926,1.84945e-06,-6.04201e-10,4.56621e-14,-36204.1,-137.591], Tmin=(1617.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCH3) + radical(CCJCO)"""),
)

species(
    label = 'CH2O(15)',
    structure = SMILES('C=O'),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=C(O)C(O)=CC(5562)',
    structure = SMILES('C=C(O)C(O)=CC'),
    E0 = (-369.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.87493,'amu*angstrom^2'), symmetry=1, barrier=(20.1164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.877997,'amu*angstrom^2'), symmetry=1, barrier=(20.1869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874595,'amu*angstrom^2'), symmetry=1, barrier=(20.1087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875946,'amu*angstrom^2'), symmetry=1, barrier=(20.1397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884099,0.0877395,-9.34163e-05,4.76671e-08,-9.0713e-12,-44294.8,25.8217], Tmin=(100,'K'), Tmax=(1473.72,'K')), NASAPolynomial(coeffs=[23.5689,0.00887539,-4.29798e-07,-1.49469e-10,1.60597e-14,-50145.5,-97.0299], Tmin=(1473.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]OC1(O)C(C)C1([CH2])O(32214)',
    structure = SMILES('[CH2]OC1(O)C(C)C1([CH2])O'),
    E0 = (-168.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53314,0.0993855,-5.33192e-05,-3.02489e-08,2.70067e-11,-20104.6,33.8748], Tmin=(100,'K'), Tmax=(919.771,'K')), NASAPolynomial(coeffs=[33.1764,0.00716958,1.28703e-06,-4.0297e-10,2.40791e-14,-28973.9,-144.181], Tmin=(919.771,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CsJOCH3) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1(O)COC1(O)[CH]C(32215)',
    structure = SMILES('[CH2]C1(O)COC1(O)[CH]C'),
    E0 = (-171.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.49682,0.116805,-0.000120167,5.9054e-08,-1.04876e-11,-20255.7,42.398], Tmin=(100,'K'), Tmax=(1700.02,'K')), NASAPolynomial(coeffs=[27.2597,0.0105629,3.46478e-06,-1.14965e-09,8.94398e-14,-25818,-107.908], Tmin=(1700.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)2O) + radical(CCJCO)"""),
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
    label = '[CH2]OC(O)(C=C)C(=C)O(32216)',
    structure = SMILES('[CH2]OC(O)(C=C)C(=C)O'),
    E0 = (-305.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,180,180,431.931,698.89,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(3.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(3.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(3.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(3.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(3.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155489,'amu*angstrom^2'), symmetry=1, barrier=(3.575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.74011,0.116405,-0.00012531,6.25157e-08,-1.14922e-11,-36457.4,40.0415], Tmin=(100,'K'), Tmax=(1549.16,'K')), NASAPolynomial(coeffs=[32.4515,0.00537274,1.7259e-06,-5.5627e-10,4.22993e-14,-44941,-137.331], Tmin=(1549.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-305.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]C(O)=C(O)[CH]C(4609)',
    structure = SMILES('[CH2]C(O)=C(O)[CH]C'),
    E0 = (-122.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,322.653],'cm^-1')),
        HinderedRotor(inertia=(0.00160507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208419,'amu*angstrom^2'), symmetry=1, barrier=(15.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20977,'amu*angstrom^2'), symmetry=1, barrier=(15.5164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209392,'amu*angstrom^2'), symmetry=1, barrier=(15.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207973,'amu*angstrom^2'), symmetry=1, barrier=(15.5093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141743,0.0688338,-2.99459e-05,-3.072e-08,2.30897e-11,-14621,29.0351], Tmin=(100,'K'), Tmax=(905.863,'K')), NASAPolynomial(coeffs=[24.2058,0.00628117,1.26066e-06,-4.23747e-10,2.91719e-14,-20774,-94.5791], Tmin=(905.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCJCO)"""),
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
    label = '[CH2]OC(O)=CC(32217)',
    structure = SMILES('[CH2]OC(O)=CC'),
    E0 = (-149.114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,203.02,203.066],'cm^-1')),
        HinderedRotor(inertia=(0.00408986,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516333,'amu*angstrom^2'), symmetry=1, barrier=(15.0959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516697,'amu*angstrom^2'), symmetry=1, barrier=(15.0963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516437,'amu*angstrom^2'), symmetry=1, barrier=(15.0944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.913578,0.0602387,-4.85617e-05,1.30832e-08,1.33132e-12,-17816.3,22.0092], Tmin=(100,'K'), Tmax=(970.649,'K')), NASAPolynomial(coeffs=[15.7672,0.0156046,-5.20359e-06,8.98498e-10,-6.21111e-14,-21480.7,-53.2313], Tmin=(970.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COCJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]OC(=CC)C(=C)O(32218)',
    structure = SMILES('[CH2]OC(=CC)C(=C)O'),
    E0 = (-159.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1086,0.0952661,-9.71941e-05,4.86991e-08,-9.30675e-12,-18948.7,30.7499], Tmin=(100,'K'), Tmax=(1393.23,'K')), NASAPolynomial(coeffs=[24.3204,0.0154245,-3.87597e-06,5.25264e-10,-3.07111e-14,-25371.1,-97.983], Tmin=(1393.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-159.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]CC(O)(O[CH2])C(=C)O(5957)',
    structure = SMILES('[CH2]CC(O)(O[CH2])C(=C)O'),
    E0 = (-226.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,350,440,435,1725,180,180,180,402.603,723.927,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154828,'amu*angstrom^2'), symmetry=1, barrier=(3.55979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.90019,0.120701,-0.000129375,6.47816e-08,-1.20025e-11,-26999,42.594], Tmin=(100,'K'), Tmax=(1526.6,'K')), NASAPolynomial(coeffs=[32.9885,0.0073484,9.82295e-07,-4.34004e-10,3.47321e-14,-35705.7,-138.358], Tmin=(1526.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CsJOCH3)"""),
)

species(
    label = '[CH2]OC([O])(CC)C(=C)O(5515)',
    structure = SMILES('[CH2]OC([O])(CC)C(=C)O'),
    E0 = (-198.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1099.37,1115.66,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.163799,'amu*angstrom^2'), symmetry=1, barrier=(3.76606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163799,'amu*angstrom^2'), symmetry=1, barrier=(3.76606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163799,'amu*angstrom^2'), symmetry=1, barrier=(3.76606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163799,'amu*angstrom^2'), symmetry=1, barrier=(3.76606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163799,'amu*angstrom^2'), symmetry=1, barrier=(3.76606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163799,'amu*angstrom^2'), symmetry=1, barrier=(3.76606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4530.1,'J/mol'), sigma=(7.54217,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=707.59 K, Pc=23.96 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.01429,0.114418,-0.000124342,6.58463e-08,-1.32313e-11,-23698.7,38.4921], Tmin=(100,'K'), Tmax=(1329.16,'K')), NASAPolynomial(coeffs=[28.0194,0.0156044,-3.31397e-06,3.70549e-10,-1.85318e-14,-30938,-112.148], Tmin=(1329.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCH3) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)OC(5956)',
    structure = SMILES('C=C(O)C([O])([CH]C)OC'),
    E0 = (-186.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,440,435,1725,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1065.22,1151.08,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.165296,'amu*angstrom^2'), symmetry=1, barrier=(3.80048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165296,'amu*angstrom^2'), symmetry=1, barrier=(3.80048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165296,'amu*angstrom^2'), symmetry=1, barrier=(3.80048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165296,'amu*angstrom^2'), symmetry=1, barrier=(3.80048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165296,'amu*angstrom^2'), symmetry=1, barrier=(3.80048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165296,'amu*angstrom^2'), symmetry=1, barrier=(3.80048,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.78443,0.108504,-0.000113186,5.8035e-08,-1.13199e-11,-22259.3,39.3765], Tmin=(100,'K'), Tmax=(1374.61,'K')), NASAPolynomial(coeffs=[26.7615,0.0168676,-3.83855e-06,4.67273e-10,-2.5132e-14,-29297.5,-104.478], Tmin=(1374.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]OC(O)(CC)C(=C)[O](5958)',
    structure = SMILES('[CH2]OC(O)(CC)C(=C)[O]'),
    E0 = (-294.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.35013,0.113423,-0.000117594,5.83098e-08,-1.08293e-11,-35134.4,40.3319], Tmin=(100,'K'), Tmax=(1495.46,'K')), NASAPolynomial(coeffs=[30.0167,0.0122621,-1.49485e-06,3.08492e-11,3.72211e-15,-43183.9,-123.392], Tmin=(1495.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-294.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CsJOCH3)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)O[CH2](5959)',
    structure = SMILES('[CH]=C(O)C(O)(CC)O[CH2]'),
    E0 = (-184.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1600,1603.76,2882.23,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15088,'amu*angstrom^2'), symmetry=1, barrier=(3.46903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.95809,0.122919,-0.000133823,6.79397e-08,-1.27617e-11,-21964.4,41.4145], Tmin=(100,'K'), Tmax=(1502.01,'K')), NASAPolynomial(coeffs=[33.4589,0.00680793,1.23591e-06,-4.85389e-10,3.84882e-14,-30746.3,-141.908], Tmin=(1502.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CsJOCH3) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C(O)(OC)C(=C)O(32219)',
    structure = SMILES('[CH2][CH]C(O)(OC)C(=C)O'),
    E0 = (-214.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,445.742,682.215,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155755,'amu*angstrom^2'), symmetry=1, barrier=(3.58111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.70995,0.115193,-0.000119389,5.8189e-08,-1.05024e-11,-25557.7,43.6251], Tmin=(100,'K'), Tmax=(1574.78,'K')), NASAPolynomial(coeffs=[31.4882,0.00892728,3.09666e-07,-3.07507e-10,2.59573e-14,-33923,-129.253], Tmin=(1574.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'C=C([O])C(O)([CH]C)OC(32220)',
    structure = SMILES('C=C([O])C(O)([CH]C)OC'),
    E0 = (-282.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.156,0.107875,-0.000107496,5.16014e-08,-9.29077e-12,-33693.2,41.3486], Tmin=(100,'K'), Tmax=(1550.72,'K')), NASAPolynomial(coeffs=[28.6301,0.0136857,-2.09157e-06,1.41545e-10,-3.86514e-15,-41464.4,-114.953], Tmin=(1550.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C(O)C(O)([CH]C)OC(32221)',
    structure = SMILES('[CH]=C(O)C(O)([CH]C)OC'),
    E0 = (-172.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180,180,180,1600,1634.41,2853.66,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151895,'amu*angstrom^2'), symmetry=1, barrier=(3.49236,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.76119,0.117345,-0.000123653,6.11639e-08,-1.1203e-11,-20523.4,42.4209], Tmin=(100,'K'), Tmax=(1547.54,'K')), NASAPolynomial(coeffs=[31.9829,0.00836008,5.73778e-07,-3.60648e-10,2.98192e-14,-28980.3,-132.95], Tmin=(1547.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = '[CH2][O](167)',
    structure = SMILES('[CH2][O]'),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88409,-0.00363885,3.28543e-05,-4.13611e-08,1.59631e-11,23210.8,7.47983], Tmin=(100,'K'), Tmax=(933.06,'K')), NASAPolynomial(coeffs=[6.69335,0.000289989,8.61416e-07,-1.56351e-10,7.33778e-15,21991.3,-9.6043], Tmin=(933.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(H3COJ) + radical(CsJOH)"""),
)

species(
    label = '[CH2]OC1(O)[C](O)CC1C(32222)',
    structure = SMILES('[CH2]OC1(O)[C](O)CC1C'),
    E0 = (-190.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.870401,0.0842029,-2.21366e-05,-5.15043e-08,3.13989e-11,-22753,34.0202], Tmin=(100,'K'), Tmax=(936.744,'K')), NASAPolynomial(coeffs=[29.1036,0.014044,-2.39923e-06,3.56359e-10,-3.12861e-14,-30905.9,-122.172], Tmin=(936.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclobutane) + radical(CsJOCH3) + radical(C2CsJOH)"""),
)

species(
    label = 'C[CH]C1(O)OCC[C]1O(32223)',
    structure = SMILES('C[CH]C1(O)OCC[C]1O'),
    E0 = (-266.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.33742,0.0704856,1.80145e-05,-9.99916e-08,5.2698e-11,-31917.6,31.3417], Tmin=(100,'K'), Tmax=(887.838,'K')), NASAPolynomial(coeffs=[28.506,0.00985057,3.35193e-06,-1.04877e-09,7.70641e-14,-39771.1,-119.764], Tmin=(887.838,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-266.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(O)(OC)C(=C)O(32224)',
    structure = SMILES('C=CC(O)(OC)C(=C)O'),
    E0 = (-493.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.57665,0.114181,-0.000117374,5.7053e-08,-1.03176e-11,-59060.1,38.6687], Tmin=(100,'K'), Tmax=(1557.79,'K')), NASAPolynomial(coeffs=[31.0961,0.0102753,-5.27233e-07,-1.40018e-10,1.44587e-14,-67434.7,-131.916], Tmin=(1557.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]O[C](O)C(C)C(=C)O(32225)',
    structure = SMILES('[CH2]O[C](O)C(C)C(=C)O'),
    E0 = (-209.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
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
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09541,0.103417,-0.000105929,5.47381e-08,-1.10809e-11,-24949.4,39.1028], Tmin=(100,'K'), Tmax=(1208.14,'K')), NASAPolynomial(coeffs=[22.2536,0.0261111,-9.94649e-06,1.77361e-09,-1.20864e-13,-30591.2,-77.9574], Tmin=(1208.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_P) + radical(CsJOCH3)"""),
)

species(
    label = 'C=C(O)C1(O)OCC1C(31894)',
    structure = SMILES('C=C(O)C1(O)OCC1C'),
    E0 = (-489.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.702902,0.0728009,3.06256e-05,-1.27487e-07,6.63106e-11,-58665.5,29.4672], Tmin=(100,'K'), Tmax=(885.342,'K')), NASAPolynomial(coeffs=[34.2423,0.00030067,8.79859e-06,-2.12059e-09,1.50508e-13,-68199.5,-153.766], Tmin=(885.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-489.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Oxetane)"""),
)

species(
    label = '[CH2]OC(O)([CH]C)C(C)=O(32226)',
    structure = SMILES('[CH2]OC(O)([CH]C)C(C)=O'),
    E0 = (-256.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.8869,0.0905553,-5.43013e-05,-8.98979e-09,1.39548e-11,-30634.8,38.4561], Tmin=(100,'K'), Tmax=(957.439,'K')), NASAPolynomial(coeffs=[25.6781,0.0196327,-5.95102e-06,1.04554e-09,-7.66843e-14,-37557.8,-98.1382], Tmin=(957.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-256.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CsJOCH3) + radical(CCJCO)"""),
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
    label = 'C=C(O)C([O])(O)[CH]C(31687)',
    structure = SMILES('C=C(O)C([O])(O)[CH]C'),
    E0 = (-207.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,440,435,1725,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180,1600,1767.23,2733.8,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156597,'amu*angstrom^2'), symmetry=1, barrier=(3.60046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156597,'amu*angstrom^2'), symmetry=1, barrier=(3.60046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156597,'amu*angstrom^2'), symmetry=1, barrier=(3.60046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156597,'amu*angstrom^2'), symmetry=1, barrier=(3.60046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156597,'amu*angstrom^2'), symmetry=1, barrier=(3.60046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58664,0.0981833,-0.000107415,5.53071e-08,-1.04939e-11,-24756.3,37.3299], Tmin=(100,'K'), Tmax=(1507.57,'K')), NASAPolynomial(coeffs=[26.5121,0.00615983,1.52859e-06,-5.55416e-10,4.43984e-14,-31243.3,-103.181], Tmin=(1507.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)OJ) + radical(CCJCO)"""),
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
    label = '[CH2]O[C](O)C(=C)O(32227)',
    structure = SMILES('[CH2]OC(O)=C([CH2])O'),
    E0 = (-204.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,364.81,364.81],'cm^-1')),
        HinderedRotor(inertia=(0.16712,'amu*angstrom^2'), symmetry=1, barrier=(15.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167121,'amu*angstrom^2'), symmetry=1, barrier=(15.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167121,'amu*angstrom^2'), symmetry=1, barrier=(15.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167121,'amu*angstrom^2'), symmetry=1, barrier=(15.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167121,'amu*angstrom^2'), symmetry=1, barrier=(15.783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.858891,0.0882372,-0.000105136,5.74744e-08,-1.15074e-11,-24435.5,29.1136], Tmin=(100,'K'), Tmax=(1424.53,'K')), NASAPolynomial(coeffs=[24.6873,0.000998281,3.05272e-06,-7.98537e-10,6.03452e-14,-30140.3,-97.6487], Tmin=(1424.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(C=C(O)CJ) + radical(C=COCJ)"""),
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
    E0 = (-232.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-168.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-169.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-90.3966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-213.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-38.0267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-118.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-114.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-126.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-89.8699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-148.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-140.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-160.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-199.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-139.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (70.0577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-107.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-117.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-207.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-49.0927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-223.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-27.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (173.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (139.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['CH2O(15)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]OC1(O)C(C)C1([CH2])O(32214)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.39513e+10,'s^-1'), n=0.560608, Ea=(63.126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd] + [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 63.1 to 63.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]C1(O)COC1(O)[CH]C(32215)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.48e+06,'s^-1'), n=1.25, Ea=(62.76,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 349 used for R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2]OC(O)(C=C)C(=C)O(32216)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.014e+08,'cm^3/(mol*s)'), n=1.733, Ea=(3.17984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2823 used for Cds-HH_Cds-Cs\O2s/H;HJ
Exact match found for rate rule [Cds-HH_Cds-Cs\O2s/H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2O(15)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2330,'cm^3/(mol*s)'), n=3.17, Ea=(28.0328,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2COH(99)', '[CH2]OC(O)=CC(32217)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsOs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(5)', '[CH2]OC(=CC)C(=C)O(32218)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.812049,'m^3/(mol*s)'), n=2.01336, Ea=(12.3541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;YJ] for rate rule [Cds-CdOs_Cds;OJ_pri]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]CC(O)(O[CH2])C(=C)O(5957)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]OC([O])(CC)C(=C)O(5515)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C(O)C([O])([CH]C)OC(5956)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.3e+09,'s^-1','*|/',2.51189), n=0.75, Ea=(97.0688,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R4H_SSS;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]OC(O)(CC)C(=C)[O](5958)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(O)C(O)(CC)O[CH2](5959)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C(O)(OC)C(=C)O(32219)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['C=C([O])C(O)([CH]C)OC(32220)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.06724e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(O)C(O)([CH]C)OC(32221)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][O](167)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]OC1(O)[C](O)CC1C(32222)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.71e+08,'s^-1'), n=0.99, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['C[CH]C1(O)OCC[C]1O(32223)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.71e+11,'s^-1'), n=0.2, Ea=(115.06,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H
Exact match found for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['C=CC(O)(OC)C(=C)O(32224)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]O[C](O)C(C)C(=C)O(32225)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['C=C(O)C1(O)OCC1C(31894)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    products = ['[CH2]OC(O)([CH]C)C(C)=O(32226)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_Cs;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(19)', 'C=C(O)C([O])(O)[CH]C(31687)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CHCH3(T)(95)', '[CH2]O[C](O)C(=C)O(32227)'],
    products = ['[CH2]OC(O)([CH]C)C(=C)O(5955)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5233',
    isomers = [
        '[CH2]OC(O)([CH]C)C(=C)O(5955)',
    ],
    reactants = [
        ('CH2O(15)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5233',
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

