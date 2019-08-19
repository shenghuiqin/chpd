species(
    label = 'S(161)(160)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,375,552.5,462.5,1710,332.668,332.668],'cm^-1')),
        HinderedRotor(inertia=(0.152191,'amu*angstrom^2'), symmetry=1, barrier=(11.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152191,'amu*angstrom^2'), symmetry=1, barrier=(11.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152191,'amu*angstrom^2'), symmetry=1, barrier=(11.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152191,'amu*angstrom^2'), symmetry=1, barrier=(11.9519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'S(168)(167)',
    structure = SMILES('C=C(O)C(O)=CC'),
    E0 = (-369.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(0.87423,'amu*angstrom^2'), symmetry=1, barrier=(20.1003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.87642,'amu*angstrom^2'), symmetry=1, barrier=(20.1506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.878567,'amu*angstrom^2'), symmetry=1, barrier=(20.2,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875855,'amu*angstrom^2'), symmetry=1, barrier=(20.1376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884099,0.0877395,-9.34163e-05,4.76671e-08,-9.0713e-12,-44294.8,25.8217], Tmin=(100,'K'), Tmax=(1473.72,'K')), NASAPolynomial(coeffs=[23.5689,0.00887539,-4.29798e-07,-1.49469e-10,1.60597e-14,-50145.5,-97.0299], Tmin=(1473.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = 'S(169)(168)',
    structure = SMILES('CC1CC(O)=C1O'),
    E0 = (-316.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4460.07,'J/mol'), sigma=(7.06575,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=696.65 K, Pc=28.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716318,0.0465129,4.18254e-05,-1.08121e-07,5.08748e-11,-37912.2,22.5083], Tmin=(100,'K'), Tmax=(918.644,'K')), NASAPolynomial(coeffs=[26.071,0.00311907,3.26936e-06,-7.39926e-10,4.4025e-14,-45398,-113.051], Tmin=(918.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene)"""),
)

species(
    label = 'S(172)(171)',
    structure = SMILES('C=CC(O)=C(C)O'),
    E0 = (-370.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(0.869774,'amu*angstrom^2'), symmetry=1, barrier=(19.9978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86994,'amu*angstrom^2'), symmetry=1, barrier=(20.0016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869436,'amu*angstrom^2'), symmetry=1, barrier=(19.99,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869771,'amu*angstrom^2'), symmetry=1, barrier=(19.9978,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1976,0.0877436,-9.10797e-05,4.46402e-08,-8.05724e-12,-44357.6,27.774], Tmin=(100,'K'), Tmax=(1598.72,'K')), NASAPolynomial(coeffs=[24.3035,0.00683599,8.79619e-07,-3.98758e-10,3.21892e-14,-50325.6,-100.383], Tmin=(1598.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'S(173)(172)',
    structure = SMILES('C=CC(O)C(C)=O'),
    E0 = (-325.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,347.492,348.832],'cm^-1')),
        HinderedRotor(inertia=(0.00139452,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00140098,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128246,'amu*angstrom^2'), symmetry=1, barrier=(10.9614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126695,'amu*angstrom^2'), symmetry=1, barrier=(10.9669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4073.27,'J/mol'), sigma=(6.50815,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.24 K, Pc=33.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0618,0.0556095,-2.99883e-05,1.34976e-09,2.75317e-12,-39000.2,28.3987], Tmin=(100,'K'), Tmax=(1110.77,'K')), NASAPolynomial(coeffs=[13.3828,0.0262972,-1.07381e-05,2.00002e-09,-1.39935e-13,-42666.2,-36.5188], Tmin=(1110.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-325.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'S(175)(174)',
    structure = SMILES('CC1(O)CC=C1O'),
    E0 = (-291.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4466.08,'J/mol'), sigma=(7.21957,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=697.59 K, Pc=26.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8346,0.0482796,2.58179e-05,-8.3966e-08,4.04937e-11,-34951.8,23.8218], Tmin=(100,'K'), Tmax=(920.94,'K')), NASAPolynomial(coeffs=[22.7811,0.008409,4.40278e-07,-2.14397e-10,9.95424e-15,-41345.6,-93.0167], Tmin=(920.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-291.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'S(176)(175)',
    structure = SMILES('C=CC(=O)C(C)O'),
    E0 = (-321.792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,180,791.468],'cm^-1')),
        HinderedRotor(inertia=(0.00890222,'amu*angstrom^2'), symmetry=1, barrier=(3.95866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.025921,'amu*angstrom^2'), symmetry=1, barrier=(11.5338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.501926,'amu*angstrom^2'), symmetry=1, barrier=(11.5403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502049,'amu*angstrom^2'), symmetry=1, barrier=(11.5431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4073.27,'J/mol'), sigma=(6.50815,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.24 K, Pc=33.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02967,0.0574933,-3.81026e-05,1.21107e-08,-1.53223e-12,-38589.2,28.2443], Tmin=(100,'K'), Tmax=(1830.19,'K')), NASAPolynomial(coeffs=[16.9171,0.02277,-9.64365e-06,1.74414e-09,-1.16172e-13,-44404.6,-58.0061], Tmin=(1830.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-321.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C3H4O(121)(120)',
    structure = SMILES('CC=C=O'),
    E0 = (-77.7156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.447917,'amu*angstrom^2'), symmetry=1, barrier=(10.2985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2689.99,'J/mol'), sigma=(5.15244,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=420.17 K, Pc=44.62 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00496,0.0183622,-5.86278e-07,-8.1995e-09,3.3864e-12,-9308.28,12.0089], Tmin=(100,'K'), Tmax=(1106.37,'K')), NASAPolynomial(coeffs=[5.96289,0.0150316,-6.0541e-06,1.11098e-09,-7.67696e-14,-10413.5,-4.59703], Tmin=(1106.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.7156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""methylketene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H4O(40)(40)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C=C=C([O])CC(8829)',
    structure = SMILES('C=C=C([O])CC'),
    E0 = (47.3032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725,220.708,222.018],'cm^-1')),
        HinderedRotor(inertia=(0.34305,'amu*angstrom^2'), symmetry=1, barrier=(11.7864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339408,'amu*angstrom^2'), symmetry=1, barrier=(11.791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36132,0.0508126,-3.08268e-05,2.54214e-09,3.1323e-12,5790.59,21.9667], Tmin=(100,'K'), Tmax=(1007,'K')), NASAPolynomial(coeffs=[12.5297,0.020834,-7.59769e-06,1.34845e-09,-9.28702e-14,2811.96,-35.6132], Tmin=(1007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.3032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
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
    label = '[CH2]C(=O)C(=O)CC(8842)',
    structure = SMILES('[CH2]C(=O)C(=O)CC'),
    E0 = (-162.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,365,385,505,600,445,480,1700,1720,352.088,2238.1],'cm^-1')),
        HinderedRotor(inertia=(0.0859868,'amu*angstrom^2'), symmetry=1, barrier=(7.56458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859642,'amu*angstrom^2'), symmetry=1, barrier=(7.56459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859615,'amu*angstrom^2'), symmetry=1, barrier=(7.56392,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37635,'amu*angstrom^2'), symmetry=1, barrier=(33.1076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15233,0.0680038,-7.99912e-05,5.95539e-08,-1.92307e-11,-19492.6,25.214], Tmin=(100,'K'), Tmax=(738.439,'K')), NASAPolynomial(coeffs=[6.798,0.0374228,-1.78728e-05,3.47419e-09,-2.45165e-13,-20326.4,-0.311287], Tmin=(738.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-162.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(CJCC=O)"""),
)

species(
    label = 'CH3(15)(16)',
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
    label = 'C=C([O])C(=C)O(8849)',
    structure = SMILES('C=C([O])C(=C)O'),
    E0 = (-196.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.916777,'amu*angstrom^2'), symmetry=1, barrier=(21.0785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920477,'amu*angstrom^2'), symmetry=1, barrier=(21.1636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284274,0.0664439,-7.46756e-05,3.95333e-08,-7.70645e-12,-23433.6,21.3173], Tmin=(100,'K'), Tmax=(1467.22,'K')), NASAPolynomial(coeffs=[18.8952,0.00393949,1.25465e-06,-4.33245e-10,3.4763e-14,-27628.4,-71.2884], Tmin=(1467.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C2H5(27)(28)',
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
    label = '[CH2]C(O)=C=O(8828)',
    structure = SMILES('C=C(O)[C]=O'),
    E0 = (-127.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2950,3100,1380,975,1025,1650,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.03693,'amu*angstrom^2'), symmetry=1, barrier=(23.8411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03759,'amu*angstrom^2'), symmetry=1, barrier=(23.8562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43271,0.0523177,-6.55172e-05,3.78375e-08,-8.27375e-12,-15211.5,15.8439], Tmin=(100,'K'), Tmax=(1138,'K')), NASAPolynomial(coeffs=[15.2759,0.00365942,-1.38033e-06,2.64519e-10,-1.95358e-14,-18362.2,-52.7311], Tmin=(1138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=C(O)C([O])=CC(8848)',
    structure = SMILES('C=C(O)C([O])=CC'),
    E0 = (-232.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.757919,'amu*angstrom^2'), symmetry=1, barrier=(17.426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756197,'amu*angstrom^2'), symmetry=1, barrier=(17.3865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757974,'amu*angstrom^2'), symmetry=1, barrier=(17.4273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145471,0.0794715,-8.63325e-05,4.63576e-08,-9.46463e-12,-27754.2,24.6653], Tmin=(100,'K'), Tmax=(1313.65,'K')), NASAPolynomial(coeffs=[19.8698,0.0123363,-2.60601e-06,2.80197e-10,-1.3067e-14,-32478.7,-75.3243], Tmin=(1313.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(=C)O(8855)',
    structure = SMILES('[CH2]CC(=O)C(=C)O'),
    E0 = (-146.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,382.141,382.144],'cm^-1')),
        HinderedRotor(inertia=(0.116187,'amu*angstrom^2'), symmetry=1, barrier=(12.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116191,'amu*angstrom^2'), symmetry=1, barrier=(12.0406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116189,'amu*angstrom^2'), symmetry=1, barrier=(12.0407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116187,'amu*angstrom^2'), symmetry=1, barrier=(12.0406,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661955,0.0750645,-8.3545e-05,4.91082e-08,-1.16253e-11,-17546.5,24.0304], Tmin=(100,'K'), Tmax=(1021.82,'K')), NASAPolynomial(coeffs=[13.1639,0.0261243,-1.1702e-05,2.23537e-09,-1.57235e-13,-20101.5,-36.554], Tmin=(1021.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(CJCC=O)"""),
)

species(
    label = 'C3H5O(279)(278)',
    structure = SMILES('CC[C]=O'),
    E0 = (-44.1874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.192008,'amu*angstrom^2'), symmetry=1, barrier=(4.41465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191563,'amu*angstrom^2'), symmetry=1, barrier=(4.40441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2COH(1431)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.989115,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C(O)C(=O)CC(8856)',
    structure = SMILES('[CH]=C(O)C(=O)CC'),
    E0 = (-111.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,252.473],'cm^-1')),
        HinderedRotor(inertia=(0.270183,'amu*angstrom^2'), symmetry=1, barrier=(12.2712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270436,'amu*angstrom^2'), symmetry=1, barrier=(12.2638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26826,'amu*angstrom^2'), symmetry=1, barrier=(12.2637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268894,'amu*angstrom^2'), symmetry=1, barrier=(12.2737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7613,0.0731453,-7.95332e-05,4.60811e-08,-1.08237e-11,-13279.8,23.4437], Tmin=(100,'K'), Tmax=(1025.9,'K')), NASAPolynomial(coeffs=[12.568,0.027111,-1.22252e-05,2.34201e-09,-1.64964e-13,-15702.3,-33.8185], Tmin=(1025.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C(=O)CC(8843)',
    structure = SMILES('[CH2]C([O])C(=O)CC'),
    E0 = (1.82978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,253.256,253.258,4000],'cm^-1')),
        HinderedRotor(inertia=(0.550332,'amu*angstrom^2'), symmetry=1, barrier=(25.0453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000991364,'amu*angstrom^2'), symmetry=1, barrier=(11.256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0112888,'amu*angstrom^2'), symmetry=1, barrier=(8.80378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550306,'amu*angstrom^2'), symmetry=1, barrier=(25.0455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15704,0.0648368,-5.84199e-05,2.96528e-08,-6.3845e-12,320.515,28.752], Tmin=(100,'K'), Tmax=(1085.18,'K')), NASAPolynomial(coeffs=[9.74089,0.0331967,-1.46854e-05,2.78524e-09,-1.94876e-13,-1542.5,-13.3619], Tmin=(1085.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.82978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCO) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C(O)C([O])[CH]C(8850)',
    structure = SMILES('C=C(O)C([O])[CH]C'),
    E0 = (2.59554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,334.186,336.32,336.519],'cm^-1')),
        HinderedRotor(inertia=(0.0014928,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00151766,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216028,'amu*angstrom^2'), symmetry=1, barrier=(16.9563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211367,'amu*angstrom^2'), symmetry=1, barrier=(16.9167,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582894,0.0601292,-1.74707e-05,-3.07624e-08,1.89857e-11,448.888,30.2159], Tmin=(100,'K'), Tmax=(949.697,'K')), NASAPolynomial(coeffs=[20.0817,0.0148432,-4.13141e-06,7.20283e-10,-5.43182e-14,-4916.08,-71.5952], Tmin=(949.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.59554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC([O])C(=C)O(8851)',
    structure = SMILES('[CH2]CC([O])C(=C)O'),
    E0 = (7.93984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,390.231,390.372,4000],'cm^-1')),
        HinderedRotor(inertia=(0.174374,'amu*angstrom^2'), symmetry=1, barrier=(18.8524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174465,'amu*angstrom^2'), symmetry=1, barrier=(18.8534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174604,'amu*angstrom^2'), symmetry=1, barrier=(18.8551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174453,'amu*angstrom^2'), symmetry=1, barrier=(18.852,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.425828,0.066073,-3.65313e-05,-9.05001e-09,1.07201e-11,1094.94,30.0697], Tmin=(100,'K'), Tmax=(963.282,'K')), NASAPolynomial(coeffs=[19.3512,0.0167769,-5.37961e-06,9.56946e-10,-6.94138e-14,-3910.11,-67.5799], Tmin=(963.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.93984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(O)C(=O)[CH]C(8844)',
    structure = SMILES('[CH2]C(O)C([O])=CC'),
    E0 = (-90.5188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,299.833,312.963],'cm^-1')),
        HinderedRotor(inertia=(0.00185591,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123576,'amu*angstrom^2'), symmetry=1, barrier=(8.80827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130752,'amu*angstrom^2'), symmetry=1, barrier=(8.86809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12886,'amu*angstrom^2'), symmetry=1, barrier=(8.83192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345407,0.0723962,-6.93553e-05,3.42305e-08,-6.65407e-12,-10748.6,30.4624], Tmin=(100,'K'), Tmax=(1255.15,'K')), NASAPolynomial(coeffs=[16.6574,0.0204121,-7.2304e-06,1.23314e-09,-8.16962e-14,-14843.4,-51.9406], Tmin=(1255.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.5188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(=O)C([O])CC(8770)',
    structure = SMILES('C=C([O])C([O])CC'),
    E0 = (-59.5017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,411.226,412.428,412.619,412.693],'cm^-1')),
        HinderedRotor(inertia=(0.000990404,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0846735,'amu*angstrom^2'), symmetry=1, barrier=(10.0852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0848459,'amu*angstrom^2'), symmetry=1, barrier=(10.0874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728515,0.0616743,-3.46694e-05,-2.82684e-09,6.53814e-12,-7029.53,28.6979], Tmin=(100,'K'), Tmax=(994.485,'K')), NASAPolynomial(coeffs=[16.0073,0.0221145,-8.02439e-06,1.44923e-09,-1.01988e-13,-11151.1,-50.3725], Tmin=(994.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.5017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C([O])CC(8852)',
    structure = SMILES('[CH]=C(O)C([O])CC'),
    E0 = (49.7896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,366.46,366.615],'cm^-1')),
        HinderedRotor(inertia=(0.148496,'amu*angstrom^2'), symmetry=1, barrier=(14.1347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14856,'amu*angstrom^2'), symmetry=1, barrier=(14.1332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148272,'amu*angstrom^2'), symmetry=1, barrier=(14.1345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147862,'amu*angstrom^2'), symmetry=1, barrier=(14.1341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.312838,0.068886,-4.28217e-05,-3.80964e-09,9.19362e-12,6132.1,29.0919], Tmin=(100,'K'), Tmax=(960.921,'K')), NASAPolynomial(coeffs=[19.819,0.0162025,-5.09347e-06,8.95974e-10,-6.47584e-14,1066.86,-71.0865], Tmin=(960.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.7896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC(O)=C([CH2])O(8832)',
    structure = SMILES('[CH2]CC(O)=C([CH2])O'),
    E0 = (-117.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,424.662],'cm^-1')),
        HinderedRotor(inertia=(0.112579,'amu*angstrom^2'), symmetry=1, barrier=(14.3919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112434,'amu*angstrom^2'), symmetry=1, barrier=(14.3917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112485,'amu*angstrom^2'), symmetry=1, barrier=(14.3926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112501,'amu*angstrom^2'), symmetry=1, barrier=(14.3922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11241,'amu*angstrom^2'), symmetry=1, barrier=(14.3919,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06791,0.0870529,-9.14386e-05,4.55824e-08,-8.38467e-12,-13928.7,32.6758], Tmin=(100,'K'), Tmax=(1562.56,'K')), NASAPolynomial(coeffs=[23.7823,0.00723112,7.46054e-07,-3.86281e-10,3.20893e-14,-19716.1,-91.9721], Tmin=(1562.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(RCCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C[CH]C([O])=C(C)O(8833)',
    structure = SMILES('C[CH]C([O])=C(C)O'),
    E0 = (-143.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,363.148,363.217],'cm^-1')),
        HinderedRotor(inertia=(0.0012773,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111173,'amu*angstrom^2'), symmetry=1, barrier=(10.4121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111176,'amu*angstrom^2'), symmetry=1, barrier=(10.4118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111176,'amu*angstrom^2'), symmetry=1, barrier=(10.4116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539627,0.0644069,-3.47983e-05,-1.16762e-08,1.24888e-11,-17178.7,28.4164], Tmin=(100,'K'), Tmax=(930.903,'K')), NASAPolynomial(coeffs=[18.8584,0.0154804,-3.95891e-06,6.11722e-10,-4.24718e-14,-21880,-65.5818], Tmin=(930.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C([O])=C(O)CC(8834)',
    structure = SMILES('[CH2]C([O])=C(O)CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,403.683,407.159],'cm^-1')),
        HinderedRotor(inertia=(0.0985501,'amu*angstrom^2'), symmetry=1, barrier=(11.483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101351,'amu*angstrom^2'), symmetry=1, barrier=(11.5221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101048,'amu*angstrom^2'), symmetry=1, barrier=(11.521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10204,'amu*angstrom^2'), symmetry=1, barrier=(11.505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(O)[C](O)CC(8853)',
    structure = SMILES('[CH]C(O)=C(O)CC'),
    E0 = (-110.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00203323,0.0732248,-3.43625e-05,-2.16115e-08,1.78208e-11,-13190.3,27.458], Tmin=(100,'K'), Tmax=(926.829,'K')), NASAPolynomial(coeffs=[21.8543,0.0167143,-4.08005e-06,6.09833e-10,-4.25184e-14,-18864.5,-85.0647], Tmin=(926.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])O(8845)',
    structure = SMILES('[CH2]CC(=O)C([CH2])O'),
    E0 = (-30.3142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,340.066,342.119],'cm^-1')),
        HinderedRotor(inertia=(0.00147515,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144597,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00144665,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115757,'amu*angstrom^2'), symmetry=1, barrier=(9.4902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113377,'amu*angstrom^2'), symmetry=1, barrier=(9.49143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735022,0.0732832,-7.55104e-05,4.2154e-08,-9.62104e-12,-3529.57,30.1206], Tmin=(100,'K'), Tmax=(1050.32,'K')), NASAPolynomial(coeffs=[12.2345,0.029489,-1.29662e-05,2.45533e-09,-1.71838e-13,-5945.18,-25.9221], Tmin=(1050.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.3142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(CJCO)"""),
)

species(
    label = '[CH2]CC([O])=C(C)O(8835)',
    structure = SMILES('[CH2]CC([O])=C(C)O'),
    E0 = (-138.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,318.087,318.17],'cm^-1')),
        HinderedRotor(inertia=(0.0016365,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143797,'amu*angstrom^2'), symmetry=1, barrier=(10.5463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141019,'amu*angstrom^2'), symmetry=1, barrier=(10.5785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137438,'amu*angstrom^2'), symmetry=1, barrier=(10.5767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.149023,0.0765103,-7.49082e-05,3.66676e-08,-6.8673e-12,-16509.4,30.1846], Tmin=(100,'K'), Tmax=(1435.78,'K')), NASAPolynomial(coeffs=[19.7458,0.0146968,-3.65663e-06,4.85292e-10,-2.76616e-14,-21563.9,-70.6997], Tmin=(1435.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CH2(S)(21)(22)',
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
    label = 'C=C(O)C(C)=O(8857)',
    structure = SMILES('C=C(O)C(C)=O'),
    E0 = (-333.051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,350,440,435,1725,336.509],'cm^-1')),
        HinderedRotor(inertia=(0.168454,'amu*angstrom^2'), symmetry=1, barrier=(13.5274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168467,'amu*angstrom^2'), symmetry=1, barrier=(13.5272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168491,'amu*angstrom^2'), symmetry=1, barrier=(13.5287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34048,0.0507425,-3.55578e-05,7.13228e-09,1.49089e-12,-39954.3,19.6241], Tmin=(100,'K'), Tmax=(1048.53,'K')), NASAPolynomial(coeffs=[13.8723,0.0167056,-6.56509e-06,1.2237e-09,-8.67531e-14,-43339.3,-45.0381], Tmin=(1048.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]C(O)C(=O)CC(8861)',
    structure = SMILES('[CH]C(O)C(=O)CC'),
    E0 = (-3.19103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725856,0.0785881,-0.000108626,9.24897e-08,-3.25033e-11,-272.212,27.444], Tmin=(100,'K'), Tmax=(789.086,'K')), NASAPolynomial(coeffs=[6.89798,0.0392746,-1.86372e-05,3.57145e-09,-2.48192e-13,-996.407,0.71246], Tmin=(789.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-3.19103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CCC(=O)C(C)=O(8858)',
    structure = SMILES('CCC(=O)C(C)=O'),
    E0 = (-374.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51646,0.0601905,-5.04133e-05,2.65316e-08,-6.55815e-12,-44953.8,22.5887], Tmin=(100,'K'), Tmax=(902.454,'K')), NASAPolynomial(coeffs=[5.68542,0.0417124,-1.97003e-05,3.84325e-09,-2.73025e-13,-45706.3,2.90371], Tmin=(902.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-374.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs)"""),
)

species(
    label = 'C=C(O)C(=C)OC(8859)',
    structure = SMILES('C=C(O)C(=C)OC'),
    E0 = (-315.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.852525,'amu*angstrom^2'), symmetry=1, barrier=(19.6012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854974,'amu*angstrom^2'), symmetry=1, barrier=(19.6575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853489,'amu*angstrom^2'), symmetry=1, barrier=(19.6234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.851534,'amu*angstrom^2'), symmetry=1, barrier=(19.5784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.879363,0.085231,-8.65755e-05,4.21603e-08,-7.66487e-12,-37712.8,25.643], Tmin=(100,'K'), Tmax=(1548.07,'K')), NASAPolynomial(coeffs=[23.4666,0.00958852,-9.41414e-07,-3.17674e-11,7.00573e-15,-43724.6,-97.5224], Tmin=(1548.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(CC)OO(8860)',
    structure = SMILES('C=C=C(CC)OO'),
    E0 = (33.9932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,301.12],'cm^-1')),
        HinderedRotor(inertia=(0.232477,'amu*angstrom^2'), symmetry=1, barrier=(5.34511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691181,'amu*angstrom^2'), symmetry=1, barrier=(15.8916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691199,'amu*angstrom^2'), symmetry=1, barrier=(15.892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691106,'amu*angstrom^2'), symmetry=1, barrier=(15.8899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09237,0.0654842,-5.60116e-05,2.59736e-08,-5.06682e-12,4191.76,25.4768], Tmin=(100,'K'), Tmax=(1189.76,'K')), NASAPolynomial(coeffs=[10.852,0.0326716,-1.46425e-05,2.79265e-09,-1.95839e-13,1869.45,-23.3036], Tmin=(1189.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.9932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=[C]C(O)=CC(8862)',
    structure = SMILES('C=[C]C(O)=CC'),
    E0 = (42.9982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02074,'amu*angstrom^2'), symmetry=1, barrier=(23.4687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02041,'amu*angstrom^2'), symmetry=1, barrier=(23.4612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01954,'amu*angstrom^2'), symmetry=1, barrier=(23.4413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502882,0.0666463,-6.75e-05,3.49508e-08,-6.96103e-12,5306.16,21.1922], Tmin=(100,'K'), Tmax=(1342.26,'K')), NASAPolynomial(coeffs=[16.3641,0.0144229,-3.60075e-06,4.62785e-10,-2.5161e-14,1494.65,-58.3348], Tmin=(1342.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C([O])C(O)=CC(8863)',
    structure = SMILES('C=C([O])C(O)=CC'),
    E0 = (-232.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.75831,'amu*angstrom^2'), symmetry=1, barrier=(17.435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755543,'amu*angstrom^2'), symmetry=1, barrier=(17.3714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.758681,'amu*angstrom^2'), symmetry=1, barrier=(17.4436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145462,0.0794714,-8.63323e-05,4.63573e-08,-9.46453e-12,-27754.2,24.6653], Tmin=(100,'K'), Tmax=(1313.67,'K')), NASAPolynomial(coeffs=[19.8697,0.0123364,-2.60611e-06,2.80219e-10,-1.30687e-14,-32478.6,-75.3236], Tmin=(1313.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)[C]=CC(8864)',
    structure = SMILES('C=C(O)[C]=CC'),
    E0 = (42.9982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02074,'amu*angstrom^2'), symmetry=1, barrier=(23.4687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02041,'amu*angstrom^2'), symmetry=1, barrier=(23.4612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01954,'amu*angstrom^2'), symmetry=1, barrier=(23.4413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502882,0.0666463,-6.75e-05,3.49508e-08,-6.96103e-12,5306.16,21.1922], Tmin=(100,'K'), Tmax=(1342.26,'K')), NASAPolynomial(coeffs=[16.3641,0.0144229,-3.60075e-06,4.62785e-10,-2.5161e-14,1494.65,-58.3348], Tmin=(1342.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C(O)C(=C)O(8865)',
    structure = SMILES('[CH]=C(O)C(=C)O'),
    E0 = (-86.7684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,325,375,415,465,420,450,1700,1750],'cm^-1')),
        HinderedRotor(inertia=(1.02704,'amu*angstrom^2'), symmetry=1, barrier=(23.6138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0269,'amu*angstrom^2'), symmetry=1, barrier=(23.6104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02731,'amu*angstrom^2'), symmetry=1, barrier=(23.6199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.322663,0.0759296,-9.08749e-05,4.91324e-08,-9.62861e-12,-10263.7,22.3962], Tmin=(100,'K'), Tmax=(1482.05,'K')), NASAPolynomial(coeffs=[22.3441,-0.00152328,3.98935e-06,-9.50257e-10,6.95844e-14,-15194.9,-89.8453], Tmin=(1482.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.7684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C(O)C(=C)O(8866)',
    structure = SMILES('[CH2]C=C(O)C(=C)O'),
    E0 = (-251.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(1.12227,'amu*angstrom^2'), symmetry=1, barrier=(25.8033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1216,'amu*angstrom^2'), symmetry=1, barrier=(25.7878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12387,'amu*angstrom^2'), symmetry=1, barrier=(25.8401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12013,'amu*angstrom^2'), symmetry=1, barrier=(25.754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04714,0.0838822,-8.70686e-05,4.2758e-08,-7.70112e-12,-30083.6,27.9774], Tmin=(100,'K'), Tmax=(1619.92,'K')), NASAPolynomial(coeffs=[22.8798,0.00666956,1.21675e-06,-4.84385e-10,3.87218e-14,-35456.7,-91.6557], Tmin=(1619.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'CC=[C]O(8867)',
    structure = SMILES('CC=[C]O'),
    E0 = (72.7982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.344811,'amu*angstrom^2'), symmetry=1, barrier=(7.92788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346062,'amu*angstrom^2'), symmetry=1, barrier=(7.95665,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.692,0.0334424,-4.69125e-05,4.57787e-08,-1.77014e-11,8798.16,15.2755], Tmin=(100,'K'), Tmax=(837.57,'K')), NASAPolynomial(coeffs=[2.44779,0.0239997,-1.10022e-05,2.07309e-09,-1.42159e-13,9211.19,18.6318], Tmin=(837.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.7982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=CJO)"""),
)

species(
    label = 'C=C(O)C(O)=[C]C(8868)',
    structure = SMILES('C=C(O)C(O)=[C]C'),
    E0 = (-132.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(0.845835,'amu*angstrom^2'), symmetry=1, barrier=(19.4474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845982,'amu*angstrom^2'), symmetry=1, barrier=(19.4508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846479,'amu*angstrom^2'), symmetry=1, barrier=(19.4622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.847279,'amu*angstrom^2'), symmetry=1, barrier=(19.4806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.524415,0.0878455,-0.000103087,5.81171e-08,-1.23276e-11,-15708.8,24.9954], Tmin=(100,'K'), Tmax=(1274.34,'K')), NASAPolynomial(coeffs=[22.1925,0.00883947,-1.02661e-06,-1.71752e-11,7.3761e-15,-20873.4,-87.6545], Tmin=(1274.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)C(O)=CC(8869)',
    structure = SMILES('[CH]=C(O)C(O)=CC'),
    E0 = (-122.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750],'cm^-1')),
        HinderedRotor(inertia=(0.89411,'amu*angstrom^2'), symmetry=1, barrier=(20.5573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891563,'amu*angstrom^2'), symmetry=1, barrier=(20.4988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891707,'amu*angstrom^2'), symmetry=1, barrier=(20.5021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895086,'amu*angstrom^2'), symmetry=1, barrier=(20.5798,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750107,0.088929,-0.000102429,5.58194e-08,-1.13288e-11,-14584.3,25.736], Tmin=(100,'K'), Tmax=(1366.7,'K')), NASAPolynomial(coeffs=[23.6362,0.00641342,3.64202e-07,-2.87664e-10,2.56898e-14,-20209.4,-95.7252], Tmin=(1366.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C(O)=CC(8870)',
    structure = SMILES('[CH2]C([O])C(O)=CC'),
    E0 = (2.03722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,366.378,366.418],'cm^-1')),
        HinderedRotor(inertia=(0.141215,'amu*angstrom^2'), symmetry=1, barrier=(13.4665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14135,'amu*angstrom^2'), symmetry=1, barrier=(13.4659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141244,'amu*angstrom^2'), symmetry=1, barrier=(13.4656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141248,'amu*angstrom^2'), symmetry=1, barrier=(13.4664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283774,0.0707516,-5.05058e-05,5.49106e-09,5.59804e-12,388.701,28.9881], Tmin=(100,'K'), Tmax=(967.788,'K')), NASAPolynomial(coeffs=[19.2911,0.0169688,-5.54956e-06,9.76858e-10,-6.95203e-14,-4450.66,-68.0839], Tmin=(967.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.03722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C(O)=[C]C(8871)',
    structure = SMILES('[CH2]C(O)C(O)=[C]C'),
    E0 = (9.51822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,343.879],'cm^-1')),
        HinderedRotor(inertia=(0.140603,'amu*angstrom^2'), symmetry=1, barrier=(11.7962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140581,'amu*angstrom^2'), symmetry=1, barrier=(11.7958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140604,'amu*angstrom^2'), symmetry=1, barrier=(11.796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140591,'amu*angstrom^2'), symmetry=1, barrier=(11.7954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140592,'amu*angstrom^2'), symmetry=1, barrier=(11.7954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0392474,0.0808275,-8.62719e-05,4.61606e-08,-9.57735e-12,1297.04,30.8137], Tmin=(100,'K'), Tmax=(1187.87,'K')), NASAPolynomial(coeffs=[18.9048,0.0170356,-5.7176e-06,9.51058e-10,-6.2495e-14,-3203.56,-63.8418], Tmin=(1187.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.51822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=C([O])C(O)[CH]C(8872)',
    structure = SMILES('C=C([O])C(O)[CH]C'),
    E0 = (-89.9604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,297.802,301.765,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00184902,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00187391,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280963,'amu*angstrom^2'), symmetry=1, barrier=(17.9338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275302,'amu*angstrom^2'), symmetry=1, barrier=(17.9131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773408,0.0603072,-3.14832e-05,-7.79839e-09,8.95214e-12,-10694,31.2248], Tmin=(100,'K'), Tmax=(968.867,'K')), NASAPolynomial(coeffs=[16.492,0.0198534,-6.69236e-06,1.18047e-09,-8.316e-14,-14887,-50.0315], Tmin=(968.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C(O)C(O)[CH]C(8873)',
    structure = SMILES('[CH]=C(O)C(O)[CH]C'),
    E0 = (19.3309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,3120,650,792.5,1650,402.583],'cm^-1')),
        HinderedRotor(inertia=(0.00104012,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115409,'amu*angstrom^2'), symmetry=1, barrier=(13.2738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11539,'amu*angstrom^2'), symmetry=1, barrier=(13.2731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11539,'amu*angstrom^2'), symmetry=1, barrier=(13.2743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115369,'amu*angstrom^2'), symmetry=1, barrier=(13.2728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35222,0.0675747,-3.97791e-05,-8.68595e-09,1.16147e-11,2467.87,31.639], Tmin=(100,'K'), Tmax=(944.712,'K')), NASAPolynomial(coeffs=[20.3737,0.0138257,-3.69596e-06,6.11973e-10,-4.46799e-14,-2699.43,-71.1412], Tmin=(944.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.3309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_P)"""),
)

species(
    label = 'C[C]=C(O)[C](C)O(8874)',
    structure = SMILES('C[C]C(O)=C(C)O'),
    E0 = (-33.4494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.743648,'amu*angstrom^2'), symmetry=1, barrier=(17.0979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743631,'amu*angstrom^2'), symmetry=1, barrier=(17.0976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743388,'amu*angstrom^2'), symmetry=1, barrier=(17.0919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743426,'amu*angstrom^2'), symmetry=1, barrier=(17.0928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743388,'amu*angstrom^2'), symmetry=1, barrier=(17.092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.956428,0.0868802,-9.04217e-05,4.47944e-08,-8.25531e-12,-3825.61,29.7053], Tmin=(100,'K'), Tmax=(1531.6,'K')), NASAPolynomial(coeffs=[24.0031,0.00827556,-2.96485e-07,-1.55319e-10,1.54736e-14,-9897.26,-96.2123], Tmin=(1531.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.4494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=C(O)C([CH2])O(8875)',
    structure = SMILES('[CH2]C=C(O)C([CH2])O'),
    E0 = (-76.8244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,311.868],'cm^-1')),
        HinderedRotor(inertia=(0.0017274,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.255559,'amu*angstrom^2'), symmetry=1, barrier=(17.4353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250908,'amu*angstrom^2'), symmetry=1, barrier=(17.4682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251704,'amu*angstrom^2'), symmetry=1, barrier=(17.4463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253418,'amu*angstrom^2'), symmetry=1, barrier=(17.4366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.144358,0.0719598,-4.67633e-05,-4.40047e-09,1.07493e-11,-9089.27,29.9615], Tmin=(100,'K'), Tmax=(941.806,'K')), NASAPolynomial(coeffs=[21.4901,0.0129777,-3.27468e-06,5.28015e-10,-3.87296e-14,-14514.8,-79.1981], Tmin=(941.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.8244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C(O)[C](C)O(8876)',
    structure = SMILES('[CH2]C=C(O)[C](C)O'),
    E0 = (-111.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,400.14],'cm^-1')),
        HinderedRotor(inertia=(0.125207,'amu*angstrom^2'), symmetry=1, barrier=(14.0912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125098,'amu*angstrom^2'), symmetry=1, barrier=(14.0963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123551,'amu*angstrom^2'), symmetry=1, barrier=(14.106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12439,'amu*angstrom^2'), symmetry=1, barrier=(14.0841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123248,'amu*angstrom^2'), symmetry=1, barrier=(14.0803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.303811,0.0806589,-8.15741e-05,4.06288e-08,-7.74696e-12,-13277.8,30.367], Tmin=(100,'K'), Tmax=(1379.32,'K')), NASAPolynomial(coeffs=[21.1156,0.0138305,-3.77368e-06,5.48421e-10,-3.34848e-14,-18738.4,-78.2323], Tmin=(1379.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(O)C(=C)O(8877)',
    structure = SMILES('C=C(O)C(=C)O'),
    E0 = (-333.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(1.01335,'amu*angstrom^2'), symmetry=1, barrier=(23.299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01377,'amu*angstrom^2'), symmetry=1, barrier=(23.3085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01324,'amu*angstrom^2'), symmetry=1, barrier=(23.2963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.498033,0.0751665,-8.30988e-05,4.22746e-08,-7.8095e-12,-39972.2,21.9417], Tmin=(100,'K'), Tmax=(1586.91,'K')), NASAPolynomial(coeffs=[21.8153,0.0015529,2.90201e-06,-7.522e-10,5.55296e-14,-44866.9,-89.1201], Tmin=(1586.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=C(O)C(C)=O(8878)',
    structure = SMILES('CC=C(O)C(C)=O'),
    E0 = (-369.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706729,0.0660995,-5.49688e-05,2.3391e-08,-4.00086e-12,-44265.9,23.7086], Tmin=(100,'K'), Tmax=(1391.19,'K')), NASAPolynomial(coeffs=[15.2261,0.0243525,-9.95644e-06,1.82065e-09,-1.24604e-13,-48305.7,-51.1328], Tmin=(1391.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs)"""),
)

species(
    label = 'C=C(O)C(O)[C]C(8879)',
    structure = SMILES('C=C(O)C(O)[C]C'),
    E0 = (0.198124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201235,0.0823602,-9.02139e-05,5.21228e-08,-1.20015e-11,161.704,37.0062], Tmin=(100,'K'), Tmax=(1057.37,'K')), NASAPolynomial(coeffs=[15.0406,0.0262241,-1.05795e-05,1.91443e-09,-1.3065e-13,-2976.48,-35.4132], Tmin=(1057.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.198124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]C(O)C(O)=CC(8880)',
    structure = SMILES('[CH]C(O)C(O)=CC'),
    E0 = (10.3887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.361348,0.0803021,-9.09026e-05,5.41511e-08,-1.28617e-11,1380.35,26.8385], Tmin=(100,'K'), Tmax=(1025.05,'K')), NASAPolynomial(coeffs=[14.4094,0.0254823,-1.06811e-05,1.97638e-09,-1.36615e-13,-1499.6,-41.2826], Tmin=(1025.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.3887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CC1C[C]=C1O(8881)',
    structure = SMILES('CC1C[C]=C1O'),
    E0 = (150.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68823,0.0323573,3.95721e-05,-8.33617e-08,3.68694e-11,18218.4,19.435], Tmin=(100,'K'), Tmax=(933.41,'K')), NASAPolynomial(coeffs=[17.7029,0.0108886,-1.71401e-06,2.54769e-10,-2.36683e-14,13174.3,-67.7277], Tmin=(933.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'CC1CC([O])=C1O(8882)',
    structure = SMILES('CC1CC([O])=C1O'),
    E0 = (-178.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.079,0.042591,3.40604e-05,-9.05711e-08,4.25705e-11,-21355,22.7072], Tmin=(100,'K'), Tmax=(920.381,'K')), NASAPolynomial(coeffs=[22.3403,0.00646515,1.22011e-06,-3.49733e-10,1.87182e-14,-27652.2,-91.0512], Tmin=(920.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC1[C]=C(O)C1(8883)',
    structure = SMILES('CC1[C]=C(O)C1'),
    E0 = (150.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68823,0.0323573,3.95721e-05,-8.33617e-08,3.68694e-11,18218.4,19.435], Tmin=(100,'K'), Tmax=(933.41,'K')), NASAPolynomial(coeffs=[17.7029,0.0108886,-1.71401e-06,2.54769e-10,-2.36683e-14,13174.3,-67.7277], Tmin=(933.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'CC1CC(O)=C1[O](8884)',
    structure = SMILES('CC1CC(O)=C1[O]'),
    E0 = (-178.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.079,0.042591,3.40604e-05,-9.05711e-08,4.25705e-11,-21355,22.7072], Tmin=(100,'K'), Tmax=(920.381,'K')), NASAPolynomial(coeffs=[22.3403,0.00646515,1.22011e-06,-3.49733e-10,1.87182e-14,-27652.2,-91.0512], Tmin=(920.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'OC1[CH]CC=1O(8885)',
    structure = SMILES('OC1[CH]CC=1O'),
    E0 = (-84.7475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57713,0.0289172,5.42282e-05,-1.13251e-07,5.20534e-11,-10082.7,20.6096], Tmin=(100,'K'), Tmax=(910.304,'K')), NASAPolynomial(coeffs=[24.3088,-0.0068051,7.36182e-06,-1.49648e-09,9.64789e-14,-16879.8,-101.524], Tmin=(910.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.7475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CCJCO)"""),
)

species(
    label = 'C[C]1CC(O)=C1O(8886)',
    structure = SMILES('C[C]1CC(O)=C1O'),
    E0 = (-163.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986748,0.0398836,5.3074e-05,-1.17668e-07,5.38848e-11,-19570.9,22.4258], Tmin=(100,'K'), Tmax=(919.458,'K')), NASAPolynomial(coeffs=[25.6013,0.00141814,3.88401e-06,-8.36793e-10,4.96988e-14,-26997.8,-110.031], Tmin=(919.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-163.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC1[CH]C(O)=C1O(8887)',
    structure = SMILES('CC1[CH]C(O)=C1O'),
    E0 = (-116.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904468,0.0413704,5.02283e-05,-1.17375e-07,5.46521e-11,-13875.3,24.3297], Tmin=(100,'K'), Tmax=(915.057,'K')), NASAPolynomial(coeffs=[26.6414,-0.00118709,5.33095e-06,-1.12982e-09,7.0611e-14,-21513.9,-113.553], Tmin=(915.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1CC(O)=C1O(8888)',
    structure = SMILES('[CH2]C1CC(O)=C1O'),
    E0 = (-111.319,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741935,0.0480724,2.99777e-05,-9.6883e-08,4.80368e-11,-13249.2,24.2], Tmin=(100,'K'), Tmax=(902.955,'K')), NASAPolynomial(coeffs=[25.9647,-0.000424585,5.49103e-06,-1.24364e-09,8.32039e-14,-20382.1,-109.186], Tmin=(902.955,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(Isobutyl)"""),
)

species(
    label = 'CC1CC([O])[C]1O(8889)',
    structure = SMILES('CC1CC([O])[C]1O'),
    E0 = (43.0991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11662,0.050699,-4.51997e-06,-3.12523e-08,1.57373e-11,5298.71,24.741], Tmin=(100,'K'), Tmax=(990.093,'K')), NASAPolynomial(coeffs=[14.8726,0.0244088,-9.05624e-06,1.67553e-09,-1.20153e-13,1139.41,-48.7354], Tmin=(990.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1C[C](O)C1[O](8890)',
    structure = SMILES('CC1C[C](O)C1[O]'),
    E0 = (43.0991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11662,0.050699,-4.51997e-06,-3.12523e-08,1.57373e-11,5298.71,24.741], Tmin=(100,'K'), Tmax=(990.093,'K')), NASAPolynomial(coeffs=[14.8726,0.0244088,-9.05624e-06,1.67553e-09,-1.20153e-13,1139.41,-48.7354], Tmin=(990.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[C]1C[C](O)C1O(8891)',
    structure = SMILES('C[C]1C[C](O)C1O'),
    E0 = (-34.6878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23165,0.0479754,1.13341e-06,-3.61704e-08,1.73159e-11,-4060.83,25.4083], Tmin=(100,'K'), Tmax=(986.476,'K')), NASAPolynomial(coeffs=[14.4381,0.0245509,-9.05605e-06,1.67285e-09,-1.20017e-13,-8132.21,-45.5549], Tmin=(986.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.6878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC1[CH]C(O)[C]1O(8892)',
    structure = SMILES('CC1[CH]C(O)[C]1O'),
    E0 = (12.6403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15488,0.0494119,-1.62046e-06,-3.58439e-08,1.79863e-11,1634.52,27.2916], Tmin=(100,'K'), Tmax=(971.963,'K')), NASAPolynomial(coeffs=[15.3833,0.0221034,-7.69865e-06,1.40072e-09,-1.00823e-13,-2607.35,-48.5405], Tmin=(971.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.6403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = 'CC1[CH][C](O)C1O(8893)',
    structure = SMILES('CC1[CH][C](O)C1O'),
    E0 = (12.6403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15488,0.0494119,-1.62046e-06,-3.58439e-08,1.79863e-11,1634.52,27.2916], Tmin=(100,'K'), Tmax=(971.963,'K')), NASAPolynomial(coeffs=[15.3833,0.0221034,-7.69865e-06,1.40072e-09,-1.00823e-13,-2607.35,-48.5405], Tmin=(971.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.6403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1C[C](O)C1O(8894)',
    structure = SMILES('[CH2]C1C[C](O)C1O'),
    E0 = (17.8206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990533,0.0561639,-2.22059e-05,-1.46308e-08,1.09085e-11,2260.74,27.1667], Tmin=(100,'K'), Tmax=(955.66,'K')), NASAPolynomial(coeffs=[14.5988,0.0230445,-7.63955e-06,1.31042e-09,-9.01585e-14,-1428.83,-43.5635], Tmin=(955.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = 'C[C]1CC(O)[C]1O(8895)',
    structure = SMILES('C[C]1CC(O)[C]1O'),
    E0 = (-34.6878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23165,0.0479754,1.13341e-06,-3.61704e-08,1.73159e-11,-4060.83,25.4083], Tmin=(100,'K'), Tmax=(986.476,'K')), NASAPolynomial(coeffs=[14.4381,0.0245509,-9.05605e-06,1.67285e-09,-1.20017e-13,-8132.21,-45.5549], Tmin=(986.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.6878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C1CC(O)[C]1O(8896)',
    structure = SMILES('[CH2]C1CC(O)[C]1O'),
    E0 = (17.8206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990533,0.0561639,-2.22059e-05,-1.46308e-08,1.09085e-11,2260.74,27.1667], Tmin=(100,'K'), Tmax=(955.66,'K')), NASAPolynomial(coeffs=[14.5988,0.0230445,-7.63955e-06,1.31042e-09,-9.01585e-14,-1428.83,-43.5635], Tmin=(955.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.8206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(C2CsJOH)"""),
)

species(
    label = 'OC1CCC=1O(8897)',
    structure = SMILES('OC1CCC=1O'),
    E0 = (-284.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38981,0.0340504,4.58553e-05,-1.0403e-07,4.82865e-11,-34119.7,18.092], Tmin=(100,'K'), Tmax=(913.9,'K')), NASAPolynomial(coeffs=[23.7325,-0.00248879,5.29435e-06,-1.10519e-09,6.97775e-14,-40761.4,-101.682], Tmin=(913.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-284.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]C(C)C(O)=[C]O(8898)',
    structure = SMILES('[CH2]C(C)C(O)=[C]O'),
    E0 = (-26.7986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,353.484],'cm^-1')),
        HinderedRotor(inertia=(0.169124,'amu*angstrom^2'), symmetry=1, barrier=(14.9959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169124,'amu*angstrom^2'), symmetry=1, barrier=(14.9959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169124,'amu*angstrom^2'), symmetry=1, barrier=(14.9959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169124,'amu*angstrom^2'), symmetry=1, barrier=(14.9959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169124,'amu*angstrom^2'), symmetry=1, barrier=(14.9959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02981,0.0846184,-8.77689e-05,4.37546e-08,-8.01302e-12,-3019.65,34.9679], Tmin=(100,'K'), Tmax=(1595.72,'K')), NASAPolynomial(coeffs=[22.0851,0.0084366,9.884e-07,-4.89781e-10,4.09104e-14,-8074.46,-80.0742], Tmin=(1595.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.7986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]CC(O)=[C]O(8899)',
    structure = SMILES('C[CH]CC(O)=[C]O'),
    E0 = (-29.4628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,347.627,347.629],'cm^-1')),
        HinderedRotor(inertia=(0.00139401,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148202,'amu*angstrom^2'), symmetry=1, barrier=(12.6914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148112,'amu*angstrom^2'), symmetry=1, barrier=(12.6947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147487,'amu*angstrom^2'), symmetry=1, barrier=(12.6927,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148276,'amu*angstrom^2'), symmetry=1, barrier=(12.6938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.193066,0.078205,-7.98317e-05,4.08165e-08,-7.94537e-12,-3380.61,33.1542], Tmin=(100,'K'), Tmax=(1398.04,'K')), NASAPolynomial(coeffs=[19.558,0.0141542,-3.01978e-06,3.30357e-10,-1.57248e-14,-8166.32,-66.116], Tmin=(1398.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.4628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(O)=C(O)[CH]C(8830)',
    structure = SMILES('[CH2]C(O)=C(O)[CH]C'),
    E0 = (-122.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,325,375,415,465,420,450,1700,1750,323.464],'cm^-1')),
        HinderedRotor(inertia=(0.0016112,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208938,'amu*angstrom^2'), symmetry=1, barrier=(15.5142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208954,'amu*angstrom^2'), symmetry=1, barrier=(15.5143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208969,'amu*angstrom^2'), symmetry=1, barrier=(15.5143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20896,'amu*angstrom^2'), symmetry=1, barrier=(15.5143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141743,0.0688338,-2.99459e-05,-3.072e-08,2.30897e-11,-14621,29.0351], Tmin=(100,'K'), Tmax=(905.863,'K')), NASAPolynomial(coeffs=[24.2058,0.00628117,1.26066e-06,-4.23747e-10,2.91719e-14,-20774,-94.5791], Tmin=(905.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCJCO)"""),
)

species(
    label = 'CC1CC(O)C1=O(8900)',
    structure = SMILES('CC1CC(O)C1=O'),
    E0 = (-315.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85324,0.0259792,6.87174e-05,-1.08982e-07,4.36755e-11,-37831.6,24.6967], Tmin=(100,'K'), Tmax=(960.026,'K')), NASAPolynomial(coeffs=[15.8263,0.0210483,-6.83928e-06,1.30513e-09,-1.00844e-13,-42970.1,-54.9347], Tmin=(960.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = 'CC1CC(=O)C1O(8901)',
    structure = SMILES('CC1CC(=O)C1O'),
    E0 = (-318.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90151,0.0249762,7.0621e-05,-1.09924e-07,4.37071e-11,-38228,24.3817], Tmin=(100,'K'), Tmax=(962.458,'K')), NASAPolynomial(coeffs=[15.4426,0.0218312,-7.2843e-06,1.39648e-09,-1.07413e-13,-43295.4,-53.2124], Tmin=(962.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-318.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = 'C=CC(O)=[C]C(8902)',
    structure = SMILES('C=CC(O)=[C]C'),
    E0 = (81.8445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.828236,'amu*angstrom^2'), symmetry=1, barrier=(19.0428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82817,'amu*angstrom^2'), symmetry=1, barrier=(19.0412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825568,'amu*angstrom^2'), symmetry=1, barrier=(18.9814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.525329,0.0648927,-6.35102e-05,3.12651e-08,-5.91865e-12,9978.53,21.9937], Tmin=(100,'K'), Tmax=(1400.17,'K')), NASAPolynomial(coeffs=[17.0677,0.013284,-3.56118e-06,5.02273e-10,-2.97352e-14,5772.57,-61.8587], Tmin=(1400.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.8445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(O)=C(C)[O](8903)',
    structure = SMILES('C=CC(O)=C(C)[O]'),
    E0 = (-232.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.75227,'amu*angstrom^2'), symmetry=1, barrier=(17.2962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752021,'amu*angstrom^2'), symmetry=1, barrier=(17.2904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750681,'amu*angstrom^2'), symmetry=1, barrier=(17.2596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.402824,0.0788974,-8.23181e-05,4.15686e-08,-7.8506e-12,-27819.6,26.4099], Tmin=(100,'K'), Tmax=(1484.65,'K')), NASAPolynomial(coeffs=[21.2474,0.00942733,-8.75262e-07,-5.61749e-11,9.57068e-15,-33020.5,-82.4609], Tmin=(1484.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[C]=C(C)O(8904)',
    structure = SMILES('C=C[C]=C(C)O'),
    E0 = (42.3325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.00383,'amu*angstrom^2'), symmetry=1, barrier=(23.08,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00368,'amu*angstrom^2'), symmetry=1, barrier=(23.0766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00467,'amu*angstrom^2'), symmetry=1, barrier=(23.0993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.218844,0.0663473,-6.42839e-05,3.09993e-08,-5.63176e-12,5241.97,23.0355], Tmin=(100,'K'), Tmax=(1551.19,'K')), NASAPolynomial(coeffs=[17.6522,0.0116224,-1.91746e-06,1.3535e-10,-3.13713e-15,1008.9,-64.9354], Tmin=(1551.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.3325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC([O])=C(C)O(8905)',
    structure = SMILES('C=CC([O])=C(C)O'),
    E0 = (-232.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.75227,'amu*angstrom^2'), symmetry=1, barrier=(17.2962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752021,'amu*angstrom^2'), symmetry=1, barrier=(17.2904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.750681,'amu*angstrom^2'), symmetry=1, barrier=(17.2596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.402824,0.0788974,-8.23181e-05,4.15686e-08,-7.8506e-12,-27819.6,26.4099], Tmin=(100,'K'), Tmax=(1484.65,'K')), NASAPolynomial(coeffs=[21.2474,0.00942733,-8.75262e-07,-5.61749e-11,9.57068e-15,-33020.5,-82.4609], Tmin=(1484.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(O)=[C]O(8906)',
    structure = SMILES('C=CC(O)=[C]O'),
    E0 = (-89.0201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.999743,'amu*angstrom^2'), symmetry=1, barrier=(22.9861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.998859,'amu*angstrom^2'), symmetry=1, barrier=(22.9657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99907,'amu*angstrom^2'), symmetry=1, barrier=(22.9706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383059,0.0711318,-8.04555e-05,4.12386e-08,-7.59281e-12,-10527,25.9678], Tmin=(100,'K'), Tmax=(1617.64,'K')), NASAPolynomial(coeffs=[20.9015,-0.000781229,4.10752e-06,-9.80442e-10,7.0779e-14,-14890.4,-79.1573], Tmin=(1617.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.0201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C2H3(28)(29)',
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
    label = 'CC(O)=[C]O(8907)',
    structure = SMILES('CC(O)=[C]O'),
    E0 = (-177.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.791724,'amu*angstrom^2'), symmetry=1, barrier=(18.2033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.789048,'amu*angstrom^2'), symmetry=1, barrier=(18.1418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.79285,'amu*angstrom^2'), symmetry=1, barrier=(18.2292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467636,0.0562861,-6.0594e-05,3.01606e-08,-5.41538e-12,-21200.6,22.9034], Tmin=(100,'K'), Tmax=(1662.12,'K')), NASAPolynomial(coeffs=[16.7989,0.00173007,2.40661e-06,-6.30016e-10,4.63132e-14,-24522.5,-57.8448], Tmin=(1662.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C(O)=C(C)O(8908)',
    structure = SMILES('C=[C]C(O)=C(C)O'),
    E0 = (-171.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,180],'cm^-1')),
        HinderedRotor(inertia=(0.975068,'amu*angstrom^2'), symmetry=1, barrier=(22.4187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973779,'amu*angstrom^2'), symmetry=1, barrier=(22.3891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.974247,'amu*angstrom^2'), symmetry=1, barrier=(22.3998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975367,'amu*angstrom^2'), symmetry=1, barrier=(22.4256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.770963,0.0886743,-0.00010201,5.58683e-08,-1.13508e-11,-20448.2,25.8162], Tmin=(100,'K'), Tmax=(1386.19,'K')), NASAPolynomial(coeffs=[23.0055,0.00688787,7.49759e-07,-4.10143e-10,3.58348e-14,-25773.9,-92.0895], Tmin=(1386.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC(O)=C(C)O(8909)',
    structure = SMILES('[CH]=CC(O)=C(C)O'),
    E0 = (-123.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,325,375,415,465,420,450,1700,1750],'cm^-1')),
        HinderedRotor(inertia=(0.8863,'amu*angstrom^2'), symmetry=1, barrier=(20.3778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886581,'amu*angstrom^2'), symmetry=1, barrier=(20.3843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885786,'amu*angstrom^2'), symmetry=1, barrier=(20.366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885488,'amu*angstrom^2'), symmetry=1, barrier=(20.3591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01056,0.0883913,-9.85407e-05,5.11917e-08,-9.78072e-12,-14649.6,27.4918], Tmin=(100,'K'), Tmax=(1495.46,'K')), NASAPolynomial(coeffs=[24.6933,0.003968,1.85801e-06,-5.72936e-10,4.4376e-14,-20585.1,-101], Tmin=(1495.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-123.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=C[C](O)C(C)[O](8910)',
    structure = SMILES('[CH2]C=C(O)C(C)[O]'),
    E0 = (-58.0527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,409.132,410.009],'cm^-1')),
        HinderedRotor(inertia=(0.143988,'amu*angstrom^2'), symmetry=1, barrier=(17.0917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144195,'amu*angstrom^2'), symmetry=1, barrier=(17.0892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143403,'amu*angstrom^2'), symmetry=1, barrier=(17.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143745,'amu*angstrom^2'), symmetry=1, barrier=(17.0916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.452509,0.062859,-2.13948e-05,-2.82247e-08,1.83545e-11,-6840.59,28.0297], Tmin=(100,'K'), Tmax=(950.395,'K')), NASAPolynomial(coeffs=[20.6548,0.0148924,-4.1815e-06,7.30725e-10,-5.50816e-14,-12354.4,-77.2122], Tmin=(950.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.0527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC([O])[C](C)O(8911)',
    structure = SMILES('C=CC([O])[C](C)O'),
    E0 = (21.4095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1380,1390,370,380,2900,435,364.151,364.225,364.248],'cm^-1')),
        HinderedRotor(inertia=(0.0012708,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114397,'amu*angstrom^2'), symmetry=1, barrier=(10.7714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114442,'amu*angstrom^2'), symmetry=1, barrier=(10.7715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114452,'amu*angstrom^2'), symmetry=1, barrier=(10.7714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72099,0.0703653,-6.59376e-05,3.27333e-08,-6.59756e-12,2694.48,28.9508], Tmin=(100,'K'), Tmax=(1185.82,'K')), NASAPolynomial(coeffs=[13.4716,0.0273546,-1.15306e-05,2.14531e-09,-1.48778e-13,-329.451,-34.7362], Tmin=(1185.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C=[C]C(O)[C](C)O(8912)',
    structure = SMILES('C=[C]C(O)[C](C)O'),
    E0 = (28.8905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1380,1390,370,380,2900,435,309.689,309.72],'cm^-1')),
        HinderedRotor(inertia=(0.00175756,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133358,'amu*angstrom^2'), symmetry=1, barrier=(9.07795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133358,'amu*angstrom^2'), symmetry=1, barrier=(9.07792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133344,'amu*angstrom^2'), symmetry=1, barrier=(9.07804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133373,'amu*angstrom^2'), symmetry=1, barrier=(9.07792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622586,0.07768,-9.16449e-05,6.00393e-08,-1.60309e-11,3593.41,29.9795], Tmin=(100,'K'), Tmax=(906.89,'K')), NASAPolynomial(coeffs=[11.347,0.0303773,-1.34047e-05,2.52295e-09,-1.75248e-13,1648.27,-20.7112], Tmin=(906.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.8905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)[C](C)O(8913)',
    structure = SMILES('[CH]=CC(O)[C](C)O'),
    E0 = (38.1448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1380,1390,370,380,2900,435,3120,650,792.5,1650,373.45],'cm^-1')),
        HinderedRotor(inertia=(0.0974129,'amu*angstrom^2'), symmetry=1, barrier=(9.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974141,'amu*angstrom^2'), symmetry=1, barrier=(9.64023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974107,'amu*angstrom^2'), symmetry=1, barrier=(9.6402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974126,'amu*angstrom^2'), symmetry=1, barrier=(9.64019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0974145,'amu*angstrom^2'), symmetry=1, barrier=(9.6402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554541,0.0769845,-8.5105e-05,5.05151e-08,-1.20922e-11,4710.88,30.1489], Tmin=(100,'K'), Tmax=(1011.63,'K')), NASAPolynomial(coeffs=[13.0592,0.0275411,-1.17928e-05,2.20245e-09,-1.52906e-13,2180.85,-30.3236], Tmin=(1011.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.1448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[C](O)C(C)O(8914)',
    structure = SMILES('[CH]C=C(O)C(C)O'),
    E0 = (-69.2282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288764,0.069286,-3.55241e-05,-1.0324e-08,1.09573e-11,-8181.49,29.6375], Tmin=(100,'K'), Tmax=(963.082,'K')), NASAPolynomial(coeffs=[18.3751,0.0228204,-7.78073e-06,1.36318e-09,-9.54502e-14,-12994,-63.8372], Tmin=(963.082,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.2282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]CC(O)=C(C)[O](8915)',
    structure = SMILES('[CH2]CC(O)=C(C)[O]'),
    E0 = (-138.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,318.087,318.17],'cm^-1')),
        HinderedRotor(inertia=(0.0016365,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143797,'amu*angstrom^2'), symmetry=1, barrier=(10.5463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141019,'amu*angstrom^2'), symmetry=1, barrier=(10.5785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137438,'amu*angstrom^2'), symmetry=1, barrier=(10.5767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.149023,0.0765103,-7.49082e-05,3.66676e-08,-6.8673e-12,-16509.4,30.1846], Tmin=(100,'K'), Tmax=(1435.78,'K')), NASAPolynomial(coeffs=[19.7458,0.0146968,-3.65663e-06,4.85292e-10,-2.76616e-14,-21563.9,-70.6997], Tmin=(1435.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(RCCJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C(O)=C(C)[O](8916)',
    structure = SMILES('C[CH]C(O)=C(C)[O]'),
    E0 = (-143.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,363.148,363.217],'cm^-1')),
        HinderedRotor(inertia=(0.0012773,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111173,'amu*angstrom^2'), symmetry=1, barrier=(10.4121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111176,'amu*angstrom^2'), symmetry=1, barrier=(10.4118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111176,'amu*angstrom^2'), symmetry=1, barrier=(10.4116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539627,0.0644069,-3.47983e-05,-1.16762e-08,1.24888e-11,-17178.7,28.4164], Tmin=(100,'K'), Tmax=(930.903,'K')), NASAPolynomial(coeffs=[18.8584,0.0154804,-3.95891e-06,6.11722e-10,-4.24718e-14,-21880,-65.5818], Tmin=(930.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(O)=CO(8917)',
    structure = SMILES('C=CC(O)=CO'),
    E0 = (-328.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.10776,'amu*angstrom^2'), symmetry=1, barrier=(25.4697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10596,'amu*angstrom^2'), symmetry=1, barrier=(25.4283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10762,'amu*angstrom^2'), symmetry=1, barrier=(25.4663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.579041,0.0556429,-5.74162e-06,-5.97426e-08,3.54133e-11,-39399.5,19.3796], Tmin=(100,'K'), Tmax=(893.729,'K')), NASAPolynomial(coeffs=[26.6725,-0.00674585,7.6733e-06,-1.64806e-09,1.12979e-13,-46236,-115.728], Tmin=(893.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[C]C(O)=C(C)O(8918)',
    structure = SMILES('C[C]C(O)=C(C)O'),
    E0 = (-54.1354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.966313,'amu*angstrom^2'), symmetry=1, barrier=(22.2174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966099,'amu*angstrom^2'), symmetry=1, barrier=(22.2125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966107,'amu*angstrom^2'), symmetry=1, barrier=(22.2127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966385,'amu*angstrom^2'), symmetry=1, barrier=(22.2191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.965985,'amu*angstrom^2'), symmetry=1, barrier=(22.2099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.213524,0.0860288,-9.12235e-05,4.6811e-08,-8.72027e-12,-6353.39,35.7513], Tmin=(100,'K'), Tmax=(985.257,'K')), NASAPolynomial(coeffs=[19.0858,0.0189916,-6.39011e-06,1.0657e-09,-7.05936e-14,-10705.6,-59.8575], Tmin=(985.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.1354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC(O)=C(C)O(8919)',
    structure = SMILES('[CH]CC(O)=C(C)O'),
    E0 = (-31.6994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.169862,0.0845945,-9.03424e-05,4.54212e-08,-8.05134e-12,-3656.1,26.4378], Tmin=(100,'K'), Tmax=(984.384,'K')), NASAPolynomial(coeffs=[19.7977,0.0161185,-5.29207e-06,8.87653e-10,-5.96886e-14,-8200.69,-72.6957], Tmin=(984.384,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.6994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C=C[CH]C(C)=O(8920)',
    structure = SMILES('C=CC=C(C)[O]'),
    E0 = (-18.8581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.666427,'amu*angstrom^2'), symmetry=1, barrier=(15.3225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.670589,'amu*angstrom^2'), symmetry=1, barrier=(15.4182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37494,0.0475809,-1.44791e-05,-2.06114e-08,1.30873e-11,-2164.33,20.7811], Tmin=(100,'K'), Tmax=(943.579,'K')), NASAPolynomial(coeffs=[14.6224,0.0165059,-4.95387e-06,8.3126e-10,-5.82183e-14,-5780.96,-48.2777], Tmin=(943.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC([O])C(C)=O(8921)',
    structure = SMILES('C=CC([O])C(C)=O'),
    E0 = (-81.4752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.663848,'amu*angstrom^2'), symmetry=1, barrier=(15.2632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00963721,'amu*angstrom^2'), symmetry=1, barrier=(4.88617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300233,'amu*angstrom^2'), symmetry=1, barrier=(15.2585,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16765,0.054059,-3.76151e-05,1.31365e-08,-1.84659e-12,-9690.38,28.8141], Tmin=(100,'K'), Tmax=(1662.36,'K')), NASAPolynomial(coeffs=[14.4927,0.0219958,-8.68308e-06,1.53355e-09,-1.01622e-13,-14120.5,-42.2436], Tmin=(1662.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.4752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'CH3CO(55)(54)',
    structure = SMILES('C[C]=O'),
    E0 = (-22.2282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,2865.56],'cm^-1')),
        HinderedRotor(inertia=(0.0188671,'amu*angstrom^2'), symmetry=1, barrier=(18.7749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03587,0.000877295,3.071e-05,-3.92476e-08,1.52969e-11,-2682.07,7.86177], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.31372,0.00917378,-3.32204e-06,5.39475e-10,-3.24524e-14,-3645.04,-1.67576], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-22.2282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'hydroxyl1allyl(7170)',
    structure = SMILES('[CH2]C=CO'),
    E0 = (-25.0312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.978343,'amu*angstrom^2'), symmetry=1, barrier=(22.494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.975892,'amu*angstrom^2'), symmetry=1, barrier=(22.4377,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5477,0.0242083,5.90889e-06,-2.6575e-08,1.23966e-11,-2951.28,13.5324], Tmin=(100,'K'), Tmax=(959.506,'K')), NASAPolynomial(coeffs=[10.1178,0.011534,-3.79865e-06,6.81299e-10,-4.93305e-14,-5273.27,-27.2058], Tmin=(959.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.0312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""hydroxyl1allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC(=O)[CH]O(8922)',
    structure = SMILES('CC([O])=CO'),
    E0 = (-279.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.840713,'amu*angstrom^2'), symmetry=1, barrier=(19.3296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.836685,'amu*angstrom^2'), symmetry=1, barrier=(19.237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63463,0.0388026,-7.04439e-07,-4.36344e-08,2.49535e-11,-33509.1,17.0756], Tmin=(100,'K'), Tmax=(892.463,'K')), NASAPolynomial(coeffs=[18.6178,-0.000591679,3.7837e-06,-8.79826e-10,6.12415e-14,-38003,-71.12], Tmin=(892.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(O)[C]=O(8923)',
    structure = SMILES('C=CC(O)[C]=O'),
    E0 = (-110.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,280.901],'cm^-1')),
        HinderedRotor(inertia=(0.00212029,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218986,'amu*angstrom^2'), symmetry=1, barrier=(12.2713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217281,'amu*angstrom^2'), symmetry=1, barrier=(12.274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77602,0.0418188,-2.38139e-05,-1.1635e-09,3.96456e-12,-13189,25.0572], Tmin=(100,'K'), Tmax=(1012.54,'K')), NASAPolynomial(coeffs=[12.2184,0.0150723,-5.68005e-06,1.04547e-09,-7.41146e-14,-16047.2,-29.1232], Tmin=(1012.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=O)C(O)C=C(8924)',
    structure = SMILES('C=CC(O)C(=C)[O]'),
    E0 = (-165.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,322.947,322.955,322.963],'cm^-1')),
        HinderedRotor(inertia=(0.00161626,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161556,'amu*angstrom^2'), symmetry=1, barrier=(11.9572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161555,'amu*angstrom^2'), symmetry=1, barrier=(11.9573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865519,0.0602501,-4.21973e-05,7.89108e-09,2.32314e-12,-19825.1,28.7098], Tmin=(100,'K'), Tmax=(1013.66,'K')), NASAPolynomial(coeffs=[15.3531,0.0198897,-7.34596e-06,1.32888e-09,-9.30576e-14,-23625.8,-45.6403], Tmin=(1013.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=[C]C(O)C(C)=O(8925)',
    structure = SMILES('C=[C]C(O)C(C)=O'),
    E0 = (-87.3665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,334.732,337.072],'cm^-1')),
        HinderedRotor(inertia=(0.00149913,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0015049,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115552,'amu*angstrom^2'), symmetry=1, barrier=(9.27096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116151,'amu*angstrom^2'), symmetry=1, barrier=(9.27694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26888,0.0574432,-4.53322e-05,1.85772e-08,-3.1259e-12,-10407.5,28.1243], Tmin=(100,'K'), Tmax=(1384.49,'K')), NASAPolynomial(coeffs=[12.2349,0.0257608,-1.10067e-05,2.04874e-09,-1.41343e-13,-13444,-28.3484], Tmin=(1384.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.3665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(O)C(C)=O(8926)',
    structure = SMILES('[CH]=CC(O)C(C)=O'),
    E0 = (-78.1121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3120,650,792.5,1650,407.914],'cm^-1')),
        HinderedRotor(inertia=(0.00101314,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853833,'amu*angstrom^2'), symmetry=1, barrier=(10.0818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853837,'amu*angstrom^2'), symmetry=1, barrier=(10.0818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0853826,'amu*angstrom^2'), symmetry=1, barrier=(10.0818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.95366,0.0594598,-4.74076e-05,1.91507e-08,-3.0979e-12,-9278.85,29.1948], Tmin=(100,'K'), Tmax=(1469.96,'K')), NASAPolynomial(coeffs=[14.9257,0.0214398,-8.61068e-06,1.55528e-09,-1.05409e-13,-13386.5,-43.5945], Tmin=(1469.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.1121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'CH3OH(20)(21)',
    structure = SMILES('CO'),
    E0 = (-211.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1964.49,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0117392,'amu*angstrom^2'), symmetry=1, barrier=(6.09573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.0419,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.93,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.65851,-0.0162983,6.91938e-05,-7.58373e-08,2.80428e-11,-25612,-0.897331], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.52727,0.0103179,-3.62893e-06,5.77448e-10,-3.42183e-14,-26002.9,5.16759], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-211.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH3OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C=CC=C=O(6774)',
    structure = SMILES('C=CC=C=O'),
    E0 = (28.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.764513,'amu*angstrom^2'), symmetry=1, barrier=(17.5777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.70095,0.0290788,-1.80138e-05,5.65892e-09,-7.4916e-13,3463.32,14.9666], Tmin=(100,'K'), Tmax=(1629.64,'K')), NASAPolynomial(coeffs=[7.32083,0.0177393,-7.57635e-06,1.38912e-09,-9.41402e-14,1957.57,-9.57784], Tmin=(1629.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C([O])C(O)C=C(8927)',
    structure = SMILES('[CH2]C([O])C(O)C=C'),
    E0 = (56.3707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,316.199,316.209,316.212],'cm^-1')),
        HinderedRotor(inertia=(0.00168613,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00168628,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22323,'amu*angstrom^2'), symmetry=1, barrier=(15.8365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223213,'amu*angstrom^2'), symmetry=1, barrier=(15.8365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323032,0.071439,-6.43345e-05,2.93836e-08,-5.29416e-12,6920.11,31.5942], Tmin=(100,'K'), Tmax=(1346.61,'K')), NASAPolynomial(coeffs=[17.4302,0.0206232,-7.73003e-06,1.36026e-09,-9.1549e-14,2312.81,-56.0288], Tmin=(1346.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.3707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=CC([O])C(C)[O](8928)',
    structure = SMILES('C=CC([O])C(C)[O]'),
    E0 = (75.1424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,402.069,402.076,404.545,405.077],'cm^-1')),
        HinderedRotor(inertia=(0.0010323,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134038,'amu*angstrom^2'), symmetry=1, barrier=(15.1704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130164,'amu*angstrom^2'), symmetry=1, barrier=(15.1649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770968,0.0607713,-3.38901e-05,-4.14574e-10,4.59143e-12,9162.63,29.1558], Tmin=(100,'K'), Tmax=(1049.08,'K')), NASAPolynomial(coeffs=[15.3456,0.0245424,-9.74534e-06,1.81682e-09,-1.28467e-13,5040.25,-46.9299], Tmin=(1049.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.1424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]CC([O])C(C)=O(8929)',
    structure = SMILES('[CH2]CC([O])C(C)=O'),
    E0 = (-2.87631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,180,932.631,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00407799,'amu*angstrom^2'), symmetry=1, barrier=(2.51699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109492,'amu*angstrom^2'), symmetry=1, barrier=(2.51744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109481,'amu*angstrom^2'), symmetry=1, barrier=(2.51719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81145,'amu*angstrom^2'), symmetry=1, barrier=(18.6568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0384,0.0580543,-4.08833e-05,1.46549e-08,-2.13303e-12,-233.578,31.2511], Tmin=(100,'K'), Tmax=(1596.92,'K')), NASAPolynomial(coeffs=[14.1299,0.0252624,-1.00816e-05,1.79609e-09,-1.19972e-13,-4414.81,-38.0358], Tmin=(1596.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.87631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(C=OCOJ) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]C(O)C(C)[O](8930)',
    structure = SMILES('C=[C]C(O)C(C)[O]'),
    E0 = (82.6234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,474.168,474.184,474.188],'cm^-1')),
        HinderedRotor(inertia=(0.0667971,'amu*angstrom^2'), symmetry=1, barrier=(10.6583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0667979,'amu*angstrom^2'), symmetry=1, barrier=(10.6582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0667972,'amu*angstrom^2'), symmetry=1, barrier=(10.6584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0667959,'amu*angstrom^2'), symmetry=1, barrier=(10.6583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.697401,0.0679222,-5.9482e-05,2.71022e-08,-4.97994e-12,10060.1,30.0851], Tmin=(100,'K'), Tmax=(1299.39,'K')), NASAPolynomial(coeffs=[14.5232,0.0253612,-1.03503e-05,1.89464e-09,-1.30057e-13,6467.03,-40.2374], Tmin=(1299.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.6234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[CH]C([O])C(C)=O(8931)',
    structure = SMILES('C[CH]C([O])C(C)=O'),
    E0 = (-8.22062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,406.772,406.777,4000],'cm^-1')),
        HinderedRotor(inertia=(6.98687e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.793294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198944,'amu*angstrom^2'), symmetry=1, barrier=(23.3601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198942,'amu*angstrom^2'), symmetry=1, barrier=(23.3601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00675472,'amu*angstrom^2'), symmetry=1, barrier=(0.79313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23265,0.0518316,-2.16368e-05,-6.11638e-09,5.21456e-12,-881.397,31.2538], Tmin=(100,'K'), Tmax=(1072.27,'K')), NASAPolynomial(coeffs=[12.4592,0.0270233,-1.08137e-05,2.0024e-09,-1.40122e-13,-4270.38,-28.2676], Tmin=(1072.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.22062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)C([CH2])=O(8932)',
    structure = SMILES('[CH2]CC(O)C(=C)[O]'),
    E0 = (-84.6161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,460.044,460.103,460.198],'cm^-1')),
        HinderedRotor(inertia=(0.000796783,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0525677,'amu*angstrom^2'), symmetry=1, barrier=(7.89236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052493,'amu*angstrom^2'), symmetry=1, barrier=(7.89294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0525425,'amu*angstrom^2'), symmetry=1, barrier=(7.89147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594418,0.066502,-5.13748e-05,1.49001e-08,3.16303e-13,-10047,31.1576], Tmin=(100,'K'), Tmax=(1012.21,'K')), NASAPolynomial(coeffs=[15.9084,0.0215414,-7.80052e-06,1.38434e-09,-9.55542e-14,-13944.1,-46.8461], Tmin=(1012.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.6161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]=CC(O)C(C)[O](8933)',
    structure = SMILES('[CH]=CC(O)C(C)[O]'),
    E0 = (91.8778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,337.35,337.351],'cm^-1')),
        HinderedRotor(inertia=(0.00148125,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164985,'amu*angstrom^2'), symmetry=1, barrier=(13.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164984,'amu*angstrom^2'), symmetry=1, barrier=(13.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164983,'amu*angstrom^2'), symmetry=1, barrier=(13.324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404337,0.0697194,-6.09578e-05,2.70887e-08,-4.76601e-12,11187.6,31.073], Tmin=(100,'K'), Tmax=(1373.76,'K')), NASAPolynomial(coeffs=[17.0363,0.0212916,-8.07954e-06,1.42755e-09,-9.61147e-14,6617.96,-54.4483], Tmin=(1373.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.8778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'CO(10)(11)',
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
    label = 'C=CC(C)O(8934)',
    structure = SMILES('C=CC(C)O'),
    E0 = (-189.328,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,227.484],'cm^-1')),
        HinderedRotor(inertia=(0.00325763,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368666,'amu*angstrom^2'), symmetry=1, barrier=(13.5382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368633,'amu*angstrom^2'), symmetry=1, barrier=(13.5381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.1057,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94486,0.0349451,7.66343e-06,-3.36245e-08,1.49837e-11,-22687.9,20.7218], Tmin=(100,'K'), Tmax=(989.597,'K')), NASAPolynomial(coeffs=[11.3699,0.0209937,-7.78749e-06,1.43962e-09,-1.03046e-13,-25735.5,-30.6234], Tmin=(989.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=C(C)OO(4826)',
    structure = SMILES('C=CC=C(C)OO'),
    E0 = (-32.1681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.68108,'amu*angstrom^2'), symmetry=1, barrier=(15.6594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681472,'amu*angstrom^2'), symmetry=1, barrier=(15.6684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681081,'amu*angstrom^2'), symmetry=1, barrier=(15.6594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68167,'amu*angstrom^2'), symmetry=1, barrier=(15.6729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677001,0.0672743,-5.70167e-05,2.49782e-08,-4.40826e-12,-3744.55,25.832], Tmin=(100,'K'), Tmax=(1350.2,'K')), NASAPolynomial(coeffs=[14.9071,0.0251176,-1.01834e-05,1.85433e-09,-1.26733e-13,-7587.3,-47.093], Tmin=(1350.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.1681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=COC(C)=CO(8935)',
    structure = SMILES('C=COC(C)=CO'),
    E0 = (-280.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.859513,'amu*angstrom^2'), symmetry=1, barrier=(19.7619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860787,'amu*angstrom^2'), symmetry=1, barrier=(19.7912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860305,'amu*angstrom^2'), symmetry=1, barrier=(19.7801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858128,'amu*angstrom^2'), symmetry=1, barrier=(19.7301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349842,0.0641618,-2.02372e-05,-3.49558e-08,2.25901e-11,-33622.5,26.0046], Tmin=(100,'K'), Tmax=(925.816,'K')), NASAPolynomial(coeffs=[22.2988,0.0109765,-1.54004e-06,1.67265e-10,-1.41725e-14,-39471.5,-87.8335], Tmin=(925.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-280.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(O)C(=C)O(8936)',
    structure = SMILES('C=CC(O)C(=C)O'),
    E0 = (-303.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,310.711,310.711],'cm^-1')),
        HinderedRotor(inertia=(0.217041,'amu*angstrom^2'), symmetry=1, barrier=(14.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217041,'amu*angstrom^2'), symmetry=1, barrier=(14.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217043,'amu*angstrom^2'), symmetry=1, barrier=(14.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217041,'amu*angstrom^2'), symmetry=1, barrier=(14.8691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530509,0.0638616,-3.3443e-05,-1.0758e-08,1.09951e-11,-36383.6,28.4108], Tmin=(100,'K'), Tmax=(966.149,'K')), NASAPolynomial(coeffs=[18.8665,0.016907,-5.50377e-06,9.87169e-10,-7.17454e-14,-41278.2,-66.4129], Tmin=(966.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C[C]C(O)C(C)=O(8937)',
    structure = SMILES('C[C]C(O)C(C)=O'),
    E0 = (-23.9904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06667,0.069392,-7.15337e-05,4.60184e-08,-1.29709e-11,-2784.04,36.3083], Tmin=(100,'K'), Tmax=(835.344,'K')), NASAPolynomial(coeffs=[7.2629,0.0397216,-1.82554e-05,3.49824e-09,-2.45489e-13,-3819.23,7.52985], Tmin=(835.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.9904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC(O)C(C)=O(8938)',
    structure = SMILES('[CH]CC(O)C(C)=O'),
    E0 = (-0.755936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11555,0.0625446,-5.01717e-05,2.10022e-08,-3.63937e-12,13.4123,28.5346], Tmin=(100,'K'), Tmax=(1336.94,'K')), NASAPolynomial(coeffs=[12.1966,0.029391,-1.29743e-05,2.4537e-09,-1.70894e-13,-2949.53,-28.143], Tmin=(1336.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.755936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'C[C]1OC[CH]C1O(8939)',
    structure = SMILES('C[C]1OC[CH]C1O'),
    E0 = (-26.8497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63152,0.0429715,8.43969e-06,-4.51968e-08,2.31302e-11,-3135.29,22.7227], Tmin=(100,'K'), Tmax=(872.194,'K')), NASAPolynomial(coeffs=[11.3323,0.0249844,-6.20406e-06,8.34061e-10,-4.9469e-14,-5835.51,-28.5301], Tmin=(872.194,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.8497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOCs) + radical(CCJCO)"""),
)

species(
    label = 'C[C]1CC=C1O(8940)',
    structure = SMILES('C[C]1CC=C1O'),
    E0 = (50.7307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98811,0.0211655,7.55084e-05,-1.2211e-07,5.08511e-11,6195.2,18.7733], Tmin=(100,'K'), Tmax=(933.439,'K')), NASAPolynomial(coeffs=[18.4417,0.00965599,-8.03795e-07,1.04724e-10,-1.63241e-14,553.238,-73.2402], Tmin=(933.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.7307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CCJ(C)CO)"""),
)

species(
    label = 'CC1([O])CC=C1O(8941)',
    structure = SMILES('CC1([O])CC=C1O'),
    E0 = (-62.6152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20116,0.042172,2.83038e-05,-7.73804e-08,3.56549e-11,-7412.29,23.7194], Tmin=(100,'K'), Tmax=(936.973,'K')), NASAPolynomial(coeffs=[19.7686,0.0120777,-2.23705e-06,3.59015e-10,-3.14129e-14,-13050.2,-76.1673], Tmin=(936.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.6152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'CC1(O)[C]=CC1(8942)',
    structure = SMILES('CC1(O)[C]=CC1'),
    E0 = (175.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80291,0.034165,2.34317e-05,-5.9057e-08,2.64381e-11,21178.9,20.7615], Tmin=(100,'K'), Tmax=(942.566,'K')), NASAPolynomial(coeffs=[14.4388,0.0161346,-4.51783e-06,7.74338e-10,-5.72451e-14,17215.8,-47.8393], Tmin=(942.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'CC1(O)CC=C1[O](8943)',
    structure = SMILES('CC1(O)C[CH]C1=O'),
    E0 = (-157.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40762,0.0375501,3.82811e-05,-8.37446e-08,3.67188e-11,-18842.6,22.497], Tmin=(100,'K'), Tmax=(945.848,'K')), NASAPolynomial(coeffs=[18.2446,0.0152662,-3.96003e-06,7.09974e-10,-5.67252e-14,-24215.9,-69.3617], Tmin=(945.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-157.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'CC1(O)[CH]C=C1O(8944)',
    structure = SMILES('CC1(O)[CH]C=C1O'),
    E0 = (-174.801,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,3150,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617964,0.0522113,1.81717e-05,-8.02961e-08,4.01966e-11,-20881.4,21.7987], Tmin=(100,'K'), Tmax=(920.595,'K')), NASAPolynomial(coeffs=[24.9084,0.00439381,2.02878e-06,-4.9356e-10,2.82717e-14,-27799.8,-106.664], Tmin=(920.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-174.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1(O)CC=C1O(8945)',
    structure = SMILES('[CH2]C1(O)CC=C1O'),
    E0 = (-78.2744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71788,0.0541692,5.50949e-07,-5.62018e-08,3.05029e-11,-9279.31,25.1163], Tmin=(100,'K'), Tmax=(918.15,'K')), NASAPolynomial(coeffs=[22.4151,0.00687297,6.60039e-07,-2.55393e-10,1.44153e-14,-15254.3,-88.5483], Tmin=(918.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.2744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'CC1(O)C[C]=C1O(8946)',
    structure = SMILES('CC1(O)C[C]=C1O'),
    E0 = (-39.2321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80837,0.0528011,1.28881e-06,-5.49954e-08,2.96325e-11,-4587.44,24.3897], Tmin=(100,'K'), Tmax=(918.473,'K')), NASAPolynomial(coeffs=[21.5671,0.00795054,1.38801e-07,-1.5978e-10,8.16221e-15,-10322.2,-84.4541], Tmin=(918.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.2321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'H2O(7)(7)',
    structure = SMILES('O'),
    E0 = (-251.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1591.57,3573.8,3573.81],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(6727.26,'J/mol'), sigma=(2.641,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.99882,-0.000554831,2.76774e-06,-1.55665e-09,3.0233e-13,-30274.6,-0.0308942], Tmin=(100,'K'), Tmax=(1281.44,'K')), NASAPolynomial(coeffs=[3.1956,0.0019524,-1.67118e-07,-2.97938e-11,4.45139e-15,-30068.7,4.04335], Tmin=(1281.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CC1C=CC=1O(6799)',
    structure = SMILES('CC1C=CC=1O'),
    E0 = (172.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,3150,900,1100,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27474,0.0480772,-1.42986e-05,-2.56106e-08,1.61739e-11,20800.7,16.3889], Tmin=(100,'K'), Tmax=(933.508,'K')), NASAPolynomial(coeffs=[17.0647,0.010512,-2.29256e-06,3.48393e-10,-2.62771e-14,16541.5,-65.7253], Tmin=(933.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(cyclobutadiene_13)"""),
)

species(
    label = 'C=C1CC=C1O(8947)',
    structure = SMILES('C=C1CC=C1O'),
    E0 = (5.83132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2950,3150,900,1000,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58642,0.0323418,4.28909e-05,-9.35281e-08,4.25332e-11,807.514,15.7267], Tmin=(100,'K'), Tmax=(920.693,'K')), NASAPolynomial(coeffs=[20.5985,0.00356641,2.08245e-06,-4.83818e-10,2.73233e-14,-4974.58,-86.8132], Tmin=(920.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.83132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'CC1(O)C[CH]C1[O](8948)',
    structure = SMILES('CC1(O)C[CH]C1[O]'),
    E0 = (57.4312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18334,0.0442299,2.33567e-05,-6.77241e-08,3.06876e-11,7024.89,25.6204], Tmin=(100,'K'), Tmax=(952.476,'K')), NASAPolynomial(coeffs=[17.7508,0.0186305,-5.58509e-06,1.00806e-09,-7.61831e-14,1874.05,-63.9731], Tmin=(952.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.4312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1(O)[CH]C[C]1O(8949)',
    structure = SMILES('CC1(O)[CH]C[C]1O'),
    E0 = (3.6982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865323,0.0569746,-1.96833e-05,-2.03113e-08,1.33894e-11,568.437,27.071], Tmin=(100,'K'), Tmax=(963.687,'K')), NASAPolynomial(coeffs=[16.4763,0.0204367,-6.79645e-06,1.20217e-09,-8.5426e-14,-3752.56,-54.4734], Tmin=(963.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.6982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(C2CsJOH)"""),
)

species(
    label = 'CC1([O])C[CH]C1O(8950)',
    structure = SMILES('CC1([O])C[CH]C1O'),
    E0 = (55.9863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2339,0.0443007,1.94604e-05,-6.0549e-08,2.72515e-11,6848.12,25.6924], Tmin=(100,'K'), Tmax=(963.525,'K')), NASAPolynomial(coeffs=[16.635,0.0208102,-6.93601e-06,1.28114e-09,-9.51502e-14,2002.77,-57.7797], Tmin=(963.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.9863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = 'CC1(O)[CH][CH]C1O(8951)',
    structure = SMILES('CC1(O)[CH][CH]C1O'),
    E0 = (26.9724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21666,0.0429913,2.61417e-05,-7.2265e-08,3.29667e-11,3360.92,28.1892], Tmin=(100,'K'), Tmax=(942.843,'K')), NASAPolynomial(coeffs=[18.3331,0.0162066,-4.16035e-06,7.17616e-10,-5.55695e-14,-1903.81,-64.1833], Tmin=(942.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.9724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C1(O)C[CH]C1O(8952)',
    structure = SMILES('[CH2]C1(O)C[CH]C1O'),
    E0 = (39.0865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729922,0.0585644,-1.84114e-05,-2.53662e-08,1.60311e-11,4830.86,27.2784], Tmin=(100,'K'), Tmax=(956.128,'K')), NASAPolynomial(coeffs=[18.0526,0.0182878,-5.73056e-06,1.0078e-09,-7.28978e-14,46.7886,-63.2116], Tmin=(956.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.0865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CCJCO)"""),
)

species(
    label = 'CC1([O])CC[C]1O(8953)',
    structure = SMILES('CC1([O])CC[C]1O'),
    E0 = (32.7121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862376,0.058503,-2.70218e-05,-7.9399e-09,7.49874e-12,4056.51,24.6477], Tmin=(100,'K'), Tmax=(1014.27,'K')), NASAPolynomial(coeffs=[14.9737,0.0247156,-9.388e-06,1.72274e-09,-1.21478e-13,69.3765,-49.1751], Tmin=(1014.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.7121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)2OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1(O)CC[C]1O(8954)',
    structure = SMILES('[CH2]C1(O)CC[C]1O'),
    E0 = (15.8123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.20329,0.0745427,-7.08267e-05,3.44907e-08,-6.60524e-12,2046.03,26.7934], Tmin=(100,'K'), Tmax=(1274.51,'K')), NASAPolynomial(coeffs=[17.4165,0.0205196,-7.2455e-06,1.23282e-09,-8.15536e-14,-2341.64,-60.4256], Tmin=(1274.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.8123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(C2CsJOH)"""),
)

species(
    label = 'OC1=CCC1O(8955)',
    structure = SMILES('OC1=CCC1O'),
    E0 = (-242.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70002,0.0298393,4.56077e-05,-9.38361e-08,4.18318e-11,-29111.9,20.1871], Tmin=(100,'K'), Tmax=(928.648,'K')), NASAPolynomial(coeffs=[20.1069,0.00386597,1.45045e-06,-3.18036e-10,1.40431e-14,-34829.3,-79.6294], Tmin=(928.648,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-242.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=C(O)C([CH2])(C)O(8956)',
    structure = SMILES('[CH]=C(O)C([CH2])(C)O'),
    E0 = (7.83241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.576215,0.0870349,-9.60557e-05,5.16868e-08,-1.05082e-11,1118.45,31.6321], Tmin=(100,'K'), Tmax=(1333.18,'K')), NASAPolynomial(coeffs=[22.0148,0.0113067,-1.90974e-06,1.36956e-10,-3.03866e-15,-4198.86,-81.2042], Tmin=(1333.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.83241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(Cds_P)"""),
)

species(
    label = 'C[C](O)CC=[C]O(8957)',
    structure = SMILES('C[C](O)CC=[C]O'),
    E0 = (-5.19306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2750,2850,1437.5,1250,1305,750,350,325.182,325.183],'cm^-1')),
        HinderedRotor(inertia=(0.147955,'amu*angstrom^2'), symmetry=1, barrier=(11.1022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147955,'amu*angstrom^2'), symmetry=1, barrier=(11.1022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147955,'amu*angstrom^2'), symmetry=1, barrier=(11.1022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147955,'amu*angstrom^2'), symmetry=1, barrier=(11.1022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147955,'amu*angstrom^2'), symmetry=1, barrier=(11.1022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.245745,0.0794806,-8.81305e-05,5.10007e-08,-1.16232e-11,-486.708,30.5679], Tmin=(100,'K'), Tmax=(1075.5,'K')), NASAPolynomial(coeffs=[15.6728,0.0221043,-8.10769e-06,1.39707e-09,-9.27992e-14,-3805.05,-44.9813], Tmin=(1075.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.19306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C2CsJOH) + radical(C=CJO)"""),
)

species(
    label = 'CC1(O)CCC1=O(8958)',
    structure = SMILES('CC1(O)CCC1=O'),
    E0 = (-319.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29839,0.0448333,1.27489e-05,-5.07676e-08,2.32795e-11,-38342.9,22.5038], Tmin=(100,'K'), Tmax=(964.706,'K')), NASAPolynomial(coeffs=[15.2243,0.0225177,-7.63645e-06,1.38558e-09,-1.00423e-13,-42678.2,-52.7244], Tmin=(964.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-319.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = 'CC1(O)C[C]C1O(8959)',
    structure = SMILES('CC1(O)C[C]C1O'),
    E0 = (55.8322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.701585,0.0620931,-3.18409e-05,-5.27136e-09,7.10484e-12,6843.02,34.3436], Tmin=(100,'K'), Tmax=(1001,'K')), NASAPolynomial(coeffs=[15.4844,0.0250327,-9.29069e-06,1.6783e-09,-1.174e-13,2780.71,-42.4983], Tmin=(1001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.8322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(CsJ2_singlet-CsH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=CC(=O)[CH]C(8960)',
    structure = SMILES('C=CC([O])=CC'),
    E0 = (-18.1925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.679014,'amu*angstrom^2'), symmetry=1, barrier=(15.6119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.681937,'amu*angstrom^2'), symmetry=1, barrier=(15.6791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.301,0.0519692,-3.14497e-05,5.18841e-10,4.67222e-12,-2084.25,20.2318], Tmin=(100,'K'), Tmax=(967.226,'K')), NASAPolynomial(coeffs=[13.2178,0.0193052,-6.56546e-06,1.13077e-09,-7.72981e-14,-5166.83,-40.8811], Tmin=(967.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.1925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=O)C(C)[O](8961)',
    structure = SMILES('C=CC(=O)C(C)[O]'),
    E0 = (-78.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0148917,'amu*angstrom^2'), symmetry=1, barrier=(6.10245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199027,'amu*angstrom^2'), symmetry=1, barrier=(17.5029,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.760944,'amu*angstrom^2'), symmetry=1, barrier=(17.4956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54456,0.0516492,-3.29251e-05,1.002e-08,-1.22776e-12,-9298.6,27.1531], Tmin=(100,'K'), Tmax=(1834.64,'K')), NASAPolynomial(coeffs=[13.9292,0.0246474,-1.08483e-05,1.99781e-09,-1.34599e-13,-13842.9,-40.1112], Tmin=(1834.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC(=O)[CH]O(8962)',
    structure = SMILES('C=CC([O])=CO'),
    E0 = (-190.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10149,'amu*angstrom^2'), symmetry=1, barrier=(25.3253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09938,'amu*angstrom^2'), symmetry=1, barrier=(25.277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.629641,0.0698866,-7.57188e-05,3.72983e-08,-6.55468e-12,-22772.7,25.2448], Tmin=(100,'K'), Tmax=(1727.13,'K')), NASAPolynomial(coeffs=[19.7192,0.000412813,4.02616e-06,-9.7421e-10,6.95552e-14,-26468.8,-74.3978], Tmin=(1727.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-190.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C2H5O(53)(52)',
    structure = SMILES('C[CH]O'),
    E0 = (-69.1109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.000744213,'amu*angstrom^2'), symmetry=1, barrier=(8.44984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226097,'amu*angstrom^2'), symmetry=1, barrier=(24.137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0605,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.22283,0.00512175,3.48387e-05,-4.91944e-08,2.01184e-11,-8356.22,8.01676], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.35842,0.0124356,-4.33097e-06,6.8453e-10,-4.03713e-14,-9530.19,-6.05106], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-69.1109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CH3CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH2CHCO(3566)',
    structure = SMILES('C=C[C]=O'),
    E0 = (83.3963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.27992,'amu*angstrom^2'), symmetry=1, barrier=(29.4278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08334,0.0153204,6.54286e-06,-1.77538e-08,7.39385e-12,10067.5,11.895], Tmin=(100,'K'), Tmax=(1002.99,'K')), NASAPolynomial(coeffs=[6.97837,0.0111885,-4.32951e-06,8.06755e-10,-5.75269e-14,8712.68,-9.76679], Tmin=(1002.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(O)C(=O)C=C(8963)',
    structure = SMILES('[CH2]C(O)C(=O)C=C'),
    E0 = (-110.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,641.153,641.159],'cm^-1')),
        HinderedRotor(inertia=(0.00756797,'amu*angstrom^2'), symmetry=1, barrier=(2.20759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0960008,'amu*angstrom^2'), symmetry=1, barrier=(2.20725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575002,'amu*angstrom^2'), symmetry=1, barrier=(13.2204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575003,'amu*angstrom^2'), symmetry=1, barrier=(13.2204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24388,0.0587792,-4.59841e-05,1.81134e-08,-2.91556e-12,-13154,28.0796], Tmin=(100,'K'), Tmax=(1438.26,'K')), NASAPolynomial(coeffs=[13.051,0.0259419,-1.17371e-05,2.23905e-09,-1.56261e-13,-16550.3,-33.1742], Tmin=(1438.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(CJCO)"""),
)

species(
    label = 'CC(O)[C]=O(8964)',
    structure = SMILES('CC(O)[C]=O'),
    E0 = (-213.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00163274,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146615,'amu*angstrom^2'), symmetry=1, barrier=(10.6834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146646,'amu*angstrom^2'), symmetry=1, barrier=(10.7122,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26172,0.0301261,-1.75819e-06,-2.18079e-08,1.12305e-11,-25577.5,21.2628], Tmin=(100,'K'), Tmax=(963.868,'K')), NASAPolynomial(coeffs=[11.3006,0.0120966,-4.01738e-06,7.23927e-10,-5.24577e-14,-28224.9,-26.7063], Tmin=(963.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C(=O)C(C)O(8965)',
    structure = SMILES('C=C=C([O])C(C)O'),
    E0 = (-125.504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,293.612,293.612],'cm^-1')),
        HinderedRotor(inertia=(0.194663,'amu*angstrom^2'), symmetry=1, barrier=(11.9085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194663,'amu*angstrom^2'), symmetry=1, barrier=(11.9085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194663,'amu*angstrom^2'), symmetry=1, barrier=(11.9085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.664821,0.0653836,-5.37331e-05,1.74633e-08,-4.04936e-13,-14967.5,27.7334], Tmin=(100,'K'), Tmax=(1002.57,'K')), NASAPolynomial(coeffs=[16.0869,0.0189106,-6.73027e-06,1.18849e-09,-8.2114e-14,-18816.6,-50.4828], Tmin=(1002.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC(=O)C(C)O(8966)',
    structure = SMILES('[CH]=CC(=O)C(C)O'),
    E0 = (-74.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3120,650,792.5,1650,518.383],'cm^-1')),
        HinderedRotor(inertia=(0.191477,'amu*angstrom^2'), symmetry=1, barrier=(4.40244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84526,'amu*angstrom^2'), symmetry=1, barrier=(19.4342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00171167,'amu*angstrom^2'), symmetry=1, barrier=(19.4344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845455,'amu*angstrom^2'), symmetry=1, barrier=(19.4387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32634,0.0570614,-4.26569e-05,1.59072e-08,-2.4258e-12,-8886.62,27.5528], Tmin=(100,'K'), Tmax=(1504.39,'K')), NASAPolynomial(coeffs=[12.9506,0.0261539,-1.18398e-05,2.25072e-09,-1.56375e-13,-12384.1,-33.2748], Tmin=(1504.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(=O)C=C(6849)',
    structure = SMILES('C=CC(=O)C=C'),
    E0 = (-23.1554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.439253,'amu*angstrom^2'), symmetry=1, barrier=(10.0993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114097,'amu*angstrom^2'), symmetry=1, barrier=(20.0775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.40946,0.0377612,-1.86611e-05,3.53817e-09,-1.87846e-13,-2730.31,18.3203], Tmin=(100,'K'), Tmax=(2103.34,'K')), NASAPolynomial(coeffs=[16.8425,0.0162284,-7.52334e-06,1.345e-09,-8.60854e-14,-10110.2,-65.1523], Tmin=(2103.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.1554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=[C]C([O])C(C)O(8967)',
    structure = SMILES('C=[C]C([O])C(C)O'),
    E0 = (82.6234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,474.168,474.184,474.188],'cm^-1')),
        HinderedRotor(inertia=(0.0667971,'amu*angstrom^2'), symmetry=1, barrier=(10.6583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0667979,'amu*angstrom^2'), symmetry=1, barrier=(10.6582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0667972,'amu*angstrom^2'), symmetry=1, barrier=(10.6584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0667959,'amu*angstrom^2'), symmetry=1, barrier=(10.6583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.697401,0.0679222,-5.9482e-05,2.71022e-08,-4.97994e-12,10060.1,30.0851], Tmin=(100,'K'), Tmax=(1299.39,'K')), NASAPolynomial(coeffs=[14.5232,0.0253612,-1.03503e-05,1.89464e-09,-1.30057e-13,6467.03,-40.2374], Tmin=(1299.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.6234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C(O)C([O])C=C(8968)',
    structure = SMILES('[CH2]C(O)C([O])C=C'),
    E0 = (56.3707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,316.199,316.209,316.212],'cm^-1')),
        HinderedRotor(inertia=(0.00168613,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00168628,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22323,'amu*angstrom^2'), symmetry=1, barrier=(15.8365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223213,'amu*angstrom^2'), symmetry=1, barrier=(15.8365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.323032,0.071439,-6.43345e-05,2.93836e-08,-5.29416e-12,6920.11,31.5942], Tmin=(100,'K'), Tmax=(1346.61,'K')), NASAPolynomial(coeffs=[17.4302,0.0206232,-7.73003e-06,1.36026e-09,-9.1549e-14,2312.81,-56.0288], Tmin=(1346.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.3707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=CC([O])C(C)O(8969)',
    structure = SMILES('[CH]=CC([O])C(C)O'),
    E0 = (91.8778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,337.35,337.351],'cm^-1')),
        HinderedRotor(inertia=(0.00148125,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164985,'amu*angstrom^2'), symmetry=1, barrier=(13.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164984,'amu*angstrom^2'), symmetry=1, barrier=(13.324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164983,'amu*angstrom^2'), symmetry=1, barrier=(13.324,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404337,0.0697194,-6.09578e-05,2.70887e-08,-4.76601e-12,11187.6,31.073], Tmin=(100,'K'), Tmax=(1373.76,'K')), NASAPolynomial(coeffs=[17.0363,0.0212916,-8.07954e-06,1.42755e-09,-9.61147e-14,6617.96,-54.4483], Tmin=(1373.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.8778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC(=O)C(C)[O](8970)',
    structure = SMILES('[CH2]CC(=O)C(C)[O]'),
    E0 = (1.82978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,253.256,253.258,4000],'cm^-1')),
        HinderedRotor(inertia=(0.550332,'amu*angstrom^2'), symmetry=1, barrier=(25.0453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000991364,'amu*angstrom^2'), symmetry=1, barrier=(11.256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0112888,'amu*angstrom^2'), symmetry=1, barrier=(8.80378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.550306,'amu*angstrom^2'), symmetry=1, barrier=(25.0455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15704,0.0648368,-5.84199e-05,2.96528e-08,-6.3845e-12,320.515,28.752], Tmin=(100,'K'), Tmax=(1085.18,'K')), NASAPolynomial(coeffs=[9.74089,0.0331967,-1.46854e-05,2.78524e-09,-1.94876e-13,-1542.5,-13.3619], Tmin=(1085.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.82978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(CJCC=O) + radical(C=OCOJ)"""),
)

species(
    label = 'C[CH]C(=O)C(C)[O](8971)',
    structure = SMILES('CC=C([O])C(C)[O]'),
    E0 = (-71.7471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,398.959,399.32,399.613],'cm^-1')),
        HinderedRotor(inertia=(0.0740137,'amu*angstrom^2'), symmetry=1, barrier=(8.39472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0741932,'amu*angstrom^2'), symmetry=1, barrier=(8.4121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0736857,'amu*angstrom^2'), symmetry=1, barrier=(8.38907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7656,0.0620058,-3.96494e-05,5.06881e-09,3.09157e-12,-8504.78,28.1269], Tmin=(100,'K'), Tmax=(1023.92,'K')), NASAPolynomial(coeffs=[15.0517,0.023587,-8.8441e-06,1.59946e-09,-1.11411e-13,-12342,-45.5848], Tmin=(1023.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.7471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=O)CO(8972)',
    structure = SMILES('C=CC(=O)CO'),
    E0 = (-271.747,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,461.431,462.172],'cm^-1')),
        HinderedRotor(inertia=(0.0594077,'amu*angstrom^2'), symmetry=1, barrier=(8.52855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0557732,'amu*angstrom^2'), symmetry=1, barrier=(8.47025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0606532,'amu*angstrom^2'), symmetry=1, barrier=(8.57595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17634,0.0439245,-3.05575e-05,1.0846e-08,-1.65296e-12,-32621.7,19.9132], Tmin=(100,'K'), Tmax=(1415.25,'K')), NASAPolynomial(coeffs=[8.11206,0.027148,-1.27764e-05,2.47004e-09,-1.73372e-13,-34301.8,-10.7848], Tmin=(1415.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-271.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=CC)OO(8973)',
    structure = SMILES('C=CC(=CC)OO'),
    E0 = (-31.5024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.683899,'amu*angstrom^2'), symmetry=1, barrier=(15.7242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.684139,'amu*angstrom^2'), symmetry=1, barrier=(15.7297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683793,'amu*angstrom^2'), symmetry=1, barrier=(15.7217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683817,'amu*angstrom^2'), symmetry=1, barrier=(15.7223,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921293,0.0679381,-6.11397e-05,2.97657e-08,-6.00315e-12,-3678.24,24.1399], Tmin=(100,'K'), Tmax=(1168.95,'K')), NASAPolynomial(coeffs=[11.8029,0.0307028,-1.33596e-05,2.51632e-09,-1.75439e-13,-6222.27,-30.0564], Tmin=(1168.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.5024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(=CO)OC(8974)',
    structure = SMILES('C=CC(=CO)OC'),
    E0 = (-310.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.966386,'amu*angstrom^2'), symmetry=1, barrier=(22.2191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96698,'amu*angstrom^2'), symmetry=1, barrier=(22.2328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968696,'amu*angstrom^2'), symmetry=1, barrier=(22.2722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970924,'amu*angstrom^2'), symmetry=1, barrier=(22.3235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0381864,0.0676991,-1.6784e-05,-4.90185e-08,3.04164e-11,-37133.2,22.9535], Tmin=(100,'K'), Tmax=(911.896,'K')), NASAPolynomial(coeffs=[26.3376,0.00421302,2.31441e-06,-5.97094e-10,3.86644e-14,-44086.6,-113.327], Tmin=(911.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-310.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C(O)C(C)O(8975)',
    structure = SMILES('C=C=C(O)C(C)O'),
    E0 = (-263.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.707698,'amu*angstrom^2'), symmetry=1, barrier=(16.2714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706982,'amu*angstrom^2'), symmetry=1, barrier=(16.2549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708,'amu*angstrom^2'), symmetry=1, barrier=(16.2783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707759,'amu*angstrom^2'), symmetry=1, barrier=(16.2728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.33177,0.0689642,-4.48301e-05,-1.44485e-09,8.41019e-12,-31526,27.4279], Tmin=(100,'K'), Tmax=(956.809,'K')), NASAPolynomial(coeffs=[19.6269,0.0158855,-4.8648e-06,8.41474e-10,-6.03743e-14,-36481.1,-71.4067], Tmin=(956.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-263.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C[C]C(=O)C(C)O(8976)',
    structure = SMILES('C[C]C(=O)C(C)O'),
    E0 = (-23.0382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08259,0.0639423,-4.98741e-05,2.06494e-08,-3.58759e-12,-2665.9,37.7734], Tmin=(100,'K'), Tmax=(1320.49,'K')), NASAPolynomial(coeffs=[11.4258,0.0326111,-1.4284e-05,2.68141e-09,-1.8585e-13,-5397.55,-15.0022], Tmin=(1320.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.0382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]CC(=O)C(C)O(8977)',
    structure = SMILES('[CH]CC(=O)C(C)O'),
    E0 = (-0.602226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,375,552.5,462.5,1710,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0889,0.0629364,-5.04251e-05,2.10124e-08,-3.61873e-12,33.0103,28.5947], Tmin=(100,'K'), Tmax=(1345.88,'K')), NASAPolynomial(coeffs=[12.4029,0.0293107,-1.29489e-05,2.44898e-09,-1.70533e-13,-3012.45,-29.3497], Tmin=(1345.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.602226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (79.8717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (48.9099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-59.8722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-19.4004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-17.3313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (64.9136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (59.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (100.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (24.6912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (24.8964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (71.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-27.1186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (3.89844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (113.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-99.7195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-135.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-150.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-76.4621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-5.34099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-104.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (86.8112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (53.6573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-154.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-1.38046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-230.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (177.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (75.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-20.2932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (75.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-20.2932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (49.4191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-40.0423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (176.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (79.7437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (88.9981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (24.8986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (25.4569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-94.6401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (-27.1186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (72.9184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-26.5603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (82.7311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-109.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (-25.0814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (-150.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (-76.4621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (-51.8511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-77.2676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (85.9973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-206.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-165.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (57.0465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (67.2371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (-247.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (179.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (33.1952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (179.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (33.1952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (57.2038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (48.4831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (101.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (100.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (65.9604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (65.9604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (-11.5419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (35.7862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (76.0405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (81.2208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (28.7124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (52.3386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (135.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (-18.5143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (-21.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (-114.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (-173.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (-173.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (110.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (-20.9589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (74.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (-20.9589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (47.1674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (-40.0423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (108.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (44.4283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (88.3324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (-35.1913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (44.2709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (-53.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (51.7519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (-75.2132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (101.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (-118.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (-34.7102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (-113.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (-92.5282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (-118.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (-97.8726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (91.0976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (-288.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (-230.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (-27.8778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (-14.4718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (-248.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (9.51383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS105',
    E0 = (130.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS106',
    E0 = (-47.2594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS107',
    E0 = (9.89984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS108',
    E0 = (-20.9589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS109',
    E0 = (25.8106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS110',
    E0 = (45.9549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS111',
    E0 = (124.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS112',
    E0 = (133.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS113',
    E0 = (145.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS114',
    E0 = (-35.7518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS115',
    E0 = (-116.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS116',
    E0 = (79.2321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS117',
    E0 = (138.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS118',
    E0 = (60.5238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS119',
    E0 = (171.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS120',
    E0 = (55.9275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS121',
    E0 = (16.7526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS122',
    E0 = (89.5585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS123',
    E0 = (-59.6429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS124',
    E0 = (116.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS125',
    E0 = (-64.9872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS126',
    E0 = (72.6628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS127',
    E0 = (128.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS128',
    E0 = (49.8801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS129',
    E0 = (-82.8463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS130',
    E0 = (-99.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS131',
    E0 = (63.4488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS132',
    E0 = (36.2821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS133',
    E0 = (123.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS134',
    E0 = (79.6215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS135',
    E0 = (149.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS136',
    E0 = (203.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS137',
    E0 = (54.2021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS138',
    E0 = (51.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS139',
    E0 = (42.7548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS140',
    E0 = (133.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS141',
    E0 = (172.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS142',
    E0 = (152.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS143',
    E0 = (-42.9661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS144',
    E0 = (80.2925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS145',
    E0 = (26.8441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS146',
    E0 = (119.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS147',
    E0 = (90.3726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS148',
    E0 = (102.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS149',
    E0 = (67.2301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS150',
    E0 = (50.3303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS151',
    E0 = (176.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS152',
    E0 = (16.1167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS153',
    E0 = (3.09126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS154',
    E0 = (-103.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS155',
    E0 = (-148.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS156',
    E0 = (143.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS157',
    E0 = (13.1415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS158',
    E0 = (133.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS159',
    E0 = (-54.7678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS160',
    E0 = (14.2854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS161',
    E0 = (-20.4401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS162',
    E0 = (101.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS163',
    E0 = (73.1164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS164',
    E0 = (90.4851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS165',
    E0 = (137.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS166',
    E0 = (-62.7479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS167',
    E0 = (112.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS168',
    E0 = (43.7104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS169',
    E0 = (105.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS170',
    E0 = (138.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS171',
    E0 = (119.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS172',
    E0 = (-75.2132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS173',
    E0 = (155.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS174',
    E0 = (-40.2707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS175',
    E0 = (-59.0424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS176',
    E0 = (-135.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS177',
    E0 = (-34.7102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS178',
    E0 = (26.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS179',
    E0 = (-5.34099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS180',
    E0 = (-46.7738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS181',
    E0 = (-65.5455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS182',
    E0 = (148.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS183',
    E0 = (282.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS184',
    E0 = (3.71981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS185',
    E0 = (-120.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS186',
    E0 = (33.8102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS187',
    E0 = (36.4358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction74',
    reactants = ['OH(5)(5)', 'C=C=C([O])CC(8829)'],
    products = ['S(161)(160)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.117e+08,'m^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/OneDe] for rate rule [O_pri_rad;Cd_rad/OneDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction75',
    reactants = ['H(3)(3)', '[CH2]C(=O)C(=O)CC(8842)'],
    products = ['S(161)(160)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction76',
    reactants = ['CH3(15)(16)', 'C=C([O])C(=C)O(8849)'],
    products = ['S(161)(160)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.54203e+08,'m^3/(mol*s)'), n=-0.441, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;C_methyl] for rate rule [C_rad/H2/CO;C_methyl]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction77',
    reactants = ['C2H5(27)(28)', '[CH2]C(O)=C=O(8828)'],
    products = ['S(161)(160)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.12e+14,'cm^3/(mol*s)','*|/',3), n=-0.5, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [CO_sec_rad;C_rad/H2/Cs] for rate rule [CO_rad/OneDe;C_rad/H2/Cs]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction78',
    reactants = ['H(3)(3)', 'C=C(O)C([O])=CC(8848)'],
    products = ['S(161)(160)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H/OneDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction79',
    reactants = ['H(3)(3)', '[CH2]CC(=O)C(=C)O(8855)'],
    products = ['S(161)(160)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction80',
    reactants = ['C3H5O(279)(278)', 'CH2COH(1431)'],
    products = ['S(161)(160)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.81e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;CO_rad] for rate rule [Cd_rad/NonDe;CO_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Recombination"""),
)

reaction(
    label = 'reaction81',
    reactants = ['H(3)(3)', '[CH]=C(O)C(=O)CC(8856)'],
    products = ['S(161)(160)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction82',
    reactants = ['[CH2]C([O])C(=O)CC(8843)'],
    products = ['S(161)(160)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction83',
    reactants = ['C=C(O)C([O])[CH]C(8850)'],
    products = ['S(161)(160)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction84',
    reactants = ['[CH2]CC([O])C(=C)O(8851)'],
    products = ['S(161)(160)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction85',
    reactants = ['[CH2]C(O)C(=O)[CH]C(8844)'],
    products = ['S(161)(160)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH2]C(=O)C([O])CC(8770)'],
    products = ['S(161)(160)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH]=C(O)C([O])CC(8852)'],
    products = ['S(161)(160)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction88',
    reactants = ['[CH2]CC(O)=C([CH2])O(8832)'],
    products = ['S(161)(160)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction89',
    reactants = ['C[CH]C([O])=C(C)O(8833)'],
    products = ['S(161)(160)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction90',
    reactants = ['[CH2]C([O])=C(O)CC(8834)'],
    products = ['S(161)(160)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction91',
    reactants = ['[CH]=C(O)[C](O)CC(8853)'],
    products = ['S(161)(160)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction92',
    reactants = ['[CH2]CC(=O)C([CH2])O(8845)'],
    products = ['S(161)(160)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction93',
    reactants = ['[CH2]CC([O])=C(C)O(8835)'],
    products = ['S(161)(160)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction94',
    reactants = ['CH2(S)(21)(22)', 'C=C(O)C(C)=O(8857)'],
    products = ['S(161)(160)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction99',
    reactants = ['[CH]C(O)C(=O)CC(8861)'],
    products = ['S(161)(160)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH]
Euclidian distance = 1.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction95',
    reactants = ['S(161)(160)'],
    products = ['CCC(=O)C(C)=O(8858)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction96',
    reactants = ['C=C(O)C(=C)OC(8859)'],
    products = ['S(161)(160)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction97',
    reactants = ['S(168)(167)'],
    products = ['S(161)(160)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction98',
    reactants = ['C=C=C(CC)OO(8860)'],
    products = ['S(161)(160)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction100',
    reactants = ['OH(5)(5)', 'C=[C]C(O)=CC(8862)'],
    products = ['S(168)(167)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/Cd] for rate rule [O_pri_rad;Cd_rad/Cd]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction101',
    reactants = ['H(3)(3)', 'C=C([O])C(O)=CC(8863)'],
    products = ['S(168)(167)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction102',
    reactants = ['OH(5)(5)', 'C=C(O)[C]=CC(8864)'],
    products = ['S(168)(167)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/Cd;Y_rad] for rate rule [Cd_rad/Cd;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction103',
    reactants = ['H(3)(3)', 'C=C(O)C([O])=CC(8848)'],
    products = ['S(168)(167)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction104',
    reactants = ['CH3(15)(16)', '[CH]=C(O)C(=C)O(8865)'],
    products = ['S(168)(167)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.81e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 69 used for C_methyl;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction105',
    reactants = ['H(3)(3)', '[CH2]C=C(O)C(=C)O(8866)'],
    products = ['S(168)(167)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.156e+12,'cm^3/(mol*s)'), n=0.461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 15 used for C_rad/H2/Cd;H_rad
Exact match found for rate rule [C_rad/H2/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction106',
    reactants = ['CH2COH(1431)', 'CC=[C]O(8867)'],
    products = ['S(168)(167)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(556926,'m^3/(mol*s)'), n=0.4, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;Cd_rad] for rate rule [Cd_rad/NonDe;Cd_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Recombination
Ea raised from -2.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction107',
    reactants = ['H(3)(3)', 'C=C(O)C(O)=[C]C(8868)'],
    products = ['S(168)(167)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction108',
    reactants = ['H(3)(3)', '[CH]=C(O)C(O)=CC(8869)'],
    products = ['S(168)(167)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction109',
    reactants = ['[CH2]C([O])C(O)=CC(8870)'],
    products = ['S(168)(167)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction110',
    reactants = ['C=C(O)C([O])[CH]C(8850)'],
    products = ['S(168)(167)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction111',
    reactants = ['[CH2]CC(O)=C([CH2])O(8832)'],
    products = ['S(168)(167)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction112',
    reactants = ['[CH2]C(O)C(=O)[CH]C(8844)'],
    products = ['S(168)(167)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction113',
    reactants = ['[CH2]C(O)C(O)=[C]C(8871)'],
    products = ['S(168)(167)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction114',
    reactants = ['C=C([O])C(O)[CH]C(8872)'],
    products = ['S(168)(167)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction115',
    reactants = ['[CH]=C(O)C(O)[CH]C(8873)'],
    products = ['S(168)(167)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction116',
    reactants = ['C[CH]C([O])=C(C)O(8833)'],
    products = ['S(168)(167)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction117',
    reactants = ['C[C]=C(O)[C](C)O(8874)'],
    products = ['S(168)(167)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction118',
    reactants = ['[CH2]C([O])=C(O)CC(8834)'],
    products = ['S(168)(167)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction119',
    reactants = ['[CH]=C(O)[C](O)CC(8853)'],
    products = ['S(168)(167)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction120',
    reactants = ['[CH2]C=C(O)C([CH2])O(8875)'],
    products = ['S(168)(167)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction121',
    reactants = ['[CH2]C=C(O)[C](C)O(8876)'],
    products = ['S(168)(167)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction122',
    reactants = ['CH2(S)(21)(22)', 'C=C(O)C(=C)O(8877)'],
    products = ['S(168)(167)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction123',
    reactants = ['S(168)(167)'],
    products = ['S(172)(171)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CdH2_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction124',
    reactants = ['S(168)(167)'],
    products = ['CC=C(O)C(C)=O(8878)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction125',
    reactants = ['C=C(O)C(O)[C]C(8879)'],
    products = ['S(168)(167)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2C;CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction126',
    reactants = ['[CH]C(O)C(O)=CC(8880)'],
    products = ['S(168)(167)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH]
Euclidian distance = 1.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction127',
    reactants = ['S(168)(167)'],
    products = ['S(169)(168)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction128',
    reactants = ['OH(5)(5)', 'CC1C[C]=C1O(8881)'],
    products = ['S(169)(168)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [O_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction129',
    reactants = ['H(3)(3)', 'CC1CC([O])=C1O(8882)'],
    products = ['S(169)(168)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction130',
    reactants = ['OH(5)(5)', 'CC1[C]=C(O)C1(8883)'],
    products = ['S(169)(168)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction131',
    reactants = ['H(3)(3)', 'CC1CC(O)=C1[O](8884)'],
    products = ['S(169)(168)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction132',
    reactants = ['CH3(15)(16)', 'OC1[CH]CC=1O(8885)'],
    products = ['S(169)(168)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/H/CdCs;Y_rad] for rate rule [C_rad/H/CdCs;C_methyl]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction133',
    reactants = ['H(3)(3)', 'C[C]1CC(O)=C1O(8886)'],
    products = ['S(169)(168)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction134',
    reactants = ['H(3)(3)', 'CC1[CH]C(O)=C1O(8887)'],
    products = ['S(169)(168)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction135',
    reactants = ['H(3)(3)', '[CH2]C1CC(O)=C1O(8888)'],
    products = ['S(169)(168)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction136',
    reactants = ['CC1CC([O])[C]1O(8889)'],
    products = ['S(169)(168)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction137',
    reactants = ['CC1C[C](O)C1[O](8890)'],
    products = ['S(169)(168)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction138',
    reactants = ['C[C]1C[C](O)C1O(8891)'],
    products = ['S(169)(168)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction139',
    reactants = ['CC1[CH]C(O)[C]1O(8892)'],
    products = ['S(169)(168)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(2.5515e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction140',
    reactants = ['CC1[CH][C](O)C1O(8893)'],
    products = ['S(169)(168)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction141',
    reactants = ['[CH2]C1C[C](O)C1O(8894)'],
    products = ['S(169)(168)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction142',
    reactants = ['C[C]1CC(O)[C]1O(8895)'],
    products = ['S(169)(168)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction143',
    reactants = ['[CH2]C1CC(O)[C]1O(8896)'],
    products = ['S(169)(168)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction144',
    reactants = ['CH2(S)(21)(22)', 'OC1CCC=1O(8897)'],
    products = ['S(169)(168)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.74695e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/OneDeC]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction145',
    reactants = ['[CH2]C(C)C(O)=[C]O(8898)'],
    products = ['S(169)(168)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriND_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction146',
    reactants = ['C[CH]CC(O)=[C]O(8899)'],
    products = ['S(169)(168)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriND_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction147',
    reactants = ['[CH2]C(O)=C(O)[CH]C(8830)'],
    products = ['S(169)(168)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Cpri_rad_out_2H] + [R4;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction148',
    reactants = ['S(169)(168)'],
    products = ['CC1CC(O)C1=O(8900)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction149',
    reactants = ['S(169)(168)'],
    products = ['CC1CC(=O)C1O(8901)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction150',
    reactants = ['OH(5)(5)', 'C=CC(O)=[C]C(8902)'],
    products = ['S(172)(171)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [O_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction151',
    reactants = ['H(3)(3)', 'C=CC(O)=C(C)[O](8903)'],
    products = ['S(172)(171)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction152',
    reactants = ['OH(5)(5)', 'C=C[C]=C(C)O(8904)'],
    products = ['S(172)(171)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/Cd;Y_rad] for rate rule [Cd_rad/Cd;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction153',
    reactants = ['H(3)(3)', 'C=CC([O])=C(C)O(8905)'],
    products = ['S(172)(171)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction154',
    reactants = ['CH3(15)(16)', 'C=CC(O)=[C]O(8906)'],
    products = ['S(172)(171)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(1.87377e+06,'m^3/(mol*s)'), n=0.2675, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_rad;C_methyl] + [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;C_methyl]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction155',
    reactants = ['H(3)(3)', '[CH2]C=C(O)C(=C)O(8866)'],
    products = ['S(172)(171)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(3.156e+12,'cm^3/(mol*s)'), n=0.461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 15 used for C_rad/H2/Cd;H_rad
Exact match found for rate rule [C_rad/H2/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction156',
    reactants = ['C2H3(28)(29)', 'CC(O)=[C]O(8907)'],
    products = ['S(172)(171)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(2.68887e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_rad;Cd_pri_rad] + [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction157',
    reactants = ['H(3)(3)', 'C=[C]C(O)=C(C)O(8908)'],
    products = ['S(172)(171)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [Cd_rad/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction158',
    reactants = ['H(3)(3)', '[CH]=CC(O)=C(C)O(8909)'],
    products = ['S(172)(171)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction159',
    reactants = ['C=C[C](O)C(C)[O](8910)'],
    products = ['S(172)(171)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction160',
    reactants = ['C=CC([O])[C](C)O(8911)'],
    products = ['S(172)(171)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction161',
    reactants = ['[CH2]C=C(O)C([CH2])O(8875)'],
    products = ['S(172)(171)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction162',
    reactants = ['C=[C]C(O)[C](C)O(8912)'],
    products = ['S(172)(171)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction163',
    reactants = ['[CH2]CC([O])=C(C)O(8835)'],
    products = ['S(172)(171)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction164',
    reactants = ['[CH]=CC(O)[C](C)O(8913)'],
    products = ['S(172)(171)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction165',
    reactants = ['C[CH]C([O])=C(C)O(8833)'],
    products = ['S(172)(171)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction166',
    reactants = ['[CH]=C[C](O)C(C)O(8914)'],
    products = ['S(172)(171)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction167',
    reactants = ['[CH2]CC(O)=C(C)[O](8915)'],
    products = ['S(172)(171)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction168',
    reactants = ['[CH2]CC(O)=C([CH2])O(8832)'],
    products = ['S(172)(171)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction169',
    reactants = ['C[CH]C(O)=C(C)[O](8916)'],
    products = ['S(172)(171)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction170',
    reactants = ['[CH2]C(O)=C(O)[CH]C(8830)'],
    products = ['S(172)(171)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction171',
    reactants = ['CH2(S)(21)(22)', 'C=CC(O)=CO(8917)'],
    products = ['S(172)(171)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(884144,'m^3/(mol*s)'), n=0.112, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cd_sec] for rate rule [carbene;Cd/H/NonDeO]
Euclidian distance = 1.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction172',
    reactants = ['S(172)(171)'],
    products = ['S(173)(172)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction173',
    reactants = ['S(172)(171)'],
    products = ['S(176)(175)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction174',
    reactants = ['C[C]C(O)=C(C)O(8918)'],
    products = ['S(172)(171)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(7.00341e+13,'s^-1'), n=-1.27142, Ea=(26.2576,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(C=C);CH] for rate rule [CsJ2-C;CsJ2(C=C);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction175',
    reactants = ['[CH]CC(O)=C(C)O(8919)'],
    products = ['S(172)(171)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(6.84965e+11,'s^-1'), n=0.4135, Ea=(17.2276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [singletcarbene_CH;singletcarbene;CH2(C=C)] for rate rule [CsJ2-C;CsJ2H;CH2(C=C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction176',
    reactants = ['S(172)(171)'],
    products = ['S(175)(174)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;C=C_1;CdH2_2]
Euclidian distance = 1.0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction177',
    reactants = ['OH(5)(5)', 'C=C[CH]C(C)=O(8920)'],
    products = ['S(173)(172)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad;C_sec_rad] for rate rule [O_pri_rad;C_rad/H/TwoDe]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction178',
    reactants = ['H(3)(3)', 'C=CC([O])C(C)=O(8921)'],
    products = ['S(173)(172)'],
    transitionState = 'TS105',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction179',
    reactants = ['CH3CO(55)(54)', 'hydroxyl1allyl(7170)'],
    products = ['S(173)(172)'],
    transitionState = 'TS106',
    kinetics = Arrhenius(A=(6.64e+13,'cm^3/(mol*s)','*|/',2), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_sec_rad;CO_rad/NonDe] for rate rule [C_rad/H/OneDeO;CO_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction180',
    reactants = ['C2H3(28)(29)', 'CC(=O)[CH]O(8922)'],
    products = ['S(173)(172)'],
    transitionState = 'TS107',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction181',
    reactants = ['H(3)(3)', 'C=CC(O)=C(C)[O](8903)'],
    products = ['S(173)(172)'],
    transitionState = 'TS108',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction182',
    reactants = ['CH3(15)(16)', 'C=CC(O)[C]=O(8923)'],
    products = ['S(173)(172)'],
    transitionState = 'TS109',
    kinetics = Arrhenius(A=(4.2e+13,'cm^3/(mol*s)','+|-',8.4e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 72 used for C_methyl;CO_rad/NonDe
Exact match found for rate rule [CO_rad/NonDe;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction183',
    reactants = ['H(3)(3)', '[CH2]C(=O)C(O)C=C(8924)'],
    products = ['S(173)(172)'],
    transitionState = 'TS110',
    kinetics = Arrhenius(A=(8.83648e+06,'m^3/(mol*s)'), n=0.3308, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;H_rad] for rate rule [C_rad/H2/CO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction184',
    reactants = ['H(3)(3)', 'C=[C]C(O)C(C)=O(8925)'],
    products = ['S(173)(172)'],
    transitionState = 'TS111',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction185',
    reactants = ['H(3)(3)', '[CH]=CC(O)C(C)=O(8926)'],
    products = ['S(173)(172)'],
    transitionState = 'TS112',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction186',
    reactants = ['CH3OH(20)(21)', 'C=CC=C=O(6774)'],
    products = ['S(173)(172)'],
    transitionState = 'TS113',
    kinetics = Arrhenius(A=(1.79e-05,'cm^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [doublebond;CH3OH] for rate rule [Cdd_Cd_HDe;CH3OH]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction187',
    reactants = ['C=C[C](O)C(C)[O](8910)'],
    products = ['S(173)(172)'],
    transitionState = 'TS114',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction188',
    reactants = ['[CH2]CC(O)=C(C)[O](8915)'],
    products = ['S(173)(172)'],
    transitionState = 'TS115',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction189',
    reactants = ['[CH2]C([O])C(O)C=C(8927)'],
    products = ['S(173)(172)'],
    transitionState = 'TS116',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction190',
    reactants = ['C=CC([O])C(C)[O](8928)'],
    products = ['S(173)(172)'],
    transitionState = 'TS117',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction191',
    reactants = ['[CH2]CC([O])C(C)=O(8929)'],
    products = ['S(173)(172)'],
    transitionState = 'TS118',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction192',
    reactants = ['C=[C]C(O)C(C)[O](8930)'],
    products = ['S(173)(172)'],
    transitionState = 'TS119',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction193',
    reactants = ['C=CC([O])[C](C)O(8911)'],
    products = ['S(173)(172)'],
    transitionState = 'TS120',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction194',
    reactants = ['C[CH]C([O])C(C)=O(8931)'],
    products = ['S(173)(172)'],
    transitionState = 'TS121',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction195',
    reactants = ['C=[C]C(O)[C](C)O(8912)'],
    products = ['S(173)(172)'],
    transitionState = 'TS122',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction196',
    reactants = ['[CH2]CC(O)C([CH2])=O(8932)'],
    products = ['S(173)(172)'],
    transitionState = 'TS123',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction197',
    reactants = ['[CH]=CC(O)C(C)[O](8933)'],
    products = ['S(173)(172)'],
    transitionState = 'TS124',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction198',
    reactants = ['C=C([O])C(O)[CH]C(8872)'],
    products = ['S(173)(172)'],
    transitionState = 'TS125',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction199',
    reactants = ['[CH]=CC(O)[C](C)O(8913)'],
    products = ['S(173)(172)'],
    transitionState = 'TS126',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction200',
    reactants = ['CO(10)(11)', 'C=CC(C)O(8934)'],
    products = ['S(173)(172)'],
    transitionState = 'TS127',
    kinetics = Arrhenius(A=(538,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [CO;Cs_Cs]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction201',
    reactants = ['C=CC=C(C)OO(4826)'],
    products = ['S(173)(172)'],
    transitionState = 'TS128',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CH3;R_O_R]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction202',
    reactants = ['C=COC(C)=CO(8935)'],
    products = ['S(173)(172)'],
    transitionState = 'TS129',
    kinetics = Arrhenius(A=(855.663,'s^-1'), n=2.935, Ea=(197.924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CH3;R_O_C]
Euclidian distance = 2.2360679775
family: ketoenol"""),
)

reaction(
    label = 'reaction203',
    reactants = ['C=CC(O)C(=C)O(8936)'],
    products = ['S(173)(172)'],
    transitionState = 'TS130',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction204',
    reactants = ['C[C]C(O)C(C)=O(8937)'],
    products = ['S(173)(172)'],
    transitionState = 'TS131',
    kinetics = Arrhenius(A=(4.85495e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [CsJ2-C;CsJ2(CsC);CH] for rate rule [CsJ2-C;CsJ2(CsC);CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction205',
    reactants = ['[CH]CC(O)C(C)=O(8938)'],
    products = ['S(173)(172)'],
    transitionState = 'TS132',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction206',
    reactants = ['C[C]1OC[CH]C1O(8939)'],
    products = ['S(173)(172)'],
    transitionState = 'TS133',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction207',
    reactants = ['OH(5)(5)', 'C[C]1CC=C1O(8940)'],
    products = ['S(175)(174)'],
    transitionState = 'TS134',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDeCs2;Y_rad] for rate rule [C_rad/OneDeCs2;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction208',
    reactants = ['H(3)(3)', 'CC1([O])CC=C1O(8941)'],
    products = ['S(175)(174)'],
    transitionState = 'TS135',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [O_sec_rad;H_rad] + [O_rad/NonDe;Y_rad] for rate rule [O_rad/NonDe;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction209',
    reactants = ['OH(5)(5)', 'CC1(O)[C]=CC1(8942)'],
    products = ['S(175)(174)'],
    transitionState = 'TS136',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/NonDe] for rate rule [O_pri_rad;Cd_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction210',
    reactants = ['H(3)(3)', 'CC1(O)CC=C1[O](8943)'],
    products = ['S(175)(174)'],
    transitionState = 'TS137',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction211',
    reactants = ['CH3(15)(16)', 'OC1[CH]CC=1O(8885)'],
    products = ['S(175)(174)'],
    transitionState = 'TS138',
    kinetics = Arrhenius(A=(4.88e+15,'cm^3/(mol*s)','*|/',2), n=-1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [C_ter_rad;C_methyl] for rate rule [C_rad/OneDeO;C_methyl]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction212',
    reactants = ['H(3)(3)', 'CC1(O)[CH]C=C1O(8944)'],
    products = ['S(175)(174)'],
    transitionState = 'TS139',
    kinetics = Arrhenius(A=(2.71464e+07,'m^3/(mol*s)'), n=0.107721, Ea=(5.76381,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 36 used for C_rad/H/CdCs;H_rad
Exact match found for rate rule [C_rad/H/CdCs;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction213',
    reactants = ['H(3)(3)', '[CH2]C1(O)CC=C1O(8945)'],
    products = ['S(175)(174)'],
    transitionState = 'TS140',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction214',
    reactants = ['H(3)(3)', 'CC1(O)C[C]=C1O(8946)'],
    products = ['S(175)(174)'],
    transitionState = 'TS141',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction215',
    reactants = ['H2O(7)(7)', 'CC1C=CC=1O(6799)'],
    products = ['S(175)(174)'],
    transitionState = 'TS142',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/H/De_Cd/Nd/De;H_OH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction216',
    reactants = ['H2O(7)(7)', 'C=C1CC=C1O(8947)'],
    products = ['S(175)(174)'],
    transitionState = 'TS143',
    kinetics = Arrhenius(A=(2500,'cm^3/(mol*s)'), n=2.76, Ea=(202.924,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd/unsub_Cd/disub;H_OH] for rate rule [Cd/H2_Cd/Nd/De;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction217',
    reactants = ['CC1(O)C[CH]C1[O](8948)'],
    products = ['S(175)(174)'],
    transitionState = 'TS144',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction218',
    reactants = ['CC1(O)[CH]C[C]1O(8949)'],
    products = ['S(175)(174)'],
    transitionState = 'TS145',
    kinetics = Arrhenius(A=(5.10299e+10,'s^-1'), n=0.2847, Ea=(23.1459,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad_NDe] + [R2radExo;Y_rad;XH_Rrad_NDe] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction219',
    reactants = ['CC1([O])C[CH]C1O(8950)'],
    products = ['S(175)(174)'],
    transitionState = 'TS146',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction220',
    reactants = ['CC1(O)[CH][CH]C1O(8951)'],
    products = ['S(175)(174)'],
    transitionState = 'TS147',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction221',
    reactants = ['[CH2]C1(O)C[CH]C1O(8952)'],
    products = ['S(175)(174)'],
    transitionState = 'TS148',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction222',
    reactants = ['CC1([O])CC[C]1O(8953)'],
    products = ['S(175)(174)'],
    transitionState = 'TS149',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction223',
    reactants = ['[CH2]C1(O)CC[C]1O(8954)'],
    products = ['S(175)(174)'],
    transitionState = 'TS150',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction224',
    reactants = ['CH2(S)(21)(22)', 'OC1=CCC1O(8955)'],
    products = ['S(175)(174)'],
    transitionState = 'TS151',
    kinetics = Arrhenius(A=(436738,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H/ODMustO]
Euclidian distance = 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction225',
    reactants = ['[CH]=C(O)C([CH2])(C)O(8956)'],
    products = ['S(175)(174)'],
    transitionState = 'TS152',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction226',
    reactants = ['C[C](O)CC=[C]O(8957)'],
    products = ['S(175)(174)'],
    transitionState = 'TS153',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_NDMustO;CdsinglepriND_rad_out]
Euclidian distance = 3.74165738677
family: Birad_recombination"""),
)

reaction(
    label = 'reaction227',
    reactants = ['[CH2]C=C(O)[C](C)O(8876)'],
    products = ['S(175)(174)'],
    transitionState = 'TS154',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SDS;C_rad_out_NDMustO;Cpri_rad_out_2H]
Euclidian distance = 3.16227766017
family: Birad_recombination"""),
)

reaction(
    label = 'reaction228',
    reactants = ['S(175)(174)'],
    products = ['CC1(O)CCC1=O(8958)'],
    transitionState = 'TS155',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_H]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction229',
    reactants = ['CC1(O)C[C]C1O(8959)'],
    products = ['S(175)(174)'],
    transitionState = 'TS156',
    kinetics = Arrhenius(A=(1.61832e+16,'s^-1'), n=-0.885455, Ea=(87.4392,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2(CsC);CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction230',
    reactants = ['OH(5)(5)', 'C=CC(=O)[CH]C(8960)'],
    products = ['S(176)(175)'],
    transitionState = 'TS157',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/OneDeC] for rate rule [O_pri_rad;C_rad/H/OneDeC]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction231',
    reactants = ['H(3)(3)', 'C=CC(=O)C(C)[O](8961)'],
    products = ['S(176)(175)'],
    transitionState = 'TS158',
    kinetics = Arrhenius(A=(5.21063e+06,'m^3/(mol*s)'), n=0.156446, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [H_rad;O_sec_rad] for rate rule [H_rad;O_rad/NonDe]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction232',
    reactants = ['CH3(15)(16)', 'C=CC(=O)[CH]O(8962)'],
    products = ['S(176)(175)'],
    transitionState = 'TS159',
    kinetics = Arrhenius(A=(4.64257e+09,'m^3/(mol*s)'), n=-0.82, Ea=(0.004184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_sec_rad;C_methyl] for rate rule [C_rad/H/OneDeO;C_methyl]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction233',
    reactants = ['C2H5O(53)(52)', 'CH2CHCO(3566)'],
    products = ['S(176)(175)'],
    transitionState = 'TS160',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/NonDe;CO_sec_rad] for rate rule [C_rad/H/CsO;CO_rad/OneDe]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction234',
    reactants = ['H(3)(3)', 'C=CC([O])=C(C)O(8905)'],
    products = ['S(176)(175)'],
    transitionState = 'TS161',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction235',
    reactants = ['H(3)(3)', '[CH2]C(O)C(=O)C=C(8963)'],
    products = ['S(176)(175)'],
    transitionState = 'TS162',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction236',
    reactants = ['C2H3(28)(29)', 'CC(O)[C]=O(8964)'],
    products = ['S(176)(175)'],
    transitionState = 'TS163',
    kinetics = Arrhenius(A=(1.6647e+07,'m^3/(mol*s)'), n=-0.0666667, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_rad;Cd_pri_rad] + [CO_rad/NonDe;Y_rad] for rate rule [CO_rad/NonDe;Cd_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction237',
    reactants = ['H(3)(3)', 'C=[C]C(=O)C(C)O(8965)'],
    products = ['S(176)(175)'],
    transitionState = 'TS164',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction238',
    reactants = ['H(3)(3)', '[CH]=CC(=O)C(C)O(8966)'],
    products = ['S(176)(175)'],
    transitionState = 'TS165',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction239',
    reactants = ['H2O(7)(7)', 'C=CC(=O)C=C(6849)'],
    products = ['S(176)(175)'],
    transitionState = 'TS166',
    kinetics = Arrhenius(A=(260.8,'cm^3/(mol*s)'), n=2.92, Ea=(212.129,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd/unsub_Cd/monosub;H_OH] for rate rule [Cd/H2_Cd/H/De;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction240',
    reactants = ['C3H4O(121)(120)', 'C2H4O(40)(40)'],
    products = ['S(176)(175)'],
    transitionState = 'TS167',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;R_OH] for rate rule [Cdd_Cd_HNd;Cd_pri_OH]
Euclidian distance = 2.82842712475
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction241',
    reactants = ['C=CC([O])[C](C)O(8911)'],
    products = ['S(176)(175)'],
    transitionState = 'TS168',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction242',
    reactants = ['C=[C]C([O])C(C)O(8967)'],
    products = ['S(176)(175)'],
    transitionState = 'TS169',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction243',
    reactants = ['C=CC([O])C(C)[O](8928)'],
    products = ['S(176)(175)'],
    transitionState = 'TS170',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction244',
    reactants = ['[CH2]C(O)C([O])C=C(8968)'],
    products = ['S(176)(175)'],
    transitionState = 'TS171',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction245',
    reactants = ['[CH2]CC([O])=C(C)O(8835)'],
    products = ['S(176)(175)'],
    transitionState = 'TS172',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction246',
    reactants = ['[CH]=CC([O])C(C)O(8969)'],
    products = ['S(176)(175)'],
    transitionState = 'TS173',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction247',
    reactants = ['C=C[C](O)C(C)[O](8910)'],
    products = ['S(176)(175)'],
    transitionState = 'TS174',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction248',
    reactants = ['[CH2]C=C(O)C([CH2])O(8875)'],
    products = ['S(176)(175)'],
    transitionState = 'TS175',
    kinetics = Arrhenius(A=(2.59e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction249',
    reactants = ['C[CH]C([O])=C(C)O(8833)'],
    products = ['S(176)(175)'],
    transitionState = 'TS176',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction250',
    reactants = ['[CH]=C[C](O)C(C)O(8914)'],
    products = ['S(176)(175)'],
    transitionState = 'TS177',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction251',
    reactants = ['[CH2]CC(=O)C(C)[O](8970)'],
    products = ['S(176)(175)'],
    transitionState = 'TS178',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction252',
    reactants = ['[CH2]CC(=O)C([CH2])O(8845)'],
    products = ['S(176)(175)'],
    transitionState = 'TS179',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction253',
    reactants = ['C[CH]C(=O)C(C)[O](8971)'],
    products = ['S(176)(175)'],
    transitionState = 'TS180',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction254',
    reactants = ['[CH2]C(O)C(=O)[CH]C(8844)'],
    products = ['S(176)(175)'],
    transitionState = 'TS181',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction255',
    reactants = ['CH2(S)(21)(22)', 'C=CC(=O)CO(8972)'],
    products = ['S(176)(175)'],
    transitionState = 'TS182',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/OneDeO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction256',
    reactants = ['C=CC(=CC)OO(8973)'],
    products = ['S(176)(175)'],
    transitionState = 'TS183',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_R] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond;R_O_R]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction257',
    reactants = ['C=CC(=CO)OC(8974)'],
    products = ['S(176)(175)'],
    transitionState = 'TS184',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction258',
    reactants = ['C=C=C(O)C(C)O(8975)'],
    products = ['S(176)(175)'],
    transitionState = 'TS185',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction259',
    reactants = ['C[C]C(=O)C(C)O(8976)'],
    products = ['S(176)(175)'],
    transitionState = 'TS186',
    kinetics = Arrhenius(A=(1.84394e+15,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction260',
    reactants = ['[CH]CC(=O)C(C)O(8977)'],
    products = ['S(176)(175)'],
    transitionState = 'TS187',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;singletcarbene;CH2(C)] + [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH2(C)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = '831',
    isomers = [
        'S(161)(160)',
        'S(168)(167)',
        'S(169)(168)',
        'S(172)(171)',
        'S(173)(172)',
        'S(175)(174)',
        'S(176)(175)',
    ],
    reactants = [
        ('C3H4O(121)(120)', 'C2H4O(40)(40)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '831',
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

