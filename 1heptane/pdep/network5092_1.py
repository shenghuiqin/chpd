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
    label = 'C=[C]C(O)=CC(31155)',
    structure = SMILES('C=[C]C(O)=CC'),
    E0 = (42.9982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02093,'amu*angstrom^2'), symmetry=1, barrier=(23.4731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01944,'amu*angstrom^2'), symmetry=1, barrier=(23.4389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0196,'amu*angstrom^2'), symmetry=1, barrier=(23.4427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502882,0.0666463,-6.75e-05,3.49508e-08,-6.96103e-12,5306.16,21.1922], Tmin=(100,'K'), Tmax=(1342.26,'K')), NASAPolynomial(coeffs=[16.3641,0.0144229,-3.60075e-06,4.62785e-10,-2.5161e-14,1494.65,-58.3348], Tmin=(1342.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
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
    label = 'C=C([O])C(O)=CC(31156)',
    structure = SMILES('C=C([O])C(O)=CC'),
    E0 = (-232.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.759458,'amu*angstrom^2'), symmetry=1, barrier=(17.4614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756943,'amu*angstrom^2'), symmetry=1, barrier=(17.4036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756852,'amu*angstrom^2'), symmetry=1, barrier=(17.4015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145462,0.0794714,-8.63323e-05,4.63573e-08,-9.46453e-12,-27754.2,24.6653], Tmin=(100,'K'), Tmax=(1313.67,'K')), NASAPolynomial(coeffs=[19.8697,0.0123364,-2.60611e-06,2.80219e-10,-1.30687e-14,-32478.6,-75.3236], Tmin=(1313.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)[C]=CC(31157)',
    structure = SMILES('C=C(O)[C]=CC'),
    E0 = (42.9982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.02093,'amu*angstrom^2'), symmetry=1, barrier=(23.4731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01944,'amu*angstrom^2'), symmetry=1, barrier=(23.4389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0196,'amu*angstrom^2'), symmetry=1, barrier=(23.4427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502882,0.0666463,-6.75e-05,3.49508e-08,-6.96103e-12,5306.16,21.1922], Tmin=(100,'K'), Tmax=(1342.26,'K')), NASAPolynomial(coeffs=[16.3641,0.0144229,-3.60075e-06,4.62785e-10,-2.5161e-14,1494.65,-58.3348], Tmin=(1342.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.9982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C(O)C([O])=CC(4627)',
    structure = SMILES('C=C(O)C([O])=CC'),
    E0 = (-232.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.757823,'amu*angstrom^2'), symmetry=1, barrier=(17.4238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756881,'amu*angstrom^2'), symmetry=1, barrier=(17.4022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759059,'amu*angstrom^2'), symmetry=1, barrier=(17.4522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145471,0.0794715,-8.63325e-05,4.63576e-08,-9.46463e-12,-27754.2,24.6653], Tmin=(100,'K'), Tmax=(1313.65,'K')), NASAPolynomial(coeffs=[19.8698,0.0123363,-2.60601e-06,2.80197e-10,-1.3067e-14,-32478.7,-75.3243], Tmin=(1313.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = '[CH]=C(O)C(=C)O(8478)',
    structure = SMILES('[CH]=C(O)C(=C)O'),
    E0 = (-86.7684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(1.02953,'amu*angstrom^2'), symmetry=1, barrier=(23.671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02396,'amu*angstrom^2'), symmetry=1, barrier=(23.543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02862,'amu*angstrom^2'), symmetry=1, barrier=(23.6499,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.322663,0.0759296,-9.08749e-05,4.91324e-08,-9.62861e-12,-10263.7,22.3962], Tmin=(100,'K'), Tmax=(1482.05,'K')), NASAPolynomial(coeffs=[22.3441,-0.00152328,3.98935e-06,-9.50257e-10,6.95844e-14,-15194.9,-89.8453], Tmin=(1482.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.7684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C(O)C(=C)O(31158)',
    structure = SMILES('[CH2]C=C(O)C(=C)O'),
    E0 = (-251.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(1.12181,'amu*angstrom^2'), symmetry=1, barrier=(25.7926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12332,'amu*angstrom^2'), symmetry=1, barrier=(25.8274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1216,'amu*angstrom^2'), symmetry=1, barrier=(25.7877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12186,'amu*angstrom^2'), symmetry=1, barrier=(25.7938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.04714,0.0838822,-8.70686e-05,4.2758e-08,-7.70112e-12,-30083.6,27.9774], Tmin=(100,'K'), Tmax=(1619.92,'K')), NASAPolynomial(coeffs=[22.8798,0.00666956,1.21675e-06,-4.84385e-10,3.87218e-14,-35456.7,-91.6557], Tmin=(1619.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
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
    label = 'CC=[C]O(31159)',
    structure = SMILES('CC=[C]O'),
    E0 = (72.7982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
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
    label = 'C=C(O)C(O)=[C]C(31160)',
    structure = SMILES('C=C(O)C(O)=[C]C'),
    E0 = (-132.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.846386,'amu*angstrom^2'), symmetry=1, barrier=(19.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846124,'amu*angstrom^2'), symmetry=1, barrier=(19.4541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846658,'amu*angstrom^2'), symmetry=1, barrier=(19.4663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.846444,'amu*angstrom^2'), symmetry=1, barrier=(19.4614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.524415,0.0878455,-0.000103087,5.81171e-08,-1.23276e-11,-15708.8,24.9954], Tmin=(100,'K'), Tmax=(1274.34,'K')), NASAPolynomial(coeffs=[22.1925,0.00883947,-1.02661e-06,-1.71752e-11,7.3761e-15,-20873.4,-87.6545], Tmin=(1274.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C(O)C(O)=CC(31161)',
    structure = SMILES('[CH]=C(O)C(O)=CC'),
    E0 = (-122.794,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.892967,'amu*angstrom^2'), symmetry=1, barrier=(20.5311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892677,'amu*angstrom^2'), symmetry=1, barrier=(20.5244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.892566,'amu*angstrom^2'), symmetry=1, barrier=(20.5219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.893372,'amu*angstrom^2'), symmetry=1, barrier=(20.5404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.750107,0.088929,-0.000102429,5.58194e-08,-1.13288e-11,-14584.3,25.736], Tmin=(100,'K'), Tmax=(1366.7,'K')), NASAPolynomial(coeffs=[23.6362,0.00641342,3.64202e-07,-2.87664e-10,2.56898e-14,-20209.4,-95.7252], Tmin=(1366.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])C(O)=CC(31162)',
    structure = SMILES('[CH2]C([O])C(O)=CC'),
    E0 = (2.03722,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,366.473,366.476],'cm^-1')),
        HinderedRotor(inertia=(0.141294,'amu*angstrom^2'), symmetry=1, barrier=(13.4661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141293,'amu*angstrom^2'), symmetry=1, barrier=(13.466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141295,'amu*angstrom^2'), symmetry=1, barrier=(13.4661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141293,'amu*angstrom^2'), symmetry=1, barrier=(13.4661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283774,0.0707516,-5.05058e-05,5.49106e-09,5.59804e-12,388.701,28.9881], Tmin=(100,'K'), Tmax=(967.788,'K')), NASAPolynomial(coeffs=[19.2911,0.0169688,-5.54956e-06,9.76858e-10,-6.95203e-14,-4450.66,-68.0839], Tmin=(967.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.03722,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C(O)C([O])[CH]C(4629)',
    structure = SMILES('C=C(O)C([O])[CH]C'),
    E0 = (2.59554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,334.962,334.964,334.965],'cm^-1')),
        HinderedRotor(inertia=(0.00150265,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00150216,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212724,'amu*angstrom^2'), symmetry=1, barrier=(16.9365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212695,'amu*angstrom^2'), symmetry=1, barrier=(16.9363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582894,0.0601292,-1.74707e-05,-3.07624e-08,1.89857e-11,448.888,30.2159], Tmin=(100,'K'), Tmax=(949.697,'K')), NASAPolynomial(coeffs=[20.0817,0.0148432,-4.13141e-06,7.20283e-10,-5.43182e-14,-4916.08,-71.5952], Tmin=(949.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.59554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)=C([CH2])O(4611)',
    structure = SMILES('[CH2]CC(O)=C([CH2])O'),
    E0 = (-117.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,325,375,415,465,420,450,1700,1750,424.855],'cm^-1')),
        HinderedRotor(inertia=(0.112797,'amu*angstrom^2'), symmetry=1, barrier=(14.3925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112269,'amu*angstrom^2'), symmetry=1, barrier=(14.3897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112576,'amu*angstrom^2'), symmetry=1, barrier=(14.3902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112523,'amu*angstrom^2'), symmetry=1, barrier=(14.3941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112345,'amu*angstrom^2'), symmetry=1, barrier=(14.394,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06791,0.0870529,-9.14386e-05,4.55824e-08,-8.38467e-12,-13928.7,32.6758], Tmin=(100,'K'), Tmax=(1562.56,'K')), NASAPolynomial(coeffs=[23.7823,0.00723112,7.46054e-07,-3.86281e-10,3.20893e-14,-19716.1,-91.9721], Tmin=(1562.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(O)C(=O)[CH]C(4622)',
    structure = SMILES('[CH2]C(O)C([O])=CC'),
    E0 = (-90.5188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,307.36,309.933],'cm^-1')),
        HinderedRotor(inertia=(0.00176934,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131214,'amu*angstrom^2'), symmetry=1, barrier=(8.83614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131255,'amu*angstrom^2'), symmetry=1, barrier=(8.8339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130366,'amu*angstrom^2'), symmetry=1, barrier=(8.83892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.345407,0.0723962,-6.93553e-05,3.42305e-08,-6.65407e-12,-10748.6,30.4624], Tmin=(100,'K'), Tmax=(1255.15,'K')), NASAPolynomial(coeffs=[16.6574,0.0204121,-7.2304e-06,1.23314e-09,-8.16962e-14,-14843.4,-51.9406], Tmin=(1255.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.5188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(O)C(O)=[C]C(31163)',
    structure = SMILES('[CH2]C(O)C(O)=[C]C'),
    E0 = (9.51822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200.168],'cm^-1')),
        HinderedRotor(inertia=(0.00386889,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449267,'amu*angstrom^2'), symmetry=1, barrier=(12.9706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426874,'amu*angstrom^2'), symmetry=1, barrier=(12.9173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449783,'amu*angstrom^2'), symmetry=1, barrier=(12.9095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.445007,'amu*angstrom^2'), symmetry=1, barrier=(12.9281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0392474,0.0808275,-8.62719e-05,4.61606e-08,-9.57735e-12,1297.04,30.8137], Tmin=(100,'K'), Tmax=(1187.87,'K')), NASAPolynomial(coeffs=[18.9048,0.0170356,-5.7176e-06,9.51058e-10,-6.2495e-14,-3203.56,-63.8418], Tmin=(1187.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.51822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'C=C([O])C(O)[CH]C(11048)',
    structure = SMILES('C=C([O])C(O)[CH]C'),
    E0 = (-89.9604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,300.618,300.7,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00186589,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186341,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279062,'amu*angstrom^2'), symmetry=1, barrier=(17.9225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279467,'amu*angstrom^2'), symmetry=1, barrier=(17.9245,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773408,0.0603072,-3.14832e-05,-7.79839e-09,8.95214e-12,-10694,31.2248], Tmin=(100,'K'), Tmax=(968.867,'K')), NASAPolynomial(coeffs=[16.492,0.0198534,-6.69236e-06,1.18047e-09,-8.316e-14,-14887,-50.0315], Tmin=(968.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)[CH]C(16341)',
    structure = SMILES('[CH]=C(O)C(O)[CH]C'),
    E0 = (19.3309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,350,440,435,1725,402.616],'cm^-1')),
        HinderedRotor(inertia=(0.00104002,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115404,'amu*angstrom^2'), symmetry=1, barrier=(13.2734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115396,'amu*angstrom^2'), symmetry=1, barrier=(13.2735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115397,'amu*angstrom^2'), symmetry=1, barrier=(13.2735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115395,'amu*angstrom^2'), symmetry=1, barrier=(13.2736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35222,0.0675747,-3.97791e-05,-8.68595e-09,1.16147e-11,2467.87,31.639], Tmin=(100,'K'), Tmax=(944.712,'K')), NASAPolynomial(coeffs=[20.3737,0.0138257,-3.69596e-06,6.11973e-10,-4.46799e-14,-2699.43,-71.1412], Tmin=(944.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.3309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = 'C[CH]C([O])=C(C)O(4612)',
    structure = SMILES('C[CH]C([O])=C(C)O'),
    E0 = (-143.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,363.145,364.572],'cm^-1')),
        HinderedRotor(inertia=(0.00128459,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110513,'amu*angstrom^2'), symmetry=1, barrier=(10.4043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111849,'amu*angstrom^2'), symmetry=1, barrier=(10.4087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111367,'amu*angstrom^2'), symmetry=1, barrier=(10.4225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.539627,0.0644069,-3.47983e-05,-1.16762e-08,1.24888e-11,-17178.7,28.4164], Tmin=(100,'K'), Tmax=(930.903,'K')), NASAPolynomial(coeffs=[18.8584,0.0154804,-3.95891e-06,6.11722e-10,-4.24718e-14,-21880,-65.5818], Tmin=(930.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = 'C[C]=C(O)[C](C)O(31164)',
    structure = SMILES('C[C]C(O)=C(C)O'),
    E0 = (-33.4494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.743435,'amu*angstrom^2'), symmetry=1, barrier=(17.093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744187,'amu*angstrom^2'), symmetry=1, barrier=(17.1103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.743187,'amu*angstrom^2'), symmetry=1, barrier=(17.0873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744501,'amu*angstrom^2'), symmetry=1, barrier=(17.1175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744324,'amu*angstrom^2'), symmetry=1, barrier=(17.1135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.956428,0.0868802,-9.04217e-05,4.47944e-08,-8.25531e-12,-3825.61,29.7053], Tmin=(100,'K'), Tmax=(1531.6,'K')), NASAPolynomial(coeffs=[24.0031,0.00827556,-2.96485e-07,-1.55319e-10,1.54736e-14,-9897.26,-96.2123], Tmin=(1531.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.4494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])=C(O)CC(4613)',
    structure = SMILES('[CH2]C([O])=C(O)CC'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(O)[C](O)CC(4632)',
    structure = SMILES('[CH]C(O)=C(O)CC'),
    E0 = (-110.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1200,1600],'cm^-1')),
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
    label = '[CH2]C=C(O)C([CH2])O(31165)',
    structure = SMILES('[CH2]C=C(O)C([CH2])O'),
    E0 = (-76.8244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,311.91],'cm^-1')),
        HinderedRotor(inertia=(0.0017328,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25271,'amu*angstrom^2'), symmetry=1, barrier=(17.4465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252707,'amu*angstrom^2'), symmetry=1, barrier=(17.4466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25271,'amu*angstrom^2'), symmetry=1, barrier=(17.4465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252711,'amu*angstrom^2'), symmetry=1, barrier=(17.4465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.14439,0.0719594,-4.67619e-05,-4.40244e-09,1.07502e-11,-9089.27,29.9614], Tmin=(100,'K'), Tmax=(941.8,'K')), NASAPolynomial(coeffs=[21.49,0.0129779,-3.27478e-06,5.2804e-10,-3.87317e-14,-14514.8,-79.1975], Tmin=(941.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.8244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C=C(O)[C](C)O(31166)',
    structure = SMILES('[CH2]C=C(O)[C](C)O'),
    E0 = (-111.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,400.151],'cm^-1')),
        HinderedRotor(inertia=(0.124318,'amu*angstrom^2'), symmetry=1, barrier=(14.0911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124327,'amu*angstrom^2'), symmetry=1, barrier=(14.0898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124332,'amu*angstrom^2'), symmetry=1, barrier=(14.0949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124279,'amu*angstrom^2'), symmetry=1, barrier=(14.0919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124231,'amu*angstrom^2'), symmetry=1, barrier=(14.0903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.303811,0.0806589,-8.15741e-05,4.06288e-08,-7.74696e-12,-13277.8,30.367], Tmin=(100,'K'), Tmax=(1379.32,'K')), NASAPolynomial(coeffs=[21.1156,0.0138305,-3.77368e-06,5.48421e-10,-3.34848e-14,-18738.4,-78.2323], Tmin=(1379.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C2CsJOH)"""),
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
    label = 'C=C(O)C(=C)O(5772)',
    structure = SMILES('C=C(O)C(=C)O'),
    E0 = (-333.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.01466,'amu*angstrom^2'), symmetry=1, barrier=(23.3291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01413,'amu*angstrom^2'), symmetry=1, barrier=(23.3169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01166,'amu*angstrom^2'), symmetry=1, barrier=(23.2599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.498033,0.0751665,-8.30988e-05,4.22746e-08,-7.8095e-12,-39972.2,21.9417], Tmin=(100,'K'), Tmax=(1586.91,'K')), NASAPolynomial(coeffs=[21.8153,0.0015529,2.90201e-06,-7.522e-10,5.55296e-14,-44866.9,-89.1201], Tmin=(1586.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-333.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(O)=C(C)O(31167)',
    structure = SMILES('C=CC(O)=C(C)O'),
    E0 = (-370.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.1976,0.0877436,-9.10797e-05,4.46402e-08,-8.05724e-12,-44357.6,27.774], Tmin=(100,'K'), Tmax=(1598.72,'K')), NASAPolynomial(coeffs=[24.3035,0.00683599,8.79619e-07,-3.98758e-10,3.21892e-14,-50325.6,-100.383], Tmin=(1598.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(O)C(O)[C]C(31169)',
    structure = SMILES('C=C(O)C(O)[C]C'),
    E0 = (0.198124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    label = '[CH]C(O)C(O)=CC(31170)',
    structure = SMILES('[CH]C(O)C(O)=CC'),
    E0 = (10.3887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
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
    label = 'CC1CC(O)=C1O(31171)',
    structure = SMILES('CC1CC(O)=C1O'),
    E0 = (-316.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716318,0.0465129,4.18254e-05,-1.08121e-07,5.08748e-11,-37912.2,22.5083], Tmin=(100,'K'), Tmax=(918.644,'K')), NASAPolynomial(coeffs=[26.071,0.00311907,3.26936e-06,-7.39926e-10,4.4025e-14,-45398,-113.051], Tmin=(918.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene)"""),
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
    label = 'CC=C(O)C(C)=O(31168)',
    structure = SMILES('CC=C(O)C(C)=O'),
    E0 = (-369.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706729,0.0660995,-5.49688e-05,2.3391e-08,-4.00086e-12,-44265.9,23.7086], Tmin=(100,'K'), Tmax=(1391.19,'K')), NASAPolynomial(coeffs=[15.2261,0.0243525,-9.95644e-06,1.82065e-09,-1.24604e-13,-48305.7,-51.1328], Tmin=(1391.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs)"""),
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
    E0 = (75.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-20.2932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (75.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-20.2932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (49.4191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-40.0423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (176.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (79.7437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (88.9981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (24.8986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (25.4569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-94.6401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-27.1186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (72.9184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-26.5603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (82.7311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-109.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-25.0814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-150.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-76.4621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-51.8511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-77.2676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (85.9973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-206.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (57.0465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (67.2371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-247.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-230.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-165.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction27',
    reactants = ['OH(5)', 'C=[C]C(O)=CC(31155)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Cd_rad/Cd] for rate rule [O_pri_rad;Cd_rad/Cd]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C=C([O])C(O)=CC(31156)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['OH(5)', 'C=C(O)[C]=CC(31157)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad/Cd;Y_rad] for rate rule [Cd_rad/Cd;O_pri_rad]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', 'C=C(O)C([O])=CC(4627)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [O_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH3(17)', '[CH]=C(O)C(=C)O(8478)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.23e+13,'cm^3/(mol*s)','+|-',1.81e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 69 used for C_methyl;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;C_methyl]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]C=C(O)C(=C)O(31158)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.156e+12,'cm^3/(mol*s)'), n=0.461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 15 used for C_rad/H2/Cd;H_rad
Exact match found for rate rule [C_rad/H2/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2COH(99)', 'CC=[C]O(31159)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(556926,'m^3/(mol*s)'), n=0.4, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;Cd_rad] for rate rule [Cd_rad/NonDe;Cd_rad/NonDe]
Euclidian distance = 2.82842712475
family: R_Recombination
Ea raised from -2.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', 'C=C(O)C(O)=[C]C(31160)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH]=C(O)C(O)=CC(31161)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([O])C(O)=CC(31162)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C(O)C([O])[CH]C(4629)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]CC(O)=C([CH2])O(4611)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(O)C(=O)[CH]C(4622)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(O)C(O)=[C]C(31163)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])C(O)[CH]C(11048)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C(O)C(O)[CH]C(16341)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[CH]C([O])=C(C)O(4612)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C[C]=C(O)[C](C)O(31164)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([O])=C(O)CC(4613)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C(O)[C](O)CC(4632)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C=C(O)C([CH2])O(31165)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C=C(O)[C](C)O(31166)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CH2(S)(23)', 'C=C(O)C(=C)O(5772)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.94e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['C=CC(O)=C(C)O(31167)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.02873e+09,'s^-1'), n=1.23767, Ea=(163.714,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH3_1;CdH2_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C(O)C(O)[C]C(31169)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2C;CH]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]C(O)C(O)=CC(31170)'],
    products = ['C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CsJ2-C;singletcarbene;CH] for rate rule [CsJ2-C;CsJ2H;CH]
Euclidian distance = 1.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['CC1CC(O)=C1O(31171)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction97',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['C=C(O)C(=O)CC(4626)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C(O)C(O)=CC(5562)'],
    products = ['CC=C(O)C(C)=O(31168)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

network(
    label = '5092',
    isomers = [
        'C=C(O)C(O)=CC(5562)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5092',
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

