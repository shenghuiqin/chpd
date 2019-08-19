species(
    label = 'C=[C]C(C)O[CH]O(21255)',
    structure = SMILES('C=[C]C(C)O[CH]O'),
    E0 = (66.2627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.549696,0.0757718,-7.90909e-05,4.36831e-08,-9.72106e-12,8093.97,29.245], Tmin=(100,'K'), Tmax=(1084.57,'K')), NASAPolynomial(coeffs=[13.7911,0.0269358,-1.15484e-05,2.16548e-09,-1.50918e-13,5221.74,-35.712], Tmin=(1084.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.2627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(OCJO) + radical(Cds_S)"""),
)

species(
    label = 'HOCHO(59)',
    structure = SMILES('O=CO'),
    E0 = (-389.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1051.41,1051.94,1052.55,2787.25,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00861396,'amu*angstrom^2'), symmetry=1, barrier=(6.76824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4005.91,'J/mol'), sigma=(3.626,'angstroms'), dipoleMoment=(1.7,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46778.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-45330.3,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-389.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HOCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH3CHCCH2(18175)',
    structure = SMILES('C=C=CC'),
    E0 = (145.615,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.759584,'amu*angstrom^2'), symmetry=1, barrier=(17.4643,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2996.71,'J/mol'), sigma=(5.18551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=468.08 K, Pc=48.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74635,0.0218189,8.22353e-06,-2.14768e-08,8.55624e-12,17563.6,12.7381], Tmin=(100,'K'), Tmax=(1025.6,'K')), NASAPolynomial(coeffs=[6.82078,0.0192338,-7.45622e-06,1.36536e-09,-9.53195e-14,16028,-10.4333], Tmin=(1025.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""CH3CHCCH2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C]C(C)OC=O(25315)',
    structure = SMILES('C=[C]C(C)OC=O'),
    E0 = (-102.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,366.954,367.075,367.232],'cm^-1')),
        HinderedRotor(inertia=(0.00125125,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0012508,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310625,'amu*angstrom^2'), symmetry=1, barrier=(29.7002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31061,'amu*angstrom^2'), symmetry=1, barrier=(29.7022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33263,0.0493572,-1.83577e-05,-7.07195e-09,4.7321e-12,-12228.2,27.1109], Tmin=(100,'K'), Tmax=(1146.68,'K')), NASAPolynomial(coeffs=[12.6009,0.0272556,-1.19538e-05,2.29056e-09,-1.62101e-13,-15943.6,-33.7266], Tmin=(1146.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C(C)O[CH]O(25316)',
    structure = SMILES('C=C=C(C)O[CH]O'),
    E0 = (-52.1685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.70987,'amu*angstrom^2'), symmetry=1, barrier=(16.3213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709909,'amu*angstrom^2'), symmetry=1, barrier=(16.3222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.71012,'amu*angstrom^2'), symmetry=1, barrier=(16.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709676,'amu*angstrom^2'), symmetry=1, barrier=(16.3169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0548679,0.0795368,-8.34418e-05,4.07662e-08,-6.90037e-12,-6125.97,25.5769], Tmin=(100,'K'), Tmax=(992.587,'K')), NASAPolynomial(coeffs=[19.1539,0.0151212,-5.06404e-06,8.63236e-10,-5.87077e-14,-10535.7,-69.5371], Tmin=(992.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.1685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(OCJO)"""),
)

species(
    label = 'C#CC(C)O[CH]O(25317)',
    structure = SMILES('C#CC(C)O[CH]O'),
    E0 = (-5.62903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.232045,0.0779753,-8.73943e-05,5.01659e-08,-1.12247e-11,-537.121,27.5221], Tmin=(100,'K'), Tmax=(1101.04,'K')), NASAPolynomial(coeffs=[16.632,0.0183957,-6.2267e-06,1.02028e-09,-6.58124e-14,-4148.56,-53.177], Tmin=(1101.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.62903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(OCJO)"""),
)

species(
    label = 'C=[C][CH]C(18176)',
    structure = SMILES('[CH2][C]=CC'),
    E0 = (361.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.352622,'amu*angstrom^2'), symmetry=1, barrier=(8.10748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828631,'amu*angstrom^2'), symmetry=1, barrier=(19.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42015,0.030446,-1.69076e-05,4.64684e-09,-5.12013e-13,43485.7,14.8304], Tmin=(100,'K'), Tmax=(2065.83,'K')), NASAPolynomial(coeffs=[10.7464,0.014324,-5.20136e-06,8.69079e-10,-5.48385e-14,40045.6,-31.3799], Tmin=(2065.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(Allyl_P)"""),
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
    label = 'C=C=CO[CH]O(25318)',
    structure = SMILES('C=C=CO[CH]O'),
    E0 = (-10.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05,'amu*angstrom^2'), symmetry=1, barrier=(24.1416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0484,'amu*angstrom^2'), symmetry=1, barrier=(24.1047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04819,'amu*angstrom^2'), symmetry=1, barrier=(24.1,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194548,0.0709862,-7.91838e-05,4.03018e-08,-7.46796e-12,-1078.75,24.4774], Tmin=(100,'K'), Tmax=(1567.3,'K')), NASAPolynomial(coeffs=[21.434,0.000652077,2.61506e-06,-6.53439e-10,4.76096e-14,-5999.57,-83.6565], Tmin=(1567.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(OCJO)"""),
)

species(
    label = 'C=[C]C(C)OC[O](25319)',
    structure = SMILES('C=[C]C(C)OC[O]'),
    E0 = (104.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790347,0.0607237,-3.67548e-05,5.71649e-09,1.50056e-12,12716.3,29.0714], Tmin=(100,'K'), Tmax=(1152.8,'K')), NASAPolynomial(coeffs=[15.0578,0.0261204,-1.11202e-05,2.10561e-09,-1.4823e-13,8436.63,-46.0845], Tmin=(1152.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(OCOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C](C)O[CH]O(21253)',
    structure = SMILES('[CH2]C=C(C)O[CH]O'),
    E0 = (-77.2736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.455408,0.0828168,-8.41635e-05,4.18356e-08,-7.93018e-12,-9120.54,28.6433], Tmin=(100,'K'), Tmax=(1402.27,'K')), NASAPolynomial(coeffs=[21.9423,0.0131103,-3.37684e-06,4.69968e-10,-2.80616e-14,-14830.2,-84.9462], Tmin=(1402.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.2736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(OCJO)"""),
)

species(
    label = '[CH]=CC(C)O[CH]O(21257)',
    structure = SMILES('[CH]=CC(C)O[CH]O'),
    E0 = (75.5171,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355831,0.0765347,-7.7517e-05,4.04005e-08,-8.36095e-12,9216.89,29.8672], Tmin=(100,'K'), Tmax=(1172.53,'K')), NASAPolynomial(coeffs=[15.9294,0.0234062,-9.54983e-06,1.75602e-09,-1.21346e-13,5564.82,-47.745], Tmin=(1172.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.5171,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(OCJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][C](C)OCO(25320)',
    structure = SMILES('[CH2][C]=C(C)OCO'),
    E0 = (-28.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0841719,0.065004,-7.74544e-06,-5.87057e-08,3.36518e-11,-3254.13,27.384], Tmin=(100,'K'), Tmax=(919.289,'K')), NASAPolynomial(coeffs=[26.9426,0.00347141,2.36922e-06,-5.64242e-10,3.39813e-14,-10530.4,-112.649], Tmin=(919.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(C=C)O[CH]O(21218)',
    structure = SMILES('[CH2]C(C=C)O[CH]O'),
    E0 = (38.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3800.23,'J/mol'), sigma=(6.47665,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.59 K, Pc=31.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.197079,0.0843085,-9.8998e-05,6.05916e-08,-1.47169e-11,4818.78,29.7911], Tmin=(100,'K'), Tmax=(1005.09,'K')), NASAPolynomial(coeffs=[15.0784,0.0250865,-1.06174e-05,1.97133e-09,-1.36454e-13,1827.27,-42.0786], Tmin=(1005.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(OCJO)"""),
)

species(
    label = '[CH2]C([C]=C)OCO(21256)',
    structure = SMILES('[CH2]C([C]=C)OCO'),
    E0 = (87.8126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.875591,0.0852489,-8.7039e-05,4.22065e-08,-7.63924e-12,10755.8,34.3339], Tmin=(100,'K'), Tmax=(1546.97,'K')), NASAPolynomial(coeffs=[24.1098,0.00849619,-8.37412e-07,-1.84198e-11,4.93294e-15,4479.06,-92.4089], Tmin=(1546.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.8126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CC(C)O[CH][O](21258)',
    structure = SMILES('C=CC(C)O[CH][O]'),
    E0 = (55.8188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89367,0.0560825,-1.31518e-05,-7.7301e-08,8.2709e-11,6779.24,24.5982], Tmin=(100,'K'), Tmax=(452.439,'K')), NASAPolynomial(coeffs=[4.76733,0.0447169,-2.20198e-05,4.35547e-09,-3.11016e-13,6375.5,11.4255], Tmin=(452.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.8188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(OCJO) + radical(OCOJ)"""),
)

species(
    label = '[CH]=[C]C(C)OCO(25321)',
    structure = SMILES('[CH]=[C]C(C)OCO'),
    E0 = (124.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272817,0.0659899,-2.62381e-05,-2.77906e-08,1.94596e-11,15110.6,30.847], Tmin=(100,'K'), Tmax=(940.407,'K')), NASAPolynomial(coeffs=[22.5914,0.0109992,-2.23317e-06,3.55461e-10,-2.92185e-14,9146.72,-84.8464], Tmin=(940.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[O][CH]O(171)',
    structure = SMILES('[O][CH]O'),
    E0 = (5.37818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,1989.68],'cm^-1')),
        HinderedRotor(inertia=(0.0937094,'amu*angstrom^2'), symmetry=1, barrier=(2.15456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40204,0.0209571,-4.96202e-05,5.95355e-08,-2.44634e-11,660.882,10.3084], Tmin=(100,'K'), Tmax=(868.739,'K')), NASAPolynomial(coeffs=[-0.343229,0.0179325,-9.40003e-06,1.81361e-09,-1.23835e-13,2076.48,32.2523], Tmin=(868.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.37818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(OCJO) + radical(OCOJ)"""),
)

species(
    label = 'C=C=C(C)OCO(25322)',
    structure = SMILES('C=C=C(C)OCO'),
    E0 = (-241.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.155619,0.0625561,-6.12712e-07,-6.60992e-08,3.62413e-11,-28842.2,25.4269], Tmin=(100,'K'), Tmax=(920.483,'K')), NASAPolynomial(coeffs=[27.0416,0.00316831,2.55093e-06,-5.90271e-10,3.50206e-14,-36225.5,-115.274], Tmin=(920.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC(C)OC=O(876)',
    structure = SMILES('C=CC(C)OC=O'),
    E0 = (-340.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39287,0.0445223,6.82342e-06,-3.61916e-08,1.53325e-11,-40832.8,26.4155], Tmin=(100,'K'), Tmax=(1047.81,'K')), NASAPolynomial(coeffs=[13.0834,0.0288876,-1.23008e-05,2.38437e-09,-1.72327e-13,-44874.3,-38.1254], Tmin=(1047.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-OdOsH)"""),
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
    label = 'C=[C]CO[CH]O(25323)',
    structure = SMILES('C=[C]CO[CH]O'),
    E0 = (108.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,3025,407.5,1350,352.5,368.62,369.366,369.966,370.456],'cm^-1')),
        HinderedRotor(inertia=(0.00123094,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.086085,'amu*angstrom^2'), symmetry=1, barrier=(8.36671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.085679,'amu*angstrom^2'), symmetry=1, barrier=(8.36367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859646,'amu*angstrom^2'), symmetry=1, barrier=(8.36095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44516,0.0620532,-7.41096e-05,4.08391e-08,-3.894e-12,13145.6,23.6439], Tmin=(100,'K'), Tmax=(612.808,'K')), NASAPolynomial(coeffs=[8.76727,0.0253454,-1.13938e-05,2.13202e-09,-1.4637e-13,12040,-9.7938], Tmin=(612.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(OCJO)"""),
)

species(
    label = 'C=C(C)[CH]O[CH]O(25324)',
    structure = SMILES('[CH2]C(C)=CO[CH]O'),
    E0 = (-74.5371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10035,0.0890558,-9.25549e-05,4.56747e-08,-8.38048e-12,-8761.29,29.7161], Tmin=(100,'K'), Tmax=(1539.95,'K')), NASAPolynomial(coeffs=[24.6365,0.00825443,-2.61303e-07,-1.63062e-10,1.5953e-14,-15033.9,-100.191], Tmin=(1539.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.5371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(OCJO) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1C(C)OC1O(25325)',
    structure = SMILES('C=C1C(C)OC1O'),
    E0 = (-239.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40964,0.0398403,2.85601e-05,-6.5621e-08,2.74821e-11,-28718.2,23.42], Tmin=(100,'K'), Tmax=(991.121,'K')), NASAPolynomial(coeffs=[15.7326,0.0233725,-9.07842e-06,1.77748e-09,-1.33084e-13,-33587.7,-55.795], Tmin=(991.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]C(C)[O](3245)',
    structure = SMILES('C=[C]C(C)[O]'),
    E0 = (278.874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,461.597,461.616],'cm^-1')),
        HinderedRotor(inertia=(0.0767265,'amu*angstrom^2'), symmetry=1, barrier=(11.6024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0767264,'amu*angstrom^2'), symmetry=1, barrier=(11.6031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05953,0.035711,-1.16803e-05,-8.94971e-09,5.51275e-12,33616.7,20.5932], Tmin=(100,'K'), Tmax=(1043.68,'K')), NASAPolynomial(coeffs=[10.303,0.0183923,-7.30569e-06,1.36087e-09,-9.60942e-14,31118.5,-23.2536], Tmin=(1043.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_S)"""),
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
    label = 'C[CH]O[CH]O(2207)',
    structure = SMILES('C[CH]O[CH]O'),
    E0 = (-58.3656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,335.907,336.099,336.624],'cm^-1')),
        HinderedRotor(inertia=(0.00149336,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145994,'amu*angstrom^2'), symmetry=1, barrier=(11.7512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146216,'amu*angstrom^2'), symmetry=1, barrier=(11.7502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146642,'amu*angstrom^2'), symmetry=1, barrier=(11.7484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16954,0.062401,-8.07621e-05,5.38881e-08,-1.4043e-11,-6917.97,20.4625], Tmin=(100,'K'), Tmax=(946.379,'K')), NASAPolynomial(coeffs=[12.3848,0.0149995,-5.63345e-06,9.66022e-10,-6.32765e-14,-9040.81,-33.0269], Tmin=(946.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.3656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsOsHH) + radical(OCJO) + radical(CCsJOCs)"""),
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
    E0 = (66.2627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (145.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (172.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (219.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (66.2627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (159.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (228.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (248.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (180.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (208.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (218.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (150.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (217.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (188.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (366.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (129.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (105.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (528.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (210.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (74.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (484.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (342.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['HOCHO(59)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=[C]C(C)OC=O(25315)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Od_CO-NdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C=C(C)O[CH]O(25316)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#CC(C)O[CH]O(25317)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['HOCHO(59)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4000,'m^3/(mol*s)'), n=1.39, Ea=(94.4171,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_CO-NdH;YJ] for rate rule [Od_CO-NdH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 92.5 to 94.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', 'C=C=CO[CH]O(25318)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0105855,'m^3/(mol*s)'), n=2.41, Ea=(33.7021,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]C(C)OC[O](25319)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(779620,'s^-1'), n=2.1225, Ea=(123.773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=C[C](C)O[CH]O(21253)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC(C)O[CH]O(21257)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=[C][C](C)OCO(25320)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeO;XH_out] for rate rule [R3H_SS_O;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['[CH2]C(C=C)O[CH]O(21218)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([C]=C)OCO(21256)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(304,'s^-1'), n=2.77, Ea=(62.3834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 324 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=CC(C)O[CH][O](21258)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]C(C)OCO(25321)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.1728e+06,'s^-1'), n=1.70245, Ea=(63.8935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cs_H_out_1H] + [R5Hall;Cd_rad_out_singleH;Cs_H_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]O(171)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=C=C(C)OCO(25322)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3radExo;Y_rad_NDe;XH_Rrad] for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=CC(C)OC=O(876)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2(S)(23)', 'C=[C]CO[CH]O(25323)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=C(C)[CH]O[CH]O(25324)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HR!H)CJ;CJ;CH3] + [cCs(-HR!H)CJ;CdsJ;C] for rate rule [cCs(-HR!H)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]C(C)O[CH]O(21255)'],
    products = ['C=C1C(C)OC1O(25325)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HCOH(T)(285)', 'C=[C]C(C)[O](3245)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H2CC(41)', 'C[CH]O[CH]O(2207)'],
    products = ['C=[C]C(C)O[CH]O(21255)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4227',
    isomers = [
        'C=[C]C(C)O[CH]O(21255)',
    ],
    reactants = [
        ('HOCHO(59)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4227',
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

