species(
    label = '[CH2]OCC(O)=C(O)[CH]C(5903)',
    structure = SMILES('[CH2]OCC(O)=C(O)[CH]C'),
    E0 = (-226.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73665,0.104275,-0.000104635,5.14809e-08,-9.5976e-12,-27002.1,41.8774], Tmin=(100,'K'), Tmax=(1465.31,'K')), NASAPolynomial(coeffs=[26.5599,0.0162077,-3.40186e-06,3.82299e-10,-1.94725e-14,-34132.7,-101.483], Tmin=(1465.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC) + radical(CCJCO)"""),
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
    label = 'C[CH][C](O)C1(O)COC1(32242)',
    structure = SMILES('C[CH][C](O)C1(O)COC1'),
    E0 = (-154.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.388932,0.0918733,-7.80707e-05,2.3224e-08,3.50841e-12,-18365.8,33.6819], Tmin=(100,'K'), Tmax=(831.641,'K')), NASAPolynomial(coeffs=[17.4491,0.0295811,-8.11106e-06,1.12715e-09,-6.5215e-14,-22145.6,-53.9749], Tmin=(831.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(C2CsJOH) + radical(CCJCO)"""),
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
    label = '[CH2]OCC(O)=C(O)C=C(32243)',
    structure = SMILES('[CH2]OCC(O)=C(O)C=C'),
    E0 = (-315.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82099,0.108873,-0.000118536,6.23886e-08,-1.23688e-11,-37680.7,36.7765], Tmin=(100,'K'), Tmax=(1376.14,'K')), NASAPolynomial(coeffs=[27.4258,0.0128806,-1.9345e-06,1.02479e-10,5.11711e-17,-44690.5,-109.882], Tmin=(1376.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CsJOCC)"""),
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
    label = '[CH2]OCC(O)=C=CC(32244)',
    structure = SMILES('[CH2]OCC(O)=C=CC'),
    E0 = (-48.5002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0228084,0.0835091,-8.14288e-05,4.20667e-08,-8.71579e-12,-5686.66,30.9462], Tmin=(100,'K'), Tmax=(1166.62,'K')), NASAPolynomial(coeffs=[15.9107,0.0290347,-1.13882e-05,2.04231e-09,-1.38871e-13,-9393.71,-48.1522], Tmin=(1166.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.5002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]CC(O)=C(O)CO[CH2](5904)',
    structure = SMILES('[CH2]CC(O)=C(O)CO[CH2]'),
    E0 = (-221.037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,325,375,415,465,420,450,1700,1750,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57065,0.106525,-0.000111257,5.76506e-08,-1.14517e-11,-26370.4,40.5637], Tmin=(100,'K'), Tmax=(1316.15,'K')), NASAPolynomial(coeffs=[25.5478,0.0186985,-4.99788e-06,7.05067e-10,-4.19267e-14,-33040.3,-95.9371], Tmin=(1316.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C(O)=C(O)[CH]OC(32245)',
    structure = SMILES('C[CH]C(O)=C(O)[CH]OC'),
    E0 = (-303.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03201,0.0987035,-9.56387e-05,4.70416e-08,-9.06385e-12,-36276.8,38.5883], Tmin=(100,'K'), Tmax=(1269.48,'K')), NASAPolynomial(coeffs=[22.2451,0.0253601,-8.97775e-06,1.53191e-09,-1.01645e-13,-42186.9,-79.2647], Tmin=(1269.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-303.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=CCJ(O)C) + radical(CCJCO)"""),
)

species(
    label = '[CH2]OC[C](O)C(=O)CC(5869)',
    structure = SMILES('[CH2]OCC(O)=C([O])CC'),
    E0 = (-288.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4583.77,'J/mol'), sigma=(7.56195,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=715.97 K, Pc=24.05 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02323,0.0992727,-9.95574e-05,5.12942e-08,-1.0335e-11,-34505.5,38.3117], Tmin=(100,'K'), Tmax=(1218.57,'K')), NASAPolynomial(coeffs=[21.5741,0.0250953,-8.24785e-06,1.3392e-09,-8.61932e-14,-40012.8,-75.1742], Tmin=(1218.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]OCC([O])=C(O)CC(5905)',
    structure = SMILES('[CH2]OCC([O])=C(O)CC'),
    E0 = (-288.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02323,0.0992727,-9.95574e-05,5.12942e-08,-1.0335e-11,-34505.5,38.3117], Tmin=(100,'K'), Tmax=(1218.57,'K')), NASAPolynomial(coeffs=[21.5741,0.0250953,-8.24785e-06,1.3392e-09,-8.61932e-14,-40012.8,-75.1742], Tmin=(1218.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]O[CH]C(O)=C(O)CC(5906)',
    structure = SMILES('[CH2]O[CH]C(O)=C(O)CC'),
    E0 = (-315.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32566,0.10256,-0.000102539,5.15904e-08,-1.00664e-11,-37723,38.8356], Tmin=(100,'K'), Tmax=(1289.53,'K')), NASAPolynomial(coeffs=[24.018,0.0221635,-6.94729e-06,1.09914e-09,-6.98952e-14,-44111,-89.3024], Tmin=(1289.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C[CH]C(O)=C([O])COC(32246)',
    structure = SMILES('C[CH]C(O)=C([O])COC'),
    E0 = (-276.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729685,0.0954129,-9.26286e-05,4.66971e-08,-9.31075e-12,-33059.4,38.0652], Tmin=(100,'K'), Tmax=(1221.52,'K')), NASAPolynomial(coeffs=[19.7341,0.0284016,-1.03398e-05,1.78621e-09,-1.19107e-13,-38058.7,-64.7555], Tmin=(1221.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH]C([O])=C(O)COC(5908)',
    structure = SMILES('C[CH]C([O])=C(O)COC'),
    E0 = (-276.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.729685,0.0954129,-9.26286e-05,4.66971e-08,-9.31075e-12,-33059.4,38.0652], Tmin=(100,'K'), Tmax=(1221.52,'K')), NASAPolynomial(coeffs=[19.7341,0.0284016,-1.03398e-05,1.78621e-09,-1.19107e-13,-38058.7,-64.7555], Tmin=(1221.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-276.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2][CH]C(O)=C(O)COC(32247)',
    structure = SMILES('[CH2]C=C(O)[C](O)COC'),
    E0 = (-250.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.872169,0.106563,-0.000118385,6.85936e-08,-1.58185e-11,-29996.3,36.4475], Tmin=(100,'K'), Tmax=(1055.28,'K')), NASAPolynomial(coeffs=[18.5239,0.0330426,-1.3882e-05,2.57468e-09,-1.7838e-13,-34090,-58.171], Tmin=(1055.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-250.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(C2CsJOH)"""),
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
    label = '[CH2]O[CH2](103)',
    structure = SMILES('[CH2]O[CH2]'),
    E0 = (178.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,672.809],'cm^-1')),
        HinderedRotor(inertia=(0.0102546,'amu*angstrom^2'), symmetry=1, barrier=(3.29362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14325,'amu*angstrom^2'), symmetry=1, barrier=(3.29359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05226,0.0166761,-4.43211e-07,-1.08823e-08,5.26552e-12,21459,11.4763], Tmin=(100,'K'), Tmax=(996.635,'K')), NASAPolynomial(coeffs=[7.42989,0.00875132,-3.03205e-06,5.59498e-10,-4.03289e-14,20107.5,-12.0316], Tmin=(996.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CsJOCH3) + radical(CsJOCH3)"""),
)

species(
    label = 'C[CH]C(O)=[C]O(31655)',
    structure = SMILES('C[CH]C(O)=[C]O'),
    E0 = (-0.226743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,1685,370,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,357.597],'cm^-1')),
        HinderedRotor(inertia=(0.161034,'amu*angstrom^2'), symmetry=1, barrier=(14.6141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161061,'amu*angstrom^2'), symmetry=1, barrier=(14.6129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161086,'amu*angstrom^2'), symmetry=1, barrier=(14.6132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161124,'amu*angstrom^2'), symmetry=1, barrier=(14.6146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01597,0.0512908,-1.42347e-05,-3.65113e-08,2.34132e-11,93.6637,26.3323], Tmin=(100,'K'), Tmax=(901.306,'K')), NASAPolynomial(coeffs=[21.0416,0.00174384,2.77393e-06,-6.8075e-10,4.66332e-14,-5113.55,-77.0604], Tmin=(901.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.226743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCJCO)"""),
)

species(
    label = '[CH2]OCC1(O)[C](O)C1C(32248)',
    structure = SMILES('[CH2]OCC1(O)[C](O)C1C'),
    E0 = (-152.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.835198,0.0998388,-0.000102429,5.47211e-08,-1.15514e-11,-18118.1,35.4215], Tmin=(100,'K'), Tmax=(1156.42,'K')), NASAPolynomial(coeffs=[19.5996,0.0291552,-1.07447e-05,1.86529e-09,-1.24709e-13,-22844.3,-66.1343], Tmin=(1156.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + ring(Cyclopropane) + radical(CsJOCC) + radical(C2CsJOH)"""),
)

species(
    label = 'C[CH]C1(O)COC[C]1O(32249)',
    structure = SMILES('C[CH]C1(O)COC[C]1O'),
    E0 = (-232.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.978735,0.0913393,-8.33678e-05,4.08527e-08,-7.73462e-12,-27770.6,34.2792], Tmin=(100,'K'), Tmax=(1467.4,'K')), NASAPolynomial(coeffs=[18.7932,0.0265823,-6.0702e-06,6.91141e-10,-3.3011e-14,-32404,-64.707], Tmin=(1467.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-232.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(O)=C(O)COC(32250)',
    structure = SMILES('C=CC(O)=C(O)COC'),
    E0 = (-502.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60513,0.108836,-0.000115283,6.06919e-08,-1.22925e-11,-60276.3,34.313], Tmin=(100,'K'), Tmax=(1261.39,'K')), NASAPolynomial(coeffs=[25.5338,0.0197255,-5.68867e-06,8.52069e-10,-5.26242e-14,-66880.2,-101.957], Tmin=(1261.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-502.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC1COCC(O)=C1O(31892)',
    structure = SMILES('CC1COCC(O)=C1O'),
    E0 = (-581.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.24261,0.0685026,1.92489e-05,-9.46054e-08,4.80715e-11,-69716.1,26.111], Tmin=(100,'K'), Tmax=(907.85,'K')), NASAPolynomial(coeffs=[27.1019,0.0140979,-3.41928e-08,-2.75012e-10,1.82965e-14,-77403.9,-118.163], Tmin=(907.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-581.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(36dihydro2hpyran)"""),
)

species(
    label = '[CH2]OCC(=O)C(O)[CH]C(32251)',
    structure = SMILES('[CH2]OCC(=O)C(O)[CH]C'),
    E0 = (-188.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.478176,0.0874897,-7.37404e-05,3.16991e-08,-5.43283e-12,-22523.8,41.6148], Tmin=(100,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[19.8599,0.0293812,-1.14815e-05,2.05206e-09,-1.38723e-13,-28218.5,-63.3477], Tmin=(1400,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-OdCsCs) + radical(CsJOCC) + radical(CCJCO)"""),
)

species(
    label = '[CH2]OCC(O)C(=O)[CH]C(5921)',
    structure = SMILES('[CH2]OCC(O)C([O])=CC'),
    E0 = (-253.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399044,0.0944428,-9.40138e-05,5.00284e-08,-1.0736e-11,-30318.7,38.2312], Tmin=(100,'K'), Tmax=(1124.05,'K')), NASAPolynomial(coeffs=[16.4619,0.0344421,-1.39454e-05,2.54049e-09,-1.74199e-13,-34109.3,-45.085], Tmin=(1124.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCC) + radical(C=C(C)OJ)"""),
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
    label = 'C[CH]C(O)=C(O)C[O](31659)',
    structure = SMILES('C[CH]C(O)=C(O)C[O]'),
    E0 = (-209.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.397705,0.0846344,-8.53876e-05,4.31776e-08,-8.44867e-12,-24991.5,34.8991], Tmin=(100,'K'), Tmax=(1290.94,'K')), NASAPolynomial(coeffs=[20.7241,0.0174893,-5.39474e-06,8.48309e-10,-5.38686e-14,-30303.3,-71.8468], Tmin=(1290.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCOJ) + radical(CCJCO)"""),
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
    label = '[CH2]OCC(O)=[C]O(32252)',
    structure = SMILES('[CH2]OCC(O)=[C]O'),
    E0 = (-122.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,350,440,435,1725,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (102.089,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.146746,0.0773223,-8.77774e-05,4.7616e-08,-9.62493e-12,-14524.2,31.8725], Tmin=(100,'K'), Tmax=(1379.9,'K')), NASAPolynomial(coeffs=[20.6454,0.00677118,8.84724e-08,-2.32952e-10,2.20628e-14,-19283.7,-71.5865], Tmin=(1379.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]OCC(=O)C(O)=CC(32253)',
    structure = SMILES('[CH2]OCC(=O)C(O)=CC'),
    E0 = (-305.813,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.328533,0.0938792,-9.44696e-05,4.956e-08,-1.04588e-11,-36623.8,32.4161], Tmin=(100,'K'), Tmax=(1140.76,'K')), NASAPolynomial(coeffs=[17.0636,0.032895,-1.4281e-05,2.69745e-09,-1.88824e-13,-40591.8,-53.7815], Tmin=(1140.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-305.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]OC=C(O)C(O)=CC(32254)',
    structure = SMILES('[CH2]OC=C(O)C(O)=CC'),
    E0 = (-367.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.90214,0.119923,-0.000132821,6.81365e-08,-1.28232e-11,-43983.2,37.0485], Tmin=(100,'K'), Tmax=(1525.19,'K')), NASAPolynomial(coeffs=[32.6055,0.00436037,2.9019e-06,-8.34736e-10,6.32961e-14,-52204.5,-140.688], Tmin=(1525.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COCJ)"""),
)

species(
    label = '[CH2]OCC([O])C(O)=CC(32255)',
    structure = SMILES('[CH2]OCC([O])C(O)=CC'),
    E0 = (-160.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979578,0.098832,-9.57853e-05,4.72083e-08,-9.14118e-12,-19158.9,38.6231], Tmin=(100,'K'), Tmax=(1260.31,'K')), NASAPolynomial(coeffs=[21.7958,0.0265476,-9.75445e-06,1.70098e-09,-1.14249e-13,-24899.8,-76.5246], Tmin=(1260.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-160.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CsJOCC)"""),
)

species(
    label = '[CH2]O[CH]C(O)C(O)=CC(32256)',
    structure = SMILES('[CH2]O[CH]C(O)C(O)=CC'),
    E0 = (-210.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43046,0.109765,-0.00012098,6.67288e-08,-1.42462e-11,-25145.6,39.8743], Tmin=(100,'K'), Tmax=(1156.21,'K')), NASAPolynomial(coeffs=[23.8868,0.0221767,-7.34713e-06,1.20728e-09,-7.86869e-14,-30999.9,-85.9411], Tmin=(1156.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCC) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]OCC(O)C(O)=[C]C(32257)',
    structure = SMILES('[CH2]OCC(O)C(O)=[C]C'),
    E0 = (-153.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3580,3650,1210,1345,900,1100,1685,370,1380,1390,370,380,2900,435,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.822055,0.103293,-0.000112242,6.34542e-08,-1.42151e-11,-18271.4,38.7225], Tmin=(100,'K'), Tmax=(1088.61,'K')), NASAPolynomial(coeffs=[18.9274,0.0307241,-1.22475e-05,2.21668e-09,-1.51661e-13,-22571.2,-58.2335], Tmin=(1088.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCC) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C(O)C(O)CO[CH2](32258)',
    structure = SMILES('[CH2]C=C(O)C(O)CO[CH2]'),
    E0 = (-239.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00682,0.0988014,-8.8165e-05,3.305e-08,-2.54242e-12,-28641.8,39.1891], Tmin=(100,'K'), Tmax=(993.915,'K')), NASAPolynomial(coeffs=[22.8339,0.0244467,-8.536e-06,1.49599e-09,-1.03335e-13,-34447.4,-81.0479], Tmin=(993.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-239.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CsJOCC) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C(O)[C](O)COC(32259)',
    structure = SMILES('C[C]C(O)=C(O)COC'),
    E0 = (-165.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.39923,0.108329,-0.000115646,6.19322e-08,-1.28752e-11,-19742.6,36.3755], Tmin=(100,'K'), Tmax=(1184.76,'K')), NASAPolynomial(coeffs=[23.796,0.0232651,-7.94908e-06,1.3315e-09,-8.77352e-14,-25712.7,-89.4486], Tmin=(1184.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC1COC[C](O)[C]1O(32260)',
    structure = SMILES('CC1COC[C](O)[C]1O'),
    E0 = (-268.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0785146,0.0846275,-6.3319e-05,1.66608e-08,2.82473e-12,-32161.7,28.9871], Tmin=(100,'K'), Tmax=(867.276,'K')), NASAPolynomial(coeffs=[14.6085,0.03574,-1.13695e-05,1.79007e-09,-1.13129e-13,-35418.2,-43.8653], Tmin=(867.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxane) + radical(C2CsJOH) + radical(C2CsJOH)"""),
)

species(
    label = 'CC=C(O)C(O)=COC(32261)',
    structure = SMILES('CC=C(O)C(O)=COC'),
    E0 = (-559.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.06556,0.122533,-0.000132748,6.7256e-08,-1.25429e-11,-67069.7,35.1525], Tmin=(100,'K'), Tmax=(1536,'K')), NASAPolynomial(coeffs=[32.79,0.0068548,1.99986e-06,-6.81913e-10,5.33654e-14,-75453.4,-144.654], Tmin=(1536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-559.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'CC=C(O)C(=O)COC(32262)',
    structure = SMILES('CC=C(O)C(=O)COC'),
    E0 = (-493.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1392,0.0940875,-9.18296e-05,4.8441e-08,-1.05769e-11,-59217.9,30.0533], Tmin=(100,'K'), Tmax=(1085.16,'K')), NASAPolynomial(coeffs=[14.1694,0.0413443,-1.89233e-05,3.6509e-09,-2.58058e-13,-62323.3,-40.147], Tmin=(1085.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH)"""),
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
    label = '[CH2]OC[C](O)C(=C)O(32263)',
    structure = SMILES('[CH2]OCC(O)=C([CH2])O'),
    E0 = (-244.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3580,3650,1210,1345,900,1100,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,325,375,415,465,420,450,1700,1750,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10028,0.0957872,-0.000106641,5.73741e-08,-1.15814e-11,-29235.4,34.8607], Tmin=(100,'K'), Tmax=(1365.41,'K')), NASAPolynomial(coeffs=[24.2283,0.010625,-1.04184e-06,-6.45313e-11,1.18343e-14,-35130.4,-91.4818], Tmin=(1365.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-244.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CsJOCC) + radical(C=C(O)CJ)"""),
)

species(
    label = 'CC=C(O)C1(O)COC1(31895)',
    structure = SMILES('CC=C(O)C1(O)COC1'),
    E0 = (-464.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.02603,0.0982678,-9.13231e-05,4.32611e-08,-7.61469e-12,-55623.3,36.8082], Tmin=(100,'K'), Tmax=(1678.46,'K')), NASAPolynomial(coeffs=[21.3388,0.0203934,-1.89548e-06,-1.36148e-10,2.2495e-14,-60340.6,-78.7015], Tmin=(1678.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-464.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(440.667,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Oxetane)"""),
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
    E0 = (-226.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-154.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-103.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-196.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-20.1282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-107.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-106.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-68.9999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-174.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (192.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-193.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-145.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-163.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (70.0577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (178.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (4.83535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-140.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-201.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-218.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-98.5336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-86.7073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (172.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (221.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-68.4981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-149.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1.71016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-64.1959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-13.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-76.5937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-3.91829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-200.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-162.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-208.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (175.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-218.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['CH2O(15)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['C[CH][C](O)C1(O)COC1(32242)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.20551e+07,'s^-1'), n=1.225, Ea=(72.3261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 72.1 to 72.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', '[CH2]OCC(O)=C(O)C=C(32243)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2544 used for Cds-HH_Cds-CdH;HJ
Exact match found for rate rule [Cds-HH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2O(15)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.000198678,'m^3/(mol*s)'), n=2.77646, Ea=(45.4073,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_R;CsJ-CdHH] for rate rule [Od_CO-HH;CsJ-CdHH]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OH(5)', '[CH2]OCC(O)=C=CC(32244)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(41610,'cm^3/(mol*s)'), n=2.487, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds;OJ_pri] for rate rule [Ca_Cds-CsH;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -7.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]CC(O)=C(O)CO[CH2](5904)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.99, Ea=(113.805,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 84 used for R2H_S;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['C[CH]C(O)=C(O)[CH]OC(32245)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00568695,'s^-1'), n=4.30267, Ea=(120.197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R3H_SS_O;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2]OC[C](O)C(=O)CC(5869)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(187863,'s^-1'), n=2.03383, Ea=(157.381,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;C_rad_out_H/NonDeC;O_H_out] + [R3H_SS_2Cd;C_rad_out_H/NonDeC;XH_out] for rate rule [R3H_SS_2Cd;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OCC([O])=C(O)CC(5905)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2]O[CH]C(O)=C(O)CC(5906)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+10,'s^-1'), n=0, Ea=(418.4,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['C[CH]C(O)=C([O])COC(32246)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.06724e+06,'s^-1'), n=1.18977, Ea=(32.9114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_2H;XH_out] for rate rule [R5H_SSSS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['C[CH]C([O])=C(O)COC(5908)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(12584.6,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;C_rad_out_2H;XH_out] for rate rule [R6H_RSSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2][CH]C(O)=C(O)COC(32247)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R7HJ_5;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][O](167)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.56662e+07,'m^3/(mol*s)'), n=0.074875, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [C_rad/H2/Cd;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]O[CH2](103)', 'C[CH]C(O)=[C]O(31655)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.28206e+10,'m^3/(mol*s)'), n=-0.9855, Ea=(1.09644,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_rad;C_pri_rad] for rate rule [Cd_rad;C_rad/H2/O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2]OCC1(O)[C](O)C1C(32248)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_NdNd;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['C[CH]C1(O)COC[C]1O(32249)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.95277e+09,'s^-1'), n=0.6, Ea=(85.772,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_secNd;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['C=CC(O)=C(O)COC(32250)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['CC1COCC(O)=C1O(31892)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;C_rad_out_2H;Cpri_rad_out_single] + [R6_SSSDS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R6_SSSDS;C_rad_out_2H;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2]OCC(=O)C(O)[CH]C(32251)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1713.16,'s^-1'), n=2.9, Ea=(127.847,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_Cs;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2]OCC(O)C(=O)[CH]C(5921)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1290.48,'s^-1'), n=2.90375, Ea=(139.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH2(19)', 'C[CH]C(O)=C(O)C[O](31659)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['CHCH3(T)(95)', '[CH2]OCC(O)=[C]O(32252)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]OCC(=O)C(O)=CC(32253)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.39e+06,'cm^3/(mol*s)'), n=2.09, Ea=(25.5224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2835 used for Od_CO-CdCs;HJ
Exact match found for rate rule [Od_CO-CdCs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]OC=C(O)C(O)=CC(32254)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.182e+10,'cm^3/(mol*s)'), n=0.859, Ea=(6.76971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-OsH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]OCC([O])C(O)=CC(32255)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.52488e+09,'s^-1'), n=1.21745, Ea=(162.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cd] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]O[CH]C(O)C(O)=CC(32256)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.96379e+08,'s^-1'), n=1.46437, Ea=(146.571,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_1H;Cs_H_out_Cd] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]OCC(O)C(O)=[C]C(32257)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.05e+09,'s^-1'), n=0.86, Ea=(139.746,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_NonDe] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['[CH2]C=C(O)C(O)CO[CH2](32258)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(800000,'s^-1'), n=1.81, Ea=(149.787,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_2H] for rate rule [R4H_SDS;C_rad_out_NDMustO;Cs_H_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[C]=C(O)[C](O)COC(32259)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['CC1COC[C](O)[C]1O(32260)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.46e+06,'s^-1'), n=1.02, Ea=(25.5224,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_secNd;radadd_intra_cs2H] for rate rule [R6_linear;doublebond_intra_secNd_HNd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['CC=C(O)C(O)=COC(32261)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['CC=C(O)C(=O)COC(32262)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radEndo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CH2(S)(23)', '[CH2]OC[C](O)C(=C)O(32263)'],
    products = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.97e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.324, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 4 used for carbene;Cd_pri
Exact match found for rate rule [carbene;Cd_pri]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -3.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]OCC(O)=C(O)[CH]C(5903)'],
    products = ['CC=C(O)C1(O)COC1(31895)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_OneDe/O]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

network(
    label = '5234',
    isomers = [
        '[CH2]OCC(O)=C(O)[CH]C(5903)',
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
    label = '5234',
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

