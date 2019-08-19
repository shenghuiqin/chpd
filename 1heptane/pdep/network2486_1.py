species(
    label = '[CH2]CC[CH]OO(2083)',
    structure = SMILES('[CH2]CC[CH]OO'),
    E0 = (161.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,326.383,1690.44],'cm^-1')),
        HinderedRotor(inertia=(0.00158631,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629709,'amu*angstrom^2'), symmetry=1, barrier=(47.5617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103571,'amu*angstrom^2'), symmetry=1, barrier=(7.81773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103423,'amu*angstrom^2'), symmetry=1, barrier=(7.8149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629421,'amu*angstrom^2'), symmetry=1, barrier=(47.5599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29161,0.0643247,-6.79004e-05,4.43562e-08,-1.26393e-11,19481.5,24.7429], Tmin=(100,'K'), Tmax=(828.046,'K')), NASAPolynomial(coeffs=[7.12676,0.0361379,-1.68416e-05,3.2495e-09,-2.28917e-13,18515.1,-2.30744], Tmin=(828.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(RCCJ)"""),
)

species(
    label = 'CH2CHOOH(64)',
    structure = SMILES('C=COO'),
    E0 = (-53.0705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.754187,'amu*angstrom^2'), symmetry=1, barrier=(17.3402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754176,'amu*angstrom^2'), symmetry=1, barrier=(17.34,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79821,0.0274377,-2.03468e-05,7.62127e-09,-1.12671e-12,-6315.53,9.11829], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.82519,0.0274167,-2.04456e-05,7.77399e-09,-1.18661e-12,-6325.27,8.96641], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), E0=(-53.0705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CH2CHOOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'C2H4(21)',
    structure = SMILES('C=C'),
    E0 = (42.0619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]CC=COO(14318)',
    structure = SMILES('[CH2]CC=COO'),
    E0 = (104.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.843123,'amu*angstrom^2'), symmetry=1, barrier=(19.3851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.844362,'amu*angstrom^2'), symmetry=1, barrier=(19.4135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843517,'amu*angstrom^2'), symmetry=1, barrier=(19.3941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00548084,'amu*angstrom^2'), symmetry=1, barrier=(7.51602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22354,0.0536522,-3.84767e-05,1.1081e-08,-2.90485e-13,12665.9,25.7347], Tmin=(100,'K'), Tmax=(1100.75,'K')), NASAPolynomial(coeffs=[13.1919,0.0207835,-8.16196e-06,1.48829e-09,-1.03014e-13,9387.5,-36.0776], Tmin=(1100.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJ)"""),
)

species(
    label = 'C=CC[CH]OO(1931)',
    structure = SMILES('C=CC[CH]OO'),
    E0 = (83.7382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,180,1544.02],'cm^-1')),
        HinderedRotor(inertia=(0.68936,'amu*angstrom^2'), symmetry=1, barrier=(15.8497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.685354,'amu*angstrom^2'), symmetry=1, barrier=(15.7576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.818608,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60568,0.0554494,-4.7278e-05,2.24922e-08,-4.60428e-12,10155.2,22.2921], Tmin=(100,'K'), Tmax=(1121.97,'K')), NASAPolynomial(coeffs=[8.60356,0.0305008,-1.39231e-05,2.67291e-09,-1.88084e-13,8584.96,-12.2739], Tmin=(1121.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.7382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOOH)"""),
)

species(
    label = 'C2H4(T)(94)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1436.54,1437.15,2688.96,2689.16],'cm^-1')),
        HinderedRotor(inertia=(0.0257549,'amu*angstrom^2'), symmetry=1, barrier=(17.2441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40736,0.0100312,6.40927e-06,-1.41291e-08,5.92671e-12,38288.2,6.11703], Tmin=(100,'K'), Tmax=(954.26,'K')), NASAPolynomial(coeffs=[5.52249,0.00856173,-2.90743e-06,5.02353e-10,-3.44572e-14,37547.8,-5.75276], Tmin=(954.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]OO(104)',
    structure = SMILES('[CH2][CH]OO'),
    E0 = (224.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00920734,'amu*angstrom^2'), symmetry=1, barrier=(3.53679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00921023,'amu*angstrom^2'), symmetry=1, barrier=(3.53685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43223,'amu*angstrom^2'), symmetry=1, barrier=(32.9297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52392,0.0279389,-1.79903e-05,3.58397e-09,4.58838e-13,27095.5,18.6054], Tmin=(100,'K'), Tmax=(1150.47,'K')), NASAPolynomial(coeffs=[9.41961,0.0107482,-4.42273e-06,8.47943e-10,-6.05179e-14,25059.9,-17.5802], Tmin=(1150.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(CJCOOH)"""),
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
    label = '[CH2]CCC=O(1026)',
    structure = SMILES('[CH2]CCC=O'),
    E0 = (-21.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.213965,'amu*angstrom^2'), symmetry=1, barrier=(4.91947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213194,'amu*angstrom^2'), symmetry=1, barrier=(4.90175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00471548,'amu*angstrom^2'), symmetry=1, barrier=(4.90985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49762,0.0360988,-1.83384e-05,4.20693e-09,-3.76081e-13,-2578.43,18.7989], Tmin=(100,'K'), Tmax=(2478.97,'K')), NASAPolynomial(coeffs=[13.6768,0.0180607,-7.42392e-06,1.27176e-09,-8.00786e-14,-8121.12,-45.2836], Tmin=(2478.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]COO(2082)',
    structure = SMILES('[CH2]C[CH]COO'),
    E0 = (173.035,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,359.746,2020.54],'cm^-1')),
        HinderedRotor(inertia=(0.00355485,'amu*angstrom^2'), symmetry=1, barrier=(10.3,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412618,'amu*angstrom^2'), symmetry=1, barrier=(37.9518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112283,'amu*angstrom^2'), symmetry=1, barrier=(10.2948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111948,'amu*angstrom^2'), symmetry=1, barrier=(10.298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411175,'amu*angstrom^2'), symmetry=1, barrier=(37.9524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33069,0.0634044,-6.62501e-05,4.25938e-08,-1.19396e-11,20903.3,26.7584], Tmin=(100,'K'), Tmax=(840.727,'K')), NASAPolynomial(coeffs=[7.18344,0.0355588,-1.65699e-05,3.20006e-09,-2.25637e-13,19919.2,-0.462476], Tmin=(840.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C[CH]OO(309)',
    structure = SMILES('C[CH]C[CH]OO'),
    E0 = (150.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,215.104,1603.04],'cm^-1')),
        HinderedRotor(inertia=(0.168488,'amu*angstrom^2'), symmetry=1, barrier=(5.53204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303367,'amu*angstrom^2'), symmetry=1, barrier=(5.53202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.7442,'amu*angstrom^2'), symmetry=1, barrier=(57.2689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168485,'amu*angstrom^2'), symmetry=1, barrier=(5.53203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74446,'amu*angstrom^2'), symmetry=1, barrier=(57.2689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26598,0.0662667,-8.68783e-05,7.58369e-08,-2.75572e-11,18181.7,25.7673], Tmin=(100,'K'), Tmax=(806.273,'K')), NASAPolynomial(coeffs=[4.55228,0.0397216,-1.84402e-05,3.49454e-09,-2.41229e-13,17984.7,12.6848], Tmin=(806.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]CCCO[O](272)',
    structure = SMILES('[CH2]CCCO[O]'),
    E0 = (124.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,492.5,1135,1000,3000,3100,440,815,1455,1000,211.819,1654.03],'cm^-1')),
        HinderedRotor(inertia=(0.262171,'amu*angstrom^2'), symmetry=1, barrier=(8.34708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00375289,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261522,'amu*angstrom^2'), symmetry=1, barrier=(8.3433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262083,'amu*angstrom^2'), symmetry=1, barrier=(8.34519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51474,0.0584856,-5.77493e-05,3.6598e-08,-1.02894e-11,15075,24.5481], Tmin=(100,'K'), Tmax=(833.732,'K')), NASAPolynomial(coeffs=[6.38667,0.0351122,-1.56986e-05,2.97451e-09,-2.07471e-13,14262.6,1.92965], Tmin=(833.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2][CH]CCOO(310)',
    structure = SMILES('[CH2][CH]CCOO'),
    E0 = (167.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1842.41,1843.45],'cm^-1')),
        HinderedRotor(inertia=(0.204753,'amu*angstrom^2'), symmetry=1, barrier=(9.22363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208519,'amu*angstrom^2'), symmetry=1, barrier=(9.26148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00383821,'amu*angstrom^2'), symmetry=1, barrier=(9.25006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12134,'amu*angstrom^2'), symmetry=1, barrier=(49.1916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07868,'amu*angstrom^2'), symmetry=1, barrier=(49.2764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33193,0.063861,-7.90196e-05,6.60931e-08,-2.3602e-11,20184.5,27.1691], Tmin=(100,'K'), Tmax=(784.709,'K')), NASAPolynomial(coeffs=[5.06473,0.038298,-1.76625e-05,3.35254e-09,-2.32269e-13,19799.9,11.3475], Tmin=(784.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'CC[CH][CH]OO(5724)',
    structure = SMILES('CC[CH][CH]OO'),
    E0 = (156.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,308.52,1927.19],'cm^-1')),
        HinderedRotor(inertia=(0.134124,'amu*angstrom^2'), symmetry=1, barrier=(9.0915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672494,'amu*angstrom^2'), symmetry=1, barrier=(45.2172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00176209,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134794,'amu*angstrom^2'), symmetry=1, barrier=(9.09234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.668937,'amu*angstrom^2'), symmetry=1, barrier=(45.2136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19155,0.0667438,-7.77602e-05,5.7783e-08,-1.86258e-11,18903.6,25.6148], Tmin=(100,'K'), Tmax=(740.163,'K')), NASAPolynomial(coeffs=[6.71054,0.0369166,-1.73105e-05,3.33349e-09,-2.339e-13,18086.7,0.649718], Tmin=(740.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOOH) + radical(CCJCOOH)"""),
)

species(
    label = 'CCC[CH]O[O](257)',
    structure = SMILES('CCC[CH]O[O]'),
    E0 = (107.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.119013,'amu*angstrom^2'), symmetry=1, barrier=(2.73635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120007,'amu*angstrom^2'), symmetry=1, barrier=(2.75919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121685,'amu*angstrom^2'), symmetry=1, barrier=(2.79779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481289,'amu*angstrom^2'), symmetry=1, barrier=(11.0658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80219,0.0553016,-3.67985e-05,-1.19632e-08,2.59165e-11,13057,21.9547], Tmin=(100,'K'), Tmax=(531.462,'K')), NASAPolynomial(coeffs=[5.78891,0.0367206,-1.66008e-05,3.14923e-09,-2.19362e-13,12471.9,3.72319], Tmin=(531.462,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = 'CCC=COO(5725)',
    structure = SMILES('CCC=COO'),
    E0 = (-100.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2109,0.052651,-2.72039e-05,-2.05256e-09,4.43529e-12,-12017.9,23.9923], Tmin=(100,'K'), Tmax=(1043.51,'K')), NASAPolynomial(coeffs=[13.1657,0.0234554,-9.1411e-06,1.67948e-09,-1.17579e-13,-15418.3,-38.5302], Tmin=(1043.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH)"""),
)

species(
    label = 'C=CCCOO(306)',
    structure = SMILES('C=CCCOO'),
    E0 = (-104.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655,323.996,323.996],'cm^-1')),
        HinderedRotor(inertia=(0.111179,'amu*angstrom^2'), symmetry=1, barrier=(8.28185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331823,'amu*angstrom^2'), symmetry=1, barrier=(24.718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111178,'amu*angstrom^2'), symmetry=1, barrier=(8.28185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.419395,'amu*angstrom^2'), symmetry=1, barrier=(31.2413,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37317,0.05546,-4.0154e-05,1.50002e-08,-2.32081e-12,-12513.7,22.9715], Tmin=(100,'K'), Tmax=(1476.88,'K')), NASAPolynomial(coeffs=[11.7098,0.0274642,-1.17202e-05,2.16524e-09,-1.48175e-13,-15566.9,-30.9276], Tmin=(1476.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-104.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]CC([CH2])OO(2134)',
    structure = SMILES('[CH2]CC([CH2])OO'),
    E0 = (173.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1380,1390,370,380,2900,435,259.801],'cm^-1')),
        HinderedRotor(inertia=(0.734078,'amu*angstrom^2'), symmetry=1, barrier=(16.8779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00276332,'amu*angstrom^2'), symmetry=1, barrier=(31.3749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10436,'amu*angstrom^2'), symmetry=1, barrier=(48.3833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00276351,'amu*angstrom^2'), symmetry=1, barrier=(31.3771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36447,'amu*angstrom^2'), symmetry=1, barrier=(31.3717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1683,0.0629297,-5.69014e-05,2.78651e-08,-5.66751e-12,21021.5,26.6726], Tmin=(100,'K'), Tmax=(1157.5,'K')), NASAPolynomial(coeffs=[11.0409,0.028813,-1.269e-05,2.40154e-09,-1.6784e-13,18736,-22.4009], Tmin=(1157.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CJCOOH) + radical(RCCJ)"""),
)

species(
    label = 'OOC1CCC1(13167)',
    structure = SMILES('OOC1CCC1'),
    E0 = (-87.0924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80286,0.0329153,3.44181e-05,-6.5822e-08,2.65312e-11,-10382,20.1791], Tmin=(100,'K'), Tmax=(991.679,'K')), NASAPolynomial(coeffs=[13.1921,0.0244671,-9.51187e-06,1.83336e-09,-1.35215e-13,-14484.4,-43.9673], Tmin=(991.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.0924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]CCC([O])O(4408)',
    structure = SMILES('[CH2]CCC([O])O'),
    E0 = (-73.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,211.106,1680.17,1680.21],'cm^-1')),
        HinderedRotor(inertia=(0.00377893,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00377916,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250357,'amu*angstrom^2'), symmetry=1, barrier=(7.92332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250121,'amu*angstrom^2'), symmetry=1, barrier=(7.92218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.1051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33272,0.0453153,6.15728e-06,-1.07263e-07,1.06034e-10,-8778.37,24.2952], Tmin=(100,'K'), Tmax=(432.868,'K')), NASAPolynomial(coeffs=[4.54075,0.038633,-1.82347e-05,3.53317e-09,-2.49027e-13,-9098.08,14.0067], Tmin=(432.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH2](128)',
    structure = SMILES('[CH2]C[CH2]'),
    E0 = (290.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.000703322,'amu*angstrom^2'), symmetry=1, barrier=(7.98556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0167612,'amu*angstrom^2'), symmetry=1, barrier=(29.9599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.12863,0.0145638,1.31003e-05,-2.16419e-08,7.8594e-12,34974.2,14.113], Tmin=(100,'K'), Tmax=(1043.7,'K')), NASAPolynomial(coeffs=[5.43609,0.0169627,-6.50459e-06,1.20126e-09,-8.4342e-14,33880.3,-0.0512239], Tmin=(1043.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH]OO(168)',
    structure = SMILES('[CH]OO'),
    E0 = (320.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00655552,'amu*angstrom^2'), symmetry=1, barrier=(3.39708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148635,'amu*angstrom^2'), symmetry=1, barrier=(16.5265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06576,0.0188176,-2.02137e-05,1.03271e-08,-2.01695e-12,38527.3,11.7484], Tmin=(100,'K'), Tmax=(1262.94,'K')), NASAPolynomial(coeffs=[8.15329,0.00270404,-1.0753e-06,2.24439e-10,-1.70827e-14,37242.2,-13.9836], Tmin=(1262.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CH2_triplet)"""),
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
    label = '[CH2]C[CH]OO(468)',
    structure = SMILES('[CH2]C[CH]OO'),
    E0 = (192.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2231.58],'cm^-1')),
        HinderedRotor(inertia=(0.046424,'amu*angstrom^2'), symmetry=1, barrier=(9.87208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218572,'amu*angstrom^2'), symmetry=1, barrier=(3.32054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154635,'amu*angstrom^2'), symmetry=1, barrier=(3.55536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99023,'amu*angstrom^2'), symmetry=1, barrier=(45.7599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09557,0.0449449,-4.38708e-05,2.57348e-08,-6.52716e-12,23271.7,20.8749], Tmin=(100,'K'), Tmax=(925.066,'K')), NASAPolynomial(coeffs=[6.76412,0.0247575,-1.1136e-05,2.14313e-09,-1.51319e-13,22407.9,-1.28434], Tmin=(925.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCsJOOH) + radical(RCCJ)"""),
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
    E0 = (161.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (322.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (305.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (296.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (288.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (161.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (286.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (313.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (290.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (303.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (279.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (234.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (542.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (224.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (224.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (332.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (169.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (253.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (610.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (574.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['CH2CHOOH(64)', 'C2H4(21)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', '[CH2]CC=COO(14318)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=CC[CH]OO(1931)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2CHOOH(64)', 'C2H4(T)(94)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.0133609,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]OO(104)', 'C2H4(21)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00310016,'m^3/(mol*s)'), n=2.49247, Ea=(21.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OH(5)', '[CH2]CCC=O(1026)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(154.687,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 150.6 to 154.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['[CH2]C[CH]COO(2082)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.4e-19,'s^-1'), n=8.84, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 342 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['C[CH]C[CH]OO(309)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]CCCO[O](272)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.44e+09,'s^-1'), n=1.17, Ea=(165.937,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 242 used for R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['[CH2][CH]CCOO(310)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['CC[CH][CH]OO(5724)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['CCC[CH]O[O](257)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]OO(104)', 'C2H4(T)(94)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['CCC=COO(5725)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['C=CCCOO(306)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC([CH2])OO(2134)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['OOC1CCC1(13167)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC[CH]OO(2083)'],
    products = ['[CH2]CCC([O])O(4408)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.263e+10,'s^-1'), n=0, Ea=(92.7451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnOOH;C_rad_out_H/NonDeC] for rate rule [ROOH;C_rad_out_H/NonDeC]
Euclidian distance = 1.0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C[CH2](128)', '[CH]OO(168)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(19)', '[CH2]C[CH]OO(468)'],
    products = ['[CH2]CC[CH]OO(2083)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2486',
    isomers = [
        '[CH2]CC[CH]OO(2083)',
    ],
    reactants = [
        ('CH2CHOOH(64)', 'C2H4(21)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2486',
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

