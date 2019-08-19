species(
    label = 'C#C[CH]CC[C]=CC(26584)',
    structure = SMILES('C#C[CH]CC[C]=CC'),
    E0 = (554.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346599,0.0809579,-7.75396e-05,4.41698e-08,-1.05714e-11,66767.8,30.364], Tmin=(100,'K'), Tmax=(996.474,'K')), NASAPolynomial(coeffs=[10.5847,0.0398595,-1.5672e-05,2.77757e-09,-1.86406e-13,64727.5,-18.9922], Tmin=(996.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Sec_Propargyl)"""),
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
    label = '[CH]=C1[CH]CCC1=CC(28097)',
    structure = SMILES('[CH]C1=CCCC1=CC'),
    E0 = (394.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.811404,0.0529502,2.62969e-05,-6.56743e-08,2.70874e-11,47520.2,26.1971], Tmin=(100,'K'), Tmax=(1000.14,'K')), NASAPolynomial(coeffs=[14.1263,0.0403174,-1.56765e-05,2.91134e-09,-2.0695e-13,42825.3,-48.1979], Tmin=(1000.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(3-Methylenecyclopentene) + radical(AllylJ2_triplet)"""),
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
    label = 'C#CC=CC[C]=CC(28098)',
    structure = SMILES('C#CC=CC[C]=CC'),
    E0 = (519.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.845525,'amu*angstrom^2'), symmetry=1, barrier=(19.4403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845679,'amu*angstrom^2'), symmetry=1, barrier=(19.4438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845685,'amu*angstrom^2'), symmetry=1, barrier=(19.444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845655,'amu*angstrom^2'), symmetry=1, barrier=(19.4433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.36666,0.0726553,-5.51667e-05,2.14891e-08,-3.39903e-12,62572.8,29.6821], Tmin=(100,'K'), Tmax=(1481.68,'K')), NASAPolynomial(coeffs=[16.0159,0.0304082,-1.23975e-05,2.24559e-09,-1.5214e-13,57935.3,-51.9695], Tmin=(1481.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]CC=C=CC(28099)',
    structure = SMILES('C#C[CH]CC=C=CC'),
    E0 = (480.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,206.397,206.409],'cm^-1')),
        HinderedRotor(inertia=(0.40362,'amu*angstrom^2'), symmetry=1, barrier=(12.2049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.60923,'amu*angstrom^2'), symmetry=1, barrier=(78.8989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403702,'amu*angstrom^2'), symmetry=1, barrier=(12.2049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.60963,'amu*angstrom^2'), symmetry=1, barrier=(78.8992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.364168,0.0774807,-7.10577e-05,3.70863e-08,-7.98397e-12,57931.8,28.4379], Tmin=(100,'K'), Tmax=(1111.54,'K')), NASAPolynomial(coeffs=[12.3237,0.0344434,-1.29802e-05,2.25342e-09,-1.49646e-13,55273.1,-30.5245], Tmin=(1111.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[CH]CCC#CC(28100)',
    structure = SMILES('C#C[CH]CCC#CC'),
    E0 = (476.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,750,770,3400,2100,2100,2175,2250,500,525,550,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631822,0.0717014,-5.16779e-05,9.32997e-09,5.89761e-12,57392.1,28.3344], Tmin=(100,'K'), Tmax=(801.273,'K')), NASAPolynomial(coeffs=[12.0883,0.0318869,-9.67434e-06,1.44786e-09,-8.75999e-14,54998.4,-27.8791], Tmin=(801.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
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
    label = 'C#C[CH]CCC#C(28101)',
    structure = SMILES('C#C[CH]CCC#C'),
    E0 = (518.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2100,2250,500,550,750,756.667,763.333,770,3350,3450,2000,2200,3025,407.5,1350,352.5,180,1773.89],'cm^-1')),
        HinderedRotor(inertia=(3.87909,'amu*angstrom^2'), symmetry=1, barrier=(89.188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.87839,'amu*angstrom^2'), symmetry=1, barrier=(89.1718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.87904,'amu*angstrom^2'), symmetry=1, barrier=(89.1869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.87884,'amu*angstrom^2'), symmetry=1, barrier=(89.1822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00545,0.0640795,-5.44554e-05,1.68941e-08,2.39826e-12,62451.3,23.5704], Tmin=(100,'K'), Tmax=(799.795,'K')), NASAPolynomial(coeffs=[12.1558,0.0235118,-6.87589e-06,9.94077e-10,-5.85116e-14,60181.6,-30.7715], Tmin=(799.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C[C]=CC(27995)',
    structure = SMILES('C#CC[CH]C[C]=CC'),
    E0 = (602.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568446,0.0810105,-9.46252e-05,7.51152e-08,-2.54608e-11,72557.9,32.2917], Tmin=(100,'K'), Tmax=(818.037,'K')), NASAPolynomial(coeffs=[5.72931,0.0482894,-2.08999e-05,3.84578e-09,-2.6158e-13,71964,9.96092], Tmin=(818.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJCC) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C[CH]C=CC(27671)',
    structure = SMILES('C#C[CH]C[CH]C=CC'),
    E0 = (457.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.228032,0.0752685,-5.82293e-05,2.49434e-08,-4.411e-12,55145.1,29.2069], Tmin=(100,'K'), Tmax=(1335.68,'K')), NASAPolynomial(coeffs=[13.8286,0.0345384,-1.24883e-05,2.11294e-09,-1.37806e-13,51511.9,-40.3445], Tmin=(1335.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[CH]CCC=[C]C(27667)',
    structure = SMILES('C#C[CH]CCC=[C]C'),
    E0 = (554.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346599,0.0809579,-7.75396e-05,4.41698e-08,-1.05714e-11,66767.8,30.364], Tmin=(100,'K'), Tmax=(996.474,'K')), NASAPolynomial(coeffs=[10.5847,0.0398595,-1.5672e-05,2.77757e-09,-1.86406e-13,64727.5,-18.9922], Tmin=(996.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CCC[CH][C]=CC(28102)',
    structure = SMILES('C#CCC[CH][C]=CC'),
    E0 = (548.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.584114,0.0744971,-5.88372e-05,2.62785e-08,-5.00292e-12,66147.1,28.5618], Tmin=(100,'K'), Tmax=(1213.02,'K')), NASAPolynomial(coeffs=[10.9733,0.040238,-1.64726e-05,2.99499e-09,-2.04213e-13,63626.6,-23.566], Tmin=(1213.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[CH][CH]CC=CC(27673)',
    structure = SMILES('[CH]=C=C[CH]CC=CC'),
    E0 = (463.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194458,0.0734606,-5.22973e-05,1.95026e-08,-2.96779e-12,55932.8,30.9012], Tmin=(100,'K'), Tmax=(1537.53,'K')), NASAPolynomial(coeffs=[16.0096,0.032316,-1.21566e-05,2.09754e-09,-1.37726e-13,51069.7,-52.2006], Tmin=(1537.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_S)"""),
)

species(
    label = 'C#C[CH]CC[CH]C=C(26564)',
    structure = SMILES('C#C[CH]CCC=C[CH2]'),
    E0 = (467.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3593.44,'J/mol'), sigma=(6.30673,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.29 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.158128,0.0765436,-5.38631e-05,1.44172e-08,7.8694e-13,56397.5,30.8415], Tmin=(100,'K'), Tmax=(986.639,'K')), NASAPolynomial(coeffs=[14.6957,0.0332245,-1.175e-05,2.00637e-09,-1.33839e-13,52768.6,-42.951], Tmin=(986.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Allyl_P)"""),
)

species(
    label = '[C]#CCCC[C]=CC(28103)',
    structure = SMILES('[C]#CCCC[C]=CC'),
    E0 = (744.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2175,525,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.343392,0.0839569,-9.2539e-05,6.43861e-08,-1.90129e-11,89729.1,31.2811], Tmin=(100,'K'), Tmax=(813.456,'K')), NASAPolynomial(coeffs=[8.5746,0.0434834,-1.79096e-05,3.22633e-09,-2.17401e-13,88389.9,-6.73059], Tmin=(813.456,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Acetyl)"""),
)

species(
    label = 'C#CCCC[C]=[C]C(28104)',
    structure = SMILES('C#CCCC[C]=[C]C'),
    E0 = (645.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,750,770,3400,2100,2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376547,0.0838892,-9.04742e-05,6.07425e-08,-1.74078e-11,77784.2,30.8988], Tmin=(100,'K'), Tmax=(833.8,'K')), NASAPolynomial(coeffs=[8.67484,0.0440794,-1.88563e-05,3.47992e-09,-2.38488e-13,76400.4,-7.62736], Tmin=(833.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[C]#C[CH]CCC=CC(27674)',
    structure = SMILES('[C]#C[CH]CCC=CC'),
    E0 = (653.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32405,0.0808561,-7.88031e-05,4.64851e-08,-1.14853e-11,78712.3,30.7112], Tmin=(100,'K'), Tmax=(970.374,'K')), NASAPolynomial(coeffs=[10.3583,0.0394935,-1.48648e-05,2.55815e-09,-1.68232e-13,76764.9,-17.3966], Tmin=(970.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CCCC[C]=C[CH2](27670)',
    structure = SMILES('C#CCCC[C]=C[CH2]'),
    E0 = (559.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.208446,0.0793054,-6.64294e-05,3.05949e-08,-5.8189e-12,67412.7,31.2974], Tmin=(100,'K'), Tmax=(1243.31,'K')), NASAPolynomial(coeffs=[13.711,0.035864,-1.40181e-05,2.49123e-09,-1.67807e-13,64055.2,-36.7847], Tmin=(1243.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=[C]CCC1[C]=C1(28105)',
    structure = SMILES('CC=[C]CCC1[C]=C1'),
    E0 = (727.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637542,0.0790637,-7.39885e-05,4.20344e-08,-1.04886e-11,87621.7,27.6956], Tmin=(100,'K'), Tmax=(932.137,'K')), NASAPolynomial(coeffs=[8.34579,0.0459858,-2.07591e-05,3.96448e-09,-2.78211e-13,86184.7,-8.95048], Tmin=(932.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Cds_S) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = 'CC=C1C=[C][CH]CC1(28106)',
    structure = SMILES('CC=C1[CH][C]=CCC1'),
    E0 = (369.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73933,0.0239625,0.000104374,-1.44144e-07,5.3867e-11,44542.1,22.3173], Tmin=(100,'K'), Tmax=(987.575,'K')), NASAPolynomial(coeffs=[14.4487,0.037447,-1.4776e-05,2.89494e-09,-2.16226e-13,38864,-54.8777], Tmin=(987.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'C#CC=CCC=CC(27677)',
    structure = SMILES('C#CC=CCC=CC'),
    E0 = (281.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.346027,0.068757,-3.32281e-05,-3.41088e-09,5.38566e-12,33971.8,29.2787], Tmin=(100,'K'), Tmax=(1079.74,'K')), NASAPolynomial(coeffs=[15.5741,0.0334471,-1.34927e-05,2.5058e-09,-1.75548e-13,29453.2,-51.0528], Tmin=(1079.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CCCC=C=CC(28107)',
    structure = SMILES('C#CCCC=C=CC'),
    E0 = (334.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394085,0.0762943,-6.12126e-05,2.72121e-08,-5.05596e-12,40346,28.4916], Tmin=(100,'K'), Tmax=(1258.63,'K')), NASAPolynomial(coeffs=[12.6304,0.0374072,-1.48686e-05,2.66508e-09,-1.80274e-13,37265.8,-33.3562], Tmin=(1258.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH)"""),
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
    label = 'C#C[CH]CC[C]=C(28108)',
    structure = SMILES('C#C[CH]CC[C]=C'),
    E0 = (590.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,327.222,327.266,327.275],'cm^-1')),
        HinderedRotor(inertia=(0.00157462,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0829663,'amu*angstrom^2'), symmetry=1, barrier=(6.30518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06129,'amu*angstrom^2'), symmetry=1, barrier=(80.6733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06158,'amu*angstrom^2'), symmetry=1, barrier=(80.6717,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844363,0.0672034,-6.36282e-05,3.47187e-08,-7.78691e-12,71085.1,26.7662], Tmin=(100,'K'), Tmax=(1072.93,'K')), NASAPolynomial(coeffs=[11.0084,0.0293116,-1.06551e-05,1.80448e-09,-1.17833e-13,68904,-22.9851], Tmin=(1072.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC([CH2])C[C]=CC(26585)',
    structure = SMILES('C#CC([CH2])C[C]=CC'),
    E0 = (603.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,750,770,3400,2100,218.855,219.126],'cm^-1')),
        HinderedRotor(inertia=(0.00350253,'amu*angstrom^2'), symmetry=1, barrier=(0.119786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161609,'amu*angstrom^2'), symmetry=1, barrier=(5.50317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.9251,'amu*angstrom^2'), symmetry=1, barrier=(101.848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170685,'amu*angstrom^2'), symmetry=1, barrier=(5.50478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.01599,'amu*angstrom^2'), symmetry=1, barrier=(102.07,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3619.53,'J/mol'), sigma=(6.34117,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=565.36 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301188,0.0802179,-7.47905e-05,4.08001e-08,-9.26938e-12,72754.1,32.1875], Tmin=(100,'K'), Tmax=(1051.89,'K')), NASAPolynomial(coeffs=[11.4468,0.0378349,-1.43523e-05,2.4957e-09,-1.65688e-13,70409.3,-22.1473], Tmin=(1051.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]CC(=C)[CH]C(26576)',
    structure = SMILES('C#C[CH]CC([CH2])=CC'),
    E0 = (452.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0773253,0.0821768,-7.17924e-05,3.45969e-08,-6.76576e-12,54568.1,30.1996], Tmin=(100,'K'), Tmax=(1230.35,'K')), NASAPolynomial(coeffs=[15.1429,0.032694,-1.14643e-05,1.90791e-09,-1.23502e-13,50822.9,-46.3841], Tmin=(1230.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC1CCC1=CC(28109)',
    structure = SMILES('C#CC1CCC1=CC'),
    E0 = (277.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.7388,-0.000635209,0.000116702,-1.17013e-07,2.68619e-11,33013.3,-20.3185], Tmin=(100,'K'), Tmax=(1749.67,'K')), NASAPolynomial(coeffs=[82.4395,0.0271591,-7.14814e-05,1.73124e-08,-1.27893e-12,-21422.2,-490.203], Tmin=(1749.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane)"""),
)

species(
    label = 'C3H2(81)',
    structure = SMILES('[CH]C#C'),
    E0 = (530.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.621997,'amu*angstrom^2'), symmetry=1, barrier=(14.3009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6874,0.0353211,-7.31245e-05,7.03801e-08,-2.45313e-11,63833.3,8.93418], Tmin=(100,'K'), Tmax=(908.454,'K')), NASAPolynomial(coeffs=[4.78002,0.0100611,-4.92178e-06,8.868e-10,-5.66859e-14,64115.2,2.68365], Tmin=(908.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""HCCCH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[C]=CC(24198)',
    structure = SMILES('[CH2]C[C]=CC'),
    E0 = (390.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,218.707],'cm^-1')),
        HinderedRotor(inertia=(0.137812,'amu*angstrom^2'), symmetry=1, barrier=(4.67652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137786,'amu*angstrom^2'), symmetry=1, barrier=(4.67628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137779,'amu*angstrom^2'), symmetry=1, barrier=(4.67636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04972,0.0405227,-2.17348e-05,5.57612e-09,-5.73511e-13,47043.1,21.1549], Tmin=(100,'K'), Tmax=(2150.78,'K')), NASAPolynomial(coeffs=[12.3921,0.021288,-8.32003e-06,1.41801e-09,-9.01841e-14,42594.3,-36.6615], Tmin=(2150.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[C]=CC(24199)',
    structure = SMILES('[C]=CC'),
    E0 = (564.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.40488,'amu*angstrom^2'), symmetry=1, barrier=(9.30899,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28451,0.0144722,-2.9151e-06,-2.04635e-09,7.91551e-13,67868.8,10.0416], Tmin=(100,'K'), Tmax=(1455.86,'K')), NASAPolynomial(coeffs=[5.67821,0.0121697,-4.94658e-06,9.00486e-10,-6.07653e-14,66718.9,-3.96135], Tmin=(1455.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C=CC[CH2](26618)',
    structure = SMILES('C#C[CH]C[CH2]'),
    E0 = (477.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0141456,'amu*angstrom^2'), symmetry=1, barrier=(6.90141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.94797,'amu*angstrom^2'), symmetry=1, barrier=(67.7796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0405762,'amu*angstrom^2'), symmetry=1, barrier=(67.6916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76795,0.0430915,-3.64272e-05,1.77221e-08,-3.47093e-12,57486.3,19.7996], Tmin=(100,'K'), Tmax=(1346.56,'K')), NASAPolynomial(coeffs=[9.46877,0.0175599,-5.02768e-06,7.11764e-10,-4.08662e-14,55653.2,-18.7496], Tmin=(1346.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]=[C]C1CCC1=CC(28110)',
    structure = SMILES('[CH]=[C]C1CCC1=CC'),
    E0 = (620.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.4585,-0.00407591,0.000119164,-1.16947e-07,2.65163e-11,74276.4,-19.3863], Tmin=(100,'K'), Tmax=(1769.07,'K')), NASAPolynomial(coeffs=[84.8065,0.0254037,-7.14488e-05,1.72963e-08,-1.27453e-12,17760.4,-501.475], Tmin=(1769.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=C=[C]CC[C]=CC(28111)',
    structure = SMILES('[CH2]C#CCC[C]=CC'),
    E0 = (545.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.786799,0.0711967,-5.35741e-05,2.30143e-08,-4.29003e-12,65765.1,31.1734], Tmin=(100,'K'), Tmax=(1218.2,'K')), NASAPolynomial(coeffs=[9.75861,0.0417377,-1.73008e-05,3.16369e-09,-2.163e-13,63579.2,-13.8813], Tmin=(1218.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C=[C]CCC=CC(27690)',
    structure = SMILES('[CH]=C=[C]CCC=CC'),
    E0 = (560.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,210.536,210.54],'cm^-1')),
        HinderedRotor(inertia=(0.477152,'amu*angstrom^2'), symmetry=1, barrier=(15.0099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364056,'amu*angstrom^2'), symmetry=1, barrier=(11.4515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364063,'amu*angstrom^2'), symmetry=1, barrier=(11.4515,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364058,'amu*angstrom^2'), symmetry=1, barrier=(11.4514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.526638,0.0768565,-6.45811e-05,3.09018e-08,-6.28367e-12,67545.8,31.2763], Tmin=(100,'K'), Tmax=(1145.37,'K')), NASAPolynomial(coeffs=[10.9767,0.0403615,-1.67866e-05,3.08297e-09,-2.11667e-13,65152,-20.5579], Tmin=(1145.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(560.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C=CC[C]=CC(27958)',
    structure = SMILES('C=[C]C=CC[C]=CC'),
    E0 = (541.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.713501,'amu*angstrom^2'), symmetry=1, barrier=(16.4048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714239,'amu*angstrom^2'), symmetry=1, barrier=(16.4218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714825,'amu*angstrom^2'), symmetry=1, barrier=(16.4352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714532,'amu*angstrom^2'), symmetry=1, barrier=(16.4285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283233,0.0780811,-6.41802e-05,2.88981e-08,-5.38961e-12,65244.1,30.6758], Tmin=(100,'K'), Tmax=(1261.83,'K')), NASAPolynomial(coeffs=[13.5051,0.036168,-1.43563e-05,2.5746e-09,-1.74308e-13,61907.4,-36.1871], Tmin=(1261.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC=C=CC(27658)',
    structure = SMILES('[CH2][C]=CCC=C=CC'),
    E0 = (546.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.328597,'amu*angstrom^2'), symmetry=1, barrier=(7.55509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243112,'amu*angstrom^2'), symmetry=1, barrier=(17.4589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759703,'amu*angstrom^2'), symmetry=1, barrier=(17.4671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75971,'amu*angstrom^2'), symmetry=1, barrier=(17.4672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503535,0.0741618,-5.45351e-05,2.10444e-08,-3.38084e-12,65870.7,30.5981], Tmin=(100,'K'), Tmax=(1425.5,'K')), NASAPolynomial(coeffs=[13.6464,0.0372832,-1.57298e-05,2.89665e-09,-1.982e-13,62123.6,-37.4685], Tmin=(1425.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=CCC[C]=[C]C(28112)',
    structure = SMILES('C=C=CCC[C]=[C]C'),
    E0 = (643.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,246.83,246.83,246.83],'cm^-1')),
        HinderedRotor(inertia=(0.247481,'amu*angstrom^2'), symmetry=1, barrier=(10.6995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247481,'amu*angstrom^2'), symmetry=1, barrier=(10.6995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24748,'amu*angstrom^2'), symmetry=1, barrier=(10.6995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24748,'amu*angstrom^2'), symmetry=1, barrier=(10.6995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.724632,0.0773667,-6.9723e-05,3.80762e-08,-9.22392e-12,77561.1,30.4288], Tmin=(100,'K'), Tmax=(951.371,'K')), NASAPolynomial(coeffs=[8.04726,0.0465794,-2.1182e-05,4.06177e-09,-2.85759e-13,76167.8,-4.53372], Tmin=(951.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CCC=C=C(27691)',
    structure = SMILES('[CH2]C=[C]CCC=C=C'),
    E0 = (557.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,207.071,207.843,937.392],'cm^-1')),
        HinderedRotor(inertia=(0.255796,'amu*angstrom^2'), symmetry=1, barrier=(7.76236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.634246,'amu*angstrom^2'), symmetry=1, barrier=(19.5886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.641376,'amu*angstrom^2'), symmetry=1, barrier=(19.5882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.642606,'amu*angstrom^2'), symmetry=1, barrier=(19.5885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.290965,0.0759248,-5.6612e-05,2.18964e-08,-3.47566e-12,67201,31.7779], Tmin=(100,'K'), Tmax=(1462.23,'K')), NASAPolynomial(coeffs=[15.3725,0.0346687,-1.42904e-05,2.60102e-09,-1.76714e-13,62790.4,-46.7124], Tmin=(1462.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C=CCC=C=CC(27660)',
    structure = SMILES('C=C=CCC=C=CC'),
    E0 = (333.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449316,0.0731152,-5.19279e-05,1.90028e-08,-2.86436e-12,40288.2,29.0972], Tmin=(100,'K'), Tmax=(1522.24,'K')), NASAPolynomial(coeffs=[14.861,0.0352455,-1.46114e-05,2.65992e-09,-1.80324e-13,35900.6,-46.4862], Tmin=(1522.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    E0 = (554.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (679.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (733.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (700.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (706.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (638.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (613.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (684.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (727.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (711.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (788.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (755.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (696.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (715.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (892.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (795.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (822.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (641.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (779.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (585.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (647.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (647.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1009.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (763.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (724.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (561.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (920.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1041.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (628.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (686.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (706.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (598.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (587.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (825.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (769.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (571.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['[CH]=C1[CH]CCC1=CC(28097)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.03297e+11,'s^-1'), n=0.45637, Ea=(125.756,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;triplebond_intra_H;radadd_intra] for rate rule [R6;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C#CC=CC[C]=CC(28098)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.31e+08,'cm^3/(mol*s)'), n=1.64, Ea=(2.76144,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2560 used for Cds-CsH_Cds-CtH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CtH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#C[CH]CC=C=CC(28099)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#C[CH]CCC#CC(28100)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0132398,'m^3/(mol*s)'), n=2.333, Ea=(2.89605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CtH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', 'C#C[CH]CCC#C(28101)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC[CH]C[C]=CC(27995)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Ct]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#C[CH]C[CH]C=CC(27671)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#C[CH]CCC=[C]C(27667)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#CCC[CH][C]=CC(28102)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.14401e+08,'s^-1'), n=1.25033, Ea=(201.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_H/OneDe;XH_out] for rate rule [R3H_SS_Cs;C_rad_out_H/Ct;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#C[CH][CH]CC=CC(27673)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#C[CH]CC[CH]C=C(26564)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[C]#CCCC[C]=CC(28103)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CCCC[C]=[C]C(28104)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R5HJ_1;Cd_rad_out_Cs;Cs_H_out_H/Ct]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]#C[CH]CCC=CC(27674)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.48073e+07,'s^-1'), n=1.57892, Ea=(169.446,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_doubleC] for rate rule [R6HJ_2;Ct_rad_out;Cd_H_out_doubleC]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CCCC[C]=C[CH2](27670)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(54838.5,'s^-1'), n=1.94766, Ea=(82.2091,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R6HJ_2;C_rad_out_2H;Cs_H_out_H/Ct]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['CC=[C]CCC1[C]=C1(28105)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_T;triplebond_intra_H;radadd_intra_cs] for rate rule [R3_T;triplebond_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['CC=C1C=[C][CH]CC1(28106)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#CC=CCC=CC(27677)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#CCCC=C=CC(28107)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(S)(23)', 'C#C[CH]CC[C]=C(28108)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C[C]=CC(26585)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#C[CH]CC(=C)[CH]C(26576)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C#CC1CCC1=CC(28109)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C3H2(81)', '[CH2]C[C]=CC(24198)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[C]=CC(24199)', '[CH]=C=CC[CH2](26618)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['[CH]=[C]C1CCC1=CC(28110)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.48517e+07,'s^-1'), n=1.03851, Ea=(74.6464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C=C=[C]CC[C]=CC(28111)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(443913,'s^-1'), n=2.07262, Ea=(132.134,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;Cd_H_out_singleNd] + [R3H;Cd_rad_out;Cd_H_out_singleNd] + [R3H;Cd_rad_out_singleH;XH_out] for rate rule [R3H;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C=[C]CCC=CC(27690)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cd_H_out_doubleC] for rate rule [R4H_SSS;Cd_rad_out_double;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C=[C]C=CC[C]=CC(27958)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['[CH2][C]=CCC=C=CC(27658)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C=CCC[C]=[C]C(28112)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cd_H_out_singleH] for rate rule [R7HJ_1;Cd_rad_out_Cs;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C=[C]CCC=C=C(27691)'],
    products = ['C#C[CH]CC[C]=CC(26584)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(291.349,'s^-1'), n=2.95077, Ea=(212.39,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_singleH] for rate rule [R8Hall;C_rad_out_2H;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#C[CH]CC[C]=CC(26584)'],
    products = ['C=C=CCC=C=CC(27660)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

network(
    label = '4556',
    isomers = [
        'C#C[CH]CC[C]=CC(26584)',
    ],
    reactants = [
        ('CH3CHCCH2(18175)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4556',
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

