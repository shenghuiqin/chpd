species(
    label = 'C#CC(C)O[O](26146)',
    structure = SMILES('C#CC(C)O[O]'),
    E0 = (200.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,492.5,1135,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.72109,'amu*angstrom^2'), symmetry=1, barrier=(16.5793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722427,'amu*angstrom^2'), symmetry=1, barrier=(16.61,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19122,'amu*angstrom^2'), symmetry=1, barrier=(50.3804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3669.66,'J/mol'), sigma=(6.11465,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=573.19 K, Pc=36.42 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39394,0.052523,-5.20004e-05,2.72774e-08,-5.66058e-12,24168.1,21.8356], Tmin=(100,'K'), Tmax=(1177.97,'K')), NASAPolynomial(coeffs=[12.1857,0.0158771,-5.33574e-06,8.67386e-10,-5.552e-14,21625.7,-31.9959], Tmin=(1177.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = 'CC1[C]=COO1(26431)',
    structure = SMILES('CC1[C]=COO1'),
    E0 = (201.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3594.27,'J/mol'), sigma=(5.99159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.42 K, Pc=37.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14932,0.0249527,3.74048e-05,-7.02522e-08,2.94238e-11,24276.8,20.3516], Tmin=(100,'K'), Tmax=(961.705,'K')), NASAPolynomial(coeffs=[14.6463,0.012195,-3.87025e-06,7.66602e-10,-6.16998e-14,20059.4,-48.8808], Tmin=(961.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(887.157,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36759e-06,-4.93128e-09,1.45955e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.71,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00335e-07,1.5903e-10,-1.14891e-14,-1048.45,6.08303], Tmin=(1087.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CH3CHCCH(24127)',
    structure = SMILES('C#C[CH]C'),
    E0 = (305.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525],'cm^-1')),
        HinderedRotor(inertia=(0.473568,'amu*angstrom^2'), symmetry=1, barrier=(10.8883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0735085,'amu*angstrom^2'), symmetry=1, barrier=(26.973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2908.77,'J/mol'), sigma=(5.16769,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=454.34 K, Pc=47.83 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6542,0.0250419,-5.9594e-06,-7.9258e-09,4.31635e-12,36808.5,13.9197], Tmin=(100,'K'), Tmax=(1027.31,'K')), NASAPolynomial(coeffs=[7.5854,0.0154004,-5.83859e-06,1.0531e-09,-7.28218e-14,35290.9,-12.4583], Tmin=(1027.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH3CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C#C[CH]O[O](26427)',
    structure = SMILES('[CH]=C=CO[O]'),
    E0 = (404.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,492.5,1135,1000,180,180,180,2867.76],'cm^-1')),
        HinderedRotor(inertia=(0.010518,'amu*angstrom^2'), symmetry=1, barrier=(9.12083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13751,0.0383774,-4.88548e-05,3.09451e-08,-7.37348e-12,48769,18.1869], Tmin=(100,'K'), Tmax=(1167.31,'K')), NASAPolynomial(coeffs=[10.2545,0.00577971,-8.1986e-07,1.20449e-12,5.54427e-15,47199.8,-20.8326], Tmin=(1167.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ)"""),
)

species(
    label = 'C[CH]O[O](97)',
    structure = SMILES('C[CH]O[O]'),
    E0 = (162.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.659845,'amu*angstrom^2'), symmetry=1, barrier=(15.1711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000844769,'amu*angstrom^2'), symmetry=1, barrier=(5.4037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88063,0.0212973,-8.38042e-06,-1.15683e-09,1.18565e-12,19629.9,16.0899], Tmin=(100,'K'), Tmax=(1195.85,'K')), NASAPolynomial(coeffs=[6.93453,0.0132692,-5.24916e-06,9.65724e-10,-6.67381e-14,18264.7,-5.84667], Tmin=(1195.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'C2H(33)',
    structure = SMILES('[C]#C'),
    E0 = (557.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(557.301,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""C2H""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C#C[C](C)O[O](26428)',
    structure = SMILES('[CH]=C=C(C)O[O]'),
    E0 = (363.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3120,650,792.5,1650,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.433796,'amu*angstrom^2'), symmetry=1, barrier=(9.97381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67616,0.0550901,-8.07409e-05,6.61287e-08,-2.12061e-11,43753.2,21.8523], Tmin=(100,'K'), Tmax=(891.7,'K')), NASAPolynomial(coeffs=[7.45622,0.0206786,-8.5842e-06,1.51268e-09,-9.88831e-14,43059.6,-3.47947], Tmin=(891.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([CH2])O[O](26429)',
    structure = SMILES('C#CC([CH2])O[O]'),
    E0 = (414.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,492.5,1135,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(2.74685,'amu*angstrom^2'), symmetry=1, barrier=(71.8847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46149,'amu*angstrom^2'), symmetry=1, barrier=(12.3367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.730602,'amu*angstrom^2'), symmetry=1, barrier=(19.3701,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44815,0.0538187,-6.01333e-05,3.32997e-08,-6.73708e-12,49898,23.6377], Tmin=(100,'K'), Tmax=(913.495,'K')), NASAPolynomial(coeffs=[12.329,0.0129826,-4.25889e-06,6.81984e-10,-4.34646e-14,47626,-29.4265], Tmin=(913.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = '[C]#CC(C)O[O](26430)',
    structure = SMILES('[C]#CC(C)O[O]'),
    E0 = (537.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,1380,1390,370,380,2900,435,193.597],'cm^-1')),
        HinderedRotor(inertia=(0.388945,'amu*angstrom^2'), symmetry=1, barrier=(10.2201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387159,'amu*angstrom^2'), symmetry=1, barrier=(10.2226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.49605,'amu*angstrom^2'), symmetry=1, barrier=(66.181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70966,0.0522878,-5.95542e-05,3.13581e-08,-3.78944e-12,64700.1,21.4645], Tmin=(100,'K'), Tmax=(715.35,'K')), NASAPolynomial(coeffs=[9.92681,0.0164602,-5.64853e-06,8.97155e-10,-5.53538e-14,63265.6,-17.2357], Tmin=(715.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(Acetyl)"""),
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
    label = 'C#CCO[O](26153)',
    structure = SMILES('C#CCO[O]'),
    E0 = (242.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,492.5,1135,1000,2175,525],'cm^-1')),
        HinderedRotor(inertia=(0.595333,'amu*angstrom^2'), symmetry=1, barrier=(13.6879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74398,'amu*angstrom^2'), symmetry=1, barrier=(63.0895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29019,0.0394598,-5.33345e-05,4.08817e-08,-1.22873e-11,29219.5,16.1946], Tmin=(100,'K'), Tmax=(942.249,'K')), NASAPolynomial(coeffs=[6.81691,0.0148149,-5.46e-06,8.95135e-10,-5.57483e-14,28607.4,-4.09627], Tmin=(942.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C1OOC1C(26423)',
    structure = SMILES('[CH]=C1OOC1C'),
    E0 = (295.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23895,0.0220754,4.61641e-05,-7.68615e-08,3.06646e-11,35636.2,20.632], Tmin=(100,'K'), Tmax=(982.407,'K')), NASAPolynomial(coeffs=[14.1783,0.0147804,-5.78352e-06,1.20107e-09,-9.48541e-14,31296.5,-46.9039], Tmin=(982.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_P)"""),
)

species(
    label = 'C#C[C](C)OO(26424)',
    structure = SMILES('[CH]=C=C(C)OO'),
    E0 = (211.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,3120,650,792.5,1650,3615,1310,387.5,850,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.677736,'amu*angstrom^2'), symmetry=1, barrier=(15.5825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.675718,'amu*angstrom^2'), symmetry=1, barrier=(15.5361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50148,0.0575454,-7.16632e-05,5.05832e-08,-1.45384e-11,25478.9,21.726], Tmin=(100,'K'), Tmax=(846.138,'K')), NASAPolynomial(coeffs=[8.88938,0.0226184,-9.74288e-06,1.79415e-09,-1.22462e-13,24228.8,-12.6817], Tmin=(846.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([CH2])OO(26425)',
    structure = SMILES('C#CC([CH2])OO'),
    E0 = (262.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00558,0.0595106,-6.271e-05,3.32896e-08,-6.86408e-12,31635.3,24.4669], Tmin=(100,'K'), Tmax=(1194.86,'K')), NASAPolynomial(coeffs=[14.9037,0.0129835,-4.29995e-06,6.99496e-10,-4.51513e-14,28314,-45.0576], Tmin=(1194.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(262.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCOOH)"""),
)

species(
    label = '[C]#CC(C)OO(26426)',
    structure = SMILES('[C]#CC(C)OO'),
    E0 = (385.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,240.467],'cm^-1')),
        HinderedRotor(inertia=(1.79664,'amu*angstrom^2'), symmetry=1, barrier=(75.8469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265996,'amu*angstrom^2'), symmetry=1, barrier=(11.0962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644069,'amu*angstrom^2'), symmetry=1, barrier=(26.8994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83305,'amu*angstrom^2'), symmetry=1, barrier=(75.8465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25981,0.0583907,-6.5421e-05,3.90834e-08,-9.23975e-12,46437.6,22.3011], Tmin=(100,'K'), Tmax=(1036.72,'K')), NASAPolynomial(coeffs=[11.8543,0.0175133,-6.27608e-06,1.04966e-09,-6.80139e-14,44240.9,-29.1932], Tmin=(1036.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = 'HO2(10)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.87083,0.0182042,1.06711e-05,-2.72492e-08,1.19478e-11,33023.8,11.2934], Tmin=(100,'K'), Tmax=(955.249,'K')), NASAPolynomial(coeffs=[8.52653,0.0108962,-3.56564e-06,6.31243e-10,-4.51891e-14,31196.2,-19.6435], Tmin=(955.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CH2CHCCH""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C#CC(C)[O](4289)',
    structure = SMILES('C#CC(C)[O]'),
    E0 = (206.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,248.276],'cm^-1')),
        HinderedRotor(inertia=(0.59556,'amu*angstrom^2'), symmetry=1, barrier=(26.1247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.599472,'amu*angstrom^2'), symmetry=1, barrier=(26.1303,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10534,0.0335727,-4.56046e-06,-2.27819e-08,1.28075e-11,24969.9,17.5705], Tmin=(100,'K'), Tmax=(925.951,'K')), NASAPolynomial(coeffs=[11.9497,0.0118724,-3.14438e-06,4.88913e-10,-3.35953e-14,22254,-33.9865], Tmin=(925.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1=C=COO1(26432)',
    structure = SMILES('CC1=C=COO1'),
    E0 = (393.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2950,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16001,0.0314729,5.40636e-07,-2.23655e-08,1.01967e-11,47406.2,16.2598], Tmin=(100,'K'), Tmax=(1034.52,'K')), NASAPolynomial(coeffs=[11.254,0.0165027,-7.0307e-06,1.38069e-09,-1.01142e-13,44444.1,-33.1444], Tmin=(1034.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene)"""),
)

species(
    label = 'CC1C#COO1(26433)',
    structure = SMILES('CC1C#COO1'),
    E0 = (468.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85658,0.0399862,-1.96103e-05,-6.97782e-10,1.85048e-12,56417.9,16.9705], Tmin=(100,'K'), Tmax=(1307.74,'K')), NASAPolynomial(coeffs=[13.3932,0.0179758,-9.59239e-06,1.9584e-09,-1.41375e-13,52265.3,-46.1227], Tmin=(1307.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + ring(Cyclopentyne)"""),
)

species(
    label = 'C1=COOC=1(26434)',
    structure = SMILES('C1=COOC=1'),
    E0 = (435.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,750.683,751.338,751.366,751.383,751.481,751.496,751.499,751.568,751.587,751.846,751.98],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79508,0.0126821,3.97973e-05,-6.71799e-08,2.81361e-11,52414.5,11.2805], Tmin=(100,'K'), Tmax=(945.135,'K')), NASAPolynomial(coeffs=[13.4398,0.00265609,1.22671e-07,1.42216e-11,-8.78205e-15,48838.1,-47.7493], Tmin=(945.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(435.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene)"""),
)

species(
    label = 'C[C]1C=COO1(26435)',
    structure = SMILES('CC1=C[CH]OO1'),
    E0 = (81.1547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49683,0.0168717,5.50241e-05,-8.30537e-08,3.24042e-11,9829.3,19.5015], Tmin=(100,'K'), Tmax=(975.788,'K')), NASAPolynomial(coeffs=[12.6044,0.0159904,-5.95881e-06,1.19982e-09,-9.33835e-14,5926.12,-38.9064], Tmin=(975.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.1547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxolene) + radical(C=CCJO)"""),
)

species(
    label = 'CC1C=[C]OO1(26436)',
    structure = SMILES('CC1C=[C]OO1'),
    E0 = (203.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38922,0.0221283,3.53449e-05,-6.09514e-08,2.43478e-11,24494.6,22.1464], Tmin=(100,'K'), Tmax=(983.272,'K')), NASAPolynomial(coeffs=[11.8024,0.0166442,-6.33984e-06,1.24614e-09,-9.41975e-14,21057.4,-31.173], Tmin=(983.272,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1C=COO1(26437)',
    structure = SMILES('[CH2]C1C=COO1'),
    E0 = (177.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10907,0.023165,4.87186e-05,-8.65618e-08,3.63638e-11,21408.8,22.0175], Tmin=(100,'K'), Tmax=(950.281,'K')), NASAPolynomial(coeffs=[16.5314,0.00887181,-1.98459e-06,4.07473e-10,-3.80181e-14,16572.1,-57.8533], Tmin=(950.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(CJCOOH)"""),
)

species(
    label = '[C]1[CH]OOC=1(26438)',
    structure = SMILES('[C]1[CH]OOC=1'),
    E0 = (360.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85,819.85],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08524,0.00281293,6.91832e-05,-9.84561e-08,3.94476e-11,43441.6,15.165], Tmin=(100,'K'), Tmax=(940.003,'K')), NASAPolynomial(coeffs=[13.8419,0.00124637,1.14098e-06,-1.69689e-10,1.91103e-15,39466.3,-46.4528], Tmin=(940.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = 'C[C]1[C]=COO1(26439)',
    structure = SMILES('CC1=[C][CH]OO1'),
    E0 = (318.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46446,0.0214584,3.03075e-05,-5.39146e-08,2.1511e-11,38432.6,20.0918], Tmin=(100,'K'), Tmax=(992.782,'K')), NASAPolynomial(coeffs=[11.46,0.0154162,-6.19456e-06,1.23909e-09,-9.39135e-14,35158.2,-30.737], Tmin=(992.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxolene) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C1[C]=COO1(26440)',
    structure = SMILES('[CH2]C1[C]=COO1'),
    E0 = (415.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07898,0.0277328,2.40244e-05,-5.73753e-08,2.54107e-11,50012,22.5992], Tmin=(100,'K'), Tmax=(956.119,'K')), NASAPolynomial(coeffs=[15.3356,0.00838229,-2.2681e-06,4.57837e-10,-3.94568e-14,45826.6,-49.393], Tmin=(956.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S) + radical(CJCOOH)"""),
)

species(
    label = 'CC1[C]=[C]OO1(26441)',
    structure = SMILES('CC1[C]=[C]OO1'),
    E0 = (440.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35324,0.0267531,1.05207e-05,-3.17163e-08,1.34358e-11,53098.1,22.7499], Tmin=(100,'K'), Tmax=(1011.54,'K')), NASAPolynomial(coeffs=[10.7011,0.0159989,-6.53555e-06,1.27611e-09,-9.39657e-14,50270.6,-23.2477], Tmin=(1011.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'CC=C=CO[O](26151)',
    structure = SMILES('CC=C=CO[O]'),
    E0 = (214.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.690031,'amu*angstrom^2'), symmetry=1, barrier=(15.8652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.691142,'amu*angstrom^2'), symmetry=1, barrier=(15.8907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75181,0.0496222,-4.93503e-05,2.71166e-08,-6.10875e-12,25868.3,20.7585], Tmin=(100,'K'), Tmax=(1064.11,'K')), NASAPolynomial(coeffs=[9.44376,0.0207069,-8.58861e-06,1.57823e-09,-1.08556e-13,24231.4,-16.8283], Tmin=(1064.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = '[C]1=COOC1(26442)',
    structure = SMILES('[C]1=COOC1'),
    E0 = (243.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,180,180,568.465,1028.48,1028.48,1028.48,1028.48,1028.49,1028.49,1028.5,1028.51,2930.72],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93792,0.0131609,3.16359e-05,-5.09516e-08,2.04045e-11,29332.9,15.7895], Tmin=(100,'K'), Tmax=(965.366,'K')), NASAPolynomial(coeffs=[9.78787,0.0102781,-3.50705e-06,6.80132e-10,-5.25119e-14,26822.1,-23.1703], Tmin=(965.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(Cds_S)"""),
)

species(
    label = 'CC1[CH]OOC=1(26443)',
    structure = SMILES('CC1[CH]OOC=1'),
    E0 = (83.8912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28116,0.0181387,6.36799e-05,-1.00975e-07,4.11328e-11,10169.6,18.3341], Tmin=(100,'K'), Tmax=(950.27,'K')), NASAPolynomial(coeffs=[15.9308,0.0103825,-2.52771e-06,5.10912e-10,-4.58198e-14,5331.46,-58.6283], Tmin=(950.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.8912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(12dioxolene) + radical(C=CCJO)"""),
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
    E0 = (299.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (541.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (720.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (575.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (625.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (749.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (662.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (322.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (386.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (343.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (429.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (258.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (301.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (449.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (617.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (695.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (605.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (361.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (435.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (353.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (497.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (530.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (626.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (652.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (284.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (663.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (345.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', 'CH3CHCCH(24127)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.19343e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDeC;Y_rad] for rate rule [C_rad/H/CtCs;O2_birad]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', 'C#C[CH]O[O](26427)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.64257e+09,'m^3/(mol*s)'), n=-0.82, Ea=(0.004184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_sec_rad;C_methyl] for rate rule [C_rad/H/OneDeO;C_methyl]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]O[O](97)', 'C2H(33)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.18e+12,'cm^3/(mol*s)'), n=-0.085, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_rad/H/CsO] for rate rule [Ct_rad/Ct;C_rad/H/CsO]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C#C[C](C)O[O](26428)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [C_rad/OneDe;H_rad] for rate rule [C_rad/OneDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C#CC([CH2])O[O](26429)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[C]#CC(C)O[O](26430)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2(S)(23)', 'C#CCO[O](26153)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/OneDeO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC(C)O[O](26146)'],
    products = ['[CH]=C1OOC1C(26423)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#C[C](C)OO(26424)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.47099e+07,'s^-1'), n=1.47622, Ea=(175.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_OneDe/Cs;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC(C)O[O](26146)'],
    products = ['C#CC([CH2])OO(26425)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.15, Ea=(143.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 246 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]#CC(C)OO(26426)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC(C)O[O](26146)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC(C)O[O](26146)'],
    products = ['HO2(10)', 'CH2CHCCH(26391)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.4759e+06,'s^-1','*|/',5), n=1.829, Ea=(101.449,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R2OO_2H_HDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O(4)', 'C#CC(C)[O](4289)'],
    products = ['C#CC(C)O[O](26146)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'CC1=C=COO1(26432)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'CC1C#COO1(26433)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH3(17)', 'C1=COOC=1(26434)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.021171,'m^3/(mol*s)'), n=2.41, Ea=(33.7021,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;CsJ-HHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['C[C]1C=COO1(26435)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.49816e+10,'s^-1'), n=0.9445, Ea=(160.64,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out] for rate rule [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out_OOH/Cs]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['CC1C=[C]OO1(26436)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['[CH2]C1C=COO1(26437)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R3H_SS_12cy5;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', '[C]1[CH]OOC=1(26438)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.68936e+07,'m^3/(mol*s)'), n=0.0549843, Ea=(0.0928142,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;C_methyl]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'C[C]1[C]=COO1(26439)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C1[C]=COO1(26440)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.48677e-12,'cm^3/(molecule*s)'), n=0.6, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for C_rad/H2/Cs;H_rad
Exact match found for rate rule [C_rad/H2/Cs;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'CC1[C]=[C]OO1(26441)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34601e+06,'m^3/(mol*s)'), n=0.278532, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC=C=CO[O](26151)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(S)(23)', '[C]1=COOC1(26442)'],
    products = ['CC1[C]=COO1(26431)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC1[C]=COO1(26431)'],
    products = ['CC1[CH]OOC=1(26443)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HR!H)CJ;CJ;CH3] + [cCs(-HR!H)CJ;CdsJ;C] for rate rule [cCs(-HR!H)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

network(
    label = '4416',
    isomers = [
        'C#CC(C)O[O](26146)',
        'CC1[C]=COO1(26431)',
    ],
    reactants = [
        ('O2(2)', 'CH3CHCCH(24127)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4416',
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

