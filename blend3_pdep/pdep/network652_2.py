species(
    label = 'S(247)(246)',
    structure = SMILES('CC(=O)C=CC=O'),
    E0 = (-249.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,338.299],'cm^-1')),
        HinderedRotor(inertia=(0.102894,'amu*angstrom^2'), symmetry=1, barrier=(8.35729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1029,'amu*angstrom^2'), symmetry=1, barrier=(8.35732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102901,'amu*angstrom^2'), symmetry=1, barrier=(8.35737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3947.01,'J/mol'), sigma=(6.07703,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=616.51 K, Pc=39.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72359,0.0430426,-2.26393e-05,4.45245e-09,-2.82809e-13,-30034.5,19.0288], Tmin=(100,'K'), Tmax=(2443.83,'K')), NASAPolynomial(coeffs=[29.6812,0.00674437,-5.16283e-06,9.95166e-10,-6.31701e-14,-45547.2,-139.895], Tmin=(2443.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'S(252)(251)',
    structure = SMILES('CC(O)=C=CC=O'),
    E0 = (-191.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.833403,'amu*angstrom^2'), symmetry=1, barrier=(19.1616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.82892,'amu*angstrom^2'), symmetry=1, barrier=(19.0585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828534,'amu*angstrom^2'), symmetry=1, barrier=(19.0496,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4176.29,'J/mol'), sigma=(6.51688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=652.33 K, Pc=34.24 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986472,0.0659047,-6.4425e-05,3.18942e-08,-6.35286e-12,-22904.6,21.9789], Tmin=(100,'K'), Tmax=(1200.1,'K')), NASAPolynomial(coeffs=[13.8493,0.0230323,-1.08392e-05,2.12683e-09,-1.51862e-13,-25992,-42.4231], Tmin=(1200.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
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
    label = 'S(780)(779)',
    structure = SMILES('O=[C]C=CC=O'),
    E0 = (-44.2053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.882295,'amu*angstrom^2'), symmetry=1, barrier=(20.2857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879664,'amu*angstrom^2'), symmetry=1, barrier=(20.2252,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3741.66,'J/mol'), sigma=(5.7541,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.44 K, Pc=44.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70479,0.0574491,-9.74937e-05,9.58172e-08,-3.8062e-11,-5240.81,19.2682], Tmin=(100,'K'), Tmax=(736.114,'K')), NASAPolynomial(coeffs=[5.25854,0.0278824,-1.6346e-05,3.39827e-09,-2.46552e-13,-5486.14,5.0995], Tmin=(736.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.2053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
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
    label = 'C[C]=C=CC=O(7615)',
    structure = SMILES('CC#CC=C[O]'),
    E0 = (161.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.56284,'amu*angstrom^2'), symmetry=1, barrier=(35.9327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57121,'amu*angstrom^2'), symmetry=1, barrier=(36.1253,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80258,0.0332802,2.22527e-05,-6.26465e-08,2.95311e-11,19480.5,19.4671], Tmin=(100,'K'), Tmax=(920.512,'K')), NASAPolynomial(coeffs=[16.5497,0.00825504,-6.12625e-07,6.96751e-12,-3.27546e-15,15110.8,-59.4458], Tmin=(920.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=COJ)"""),
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
    label = 'CC(=O)[C]=CC=O(7584)',
    structure = SMILES('CC([O])=C=CC=O'),
    E0 = (-53.5394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.653958,'amu*angstrom^2'), symmetry=1, barrier=(15.0358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652529,'amu*angstrom^2'), symmetry=1, barrier=(15.0029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48531,0.0601767,-6.50487e-05,3.92374e-08,-9.97775e-12,-6352.84,21.7039], Tmin=(100,'K'), Tmax=(931.209,'K')), NASAPolynomial(coeffs=[8.83561,0.028604,-1.41916e-05,2.82842e-09,-2.03226e-13,-7721.79,-13.2332], Tmin=(931.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.5394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=CC=C=[C]O(7616)',
    structure = SMILES('[O]C=CC#CO'),
    E0 = (61.9116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44972,'amu*angstrom^2'), symmetry=1, barrier=(33.3319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45335,'amu*angstrom^2'), symmetry=1, barrier=(33.4153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8124,0.0347446,3.22938e-06,-3.94368e-08,2.01747e-11,7537.22,19.3479], Tmin=(100,'K'), Tmax=(948.051,'K')), NASAPolynomial(coeffs=[17.0656,0.00367907,-2.9063e-07,7.7097e-11,-1.21406e-14,3149,-61.3163], Tmin=(948.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)=C=CC=O(7617)',
    structure = SMILES('C=C(O)[C]=CC=O'),
    E0 = (-44.4643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,350,440,435,1725,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.08303,'amu*angstrom^2'), symmetry=1, barrier=(24.9011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08446,'amu*angstrom^2'), symmetry=1, barrier=(24.9338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08256,'amu*angstrom^2'), symmetry=1, barrier=(24.8902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.661659,0.0720618,-8.35572e-05,4.86205e-08,-1.10886e-11,-5226.36,21.8972], Tmin=(100,'K'), Tmax=(1074.06,'K')), NASAPolynomial(coeffs=[15.2881,0.0175901,-7.48361e-06,1.40174e-09,-9.79087e-14,-8368.29,-49.7119], Tmin=(1074.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.4643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C)"""),
)

species(
    label = 'HCO(14)(15)',
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
    label = '[CH]=C=C(C)O(7618)',
    structure = SMILES('[CH]=C=C(C)O'),
    E0 = (86.6206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.984925,'amu*angstrom^2'), symmetry=1, barrier=(22.6454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10635,0.0503488,-5.11114e-05,2.5565e-08,-4.75667e-12,10533.6,19.0983], Tmin=(100,'K'), Tmax=(1545.92,'K')), NASAPolynomial(coeffs=[14.256,0.00720865,-4.07516e-07,-1.15074e-10,1.30724e-14,7557.26,-46.5464], Tmin=(1545.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.6206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ)"""),
)

species(
    label = 'CC(O)=C=[C]C=O(7619)',
    structure = SMILES('CC(O)=C=C=C[O]'),
    E0 = (5.39251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.19634,'amu*angstrom^2'), symmetry=1, barrier=(27.5061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19457,'amu*angstrom^2'), symmetry=1, barrier=(27.4654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.341136,0.0756378,-8.14065e-05,4.07995e-08,-7.55309e-12,821.889,25.803], Tmin=(100,'K'), Tmax=(1532.39,'K')), NASAPolynomial(coeffs=[22.2459,0.00415139,8.31734e-07,-3.13139e-10,2.45532e-14,-4629.67,-88.0078], Tmin=(1532.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.39251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CC(O)=C=C[C]=O(7620)',
    structure = SMILES('CC(O)=C=C[C]=O'),
    E0 = (-36.6733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1855,455,950,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.849071,'amu*angstrom^2'), symmetry=1, barrier=(19.5218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.848978,'amu*angstrom^2'), symmetry=1, barrier=(19.5197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.848812,'amu*angstrom^2'), symmetry=1, barrier=(19.5159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11327,0.0653135,-7.04144e-05,3.86206e-08,-8.51992e-12,-4308.29,22.2071], Tmin=(100,'K'), Tmax=(1089.85,'K')), NASAPolynomial(coeffs=[12.9142,0.0220018,-1.08031e-05,2.15622e-09,-1.55423e-13,-6880.54,-35.7407], Tmin=(1089.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-36.6733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = 'CC([O])[C]=CC=O(7588)',
    structure = SMILES('CC([O])[C]=CC=O'),
    E0 = (155.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,340.004,340.28],'cm^-1')),
        HinderedRotor(inertia=(0.146795,'amu*angstrom^2'), symmetry=1, barrier=(12.0322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146911,'amu*angstrom^2'), symmetry=1, barrier=(12.0323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146902,'amu*angstrom^2'), symmetry=1, barrier=(12.0314,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69495,0.0553565,-4.39351e-05,1.76841e-08,-3.01337e-12,18767.3,24.9739], Tmin=(100,'K'), Tmax=(1312.37,'K')), NASAPolynomial(coeffs=[10.0004,0.0300422,-1.50016e-05,2.98625e-09,-2.13492e-13,16587.4,-17.3527], Tmin=(1312.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C(O)[C]=CC=O(7621)',
    structure = SMILES('[CH2]C(O)[C]=CC=O'),
    E0 = (136.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,351.176],'cm^-1')),
        HinderedRotor(inertia=(0.125409,'amu*angstrom^2'), symmetry=1, barrier=(10.9831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125442,'amu*angstrom^2'), symmetry=1, barrier=(10.9835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125499,'amu*angstrom^2'), symmetry=1, barrier=(10.9837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125448,'amu*angstrom^2'), symmetry=1, barrier=(10.9831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28861,0.0653089,-7.11008e-05,4.26404e-08,-1.07621e-11,16523.8,27.2814], Tmin=(100,'K'), Tmax=(937.534,'K')), NASAPolynomial(coeffs=[9.42849,0.0305799,-1.55362e-05,3.12905e-09,-2.26125e-13,14997.5,-11.4638], Tmin=(937.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(CJCO)"""),
)

species(
    label = 'CC(O)=C=[C]C[O](7622)',
    structure = SMILES('C[C](O)C#CC[O]'),
    E0 = (183.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2250,500,550,2750,2850,1437.5,1250,1305,750,350,186.269,186.368,1386.05],'cm^-1')),
        HinderedRotor(inertia=(0.247197,'amu*angstrom^2'), symmetry=1, barrier=(6.09796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247449,'amu*angstrom^2'), symmetry=1, barrier=(6.09564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00447117,'amu*angstrom^2'), symmetry=1, barrier=(6.09566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.5961,'amu*angstrom^2'), symmetry=1, barrier=(63.7238,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20723,0.0659955,-8.88176e-05,7.17049e-08,-2.30887e-11,22182.7,28.4473], Tmin=(100,'K'), Tmax=(893.986,'K')), NASAPolynomial(coeffs=[6.92393,0.0304159,-1.23387e-05,2.15905e-09,-1.41061e-13,21560.2,3.74343], Tmin=(893.986,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCOJ) + radical(C2CsJOH)"""),
)

species(
    label = 'C[C](O)C=[C]C=O(7595)',
    structure = SMILES('CC(O)=C[C]=C[O]'),
    E0 = (-24.9975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17014,'amu*angstrom^2'), symmetry=1, barrier=(26.9037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16965,'amu*angstrom^2'), symmetry=1, barrier=(26.8925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16976,'amu*angstrom^2'), symmetry=1, barrier=(26.8951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.928288,0.0811938,-8.62323e-05,4.28041e-08,-7.74497e-12,-2805.62,28.4291], Tmin=(100,'K'), Tmax=(1623.21,'K')), NASAPolynomial(coeffs=[22.5678,0.00443465,2.12776e-06,-6.43719e-10,4.90632e-14,-7948.99,-88.6548], Tmin=(1623.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.9975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CC(O)=[C]C[C]=O(7623)',
    structure = SMILES('CC(O)=[C]C[C]=O'),
    E0 = (54.0795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,291.665],'cm^-1')),
        HinderedRotor(inertia=(0.221726,'amu*angstrom^2'), symmetry=1, barrier=(13.3732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221637,'amu*angstrom^2'), symmetry=1, barrier=(13.3728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221782,'amu*angstrom^2'), symmetry=1, barrier=(13.3737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222937,'amu*angstrom^2'), symmetry=1, barrier=(13.3742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808721,0.0631018,-5.69119e-05,2.49064e-08,-4.28248e-12,6624.97,26.2834], Tmin=(100,'K'), Tmax=(1402.58,'K')), NASAPolynomial(coeffs=[16.9855,0.0169681,-7.57439e-06,1.4559e-09,-1.02641e-13,2087.07,-57.2331], Tmin=(1402.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.0795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([O])=CC=C[O](7453)',
    structure = SMILES('CC([O])=CC=C[O]'),
    E0 = (-86.1882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.912938,'amu*angstrom^2'), symmetry=1, barrier=(20.9902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.915391,'amu*angstrom^2'), symmetry=1, barrier=(21.0466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.756235,0.0562742,-1.51561e-05,-3.6247e-08,2.26938e-11,-10235.2,24.2724], Tmin=(100,'K'), Tmax=(917.521,'K')), NASAPolynomial(coeffs=[21.0957,0.00729401,3.01004e-08,-1.33448e-10,7.27657e-15,-15638.3,-81.2071], Tmin=(917.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-86.1882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(O)=CC=C[O](5193)',
    structure = SMILES('C=C(O)[CH]C=C[O]'),
    E0 = (-76.2221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35339,'amu*angstrom^2'), symmetry=1, barrier=(31.117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35269,'amu*angstrom^2'), symmetry=1, barrier=(31.1009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35312,'amu*angstrom^2'), symmetry=1, barrier=(31.1109,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400843,0.0594199,-7.08619e-06,-5.37493e-08,3.0632e-11,-9019.51,24.2686], Tmin=(100,'K'), Tmax=(924.396,'K')), NASAPolynomial(coeffs=[25.3551,0.00260988,2.06445e-06,-4.65225e-10,2.62503e-14,-15819.3,-105.985], Tmin=(924.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.2221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C[C](O)C=C[C]=O(7597)',
    structure = SMILES('CC(O)=C[CH][C]=O'),
    E0 = (-21.6336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1855,455,950,350,440,435,1725,375.438],'cm^-1')),
        HinderedRotor(inertia=(0.191574,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191577,'amu*angstrom^2'), symmetry=1, barrier=(19.162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191572,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191577,'amu*angstrom^2'), symmetry=1, barrier=(19.1619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.943088,0.0514379,-8.03709e-06,-3.43369e-08,1.82981e-11,-2477.79,25.7054], Tmin=(100,'K'), Tmax=(990.28,'K')), NASAPolynomial(coeffs=[19.613,0.0127286,-4.99884e-06,1.04521e-09,-8.29423e-14,-7975.14,-73.2705], Tmin=(990.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.6336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([O])=[C]CC=O(7624)',
    structure = SMILES('CC([O])=[C]CC=O'),
    E0 = (31.9237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,353.881,354.148],'cm^-1')),
        HinderedRotor(inertia=(0.00134351,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122626,'amu*angstrom^2'), symmetry=1, barrier=(10.8553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122651,'amu*angstrom^2'), symmetry=1, barrier=(10.8509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54122,0.0527991,-3.7905e-05,1.30063e-08,-1.80161e-12,3928.46,24.3284], Tmin=(100,'K'), Tmax=(1643.99,'K')), NASAPolynomial(coeffs=[13.6165,0.0234186,-1.10978e-05,2.13555e-09,-1.48499e-13,-41.8842,-39.931], Tmin=(1643.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=[C]CC=O(7625)',
    structure = SMILES('[CH2]C(O)=[C]CC=O'),
    E0 = (53.0355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,350,440,435,1725,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,277.581],'cm^-1')),
        HinderedRotor(inertia=(0.330374,'amu*angstrom^2'), symmetry=1, barrier=(18.0947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00219251,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.329484,'amu*angstrom^2'), symmetry=1, barrier=(18.0948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.330383,'amu*angstrom^2'), symmetry=1, barrier=(18.0955,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755888,0.0617819,-4.9183e-05,1.55815e-08,-8.70302e-13,6503.59,26.3414], Tmin=(100,'K'), Tmax=(1120.64,'K')), NASAPolynomial(coeffs=[17.0284,0.0171408,-7.42234e-06,1.44191e-09,-1.03839e-13,2012.42,-57.7836], Tmin=(1120.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.0355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(C=C(O)CJ)"""),
)

species(
    label = 'CC(O)[C]=C[C]=O(7626)',
    structure = SMILES('CC(O)[C]=C[C]=O'),
    E0 = (85.6467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.589682,'amu*angstrom^2'), symmetry=1, barrier=(13.5579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589875,'amu*angstrom^2'), symmetry=1, barrier=(13.5624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589521,'amu*angstrom^2'), symmetry=1, barrier=(13.5543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.589405,'amu*angstrom^2'), symmetry=1, barrier=(13.5516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.938233,0.0734944,-9.53377e-05,6.89363e-08,-2.05817e-11,10405.7,27.4286], Tmin=(100,'K'), Tmax=(808.472,'K')), NASAPolynomial(coeffs=[9.63505,0.0304642,-1.54984e-05,3.09815e-09,-2.22077e-13,8999.52,-12.6791], Tmin=(808.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC(=O)[C]=CC[O](7593)',
    structure = SMILES('CC([O])=C=CC[O]'),
    E0 = (106.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180,3259.98],'cm^-1')),
        HinderedRotor(inertia=(0.418015,'amu*angstrom^2'), symmetry=1, barrier=(9.61098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.418255,'amu*angstrom^2'), symmetry=1, barrier=(9.6165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29486,0.0649342,-8.72606e-05,7.4083e-08,-2.59591e-11,12900.2,25.7155], Tmin=(100,'K'), Tmax=(808.554,'K')), NASAPolynomial(coeffs=[5.85469,0.0344025,-1.58267e-05,2.98776e-09,-2.05742e-13,12423.5,6.29782], Tmin=(808.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=C=CC[O](7627)',
    structure = SMILES('C=C(O)[C]=CC[O]'),
    E0 = (115.568,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.565207,'amu*angstrom^2'), symmetry=1, barrier=(12.9952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564091,'amu*angstrom^2'), symmetry=1, barrier=(12.9696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.95879,'amu*angstrom^2'), symmetry=1, barrier=(91.0204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75812,0.0731697,-9.16994e-05,6.31982e-08,-1.74356e-11,14014.6,24.8969], Tmin=(100,'K'), Tmax=(886.799,'K')), NASAPolynomial(coeffs=[11.4675,0.0248632,-9.98885e-06,1.7701e-09,-1.1798e-13,12115.2,-25.4831], Tmin=(886.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.568,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CC([O])=[C]C=CO(7454)',
    structure = SMILES('CC([O])=[C]C=CO'),
    E0 = (-28.6553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03859,'amu*angstrom^2'), symmetry=1, barrier=(23.8792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04236,'amu*angstrom^2'), symmetry=1, barrier=(23.9658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04283,'amu*angstrom^2'), symmetry=1, barrier=(23.9768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.794495,0.082309,-9.15949e-05,4.78442e-08,-9.10871e-12,-3253.9,28.024], Tmin=(100,'K'), Tmax=(1537.86,'K')), NASAPolynomial(coeffs=[22.2025,0.0040977,2.63382e-06,-7.82469e-10,6.06784e-14,-8151.85,-85.7487], Tmin=(1537.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.6553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(O)=C=C[CH]O(7628)',
    structure = SMILES('[CH2]C(O)=[C]C=CO'),
    E0 = (-7.54347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(1.15298,'amu*angstrom^2'), symmetry=1, barrier=(26.5093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15382,'amu*angstrom^2'), symmetry=1, barrier=(26.5285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15375,'amu*angstrom^2'), symmetry=1, barrier=(26.5269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15289,'amu*angstrom^2'), symmetry=1, barrier=(26.5073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71522,0.0928755,-0.000108217,5.68818e-08,-1.06769e-11,-673.167,30.5215], Tmin=(100,'K'), Tmax=(1599.7,'K')), NASAPolynomial(coeffs=[25.5084,-0.00234577,6.52757e-06,-1.54657e-09,1.12258e-13,-5909.3,-102.749], Tmin=(1599.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.54347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(O)CJ) + radical(C=CJC=C)"""),
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
    label = 'O=CC=C=CO(7629)',
    structure = SMILES('O=CC=C=CO'),
    E0 = (-149.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(1.15186,'amu*angstrom^2'), symmetry=1, barrier=(26.4835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1532,'amu*angstrom^2'), symmetry=1, barrier=(26.5143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25519,0.05143,-4.02864e-05,6.80829e-09,3.07965e-12,-17880.3,19.0076], Tmin=(100,'K'), Tmax=(998.283,'K')), NASAPolynomial(coeffs=[16.745,0.00799049,-3.00262e-06,6.00153e-10,-4.62935e-14,-21901,-60.3438], Tmin=(998.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-149.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC(O)=C[C]C=O(7632)',
    structure = SMILES('CC(O)=C[C]C=O'),
    E0 = (79.4141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.778439,'amu*angstrom^2'), symmetry=1, barrier=(17.8978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778892,'amu*angstrom^2'), symmetry=1, barrier=(17.9083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.56374,'amu*angstrom^2'), symmetry=1, barrier=(58.9455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57311,'amu*angstrom^2'), symmetry=1, barrier=(36.169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02913,0.0671072,-7.02132e-05,4.01889e-08,-9.44086e-12,9656.86,34.1533], Tmin=(100,'K'), Tmax=(1020.45,'K')), NASAPolynomial(coeffs=[11.0908,0.0276666,-1.22373e-05,2.31247e-09,-1.61439e-13,7603.39,-14.5922], Tmin=(1020.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.4141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
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
    label = 'C=C=C(C)O(7630)',
    structure = SMILES('C=C=C(C)O'),
    E0 = (-67.8562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(0.902443,'amu*angstrom^2'), symmetry=1, barrier=(20.7489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900607,'amu*angstrom^2'), symmetry=1, barrier=(20.7067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69403,0.0423221,-1.88022e-05,-1.16794e-08,9.29029e-12,-8070.53,16.3593], Tmin=(100,'K'), Tmax=(945.589,'K')), NASAPolynomial(coeffs=[13.6968,0.0121323,-3.56472e-06,5.98969e-10,-4.24187e-14,-11260.7,-45.7418], Tmin=(945.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.8562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC(O)=C=C=CO(7631)',
    structure = SMILES('CC(O)=C=C=CO'),
    E0 = (-136.07,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725],'cm^-1')),
        HinderedRotor(inertia=(1.10911,'amu*angstrom^2'), symmetry=1, barrier=(25.5006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1103,'amu*angstrom^2'), symmetry=1, barrier=(25.528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1107,'amu*angstrom^2'), symmetry=1, barrier=(25.5371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00863,0.0856659,-9.57221e-05,4.91098e-08,-9.18995e-12,-16164.1,26.7861], Tmin=(100,'K'), Tmax=(1548.43,'K')), NASAPolynomial(coeffs=[24.6736,0.00158147,2.91807e-06,-7.58233e-10,5.60698e-14,-21990.7,-101.478], Tmin=(1548.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C(=O)C=CC=O(2803)',
    structure = SMILES('C=C([O])C=CC=O'),
    E0 = (-105.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.804583,'amu*angstrom^2'), symmetry=1, barrier=(18.4989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.805883,'amu*angstrom^2'), symmetry=1, barrier=(18.5288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16866,0.0608176,-5.94524e-05,2.9485e-08,-5.85639e-12,-12604.2,21.9805], Tmin=(100,'K'), Tmax=(1208.63,'K')), NASAPolynomial(coeffs=[13.3958,0.0203512,-9.23039e-06,1.78305e-09,-1.26336e-13,-15559.8,-39.3255], Tmin=(1208.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ)"""),
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
    label = 'CHCHCHO(2768)',
    structure = SMILES('[CH]=CC=O'),
    E0 = (171.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,688.575,691.322,692.403,693.509,695.331],'cm^-1')),
        HinderedRotor(inertia=(0.0100227,'amu*angstrom^2'), symmetry=1, barrier=(3.38495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95468,0.0145151,1.89382e-05,-3.59857e-08,1.47962e-11,20706.8,12.171], Tmin=(100,'K'), Tmax=(981.049,'K')), NASAPolynomial(coeffs=[9.63601,0.00810119,-3.09994e-06,6.30252e-10,-4.9152e-14,18393.6,-25.043], Tmin=(981.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CHCHCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=CC(C)=O(7465)',
    structure = SMILES('[CH]=CC(C)=O'),
    E0 = (120.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710],'cm^-1')),
        HinderedRotor(inertia=(0.223718,'amu*angstrom^2'), symmetry=1, barrier=(5.14373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225086,'amu*angstrom^2'), symmetry=1, barrier=(5.17517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72572,0.0297242,-1.46437e-05,2.73924e-09,-1.34393e-13,14549.4,16.4009], Tmin=(100,'K'), Tmax=(2060.72,'K')), NASAPolynomial(coeffs=[13.6528,0.0132567,-6.10919e-06,1.09503e-09,-7.04091e-14,9038.89,-46.66], Tmin=(2060.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'CC(=O)C=[C]C=O(7585)',
    structure = SMILES('CC(=O)C=C=C[O]'),
    E0 = (-53.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.863481,'amu*angstrom^2'), symmetry=1, barrier=(19.8531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.863734,'amu*angstrom^2'), symmetry=1, barrier=(19.8589,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18871,0.0546892,-4.44353e-05,1.77663e-08,-2.81993e-12,-6296.83,23.6411], Tmin=(100,'K'), Tmax=(1501.54,'K')), NASAPolynomial(coeffs=[15.0289,0.0178197,-7.60341e-06,1.41325e-09,-9.72104e-14,-10453.1,-48.7556], Tmin=(1501.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-53.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'CC(=O)C=C[C]=O(7586)',
    structure = SMILES('CC(=O)[CH]C=C=O'),
    E0 = (-90.2611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,292.839,2213.46],'cm^-1')),
        HinderedRotor(inertia=(0.00199692,'amu*angstrom^2'), symmetry=1, barrier=(0.121771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172437,'amu*angstrom^2'), symmetry=1, barrier=(10.3441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.532925,'amu*angstrom^2'), symmetry=1, barrier=(32.1618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91384,0.05033,-4.44228e-05,2.31827e-08,-5.42252e-12,-10784.5,23.5552], Tmin=(100,'K'), Tmax=(971.563,'K')), NASAPolynomial(coeffs=[6.55997,0.0312011,-1.4889e-05,2.91675e-09,-2.07633e-13,-11687.3,1.27448], Tmin=(971.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.2611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-(Cdd-O2d)CsH) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[CH2]C([O])C=CC=O(7587)',
    structure = SMILES('[CH2]C([O])C=CC=O'),
    E0 = (129.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,390.996,404.862],'cm^-1')),
        HinderedRotor(inertia=(0.117982,'amu*angstrom^2'), symmetry=1, barrier=(13.4442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120356,'amu*angstrom^2'), symmetry=1, barrier=(13.3783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116909,'amu*angstrom^2'), symmetry=1, barrier=(13.4118,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33486,0.0587495,-4.85124e-05,1.97568e-08,-3.27817e-12,15626.6,26.428], Tmin=(100,'K'), Tmax=(1396.45,'K')), NASAPolynomial(coeffs=[13.1521,0.0248999,-1.21527e-05,2.39859e-09,-1.7061e-13,12326.1,-34.5297], Tmin=(1396.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = 'CC(=O)C=[C]C[O](7589)',
    structure = SMILES('CC(=O)C=[C]C[O]'),
    E0 = (147.892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,180,2170.58,2171.12],'cm^-1')),
        HinderedRotor(inertia=(0.238462,'amu*angstrom^2'), symmetry=1, barrier=(5.48271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238624,'amu*angstrom^2'), symmetry=1, barrier=(5.48643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238709,'amu*angstrom^2'), symmetry=1, barrier=(5.4884,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7306,0.0609164,-9.78591e-05,1.04692e-07,-4.30921e-11,17858.2,26.4444], Tmin=(100,'K'), Tmax=(816.624,'K')), NASAPolynomial(coeffs=[-0.866294,0.0471679,-2.39873e-05,4.69441e-09,-3.28286e-13,19164.9,43.8506], Tmin=(816.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = 'CC(=O)[CH]C[C]=O(7590)',
    structure = SMILES('CC([O])=CC[C]=O'),
    E0 = (-45.9575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,408.961,410.831],'cm^-1')),
        HinderedRotor(inertia=(0.0772605,'amu*angstrom^2'), symmetry=1, barrier=(9.48362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0788662,'amu*angstrom^2'), symmetry=1, barrier=(9.49735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0797959,'amu*angstrom^2'), symmetry=1, barrier=(9.51697,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09955,0.0556405,-4.28036e-05,1.58893e-08,-2.33272e-12,-5416.27,26.2785], Tmin=(100,'K'), Tmax=(1614.97,'K')), NASAPolynomial(coeffs=[16.2791,0.0180429,-7.88235e-06,1.47346e-09,-1.01106e-13,-10319.1,-54.2298], Tmin=(1614.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.9575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=O)C[CH]C=O(7591)',
    structure = SMILES('C=C([O])CC=C[O]'),
    E0 = (-55.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,300.616,300.618,300.623,300.634],'cm^-1')),
        HinderedRotor(inertia=(0.2712,'amu*angstrom^2'), symmetry=1, barrier=(17.3994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2713,'amu*angstrom^2'), symmetry=1, barrier=(17.4,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.978133,0.0515904,-7.28955e-06,-3.97615e-08,2.25807e-11,-6532.62,26.4979], Tmin=(100,'K'), Tmax=(929.656,'K')), NASAPolynomial(coeffs=[19.5071,0.00995383,-1.56325e-06,2.01749e-10,-1.71757e-14,-11623.6,-70.3944], Tmin=(929.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC([O])C=[C]C=O(7592)',
    structure = SMILES('CC([O])C=C=C[O]'),
    E0 = (114.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03546,'amu*angstrom^2'), symmetry=1, barrier=(23.8072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03452,'amu*angstrom^2'), symmetry=1, barrier=(23.7857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.867567,0.05506,-1.75661e-05,-2.6047e-08,1.63774e-11,13870.1,26.5251], Tmin=(100,'K'), Tmax=(958.641,'K')), NASAPolynomial(coeffs=[18.9198,0.0130436,-3.93969e-06,7.20643e-10,-5.50768e-14,8878.54,-67.7866], Tmin=(958.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(=O)[CH]CC=O(7594)',
    structure = SMILES('[CH2]C([O])=CCC=O'),
    E0 = (-47.0015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,350,440,435,1725,459.256,466.743],'cm^-1')),
        HinderedRotor(inertia=(0.000774916,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933837,'amu*angstrom^2'), symmetry=1, barrier=(14.5098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0957514,'amu*angstrom^2'), symmetry=1, barrier=(14.4601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16182,0.0531124,-3.14975e-05,2.75342e-09,2.39156e-12,-5542.96,25.9129], Tmin=(100,'K'), Tmax=(1122.95,'K')), NASAPolynomial(coeffs=[14.5717,0.0208532,-9.12096e-06,1.76655e-09,-1.26488e-13,-9532.46,-44.6906], Tmin=(1122.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.0015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC([O])C=C[C]=O(7596)',
    structure = SMILES('CC([O])[CH]C=C=O'),
    E0 = (49.3434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,366.891,366.892,366.934],'cm^-1')),
        HinderedRotor(inertia=(0.240712,'amu*angstrom^2'), symmetry=1, barrier=(22.9862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240607,'amu*angstrom^2'), symmetry=1, barrier=(22.9869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428127,'amu*angstrom^2'), symmetry=1, barrier=(40.9155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972108,0.0639578,-5.76389e-05,2.68087e-08,-5.04426e-12,6045.78,23.2928], Tmin=(100,'K'), Tmax=(1265.17,'K')), NASAPolynomial(coeffs=[13.5436,0.024212,-1.05164e-05,1.97844e-09,-1.37825e-13,2864.72,-40.3145], Tmin=(1265.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.3434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C([O])C=CC[O](7417)',
    structure = SMILES('C=C([O])C=CC[O]'),
    E0 = (54.3769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,228.049,228.049,228.049,2964.37],'cm^-1')),
        HinderedRotor(inertia=(0.387588,'amu*angstrom^2'), symmetry=1, barrier=(14.3039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.387588,'amu*angstrom^2'), symmetry=1, barrier=(14.3039,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34185,0.0609452,-6.38806e-05,3.89658e-08,-9.95247e-12,6633.63,24.7107], Tmin=(100,'K'), Tmax=(935.551,'K')), NASAPolynomial(coeffs=[8.84339,0.0288723,-1.24575e-05,2.32235e-09,-1.60615e-13,5230,-10.9801], Tmin=(935.551,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.3769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])=CC=CO(5157)',
    structure = SMILES('C=C([O])[CH]C=CO'),
    E0 = (-79.88,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23817,'amu*angstrom^2'), symmetry=1, barrier=(28.4679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23448,'amu*angstrom^2'), symmetry=1, barrier=(28.3831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23815,'amu*angstrom^2'), symmetry=1, barrier=(28.4676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347664,0.0625043,-1.82904e-05,-4.25449e-08,2.72266e-11,-9459.18,24.5509], Tmin=(100,'K'), Tmax=(911.865,'K')), NASAPolynomial(coeffs=[24.9673,0.00211818,2.72483e-06,-6.50021e-10,4.22106e-14,-15928.6,-102.807], Tmin=(911.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-79.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=CC(C)=O(7600)',
    structure = SMILES('C=CC(C)=O'),
    E0 = (-126.494,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,180],'cm^-1')),
        HinderedRotor(inertia=(0.114507,'amu*angstrom^2'), symmetry=1, barrier=(2.63275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110505,'amu*angstrom^2'), symmetry=1, barrier=(2.54072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56497,0.0288349,-6.59621e-06,-4.30169e-09,1.71648e-12,-15159.9,16.5833], Tmin=(100,'K'), Tmax=(1456.22,'K')), NASAPolynomial(coeffs=[8.25968,0.0228095,-1.02958e-05,1.92717e-09,-1.31456e-13,-17838.1,-16.5318], Tmin=(1456.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-126.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC=CC=O(7601)',
    structure = SMILES('CC=CC=O'),
    E0 = (-117.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.626371,'amu*angstrom^2'), symmetry=1, barrier=(14.4015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625884,'amu*angstrom^2'), symmetry=1, barrier=(14.3903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4624,0.0331546,-1.4828e-05,1.43613e-09,2.81956e-13,-14059.8,15.7893], Tmin=(100,'K'), Tmax=(1739.36,'K')), NASAPolynomial(coeffs=[12.1692,0.0184725,-8.75534e-06,1.63406e-09,-1.09479e-13,-18592.3,-39.7353], Tmin=(1739.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'C=C(O)C=CC=O(5192)',
    structure = SMILES('C=C(O)C=CC=O'),
    E0 = (-243.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.939432,'amu*angstrom^2'), symmetry=1, barrier=(21.5994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941171,'amu*angstrom^2'), symmetry=1, barrier=(21.6394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937363,'amu*angstrom^2'), symmetry=1, barrier=(21.5518,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.402004,0.0694042,-6.75524e-05,3.19153e-08,-5.84626e-12,-29143.7,23.2376], Tmin=(100,'K'), Tmax=(1337.05,'K')), NASAPolynomial(coeffs=[18.8569,0.0141932,-5.61246e-06,1.03134e-09,-7.16018e-14,-34078.7,-71.1572], Tmin=(1337.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H)"""),
)

species(
    label = 'CC(=O)C=C=CO(7463)',
    structure = SMILES('CC(=O)C=C=CO'),
    E0 = (-194.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,180],'cm^-1')),
        HinderedRotor(inertia=(0.857163,'amu*angstrom^2'), symmetry=1, barrier=(19.7079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.855569,'amu*angstrom^2'), symmetry=1, barrier=(19.6712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854091,'amu*angstrom^2'), symmetry=1, barrier=(19.6372,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801099,0.0615574,-4.83841e-05,1.36454e-08,3.99662e-13,-23295.1,23.6111], Tmin=(100,'K'), Tmax=(1036.66,'K')), NASAPolynomial(coeffs=[16.3184,0.0171946,-6.63719e-06,1.23191e-09,-8.74905e-14,-27345.8,-55.8294], Tmin=(1036.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC(=O)[C]CC=O(7598)',
    structure = SMILES('CC(=O)[C]CC=O'),
    E0 = (67.471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,202.299,202.304,2083.49],'cm^-1')),
        HinderedRotor(inertia=(0.44347,'amu*angstrom^2'), symmetry=1, barrier=(12.8774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78522,'amu*angstrom^2'), symmetry=1, barrier=(51.8232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78514,'amu*angstrom^2'), symmetry=1, barrier=(51.823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168232,'amu*angstrom^2'), symmetry=1, barrier=(51.8234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02252,0.0504716,-2.55503e-05,-2.29367e-08,3.06168e-11,8179.23,33.5069], Tmin=(100,'K'), Tmax=(502.973,'K')), NASAPolynomial(coeffs=[4.56244,0.0392514,-1.88666e-05,3.69744e-09,-2.6315e-13,7810.15,21.8697], Tmin=(502.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = 'CC(=O)C[C]C=O(7599)',
    structure = SMILES('CC(=O)C[C]C=O'),
    E0 = (67.471,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,202.299,202.304,2083.49],'cm^-1')),
        HinderedRotor(inertia=(0.44347,'amu*angstrom^2'), symmetry=1, barrier=(12.8774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78522,'amu*angstrom^2'), symmetry=1, barrier=(51.8232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78514,'amu*angstrom^2'), symmetry=1, barrier=(51.823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168232,'amu*angstrom^2'), symmetry=1, barrier=(51.8234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02252,0.0504716,-2.55503e-05,-2.29367e-08,3.06168e-11,8179.23,33.5069], Tmin=(100,'K'), Tmax=(502.973,'K')), NASAPolynomial(coeffs=[4.56244,0.0392514,-1.88666e-05,3.69744e-09,-2.6315e-13,7810.15,21.8697], Tmin=(502.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (189.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (158.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (198.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (167.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (119.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (217.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (175.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (183.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (164.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (206.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-2.13613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (81.9323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-7.94108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (2.02494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (41.7666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (56.8969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (78.0088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (110.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (131.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (140.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-3.68206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (17.4298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (270.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (136.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (118.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-109.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (7.79664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (91.9822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (106.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (149.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (162.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (153.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (162.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (121.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (151.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (178.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (170.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-18.1046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (22.9131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (203.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (195.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-22.0282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (35.6705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (10.5697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (74.3167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (12.8844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (79.3501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (-54.9067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (59.5938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (200.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-39.2807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (-50.8409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (104.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (104.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction29',
    reactants = ['OH(5)(5)', 'C[C]=C=CC=O(7615)'],
    products = ['S(252)(251)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.11458e+06,'m^3/(mol*s)'), n=0.373849, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;Cd_allenic] + [O_pri_rad;Y_rad] for rate rule [O_pri_rad;Cd_allenic]
Euclidian distance = 2.0
family: R_Recombination
Ea raised from -5.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)(3)', 'CC(=O)[C]=CC=O(7584)'],
    products = ['S(252)(251)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7e+11,'cm^3/(mol*s)'), n=0.493, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2500,'K'), comment="""From training reaction 25 used for H_rad;O_rad/OneDe
Exact match found for rate rule [H_rad;O_rad/OneDe]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH3(15)(16)', 'O=CC=C=[C]O(7616)'],
    products = ['S(252)(251)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.705e+09,'cm^3/(mol*s)'), n=1.07, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 48 used for C_methyl;Cd_allenic
Exact match found for rate rule [Cd_allenic;C_methyl]
Euclidian distance = 0
family: R_Recombination
Ea raised from -9.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)(3)', '[CH2]C(O)=C=CC=O(7617)'],
    products = ['S(252)(251)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.156e+12,'cm^3/(mol*s)'), n=0.461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 15 used for C_rad/H2/Cd;H_rad
Exact match found for rate rule [C_rad/H2/Cd;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['HCO(14)(15)', '[CH]=C=C(C)O(7618)'],
    products = ['S(252)(251)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [Cd_rad;CO_pri_rad] for rate rule [Cd_allenic;CO_pri_rad]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)(3)', 'CC(O)=C=[C]C=O(7619)'],
    products = ['S(252)(251)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)(3)', 'CC(O)=C=C[C]=O(7620)'],
    products = ['S(252)(251)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CC([O])[C]=CC=O(7588)'],
    products = ['S(252)(251)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)[C]=CC=O(7621)'],
    products = ['S(252)(251)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC(O)=C=[C]C[O](7622)'],
    products = ['S(252)(251)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C[C](O)C=[C]C=O(7595)'],
    products = ['S(252)(251)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['CC(O)=[C]C[C]=O(7623)'],
    products = ['S(252)(251)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CC([O])=CC=C[O](7453)'],
    products = ['S(252)(251)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(O)=CC=C[O](5193)'],
    products = ['S(252)(251)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[C](O)C=C[C]=O(7597)'],
    products = ['S(252)(251)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CC([O])=[C]CC=O(7624)'],
    products = ['S(252)(251)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(O)=[C]CC=O(7625)'],
    products = ['S(252)(251)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC(O)[C]=C[C]=O(7626)'],
    products = ['S(252)(251)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['CC(=O)[C]=CC[O](7593)'],
    products = ['S(252)(251)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(O)=C=CC[O](7627)'],
    products = ['S(252)(251)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['CC([O])=[C]C=CO(7454)'],
    products = ['S(252)(251)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(O)=C=C[CH]O(7628)'],
    products = ['S(252)(251)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['CH2(S)(21)(22)', 'O=CC=C=CO(7629)'],
    products = ['S(252)(251)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(884144,'m^3/(mol*s)'), n=0.112, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cd_sec] for rate rule [carbene;Cd/H/NonDeO]
Euclidian distance = 1.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction54',
    reactants = ['CC(O)=C[C]C=O(7632)'],
    products = ['S(252)(251)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.14647e+14,'s^-1'), n=-1.07844, Ea=(56.8484,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CsJ2-C;CsJ2C;CH=C]
Euclidian distance = 0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CO(10)(11)', 'C=C=C(C)O(7630)'],
    products = ['S(252)(251)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction27',
    reactants = ['S(252)(251)'],
    products = ['S(247)(246)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(104,'s^-1'), n=3.21, Ea=(82.0482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH3;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction53',
    reactants = ['CC(O)=C=C=CO(7631)'],
    products = ['S(252)(251)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CH3(15)(16)', 'S(780)(779)'],
    products = ['S(247)(246)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.2e+13,'cm^3/(mol*s)','+|-',8.4e+12), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""Estimated using template [CO_sec_rad;C_methyl] for rate rule [CO_rad/OneDe;C_methyl]
Euclidian distance = 1.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)(3)', '[CH2]C(=O)C=CC=O(2803)'],
    products = ['S(247)(246)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.83648e+06,'m^3/(mol*s)'), n=0.3308, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_pri_rad;H_rad] for rate rule [C_rad/H2/CO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH3CO(55)(54)', 'CHCHCHO(2768)'],
    products = ['S(247)(246)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.6647e+07,'m^3/(mol*s)'), n=-0.0666667, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;CO_rad/NonDe] + [Cd_pri_rad;CO_rad] for rate rule [Cd_pri_rad;CO_rad/NonDe]
Euclidian distance = 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', 'CC(=O)[C]=CC=O(7584)'],
    products = ['S(247)(246)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['HCO(14)(15)', '[CH]=CC(C)=O(7465)'],
    products = ['S(247)(246)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.81e+13,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 91 used for Cd_pri_rad;CO_pri_rad
Exact match found for rate rule [CO_pri_rad;Cd_pri_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)(3)', 'CC(=O)C=[C]C=O(7585)'],
    products = ['S(247)(246)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad/OneDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)(3)', 'CC(=O)C=C[C]=O(7586)'],
    products = ['S(247)(246)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;CO_sec_rad] for rate rule [H_rad;CO_rad/OneDe]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([O])C=CC=O(7587)'],
    products = ['S(247)(246)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC([O])[C]=CC=O(7588)'],
    products = ['S(247)(246)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC(=O)C=[C]C[O](7589)'],
    products = ['S(247)(246)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CC(=O)[CH]C[C]=O(7590)'],
    products = ['S(247)(246)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(=O)C[CH]C=O(7591)'],
    products = ['S(247)(246)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CC([O])C=[C]C=O(7592)'],
    products = ['S(247)(246)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CC(=O)[C]=CC[O](7593)'],
    products = ['S(247)(246)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(=O)[CH]CC=O(7594)'],
    products = ['S(247)(246)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[C](O)C=[C]C=O(7595)'],
    products = ['S(247)(246)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(60.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC([O])=[C]C=CO(7454)'],
    products = ['S(247)(246)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])C=C[C]=O(7596)'],
    products = ['S(247)(246)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[C](O)C=C[C]=O(7597)'],
    products = ['S(247)(246)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C=CC[O](7417)'],
    products = ['S(247)(246)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([O])=CC=CO(5157)'],
    products = ['S(247)(246)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CO(10)(11)', 'C=CC(C)=O(7600)'],
    products = ['S(247)(246)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO;R_H] for rate rule [CO;Cd_pri]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CO(10)(11)', 'CC=CC=O(7601)'],
    products = ['S(247)(246)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(538,'cm^3/(mol*s)'), n=3.29, Ea=(437.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO;R_R'] for rate rule [CO;C_methyl_Cd_pri]
Euclidian distance = 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C=CC=O(5192)'],
    products = ['S(247)(246)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(204.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC(=O)C=C=CO(7463)'],
    products = ['S(247)(246)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(143.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CC(=O)[C]CC=O(7598)'],
    products = ['S(247)(246)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC(=O)C[C]C=O(7599)'],
    products = ['S(247)(246)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.90176e+13,'s^-1'), n=-0.332469, Ea=(37.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [singletcarbene_CH;CsJ2C;CH2(C)] + [CsJ2-C;CsJ2C;CH] for rate rule [CsJ2-C;CsJ2C;CH2(C)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Singlet_Carbene_Intra_Disproportionation"""),
)

network(
    label = '652',
    isomers = [
        'S(247)(246)',
        'S(252)(251)',
    ],
    reactants = [
        ('CH3(15)(16)', 'S(780)(779)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '652',
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

