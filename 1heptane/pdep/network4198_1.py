species(
    label = 'C=[C]OC[C]=CC(20247)',
    structure = SMILES('C=[C]OC[C]=CC'),
    E0 = (365.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,272.948,272.948,272.949,272.951],'cm^-1')),
        HinderedRotor(inertia=(0.00226279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272381,'amu*angstrom^2'), symmetry=1, barrier=(14.4001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272382,'amu*angstrom^2'), symmetry=1, barrier=(14.4001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272383,'amu*angstrom^2'), symmetry=1, barrier=(14.4001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784898,0.0658741,-5.59174e-05,2.5085e-08,-4.56791e-12,44092.7,28.7358], Tmin=(100,'K'), Tmax=(1306.52,'K')), NASAPolynomial(coeffs=[13.7228,0.0262639,-1.04413e-05,1.88032e-09,-1.27742e-13,40712,-37.1411], Tmin=(1306.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'CH2CO(28)',
    structure = SMILES('C=C=O'),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C=[C]OC=C=CC(25009)',
    structure = SMILES('C=[C]OC=C=CC'),
    E0 = (314.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.761513,'amu*angstrom^2'), symmetry=1, barrier=(17.5087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762266,'amu*angstrom^2'), symmetry=1, barrier=(17.526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762472,'amu*angstrom^2'), symmetry=1, barrier=(17.5307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975509,0.0621175,-5.29174e-05,2.38235e-08,-4.35509e-12,37891.3,27.691], Tmin=(100,'K'), Tmax=(1301.33,'K')), NASAPolynomial(coeffs=[13.1163,0.024799,-9.90127e-06,1.78631e-09,-1.21479e-13,34731.5,-34.0791], Tmin=(1301.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OCC#CC(25010)',
    structure = SMILES('C=[C]OCC#CC'),
    E0 = (287.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,363.373,363.374,363.374,363.374],'cm^-1')),
        HinderedRotor(inertia=(0.868631,'amu*angstrom^2'), symmetry=1, barrier=(81.3898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155668,'amu*angstrom^2'), symmetry=1, barrier=(14.5859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155668,'amu*angstrom^2'), symmetry=1, barrier=(14.5859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868628,'amu*angstrom^2'), symmetry=1, barrier=(81.3898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313725,0.0650007,-5.74069e-05,2.65847e-08,-4.73027e-12,34732.6,29.1924], Tmin=(100,'K'), Tmax=(1571.1,'K')), NASAPolynomial(coeffs=[15.5638,0.0176267,-4.0162e-06,4.66512e-10,-2.32191e-14,30995.6,-47.9127], Tmin=(1571.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=CJO)"""),
)

species(
    label = 'C#COC[C]=CC(25011)',
    structure = SMILES('C#COC[C]=CC'),
    E0 = (335.195,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.889033,'amu*angstrom^2'), symmetry=1, barrier=(20.4406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888885,'amu*angstrom^2'), symmetry=1, barrier=(20.4372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889142,'amu*angstrom^2'), symmetry=1, barrier=(20.4431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888875,'amu*angstrom^2'), symmetry=1, barrier=(20.437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.524357,0.0703522,-6.72452e-05,3.28028e-08,-6.33448e-12,40444.8,23.9993], Tmin=(100,'K'), Tmax=(1256.09,'K')), NASAPolynomial(coeffs=[16.0394,0.0209442,-8.24249e-06,1.48693e-09,-1.01611e-13,36547.2,-54.3892], Tmin=(1256.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.195,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_S)"""),
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
    label = 'C#CCO[C]=C(25012)',
    structure = SMILES('C#CCO[C]=C'),
    E0 = (329.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2175,525,1685,370,750,770,3400,2100,2950,3100,1380,975,1025,1650,246.441,248.113,249.108],'cm^-1')),
        HinderedRotor(inertia=(0.535391,'amu*angstrom^2'), symmetry=1, barrier=(23.1973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.903095,'amu*angstrom^2'), symmetry=1, barrier=(39.1059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528489,'amu*angstrom^2'), symmetry=1, barrier=(23.1845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83228,0.0556813,-5.42551e-05,2.63745e-08,-4.82784e-12,39785.3,23.9072], Tmin=(100,'K'), Tmax=(1543.83,'K')), NASAPolynomial(coeffs=[15.2661,0.00987263,-1.57459e-06,9.66613e-11,-1.06292e-15,36331,-48.7495], Tmin=(1543.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]O[CH]C=CC(20249)',
    structure = SMILES('C=[C]O[CH]C=CC'),
    E0 = (238.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76025,0.0618426,-4.52392e-05,1.68737e-08,-2.53627e-12,28834.6,29.3024], Tmin=(100,'K'), Tmax=(1567.56,'K')), NASAPolynomial(coeffs=[15.4099,0.0244606,-9.46817e-06,1.66062e-09,-1.10038e-13,24241.8,-47.9588], Tmin=(1567.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=[C]OCC=[C]C(20246)',
    structure = SMILES('C=[C]OCC=[C]C'),
    E0 = (365.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,272.948,272.948,272.949,272.951],'cm^-1')),
        HinderedRotor(inertia=(0.00226279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272381,'amu*angstrom^2'), symmetry=1, barrier=(14.4001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272382,'amu*angstrom^2'), symmetry=1, barrier=(14.4001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.272383,'amu*angstrom^2'), symmetry=1, barrier=(14.4001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784898,0.0658741,-5.59174e-05,2.5085e-08,-4.56791e-12,44092.7,28.7358], Tmin=(100,'K'), Tmax=(1306.52,'K')), NASAPolynomial(coeffs=[13.7228,0.0262639,-1.04413e-05,1.88032e-09,-1.27742e-13,40712,-37.1411], Tmin=(1306.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC[C]=CC(25013)',
    structure = SMILES('[CH]=COC[C]=CC'),
    E0 = (372.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.746414,'amu*angstrom^2'), symmetry=1, barrier=(17.1615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745338,'amu*angstrom^2'), symmetry=1, barrier=(17.1368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745314,'amu*angstrom^2'), symmetry=1, barrier=(17.1362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.744785,'amu*angstrom^2'), symmetry=1, barrier=(17.1241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44493,0.0683515,-4.84933e-05,7.91368e-09,3.61347e-12,44994,27.2283], Tmin=(100,'K'), Tmax=(986.688,'K')), NASAPolynomial(coeffs=[17.34,0.0205402,-7.24795e-06,1.28815e-09,-9.00321e-14,40653.3,-59.1553], Tmin=(986.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=CO[CH][C]=CC(25014)',
    structure = SMILES('C=CO[CH][C]=CC'),
    E0 = (236.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695502,0.0628493,-3.79876e-05,2.50338e-09,3.99556e-12,28608.9,26.8627], Tmin=(100,'K'), Tmax=(1029.77,'K')), NASAPolynomial(coeffs=[15.4447,0.0243752,-9.35418e-06,1.7108e-09,-1.19861e-13,24573.5,-49.5706], Tmin=(1029.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OC[CH]C=C(19795)',
    structure = SMILES('C=[C]OC[CH]C=C'),
    E0 = (251.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3355.68,'J/mol'), sigma=(5.84875,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=524.15 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.389317,0.0671199,-3.70254e-05,-6.52596e-09,8.98886e-12,30374.5,26.9345], Tmin=(100,'K'), Tmax=(985.206,'K')), NASAPolynomial(coeffs=[18.5341,0.0205568,-7.40134e-06,1.35417e-09,-9.71101e-14,25483.7,-67.0094], Tmin=(985.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OCC=CC(20250)',
    structure = SMILES('[CH]=[C]OCC=CC'),
    E0 = (374.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,326.515,326.526,326.533],'cm^-1')),
        HinderedRotor(inertia=(0.183354,'amu*angstrom^2'), symmetry=1, barrier=(13.8727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183364,'amu*angstrom^2'), symmetry=1, barrier=(13.8726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183367,'amu*angstrom^2'), symmetry=1, barrier=(13.8726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183359,'amu*angstrom^2'), symmetry=1, barrier=(13.8727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.488289,0.067706,-5.74859e-05,2.51609e-08,-4.38203e-12,45220.4,29.737], Tmin=(100,'K'), Tmax=(1384.49,'K')), NASAPolynomial(coeffs=[16.2441,0.0221851,-8.16727e-06,1.41275e-09,-9.37793e-14,40857.7,-51.4015], Tmin=(1384.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = 'C=COC[C]=[C]C(25015)',
    structure = SMILES('C=COC[C]=[C]C'),
    E0 = (363.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,265.163,265.165,265.171,265.245],'cm^-1')),
        HinderedRotor(inertia=(0.317608,'amu*angstrom^2'), symmetry=1, barrier=(15.8555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317537,'amu*angstrom^2'), symmetry=1, barrier=(15.8556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317705,'amu*angstrom^2'), symmetry=1, barrier=(15.8555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317637,'amu*angstrom^2'), symmetry=1, barrier=(15.8555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.447313,0.069868,-5.80844e-05,2.1559e-08,-2.11626e-12,43879.3,27.291], Tmin=(100,'K'), Tmax=(1045.68,'K')), NASAPolynomial(coeffs=[16.1831,0.022419,-8.30183e-06,1.47572e-09,-1.01277e-13,39891.5,-52.6598], Tmin=(1045.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]COC=C(20248)',
    structure = SMILES('[CH2]C=[C]COC=C'),
    E0 = (277.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,314.399,314.4,314.402,314.403],'cm^-1')),
        HinderedRotor(inertia=(0.303549,'amu*angstrom^2'), symmetry=1, barrier=(21.2923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303541,'amu*angstrom^2'), symmetry=1, barrier=(21.2923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303553,'amu*angstrom^2'), symmetry=1, barrier=(21.2923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303548,'amu*angstrom^2'), symmetry=1, barrier=(21.2923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544018,0.062026,-2.21586e-05,-2.43747e-08,1.62461e-11,33496.7,26.7503], Tmin=(100,'K'), Tmax=(954.966,'K')), NASAPolynomial(coeffs=[19.1883,0.0176648,-5.4645e-06,9.60729e-10,-6.99606e-14,28397.6,-70.3921], Tmin=(954.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C][O](173)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39563,0.0101365,2.30741e-06,-8.97566e-09,3.68242e-12,19290.3,10.0703], Tmin=(100,'K'), Tmax=(1068.9,'K')), NASAPolynomial(coeffs=[6.35055,0.00638951,-2.69368e-06,5.4221e-10,-4.02476e-14,18240.9,-6.33602], Tmin=(1068.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CsCJ=O) + radical(CJC=O)"""),
)

species(
    label = 'C=COC=C=CC(25016)',
    structure = SMILES('C=COC=C=CC'),
    E0 = (74.3669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708189,0.0610794,-2.8829e-05,-1.07915e-08,9.79334e-12,9072.89,25.5199], Tmin=(100,'K'), Tmax=(979.139,'K')), NASAPolynomial(coeffs=[16.5553,0.02181,-7.68825e-06,1.38072e-09,-9.76571e-14,4748.69,-56.8337], Tmin=(979.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.3669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
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
    label = 'C=[C]CO[C]=C(25017)',
    structure = SMILES('C=[C]CO[C]=C'),
    E0 = (401.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1670,1700,300,440,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,402.886,402.887,402.889,402.89],'cm^-1')),
        HinderedRotor(inertia=(0.13094,'amu*angstrom^2'), symmetry=1, barrier=(15.0822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130938,'amu*angstrom^2'), symmetry=1, barrier=(15.0822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130939,'amu*angstrom^2'), symmetry=1, barrier=(15.0822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38371,0.0508909,-3.76247e-05,9.99712e-09,5.44099e-13,48405.9,24.7791], Tmin=(100,'K'), Tmax=(1024.69,'K')), NASAPolynomial(coeffs=[12.9457,0.0177132,-6.55907e-06,1.17246e-09,-8.09999e-14,45408.7,-34.3456], Tmin=(1024.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OCC1=CC(25018)',
    structure = SMILES('C=C1OCC1=CC'),
    E0 = (6.94949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24845,0.0457843,1.1732e-05,-4.91138e-08,2.22783e-11,948.07,18.9856], Tmin=(100,'K'), Tmax=(977.642,'K')), NASAPolynomial(coeffs=[15.2425,0.0237671,-8.55467e-06,1.58953e-09,-1.15649e-13,-3472.23,-56.8243], Tmin=(977.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.94949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane)"""),
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
    label = 'CC=[C]C[O](24539)',
    structure = SMILES('CC=[C]C[O]'),
    E0 = (280.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,2871.59],'cm^-1')),
        HinderedRotor(inertia=(0.255953,'amu*angstrom^2'), symmetry=1, barrier=(5.88487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25515,'amu*angstrom^2'), symmetry=1, barrier=(5.86641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4896,0.0397772,-5.33783e-05,5.63765e-08,-2.35894e-11,33785.4,19.4797], Tmin=(100,'K'), Tmax=(821.312,'K')), NASAPolynomial(coeffs=[0.0204228,0.0358527,-1.70807e-05,3.26805e-09,-2.26308e-13,34729,34.181], Tmin=(821.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH2]O[C]=C(2374)',
    structure = SMILES('[CH2]O[C]=C'),
    E0 = (283.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,607.226,610.81],'cm^-1')),
        HinderedRotor(inertia=(0.0541951,'amu*angstrom^2'), symmetry=1, barrier=(14.1831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545119,'amu*angstrom^2'), symmetry=1, barrier=(14.2038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71895,0.0199743,1.11331e-05,-3.28562e-08,1.53117e-11,34186,17.0659], Tmin=(100,'K'), Tmax=(932.118,'K')), NASAPolynomial(coeffs=[10.6429,0.00688602,-1.46284e-06,2.25572e-10,-1.7517e-14,31800.2,-25.4794], Tmin=(932.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=COCJ)"""),
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
    E0 = (365.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (538.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (517.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (562.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (365.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (496.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (504.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (599.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (478.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (507.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (527.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (683.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (606.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (508.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (521.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (458.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (821.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (373.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (681.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (847.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['CH2CO(28)', 'CH3CHCCH2(18175)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=[C]OC=C=CC(25009)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=[C]OCC#CC(25010)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.02e+09,'cm^3/(mol*s)'), n=1.64, Ea=(18.4933,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2702 used for Ct-Cs_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C#COC[C]=CC(25011)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2CO(28)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(65.3742,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from 60.1 to 65.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3(17)', 'C#CCO[C]=C(25012)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(178000,'cm^3/(mol*s)'), n=2.41, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2271 used for Ct-H_Ct-Cs;CsJ-HHH
Exact match found for rate rule [Ct-H_Ct-Cs;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['C=[C]O[CH]C=CC(20249)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]OCC=[C]C(20246)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=COC[C]=CC(25013)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['C=CO[CH][C]=CC(25014)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['C=[C]OC[CH]C=C(19795)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OCC=CC(20250)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=COC[C]=[C]C(25015)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.24775e+08,'s^-1'), n=1.35208, Ea=(243.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_single;Cd_H_out_doubleC] + [RnH;Cd_rad_out_Cs;XH_out] for rate rule [R5HJ_1;Cd_rad_out_Cs;Cd_H_out_doubleC]
Euclidian distance = 2.82842712475
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['[CH2]C=[C]COC=C(20248)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][O](173)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['C=COC=C=CC(25016)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(S)(23)', 'C=[C]CO[C]=C(25017)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]OC[C]=CC(20247)'],
    products = ['C=C1OCC1=CC(25018)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H2CC(41)', 'CC=[C]C[O](24539)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]O[C]=C(2374)', '[C]=CC(24199)'],
    products = ['C=[C]OC[C]=CC(20247)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4198',
    isomers = [
        'C=[C]OC[C]=CC(20247)',
    ],
    reactants = [
        ('CH2CO(28)', 'CH3CHCCH2(18175)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4198',
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

