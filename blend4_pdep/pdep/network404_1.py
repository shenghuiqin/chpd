species(
    label = 'S(554)(553)',
    structure = SMILES('[O]OC1[CH]C2=CC=C(C1)O2'),
    E0 = (106.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4741.84,'J/mol'), sigma=(7.48196,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=740.66 K, Pc=25.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12027,0.0532238,-4.94984e-06,-2.12056e-08,7.88641e-12,12884.8,25.7164], Tmin=(100,'K'), Tmax=(1314.12,'K')), NASAPolynomial(coeffs=[16.8089,0.0359897,-2.01151e-05,4.16103e-09,-3.01515e-13,6126.14,-64.2848], Tmin=(1314.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)(2)',
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
    label = 'C7H6O(490)(489)',
    structure = SMILES('C1=CC2=CC=C(C1)O2'),
    E0 = (41.2019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4154,'J/mol'), sigma=(6.51034,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=648.85 K, Pc=34.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20282,0.0266771,3.98414e-05,-5.77521e-08,1.86939e-11,5030.4,17.6585], Tmin=(100,'K'), Tmax=(1194.95,'K')), NASAPolynomial(coeffs=[12.8261,0.0324323,-1.92456e-05,4.14704e-09,-3.09659e-13,-458.206,-47.8272], Tmin=(1194.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.2019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan)"""),
)

species(
    label = 'O2(S)(162)(161)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1857.18,'J/mol'), sigma=(4.34667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=290.09 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]1C=C2CC3[CH]C1(OO3)O2(6936)',
    structure = SMILES('[CH]1C2C[C]3C=CC1(OO2)O3'),
    E0 = (309.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58623,0.0523956,-1.89211e-06,-2.22118e-08,7.11241e-12,37335,17.8805], Tmin=(100,'K'), Tmax=(1494.74,'K')), NASAPolynomial(coeffs=[22.044,0.0363611,-2.46491e-05,5.26453e-09,-3.81016e-13,26894.6,-103.506], Tmin=(1494.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_ane) + polycyclic(s1_5_5_ene_1) + polycyclic(s3_5_6_ene_5) - ring(12dioxolane) - ring(25dihydrofuran) - ring(Oxane) + radical(C2CsJOCs) + radical(CCJCOOH)"""),
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
    label = '[O]OC1=CC2=CC=C(C1)O2(6937)',
    structure = SMILES('[O]OC1=CC2=CC=C(C1)O2'),
    E0 = (103.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48275,0.0530868,-2.06193e-05,-4.24223e-09,2.50649e-12,12497.4,25.1327], Tmin=(100,'K'), Tmax=(1556.99,'K')), NASAPolynomial(coeffs=[20.8941,0.0260891,-1.66438e-05,3.48993e-09,-2.49869e-13,3680.46,-86.0131], Tmin=(1556.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(ROOJ)"""),
)

species(
    label = '[O]O[C]1CC2=CC=C(C1)O2(6938)',
    structure = SMILES('[O]O[C]1CC2=CC=C(C1)O2'),
    E0 = (176.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86162,0.0470221,-4.66357e-06,-1.42582e-08,4.59397e-12,21266.6,26.1663], Tmin=(100,'K'), Tmax=(1529.17,'K')), NASAPolynomial(coeffs=[16.185,0.0363935,-2.05643e-05,4.15155e-09,-2.92228e-13,13748.2,-59.2794], Tmin=(1529.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C2CsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'OO[C]1[CH]C2=CC=C(C1)O2(6939)',
    structure = SMILES('OOC1=C[C]2C=C[C](C1)O2'),
    E0 = (312.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28216,0.0630576,-3.80172e-05,9.85411e-09,-9.79256e-13,37647.7,28.1041], Tmin=(100,'K'), Tmax=(2302.9,'K')), NASAPolynomial(coeffs=[24.8014,0.0222052,-1.14074e-05,2.15069e-09,-1.42965e-13,26815.4,-104.981], Tmin=(2302.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs) + radical(C2CsJOCs)"""),
)

species(
    label = '[O]OC1CC2=[C]C=C(C1)O2(6940)',
    structure = SMILES('[O]OC1CC2=[C]C=C(C1)O2'),
    E0 = (188.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51105,0.0527299,-1.5035e-05,-7.23793e-09,3.03484e-12,22735.9,26.7221], Tmin=(100,'K'), Tmax=(1545.36,'K')), NASAPolynomial(coeffs=[17.4379,0.0339747,-1.86407e-05,3.72656e-09,-2.6107e-13,15130.3,-65.7293], Tmin=(1545.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = 'OOC1[CH]C2=CC=C([CH]1)O2(6941)',
    structure = SMILES('OOC1[CH]C2=CC=C([CH]1)O2'),
    E0 = (71.1143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.559843,0.0617397,-1.1561e-05,-2.18695e-08,9.16718e-12,8688.08,24.1642], Tmin=(100,'K'), Tmax=(1259.39,'K')), NASAPolynomial(coeffs=[20.023,0.0343028,-1.98312e-05,4.18489e-09,-3.07785e-13,1059.23,-85.0475], Tmin=(1259.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.1143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CCJCO) + radical(C=CCJCO)"""),
)

species(
    label = 'OOC1[CH]C2=C[C]=C(C1)O2(6942)',
    structure = SMILES('OOC1[CH]C2=C[C]=C(C1)O2'),
    E0 = (153.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847334,0.0623773,-2.51936e-05,-3.88621e-09,2.85124e-12,18543.8,25.5462], Tmin=(100,'K'), Tmax=(1470.11,'K')), NASAPolynomial(coeffs=[20.2126,0.0327728,-1.85428e-05,3.77983e-09,-2.68945e-13,10355.3,-83.8271], Tmin=(1470.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CCJCO) + radical(C=CJC=C)"""),
)

species(
    label = 'OOC1[CH]C2=[C]C=C(C1)O2(6943)',
    structure = SMILES('OOC1[CH]C2=[C]C=C(C1)O2'),
    E0 = (153.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847334,0.0623773,-2.51936e-05,-3.88621e-09,2.85124e-12,18543.8,25.5462], Tmin=(100,'K'), Tmax=(1470.11,'K')), NASAPolynomial(coeffs=[20.2126,0.0327728,-1.85428e-05,3.77983e-09,-2.68945e-13,10355.3,-83.8271], Tmin=(1470.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CJC=C) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1C=C2C=C[C](C1)O2(5783)',
    structure = SMILES('C1=C[C]2C=C[C](C1)O2'),
    E0 = (402.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97461,0.0358572,1.03846e-05,-2.75166e-08,9.44083e-12,48466,21.5381], Tmin=(100,'K'), Tmax=(1209.37,'K')), NASAPolynomial(coeffs=[9.6968,0.0338035,-1.62004e-05,3.19765e-09,-2.28154e-13,44880.6,-24.2863], Tmin=(1209.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(C2CsJOCs) + radical(C2CsJOCs)"""),
)

species(
    label = '[O]OC1CC2=CC3[C](O2)C13(6944)',
    structure = SMILES('[O]OC1CC2=CC3[C](O2)C13'),
    E0 = (299.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22396,0.0627237,-3.02291e-05,-2.16288e-10,1.687e-12,36118.4,22.565], Tmin=(100,'K'), Tmax=(1668.2,'K')), NASAPolynomial(coeffs=[28.245,0.0232142,-1.74355e-05,3.75555e-09,-2.69658e-13,23585.4,-132.167], Tmin=(1668.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.489,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(2,3-Dihydrofuran) + polycyclic(s2_3_6_ane) + polycyclic(s2_3_5_ene_1) - ring(Oxane) - ring(2,3-Dihydrofuran) - ring(Cyclopropane) + radical(ROOJ) + radical(C2CsJOC(O))"""),
)

species(
    label = '[O]OC1C[C]2OC3=CC2C31(6945)',
    structure = SMILES('[O]OC1C[C]2OC3=CC2C31'),
    E0 = (552.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949728,0.0749133,-5.35666e-05,1.68312e-08,-2.09765e-12,66579.7,18.1706], Tmin=(100,'K'), Tmax=(1784.24,'K')), NASAPolynomial(coeffs=[19.503,0.0333192,-1.85984e-05,3.76553e-09,-2.66927e-13,59959,-82.0802], Tmin=(1784.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_5_5_ene_1) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + polycyclic(s3_4_5_ene_3) - ring(2,3-Dihydrofuran) - ring(Cyclopentane) - ring(Cyclobutene) + radical(ROOJ) + radical(C2CsJOC(O))"""),
)

species(
    label = '[O]OC1CC23[CH]C=C(O2)C13(6946)',
    structure = SMILES('[O]OC1CC23C=C[C](O2)C13'),
    E0 = (477.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50949,0.0576128,-3.05185e-05,5.37012e-09,-2.17515e-13,57495,16.2525], Tmin=(100,'K'), Tmax=(2333.05,'K')), NASAPolynomial(coeffs=[40.315,0.00758253,-7.85942e-06,1.61198e-09,-1.05917e-13,35830.2,-206.79], Tmin=(2333.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_4_5_ene_1) + polycyclic(s1_4_5_ene_1) + polycyclic(s2_4_4_ane) - ring(Cyclobutane) - ring(25dihydrofuran) - ring(Oxetane) + radical(ROOJ) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH]1C2=CC3OOC1C[C]3O2(6947)',
    structure = SMILES('[CH]1C2=CC3C[C](O2)C1OO3'),
    E0 = (188.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96088,0.0481001,1.46762e-05,-3.76207e-08,1.09735e-11,22763.4,12.7371], Tmin=(100,'K'), Tmax=(1494.27,'K')), NASAPolynomial(coeffs=[25.5163,0.0390644,-3.04802e-05,6.7189e-09,-4.9241e-13,9692.83,-130.545], Tmin=(1494.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_6_6_ene_1) + polycyclic(s2_5_6_ane) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(Tetrahydrofuran) - ring(3,4-Dihydro-2H-pyran) - ring(12dioxane) - ring(Tetrahydrofuran) + radical(C2CsJOC(O)) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1[C]2OC3=CC2OOC1C3(6948)',
    structure = SMILES('[CH]1[C]2OC3=CC2OOC1C3'),
    E0 = (303.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46322,0.0413324,2.15118e-05,-3.98177e-08,1.10571e-11,36553.6,15.3848], Tmin=(100,'K'), Tmax=(1519.2,'K')), NASAPolynomial(coeffs=[23.7554,0.0397802,-3.07761e-05,6.74547e-09,-4.9194e-13,23793.9,-116.944], Tmin=(1519.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_6_6_ane) + polycyclic(s2_5_6_ene_6) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(2,3-Dihydrofuran) - ring(Oxane) - ring(2,3-Dihydrofuran) - ring(12dioxane) + radical(C2CsJOC(O)) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]1C2=C[CH]C3(CC1OO3)O2(6949)',
    structure = SMILES('[CH]1[C]2C=CC3(CC1OO3)O2'),
    E0 = (309.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58624,0.0523955,-1.89191e-06,-2.2212e-08,7.11246e-12,37335,17.8805], Tmin=(100,'K'), Tmax=(1494.74,'K')), NASAPolynomial(coeffs=[22.0439,0.0363614,-2.46493e-05,5.26456e-09,-3.81018e-13,26894.7,-103.505], Tmin=(1494.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_ene_5) + polycyclic(s1_5_5_ene_1) + polycyclic(s3_5_6_ane) - ring(12dioxolane) - ring(25dihydrofuran) - ring(Oxane) + radical(CCJCOOH) + radical(C2CsJOCs)"""),
)

species(
    label = '[CH2]C1=CC=C(C=CO[O])O1(6950)',
    structure = SMILES('[CH2]C1=CC=C(C=CO[O])O1'),
    E0 = (184.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3100,440,815,1455,1000,2750,3150,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.10686,0.0775987,-6.74514e-05,2.97906e-08,-5.24754e-12,22313.6,33.3458], Tmin=(100,'K'), Tmax=(1362.01,'K')), NASAPolynomial(coeffs=[17.7135,0.025891,-1.0505e-05,1.91689e-09,-1.31273e-13,17517.5,-57.0359], Tmin=(1362.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Furan) + radical(C=C(O)CJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=C=CC=C([O])C1(6951)',
    structure = SMILES('[O]OC1C=[C]C=CC(=O)C1'),
    E0 = (218.616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.413631,0.0338398,0.000131277,-2.33515e-07,1.03821e-10,26464.8,31.4294], Tmin=(100,'K'), Tmax=(905.458,'K')), NASAPolynomial(coeffs=[38.3022,-0.0147732,1.50608e-05,-3.08535e-09,2.02163e-13,14735,-174.483], Tmin=(905.458,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'OOC1=CC2=CC=C(C1)O2(6952)',
    structure = SMILES('OOC1=CC2=CC=C(C1)O2'),
    E0 = (-48.8619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07328,0.0583559,-2.15449e-05,-6.67295e-09,3.55194e-12,-5766.77,25.8451], Tmin=(100,'K'), Tmax=(1470.79,'K')), NASAPolynomial(coeffs=[21.3629,0.0291571,-1.82636e-05,3.85008e-09,-2.782e-13,-14545.3,-89.4217], Tmin=(1470.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.8619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan)"""),
)

species(
    label = '[O]OC1[CH]C2C=CC(=C1)O2(6953)',
    structure = SMILES('[O]OC1[CH]C2C=CC(=C1)O2'),
    E0 = (222.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.43358,0.0576989,4.91945e-06,-4.74809e-08,2.00338e-11,26950.2,28.3546], Tmin=(100,'K'), Tmax=(1106.27,'K')), NASAPolynomial(coeffs=[22.9092,0.0242114,-1.44585e-05,3.23753e-09,-2.50491e-13,19053.7,-95.5611], Tmin=(1106.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(25dihydrofuran) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'O(4)(4)',
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
    label = 'C1C=C2OC=1CC1OC21(6954)',
    structure = SMILES('C1C=C2OC=1CC1OC21'),
    E0 = (-19.2472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50268,0.0250969,4.68143e-05,-6.00996e-08,1.76196e-11,-2255.38,21.3325], Tmin=(100,'K'), Tmax=(1308.36,'K')), NASAPolynomial(coeffs=[13.8241,0.0395352,-2.5974e-05,5.64334e-09,-4.17688e-13,-9416.16,-52.3736], Tmin=(1308.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.2472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + polycyclic(s2_3_6_ane) - ring(Oxane)"""),
)

species(
    label = '[O]O[CH]C1CC2=CC=C1O2(6955)',
    structure = SMILES('[O]O[CH]C1CC2=CC=C1O2'),
    E0 = (205.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40714,0.0532081,-1.47589e-05,-8.66142e-09,3.59967e-12,24866.7,29.0995], Tmin=(100,'K'), Tmax=(1495.78,'K')), NASAPolynomial(coeffs=[17.7742,0.0333319,-1.87862e-05,3.81232e-09,-2.69977e-13,17297.6,-65.3864], Tmin=(1495.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_5_ane) - ring(Tetrahydrofuran) - ring(Tetrahydrofuran) + ring(Tetrahydrofuran) + ring(Furan) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = 'HO2(8)(9)',
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
    label = 'C1=C[C]2C=CC(=C1)O2(5776)',
    structure = SMILES('[CH]1C=CC2=CC=C1O2'),
    E0 = (158.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97064,0.0308161,3.1354e-05,-5.28177e-08,1.77798e-11,19101.5,14.9964], Tmin=(100,'K'), Tmax=(1187.38,'K')), NASAPolynomial(coeffs=[14.8064,0.0286413,-1.77762e-05,3.89431e-09,-2.93431e-13,13158.4,-61.3232], Tmin=(1187.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(C=CCJCO)"""),
)

species(
    label = 'C1C=C2OC=1CC1OOC21(6956)',
    structure = SMILES('C1C=C2OC=1CC1OOC21'),
    E0 = (20.2103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1121,0.0301595,5.13315e-05,-6.87393e-08,2.0725e-11,2507.24,22.0853], Tmin=(100,'K'), Tmax=(1278.24,'K')), NASAPolynomial(coeffs=[15.4178,0.0431621,-2.80466e-05,6.10216e-09,-4.53268e-13,-5358.11,-62.8339], Tmin=(1278.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.2103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + polycyclic(s2_4_6_ane) - ring(Oxane)"""),
)

species(
    label = '[O]OC1[CH]C2=CC2C(=O)C1(6957)',
    structure = SMILES('[O]OC1[CH]C2=CC2C(=O)C1'),
    E0 = (297.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.291168,0.0705604,-4.48203e-05,1.10472e-08,-2.90156e-13,35922,26.57], Tmin=(100,'K'), Tmax=(1268.95,'K')), NASAPolynomial(coeffs=[16.7958,0.031307,-1.3518e-05,2.53422e-09,-1.75733e-13,30704.9,-61.0387], Tmin=(1268.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_5) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1[CH]C(=O)C2C=C2C1(6958)',
    structure = SMILES('[O]OC1[CH]C(=O)C2C=C2C1'),
    E0 = (380.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.754509,0.0609195,-2.73801e-05,-2.81952e-09,3.81131e-12,45885.7,30.1963], Tmin=(100,'K'), Tmax=(1162.41,'K')), NASAPolynomial(coeffs=[14.0545,0.0328723,-1.40532e-05,2.65126e-09,-1.85723e-13,41596.6,-41.1196], Tmin=(1162.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_5) + radical(CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[O]C1[CH]C2=CC=C(C1)O2(6959)',
    structure = SMILES('[O]C1[CH]C2=CC=C(C1)O2'),
    E0 = (113.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64673,0.036747,3.24974e-05,-5.64467e-08,1.9269e-11,13694.1,22.0947], Tmin=(100,'K'), Tmax=(1167.54,'K')), NASAPolynomial(coeffs=[14.3213,0.0353389,-1.96729e-05,4.16511e-09,-3.09469e-13,7870.88,-53.2795], Tmin=(1167.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(Oxane) + ring(Furan) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]OC1C=C2C=C=C(C1)O2(6960)',
    structure = SMILES('[O]OC1C=C2C=C=C(C1)O2'),
    E0 = (467.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41289,0.0588603,-3.20512e-05,4.12986e-09,4.20681e-13,56305.1,23.3937], Tmin=(100,'K'), Tmax=(1766.57,'K')), NASAPolynomial(coeffs=[26.4557,0.0197799,-1.3832e-05,2.90149e-09,-2.04654e-13,44707.2,-119.457], Tmin=(1766.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(1,2-Cyclopentadiene) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=C2C[C]=C(C1)O2(6961)',
    structure = SMILES('[O]OC1C=C2C[C]=C(C1)O2'),
    E0 = (303.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,281.077,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92],'cm^-1')),
        HinderedRotor(inertia=(0.150393,'amu*angstrom^2'), symmetry=1, barrier=(3.45783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75833,0.0475387,-9.94608e-06,-9.88197e-09,3.59983e-12,36640.5,27.336], Tmin=(100,'K'), Tmax=(1520.9,'K')), NASAPolynomial(coeffs=[16.0273,0.0329647,-1.82105e-05,3.66376e-09,-2.57888e-13,29645.4,-56.2139], Tmin=(1520.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC1[C]=C2CC=C(C1)O2(6962)',
    structure = SMILES('[O]OC1[C]=C2CC=C(C1)O2'),
    E0 = (303.966,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,281.077,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,827.638,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92,1700.92],'cm^-1')),
        HinderedRotor(inertia=(0.150393,'amu*angstrom^2'), symmetry=1, barrier=(3.45783,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75833,0.0475387,-9.94608e-06,-9.88197e-09,3.59983e-12,36640.5,27.336], Tmin=(100,'K'), Tmax=(1520.9,'K')), NASAPolynomial(coeffs=[16.0273,0.0329647,-1.82105e-05,3.66376e-09,-2.57888e-13,29645.4,-56.2139], Tmin=(1520.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'OOC1[C]=C2[CH]C=C(C1)O2(6963)',
    structure = SMILES('OOC1[C]=C2[CH]C=C(C1)O2'),
    E0 = (268.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.102,0.0571075,-1.98687e-05,-6.78385e-09,3.50349e-12,32448.1,26.133], Tmin=(100,'K'), Tmax=(1447.11,'K')), NASAPolynomial(coeffs=[18.8867,0.0316588,-1.80668e-05,3.7084e-09,-2.65173e-13,24818.2,-74.8182], Tmin=(1447.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC1[CH]C2=CCC(=C1)O2(6964)',
    structure = SMILES('[O]OC1[CH]C2=CCC(=C1)O2'),
    E0 = (136.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68408,0.0335604,4.56871e-05,-7.1479e-08,2.48276e-11,16487.1,27.8904], Tmin=(100,'K'), Tmax=(1119.22,'K')), NASAPolynomial(coeffs=[14.1624,0.0353392,-1.88494e-05,3.98345e-09,-2.97815e-13,10789.2,-46.6919], Tmin=(1119.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[O]O[C]1C=C2CC=C(C1)O2(6965)',
    structure = SMILES('[O]OC1=C[C]2CC=C(C1)O2'),
    E0 = (242.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702285,0.0586627,-1.79755e-05,-1.56319e-08,8.08245e-12,29308.7,28.5486], Tmin=(100,'K'), Tmax=(1183.44,'K')), NASAPolynomial(coeffs=[18.6538,0.0269248,-1.44262e-05,3.03046e-09,-2.24685e-13,23033.3,-69.6424], Tmin=(1183.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(C2CsJOC(O))"""),
)

species(
    label = '[O]OC1C=C2OC3([CH]C23)C1(6966)',
    structure = SMILES('[O]OC1C=C2OC3([CH]C23)C1'),
    E0 = (534.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.160958,0.0692008,-2.90779e-05,-1.00338e-08,6.53818e-12,64442.3,23.267], Tmin=(100,'K'), Tmax=(1246.91,'K')), NASAPolynomial(coeffs=[23.1467,0.0270891,-1.64627e-05,3.56205e-09,-2.65934e-13,56251.6,-102.557], Tmin=(1246.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_0) + polycyclic(s1_3_6_ene_2) + polycyclic(s2_3_4_ane) - ring(Cyclopropane) - ring(3,4-Dihydro-2H-pyran) - ring(Oxetane) + radical(CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[CH]1C=C2CC3OOC3[C]1O2(6967)',
    structure = SMILES('C1=C[C]2O[C]1CC1OOC21'),
    E0 = (368.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43518,0.0532954,-1.07345e-05,-1.27558e-08,4.65334e-12,44454.3,21.593], Tmin=(100,'K'), Tmax=(1491.26,'K')), NASAPolynomial(coeffs=[18.5732,0.0352535,-2.06781e-05,4.24763e-09,-3.02451e-13,36237.5,-78.3487], Tmin=(1491.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_ene_5) + polycyclic(s2_4_6_ane) - ring(Oxane) + radical(C2CsJOCs) + radical(C2CsJOCs)"""),
)

species(
    label = 'OOC1C=C2C=C=C(C1)O2(6968)',
    structure = SMILES('OOC1C=C2C=C=C(C1)O2'),
    E0 = (315.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.934868,0.064827,-3.49496e-05,3.68466e-09,8.27253e-13,38044.2,24.3601], Tmin=(100,'K'), Tmax=(1648.4,'K')), NASAPolynomial(coeffs=[26.1098,0.0237777,-1.58317e-05,3.32794e-09,-2.37176e-13,27021.8,-117.936], Tmin=(1648.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(1,2-Cyclopentadiene)"""),
)

species(
    label = 'C1=C2CC3C=C(O2)C1O3(6969)',
    structure = SMILES('C1=C2CC3C=C(O2)C1O3'),
    E0 = (11.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16189,0.0194277,5.70222e-05,-6.65268e-08,1.8173e-11,1352.88,17.5273], Tmin=(100,'K'), Tmax=(1406.08,'K')), NASAPolynomial(coeffs=[17.7783,0.0402256,-3.17103e-05,7.09562e-09,-5.26873e-13,-8923.72,-79.8967], Tmin=(1406.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(25dihydrofuran) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) + polycyclic(s2_5_5_diene_0_5) - ring(3,4-Dihydro-2H-pyran) - ring(25dihydrofuran) - ring(2,3-Dihydrofuran)"""),
)

species(
    label = '[C]1=CCC2=CC=C1O2(5822)',
    structure = SMILES('[C]1=CCC2=CC=C1O2'),
    E0 = (240.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.114,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19905,0.0322078,1.49001e-05,-3.106e-08,9.86646e-12,28959.5,17.2785], Tmin=(100,'K'), Tmax=(1304.15,'K')), NASAPolynomial(coeffs=[12.9019,0.0300912,-1.79878e-05,3.8083e-09,-2.78965e-13,23556.3,-47.2116], Tmin=(1304.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan) + radical(C=CJC=C)"""),
)

species(
    label = 'C1=C2CC3C=C(O2)C1OO3(6970)',
    structure = SMILES('C1=C2CC3C=C(O2)C1OO3'),
    E0 = (406.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38359,0.0418849,1.89546e-05,-3.59722e-08,9.95563e-12,48961.8,17.5396], Tmin=(100,'K'), Tmax=(1519.78,'K')), NASAPolynomial(coeffs=[20.5725,0.0425262,-2.95605e-05,6.31335e-09,-4.55332e-13,37830.5,-96.2571], Tmin=(1519.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_6_6_ene_4) + polycyclic(s2_5_6_diene_0_6) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(2,3-Dihydrofuran) - ring(3,4-Dihydro-2H-pyran) - ring(36dihydro12dioxin) - ring(2,3-Dihydrofuran)"""),
)

species(
    label = '[O]OC1CC2=C[CH]C(=O)C21(6971)',
    structure = SMILES('[O]OC1C[C]2C=CC(=O)C21'),
    E0 = (199.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951179,0.0600098,-3.03767e-05,4.83702e-09,2.26779e-13,24110.5,24.4402], Tmin=(100,'K'), Tmax=(1546.83,'K')), NASAPolynomial(coeffs=[15.8388,0.0326573,-1.46607e-05,2.72197e-09,-1.84275e-13,18171.4,-58.188], Tmin=(1546.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(ROOJ) + radical(Allyl_T)"""),
)

species(
    label = 'C=C1OC2([CH]C12)C=CO[O](6972)',
    structure = SMILES('C=C1OC2([CH]C12)C=CO[O]'),
    E0 = (405.074,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,2750,3150,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194457,0.0568965,3.08857e-05,-1.00073e-07,4.76741e-11,48880.8,32.3688], Tmin=(100,'K'), Tmax=(934.32,'K')), NASAPolynomial(coeffs=[27.6731,0.00755048,4.63282e-07,-1.30653e-10,-1.5729e-15,40765.1,-114.286], Tmin=(934.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.074,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + polycyclic(s2_3_4_ane) + radical(CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=C2C=CC(=C1)O2(6973)',
    structure = SMILES('[O]OC1C=C2C=CC(=C1)O2'),
    E0 = (126.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381018,0.0697065,-4.28314e-05,6.08864e-09,9.51604e-13,15306.9,25.9586], Tmin=(100,'K'), Tmax=(1472.94,'K')), NASAPolynomial(coeffs=[26.3589,0.0205928,-1.46423e-05,3.20903e-09,-2.36387e-13,5329.09,-117.322], Tmin=(1472.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(1,4-Cyclohexadiene) + ring(25dihydrofuran) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1C=C2C=[C]C(C1)O2(6974)',
    structure = SMILES('[O]OC1C=C2C=[C]C(C1)O2'),
    E0 = (260.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376885,0.0602214,-3.62246e-06,-3.80523e-08,1.66686e-11,31452.4,26.7291], Tmin=(100,'K'), Tmax=(1124.04,'K')), NASAPolynomial(coeffs=[22.7289,0.0245936,-1.468e-05,3.26277e-09,-2.50712e-13,23653.4,-96.06], Tmin=(1124.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(25dihydrofuran) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]O[C]1C=C2C=CC(C1)O2(6975)',
    structure = SMILES('[O]OC1=C[C]2C=CC(C1)O2'),
    E0 = (283.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955066,0.0566069,-1.85409e-05,-8.7494e-09,4.7842e-12,34220.2,29.4217], Tmin=(100,'K'), Tmax=(1253.97,'K')), NASAPolynomial(coeffs=[14.6318,0.0337622,-1.60736e-05,3.15526e-09,-2.24085e-13,29156.2,-46.1712], Tmin=(1253.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s3_5_6_diene_1_5) + radical(ROOJ) + radical(C2CsJOCs)"""),
)

species(
    label = '[O]OC1C=C2[C]=CC(C1)O2(6976)',
    structure = SMILES('[O]OC1C=C2[C]=CC(C1)O2'),
    E0 = (221.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.307858,0.0624776,-9.16298e-06,-3.26285e-08,1.49934e-11,26782.2,26.0979], Tmin=(100,'K'), Tmax=(1123.08,'K')), NASAPolynomial(coeffs=[22.3134,0.0252764,-1.44703e-05,3.16676e-09,-2.41596e-13,19242.7,-94.1811], Tmin=(1123.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(25dihydrofuran) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = '[O]OC1[C]=C2C=CC(C1)O2(6977)',
    structure = SMILES('[O]OC1[C]=C2C=CC(C1)O2'),
    E0 = (260.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376885,0.0602214,-3.62246e-06,-3.80523e-08,1.66686e-11,31452.4,26.7291], Tmin=(100,'K'), Tmax=(1124.04,'K')), NASAPolynomial(coeffs=[22.7289,0.0245936,-1.468e-05,3.26277e-09,-2.50712e-13,23653.4,-96.06], Tmin=(1124.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(25dihydrofuran) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC1[CH]C23C=CC2(C1)O3(6978)',
    structure = SMILES('[O]OC1[CH]C23C=CC2(C1)O3'),
    E0 = (444.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33552,0.068519,-4.50947e-05,1.21446e-08,-1.23886e-12,53605.2,22.431], Tmin=(100,'K'), Tmax=(2245.5,'K')), NASAPolynomial(coeffs=[28.1061,0.0208316,-1.32394e-05,2.68712e-09,-1.85923e-13,41582.5,-128.377], Tmin=(2245.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(Cyclopentane) + ring(Cyclobutene) + polycyclic(s2_3_5_ane) + polycyclic(s2_3_4_ene_1) - ring(Cyclobutene) - ring(Ethylene_oxide) - ring(Cyclopentane) + radical(CCJCOOH) + radical(ROOJ)"""),
)

species(
    label = 'OOC1C=C2C=CC(=C1)O2(6979)',
    structure = SMILES('OOC1C=C2C=CC(=C1)O2'),
    E0 = (-25.8855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0325553,0.0743363,-4.1871e-05,1.66434e-09,2.67377e-12,-2960.08,26.4465], Tmin=(100,'K'), Tmax=(1377.13,'K')), NASAPolynomial(coeffs=[26.6816,0.0239679,-1.64573e-05,3.61775e-09,-2.68869e-13,-12863.6,-119.956], Tmin=(1377.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.8855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(1,4-Cyclohexadiene) + ring(25dihydrofuran)"""),
)

species(
    label = '[O]OC1=CC2=CC[C](C1)O2(6980)',
    structure = SMILES('[O]OC1=CC2=CC[C](C1)O2'),
    E0 = (224.997,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.240421,0.0712418,-5.10219e-05,1.54356e-08,-1.60858e-12,27205.8,25.9211], Tmin=(100,'K'), Tmax=(1506.14,'K')), NASAPolynomial(coeffs=[22.9055,0.0211497,-1.11945e-05,2.25974e-09,-1.60721e-13,19232.7,-96.5109], Tmin=(1506.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(C2CsJOC(O))"""),
)

species(
    label = 'C1=CC23CC(C=C1O2)O3(6981)',
    structure = SMILES('C1=CC23CC(C=C1O2)O3'),
    E0 = (-21.3092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648207,0.0523338,1.81629e-05,-5.77526e-08,2.10765e-11,-2424.99,14.9157], Tmin=(100,'K'), Tmax=(1212.16,'K')), NASAPolynomial(coeffs=[27.9346,0.0230385,-2.07565e-05,4.99493e-09,-3.91325e-13,-13503,-140.384], Tmin=(1212.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.3092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(1,3-Dioxane) - ring(Tetrahydrofuran) + ring(24dihydro13dioxin) + ring(25dihydrofuran) + polycyclic(s1_4_5_ene_1) + polycyclic(s3_4_6_ene_1) - ring(25dihydrofuran) - ring(Oxetane) - ring(24dihydro13dioxin)"""),
)

species(
    label = '[CH2]C12C=CC(=CC1O[O])O2(6982)',
    structure = SMILES('[CH2]C12C=CC(=CC1O[O])O2'),
    E0 = (494.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.435689,0.0799107,-4.04921e-05,-1.86475e-08,1.63467e-11,59626.1,27.4488], Tmin=(100,'K'), Tmax=(968.132,'K')), NASAPolynomial(coeffs=[25.044,0.0158563,-5.11172e-06,9.66619e-10,-7.44652e-14,52760.9,-104.627], Tmin=(968.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s3_5_5_diene_0_4) + radical(ROOJ) + radical(CJC(C)OC)"""),
)

species(
    label = 'C1=CC23CC(C=C1O2)OO3(6983)',
    structure = SMILES('C1=CC23CC(C=C1O2)OO3'),
    E0 = (-39.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456179,0.0558557,2.22356e-05,-6.25579e-08,2.25541e-11,-4548.68,15.8014], Tmin=(100,'K'), Tmax=(1205.24,'K')), NASAPolynomial(coeffs=[26.8856,0.0302867,-2.32865e-05,5.40426e-09,-4.17223e-13,-15433.1,-135.364], Tmin=(1205.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s3_5_6_ene_1) + polycyclic(s1_5_5_ene_1) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(3,4-Dihydro-2H-pyran) + ring(25dihydrofuran) - ring(25dihydrofuran) - ring(3,4-Dihydro-2H-pyran) - ring(12dioxolane)"""),
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
    E0 = (106.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (310.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (326.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (297.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (447.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (335.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (212.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (215.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (379.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (393.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (301.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (557.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (480.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (189.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (303.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (310.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (210.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (274.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (184.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (126.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (268.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (235.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (365.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (245.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (114.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (334.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (420.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (356.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (695.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (416.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (453.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (370.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (328.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (358.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (536.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (369.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (317.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (255.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (284.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (406.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (420.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (552.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (338.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (442.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (353.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (359.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (304.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (448.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (140.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (341.552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (227.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (651.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (113.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)(2)', 'C7H6O(490)(489)'],
    products = ['S(554)(553)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(187.105,'m^3/(mol*s)'), n=1.3603, Ea=(73.6222,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Cds-CsH_Cds-(Cd-Cd-Cd)H;YJ] for rate rule [Cds-CsH_Cds-(Cd-Cd-Cd)H;O2b]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.1 to 73.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['S(554)(553)'],
    products = ['[CH]1C=C2CC3[CH]C1(OO3)O2(6936)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.36e+10,'s^-1'), n=0.44, Ea=(204.344,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;doublebond_intra_HDe_secNd;radadd_intra] for rate rule [R6;doublebond_intra_HDe_secNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)(3)', '[O]OC1=CC2=CC=C(C1)O2(6937)'],
    products = ['S(554)(553)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]O[C]1CC2=CC=C(C1)O2(6938)'],
    products = ['S(554)(553)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(40217,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OO[C]1[CH]C2=CC=C(C1)O2(6939)'],
    products = ['S(554)(553)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC1CC2=[C]C=C(C1)O2(6940)'],
    products = ['S(554)(553)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.47715e+10,'s^-1'), n=0.8, Ea=(147.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_H/NonDeC] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['S(554)(553)'],
    products = ['OOC1[CH]C2=CC=C([CH]1)O2(6941)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.56e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OOC1[CH]C2=C[C]=C(C1)O2(6942)'],
    products = ['S(554)(553)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(13085.7,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 4.58257569496
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OOC1[CH]C2=[C]C=C(C1)O2(6943)'],
    products = ['S(554)(553)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O2(2)(2)', '[CH]1C=C2C=C[C](C1)O2(5783)'],
    products = ['S(554)(553)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.08764e+06,'m^3/(mol*s)'), n=0.19382, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1CC2=CC3[C](O2)C13(6944)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(195.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_short;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c7_beta_short;doublebond_intra_secNd_HCd;radadd_intra_csHCs]
Euclidian distance = 3.60555127546
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1C[C]2OC3=CC2C31(6945)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.16118e+10,'s^-1'), n=-0.279759, Ea=(451.048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_cyclic;doublebond_intra_secNd;radadd_intra_cs] for rate rule [Rn1c5_beta_long;doublebond_intra_secNd_HCd;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1CC23[CH]C=C(O2)C13(6946)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.33004e+14,'s^-1'), n=-0.792922, Ea=(374.502,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_long_SSS_D;doublebond_intra_pri_NdNd;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['S(554)(553)'],
    products = ['[CH]1C2=CC3OOC1C[C]3O2(6947)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.0232e+08,'s^-1'), n=0.665438, Ea=(83.5089,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;doublebond_intra_secNd;radadd_intra] + [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_secNd_HCd;radadd_intra_O]
Euclidian distance = 3.74165738677
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['S(554)(553)'],
    products = ['[CH]1[C]2OC3=CC2OOC1C3(6948)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.0232e+08,'s^-1'), n=0.665438, Ea=(197.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6;doublebond_intra_secNd;radadd_intra] + [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn2c7_gamma_short;doublebond_intra_secNd_HCd;radadd_intra_O]
Euclidian distance = 3.74165738677
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['S(554)(553)'],
    products = ['[CH]1C2=C[CH]C3(CC1OO3)O2(6949)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(64782.7,'s^-1'), n=1.23, Ea=(204.344,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_NdNd;radadd_intra] for rate rule [Rn2c7_beta_long;doublebond_intra_pri_NdNd;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C1=CC=C(C=CO[O])O1(6950)'],
    products = ['S(554)(553)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.26e+06,'s^-1'), n=1.02, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra_pri_HNd;radadd_intra_cs2H] for rate rule [R6_SSS_D;doublebond_intra_pri_HNd_O;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OC1C=C=CC=C([O])C1(6951)'],
    products = ['S(554)(553)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.10864e+11,'s^-1'), n=0.310877, Ea=(56.2836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c7_gamma_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['S(554)(553)'],
    products = ['OOC1=CC2=CC=C(C1)O2(6952)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['S(554)(553)'],
    products = ['O2(S)(162)(161)', 'C7H6O(490)(489)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(20.6843,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 20.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC1[CH]C2C=CC(=C1)O2(6953)'],
    products = ['S(554)(553)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.125e+09,'s^-1'), n=0.991, Ea=(45.9529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_pentadiene;CH_end;CdHC_2] for rate rule [1_3_pentadiene;CH(CJ)_1;CdHC_2]
Euclidian distance = 1.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction22',
    reactants = ['S(554)(553)'],
    products = ['O(4)(4)', 'C1C=C2OC=1CC1OC21(6954)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2OO_S;C_sec_rad_intra;OOJ] for rate rule [R2OO_S;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]O[CH]C1CC2=CC=C1O2(6955)'],
    products = ['S(554)(553)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction24',
    reactants = ['S(554)(553)'],
    products = ['HO2(8)(9)', 'C1=C[C]2C=CC(=C1)O2(5776)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO] for rate rule [R2OO_HDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction25',
    reactants = ['S(554)(553)'],
    products = ['C1C=C2OC=1CC1OOC21(6956)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1[CH]C2=CC2C(=O)C1(6957)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5460.43,'s^-1'), n=2.73, Ea=(228.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.2360679775
family: ketoenol"""),
)

reaction(
    label = 'reaction27',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1[CH]C(=O)C2C=C2C1(6958)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)(4)', '[O]C1[CH]C2=CC=C(C1)O2(6959)'],
    products = ['S(554)(553)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)(3)', '[O]OC1C=C2C=C=C(C1)O2(6960)'],
    products = ['S(554)(553)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.149e+09,'cm^3/(mol*s)'), n=1.595, Ea=(16.7946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-OneDeH;HJ] for rate rule [Ca_Cds-CdH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC1C=C2C[C]=C(C1)O2(6961)'],
    products = ['S(554)(553)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.3e+10,'s^-1'), n=0.98, Ea=(112.55,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OC1[C]=C2CC=C(C1)O2(6962)'],
    products = ['S(554)(553)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.182e+10,'s^-1'), n=0.86, Ea=(149.369,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['OOC1[C]=C2[CH]C=C(C1)O2(6963)'],
    products = ['S(554)(553)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]OC1[CH]C2=CCC(=C1)O2(6964)'],
    products = ['S(554)(553)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(285077,'s^-1'), n=1.88833, Ea=(191.836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] + [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_1H] for rate rule [R4H_SDS;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]O[C]1C=C2CC=C(C1)O2(6965)'],
    products = ['S(554)(553)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.30755e+06,'s^-1'), n=1.805, Ea=(115.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;C_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_SDS;C_rad_out_NDMustO;Cs_H_out_H/Cd]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1C=C2OC3([CH]C23)C1(6966)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.22857e+12,'s^-1'), n=0, Ea=(430.603,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri;radadd_intra_csHCd] for rate rule [Rn0c5_beta_short;doublebond_intra_pri_NdNd;radadd_intra_csHCd]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction36',
    reactants = ['S(554)(553)'],
    products = ['[CH]1C=C2CC3OOC3[C]1O2(6967)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.44296e+10,'s^-1'), n=0.26145, Ea=(262.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra;radadd_intra] for rate rule [Rn2c6_alpha_long;doublebond_intra;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['S(554)(553)'],
    products = ['OOC1C=C2C=C=C(C1)O2(6968)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(211.191,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad_De] for rate rule [R6radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['S(554)(553)'],
    products = ['O(4)(4)', 'C1=C2CC3C=C(O2)C1O3(6969)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.63e+10,'s^-1'), n=0, Ea=(149.742,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;C_sec_rad_intra;OO] for rate rule [R4OO_SDS;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.73205080757
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['S(554)(553)'],
    products = ['HO2(8)(9)', '[C]1=CCC2=CC=C1O2(5822)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction40',
    reactants = ['S(554)(553)'],
    products = ['C1=C2CC3C=C(O2)C1OO3(6970)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(300.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination
Ea raised from 300.4 to 300.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1CC2=C[CH]C(=O)C21(6971)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond;R_O_C] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond;R_O_C]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C1OC2([CH]C12)C=CO[O](6972)'],
    products = ['S(554)(553)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6.21184e+10,'s^-1'), n=0.288169, Ea=(147.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_5_unsaturated_hexane] for rate rule [1_5_hexadiene]
Euclidian distance = 1.0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)(3)', '[O]OC1C=C2C=CC(=C1)O2(6973)'],
    products = ['S(554)(553)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(50.6821,'m^3/(mol*s)'), n=1.79231, Ea=(0.501101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDe;HJ] for rate rule [Cds-CsH_Cds-CdOs;HJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[O]OC1C=C2C=[C]C(C1)O2(6974)'],
    products = ['S(554)(553)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.677e+10,'s^-1'), n=0.839, Ea=(182.581,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out_noH] for rate rule [R2H_S_cy5;Cd_rad_out_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]O[C]1C=C2C=CC(C1)O2(6975)'],
    products = ['S(554)(553)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(8.59e-19,'s^-1'), n=8.79, Ea=(69.7054,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;Cs_H_out_OneDe] for rate rule [R3H_SS_Cs;C_rad_out_OneDe/O;Cs_H_out_(CdCdCd)]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O]OC1C=C2[C]=CC(C1)O2(6976)'],
    products = ['S(554)(553)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.20504e+10,'s^-1'), n=0.725, Ea=(137.863,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cs_H_out_NonDe] for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;Cs_H_out_NDMustO]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[O]OC1[C]=C2C=CC(C1)O2(6977)'],
    products = ['S(554)(553)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;Cs_H_out_noH]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['S(554)(553)'],
    products = ['[O]OC1[CH]C23C=CC2(C1)O3(6978)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.81521e+10,'s^-1'), n=0.581539, Ea=(342.315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c6_beta_long_SS_D;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn0c6_beta_long_SS_D;doublebond_intra_pri_NdCd;radadd_intra_csNdCd]
Euclidian distance = 2.82842712475
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['S(554)(553)'],
    products = ['OOC1C=C2C=CC(=C1)O2(6979)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O]OC1=CC2=CC[C](C1)O2(6980)'],
    products = ['S(554)(553)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;unsaturated_end]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction51',
    reactants = ['S(554)(553)'],
    products = ['O(4)(4)', 'C1=CC23CC(C=C1O2)O3(6981)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(3.31e+11,'s^-1'), n=0, Ea=(121.556,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_ter_rad_intra;OO] for rate rule [R3OO_SS;C_rad/ODMustO_intra;OOJ]
Euclidian distance = 2.2360679775
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C12C=CC(=CC1O[O])O2(6982)'],
    products = ['S(554)(553)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction53',
    reactants = ['S(554)(553)'],
    products = ['C1=CC23CC(C=C1O2)OO3(6983)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_OneDe/O]
Euclidian distance = 3.16227766017
family: Birad_recombination"""),
)

network(
    label = '404',
    isomers = [
        'S(554)(553)',
    ],
    reactants = [
        ('O2(2)(2)', 'C7H6O(490)(489)'),
        ('O2(S)(162)(161)', 'C7H6O(490)(489)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '404',
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

