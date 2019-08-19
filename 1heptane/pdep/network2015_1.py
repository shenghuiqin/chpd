species(
    label = 'C=[C]OO[CH]CC(10712)',
    structure = SMILES('C=[C]OO[CH]CC'),
    E0 = (322.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,348.669,2188.85],'cm^-1')),
        HinderedRotor(inertia=(0.106986,'amu*angstrom^2'), symmetry=1, barrier=(9.22951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106989,'amu*angstrom^2'), symmetry=1, barrier=(9.22946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498419,'amu*angstrom^2'), symmetry=1, barrier=(42.9967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106978,'amu*angstrom^2'), symmetry=1, barrier=(9.22943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498416,'amu*angstrom^2'), symmetry=1, barrier=(42.9966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13479,0.0685823,-7.08299e-05,4.64284e-08,-1.35437e-11,38902,29.32], Tmin=(100,'K'), Tmax=(803.018,'K')), NASAPolynomial(coeffs=[6.65463,0.0410863,-1.94677e-05,3.78671e-09,-2.68044e-13,38015.5,3.90098], Tmin=(803.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CCsJOOC)"""),
)

species(
    label = 'C2H5CHO(70)',
    structure = SMILES('CCC=O'),
    E0 = (-204.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.207559,'amu*angstrom^2'), symmetry=1, barrier=(4.77219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208362,'amu*angstrom^2'), symmetry=1, barrier=(4.79065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90578,0.0240644,-7.06356e-06,-9.81837e-10,5.55825e-13,-24535.9,13.5806], Tmin=(100,'K'), Tmax=(1712.49,'K')), NASAPolynomial(coeffs=[7.69109,0.0189242,-7.84934e-06,1.38273e-09,-8.99057e-14,-27060.1,-14.6647], Tmin=(1712.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""propanal""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C]OOC=CC(11022)',
    structure = SMILES('C=[C]OOC=CC'),
    E0 = (266.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.253915,'amu*angstrom^2'), symmetry=1, barrier=(5.838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728204,'amu*angstrom^2'), symmetry=1, barrier=(16.7428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663946,'amu*angstrom^2'), symmetry=1, barrier=(15.2654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727973,'amu*angstrom^2'), symmetry=1, barrier=(16.7375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25743,0.0595573,-4.84123e-05,2.07985e-08,-3.71074e-12,32195.5,28.8068], Tmin=(100,'K'), Tmax=(1301.38,'K')), NASAPolynomial(coeffs=[11.4378,0.0282661,-1.23452e-05,2.32218e-09,-1.61352e-13,29545.8,-22.9894], Tmin=(1301.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C#COO[CH]CC(11023)',
    structure = SMILES('C#COO[CH]CC'),
    E0 = (314.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,180],'cm^-1')),
        HinderedRotor(inertia=(0.384103,'amu*angstrom^2'), symmetry=1, barrier=(44.1447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0965713,'amu*angstrom^2'), symmetry=1, barrier=(10.8746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.928663,'amu*angstrom^2'), symmetry=1, barrier=(21.3518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91988,'amu*angstrom^2'), symmetry=1, barrier=(44.1417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92026,'amu*angstrom^2'), symmetry=1, barrier=(44.1506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925194,0.0730536,-7.94331e-05,4.96646e-08,-1.31722e-11,37920.3,24.5518], Tmin=(100,'K'), Tmax=(894.501,'K')), NASAPolynomial(coeffs=[9.20233,0.036039,-1.7361e-05,3.40116e-09,-2.41808e-13,36439.5,-14.4575], Tmin=(894.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCt) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOOC)"""),
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
    label = 'C=[C]OOC=C(3202)',
    structure = SMILES('C=[C]OOC=C'),
    E0 = (302.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(0.0802844,'amu*angstrom^2'), symmetry=1, barrier=(10.1627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.939515,'amu*angstrom^2'), symmetry=1, barrier=(21.6013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.93945,'amu*angstrom^2'), symmetry=1, barrier=(21.5998,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61701,0.0472731,-3.89653e-05,1.6305e-08,-2.73503e-12,36519.2,25.7161], Tmin=(100,'K'), Tmax=(1419.4,'K')), NASAPolynomial(coeffs=[12.3675,0.0169774,-6.9495e-06,1.26785e-09,-8.65498e-14,33467.3,-29.9141], Tmin=(1419.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OOC[CH]C(11024)',
    structure = SMILES('C=[C]OOC[CH]C'),
    E0 = (336.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09876,0.0694263,-7.50842e-05,5.18762e-08,-1.58454e-11,40583.9,31.4572], Tmin=(100,'K'), Tmax=(771.968,'K')), NASAPolynomial(coeffs=[6.63119,0.0407595,-1.9382e-05,3.77184e-09,-2.66846e-13,39729.7,6.19829], Tmin=(771.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]=COO[CH]CC(11025)',
    structure = SMILES('[CH]=COO[CH]CC'),
    E0 = (329.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,964.802],'cm^-1')),
        HinderedRotor(inertia=(0.984018,'amu*angstrom^2'), symmetry=1, barrier=(22.6245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.535671,'amu*angstrom^2'), symmetry=1, barrier=(12.3161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334887,'amu*angstrom^2'), symmetry=1, barrier=(7.69972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983948,'amu*angstrom^2'), symmetry=1, barrier=(22.6229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.6712,'amu*angstrom^2'), symmetry=1, barrier=(61.4161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793626,0.0712625,-6.47897e-05,3.15428e-08,-6.35782e-12,39802.8,27.8021], Tmin=(100,'K'), Tmax=(1168.49,'K')), NASAPolynomial(coeffs=[12.2903,0.0319064,-1.42674e-05,2.71762e-09,-1.9058e-13,37116,-29.4529], Tmin=(1168.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]CCOO[C]=C(11026)',
    structure = SMILES('[CH2]CCOO[C]=C'),
    E0 = (341.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19622,0.0670929,-6.58068e-05,3.9652e-08,-1.05952e-11,41161.8,30.5913], Tmin=(100,'K'), Tmax=(870.528,'K')), NASAPolynomial(coeffs=[7.13189,0.0398191,-1.88118e-05,3.66249e-09,-2.59666e-13,40128.4,2.77809], Tmin=(870.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=COO[CH][CH]C(11027)',
    structure = SMILES('C=COO[CH][CH]C'),
    E0 = (283.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.95998,0.0689547,-6.11331e-05,2.94715e-08,-5.97369e-12,34181.2,29.1082], Tmin=(100,'K'), Tmax=(1149.91,'K')), NASAPolynomial(coeffs=[11.0424,0.0338828,-1.53836e-05,2.94817e-09,-2.07333e-13,31862.4,-20.9423], Tmin=(1149.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCJCOOH) + radical(CCsJOOC)"""),
)

species(
    label = '[CH]=[C]OOCCC(11028)',
    structure = SMILES('[CH]=[C]OOCCC'),
    E0 = (383.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,330.166],'cm^-1')),
        HinderedRotor(inertia=(0.0996537,'amu*angstrom^2'), symmetry=1, barrier=(7.69195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995594,'amu*angstrom^2'), symmetry=1, barrier=(7.669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0999757,'amu*angstrom^2'), symmetry=1, barrier=(7.66901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148315,'amu*angstrom^2'), symmetry=1, barrier=(11.5152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.58783,'amu*angstrom^2'), symmetry=1, barrier=(44.7287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04949,0.0703131,-7.35662e-05,4.68707e-08,-1.30038e-11,46200.4,29.7338], Tmin=(100,'K'), Tmax=(848.681,'K')), NASAPolynomial(coeffs=[7.66446,0.0391362,-1.84639e-05,3.58709e-09,-2.53825e-13,45077.6,-1.09441], Tmin=(848.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C[CH]OOC=C(11029)',
    structure = SMILES('[CH2]C[CH]OOC=C'),
    E0 = (288.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89717,0.0685853,-5.90524e-05,2.70691e-08,-5.15811e-12,34765.9,28.8117], Tmin=(100,'K'), Tmax=(1225.49,'K')), NASAPolynomial(coeffs=[12.1036,0.0320082,-1.4283e-05,2.71498e-09,-1.89977e-13,32019.2,-27.5318], Tmin=(1225.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCsJOOC)"""),
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
    label = 'CC[CH][O](563)',
    structure = SMILES('CC[CH][O]'),
    E0 = (133.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,298.357,1774.23],'cm^-1')),
        HinderedRotor(inertia=(0.129074,'amu*angstrom^2'), symmetry=1, barrier=(8.14273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00364816,'amu*angstrom^2'), symmetry=1, barrier=(8.14268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1585,0.0245341,-8.42945e-06,1.83944e-10,2.32791e-13,16036.2,14.3859], Tmin=(100,'K'), Tmax=(2077.96,'K')), NASAPolynomial(coeffs=[11.8474,0.0146996,-6.30487e-06,1.09829e-09,-6.9226e-14,10937.4,-37.4679], Tmin=(2077.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C=COOC=CC(11030)',
    structure = SMILES('C=COOC=CC'),
    E0 = (27.1198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729695,0.0615293,-3.45834e-05,-8.64036e-10,5.05414e-12,3388.52,27.5742], Tmin=(100,'K'), Tmax=(1035.29,'K')), NASAPolynomial(coeffs=[15.6807,0.0239503,-9.38402e-06,1.74281e-09,-1.23303e-13,-789.043,-50.2992], Tmin=(1035.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.1198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
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
    label = 'C=[C]OO[CH]C(1085)',
    structure = SMILES('C=[C]OO[CH]C'),
    E0 = (346.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,260.358],'cm^-1')),
        HinderedRotor(inertia=(0.205867,'amu*angstrom^2'), symmetry=1, barrier=(9.97457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208662,'amu*angstrom^2'), symmetry=1, barrier=(9.98598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0142427,'amu*angstrom^2'), symmetry=1, barrier=(43.577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.905496,'amu*angstrom^2'), symmetry=1, barrier=(43.5734,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0892,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3210.06,'J/mol'), sigma=(5.6222,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.40 K, Pc=40.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71904,0.0546222,-6.02795e-05,4.32493e-08,-1.37781e-11,41741.9,24.9713], Tmin=(100,'K'), Tmax=(741.345,'K')), NASAPolynomial(coeffs=[5.81455,0.0325225,-1.55603e-05,3.03138e-09,-2.14437e-13,41134.7,6.43894], Tmin=(741.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CCsJOOC)"""),
)

species(
    label = '[CH2]C(C)OO[C]=C(11031)',
    structure = SMILES('[CH2]C(C)OO[C]=C'),
    E0 = (337.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,500,795,815,3000,3100,440,815,1455,1000,299.151],'cm^-1')),
        HinderedRotor(inertia=(0.315478,'amu*angstrom^2'), symmetry=1, barrier=(19.526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140708,'amu*angstrom^2'), symmetry=1, barrier=(8.50752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.089826,'amu*angstrom^2'), symmetry=1, barrier=(5.39411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00191311,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.588193,'amu*angstrom^2'), symmetry=1, barrier=(36.1214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06518,0.0674024,-6.01604e-05,2.96161e-08,-6.1795e-12,40696.5,30.9114], Tmin=(100,'K'), Tmax=(1114.42,'K')), NASAPolynomial(coeffs=[10.2675,0.0343716,-1.57003e-05,3.01865e-09,-2.1271e-13,38645.5,-14.4813], Tmin=(1114.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CJCOOH)"""),
)

species(
    label = 'C=C1OOC1CC(3475)',
    structure = SMILES('C=C1OOC1CC'),
    E0 = (24.769,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58461,0.0343002,4.55884e-05,-8.17567e-08,3.25374e-11,3082.66,24.1296], Tmin=(100,'K'), Tmax=(995.949,'K')), NASAPolynomial(coeffs=[15.3038,0.0252698,-1.01961e-05,2.02922e-09,-1.5282e-13,-1934.9,-53.4724], Tmin=(995.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.769,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C=[C]O[O](110)',
    structure = SMILES('C=[C]O[O]'),
    E0 = (339.001,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,492.5,1135,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.0942477,'amu*angstrom^2'), symmetry=1, barrier=(2.16694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10787,0.026676,-5.52143e-05,6.05409e-08,-2.36332e-11,40797.8,12.0288], Tmin=(100,'K'), Tmax=(872.797,'K')), NASAPolynomial(coeffs=[1.34334,0.0171001,-8.40162e-06,1.59781e-09,-1.08413e-13,41778.6,24.1556], Tmin=(872.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]CC(144)',
    structure = SMILES('[CH]CC'),
    E0 = (328.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,328.86,330.153,2894.73],'cm^-1')),
        HinderedRotor(inertia=(0.00155435,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00154734,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00371,0.0161454,1.44268e-05,-2.62511e-08,1.01943e-11,39516.8,12.6846], Tmin=(100,'K'), Tmax=(1003.27,'K')), NASAPolynomial(coeffs=[6.63826,0.0156478,-5.75067e-06,1.05874e-09,-7.51289e-14,38083.3,-8.37163], Tmin=(1003.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC[CH]O[O](238)',
    structure = SMILES('CC[CH]O[O]'),
    E0 = (139.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,208.348],'cm^-1')),
        HinderedRotor(inertia=(0.243143,'amu*angstrom^2'), symmetry=1, barrier=(7.00807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245015,'amu*angstrom^2'), symmetry=1, barrier=(7.01593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00088615,'amu*angstrom^2'), symmetry=1, barrier=(7.02015,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0785,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20785,0.0420296,-4.3379e-05,3.01324e-08,-9.24759e-12,16864.4,19.4406], Tmin=(100,'K'), Tmax=(771.405,'K')), NASAPolynomial(coeffs=[5.43298,0.0253065,-1.08618e-05,2.03087e-09,-1.40556e-13,16366.8,4.71818], Tmin=(771.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCsJOOH) + radical(ROOJ)"""),
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
    E0 = (322.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (485.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (541.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (467.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (448.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (435.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (465.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (475.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (528.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (465.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (322.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (361.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (766.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (513.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (330.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (667.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (540.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['C2H5CHO(70)', 'CH2CO(28)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=[C]OOC=CC(11022)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.72e+08,'cm^3/(mol*s)'), n=1.477, Ea=(6.73624,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2825 used for Cds-CsH_Cds-OsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C#COO[CH]CC(11023)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3(17)', 'C=[C]OOC=C(3202)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(41670,'cm^3/(mol*s)'), n=2.299, Ea=(28.5767,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2827 used for Cds-HH_Cds-OsH;CsJ-HHH
Exact match found for rate rule [Cds-HH_Cds-OsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['C=[C]OOC[CH]C(11024)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.4e-19,'s^-1'), n=8.84, Ea=(125.52,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 342 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=COO[CH]CC(11025)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['[CH2]CCOO[C]=C(11026)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 343 used for R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['C=COO[CH][CH]C(11027)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.54505e+10,'s^-1'), n=0.959062, Ea=(152.545,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]OOCCC(11028)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['[CH2]C[CH]OOC=C(11029)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6HJ_3;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C][O](173)', 'CC[CH][O](563)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(29.3206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 29.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['C=COOC=CC(11030)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2(S)(23)', 'C=[C]OO[CH]C(1085)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C)OO[C]=C(11031)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.62395e+08,'s^-1'), n=1.1925, Ea=(176.042,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;CH3] + [cCs(-HR!H)CJ;CsJ;CH3] for rate rule [cCs(-HR!H)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OO[CH]CC(10712)'],
    products = ['C=C1OOC1CC(3475)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]O[O](110)', '[CH]CC(144)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC[CH]O[O](238)', 'H2CC(41)'],
    products = ['C=[C]OO[CH]CC(10712)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

network(
    label = '2015',
    isomers = [
        'C=[C]OO[CH]CC(10712)',
    ],
    reactants = [
        ('C2H5CHO(70)', 'CH2CO(28)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2015',
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

