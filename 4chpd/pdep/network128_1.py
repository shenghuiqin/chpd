species(
    label = 'C=CC(C=C)O[O](423)',
    structure = SMILES('C=CC(C=C)O[O]'),
    E0 = (134.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,609.234],'cm^-1')),
        HinderedRotor(inertia=(0.58829,'amu*angstrom^2'), symmetry=1, barrier=(13.526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.588408,'amu*angstrom^2'), symmetry=1, barrier=(13.5287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.362449,'amu*angstrom^2'), symmetry=1, barrier=(13.5293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17829,0.0576802,-4.33629e-05,1.66179e-08,-2.59638e-12,16272.7,26.8607], Tmin=(100,'K'), Tmax=(1488.12,'K')), NASAPolynomial(coeffs=[13.253,0.0252236,-1.06468e-05,1.96116e-09,-1.34064e-13,12679,-36.1923], Tmin=(1488.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = '[CH2]C=CC=C(162)',
    structure = SMILES('[CH2]C=CC=C'),
    E0 = (175.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.52841,'amu*angstrom^2'), symmetry=1, barrier=(35.1412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52752,'amu*angstrom^2'), symmetry=1, barrier=(35.1208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3140.68,'J/mol'), sigma=(5.4037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=490.57 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2987,0.0238036,3.93739e-05,-6.93177e-08,2.88419e-11,21235.8,16.7723], Tmin=(100,'K'), Tmax=(942.053,'K')), NASAPolynomial(coeffs=[11.9106,0.0173605,-5.0926e-06,8.77895e-10,-6.40305e-14,17899.7,-37.1203], Tmin=(942.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C1OOC1C=C(657)',
    structure = SMILES('[CH2]C1OOC1C=C'),
    E0 = (229.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2191,0.0511756,-1.9157e-05,-1.00021e-08,6.77626e-12,27659.5,23.7932], Tmin=(100,'K'), Tmax=(1064.13,'K')), NASAPolynomial(coeffs=[13.4078,0.025146,-1.03577e-05,1.95943e-09,-1.39229e-13,23945.1,-41.032], Tmin=(1064.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CJCOOH)"""),
)

species(
    label = 'C=C[C](C=C)OO(658)',
    structure = SMILES('[CH2]C=C(C=C)OO'),
    E0 = (86.5533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.820249,0.0634873,-5.3216e-05,2.33015e-08,-4.11446e-12,10529.8,25.3688], Tmin=(100,'K'), Tmax=(1351.85,'K')), NASAPolynomial(coeffs=[14.1885,0.0239315,-9.3246e-06,1.6561e-09,-1.11484e-13,6915.44,-43.1551], Tmin=(1351.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.5533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C(C=C)OO(659)',
    structure = SMILES('C=[C]C(C=C)OO'),
    E0 = (220.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16713,0.0625803,-5.23025e-05,2.263e-08,-4.04944e-12,26593.1,27.3041], Tmin=(100,'K'), Tmax=(1297.33,'K')), NASAPolynomial(coeffs=[12.1311,0.028776,-1.32178e-05,2.54554e-09,-1.79165e-13,23748.3,-28.4451], Tmin=(1297.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C=C)OO(660)',
    structure = SMILES('[CH]=CC(C=C)OO'),
    E0 = (229.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885317,0.0642782,-5.35503e-05,2.24348e-08,-3.78936e-12,27720,28.2489], Tmin=(100,'K'), Tmax=(1397.88,'K')), NASAPolynomial(coeffs=[14.8117,0.0244282,-1.07892e-05,2.04154e-09,-1.42186e-13,23826.5,-43.6028], Tmin=(1397.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C2H3(44)',
    structure = SMILES('[CH]=C'),
    E0 = (286.361,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,677.08,1086.68,3788.01],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(286.361,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C[CH]O[O](661)',
    structure = SMILES('C=C[CH]O[O]'),
    E0 = (199.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,492.5,1135,1000,3025,407.5,1350,352.5,316.51],'cm^-1')),
        HinderedRotor(inertia=(0.157895,'amu*angstrom^2'), symmetry=1, barrier=(11.1941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157042,'amu*angstrom^2'), symmetry=1, barrier=(11.2,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21608,0.0376502,-3.45386e-05,1.68249e-08,-3.30513e-12,24056.5,18.0455], Tmin=(100,'K'), Tmax=(1220.77,'K')), NASAPolynomial(coeffs=[9.39745,0.0141196,-5.62587e-06,1.03554e-09,-7.16583e-14,22303.1,-18.0331], Tmin=(1220.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO)"""),
)

species(
    label = 'H(19)',
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
    label = 'C=C[C](C=C)O[O](662)',
    structure = SMILES('[CH2]C=C(C=C)O[O]'),
    E0 = (238.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,492.5,1135,1000,350,440,435,1725,268.262],'cm^-1')),
        HinderedRotor(inertia=(0.181066,'amu*angstrom^2'), symmetry=1, barrier=(9.35282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183789,'amu*angstrom^2'), symmetry=1, barrier=(9.34032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39122,'amu*angstrom^2'), symmetry=1, barrier=(71.1944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24575,0.0580629,-5.18992e-05,2.54273e-08,-5.1033e-12,28793.2,24.5966], Tmin=(100,'K'), Tmax=(1188.71,'K')), NASAPolynomial(coeffs=[11.2048,0.0245509,-9.61146e-06,1.71104e-09,-1.15501e-13,26425.5,-25.1717], Tmin=(1188.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C(C=C)O[O](663)',
    structure = SMILES('C=[C]C(C=C)O[O]'),
    E0 = (372.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,1685,370,492.5,1135,1000,3010,987.5,1337.5,450,1655,294.544,294.58],'cm^-1')),
        HinderedRotor(inertia=(0.194728,'amu*angstrom^2'), symmetry=1, barrier=(11.9795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194506,'amu*angstrom^2'), symmetry=1, barrier=(11.979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194532,'amu*angstrom^2'), symmetry=1, barrier=(11.9788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48711,0.0582291,-5.4043e-05,2.7923e-08,-6.11702e-12,44861.5,26.9232], Tmin=(100,'K'), Tmax=(1066.32,'K')), NASAPolynomial(coeffs=[9.15563,0.0294626,-1.35768e-05,2.6233e-09,-1.85434e-13,43226.1,-10.5654], Tmin=(1066.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C=C)O[O](664)',
    structure = SMILES('[CH]=CC(C=C)O[O]'),
    E0 = (381.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(0.602836,'amu*angstrom^2'), symmetry=1, barrier=(13.8604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603153,'amu*angstrom^2'), symmetry=1, barrier=(13.8677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603051,'amu*angstrom^2'), symmetry=1, barrier=(13.8653,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31348,0.0588054,-5.20174e-05,2.42616e-08,-4.65637e-12,45983.4,27.4687], Tmin=(100,'K'), Tmax=(1226.44,'K')), NASAPolynomial(coeffs=[11.5055,0.0255643,-1.13616e-05,2.16191e-09,-1.51506e-13,43483.5,-23.7822], Tmin=(1226.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1[CH]COO1(665)',
    structure = SMILES('C=CC1[CH]COO1'),
    E0 = (143.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6802,0.0328582,4.92418e-05,-9.68077e-08,4.32942e-11,17328,24.4193], Tmin=(100,'K'), Tmax=(901.439,'K')), NASAPolynomial(coeffs=[16.5105,0.0152161,-1.5489e-06,2.85594e-11,-4.84846e-16,12697.3,-56.4441], Tmin=(901.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(CCJCOOH)"""),
)

species(
    label = 'HO2(26)',
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
    label = 'C=C=CC=C(320)',
    structure = SMILES('C=C=CC=C'),
    E0 = (234.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.27013,'amu*angstrom^2'), symmetry=1, barrier=(29.2028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06131,0.0323124,7.18216e-06,-3.50282e-08,1.65081e-11,28283,15.8131], Tmin=(100,'K'), Tmax=(957.431,'K')), NASAPolynomial(coeffs=[12.3474,0.0146244,-4.72098e-06,8.44162e-10,-6.13759e-14,25154.4,-39.4163], Tmin=(957.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'O(20)',
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
    label = 'C=CC([O])C=C(666)',
    structure = SMILES('C=CC([O])C=C'),
    E0 = (141.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,425.747,427.177,427.601],'cm^-1')),
        HinderedRotor(inertia=(0.118091,'amu*angstrom^2'), symmetry=1, barrier=(15.3066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11825,'amu*angstrom^2'), symmetry=1, barrier=(15.2973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67575,0.0413101,-5.31841e-06,-2.05558e-08,1.00165e-11,17083.9,23.3606], Tmin=(100,'K'), Tmax=(1037.85,'K')), NASAPolynomial(coeffs=[12.0341,0.0227063,-9.2422e-06,1.75683e-09,-1.25879e-13,13785.7,-32.5286], Tmin=(1037.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CC(C)OJ)"""),
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
    E0 = (167.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (256.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (291.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (280.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (245.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (488.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (450.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (584.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (593.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (190.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (313.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (384.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', '[CH2]C=CC=C(162)'],
    products = ['C=CC(C=C)O[O](423)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.9277e+07,'m^3/(mol*s)'), n=-0.00388889, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/CdCd;Y_rad] for rate rule [C_rad/H/CdCd;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=CC(C=C)O[O](423)'],
    products = ['[CH2]C1OOC1C=C(657)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC(C=C)O[O](423)'],
    products = ['C=C[C](C=C)OO(658)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=CC(C=C)O[O](423)'],
    products = ['C=[C]C(C=C)OO(659)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(548,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=CC(C=C)O[O](423)'],
    products = ['[CH]=CC(C=C)OO(660)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.936e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2H3(44)', 'C=C[CH]O[O](661)'],
    products = ['C=CC(C=C)O[O](423)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.59671e+07,'m^3/(mol*s)'), n=0.0113737, Ea=(2.96199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/H/OneDe;Y_rad] for rate rule [C_rad/H/OneDeO;Cd_pri_rad]
Euclidian distance = 2.2360679775
family: R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(19)', 'C=C[C](C=C)O[O](662)'],
    products = ['C=CC(C=C)O[O](423)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.62e+13,'cm^3/(mol*s)'), n=0.228, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [C_rad/TwoDe;H_rad] for rate rule [C_rad/TwoDeO;H_rad]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(19)', 'C=[C]C(C=C)O[O](663)'],
    products = ['C=CC(C=C)O[O](423)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(19)', '[CH]=CC(C=C)O[O](664)'],
    products = ['C=CC(C=C)O[O](423)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.21e+14,'cm^3/(mol*s)','+|-',4.82e+13), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), comment="""From training reaction 60 used for H_rad;Cd_pri_rad
Exact match found for rate rule [Cd_pri_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=CC(C=C)O[O](423)'],
    products = ['C=CC1[CH]COO1(665)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.14356e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=CC(C=C)O[O](423)'],
    products = ['HO2(26)', 'C=C=CC=C(320)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.26e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O(20)', 'C=CC([O])C=C(666)'],
    products = ['C=CC(C=C)O[O](423)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

network(
    label = '128',
    isomers = [
        'C=CC(C=C)O[O](423)',
    ],
    reactants = [
        ('O2(2)', '[CH2]C=CC=C(162)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '128',
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

