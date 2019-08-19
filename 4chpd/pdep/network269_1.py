species(
    label = '[CH2]CC[CH]CC=CC=C(1672)',
    structure = SMILES('[CH2]CC[CH]CC=CC=C'),
    E0 = (363.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261378,0.0899612,-6.15397e-05,2.22403e-08,-3.37468e-12,43888.6,39.0906], Tmin=(100,'K'), Tmax=(1487.33,'K')), NASAPolynomial(coeffs=[15.0563,0.0487664,-1.99941e-05,3.61837e-09,-2.44611e-13,39332.1,-40.8891], Tmin=(1487.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = 'C7H10(4)',
    structure = SMILES('C=CC=CCC=C'),
    E0 = (140.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,2750,2850,1437.5,1250,1305,750,350,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.765535,'amu*angstrom^2'), symmetry=1, barrier=(17.6012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.767909,'amu*angstrom^2'), symmetry=1, barrier=(17.6557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.765812,'amu*angstrom^2'), symmetry=1, barrier=(17.6075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3438.54,'J/mol'), sigma=(5.83085,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=537.09 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928267,0.0526924,5.5285e-06,-4.52831e-08,2.12066e-11,17024.7,26.4844], Tmin=(100,'K'), Tmax=(980.344,'K')), NASAPolynomial(coeffs=[15.6612,0.0277121,-1.00058e-05,1.83677e-09,-1.31769e-13,12447.8,-52.9117], Tmin=(980.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C2H4(45)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(42.0619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]CCC1CC1[CH]C=C(1775)',
    structure = SMILES('[CH2]C=CC1CC1CC[CH2]'),
    E0 = (354.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145534,0.0697236,1.82402e-05,-7.19365e-08,3.19557e-11,42824.1,36.1442], Tmin=(100,'K'), Tmax=(986.668,'K')), NASAPolynomial(coeffs=[19.4807,0.0409388,-1.51995e-05,2.82025e-09,-2.02822e-13,36479.4,-70.8036], Tmin=(986.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC1CC=CC1[CH2](1776)',
    structure = SMILES('[CH2]CCC1CC=CC1[CH2]'),
    E0 = (320.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.206409,0.0623842,3.27591e-05,-8.33362e-08,3.54672e-11,38665.5,34.9158], Tmin=(100,'K'), Tmax=(976.866,'K')), NASAPolynomial(coeffs=[17.3246,0.0430637,-1.55381e-05,2.83152e-09,-2.01685e-13,32898.5,-59.6686], Tmin=(976.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=C[CH]C1C[CH]CCC1(1777)',
    structure = SMILES('[CH2]C=CC1C[CH]CCC1'),
    E0 = (227.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64818,0.0531183,5.18592e-05,-8.86953e-08,3.25791e-11,27556.9,30.3364], Tmin=(100,'K'), Tmax=(1053.11,'K')), NASAPolynomial(coeffs=[13.1081,0.0560609,-2.39321e-05,4.60991e-09,-3.30721e-13,22145,-43.6548], Tmin=(1053.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane) + radical(cyclohexane) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1C=CC[CH]CCC1(1745)',
    structure = SMILES('[CH2]C1C=CC[CH]CCC1'),
    E0 = (323.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896136,0.0490612,5.49855e-05,-9.36223e-08,3.60463e-11,39086.8,29.706], Tmin=(100,'K'), Tmax=(994.843,'K')), NASAPolynomial(coeffs=[11.8484,0.0517749,-1.95944e-05,3.59109e-09,-2.53168e-13,34594.2,-34.7028], Tmin=(994.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(cis-Cyclooctene) + radical(RCCJCC) + radical(Isobutyl)"""),
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
    label = '[CH2]CC=CCC=CC=C(1778)',
    structure = SMILES('[CH2]CC=CCC=CC=C'),
    E0 = (287.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.2,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.405281,0.0817332,-3.07609e-05,-1.67644e-08,1.17375e-11,34701.5,37.7775], Tmin=(100,'K'), Tmax=(1028.72,'K')), NASAPolynomial(coeffs=[18.3478,0.0400409,-1.54997e-05,2.85209e-09,-2.00408e-13,29190.9,-61.257], Tmin=(1028.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = 'C=CC=CC[CH]CC=C(1779)',
    structure = SMILES('C=CC=CC[CH]CC=C'),
    E0 = (286.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3025,407.5,1350,352.5,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.2,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.455744,0.0865599,-5.78971e-05,1.98629e-08,-2.78105e-12,34585.7,38.5019], Tmin=(100,'K'), Tmax=(1648.26,'K')), NASAPolynomial(coeffs=[18.7583,0.0399307,-1.54618e-05,2.6991e-09,-1.77702e-13,28251.8,-63.7964], Tmin=(1648.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCC=CC=CC=C(1685)',
    structure = SMILES('[CH2]CCC=CC=CC=C'),
    E0 = (268.488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3002.5,3010,3017.5,3025,975,981.25,987.5,993.75,1000,1300,1318.75,1337.5,1356.25,1375,400,425,450,475,500,1630,1642.5,1655,1667.5,1680,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (121.2,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.663285,0.0861083,-3.39134e-05,-2.04823e-08,1.48504e-11,32474.1,36.1174], Tmin=(100,'K'), Tmax=(992.866,'K')), NASAPolynomial(coeffs=[20.7262,0.0361795,-1.32383e-05,2.40173e-09,-1.69394e-13,26440.3,-75.9174], Tmin=(992.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = 'C2H4(T)(169)',
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
    label = '[CH2][CH]CC=CC=C(172)',
    structure = SMILES('[CH2][CH]CC=CC=C'),
    E0 = (411.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,395.641,395.666,1322.95],'cm^-1')),
        HinderedRotor(inertia=(0.0638666,'amu*angstrom^2'), symmetry=1, barrier=(7.08696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0638374,'amu*angstrom^2'), symmetry=1, barrier=(7.08964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0638468,'amu*angstrom^2'), symmetry=1, barrier=(7.08522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185362,'amu*angstrom^2'), symmetry=1, barrier=(20.5959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.776085,0.0623764,-3.68999e-05,8.23839e-09,2.34758e-13,49574.6,31.0578], Tmin=(100,'K'), Tmax=(1171.95,'K')), NASAPolynomial(coeffs=[12.5854,0.0327714,-1.2705e-05,2.26677e-09,-1.53515e-13,46071.7,-30.9244], Tmin=(1171.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'CH2CHCHCH(170)',
    structure = SMILES('[CH]=CC=C'),
    E0 = (346.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.31937,'amu*angstrom^2'), symmetry=1, barrier=(30.3349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64255,0.0163337,3.86225e-05,-6.71377e-08,2.83603e-11,41729.6,13.282], Tmin=(100,'K'), Tmax=(937.724,'K')), NASAPolynomial(coeffs=[12.9705,0.00669127,-1.0007e-06,1.67602e-10,-1.71436e-14,38279.7,-43.9476], Tmin=(937.724,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""CH2CHCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]CCC=C(1780)',
    structure = SMILES('[CH2]CCC=C'),
    E0 = (164.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,631.723],'cm^-1')),
        HinderedRotor(inertia=(0.00914778,'amu*angstrom^2'), symmetry=1, barrier=(2.56473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10818,'amu*angstrom^2'), symmetry=1, barrier=(2.48727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.565566,'amu*angstrom^2'), symmetry=1, barrier=(13.0035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9353,0.0370184,2.40031e-06,-2.33216e-08,9.90605e-12,19919.7,21.5221], Tmin=(100,'K'), Tmax=(1048.67,'K')), NASAPolynomial(coeffs=[9.13409,0.0280135,-1.11156e-05,2.05172e-09,-1.43458e-13,17395.2,-18.388], Tmin=(1048.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]CCC=CC=C(1771)',
    structure = SMILES('[CH2]C[CH]CCC=CC=C'),
    E0 = (363.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.261378,0.0899612,-6.15397e-05,2.22403e-08,-3.37468e-12,43888.6,39.0906], Tmin=(100,'K'), Tmax=(1487.33,'K')), NASAPolynomial(coeffs=[15.0563,0.0487664,-1.99941e-05,3.61837e-09,-2.44611e-13,39332.1,-40.8891], Tmin=(1487.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(RCCJCC)"""),
)

species(
    label = 'C=CC=CC[CH]C[CH]C(1781)',
    structure = SMILES('C=CC=CC[CH]C[CH]C'),
    E0 = (352.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28528,0.0850074,-5.56373e-05,2.01399e-08,-3.23511e-12,42564.2,38.0713], Tmin=(100,'K'), Tmax=(1338.4,'K')), NASAPolynomial(coeffs=[9.68616,0.0569111,-2.41481e-05,4.4546e-09,-3.05207e-13,40047.8,-10.0225], Tmin=(1338.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCCC=C[CH]C=C(1716)',
    structure = SMILES('[CH2]C=CC=CCCC[CH2]'),
    E0 = (274.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.507188,0.0845553,-3.22717e-05,-1.44467e-08,1.05858e-11,33245.6,38.5321], Tmin=(100,'K'), Tmax=(1036.42,'K')), NASAPolynomial(coeffs=[17.2959,0.0450518,-1.73693e-05,3.15765e-09,-2.19307e-13,27986.7,-55.562], Tmin=(1036.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.966,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCCC=CC=C(1782)',
    structure = SMILES('[CH2][CH]CCCC=CC=C'),
    E0 = (363.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.503091,0.0918059,-6.41209e-05,2.35953e-08,-3.60089e-12,43898.9,40.1368], Tmin=(100,'K'), Tmax=(1503.17,'K')), NASAPolynomial(coeffs=[16.7668,0.0458501,-1.82621e-05,3.2566e-09,-2.18273e-13,38707,-50.2194], Tmin=(1503.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(RCCJC)"""),
)

species(
    label = 'C=CC=CC[CH][CH]CC(1783)',
    structure = SMILES('C=CC=CC[CH][CH]CC'),
    E0 = (352.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44662,0.0839199,-5.49848e-05,2.05394e-08,-3.52933e-12,42558,37.3286], Tmin=(100,'K'), Tmax=(1224.72,'K')), NASAPolynomial(coeffs=[7.62247,0.0604834,-2.62808e-05,4.91481e-09,-3.39943e-13,40800.3,1.25454], Tmin=(1224.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCCC[C]=CC=C(1784)',
    structure = SMILES('[CH2]CCCC[C]=CC=C'),
    E0 = (406.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.825223,0.0972357,-7.29728e-05,2.86915e-08,-4.61987e-12,49131.3,39.0334], Tmin=(100,'K'), Tmax=(1450.66,'K')), NASAPolynomial(coeffs=[18.7108,0.0433675,-1.72721e-05,3.09348e-09,-2.08406e-13,43463.3,-62.4843], Tmin=(1450.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCCCC=[C]C=C(1785)',
    structure = SMILES('[CH2]CCCCC=[C]C=C'),
    E0 = (368.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.834007,0.0988299,-7.64391e-05,3.17813e-08,-5.45078e-12,44458.4,38.1831], Tmin=(100,'K'), Tmax=(1365.62,'K')), NASAPolynomial(coeffs=[17.3818,0.0454744,-1.78332e-05,3.17113e-09,-2.13185e-13,39483.3,-55.3737], Tmin=(1365.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = 'C=C[CH][CH]C=CCCC(1698)',
    structure = SMILES('[CH2]C=CC=C[CH]CCC'),
    E0 = (210.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.321388,0.0783178,-1.19686e-05,-3.62118e-08,1.83617e-11,25527.5,34.9714], Tmin=(100,'K'), Tmax=(1012.37,'K')), NASAPolynomial(coeffs=[17.3886,0.0451187,-1.7267e-05,3.1592e-09,-2.21624e-13,20057.1,-59.9943], Tmin=(1012.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(Allyl_S) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]CCCCC=C[C]=C(1786)',
    structure = SMILES('[CH2]CCCCC=C[C]=C'),
    E0 = (368.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.834007,0.0988299,-7.64391e-05,3.17813e-08,-5.45078e-12,44458.4,38.1831], Tmin=(100,'K'), Tmax=(1365.62,'K')), NASAPolynomial(coeffs=[17.3818,0.0454744,-1.78332e-05,3.17113e-09,-2.13185e-13,39483.3,-55.3737], Tmin=(1365.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(RCCJ)"""),
)

species(
    label = 'C=CC=[C]C[CH]CCC(1787)',
    structure = SMILES('C=CC=[C]C[CH]CCC'),
    E0 = (396.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0517867,0.0905894,-6.49186e-05,2.56668e-08,-4.3927e-12,47797.2,37.0232], Tmin=(100,'K'), Tmax=(1305.58,'K')), NASAPolynomial(coeffs=[11.8611,0.0540909,-2.29847e-05,4.25402e-09,-2.92445e-13,44686.6,-23.626], Tmin=(1305.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(RCCJCC)"""),
)

species(
    label = '[CH]=CC=CCCCC[CH2](1788)',
    structure = SMILES('[CH]=CC=CCCCC[CH2]'),
    E0 = (416.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.15462,0.0993845,-7.53806e-05,2.95625e-08,-4.67821e-12,50260.7,40.1576], Tmin=(100,'K'), Tmax=(1495.87,'K')), NASAPolynomial(coeffs=[21.4026,0.0390658,-1.48953e-05,2.60587e-09,-1.73037e-13,43512.2,-77.7519], Tmin=(1495.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]=CC[CH]CCC(1789)',
    structure = SMILES('C=C[C]=CC[CH]CCC'),
    E0 = (357.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145189,0.0930282,-7.07211e-05,3.10821e-08,-5.97728e-12,43128.4,36.488], Tmin=(100,'K'), Tmax=(1178.23,'K')), NASAPolynomial(coeffs=[10.8041,0.0558556,-2.33962e-05,4.30437e-09,-2.95448e-13,40548.3,-18.1319], Tmin=(1178.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJCC) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C=CC[CH]CCC(1790)',
    structure = SMILES('C=[C]C=CC[CH]CCC'),
    E0 = (357.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.145189,0.0930282,-7.07211e-05,3.10821e-08,-5.97728e-12,43128.4,36.488], Tmin=(100,'K'), Tmax=(1178.23,'K')), NASAPolynomial(coeffs=[10.8041,0.0558556,-2.33962e-05,4.30437e-09,-2.95448e-13,40548.3,-18.1319], Tmin=(1178.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJCC) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]CCC1C[CH]C1C=C(1791)',
    structure = SMILES('[CH2]CCC1C[CH]C1C=C'),
    E0 = (398.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.238591,0.0637753,2.17991e-05,-6.51415e-08,2.68617e-11,48119.5,37.2874], Tmin=(100,'K'), Tmax=(1021.74,'K')), NASAPolynomial(coeffs=[15.6129,0.0475563,-1.89413e-05,3.55943e-09,-2.53657e-13,42682.7,-48.4467], Tmin=(1021.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC1C[CH]C=CC1(1792)',
    structure = SMILES('[CH2]CCC1C[CH]C=CC1'),
    E0 = (242.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.489443,0.0534511,5.84425e-05,-1.05952e-07,4.17872e-11,29277,30.2518], Tmin=(100,'K'), Tmax=(996.593,'K')), NASAPolynomial(coeffs=[16.6216,0.0464966,-1.80785e-05,3.42635e-09,-2.48291e-13,23191.5,-61.9207], Tmin=(996.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(RCCJ) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C=CC1[CH]C[CH]CCC1(1793)',
    structure = SMILES('C=CC1[CH]C[CH]CCC1'),
    E0 = (275.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.626545,0.0212703,0.000206134,-3.09756e-07,1.28788e-10,33321.2,36.4218], Tmin=(100,'K'), Tmax=(916.091,'K')), NASAPolynomial(coeffs=[35.6463,0.00552588,7.32099e-06,-1.63058e-09,9.73038e-14,21149.3,-160.873], Tmin=(916.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(cycloheptane) + radical(cycloheptane)"""),
)

species(
    label = '[CH]1C=CCCCC[CH]C1(1704)',
    structure = SMILES('[CH]1C=CCCCC[CH]C1'),
    E0 = (287.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.19911,0.0858676,-4.9871e-05,1.3949e-08,-1.5696e-12,34730.3,19.7747], Tmin=(100,'K'), Tmax=(1987.44,'K')), NASAPolynomial(coeffs=[21.1878,0.0428229,-1.7383e-05,3.05113e-09,-1.98735e-13,26229.4,-98.0942], Tmin=(1987.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclononane) + radical(RCCJCC) + radical(Allyl_S)"""),
)

species(
    label = 'C=CC=CCC=CCC(1794)',
    structure = SMILES('C=CC=CCC=CCC'),
    E0 = (81.8506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.437733,0.0809491,-2.01632e-05,-2.91551e-08,1.62078e-11,10018.7,36.1071], Tmin=(100,'K'), Tmax=(1017.73,'K')), NASAPolynomial(coeffs=[18.5083,0.0424101,-1.631e-05,3.00441e-09,-2.1181e-13,4301.75,-64.77], Tmin=(1017.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.8506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=CCCCC=C(1795)',
    structure = SMILES('C=CC=CCCCC=C'),
    E0 = (91.6917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.565335,0.0847311,-3.08977e-05,-1.77733e-08,1.20851e-11,11205.9,36.3093], Tmin=(100,'K'), Tmax=(1033.76,'K')), NASAPolynomial(coeffs=[18.5555,0.0428306,-1.66554e-05,3.0654e-09,-2.15112e-13,5538.19,-64.8649], Tmin=(1033.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.6917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=CC=CCCC(1707)',
    structure = SMILES('C=CC=CC=CCCC'),
    E0 = (63.2418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.700265,0.0853725,-2.34576e-05,-3.27345e-08,1.9284e-11,7791.44,34.4636], Tmin=(100,'K'), Tmax=(988.167,'K')), NASAPolynomial(coeffs=[20.9359,0.038468,-1.40034e-05,2.54358e-09,-1.7994e-13,1529.45,-79.7094], Tmin=(988.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.2418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]CC[CH]C=CC=CC(1796)',
    structure = SMILES('[CH2]CCC=C[CH]C=CC'),
    E0 = (277.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.172363,0.0762872,-1.35884e-05,-3.04901e-08,1.54292e-11,33508.4,38.027], Tmin=(100,'K'), Tmax=(1032.88,'K')), NASAPolynomial(coeffs=[16.0408,0.0464535,-1.81213e-05,3.32606e-09,-2.32569e-13,28401.2,-49.2264], Tmin=(1032.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJC=C) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])CC=CC=C(1675)',
    structure = SMILES('[CH2]CC([CH2])CC=CC=C'),
    E0 = (368.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3743.15,'J/mol'), sigma=(6.69575,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=584.67 K, Pc=28.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.775727,0.0920908,-5.23301e-05,2.97867e-09,5.55791e-12,44508.8,39.8767], Tmin=(100,'K'), Tmax=(1022.15,'K')), NASAPolynomial(coeffs=[18.2533,0.0423809,-1.57109e-05,2.78975e-09,-1.91172e-13,39325.4,-58.6703], Tmin=(1022.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCC([CH2])C=CC=C(1797)',
    structure = SMILES('[CH2]CCC([CH2])C=CC=C'),
    E0 = (366.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.769282,0.0917391,-5.06093e-05,5.34398e-10,6.61069e-12,44234.8,39.8879], Tmin=(100,'K'), Tmax=(1012.92,'K')), NASAPolynomial(coeffs=[18.4077,0.0420072,-1.54624e-05,2.74112e-09,-1.87979e-13,39016.2,-59.4594], Tmin=(1012.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=CCC1CCC1(1679)',
    structure = SMILES('C=CC=CCC1CCC1'),
    E0 = (113.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.154891,0.0586216,5.40987e-05,-1.11174e-07,4.62934e-11,13757.9,32.5359], Tmin=(100,'K'), Tmax=(973.003,'K')), NASAPolynomial(coeffs=[20.3315,0.0392334,-1.39937e-05,2.61432e-09,-1.91826e-13,6822.89,-79.7129], Tmin=(973.003,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(544.598,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH2]CC[CH]CC1C=CC1(1798)',
    structure = SMILES('[CH2]CC[CH]CC1C=CC1'),
    E0 = (417.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.207,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0464421,0.0737579,-1.76783e-05,-1.83661e-08,9.28774e-12,50342.8,38.1346], Tmin=(100,'K'), Tmax=(1115.36,'K')), NASAPolynomial(coeffs=[13.653,0.0504114,-2.05074e-05,3.78268e-09,-2.62199e-13,45724.5,-36.0912], Tmin=(1115.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH2](171)',
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
    label = '[CH]CC=CC=C(1799)',
    structure = SMILES('[CH]CC=CC=C'),
    E0 = (483.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,410.913,410.916,410.924,410.928,410.948],'cm^-1')),
        HinderedRotor(inertia=(0.000998349,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12375,'amu*angstrom^2'), symmetry=1, barrier=(14.828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123754,'amu*angstrom^2'), symmetry=1, barrier=(14.828,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29337,0.0479245,-6.89739e-06,-2.77603e-08,1.47696e-11,58255.3,23.4321], Tmin=(100,'K'), Tmax=(973.323,'K')), NASAPolynomial(coeffs=[14.4936,0.0209449,-7.34215e-06,1.32785e-09,-9.48044e-14,54394,-46.5294], Tmin=(973.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CH2(T)(33)',
    structure = SMILES('[CH2]'),
    E0 = (381.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([971.045,2816.03,3444.23],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(381.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C[CH]CC=CC=C(1773)',
    structure = SMILES('[CH2]C[CH]CC=CC=C'),
    E0 = (387.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292654,0.0761726,-5.07746e-05,1.76972e-08,-2.56775e-12,46730.2,34.8634], Tmin=(100,'K'), Tmax=(1560.84,'K')), NASAPolynomial(coeffs=[14.3974,0.040026,-1.60369e-05,2.86e-09,-1.91282e-13,42327.2,-39.4636], Tmin=(1560.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJ)"""),
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
    E0 = (363.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (387.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (559.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (530.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (490.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (504.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (507.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (501.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (473.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (474.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (517.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (522.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (515.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (463.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (507.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (481.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (541.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (412.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (423.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (456.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (538.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (561.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (641.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (551.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (729.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (503.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (408.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (459.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (420.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (427.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (427.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (388.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (480.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (528.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (526.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (371.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (486.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (773.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (768.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C7H10(4)', 'C2H4(45)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CCC1CC1[CH]C=C(1775)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.56e+12,'s^-1'), n=0.21, Ea=(24.309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 125 used for R4_S_D;doublebond_intra_HCd_pri;radadd_intra_csHNd
Exact match found for rate rule [R4_S_D;doublebond_intra_HCd_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CCC1CC=CC1[CH2](1776)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6_SSM_D;doublebond_intra_2H_pri;radadd_intra_cs] for rate rule [R6_SSM_D;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=C[CH]C1C[CH]CCC1(1777)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.69e+11,'s^-1'), n=0.22, Ea=(167.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7;doublebond_intra_HCd_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]C1C=CC[CH]CCC1(1745)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.78819e+09,'s^-1'), n=0.5525, Ea=(127.194,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_2H_pri;radadd_intra_cs2H] for rate rule [R9;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(19)', '[CH2]CC=CCC=CC=C(1778)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(19)', 'C=CC=CC[CH]CC=C(1779)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(19)', '[CH2]CCC=CC=CC=C(1685)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0272924,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C7H10(4)', 'C2H4(T)(169)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0344575,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]CC=CC=C(172)', 'C2H4(45)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.00310016,'m^3/(mol*s)'), n=2.49247, Ea=(21.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CH2CHCHCH(170)', '[CH2]CCC=C(1780)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(14400,'cm^3/(mol*s)'), n=2.41, Ea=(6.23416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C[CH]CCC=CC=C(1771)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC=CC[CH]C[CH]C(1781)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CCCC=C[CH]C=C(1716)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.2802e+09,'s^-1'), n=0.983917, Ea=(99.768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/(Cd-Cd-Cd)] + [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2][CH]CCCC=CC=C(1782)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC=CC[CH][CH]CC(1783)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CCCC[C]=CC=C(1784)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.37e+07,'s^-1'), n=1.41, Ea=(177.82,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 207 used for R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CCCCC=[C]C=C(1785)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=C[CH][CH]C=CCCC(1698)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(62296.1,'s^-1'), n=1.86, Ea=(59.4128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/Cd] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CCCCC=C[C]=C(1786)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(252000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;Cs_H_out_1H] for rate rule [R5H_SMSS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=CC=[C]C[CH]CCC(1787)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.59786e+07,'s^-1'), n=1.58088, Ea=(142.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=CC=CCCCC[CH2](1788)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R6H_RSMSR;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 3.16227766017
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=C[C]=CC[CH]CCC(1789)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.59e+07,'s^-1'), n=1.4638, Ea=(277.467,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_singleDe] for rate rule [R7HJ_3;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=[C]C=CC[CH]CCC(1790)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.11562e+07,'s^-1'), n=1.49763, Ea=(188.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_doubleC] for rate rule [R8Hall;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]CC=CC=C(172)', 'C2H4(T)(169)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.47003e+07,'m^3/(mol*s)'), n=0.0284742, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CCC1C[CH]C1C=C(1791)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.9e+11,'s^-1'), n=0.19, Ea=(139.662,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 233 used for R4_Cs_HH_D;doublebond_intra_pri_HCd;radadd_intra_csHCs
Exact match found for rate rule [R4_Cs_HH_D;doublebond_intra_pri_HCd;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CCC1C[CH]C=CC1(1792)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.40768e+07,'s^-1'), n=0.903766, Ea=(44.7452,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_SSR;doublebond_intra_pri_2H;radadd_intra_csHCs] + [R6_SSM_D;doublebond_intra_pri_2H;radadd_intra_cs] for rate rule [R6_SSM_D;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC1[CH]C[CH]CCC1(1793)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.54269e+17,'s^-1'), n=-0.500262, Ea=(96.3509,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri_HCd;radadd_intra_cs] for rate rule [R7_linear;doublebond_intra_pri_HCd;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH]1C=CCCCC[CH]C1(1704)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.04798e+07,'s^-1'), n=0.841667, Ea=(57.3368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_2H;radadd_intra_cs2H] for rate rule [R9_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC=CCC=CCC(1794)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC=CCCCC=C(1795)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC=CC=CCCC(1707)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CC[CH]C=CC=CC(1796)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.94565e+09,'s^-1'), n=0.909333, Ea=(116.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [1_3_unsaturated_pentane_backbone;CH(CJ)_1;unsaturated_end] + [1_3_pentadiene;CH_end;unsaturated_end] for rate rule [1_3_pentadiene;CH(CJ)_1;CdH2_2]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]CC([CH2])CC=CC=C(1675)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]CCC([CH2])C=CC=C(1797)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['C=CC=CCC1CCC1(1679)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]CC[CH]CC=CC=C(1672)'],
    products = ['[CH2]CC[CH]CC1C=CC1(1798)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1,3-butadiene_backbone;C=C_1;C=C_2] for rate rule [1,3-butadiene_backbone;CdH(C)_1;CdH2_2]
Euclidian distance = 1.41421356237
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C[CH2](171)', '[CH]CC=CC=C(1799)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.4725e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(T)(33)', '[CH2]C[CH]CC=CC=C(1773)'],
    products = ['[CH2]CC[CH]CC=CC=C(1672)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.23625e+06,'m^3/(mol*s)'), n=0.36814, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

network(
    label = '269',
    isomers = [
        '[CH2]CC[CH]CC=CC=C(1672)',
    ],
    reactants = [
        ('C7H10(4)', 'C2H4(45)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '269',
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

