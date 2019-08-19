species(
    label = 'C=C[C]=CC([O])=O(26472)',
    structure = SMILES('C=C[C]=CC([O])=O'),
    E0 = (222.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180,180,1213.68,1213.69],'cm^-1')),
        HinderedRotor(inertia=(0.171321,'amu*angstrom^2'), symmetry=1, barrier=(3.939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171316,'amu*angstrom^2'), symmetry=1, barrier=(3.93889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53489,0.0631316,-0.000103752,9.98059e-08,-3.70891e-11,26856.1,23.2184], Tmin=(100,'K'), Tmax=(843.138,'K')), NASAPolynomial(coeffs=[3.55716,0.0338364,-1.65843e-05,3.16866e-09,-2.17641e-13,27215.4,17.9599], Tmin=(843.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'CO2(13)',
    structure = SMILES('O=C=O'),
    E0 = (-403.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.166,1086.67,1086.68,2300.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2779,0.00275783,7.12787e-06,-1.07855e-08,4.14228e-12,-48475.6,5.97856], Tmin=(100,'K'), Tmax=(988.185,'K')), NASAPolynomial(coeffs=[4.55071,0.00290728,-1.14643e-06,2.25798e-10,-1.69526e-14,-48986,-1.45662], Tmin=(988.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.131,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: BurkeH2O2"""),
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
    label = 'C=CC1=CC1([O])[O](27895)',
    structure = SMILES('C=CC1=CC1([O])[O]'),
    E0 = (400.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67236,0.0550864,-5.70235e-05,3.32704e-08,-8.21938e-12,48269.1,19.222], Tmin=(100,'K'), Tmax=(955.514,'K')), NASAPolynomial(coeffs=[8.37167,0.0270408,-1.29954e-05,2.55101e-09,-1.8177e-13,46988.9,-12.7933], Tmin=(955.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH2]C1[C]=CC(=O)O1(27739)',
    structure = SMILES('[CH2]C1[C]=CC(=O)O1'),
    E0 = (226.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20897,0.0315553,3.60692e-09,-1.89118e-08,8.27298e-12,27365,23.653], Tmin=(100,'K'), Tmax=(1063.97,'K')), NASAPolynomial(coeffs=[9.69765,0.0205468,-8.64815e-06,1.65482e-09,-1.18255e-13,24801,-17.5004], Tmin=(1063.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + ring(Cyclopentane) + radical(Cds_S) + radical(CJC(C)OC)"""),
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
    label = 'C=C=C=CC([O])=O(27914)',
    structure = SMILES('C=C=C=CC([O])=O'),
    E0 = (253.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,540,563.333,586.667,610,1970,2140,180,180,180,180,1739.77,1746.08],'cm^-1')),
        HinderedRotor(inertia=(0.45788,'amu*angstrom^2'), symmetry=1, barrier=(10.5276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68733,0.0626774,-0.000116836,1.21481e-07,-4.73569e-11,30502.8,22.1555], Tmin=(100,'K'), Tmax=(836.27,'K')), NASAPolynomial(coeffs=[1.43866,0.0360259,-1.90934e-05,3.75218e-09,-2.61173e-13,31518,29.1314], Tmin=(836.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = 'C=CC#CC([O])=O(27915)',
    structure = SMILES('C=CC#CC([O])=O'),
    E0 = (193.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,2100,2250,500,550,180,180,180,180,1185.39,1185.43,1185.45],'cm^-1')),
        HinderedRotor(inertia=(0.00731717,'amu*angstrom^2'), symmetry=1, barrier=(7.29679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41292,'amu*angstrom^2'), symmetry=1, barrier=(55.4778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08437,0.0467139,-6.1398e-05,5.51853e-08,-2.10157e-11,23390.8,22.2948], Tmin=(100,'K'), Tmax=(769.201,'K')), NASAPolynomial(coeffs=[4.05708,0.0295452,-1.44427e-05,2.8102e-09,-1.97329e-13,23291.8,14.624], Tmin=(769.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)O2s) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(669)',
    structure = SMILES('[O][C]=O'),
    E0 = (31.9507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90478,-0.000175995,8.15126e-06,-1.13656e-08,4.4768e-12,3848.25,8.04855], Tmin=(100,'K'), Tmax=(975.388,'K')), NASAPolynomial(coeffs=[5.59398,-0.00122084,7.11747e-07,-9.7712e-11,3.97995e-15,3238.91,-1.49318], Tmin=(975.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.9507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = 'C=[C]C=CC([O])=O(27916)',
    structure = SMILES('C=C=C[CH]C([O])=O'),
    E0 = (124.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81609,0.0448875,-2.38578e-05,3.37656e-09,2.78035e-13,15040.1,19.7756], Tmin=(100,'K'), Tmax=(1645.55,'K')), NASAPolynomial(coeffs=[16.6486,0.0197535,-1.09017e-05,2.16064e-09,-1.49951e-13,8679.98,-63.6628], Tmin=(1645.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'C=CC=[C]C([O])=O(27917)',
    structure = SMILES('C=CC=C=C([O])[O]'),
    E0 = (151.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.91501,0.060358,-5.5624e-05,2.11389e-08,-1.73774e-12,18338.4,22.3943], Tmin=(100,'K'), Tmax=(1006.78,'K')), NASAPolynomial(coeffs=[16.4865,0.0122432,-4.42632e-06,8.0386e-10,-5.70889e-14,14506.1,-56.2954], Tmin=(1006.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=CC=CC([O])=O(27918)',
    structure = SMILES('[CH]=CC=CC([O])=O'),
    E0 = (270.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,213.919,213.919,213.92,213.92,1726.42,4000],'cm^-1')),
        HinderedRotor(inertia=(0.225046,'amu*angstrom^2'), symmetry=1, barrier=(7.30803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.888843,'amu*angstrom^2'), symmetry=1, barrier=(28.8637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58208,0.0596229,-8.97746e-05,8.27475e-08,-3.08767e-11,32641.9,23.855], Tmin=(100,'K'), Tmax=(790.768,'K')), NASAPolynomial(coeffs=[4.99706,0.0315073,-1.58775e-05,3.11019e-09,-2.18297e-13,32440.8,10.3246], Tmin=(790.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C[C]=[C]C(=O)O(27919)',
    structure = SMILES('[CH2][CH]C#CC(=O)O'),
    E0 = (198.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53185,0.0585423,-8.36816e-05,6.78836e-08,-2.12972e-11,23967.5,27.4601], Tmin=(100,'K'), Tmax=(940.645,'K')), NASAPolynomial(coeffs=[7.01889,0.0233915,-8.78342e-06,1.44465e-09,-8.96363e-14,23458.1,4.10307], Tmin=(940.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)O2s) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=[C][C]=CC(=O)O(27920)',
    structure = SMILES('[CH2]C#C[CH]C(=O)O'),
    E0 = (115.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03134,0.0326345,7.37325e-06,-3.03998e-08,1.26576e-11,13933.1,25.8864], Tmin=(100,'K'), Tmax=(1056.97,'K')), NASAPolynomial(coeffs=[11.8839,0.0197745,-9.04068e-06,1.8171e-09,-1.33956e-13,10486,-28.6466], Tmin=(1056.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CtHH) + group(Cs-CtHHH) + group(Cds-OdCsOs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Propargyl) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C[C]=CC(=O)O(27921)',
    structure = SMILES('[CH]C=C=CC(=O)O'),
    E0 = (221.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35386,0.064859,-9.19574e-05,8.3056e-08,-3.07351e-11,26716.2,25.2432], Tmin=(100,'K'), Tmax=(800.177,'K')), NASAPolynomial(coeffs=[4.61778,0.0369994,-1.80928e-05,3.48683e-09,-2.42508e-13,26563.5,12.5336], Tmin=(800.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]C(=O)C=C1[CH]C1(27922)',
    structure = SMILES('[O]C(=O)[CH]C1=CC1'),
    E0 = (213.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7141,0.0408002,-7.78718e-06,-1.44439e-08,6.33181e-12,25740.4,19.6053], Tmin=(100,'K'), Tmax=(1206.25,'K')), NASAPolynomial(coeffs=[13.3116,0.0233934,-1.23191e-05,2.52872e-09,-1.84377e-13,21411,-44.8688], Tmin=(1206.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclopropene) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C[C]=C[C]1OO1(27923)',
    structure = SMILES('C=C[C]=C[C]1OO1'),
    E0 = (496.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.236503,0.0609106,-5.94472e-05,2.83051e-08,-4.96436e-12,59844.3,27.0233], Tmin=(100,'K'), Tmax=(1680.04,'K')), NASAPolynomial(coeffs=[16.3476,0.00850298,-1.12613e-07,-2.17116e-10,2.05741e-14,56413.5,-53.1617], Tmin=(1680.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cs_P) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=C[C]([O])O1(27924)',
    structure = SMILES('[CH2]C=C1[CH]C(=O)O1'),
    E0 = (-12.4008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45715,0.0363233,3.24277e-05,-7.92005e-08,3.57898e-11,-1381.95,20.3789], Tmin=(100,'K'), Tmax=(942.707,'K')), NASAPolynomial(coeffs=[19.7647,0.00806752,-1.2564e-06,2.35787e-10,-2.51241e-14,-7029.89,-78.5128], Tmin=(942.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.4008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(4-Methylene-2-oxetanone) + radical(Allyl_P) + radical(C=CCJCO)"""),
)

species(
    label = 'O=C1C=[C][CH]CO1(27912)',
    structure = SMILES('O=C1[CH][C]=CCO1'),
    E0 = (44.2624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80082,0.00683781,9.21028e-05,-1.12154e-07,3.8322e-11,5383.44,13.3297], Tmin=(100,'K'), Tmax=(1056.58,'K')), NASAPolynomial(coeffs=[10.8039,0.0296158,-1.55853e-05,3.33767e-09,-2.5411e-13,729.653,-39.7408], Tmin=(1056.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.2624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclohexane) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'C=C=C=CC(=O)O(27925)',
    structure = SMILES('C=C=C=CC(=O)O'),
    E0 = (27.3128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36661,0.0649028,-0.000102566,9.36155e-08,-3.39709e-11,3373.03,23.1415], Tmin=(100,'K'), Tmax=(808.246,'K')), NASAPolynomial(coeffs=[6.01924,0.0298049,-1.50249e-05,2.92883e-09,-2.04397e-13,3015.25,4.12508], Tmin=(808.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.3128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1=CC(=O)O1(27768)',
    structure = SMILES('C=CC1=CC(=O)O1'),
    E0 = (-32.8816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55603,0.0494598,-4.22028e-05,1.7887e-08,-3.024e-12,-3863.32,19.0962], Tmin=(100,'K'), Tmax=(1408.85,'K')), NASAPolynomial(coeffs=[13.0856,0.0167248,-7.34979e-06,1.39448e-09,-9.73854e-14,-7111.99,-40.4793], Tmin=(1408.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.8816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = 'C=C[C]=C[C]=O(27558)',
    structure = SMILES('[CH2]C=[C]C=C=O'),
    E0 = (309.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.298345,'amu*angstrom^2'), symmetry=1, barrier=(62.8412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.73562,'amu*angstrom^2'), symmetry=1, barrier=(62.8973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12659,0.0437529,-4.93794e-05,3.76269e-08,-1.22823e-11,37282,19.9963], Tmin=(100,'K'), Tmax=(833.068,'K')), NASAPolynomial(coeffs=[5.18953,0.0256982,-1.08424e-05,1.96348e-09,-1.32238e-13,36887.8,6.47605], Tmin=(833.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]C1OC1=O(27756)',
    structure = SMILES('[CH2]C=[C]C1OC1=O'),
    E0 = (219.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02412,0.0445143,-3.1469e-05,1.13819e-08,-1.73636e-12,26489.7,21.972], Tmin=(100,'K'), Tmax=(1455.3,'K')), NASAPolynomial(coeffs=[9.14844,0.0249325,-1.12857e-05,2.13599e-09,-1.48038e-13,24416.1,-15.0718], Tmin=(1455.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(2(co)oxirane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[O]C1([O])C=C=CC1(27926)',
    structure = SMILES('[O]C1([O])C=C=CC1'),
    E0 = (516.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13919,0.0536308,-5.1175e-05,2.50556e-08,-4.74759e-12,62196.4,21.5157], Tmin=(100,'K'), Tmax=(1396.71,'K')), NASAPolynomial(coeffs=[14.1919,0.0127382,-3.48727e-06,4.93709e-10,-2.90379e-14,58892.7,-44.591], Tmin=(1396.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(C=CC(C)(O)OJ) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = 'C[C]=C=CC([O])=O(27927)',
    structure = SMILES('CC#CC=C([O])[O]'),
    E0 = (145.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28877,0.0524707,-3.93188e-05,8.45344e-09,1.96637e-12,17606.2,22.5076], Tmin=(100,'K'), Tmax=(983.613,'K')), NASAPolynomial(coeffs=[14.1373,0.0152169,-5.3769e-06,9.48966e-10,-6.59353e-14,14353.2,-42.9545], Tmin=(983.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'CC=C=[C]C([O])=O(27928)',
    structure = SMILES('CC=C=C=C([O])[O]'),
    E0 = (204.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28268,0.0626832,-7.33834e-05,4.58674e-08,-1.1621e-11,24664.3,21.1159], Tmin=(100,'K'), Tmax=(954.057,'K')), NASAPolynomial(coeffs=[10.7877,0.0228313,-1.07255e-05,2.08293e-09,-1.47521e-13,22850.6,-24.2931], Tmin=(954.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(=O)C1[C]=CC1(27929)',
    structure = SMILES('[O]C(=O)C1[C]=CC1'),
    E0 = (243.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37158,0.0350868,-1.36975e-05,1.23377e-10,5.93758e-13,29345.1,19.9012], Tmin=(100,'K'), Tmax=(1695.27,'K')), NASAPolynomial(coeffs=[11.5337,0.0223872,-1.03521e-05,1.91109e-09,-1.27518e-13,24957,-32.917], Tmin=(1695.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + ring(Cyclobutene) + radical(CCOJ) + radical(cyclobutene-vinyl)"""),
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
    label = '[CH]=C=CC([O])=O(27930)',
    structure = SMILES('C#CC=C([O])[O]'),
    E0 = (187.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63825,'amu*angstrom^2'), symmetry=1, barrier=(37.6666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6892,0.044488,-4.05979e-05,1.36739e-08,-3.23674e-13,22664.3,17.6502], Tmin=(100,'K'), Tmax=(1005.71,'K')), NASAPolynomial(coeffs=[14.2148,0.00683295,-2.5767e-06,4.95366e-10,-3.68989e-14,19529.7,-45.9081], Tmin=(1005.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ) + radical(C=COJ)"""),
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
    E0 = (222.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (400.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (321.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (481.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (419.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (315.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (222.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (437.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (434.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (493.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (333.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (252.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (351.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (483.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (360.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (496.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (360.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (283.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (247.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (230.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (552.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (268.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (516.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (430.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (460.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (341.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (569.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['CO2(13)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=CC1=CC1([O])[O](27895)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.03354e+12,'s^-1'), n=0.0175, Ea=(178.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleDe] for rate rule [R4_D_CO;carbonylbond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 177.9 to 178.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['[CH2]C1[C]=CC(=O)O1(27739)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.36433e+10,'s^-1'), n=0.527281, Ea=(99.1325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_O] + [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C=C=CC([O])=O(27914)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC#CC([O])=O(27915)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C]=O(669)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CO2(13)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(174.175,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CJ] for rate rule [CO2;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 171.2 to 174.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=[C]C=CC([O])=O(27916)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=CC=[C]C([O])=O(27917)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC=CC([O])=O(27918)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=C[C]=[C]C(=O)O(27919)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.99823e+07,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=[C][C]=CC(=O)O(27920)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.316e+09,'s^-1'), n=0.699, Ea=(29.5516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_doubleC] for rate rule [R5HJ_3;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['[CH]=C[C]=CC(=O)O(27921)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.854,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=O(669)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['[O]C(=O)C=C1[CH]C1(27922)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=C[C]=C[C]1OO1(27923)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(273.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_De;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic
Ea raised from 272.4 to 273.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=CC1=C[C]([O])O1(27924)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.906e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleDe]
Euclidian distance = 1.73205080757
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['O=C1C=[C][CH]CO1(27912)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.98334e+10,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=C=C=CC(=O)O(27925)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C=CC1=CC(=O)O1(27768)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleDe_rad_out;Opri_rad]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O(4)', 'C=C[C]=C[C]=O(27558)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['[CH2]C=[C]C1OC1=O(27756)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.28489e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_(CO)_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['[O]C1([O])C=C=CC1(27926)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(293.58,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R6_SMM;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SMM;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 292.1 to 293.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['C[C]=C=CC([O])=O(27927)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['CC=C=[C]C([O])=O(27928)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC([O])=O(26472)'],
    products = ['[O]C(=O)C1[C]=CC1(27929)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(19)', '[CH]=C=CC([O])=O(27930)'],
    products = ['C=C[C]=CC([O])=O(26472)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4438',
    isomers = [
        'C=C[C]=CC([O])=O(26472)',
    ],
    reactants = [
        ('CO2(13)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4438',
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

