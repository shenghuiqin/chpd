species(
    label = 'C=C[C]=CO[CH]C(26656)',
    structure = SMILES('C=C[C]=CO[CH]C'),
    E0 = (259.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10735,'amu*angstrom^2'), symmetry=1, barrier=(25.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10718,'amu*angstrom^2'), symmetry=1, barrier=(25.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10807,'amu*angstrom^2'), symmetry=1, barrier=(25.4768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10748,'amu*angstrom^2'), symmetry=1, barrier=(25.4632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28817,0.0907966,-9.47325e-05,4.69616e-08,-8.60356e-12,31385.4,29.3827], Tmin=(100,'K'), Tmax=(1566.16,'K')), NASAPolynomial(coeffs=[24.6208,0.00800332,4.82546e-07,-3.44973e-10,2.95007e-14,25308.4,-100.728], Tmin=(1566.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJC=C)"""),
)

species(
    label = 'CH3CHO(37)',
    structure = SMILES('CC=O'),
    E0 = (-178.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1305.64,1305.66,1305.67,3976.84],'cm^-1')),
        HinderedRotor(inertia=(0.136163,'amu*angstrom^2'), symmetry=1, barrier=(3.13064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-178.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C1[C]=COC1C(28544)',
    structure = SMILES('[CH2]C1[C]=COC1C'),
    E0 = (246.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45054,0.0340436,6.11985e-05,-1.19165e-07,5.39976e-11,29726.6,20.2616], Tmin=(100,'K'), Tmax=(893.012,'K')), NASAPolynomial(coeffs=[20.003,0.0101371,1.92556e-06,-6.884e-10,4.97403e-14,24052.8,-80.3591], Tmin=(893.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(Isobutyl) + radical(Cds_S)"""),
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
    label = 'C=C[C]=COC=C(26641)',
    structure = SMILES('C=C[C]=COC=C'),
    E0 = (220.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29283,'amu*angstrom^2'), symmetry=1, barrier=(29.7247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29571,'amu*angstrom^2'), symmetry=1, barrier=(29.7909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29537,'amu*angstrom^2'), symmetry=1, barrier=(29.7831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.643595,0.0606377,-2.31174e-05,-2.5767e-08,1.82062e-11,26662.9,25.3194], Tmin=(100,'K'), Tmax=(919.334,'K')), NASAPolynomial(coeffs=[19.4156,0.0135456,-2.70938e-06,3.5352e-10,-2.44145e-14,21749.8,-71.6143], Tmin=(919.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C=C=CO[CH]C(28545)',
    structure = SMILES('C=C=C=CO[CH]C'),
    E0 = (289.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,563.333,586.667,610,1970,2140,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1327,'amu*angstrom^2'), symmetry=1, barrier=(26.0429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13298,'amu*angstrom^2'), symmetry=1, barrier=(26.0493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13205,'amu*angstrom^2'), symmetry=1, barrier=(26.0281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.714302,0.0853758,-9.02929e-05,4.53554e-08,-8.54445e-12,35013.6,26.8059], Tmin=(100,'K'), Tmax=(1458.64,'K')), NASAPolynomial(coeffs=[23.8554,0.00832432,-1.10812e-06,4.68072e-11,3.98826e-16,28875,-97.4763], Tmin=(1458.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=CC#CO[CH]C(28546)',
    structure = SMILES('C=CC#CO[CH]C'),
    E0 = (260.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,2100,2250,500,550,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15766,'amu*angstrom^2'), symmetry=1, barrier=(26.6169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1582,'amu*angstrom^2'), symmetry=1, barrier=(26.6293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15755,'amu*angstrom^2'), symmetry=1, barrier=(26.6143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15741,'amu*angstrom^2'), symmetry=1, barrier=(26.6111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.380149,0.0637644,-2.43619e-05,-2.78314e-08,1.89592e-11,31503.1,24.1062], Tmin=(100,'K'), Tmax=(946.423,'K')), NASAPolynomial(coeffs=[22.1345,0.0109095,-2.54302e-06,4.38372e-10,-3.5732e-14,25634.8,-88.8964], Tmin=(946.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(CCsJOCs)"""),
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
    label = '[CH2]COC=[C]C=C(28547)',
    structure = SMILES('[CH2]COC=[C]C=C'),
    E0 = (276.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13421,'amu*angstrom^2'), symmetry=1, barrier=(26.0778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13613,'amu*angstrom^2'), symmetry=1, barrier=(26.1218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13607,'amu*angstrom^2'), symmetry=1, barrier=(26.1204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13641,'amu*angstrom^2'), symmetry=1, barrier=(26.1283,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05788,0.0863548,-8.6551e-05,4.18202e-08,-7.5121e-12,33500.8,29.7203], Tmin=(100,'K'), Tmax=(1592.32,'K')), NASAPolynomial(coeffs=[23.1301,0.0103956,-6.79204e-07,-1.26246e-10,1.47226e-14,27724.5,-92.1748], Tmin=(1592.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CJCO)"""),
)

species(
    label = 'C=[C]C=CO[CH]C(28548)',
    structure = SMILES('C=[C]C=CO[CH]C'),
    E0 = (259.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10735,'amu*angstrom^2'), symmetry=1, barrier=(25.4601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10718,'amu*angstrom^2'), symmetry=1, barrier=(25.4564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10807,'amu*angstrom^2'), symmetry=1, barrier=(25.4768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10748,'amu*angstrom^2'), symmetry=1, barrier=(25.4632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28817,0.0907966,-9.47325e-05,4.69616e-08,-8.60356e-12,31385.4,29.3827], Tmin=(100,'K'), Tmax=(1566.16,'K')), NASAPolynomial(coeffs=[24.6208,0.00800332,4.82546e-07,-3.44973e-10,2.95007e-14,25308.4,-100.728], Tmin=(1566.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=[C]O[CH]C(28549)',
    structure = SMILES('C=CC=[C]O[CH]C'),
    E0 = (299.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.834817,'amu*angstrom^2'), symmetry=1, barrier=(19.1941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834806,'amu*angstrom^2'), symmetry=1, barrier=(19.1938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.83478,'amu*angstrom^2'), symmetry=1, barrier=(19.1932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.834829,'amu*angstrom^2'), symmetry=1, barrier=(19.1944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.775639,0.0832892,-8.26339e-05,3.94609e-08,-7.06423e-12,36264.6,31.0798], Tmin=(100,'K'), Tmax=(1562.51,'K')), NASAPolynomial(coeffs=[23.0805,0.0105484,-1.60053e-06,1.07142e-10,-2.90302e-15,30234,-90.101], Tmin=(1562.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=C[C]=[C]OCC(28550)',
    structure = SMILES('C=C[C]=[C]OCC'),
    E0 = (305.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.922854,'amu*angstrom^2'), symmetry=1, barrier=(21.2182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.927172,'amu*angstrom^2'), symmetry=1, barrier=(21.3175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925193,'amu*angstrom^2'), symmetry=1, barrier=(21.272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925911,'amu*angstrom^2'), symmetry=1, barrier=(21.2885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.381836,0.0704635,-5.11285e-05,7.47441e-09,5.00721e-12,36822.2,27.0312], Tmin=(100,'K'), Tmax=(938.04,'K')), NASAPolynomial(coeffs=[17.5144,0.0198591,-6.11125e-06,9.97122e-10,-6.70122e-14,32620.2,-59.793], Tmin=(938.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC=CO[CH]C(28551)',
    structure = SMILES('[CH]=CC=CO[CH]C'),
    E0 = (307.288,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.98068,'amu*angstrom^2'), symmetry=1, barrier=(22.5478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.980272,'amu*angstrom^2'), symmetry=1, barrier=(22.5384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.982974,'amu*angstrom^2'), symmetry=1, barrier=(22.6005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.980474,'amu*angstrom^2'), symmetry=1, barrier=(22.543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.182504,0.0749671,-3.8211e-05,-2.49321e-08,2.10651e-11,37124.7,26.209], Tmin=(100,'K'), Tmax=(920.31,'K')), NASAPolynomial(coeffs=[25.8719,0.00644485,5.84531e-07,-2.35814e-10,1.35247e-14,30435.2,-107.614], Tmin=(920.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCsJOC(O))"""),
)

species(
    label = 'C=C[CH][CH]OC=C(20237)',
    structure = SMILES('[CH2]C=C[CH]OC=C'),
    E0 = (150.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770445,0.0552698,-3.00636e-06,-4.21689e-08,2.18039e-11,18227.2,26.3997], Tmin=(100,'K'), Tmax=(962.91,'K')), NASAPolynomial(coeffs=[18.5287,0.0194874,-6.43996e-06,1.17768e-09,-8.70415e-14,13046.3,-67.7468], Tmin=(962.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C][C]=COCC(28552)',
    structure = SMILES('C=[C][C]=COCC'),
    E0 = (264.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.30683,'amu*angstrom^2'), symmetry=1, barrier=(30.0466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30642,'amu*angstrom^2'), symmetry=1, barrier=(30.0372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31143,'amu*angstrom^2'), symmetry=1, barrier=(30.1523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31024,'amu*angstrom^2'), symmetry=1, barrier=(30.1251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.691212,0.0844515,-8.53123e-05,4.28376e-08,-8.10335e-12,31967.6,27.3538], Tmin=(100,'K'), Tmax=(1482.8,'K')), NASAPolynomial(coeffs=[21.0646,0.0140147,-2.17375e-06,1.15217e-10,5.27942e-16,26807.2,-81.8205], Tmin=(1482.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[C]=COCC(28553)',
    structure = SMILES('[CH]C=C=COCC'),
    E0 = (289.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259588,0.0687834,-2.84632e-05,-1.8537e-08,1.38703e-11,34993.4,25.8833], Tmin=(100,'K'), Tmax=(964.283,'K')), NASAPolynomial(coeffs=[18.5528,0.0244015,-8.4266e-06,1.48856e-09,-1.04748e-13,30000.8,-69.2995], Tmin=(964.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C[CH][O](176)',
    structure = SMILES('C[CH][O]'),
    E0 = (157.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,1642.51],'cm^-1')),
        HinderedRotor(inertia=(0.123965,'amu*angstrom^2'), symmetry=1, barrier=(2.85019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65562,0.0114444,2.34936e-06,-4.83164e-09,1.17966e-12,18963.9,10.3625], Tmin=(100,'K'), Tmax=(1718.65,'K')), NASAPolynomial(coeffs=[6.06294,0.0136322,-6.35953e-06,1.18407e-09,-7.90642e-14,16985.9,-5.90233], Tmin=(1718.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]OC=C1[CH]C1(28554)',
    structure = SMILES('C[CH]OC=C1[CH]C1'),
    E0 = (294.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.450893,0.0568023,1.00881e-05,-7.01456e-08,3.5603e-11,35546.9,22.6063], Tmin=(100,'K'), Tmax=(931.729,'K')), NASAPolynomial(coeffs=[24.2024,0.00913151,-5.78877e-07,3.19281e-11,-8.95e-15,28764.1,-102.949], Tmin=(931.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.328,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC1C[CH][C]=CO1(28555)',
    structure = SMILES('CC1C[CH][C]=CO1'),
    E0 = (188.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73705,0.0261714,7.98961e-05,-1.30245e-07,5.45497e-11,22731.9,16.8889], Tmin=(100,'K'), Tmax=(923.852,'K')), NASAPolynomial(coeffs=[18.0664,0.0161795,-2.45025e-06,3.06788e-10,-2.64784e-14,17124,-74.6189], Tmin=(923.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=COC=C(20244)',
    structure = SMILES('C=CC=COC=C'),
    E0 = (21.5857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.715205,0.0541813,5.61043e-06,-5.80667e-08,2.96746e-11,2731.05,25.4649], Tmin=(100,'K'), Tmax=(934.322,'K')), NASAPolynomial(coeffs=[20.9328,0.0134954,-2.71129e-06,4.15908e-10,-3.34102e-14,-3049,-81.4141], Tmin=(934.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.5857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=COCC(28556)',
    structure = SMILES('C=C=C=COCC'),
    E0 = (95.6504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.302153,0.0684672,-3.77907e-05,-9.62715e-09,1.13087e-11,11648.9,23.6753], Tmin=(100,'K'), Tmax=(958.723,'K')), NASAPolynomial(coeffs=[19.8052,0.0174655,-5.50995e-06,9.66679e-10,-6.96574e-14,6513.56,-76.8726], Tmin=(958.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.6504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=CC1=COC1C(28557)',
    structure = SMILES('C=CC1=COC1C'),
    E0 = (12.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.855165,0.0373529,7.87942e-05,-1.50524e-07,6.64168e-11,1671.38,20.4954], Tmin=(100,'K'), Tmax=(925.057,'K')), NASAPolynomial(coeffs=[28.1834,0.00126434,4.2174e-06,-8.60194e-10,4.72516e-14,-6896.61,-128.201], Tmin=(925.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C[C]=C[O](26659)',
    structure = SMILES('C=C[C]=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61748,'amu*angstrom^2'), symmetry=1, barrier=(37.1891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98077,0.0322088,7.53766e-06,-4.54875e-08,2.38742e-11,27216,16.1107], Tmin=(100,'K'), Tmax=(899.949,'K')), NASAPolynomial(coeffs=[16.1069,0.00249504,1.93925e-06,-5.05324e-10,3.47549e-14,23334.1,-57.9915], Tmin=(899.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=[C]C1OC1C(28558)',
    structure = SMILES('[CH2]C=[C]C1OC1C'),
    E0 = (317.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13787,0.0497787,7.05009e-06,-5.41195e-08,2.8403e-11,38257.1,24.4293], Tmin=(100,'K'), Tmax=(888.336,'K')), NASAPolynomial(coeffs=[16.0429,0.0193677,-3.57528e-06,3.65267e-10,-1.96869e-14,34160.8,-53.8654], Tmin=(888.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C=CO[CH]C(28559)',
    structure = SMILES('CC#C[CH]O[CH]C'),
    E0 = (320.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2100,2250,500,550,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.22915,0.0808678,-9.51542e-05,6.1011e-08,-1.54027e-11,38728.9,27.5482], Tmin=(100,'K'), Tmax=(1020.08,'K')), NASAPolynomial(coeffs=[13.9236,0.0243202,-7.8142e-06,1.19336e-09,-7.17966e-14,36083.2,-38.0656], Tmin=(1020.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]O[C]=C=CC(28560)',
    structure = SMILES('C[CH]C#CO[CH]C'),
    E0 = (285.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2100,2250,500,550,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.639424,0.0813629,-8.06628e-05,3.97288e-08,-7.33595e-12,34560.1,28.6425], Tmin=(100,'K'), Tmax=(1540.8,'K')), NASAPolynomial(coeffs=[20.45,0.0134414,-1.71637e-06,2.21115e-11,6.84066e-15,29624.7,-77.1448], Tmin=(1540.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Sec_Propargyl) + radical(CCsJOCs)"""),
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
    label = '[CH2]C=C1[CH]OC1C(28561)',
    structure = SMILES('[CH2]C=C1[CH]OC1C'),
    E0 = (193.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76357,0.0311514,5.09741e-05,-8.4354e-08,3.288e-11,23388,24.451], Tmin=(100,'K'), Tmax=(994.452,'K')), NASAPolynomial(coeffs=[13.6184,0.027882,-1.10881e-05,2.16313e-09,-1.60464e-13,18834,-43.718], Tmin=(994.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(C=CCJ(O)C) + radical(Allyl_P)"""),
)

species(
    label = 'C[CH]OC1[C]=CC1(28562)',
    structure = SMILES('C[CH]OC1[C]=CC1'),
    E0 = (387.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925991,0.0550704,-1.48732e-05,-2.22658e-08,1.28585e-11,46763.4,25.9953], Tmin=(100,'K'), Tmax=(999.906,'K')), NASAPolynomial(coeffs=[15.9468,0.0229765,-8.72417e-06,1.63458e-09,-1.17892e-13,42360,-53.4679], Tmin=(999.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(CCsJOCs)"""),
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
    label = '[CH]=C=CO[CH]C(28563)',
    structure = SMILES('[CH]=C=CO[CH]C'),
    E0 = (303.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,610,2055,3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12498,'amu*angstrom^2'), symmetry=1, barrier=(25.8656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12466,'amu*angstrom^2'), symmetry=1, barrier=(25.8581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13087,'amu*angstrom^2'), symmetry=1, barrier=(26.001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399911,0.0747907,-8.15405e-05,4.15088e-08,-7.72286e-12,36677,25.4427], Tmin=(100,'K'), Tmax=(1564.34,'K')), NASAPolynomial(coeffs=[21.208,0.00361349,1.98028e-06,-5.92768e-10,4.55074e-14,31865.3,-82.2433], Tmin=(1564.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCsJOC(O)) + radical(C=C=CJ)"""),
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
    E0 = (259.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (275.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (438.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (518.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (487.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (308.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (435.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (473.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (496.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (458.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (529.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (452.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (416.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (434.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (609.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (396.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (285.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (298.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (267.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (267.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (569.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (323.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (535.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (520.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (322.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (384.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (387.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (284.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (685.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['CH3CHO(37)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['[CH2]C1[C]=COC1C(28544)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(187000,'s^-1'), n=1.48, Ea=(16.3176,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C[C]=COC=C(26641)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.67e+12,'cm^3/(mol*s)'), n=0.1, Ea=(6.4601,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2816 used for Cds-HH_Cds-OsH;HJ
Exact match found for rate rule [Cds-HH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C=C=CO[CH]C(28545)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC#CO[CH]C(28546)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH3CHO(37)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]COC=[C]C=C(28547)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.7e+13,'s^-1','+|-',2), n=-0.1, Ea=(158.364,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""From training reaction 347 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C=CO[CH]C(28548)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=CC=[C]O[CH]C(28549)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C]=[C]OCC(28550)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.4016e+09,'s^-1'), n=1.0875, Ea=(153.657,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;Cs_H_out_H/NonDeC] + [R3H_SS_O;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_O;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC=CO[CH]C(28551)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C=C[CH][CH]OC=C(20237)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C][C]=COCC(28552)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.54505e+10,'s^-1'), n=0.959062, Ea=(152.545,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C[C]=COCC(28553)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH][O](176)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C[CH]OC=C1[CH]C1(28554)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['CC1C[CH][C]=CO1(28555)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C=CC=COC=C(20244)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C=C=C=COCC(28556)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C=CC1=COC1C(28557)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CHCH3(T)(95)', 'C=C[C]=C[O](26659)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['[CH2]C=[C]C1OC1C(28558)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(64.7047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[C]=C=CO[CH]C(28559)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH]O[C]=C=CC(28560)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.992e-05,'s^-1'), n=4.805, Ea=(234.476,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R4H_MMS;Cd_rad_out_singleNd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C=CO[CH][C]=CC(25014)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R7HJ_5;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['[CH2]C=C1[CH]OC1C(28561)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.16207e+08,'s^-1'), n=0.911389, Ea=(125.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C[CH]OC1[C]=CC1(28562)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(128.613,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 125.0 to 128.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CO[CH]C(26656)'],
    products = ['C=COC=C=CC(25016)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(19)', '[CH]=C=CO[CH]C(28563)'],
    products = ['C=C[C]=CO[CH]C(26656)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4590',
    isomers = [
        'C=C[C]=CO[CH]C(26656)',
    ],
    reactants = [
        ('CH3CHO(37)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4590',
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

