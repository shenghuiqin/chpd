species(
    label = 'C#CC([CH2])C=[C]C=C(26593)',
    structure = SMILES('C#CC([CH2])C=[C]C=C'),
    E0 = (652.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14611,'amu*angstrom^2'), symmetry=1, barrier=(26.3512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99228,'amu*angstrom^2'), symmetry=1, barrier=(45.8063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13279,'amu*angstrom^2'), symmetry=1, barrier=(26.045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14769,'amu*angstrom^2'), symmetry=1, barrier=(26.3876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.106035,0.0879548,-9.62203e-05,5.49266e-08,-1.18564e-11,78637.3,28.4965], Tmin=(100,'K'), Tmax=(897.011,'K')), NASAPolynomial(coeffs=[15.6045,0.0269246,-9.25904e-06,1.51484e-09,-9.70545e-14,75455.7,-47.6144], Tmin=(897.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Isobutyl)"""),
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
    label = '[CH]=C1CC1C=[C]C=C(28250)',
    structure = SMILES('[CH]=C1CC1C=[C]C=C'),
    E0 = (756.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220129,0.0745531,-5.75751e-05,1.7652e-08,-1.00982e-13,91130.7,26.0166], Tmin=(100,'K'), Tmax=(1001.52,'K')), NASAPolynomial(coeffs=[16.3424,0.0262698,-9.38525e-06,1.63302e-09,-1.109e-13,87093.5,-55.8218], Tmin=(1001.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(756.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C(C=C)=CC1[CH2](27998)',
    structure = SMILES('[CH]=C1C(C=C)=CC1[CH2]'),
    E0 = (709.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.823671,0.0565134,-5.78462e-06,-3.63419e-08,1.92961e-11,85424.2,24.029], Tmin=(100,'K'), Tmax=(951.01,'K')), NASAPolynomial(coeffs=[15.9327,0.0255341,-8.2933e-06,1.42859e-09,-9.97137e-14,81077.6,-55.8478], Tmin=(951.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1C=[C]C([CH2])C1(28251)',
    structure = SMILES('C#CC1C=[C]C([CH2])C1'),
    E0 = (666.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813743,0.0576012,-1.14906e-05,-2.87829e-08,1.61396e-11,80237.6,24.634], Tmin=(100,'K'), Tmax=(961.369,'K')), NASAPolynomial(coeffs=[15.4956,0.0262228,-8.88541e-06,1.55465e-09,-1.08497e-13,76041.8,-52.7591], Tmin=(961.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(666.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
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
    label = 'C#CC(=C)C=[C]C=C(28252)',
    structure = SMILES('C#CC(=C)C=[C]C=C'),
    E0 = (574.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49199,'amu*angstrom^2'), symmetry=1, barrier=(34.3038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49484,'amu*angstrom^2'), symmetry=1, barrier=(34.3693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49236,'amu*angstrom^2'), symmetry=1, barrier=(34.3124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.144733,0.0727227,-4.5511e-05,-4.11023e-09,1.01744e-11,69207.2,25.6107], Tmin=(100,'K'), Tmax=(941.004,'K')), NASAPolynomial(coeffs=[20.0917,0.0176565,-5.1145e-06,8.38276e-10,-5.84141e-14,64137.1,-76.4016], Tmin=(941.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC([CH2])C#CC=C(28253)',
    structure = SMILES('C#CC([CH2])C#CC=C'),
    E0 = (662.024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2175,2250,500,525,550,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,358.709,358.709],'cm^-1')),
        HinderedRotor(inertia=(0.321591,'amu*angstrom^2'), symmetry=1, barrier=(29.364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321592,'amu*angstrom^2'), symmetry=1, barrier=(29.364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321592,'amu*angstrom^2'), symmetry=1, barrier=(29.364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.617069,'amu*angstrom^2'), symmetry=1, barrier=(56.3438,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719875,0.0614596,-2.93747e-05,-9.46421e-09,9.19196e-12,79750.7,29.0389], Tmin=(100,'K'), Tmax=(975.198,'K')), NASAPolynomial(coeffs=[15.7936,0.0237083,-8.34202e-06,1.47499e-09,-1.0276e-13,75665.9,-49.1748], Tmin=(975.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCtCsH) + group(Cs-CsHHH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC([CH2])C=C=C=C(28254)',
    structure = SMILES('C#CC([CH2])C=C=C=C'),
    E0 = (682.97,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,2175,525,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.29828,'amu*angstrom^2'), symmetry=1, barrier=(29.85,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29866,'amu*angstrom^2'), symmetry=1, barrier=(29.8588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29638,'amu*angstrom^2'), symmetry=1, barrier=(29.8063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.141,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.17262,0.085896,-0.000103093,6.75577e-08,-1.77554e-11,82278.7,26.9881], Tmin=(100,'K'), Tmax=(928.474,'K')), NASAPolynomial(coeffs=[13.2599,0.0295146,-1.20059e-05,2.15566e-09,-1.4541e-13,79848.4,-35.1792], Tmin=(928.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(682.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(Isobutyl)"""),
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
    label = 'C=C[C]=CC=C(26601)',
    structure = SMILES('C=C[C]=CC=C'),
    E0 = (344.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32866,'amu*angstrom^2'), symmetry=1, barrier=(30.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33122,'amu*angstrom^2'), symmetry=1, barrier=(30.6074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31057,0.0479268,-9.42475e-06,-2.91657e-08,1.70739e-11,41563.7,19.8807], Tmin=(100,'K'), Tmax=(926.369,'K')), NASAPolynomial(coeffs=[15.4385,0.0157821,-4.10458e-06,6.34548e-10,-4.37525e-14,37707.8,-53.8816], Tmin=(926.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C#C[C](C)C=[C]C=C(28255)',
    structure = SMILES('C#CC(C)=C[C]=C[CH2]'),
    E0 = (563.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.208575,0.0753162,-6.61471e-05,3.14929e-08,-6.03918e-12,67862.4,27.6681], Tmin=(100,'K'), Tmax=(1258.67,'K')), NASAPolynomial(coeffs=[15.1055,0.0279748,-9.72909e-06,1.61085e-09,-1.03979e-13,64112.3,-47.6282], Tmin=(1258.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC([CH2])[C]=CC=C(28256)',
    structure = SMILES('C#CC([CH2])[C]=CC=C'),
    E0 = (691.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.960941,'amu*angstrom^2'), symmetry=1, barrier=(22.0939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960605,'amu*angstrom^2'), symmetry=1, barrier=(22.0862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960243,'amu*angstrom^2'), symmetry=1, barrier=(22.0779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.960548,'amu*angstrom^2'), symmetry=1, barrier=(22.0849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0992601,0.0864987,-9.37903e-05,5.39913e-08,-1.23006e-11,83310.2,29.3466], Tmin=(100,'K'), Tmax=(1074.36,'K')), NASAPolynomial(coeffs=[16.1501,0.0259987,-9.31954e-06,1.57406e-09,-1.03023e-13,79818.7,-50.2124], Tmin=(1074.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC([CH2])C=C[C]=C(28257)',
    structure = SMILES('C#CC([CH2])C=C[C]=C'),
    E0 = (652.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14611,'amu*angstrom^2'), symmetry=1, barrier=(26.3512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99228,'amu*angstrom^2'), symmetry=1, barrier=(45.8063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13279,'amu*angstrom^2'), symmetry=1, barrier=(26.045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14769,'amu*angstrom^2'), symmetry=1, barrier=(26.3876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.106035,0.0879548,-9.62203e-05,5.49266e-08,-1.18564e-11,78637.3,28.4965], Tmin=(100,'K'), Tmax=(897.011,'K')), NASAPolynomial(coeffs=[15.6045,0.0269246,-9.25904e-06,1.51484e-09,-9.70545e-14,75455.7,-47.6144], Tmin=(897.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(C)[C]=[C]C=C(28258)',
    structure = SMILES('C#CC(C)C#C[CH][CH2]'),
    E0 = (687.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2175,2250,500,525,550,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,372.177,372.179],'cm^-1')),
        HinderedRotor(inertia=(0.0103053,'amu*angstrom^2'), symmetry=1, barrier=(21.8361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00121701,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754915,'amu*angstrom^2'), symmetry=1, barrier=(74.2046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754893,'amu*angstrom^2'), symmetry=1, barrier=(74.2043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75492,'amu*angstrom^2'), symmetry=1, barrier=(74.2046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632692,0.0674593,-4.68911e-05,1.08152e-08,1.90878e-12,82787,30.9194], Tmin=(100,'K'), Tmax=(949.495,'K')), NASAPolynomial(coeffs=[13.3964,0.0288588,-9.87581e-06,1.6524e-09,-1.09097e-13,79679.3,-33.5975], Tmin=(949.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCtCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[C]([CH2])C=CC=C(28259)',
    structure = SMILES('C#CC([CH2])=CC=C[CH2]'),
    E0 = (482.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76703,0.0575335,-5.77177e-06,-3.75025e-08,1.99709e-11,58113.1,27.6237], Tmin=(100,'K'), Tmax=(947.376,'K')), NASAPolynomial(coeffs=[16.1705,0.0258214,-8.32445e-06,1.42307e-09,-9.89511e-14,53699,-53.7492], Tmin=(947.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC([CH2])C#C(28260)',
    structure = SMILES('[CH]=CC=CC([CH2])C#C'),
    E0 = (700.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(1.02482,'amu*angstrom^2'), symmetry=1, barrier=(23.5626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02701,'amu*angstrom^2'), symmetry=1, barrier=(23.6129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02931,'amu*angstrom^2'), symmetry=1, barrier=(23.6658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02972,'amu*angstrom^2'), symmetry=1, barrier=(23.6753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.172023,0.0858053,-8.70358e-05,4.39482e-08,-8.06811e-12,84428,29.5367], Tmin=(100,'K'), Tmax=(971.704,'K')), NASAPolynomial(coeffs=[17.5881,0.0236506,-7.99884e-06,1.32412e-09,-8.66463e-14,80459.3,-58.2967], Tmin=(971.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[C]#CC(C)C=[C]C=C(28261)',
    structure = SMILES('[C]#CC(C)C=[C]C=C'),
    E0 = (784.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.942443,'amu*angstrom^2'), symmetry=1, barrier=(21.6686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09375,'amu*angstrom^2'), symmetry=1, barrier=(48.1394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.942257,'amu*angstrom^2'), symmetry=1, barrier=(21.6643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.94264,'amu*angstrom^2'), symmetry=1, barrier=(21.6731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0642801,0.0893838,-0.000104852,6.75939e-08,-1.74027e-11,94516.9,27.3072], Tmin=(100,'K'), Tmax=(950.391,'K')), NASAPolynomial(coeffs=[14.0308,0.0300592,-1.12177e-05,1.9116e-09,-1.24681e-13,91837.8,-39.9761], Tmin=(950.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(784.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=CJC=C)"""),
)

species(
    label = '[C]#CC([CH2])C=CC=C(28262)',
    structure = SMILES('[C]#CC([CH2])C=CC=C'),
    E0 = (790.729,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.742003,'amu*angstrom^2'), symmetry=1, barrier=(17.0601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.740971,'amu*angstrom^2'), symmetry=1, barrier=(17.0364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.00785,'amu*angstrom^2'), symmetry=1, barrier=(69.1564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.741514,'amu*angstrom^2'), symmetry=1, barrier=(17.0489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.000881549,0.0848936,-8.95418e-05,4.88138e-08,-9.87116e-12,95249.5,29.2576], Tmin=(100,'K'), Tmax=(906.925,'K')), NASAPolynomial(coeffs=[15.4804,0.0263954,-8.95591e-06,1.46008e-09,-9.36341e-14,92039.8,-46.1262], Tmin=(906.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(790.729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Acetyl)"""),
)

species(
    label = 'C#CC(C)C=[C][C]=C(28263)',
    structure = SMILES('C#CC(C)[CH]C#C[CH2]'),
    E0 = (623.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.224048,0.0745473,-6.89553e-05,3.56198e-08,-7.33862e-12,75125.4,30.5557], Tmin=(100,'K'), Tmax=(1257.1,'K')), NASAPolynomial(coeffs=[14.2243,0.0266921,-7.90671e-06,1.15141e-09,-6.7638e-14,71866.8,-39.1512], Tmin=(1257.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC(C)C#C(28264)',
    structure = SMILES('[CH]C=C=CC(C)C#C'),
    E0 = (671.968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07179,'amu*angstrom^2'), symmetry=1, barrier=(47.6346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07502,'amu*angstrom^2'), symmetry=1, barrier=(47.7087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06979,'amu*angstrom^2'), symmetry=1, barrier=(47.5884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08093,'amu*angstrom^2'), symmetry=1, barrier=(47.8448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.276619,0.0826163,-7.4689e-05,3.7702e-08,-7.98389e-12,80952.6,26.8864], Tmin=(100,'K'), Tmax=(1111.74,'K')), NASAPolynomial(coeffs=[12.139,0.039936,-1.71033e-05,3.17023e-09,-2.18676e-13,78315,-31.5993], Tmin=(1111.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C#CC([CH2])C=C1[CH]C1(28265)',
    structure = SMILES('C#CC([CH2])C=C1[CH]C1'),
    E0 = (687.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.475051,0.0674623,-3.80254e-05,-2.33268e-09,6.99844e-12,82849.6,25.886], Tmin=(100,'K'), Tmax=(972.213,'K')), NASAPolynomial(coeffs=[15.9175,0.0263424,-9.16657e-06,1.59312e-09,-1.0924e-13,78787.6,-53.6276], Tmin=(972.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=C[C]=CC1[C]=CC1(28026)',
    structure = SMILES('C=C[C]=CC1[C]=CC1'),
    E0 = (719.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613599,0.0646829,-3.48589e-05,-1.6223e-09,5.61018e-12,86657,25.7994], Tmin=(100,'K'), Tmax=(1007.12,'K')), NASAPolynomial(coeffs=[14.7615,0.0284016,-1.04759e-05,1.86714e-09,-1.2882e-13,82797.6,-47.5694], Tmin=(1007.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1[C]=CC(C=C)=C1(28033)',
    structure = SMILES('[CH2]C1[C]=CC(C=C)=C1'),
    E0 = (635.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18838,0.0506356,-2.03277e-06,-2.97609e-08,1.4152e-11,76545.3,23.4251], Tmin=(100,'K'), Tmax=(1000.06,'K')), NASAPolynomial(coeffs=[12.1738,0.0319869,-1.19945e-05,2.16712e-09,-1.50999e-13,73083.5,-35.8966], Tmin=(1000.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(635.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentadiene) + radical(1,3-cyclopentadiene-vinyl-1) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1C=[C][CH]CC1(28266)',
    structure = SMILES('C#CC1C=[C][CH]CC1'),
    E0 = (568.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08751,0.0487853,1.37432e-05,-5.07551e-08,2.21559e-11,68484.4,20.0027], Tmin=(100,'K'), Tmax=(997.08,'K')), NASAPolynomial(coeffs=[14.8019,0.0296376,-1.14145e-05,2.14666e-09,-1.5486e-13,63966.5,-55.0624], Tmin=(997.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(568.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclohexene) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C#CC(=C)C=CC=C(28267)',
    structure = SMILES('C#CC(=C)C=CC=C'),
    E0 = (375.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.220185,0.0662115,-1.65405e-05,-3.68123e-08,2.18604e-11,45275.2,25.743], Tmin=(100,'K'), Tmax=(951.27,'K')), NASAPolynomial(coeffs=[21.6267,0.017578,-5.10077e-06,8.97096e-10,-6.71217e-14,39330.3,-86.302], Tmin=(951.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)Ct) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(C)C=C=C=C(28268)',
    structure = SMILES('C#CC(C)C=C=C=C'),
    E0 = (477.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301294,0.0824802,-8.45299e-05,4.71692e-08,-1.07587e-11,57609,24.745], Tmin=(100,'K'), Tmax=(1051.93,'K')), NASAPolynomial(coeffs=[13.2467,0.0332551,-1.43379e-05,2.68479e-09,-1.86655e-13,54885.4,-38.3645], Tmin=(1051.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + group(Ct-CtH)"""),
)

species(
    label = 'C#C[CH]CC=[C]C=C(26590)',
    structure = SMILES('C#C[CH]CC=[C]C=C'),
    E0 = (626.781,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(3.10668,'amu*angstrom^2'), symmetry=1, barrier=(71.4287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.843074,'amu*angstrom^2'), symmetry=1, barrier=(19.3839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09628,'amu*angstrom^2'), symmetry=1, barrier=(71.1895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.09656,'amu*angstrom^2'), symmetry=1, barrier=(71.196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26904,0.0773468,-6.6166e-05,2.27803e-08,3.35569e-13,75523.2,28.3503], Tmin=(100,'K'), Tmax=(870.628,'K')), NASAPolynomial(coeffs=[15.2039,0.0261562,-7.9932e-06,1.22521e-09,-7.63792e-14,72262.2,-45.4257], Tmin=(870.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC[CH]C=[C]C=C(28269)',
    structure = SMILES('C#CCC=C[C]=C[CH2]'),
    E0 = (590.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483784,0.069306,-4.6557e-05,8.99298e-09,2.52479e-12,71205.9,28.2061], Tmin=(100,'K'), Tmax=(977.511,'K')), NASAPolynomial(coeffs=[14.4963,0.0287684,-1.01348e-05,1.73708e-09,-1.16687e-13,67663.6,-43.1838], Tmin=(977.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC1C=C(C=C)C1(28270)',
    structure = SMILES('C#CC1C=C(C=C)C1'),
    E0 = (407.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476308,0.0607852,-5.91824e-06,-4.41525e-08,2.36362e-11,49204.3,22.3493], Tmin=(100,'K'), Tmax=(952.417,'K')), NASAPolynomial(coeffs=[19.9484,0.0199387,-6.05426e-06,1.06787e-09,-7.85786e-14,43638.7,-80.3891], Tmin=(952.417,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene)"""),
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
    label = '[CH]=C=CC=[C]C=C(28029)',
    structure = SMILES('C#CC=C[C]=C[CH2]'),
    E0 = (603.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2175,525,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.8259,'amu*angstrom^2'), symmetry=1, barrier=(41.981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82498,'amu*angstrom^2'), symmetry=1, barrier=(41.96,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82933,'amu*angstrom^2'), symmetry=1, barrier=(42.0598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06706,0.0544769,-2.32027e-05,-1.48007e-08,1.16032e-11,72699.6,23.0123], Tmin=(100,'K'), Tmax=(939.868,'K')), NASAPolynomial(coeffs=[15.2958,0.018748,-5.80419e-06,9.64012e-10,-6.61249e-14,68928.4,-50.5841], Tmin=(939.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC1CC1[C]=C[CH2](28271)',
    structure = SMILES('C#CC1CC1[C]=C[CH2]'),
    E0 = (704.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862538,0.0527609,9.61622e-06,-5.6315e-08,2.7414e-11,84900.2,27.2849], Tmin=(100,'K'), Tmax=(940.505,'K')), NASAPolynomial(coeffs=[17.6436,0.0221038,-6.42227e-06,1.08082e-09,-7.72796e-14,79943,-62.2176], Tmin=(940.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(704.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclopropane) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C1CC=C=CC1[CH2](28272)',
    structure = SMILES('[CH]=C1CC=C=CC1[CH2]'),
    E0 = (650.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21752,0.045216,2.51658e-05,-5.92771e-08,2.39264e-11,78312,21.9629], Tmin=(100,'K'), Tmax=(1018.93,'K')), NASAPolynomial(coeffs=[13.6758,0.0340514,-1.39604e-05,2.67539e-09,-1.9332e-13,73813.9,-47.9892], Tmin=(1018.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(650.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC([CH2])C=C=[C]C(28273)',
    structure = SMILES('C#CC([CH2])[CH]C#CC'),
    E0 = (672.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2175,2250,500,525,550,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,221.251,221.251],'cm^-1')),
        HinderedRotor(inertia=(0.212715,'amu*angstrom^2'), symmetry=1, barrier=(7.38917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93817,'amu*angstrom^2'), symmetry=1, barrier=(67.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0444992,'amu*angstrom^2'), symmetry=1, barrier=(67.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0444992,'amu*angstrom^2'), symmetry=1, barrier=(67.3271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93816,'amu*angstrom^2'), symmetry=1, barrier=(67.3271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.132022,0.077551,-7.96946e-05,4.62698e-08,-1.04019e-11,80982.3,31.6016], Tmin=(100,'K'), Tmax=(1257.58,'K')), NASAPolynomial(coeffs=[13.0445,0.0261105,-5.96945e-06,6.29822e-10,-2.54908e-14,78554.5,-30.3928], Tmin=(1257.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#CC([CH2])[C]=C=CC(28274)',
    structure = SMILES('C#CC([CH2])C#C[CH]C'),
    E0 = (687.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2100,2175,2250,500,525,550,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,328.69,328.781],'cm^-1')),
        HinderedRotor(inertia=(0.00156087,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984979,'amu*angstrom^2'), symmetry=1, barrier=(75.549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985073,'amu*angstrom^2'), symmetry=1, barrier=(75.5502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.985383,'amu*angstrom^2'), symmetry=1, barrier=(75.5403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0285315,'amu*angstrom^2'), symmetry=1, barrier=(75.5445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.212485,0.0731479,-6.55318e-05,3.25834e-08,-6.43833e-12,82785.1,31.7281], Tmin=(100,'K'), Tmax=(1330.25,'K')), NASAPolynomial(coeffs=[14.5194,0.0258708,-7.42148e-06,1.05526e-09,-6.09826e-14,79355.4,-39.9613], Tmin=(1330.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCtCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C#C[C]([CH2])C=C=CC(28275)',
    structure = SMILES('C#CC([CH2])=C[C]=CC'),
    E0 = (563.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.208575,0.0753162,-6.61471e-05,3.14929e-08,-6.03918e-12,67862.4,27.6681], Tmin=(100,'K'), Tmax=(1258.67,'K')), NASAPolynomial(coeffs=[15.1055,0.0279748,-9.72909e-06,1.61085e-09,-1.03979e-13,64112.3,-47.6282], Tmin=(1258.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[C]#CC([CH2])C=C=CC(28276)',
    structure = SMILES('[C]#CC([CH2])C=C=CC'),
    E0 = (843.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.578872,'amu*angstrom^2'), symmetry=1, barrier=(13.3094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577836,'amu*angstrom^2'), symmetry=1, barrier=(13.2856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.57764,'amu*angstrom^2'), symmetry=1, barrier=(13.2811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7647,'amu*angstrom^2'), symmetry=1, barrier=(63.566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0567248,0.0910447,-0.000121373,9.27793e-08,-2.84045e-11,101589,29.0878], Tmin=(100,'K'), Tmax=(880.205,'K')), NASAPolynomial(coeffs=[10.7664,0.0352847,-1.42654e-05,2.50363e-09,-1.64424e-13,99978,-19.6532], Tmin=(880.205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(843.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Acetyl)"""),
)

species(
    label = 'C#CC1[CH]C(=C[CH2])C1(28277)',
    structure = SMILES('C#CC1[CH]C(=C[CH2])C1'),
    E0 = (594.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.1727,-0.00739419,0.000126318,-1.24262e-07,2.88541e-11,71123.9,-16.2095], Tmin=(100,'K'), Tmax=(1708.71,'K')), NASAPolynomial(coeffs=[74.3708,0.0302863,-7.17662e-05,1.74007e-08,-1.29163e-12,21744.8,-440.616], Tmin=(1708.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(594.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methylenecyclobutane) + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = 'C#CC([CH2])C1[C]=CC1(28278)',
    structure = SMILES('C#CC([CH2])C1[C]=CC1'),
    E0 = (783.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936154,0.0543856,-1.65318e-06,-4.16791e-08,2.2002e-11,94380,28.6296], Tmin=(100,'K'), Tmax=(925.668,'K')), NASAPolynomial(coeffs=[15.6214,0.0243033,-6.98992e-06,1.1155e-09,-7.55166e-14,90231.4,-48.8078], Tmin=(925.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(783.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[C]=CCC=C=C1(28279)',
    structure = SMILES('[CH2]C1[C]=CCC=C=C1'),
    E0 = (724.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15142,0.054436,-2.16257e-05,-2.82243e-09,3.0322e-12,87304.6,21.8676], Tmin=(100,'K'), Tmax=(1170.98,'K')), NASAPolynomial(coeffs=[10.5941,0.0353934,-1.4158e-05,2.56229e-09,-1.74734e-13,84187.3,-29.0462], Tmin=(1170.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(724.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(=C)C=C=CC(28280)',
    structure = SMILES('C#CC(=C)C=C=CC'),
    E0 = (427.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.149,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.17629,0.0735085,-5.21729e-05,1.16511e-08,1.65942e-12,51618.6,25.9324], Tmin=(100,'K'), Tmax=(1041.26,'K')), NASAPolynomial(coeffs=[17.6094,0.0252962,-9.74045e-06,1.78338e-09,-1.24925e-13,46971.3,-63.7598], Tmin=(1041.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Ct) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C=CC([CH2])C#C(28281)',
    structure = SMILES('[CH]=C=CC([CH2])C#C'),
    E0 = (696.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(1.85789,'amu*angstrom^2'), symmetry=1, barrier=(42.7166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11824,'amu*angstrom^2'), symmetry=1, barrier=(25.7106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45542,'amu*angstrom^2'), symmetry=1, barrier=(33.463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56098,0.074528,-9.1924e-05,6.08642e-08,-1.57901e-11,83938.6,25.3523], Tmin=(100,'K'), Tmax=(997.796,'K')), NASAPolynomial(coeffs=[13.2156,0.0208572,-6.81938e-06,1.04885e-09,-6.3214e-14,81559.7,-34.9374], Tmin=(997.796,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(C=C=CJ) + radical(Isobutyl)"""),
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
    E0 = (652.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (756.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (709.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (772.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (803.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (886.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (911.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (743.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (915.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (751.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (887.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (867.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (839.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (845.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (923.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (880.683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (864.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (840.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (804.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (903.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (790.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (783.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (708.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (740.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (741.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (677.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (812.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (812.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (660.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (985.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (704.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (840.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (886.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (921.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (728.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (924.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (778.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (783.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (759.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (677.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (1078.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['CH2CHCCH(26391)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['[CH]=C1CC1C=[C]C=C(28250)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(103.929,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 102.0 to 103.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['[CH]=C1C(C=C)=CC1[CH2](27998)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(56.6251,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 53.6 to 56.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC1C=[C]C([CH2])C1(28251)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C#CC(=C)C=[C]C=C(28252)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(454.257,'m^3/(mol*s)'), n=1.51715, Ea=(17.7318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-TwoDe_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#CC([CH2])C#CC=C(28253)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C#CC([CH2])C=C=C=C(28254)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C(4699)', 'CH2CHCCH(26391)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00294841,'m^3/(mol*s)'), n=2.48333, Ea=(17.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H(33)', 'C=C[C]=CC=C(26601)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00505363,'m^3/(mol*s)'), n=2.41967, Ea=(13.2813,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDeH_Cds;CJ] for rate rule [Cds-OneDeH_Cds;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#C[C](C)C=[C]C=C(28255)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC([CH2])[C]=CC=C(28256)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([CH2])C=C[C]=C(28257)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC(C)[C]=[C]C=C(28258)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#C[C]([CH2])C=CC=C(28259)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC=CC([CH2])C#C(28260)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]#CC(C)C=[C]C=C(28261)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]#CC([CH2])C=CC=C(28262)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(272058,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_single] for rate rule [R5H_TSSD;Ct_rad_out;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC(C)C=[C][C]=C(28263)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.11562e+07,'s^-1'), n=1.49763, Ea=(188.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_doubleC] for rate rule [R5HJ_3;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C[C]=CC(C)C#C(28264)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]C=C(4699)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC([CH2])C=C1[CH]C1(28265)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C=C[C]=CC1[C]=CC1(28026)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['[CH2]C1[C]=CC(C=C)=C1(28033)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.72e+09,'s^-1'), n=0.186, Ea=(55.7727,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 831 used for R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleDe
Exact match found for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC1C=[C][CH]CC1(28266)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC(=C)C=CC=C(28267)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC(C)C=C=C=C(28268)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#C[CH]CC=[C]C=C(26590)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC[CH]C=[C]C=C(28269)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC1C=C(C=C)C1(28270)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(19)', '[CH]=C=CC=[C]C=C(28029)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC1CC1[C]=C[CH2](28271)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(52.2574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 47.6 to 52.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['[CH]=C1CC=C=CC1[CH2](28272)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7_SMMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_SMMS;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([CH2])C=C=[C]C(28273)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#CC([CH2])[C]=C=CC(28274)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.992e-05,'s^-1'), n=4.805, Ea=(234.476,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R4H_MMS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#C[C]([CH2])C=C=CC(28275)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.21176e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#CC([CH2])C=C=CC(28276)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(244756,'s^-1'), n=1.235, Ea=(80.8557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;Cs_H_out_2H] for rate rule [R7H;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC1[CH]C(=C[CH2])C1(28277)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC([CH2])C1[C]=CC1(28278)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(131.123,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 126.6 to 131.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['[CH2]C1[C]=CCC=C=C1(28279)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.8036e+08,'s^-1'), n=0.568448, Ea=(107.224,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C#CC([CH2])C=[C]C=C(26593)'],
    products = ['C#CC(=C)C=C=CC(28280)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(19)', '[CH]=C=CC([CH2])C#C(28281)'],
    products = ['C#CC([CH2])C=[C]C=C(26593)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4565',
    isomers = [
        'C#CC([CH2])C=[C]C=C(26593)',
    ],
    reactants = [
        ('CH2CHCCH(26391)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4565',
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

