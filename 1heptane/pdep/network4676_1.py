species(
    label = 'C=C[C]=CC([O])C=O(29458)',
    structure = SMILES('C=C[C]=CC([O])C=O'),
    E0 = (224.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1685,370,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,343.134,346.746,348.461],'cm^-1')),
        HinderedRotor(inertia=(0.247155,'amu*angstrom^2'), symmetry=1, barrier=(20.4209,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155577,'amu*angstrom^2'), symmetry=1, barrier=(13.2131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158668,'amu*angstrom^2'), symmetry=1, barrier=(13.1681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.63554,0.0678828,-6.26057e-05,3.01866e-08,-5.80801e-12,27086.8,31.3455], Tmin=(100,'K'), Tmax=(1256.7,'K')), NASAPolynomial(coeffs=[14.8715,0.0225707,-8.52123e-06,1.49538e-09,-1.00369e-13,23508.8,-40.5875], Tmin=(1256.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'OCHCHO(48)',
    structure = SMILES('O=CC=O'),
    E0 = (-225.3,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,857.339,857.946,861.419,1717.87],'cm^-1')),
        HinderedRotor(inertia=(0.00964493,'amu*angstrom^2'), symmetry=1, barrier=(5.0669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-225.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OCHCHO""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C=C[C]=CC1OC1[O](29560)',
    structure = SMILES('C=C[C]=CC1OC1[O]'),
    E0 = (276.477,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.186822,0.0697532,-6.45394e-05,3.19724e-08,-6.09572e-12,33401.8,29.6301], Tmin=(100,'K'), Tmax=(1469.58,'K')), NASAPolynomial(coeffs=[15.1971,0.0198275,-4.32259e-06,4.55687e-10,-1.9763e-14,29969.5,-45.232], Tmin=(1469.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC([O])C1[O](29561)',
    structure = SMILES('C=CC1=CC([O])C1[O]'),
    E0 = (309.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862402,0.0463097,3.02595e-05,-8.55523e-08,3.93718e-11,37391.9,26.7559], Tmin=(100,'K'), Tmax=(947.238,'K')), NASAPolynomial(coeffs=[23.0421,0.00955602,-1.6541e-06,3.31538e-10,-3.40529e-14,30637,-92.5223], Tmin=(947.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1[C]=CC(C=O)O1(29535)',
    structure = SMILES('[CH2]C1[C]=CC(C=O)O1'),
    E0 = (203.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99728,0.050739,-1.29748e-06,-3.89668e-08,1.92644e-11,24552,27.5685], Tmin=(100,'K'), Tmax=(990.56,'K')), NASAPolynomial(coeffs=[17.6907,0.0190095,-7.28041e-06,1.42366e-09,-1.06986e-13,19494.3,-61.6451], Tmin=(990.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(25dihydrofuran) + radical(Cds_S) + radical(CJC(C)OC)"""),
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
    label = 'C=C[C]=CC(=O)C=O(29562)',
    structure = SMILES('C=C[C]=CC(=O)C=O'),
    E0 = (68.2075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1685,370,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13285,'amu*angstrom^2'), symmetry=1, barrier=(26.0464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.87887,'amu*angstrom^2'), symmetry=1, barrier=(43.199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1325,'amu*angstrom^2'), symmetry=1, barrier=(26.0385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04847,0.0677359,-7.13244e-05,4.18421e-08,-1.01799e-11,8307.41,25.2107], Tmin=(100,'K'), Tmax=(981.822,'K')), NASAPolynomial(coeffs=[10.3416,0.0298758,-1.34841e-05,2.56896e-09,-1.80029e-13,6482.52,-19.4532], Tmin=(981.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.2075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC#CC([O])C=O(29563)',
    structure = SMILES('C=CC#CC([O])C=O'),
    E0 = (195.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,2100,2250,500,550,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,486.677,486.892,486.969],'cm^-1')),
        HinderedRotor(inertia=(0.0983721,'amu*angstrom^2'), symmetry=1, barrier=(16.5835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0988586,'amu*angstrom^2'), symmetry=1, barrier=(16.5923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238985,'amu*angstrom^2'), symmetry=1, barrier=(40.2057,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16767,0.0526349,-2.4783e-05,-8.59566e-09,7.91472e-12,23659.6,30.3461], Tmin=(100,'K'), Tmax=(988.345,'K')), NASAPolynomial(coeffs=[14.5286,0.0198098,-7.21367e-06,1.3061e-09,-9.22322e-14,19980.8,-39.2061], Tmin=(988.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C=C=CC([O])C=O(29564)',
    structure = SMILES('C=C=C=CC([O])C=O'),
    E0 = (254.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,540,563.333,586.667,610,1970,2140,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.646151,'amu*angstrom^2'), symmetry=1, barrier=(14.8563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.645158,'amu*angstrom^2'), symmetry=1, barrier=(14.8335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18306,0.0626998,-5.87523e-05,2.9137e-08,-5.93912e-12,30716.4,28.8695], Tmin=(100,'K'), Tmax=(1160.97,'K')), NASAPolynomial(coeffs=[11.6756,0.0265484,-1.20433e-05,2.31496e-09,-1.63279e-13,28280.1,-23.3172], Tmin=(1160.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
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
    label = 'HCO(16)',
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
    label = 'C=C[C]=CC=O(26908)',
    structure = SMILES('C=C[C]=CC=O'),
    E0 = (169.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.18163,'amu*angstrom^2'), symmetry=1, barrier=(27.168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17982,'amu*angstrom^2'), symmetry=1, barrier=(27.1265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83685,0.0478006,-4.01103e-05,1.75891e-08,-3.19807e-12,20455.1,18.4329], Tmin=(100,'K'), Tmax=(1277.85,'K')), NASAPolynomial(coeffs=[10.001,0.0222446,-1.01115e-05,1.93838e-09,-1.36142e-13,18368.6,-22.9561], Tmin=(1277.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=C[O](2537)',
    structure = SMILES('[O]C=C[O]'),
    E0 = (-18.8461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,714.947,715.113],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7018,-0.0027881,5.18222e-05,-5.96735e-08,2.03775e-11,-2247.83,13.3237], Tmin=(100,'K'), Tmax=(1035.3,'K')), NASAPolynomial(coeffs=[6.65087,0.0113883,-5.76543e-06,1.26602e-09,-9.87778e-14,-4228.83,-7.62441], Tmin=(1035.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C=C[C]=C[C](O)C=O(29565)',
    structure = SMILES('[CH2]C=[C]C=C(O)C=O'),
    E0 = (44.9017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.282968,0.0822336,-8.51045e-05,4.33983e-08,-8.52139e-12,5564.8,28.0766], Tmin=(100,'K'), Tmax=(1292.4,'K')), NASAPolynomial(coeffs=[20.9068,0.0147511,-4.5769e-06,7.21708e-10,-4.6059e-14,246.349,-78.9729], Tmin=(1292.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.9017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=[C]C([O])C=O(29566)',
    structure = SMILES('C=CC=[C]C([O])C=O'),
    E0 = (263.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1685,370,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,376.661,376.677,376.71],'cm^-1')),
        HinderedRotor(inertia=(0.113578,'amu*angstrom^2'), symmetry=1, barrier=(11.4387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113567,'amu*angstrom^2'), symmetry=1, barrier=(11.4384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113619,'amu*angstrom^2'), symmetry=1, barrier=(11.4379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.664882,0.0660674,-5.84701e-05,2.63771e-08,-4.73037e-12,31758.8,32.1207], Tmin=(100,'K'), Tmax=(1343.84,'K')), NASAPolynomial(coeffs=[15.7588,0.0211391,-8.32026e-06,1.49793e-09,-1.0195e-13,27702.1,-45.1594], Tmin=(1343.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=OCOJ)"""),
)

species(
    label = 'C=[C]C=CC([O])C=O(29567)',
    structure = SMILES('C=C=C[CH]C([O])C=O'),
    E0 = (208.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410395,0.0691934,-5.91236e-05,2.50845e-08,-4.20575e-12,25183.5,30.9039], Tmin=(100,'K'), Tmax=(1435.83,'K')), NASAPolynomial(coeffs=[17.8413,0.0206334,-8.39302e-06,1.52972e-09,-1.04478e-13,20178,-59.4954], Tmin=(1435.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C[C]=[C]C(O)C=O(29568)',
    structure = SMILES('[CH2][CH]C#CC(O)C=O'),
    E0 = (182.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.545525,0.0669439,-6.4013e-05,3.29981e-08,-6.67045e-12,22069.2,35.1458], Tmin=(100,'K'), Tmax=(1309.31,'K')), NASAPolynomial(coeffs=[14.4526,0.0202051,-5.59569e-06,7.73157e-10,-4.3801e-14,18792,-34.3036], Tmin=(1309.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=C[C]=CC(O)[C]=O(29569)',
    structure = SMILES('C=C[C]=CC(O)[C]=O'),
    E0 = (140.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.255564,0.0745709,-7.67392e-05,4.02354e-08,-8.23121e-12,17026.5,32.7465], Tmin=(100,'K'), Tmax=(1201.37,'K')), NASAPolynomial(coeffs=[17.2675,0.0179287,-6.01657e-06,9.89605e-10,-6.42542e-14,12939,-52.4472], Tmin=(1201.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(140.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC=C[C]([O])C=O(29570)',
    structure = SMILES('C=CC=CC([O])=C[O]'),
    E0 = (2.27503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0667775,0.069152,-2.78199e-05,-3.49737e-08,2.49131e-11,431.467,26.7628], Tmin=(100,'K'), Tmax=(910.876,'K')), NASAPolynomial(coeffs=[25.3998,0.00451994,1.85071e-06,-5.068e-10,3.33175e-14,-6117.37,-103.704], Tmin=(910.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.27503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=CC=CC([O])C=O(29571)',
    structure = SMILES('[CH]=CC=CC([O])C=O'),
    E0 = (272.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,313.764,316.363],'cm^-1')),
        HinderedRotor(inertia=(0.180933,'amu*angstrom^2'), symmetry=1, barrier=(12.9523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184866,'amu*angstrom^2'), symmetry=1, barrier=(12.9513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182773,'amu*angstrom^2'), symmetry=1, barrier=(12.949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646328,0.0647935,-4.99959e-05,1.46049e-08,7.6074e-15,32874.2,32.1129], Tmin=(100,'K'), Tmax=(1049.3,'K')), NASAPolynomial(coeffs=[16.2529,0.0202986,-7.83038e-06,1.43777e-09,-1.00919e-13,28773.4,-47.8654], Tmin=(1049.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC=CC([O])[C]=O(29572)',
    structure = SMILES('C=CC=CC([O])[C]=O'),
    E0 = (185.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688494,0.0637742,-4.69009e-05,9.78896e-09,2.31241e-12,22393,32.3757], Tmin=(100,'K'), Tmax=(999.945,'K')), NASAPolynomial(coeffs=[16.4416,0.0189106,-6.83161e-06,1.22883e-09,-8.62376e-14,18335,-48.1607], Tmin=(999.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = 'C=[C][C]=CC(O)C=O(29573)',
    structure = SMILES('C=[C][C]=CC(O)C=O'),
    E0 = (179.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471593,0.0755404,-8.16923e-05,4.7141e-08,-1.08355e-11,21708.7,30.75], Tmin=(100,'K'), Tmax=(1061,'K')), NASAPolynomial(coeffs=[14.0666,0.0242868,-9.23192e-06,1.61142e-09,-1.07536e-13,18823.9,-35.6431], Tmin=(1061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[C]=CC(O)C=O(29574)',
    structure = SMILES('[CH]C=C=CC(O)C=O'),
    E0 = (204.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716162,0.0680954,-5.30057e-05,2.12635e-08,-3.50113e-12,24765.2,31.8187], Tmin=(100,'K'), Tmax=(1413.17,'K')), NASAPolynomial(coeffs=[14.0236,0.0304287,-1.30247e-05,2.40243e-09,-1.64473e-13,21004,-36.9841], Tmin=(1413.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]C(C=O)C=C1[CH]C1(29575)',
    structure = SMILES('[O]C(C=O)C=C1[CH]C1'),
    E0 = (259.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32609,0.0460734,2.81292e-07,-3.32218e-08,1.56825e-11,31294.5,28.3445], Tmin=(100,'K'), Tmax=(1006.89,'K')), NASAPolynomial(coeffs=[14.3738,0.023337,-9.19478e-06,1.75267e-09,-1.27286e-13,27191.9,-42.0171], Tmin=(1006.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C[C]=CC1[CH]OO1(29576)',
    structure = SMILES('C=C[C]=CC1[CH]OO1'),
    E0 = (489.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753582,0.0669721,-5.78162e-05,2.62764e-08,-4.84703e-12,59023.3,24.7155], Tmin=(100,'K'), Tmax=(1290.23,'K')), NASAPolynomial(coeffs=[13.8173,0.0264712,-1.07302e-05,1.94672e-09,-1.32777e-13,55652.3,-41.6384], Tmin=(1290.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC([O])[CH]O1(29551)',
    structure = SMILES('C=CC1=CC([O])[CH]O1'),
    E0 = (169.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64627,0.0534209,1.23953e-05,-7.44384e-08,3.8802e-11,20541.5,22.2536], Tmin=(100,'K'), Tmax=(905.175,'K')), NASAPolynomial(coeffs=[24.0188,0.0046999,2.71432e-06,-7.14724e-10,4.78844e-14,14075,-100.524], Tmin=(905.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(CCsJOC(O)) + radical(CC(C)OJ)"""),
)

species(
    label = 'O=CC1C=[C][CH]CO1(29577)',
    structure = SMILES('O=CC1C=[C][CH]CO1'),
    E0 = (111.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28517,0.0430298,2.42152e-05,-6.61951e-08,2.96374e-11,13485.6,21.5188], Tmin=(100,'K'), Tmax=(951.318,'K')), NASAPolynomial(coeffs=[16.3075,0.0212479,-6.68988e-06,1.18865e-09,-8.70377e-14,8754.86,-60.0474], Tmin=(951.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(36dihydro2hpyran) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=CC(=O)C=O(29578)',
    structure = SMILES('C=CC=CC(=O)C=O'),
    E0 = (-130.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.905702,0.0638732,-5.17934e-05,2.1322e-08,-3.55855e-12,-15615.5,26.1186], Tmin=(100,'K'), Tmax=(1408.2,'K')), NASAPolynomial(coeffs=[14.3295,0.0257427,-1.1177e-05,2.09341e-09,-1.44858e-13,-19396.2,-43.2387], Tmin=(1408.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=CC(O)C=O(29579)',
    structure = SMILES('C=C=C=CC(O)C=O'),
    E0 = (10.8199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816569,0.0671165,-6.01595e-05,2.76379e-08,-5.12527e-12,1418.26,29.4026], Tmin=(100,'K'), Tmax=(1283.31,'K')), NASAPolynomial(coeffs=[14.3253,0.0250111,-1.09451e-05,2.07183e-09,-1.44856e-13,-2048.98,-39.1392], Tmin=(1283.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.8199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CO(12)',
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
    label = 'C=C[C]=CC[O](26480)',
    structure = SMILES('C=C[C]=CC[O]'),
    E0 = (329.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.156902,'amu*angstrom^2'), symmetry=1, barrier=(3.6075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.08362,'amu*angstrom^2'), symmetry=1, barrier=(93.8904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3714.19,'J/mol'), sigma=(6.10713,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.15 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15325,0.0453843,-3.03584e-05,-2.52903e-09,1.30214e-11,39686.9,20.6968], Tmin=(100,'K'), Tmax=(566.82,'K')), NASAPolynomial(coeffs=[5.31978,0.0310279,-1.35096e-05,2.52177e-09,-1.74319e-13,39199.6,6.08568], Tmin=(566.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC(C=O)O1(29559)',
    structure = SMILES('C=CC1=CC(C=O)O1'),
    E0 = (-71.8343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.909895,0.0443021,3.68613e-05,-9.31426e-08,4.22548e-11,-8506.76,23.798], Tmin=(100,'K'), Tmax=(945.627,'K')), NASAPolynomial(coeffs=[23.3683,0.00880788,-1.22507e-06,2.52467e-10,-2.91039e-14,-15414.7,-97.3623], Tmin=(945.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.8343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'C=C[C]=CO[CH]C=O(29462)',
    structure = SMILES('C=C[C]=COC=C[O]'),
    E0 = (153.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39521,'amu*angstrom^2'), symmetry=1, barrier=(32.0786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39529,'amu*angstrom^2'), symmetry=1, barrier=(32.0804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39492,'amu*angstrom^2'), symmetry=1, barrier=(32.0719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4075.33,'J/mol'), sigma=(6.46898,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.56 K, Pc=34.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0187316,0.0693981,-2.39949e-05,-4.1212e-08,2.77723e-11,18592.2,28.833], Tmin=(100,'K'), Tmax=(906.39,'K')), NASAPolynomial(coeffs=[25.943,0.00424209,2.32726e-06,-6.23601e-10,4.21083e-14,11869.7,-104.849], Tmin=(906.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CC([O])=CO(29580)',
    structure = SMILES('C=C[C]=CC([O])=CO'),
    E0 = (59.8079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.22014,'amu*angstrom^2'), symmetry=1, barrier=(28.0534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20906,'amu*angstrom^2'), symmetry=1, barrier=(27.7986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21831,'amu*angstrom^2'), symmetry=1, barrier=(28.0114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.65797,0.0972808,-0.0001118,5.92559e-08,-1.13904e-11,7420.33,31.136], Tmin=(100,'K'), Tmax=(1525.95,'K')), NASAPolynomial(coeffs=[26.2443,0.00165114,4.31031e-06,-1.12902e-09,8.49263e-14,1523.11,-106.689], Tmin=(1525.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.8079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
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
    label = 'C=C[C]=C[CH]C=O(27341)',
    structure = SMILES('[CH2]C=[C]C=CC=O'),
    E0 = (251.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(1.48426,'amu*angstrom^2'), symmetry=1, barrier=(34.1261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48843,'amu*angstrom^2'), symmetry=1, barrier=(34.222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4825,'amu*angstrom^2'), symmetry=1, barrier=(34.0857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32441,0.0561728,-4.32064e-05,1.72705e-08,-2.8472e-12,30341.8,23.6367], Tmin=(100,'K'), Tmax=(1405.2,'K')), NASAPolynomial(coeffs=[11.8733,0.0261443,-1.11519e-05,2.06282e-09,-1.41578e-13,27377.2,-30.8443], Tmin=(1405.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C1OC1C=O(29511)',
    structure = SMILES('[CH2]C=[C]C1OC1C=O'),
    E0 = (236.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12257,0.0471712,1.57205e-05,-6.61625e-08,3.31727e-11,28613.2,28.9703], Tmin=(100,'K'), Tmax=(901.336,'K')), NASAPolynomial(coeffs=[18.316,0.0142372,-1.6436e-06,6.27266e-11,-2.34126e-15,23752.2,-61.9641], Tmin=(901.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[O]C1C=C=CCC1[O](29581)',
    structure = SMILES('[O]C1C=C=CCC1[O]'),
    E0 = (284.282,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.760904,0.110313,-0.000143117,9.68426e-08,-2.62052e-11,34357.9,14.6815], Tmin=(100,'K'), Tmax=(900.895,'K')), NASAPolynomial(coeffs=[16.3498,0.0343404,-1.66206e-05,3.23366e-09,-2.2835e-13,31274.9,-66.082], Tmin=(900.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C[C]=C=CC([O])C=O(29582)',
    structure = SMILES('CC#C[CH]C([O])C=O'),
    E0 = (274.574,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,2100,2250,500,550,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,434.779,434.78,434.782],'cm^-1')),
        HinderedRotor(inertia=(0.0723597,'amu*angstrom^2'), symmetry=1, barrier=(9.70653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000891769,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(6.25165e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529399,'amu*angstrom^2'), symmetry=1, barrier=(71.0155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22461,0.0539051,-2.84467e-05,-3.9214e-09,6.61879e-12,33130.1,33.6697], Tmin=(100,'K'), Tmax=(937.919,'K')), NASAPolynomial(coeffs=[12.499,0.0233742,-7.6892e-06,1.27648e-09,-8.49007e-14,30243.2,-24.1158], Tmin=(937.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'CC=C=[C]C([O])C=O(29583)',
    structure = SMILES('C[CH]C#CC([O])C=O'),
    E0 = (220.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74872,0.0632882,-5.73571e-05,2.88562e-08,-5.77663e-12,26690.1,32.7175], Tmin=(100,'K'), Tmax=(1302.95,'K')), NASAPolynomial(coeffs=[13.029,0.0223567,-6.51489e-06,9.38722e-10,-5.47922e-14,23764.3,-28.7247], Tmin=(1302.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=OCOJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'CC=C=C[C]([O])C=O(29584)',
    structure = SMILES('CC=C=CC([O])=C[O]'),
    E0 = (55.0562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.548574,0.0831534,-8.68204e-05,4.38557e-08,-8.3492e-12,6799.93,29.0053], Tmin=(100,'K'), Tmax=(1440.14,'K')), NASAPolynomial(coeffs=[22.3269,0.0105781,-1.81426e-06,1.46819e-10,-5.06346e-15,1148.43,-86.4453], Tmin=(1440.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.0562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC=C=CC([O])[C]=O(29585)',
    structure = SMILES('CC=C=CC([O])[C]=O'),
    E0 = (237.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00136081,'amu*angstrom^2'), symmetry=1, barrier=(15.4507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672518,'amu*angstrom^2'), symmetry=1, barrier=(15.4625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.672274,'amu*angstrom^2'), symmetry=1, barrier=(15.4569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07552,0.0659331,-6.43342e-05,3.43799e-08,-7.60353e-12,28717.9,31.0235], Tmin=(100,'K'), Tmax=(1074.13,'K')), NASAPolynomial(coeffs=[10.9691,0.0290905,-1.2885e-05,2.44799e-09,-1.71581e-13,26592.5,-17.4147], Tmin=(1074.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC1=CC([CH][O])O1(29539)',
    structure = SMILES('C=CC1=CC([CH][O])O1'),
    E0 = (255.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.373795,0.0661071,-3.53881e-05,-1.1082e-08,1.13639e-11,30854.3,25.9253], Tmin=(100,'K'), Tmax=(976.654,'K')), NASAPolynomial(coeffs=[20.3393,0.015435,-5.32646e-06,1.00128e-09,-7.48655e-14,25471.3,-77.5182], Tmin=(976.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C(C=O)C1[C]=CC1(29586)',
    structure = SMILES('[O]C(C=O)C1[C]=CC1'),
    E0 = (332.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44079,0.0416142,1.4638e-05,-4.99909e-08,2.21657e-11,40092.5,30.1365], Tmin=(100,'K'), Tmax=(984.224,'K')), NASAPolynomial(coeffs=[15.0446,0.0213888,-7.97371e-06,1.52025e-09,-1.12299e-13,35716.4,-43.9048], Tmin=(984.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C1[CH]OCC=C=C1(29587)',
    structure = SMILES('[O]C1[CH]OCC=C=C1'),
    E0 = (334.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.984707,0.056163,-2.7075e-05,-3.59307e-09,4.71129e-12,40371.7,21.9046], Tmin=(100,'K'), Tmax=(1092.06,'K')), NASAPolynomial(coeffs=[14.35,0.0257733,-1.08329e-05,2.05836e-09,-1.46087e-13,36345.6,-48.8206], Tmin=(1092.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'CC=C=CC(=O)C=O(29588)',
    structure = SMILES('CC=C=CC(=O)C=O'),
    E0 = (-78.0068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25934,0.0662704,-6.9558e-05,4.61257e-08,-1.36729e-11,-9288.67,24.8988], Tmin=(100,'K'), Tmax=(788.743,'K')), NASAPolynomial(coeffs=[6.44317,0.0399797,-1.95562e-05,3.86024e-09,-2.7555e-13,-10106.4,1.12034], Tmin=(788.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.0068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)H) + group(Cdd-CdsCds)"""),
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
    label = '[CH]=C=CC([O])C=O(29589)',
    structure = SMILES('[CH]=C=CC([O])C=O'),
    E0 = (268.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,322.926],'cm^-1')),
        HinderedRotor(inertia=(0.102337,'amu*angstrom^2'), symmetry=1, barrier=(7.79505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194462,'amu*angstrom^2'), symmetry=1, barrier=(14.77,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48638,0.0522564,-5.04907e-05,2.58391e-08,-5.29696e-12,32380.3,27.5444], Tmin=(100,'K'), Tmax=(1180.07,'K')), NASAPolynomial(coeffs=[11.5981,0.0179817,-6.92391e-06,1.22665e-09,-8.27843e-14,29993.8,-22.9129], Tmin=(1180.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=OCOJ)"""),
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
    E0 = (224.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (326.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (309.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (323.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (301.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (419.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (483.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (265.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (295.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (265.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (349.623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (459.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (438.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (299.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (299.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (417.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (494.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (393.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (253.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (352.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (432.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (361.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (489.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (279.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (285.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (313.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (249.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (515.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (232.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (467.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (224.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (494.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (270.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (412.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (489.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (462.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (299.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (314.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (350.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (342.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (334.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (249.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (650.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['OCHCHO(48)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=C[C]=CC1OC1[O](29560)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC1=CC([O])C1[O](29561)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(85.6196,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 83.4 to 85.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[CH2]C1[C]=CC(C=O)O1(29535)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.82166e+09,'s^-1'), n=0.527281, Ea=(99.1325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_O] + [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=C[C]=CC(=O)C=O(29562)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;HJ] for rate rule [CO-DeDe_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=CC#CC([O])C=O(29563)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C=C=C=CC([O])C=O(29564)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OCHCHO(48)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(38.9443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['HCO(16)', 'C=C[C]=CC=O(26908)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CdH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=C[O](2537)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0234394,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=C[C]=C[C](O)C=O(29565)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=CC=[C]C([O])C=O(29566)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=[C]C=CC([O])C=O(29567)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=C[C]=[C]C(O)C=O(29568)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=C[C]=CC(O)[C]=O(29569)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC=C[C]([O])C=O(29570)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC=CC([O])C=O(29571)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC=CC([O])[C]=O(29572)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.66043e+07,'s^-1'), n=1.66125, Ea=(169.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSR;Cd_rad_out_singleDe_Cd;XH_out] + [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;CO_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=[C][C]=CC(O)C=O(29573)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.658e+09,'s^-1'), n=0.699, Ea=(29.5516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_doubleC] for rate rule [R5HJ_3;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[CH]=C[C]=CC(O)C=O(29574)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.427,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_3;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C=C[O](2537)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[O]C(C=O)C=C1[CH]C1(29575)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=C[C]=CC1[CH]OO1(29576)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC1=CC([O])[CH]O1(29551)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.72e+09,'s^-1'), n=0.186, Ea=(55.7727,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleDe] for rate rule [R5_DS_CO;carbonyl_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['O=CC1C=[C][CH]CO1(29577)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC=CC(=O)C=O(29578)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=C=C=CC(O)C=O(29579)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CO(12)', 'C=C[C]=CC[O](26480)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC1=CC(C=O)O1(29559)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[C]=CO[CH]C=O(29462)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[C]=CC([O])=CO(29580)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(164.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 161.5 to 164.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', 'C=C[C]=C[CH]C=O(27341)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[CH2]C=[C]C1OC1C=O(29511)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[O]C1C=C=CCC1[O](29581)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7_SMMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_SMMS;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[C]=C=CC([O])C=O(29582)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['CC=C=[C]C([O])C=O(29583)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['CC=C=C[C]([O])C=O(29584)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.21176e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC=C=CC([O])[C]=O(29585)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.27582e+06,'s^-1'), n=1.20683, Ea=(76.1033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;Y_rad_out;Cs_H_out_2H] for rate rule [R6H;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['C=CC1=CC([CH][O])O1(29539)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[O]C(C=O)C1[C]=CC1(29586)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['[O]C1[CH]OCC=C=C1(29587)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.97949e+09,'s^-1'), n=0.649948, Ea=(110.532,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cs2H] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 109.1 to 110.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C[C]=CC([O])C=O(29458)'],
    products = ['CC=C=CC(=O)C=O(29588)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CH2(19)', '[CH]=C=CC([O])C=O(29589)'],
    products = ['C=C[C]=CC([O])C=O(29458)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4676',
    isomers = [
        'C=C[C]=CC([O])C=O(29458)',
    ],
    reactants = [
        ('OCHCHO(48)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4676',
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

