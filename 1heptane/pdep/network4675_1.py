species(
    label = 'C#CC([CH2])C([O])C=O(29457)',
    structure = SMILES('C#CC([CH2])C([O])C=O'),
    E0 = (289.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,482.678,487.449],'cm^-1')),
        HinderedRotor(inertia=(0.54176,'amu*angstrom^2'), symmetry=1, barrier=(12.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00035764,'amu*angstrom^2'), symmetry=1, barrier=(1.71898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481135,'amu*angstrom^2'), symmetry=1, barrier=(76.4926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47078,'amu*angstrom^2'), symmetry=1, barrier=(76.5333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.527693,0.0684877,-6.66024e-05,3.45847e-08,-7.08117e-12,34894,33.9991], Tmin=(100,'K'), Tmax=(1234.58,'K')), NASAPolynomial(coeffs=[14.8807,0.0203933,-6.235e-06,9.42737e-10,-5.73343e-14,31471.3,-37.7793], Tmin=(1234.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ) + radical(Isobutyl)"""),
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
    label = '[CH]=C1CC1C([O])C=O(29653)',
    structure = SMILES('[CH]=C1CC1C([O])C=O'),
    E0 = (369.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04789,0.0514753,-8.03429e-06,-3.07931e-08,1.64968e-11,44566.2,30.3517], Tmin=(100,'K'), Tmax=(978.611,'K')), NASAPolynomial(coeffs=[16.6343,0.0192431,-6.87534e-06,1.28438e-09,-9.42373e-14,40008.3,-52.2075], Tmin=(978.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CC([CH2])C1OC1[O](29654)',
    structure = SMILES('C#CC([CH2])C1OC1[O]'),
    E0 = (340.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438912,0.071842,-7.67812e-05,4.57636e-08,-1.04747e-11,41092,29.4059], Tmin=(100,'K'), Tmax=(1246.97,'K')), NASAPolynomial(coeffs=[12.4202,0.0227345,-4.86889e-06,4.52477e-10,-1.41921e-14,38933.9,-27.7137], Tmin=(1246.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.546,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CtCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C1OC(C=O)C1[CH2](29554)',
    structure = SMILES('[CH]=C1OC(C=O)C1[CH2]'),
    E0 = (248.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01831,0.0418268,4.65639e-05,-1.10993e-07,5.24403e-11,30024,26.1919], Tmin=(100,'K'), Tmax=(899.881,'K')), NASAPolynomial(coeffs=[23.9218,0.0036337,4.19062e-06,-1.04503e-09,7.11218e-14,23326.2,-96.1994], Tmin=(899.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cds_P) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC1CC([O])C1[O](29632)',
    structure = SMILES('C#CC1CC([O])C1[O]'),
    E0 = (366.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05681,0.0486461,9.53648e-06,-5.45564e-08,2.65947e-11,44243,25.2073], Tmin=(100,'K'), Tmax=(943.287,'K')), NASAPolynomial(coeffs=[17.8754,0.0175306,-4.91546e-06,8.40781e-10,-6.22088e-14,39281.4,-64.4317], Tmin=(943.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = 'C#CC(=C)C([O])C=O(29655)',
    structure = SMILES('C#CC(=C)C([O])C=O'),
    E0 = (197.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2950,3100,1380,975,1025,1650,750,770,3400,2100,350,440,435,1725,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0123294,'amu*angstrom^2'), symmetry=1, barrier=(3.94861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0535792,'amu*angstrom^2'), symmetry=1, barrier=(17.157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106479,'amu*angstrom^2'), symmetry=1, barrier=(17.157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01235,0.0619287,-5.65565e-05,2.67389e-08,-5.07866e-12,23867.2,29.0578], Tmin=(100,'K'), Tmax=(1262.44,'K')), NASAPolynomial(coeffs=[13.6186,0.0219861,-9.0978e-06,1.67701e-09,-1.15677e-13,20684.3,-34.6981], Tmin=(1262.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCtCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CC([CH2])C(=O)C=O(29656)',
    structure = SMILES('C#CC([CH2])C(=O)C=O'),
    E0 = (138.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,375,552.5,462.5,1710,750,770,3400,2100,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.311908,'amu*angstrom^2'), symmetry=1, barrier=(7.17139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658229,'amu*angstrom^2'), symmetry=1, barrier=(15.134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.57138,'amu*angstrom^2'), symmetry=1, barrier=(59.1212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.57223,'amu*angstrom^2'), symmetry=1, barrier=(59.1407,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.310475,0.0882209,-0.000136742,1.15118e-07,-3.80048e-11,16780.5,27.1933], Tmin=(100,'K'), Tmax=(854.921,'K')), NASAPolynomial(coeffs=[9.98913,0.0306199,-1.40683e-05,2.60542e-09,-1.75685e-13,15575.7,-15.3511], Tmin=(854.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)C=O)"""),
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
    label = 'C=CC([O])C=O(2572)',
    structure = SMILES('C=CC([O])C=O'),
    E0 = (-26.6044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.796297,'amu*angstrom^2'), symmetry=1, barrier=(18.3084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0379626,'amu*angstrom^2'), symmetry=1, barrier=(18.3123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03604,0.0365188,-1.44093e-05,-5.26523e-09,3.93449e-12,-3123.41,24.0883], Tmin=(100,'K'), Tmax=(1084.25,'K')), NASAPolynomial(coeffs=[10.1946,0.0190919,-7.83047e-06,1.46841e-09,-1.03414e-13,-5637.41,-19.3665], Tmin=(1084.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-26.6044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
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
    label = 'C#CC([CH2])C=O(26979)',
    structure = SMILES('C#CC([CH2])C=O'),
    E0 = (246.464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.718944,'amu*angstrom^2'), symmetry=1, barrier=(16.5299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.476468,'amu*angstrom^2'), symmetry=1, barrier=(10.9549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.70346,'amu*angstrom^2'), symmetry=1, barrier=(62.1579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38307,0.0611577,-8.56871e-05,6.78678e-08,-2.13318e-11,29733.6,20.9105], Tmin=(100,'K'), Tmax=(879.553,'K')), NASAPolynomial(coeffs=[8.30862,0.023097,-9.5821e-06,1.6972e-09,-1.1177e-13,28769.2,-10.1689], Tmin=(879.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C#C[C](C)C([O])C=O(29657)',
    structure = SMILES('[CH]=C=C(C)C([O])C=O'),
    E0 = (229.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,355.682],'cm^-1')),
        HinderedRotor(inertia=(0.111529,'amu*angstrom^2'), symmetry=1, barrier=(9.85447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110821,'amu*angstrom^2'), symmetry=1, barrier=(9.85504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110774,'amu*angstrom^2'), symmetry=1, barrier=(9.84021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.853183,0.0697191,-7.19308e-05,4.08344e-08,-9.46218e-12,27702.9,30.5785], Tmin=(100,'K'), Tmax=(1038.39,'K')), NASAPolynomial(coeffs=[11.6837,0.0279995,-1.16662e-05,2.14413e-09,-1.47417e-13,25453.6,-22.0807], Tmin=(1038.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=OCOJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC([CH2])[C](O)C=O(29658)',
    structure = SMILES('C#CC([CH2])C(O)=C[O]'),
    E0 = (119.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.22285,0.0947262,-0.000109378,5.90852e-08,-1.17473e-11,14628.9,31.5928], Tmin=(100,'K'), Tmax=(1433.92,'K')), NASAPolynomial(coeffs=[25.2834,0.00447957,2.08506e-06,-6.67334e-10,5.30767e-14,8703.74,-99.9923], Tmin=(1433.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = 'C#CC(C)[C]([O])C=O(29659)',
    structure = SMILES('C#CC(C)C([O])=C[O]'),
    E0 = (52.6468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.623668,0.0860305,-9.34009e-05,4.88815e-08,-9.6265e-12,6511.82,29.1671], Tmin=(100,'K'), Tmax=(1389.15,'K')), NASAPolynomial(coeffs=[22.6852,0.00988391,-1.42807e-06,6.3659e-11,1.1438e-15,907.141,-87.8106], Tmin=(1389.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.6468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C#C[C]([CH2])C(O)C=O(29660)',
    structure = SMILES('[CH]=C=C([CH2])C(O)C=O'),
    E0 = (137.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28564,0.0740214,-7.31545e-05,3.68724e-08,-7.2971e-12,16637,32.1125], Tmin=(100,'K'), Tmax=(1235.19,'K')), NASAPolynomial(coeffs=[17.0815,0.0196298,-7.1014e-06,1.22136e-09,-8.13228e-14,12487.8,-52.4653], Tmin=(1235.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = 'C#CC([CH2])C(O)[C]=O(29661)',
    structure = SMILES('C#CC([CH2])C(O)[C]=O'),
    E0 = (205.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.167296,0.0749685,-8.01135e-05,4.39584e-08,-9.267e-12,24832.8,35.3282], Tmin=(100,'K'), Tmax=(1270.19,'K')), NASAPolynomial(coeffs=[17.1122,0.0160047,-3.86606e-06,4.6728e-10,-2.36282e-14,20980,-48.6952], Tmin=(1270.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C#CC(C)C([O])[C]=O(29662)',
    structure = SMILES('C#CC(C)C([O])[C]=O'),
    E0 = (243.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687783,0.0668112,-5.88474e-05,2.31123e-08,-2.12015e-12,29460.1,32.4622], Tmin=(100,'K'), Tmax=(959.712,'K')), NASAPolynomial(coeffs=[14.8942,0.0205967,-6.92878e-06,1.15746e-09,-7.67767e-14,26134.8,-38.6095], Tmin=(959.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O) + radical(C=OCOJ)"""),
)

species(
    label = '[C]#CC(C)C([O])C=O(29663)',
    structure = SMILES('[C]#CC(C)C([O])C=O'),
    E0 = (421.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2782.5,750,1395,475,1775,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,258.664,258.743,2025.91],'cm^-1')),
        HinderedRotor(inertia=(0.00251592,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0025212,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.306224,'amu*angstrom^2'), symmetry=1, barrier=(14.5598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50519,'amu*angstrom^2'), symmetry=1, barrier=(71.5399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771468,0.0675355,-6.69023e-05,3.6413e-08,-7.97108e-12,50764.8,32.085], Tmin=(100,'K'), Tmax=(1109.14,'K')), NASAPolynomial(coeffs=[12.6887,0.024557,-8.77756e-06,1.47576e-09,-9.61568e-14,48121.3,-26.6429], Tmin=(1109.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(C=OCOJ)"""),
)

species(
    label = '[C]#CC([CH2])C(O)C=O(29664)',
    structure = SMILES('[C]#CC([CH2])C(O)C=O'),
    E0 = (382.441,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,2782.5,750,1395,475,1775,1000,3615,1277.5,1000,3000,3100,440,815,1455,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,290.672,1918.46],'cm^-1')),
        HinderedRotor(inertia=(0.176831,'amu*angstrom^2'), symmetry=1, barrier=(11.5202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284538,'amu*angstrom^2'), symmetry=1, barrier=(0.160954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00191367,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14898,'amu*angstrom^2'), symmetry=1, barrier=(69.7849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06918,'amu*angstrom^2'), symmetry=1, barrier=(69.5577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.681137,0.0704691,-6.91567e-05,3.14476e-08,-3.54154e-12,46119,33.4177], Tmin=(100,'K'), Tmax=(848.281,'K')), NASAPolynomial(coeffs=[13.8396,0.0218107,-6.7913e-06,1.04182e-09,-6.43796e-14,43404.9,-30.7389], Tmin=(848.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.441,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(Isobutyl)"""),
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
    label = 'C#CC([CH2])C1[CH]OO1(29665)',
    structure = SMILES('C#CC([CH2])C1[CH]OO1'),
    E0 = (553.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.786287,0.0714472,-7.75737e-05,4.88277e-08,-1.26041e-11,66723.6,25.2931], Tmin=(100,'K'), Tmax=(937.851,'K')), NASAPolynomial(coeffs=[10.431,0.030311,-1.17791e-05,2.05697e-09,-1.36326e-13,64914.5,-20.6181], Tmin=(937.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxetane) + radical(CCsJOO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[C]=COC1C=O(29504)',
    structure = SMILES('[CH2]C1[C]=COC1C=O'),
    E0 = (166.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44109,0.0313614,7.01584e-05,-1.31629e-07,5.89695e-11,20082.5,24.782], Tmin=(100,'K'), Tmax=(900.615,'K')), NASAPolynomial(coeffs=[22.2636,0.00502872,3.84408e-06,-9.87759e-10,6.68174e-14,13649.2,-88.3886], Tmin=(900.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(2,3-Dihydrofuran) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1CO[CH]C1[O](29623)',
    structure = SMILES('C#CC1CO[CH]C1[O]'),
    E0 = (276.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06062,0.0504285,9.54317e-06,-6.6139e-08,3.64543e-11,33433.1,21.5217], Tmin=(100,'K'), Tmax=(855.502,'K')), NASAPolynomial(coeffs=[18.5896,0.0117459,1.48895e-06,-7.3278e-10,6.15529e-14,28850.3,-69.5653], Tmin=(855.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CCsJOCs)"""),
)

species(
    label = 'C#CC(=C)C(O)C=O(29666)',
    structure = SMILES('C#CC(=C)C(O)C=O'),
    E0 = (-46.2098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.591777,0.0668892,-5.94869e-05,2.67827e-08,-4.77591e-12,-5428.33,29.792], Tmin=(100,'K'), Tmax=(1354.6,'K')), NASAPolynomial(coeffs=[16.3449,0.0203719,-7.97681e-06,1.43216e-09,-9.73363e-14,-9696.21,-50.9892], Tmin=(1354.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.2098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCtCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(C)C(=O)C=O(29667)',
    structure = SMILES('C#CC(C)C(=O)C=O'),
    E0 = (-72.0385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708334,0.0767667,-9.53964e-05,6.93146e-08,-2.08258e-11,-8549.71,25.4225], Tmin=(100,'K'), Tmax=(805.794,'K')), NASAPolynomial(coeffs=[9.40686,0.033587,-1.50173e-05,2.81401e-09,-1.93968e-13,-9951.56,-14.6646], Tmin=(805.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-72.0385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
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
    label = 'C#CC([CH2])C[O](26479)',
    structure = SMILES('C#CC([CH2])C[O]'),
    E0 = (387.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.296829,'amu*angstrom^2'), symmetry=1, barrier=(6.82467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259015,'amu*angstrom^2'), symmetry=1, barrier=(65.2663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0259307,'amu*angstrom^2'), symmetry=1, barrier=(65.2878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3666,'J/mol'), sigma=(6.23254,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=572.62 K, Pc=34.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65594,0.054597,-6.64663e-05,5.12778e-08,-1.60921e-11,46703.4,23.6456], Tmin=(100,'K'), Tmax=(908.831,'K')), NASAPolynomial(coeffs=[6.05667,0.0280814,-1.09074e-05,1.87035e-09,-1.20976e-13,46198.6,4.45916], Tmin=(908.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(CCOJ)"""),
)

species(
    label = 'C#C[CH]CC([O])C=O(29456)',
    structure = SMILES('C#C[CH]CC([O])C=O'),
    E0 = (239.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,414.812,415.46,417.033],'cm^-1')),
        HinderedRotor(inertia=(0.000964841,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0025394,'amu*angstrom^2'), symmetry=1, barrier=(13.55,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624137,'amu*angstrom^2'), symmetry=1, barrier=(76.3969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628933,'amu*angstrom^2'), symmetry=1, barrier=(76.3932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4157.55,'J/mol'), sigma=(6.76601,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=649.40 K, Pc=30.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723544,0.0675002,-6.3529e-05,3.07116e-08,-5.42362e-12,28901.1,31.6329], Tmin=(100,'K'), Tmax=(975.151,'K')), NASAPolynomial(coeffs=[13.4857,0.0232898,-8.04314e-06,1.33747e-09,-8.72547e-14,26025.2,-31.6001], Tmin=(975.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CC[CH]C([O])C=O(29635)',
    structure = SMILES('C#CC[CH]C([O])C=O'),
    E0 = (292.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2175,525,750,770,3400,2100,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,180,855.535,3311.6],'cm^-1')),
        HinderedRotor(inertia=(0.0350033,'amu*angstrom^2'), symmetry=1, barrier=(18.2141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0980158,'amu*angstrom^2'), symmetry=1, barrier=(2.25358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.792727,'amu*angstrom^2'), symmetry=1, barrier=(18.2263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157492,'amu*angstrom^2'), symmetry=1, barrier=(81.9043,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851198,0.0622344,-4.90523e-05,1.66783e-08,-1.02986e-12,35356.3,33.8334], Tmin=(100,'K'), Tmax=(1021.57,'K')), NASAPolynomial(coeffs=[14.0662,0.0224313,-8.14135e-06,1.42213e-09,-9.64419e-14,32033.2,-33.2529], Tmin=(1021.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'C#CC1COC1C=O(29465)',
    structure = SMILES('C#CC1COC1C=O'),
    E0 = (24.2467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54952,0.0370236,4.07892e-05,-9.16829e-08,4.31045e-11,3020.49,24.493], Tmin=(100,'K'), Tmax=(879.973,'K')), NASAPolynomial(coeffs=[16.7803,0.0144513,-2.72212e-07,-3.1676e-10,2.81283e-14,-1466.65,-57.3045], Tmin=(879.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.2467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CtCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane)"""),
)

species(
    label = 'C#CC([CH2])O[CH]C=O(29461)',
    structure = SMILES('C#CC([CH2])OC=C[O]'),
    E0 = (190.866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29346,'amu*angstrom^2'), symmetry=1, barrier=(29.7393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29354,'amu*angstrom^2'), symmetry=1, barrier=(29.7411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29345,'amu*angstrom^2'), symmetry=1, barrier=(29.7389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29334,'amu*angstrom^2'), symmetry=1, barrier=(29.7364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4029.45,'J/mol'), sigma=(6.59937,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=629.39 K, Pc=31.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.29245,0.100556,-0.000111661,5.59707e-08,-1.0091e-11,23215.1,35.4589], Tmin=(100,'K'), Tmax=(1650.12,'K')), NASAPolynomial(coeffs=[28.382,-0.000953895,5.29641e-06,-1.25364e-09,8.95498e-14,16788.5,-116.689], Tmin=(1650.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJC(C)OC) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([CH2])C([O])=CO(29668)',
    structure = SMILES('C#CC([CH2])C([O])=CO'),
    E0 = (116.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.983551,'amu*angstrom^2'), symmetry=1, barrier=(22.6138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983466,'amu*angstrom^2'), symmetry=1, barrier=(22.6118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984956,'amu*angstrom^2'), symmetry=1, barrier=(22.6461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.984672,'amu*angstrom^2'), symmetry=1, barrier=(22.6395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.16017,0.0965669,-0.000116813,6.62553e-08,-1.38178e-11,14184.1,31.4512], Tmin=(100,'K'), Tmax=(1367.37,'K')), NASAPolynomial(coeffs=[24.6591,0.00443065,2.47709e-06,-7.86932e-10,6.35357e-14,8675.61,-95.5138], Tmin=(1367.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(Isobutyl)"""),
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
    label = 'C#CC([CH2])[CH]C=O(27326)',
    structure = SMILES('C#CC([CH2])C=C[O]'),
    E0 = (334.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28068,'amu*angstrom^2'), symmetry=1, barrier=(29.4453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28073,'amu*angstrom^2'), symmetry=1, barrier=(29.4465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28038,'amu*angstrom^2'), symmetry=1, barrier=(29.4384,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0458433,0.0727493,-7.50059e-05,3.82979e-08,-7.3957e-12,40383.5,26.988], Tmin=(100,'K'), Tmax=(1429.16,'K')), NASAPolynomial(coeffs=[18.8441,0.0117163,-2.11059e-06,1.71804e-10,-5.28237e-15,35870.2,-67.4071], Tmin=(1429.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(C=COJ)"""),
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
    E0 = (289.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (369.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (391.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (411.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (366.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (424.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (382.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (289.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (550.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (289.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (372.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (432.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (442.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (406.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (364.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (364.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (365.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (517.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (426.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (432.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (419.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (554.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (347.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (326.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (352.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (352.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (573.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (448.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (538.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (297.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (504.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (289.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (650.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (577.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['OCHCHO(48)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['[CH]=C1CC1C([O])C=O(29653)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.881e+08,'s^-1'), n=1.062, Ea=(80.5299,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 18 used for R4_S_T;triplebond_intra_H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 78.7 to 80.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC([CH2])C1OC1[O](29654)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['[CH]=C1OC(C=O)C1[CH2](29554)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC1CC([O])C1[O](29632)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.46159e+06,'s^-1'), n=1.55572, Ea=(77.8233,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 76.0 to 77.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#CC(=C)C([O])C=O(29655)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(15.6063,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2632 used for Cds-CtCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CtCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C#CC([CH2])C(=O)C=O(29656)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C[O](2537)', 'CH2CHCCH(26391)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00589682,'m^3/(mol*s)'), n=2.48333, Ea=(33.6885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CtH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 31.2 to 33.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H(33)', 'C=CC([O])C=O(2572)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CJ] for rate rule [Cds-CsH_Cds-HH;CtJ_Ct]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OCHCHO(48)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(210.531,'m^3/(mol*s)'), n=1.35673, Ea=(62.7465,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-DeH_O;YJ] for rate rule [CO-DeH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 57.9 to 62.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['HCO(16)', 'C#CC([CH2])C=O(26979)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO_O;CO_pri_rad] for rate rule [CO-CsH_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#C[C](C)C([O])C=O(29657)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC([CH2])[C](O)C=O(29658)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC(C)[C]([O])C=O(29659)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#C[C]([CH2])C(O)C=O(29660)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC([CH2])C(O)[C]=O(29661)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC(C)C([O])[C]=O(29662)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(463.959,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#CC(C)C([O])C=O(29663)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.39293e+07,'s^-1'), n=1.32074, Ea=(96.0416,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_TSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#CC([CH2])C(O)C=O(29664)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C=C[O](2537)', '[CH]=[C]C=C(4699)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['[O]C(C=O)C1[C]=CC1(29586)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(130.68,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC([CH2])C1[CH]OO1(29665)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(265.786,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['[CH2]C1[C]=COC1C=O(29504)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC1CO[CH]C1[O](29623)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.64784e+07,'s^-1'), n=0.990488, Ea=(37.8683,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC(=C)C(O)C=O(29666)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC(C)C(=O)C=O(29667)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CO(12)', 'C#CC([CH2])C[O](26479)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.118397,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#C[CH]CC([O])C=O(29456)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC[CH]C([O])C=O(29635)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([CH2])C([O])C=O(29457)'],
    products = ['C#CC1COC1C=O(29465)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([CH2])O[CH]C=O(29461)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([CH2])C([O])=CO(29668)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(605.045,'s^-1'), n=2.96, Ea=(172.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_H]
Euclidian distance = 0
family: ketoenol
Ea raised from 169.2 to 172.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH2(19)', '[CH]=C=CC([O])C=O(29589)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(4)', 'C#CC([CH2])[CH]C=O(27326)'],
    products = ['C#CC([CH2])C([O])C=O(29457)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '4675',
    isomers = [
        'C#CC([CH2])C([O])C=O(29457)',
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
    label = '4675',
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

