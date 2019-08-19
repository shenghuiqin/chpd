species(
    label = 'C=[C]C(C)C=[C]C=C(26582)',
    structure = SMILES('C=[C]C(C)C=[C]C=C'),
    E0 = (543.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.7818,'amu*angstrom^2'), symmetry=1, barrier=(17.9751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781839,'amu*angstrom^2'), symmetry=1, barrier=(17.976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781771,'amu*angstrom^2'), symmetry=1, barrier=(17.9745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781781,'amu*angstrom^2'), symmetry=1, barrier=(17.9747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194054,0.0844721,-8.08726e-05,4.40912e-08,-1.00072e-11,65539.2,28.7822], Tmin=(100,'K'), Tmax=(1049.12,'K')), NASAPolynomial(coeffs=[12.0739,0.039178,-1.61131e-05,2.93984e-09,-2.0115e-13,63046.5,-29.1009], Tmin=(1049.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
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
    label = '[CH2]C1[C]=CC(C)C1=C(28074)',
    structure = SMILES('[CH2]C1[C]=CC(C)C1=C'),
    E0 = (510.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890564,0.0655118,-3.97256e-05,1.25488e-08,-1.67154e-12,61551.2,24.5099], Tmin=(100,'K'), Tmax=(1645.37,'K')), NASAPolynomial(coeffs=[11.8064,0.0389747,-1.55333e-05,2.74659e-09,-1.82185e-13,57959,-33.5884], Tmin=(1645.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
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
    label = 'C=C=C(C)C=[C]C=C(28075)',
    structure = SMILES('C=C=C(C)C=[C]C=C'),
    E0 = (447.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2031,'amu*angstrom^2'), symmetry=1, barrier=(27.6616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20132,'amu*angstrom^2'), symmetry=1, barrier=(27.6207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2029,'amu*angstrom^2'), symmetry=1, barrier=(27.6569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.446284,0.0867037,-8.17487e-05,4.01044e-08,-7.75109e-12,54012.5,27.6494], Tmin=(100,'K'), Tmax=(1264.38,'K')), NASAPolynomial(coeffs=[19.1421,0.0247325,-8.22758e-06,1.3383e-09,-8.59097e-14,49059.2,-71.448], Tmin=(1264.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C(C)C#CC=C(28076)',
    structure = SMILES('C=[C]C(C)C#CC=C'),
    E0 = (491.025,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.911558,'amu*angstrom^2'), symmetry=1, barrier=(20.9585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910296,'amu*angstrom^2'), symmetry=1, barrier=(20.9295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910574,'amu*angstrom^2'), symmetry=1, barrier=(20.9359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.185817,0.075093,-5.89663e-05,2.38428e-08,-3.89122e-12,59201.2,29.4442], Tmin=(100,'K'), Tmax=(1450.73,'K')), NASAPolynomial(coeffs=[16.8034,0.0292745,-1.15917e-05,2.07227e-09,-1.39577e-13,54379.7,-56.9087], Tmin=(1450.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)C=C=C=C(28077)',
    structure = SMILES('C=[C]C(C)C=C=C=C'),
    E0 = (574.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,1685,370,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.809551,'amu*angstrom^2'), symmetry=1, barrier=(18.6132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810178,'amu*angstrom^2'), symmetry=1, barrier=(18.6276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.810784,'amu*angstrom^2'), symmetry=1, barrier=(18.6415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.509789,0.0818615,-8.53057e-05,5.29075e-08,-1.40248e-11,69179.1,27.1486], Tmin=(100,'K'), Tmax=(893.763,'K')), NASAPolynomial(coeffs=[9.2939,0.0425488,-1.93278e-05,3.69436e-09,-2.5914e-13,67608.9,-14.2431], Tmin=(893.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C=[C]C=C(28078)',
    structure = SMILES('C#CC(C)C=[C]C=C'),
    E0 = (447.498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1564,'amu*angstrom^2'), symmetry=1, barrier=(26.588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1546,'amu*angstrom^2'), symmetry=1, barrier=(26.5464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1555,'amu*angstrom^2'), symmetry=1, barrier=(26.5672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15537,'amu*angstrom^2'), symmetry=1, barrier=(26.5642,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.157,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.172931,0.0868791,-8.59886e-05,4.55458e-08,-9.649e-12,53976,26.9524], Tmin=(100,'K'), Tmax=(1145.94,'K')), NASAPolynomial(coeffs=[16.2479,0.02956,-1.09587e-05,1.89549e-09,-1.2604e-13,50212.6,-54.5053], Tmin=(1145.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJC=C)"""),
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
    label = 'C=C=CC=[C]C=C(28079)',
    structure = SMILES('C=C=CC=[C]C=C'),
    E0 = (485.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.39542,'amu*angstrom^2'), symmetry=1, barrier=(32.0835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39881,'amu*angstrom^2'), symmetry=1, barrier=(32.1614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687433,0.0633415,-3.97624e-05,-1.28671e-09,7.52642e-12,58492.1,22.9764], Tmin=(100,'K'), Tmax=(941.087,'K')), NASAPolynomial(coeffs=[16.6902,0.0186545,-5.72336e-06,9.43693e-10,-6.42868e-14,54447,-58.745], Tmin=(941.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C)"""),
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
    label = 'C=C[C]=C[C](C)C=C(27245)',
    structure = SMILES('[CH2]C=C(C)C=[C]C=C'),
    E0 = (389.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18545,'amu*angstrom^2'), symmetry=1, barrier=(27.2558,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18514,'amu*angstrom^2'), symmetry=1, barrier=(27.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66756,'amu*angstrom^2'), symmetry=1, barrier=(38.3405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18422,'amu*angstrom^2'), symmetry=1, barrier=(27.2275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0776901,0.074885,-3.83502e-05,-8.11734e-09,1.02395e-11,46952.8,28.2695], Tmin=(100,'K'), Tmax=(952.817,'K')), NASAPolynomial(coeffs=[17.1906,0.0299609,-1.00027e-05,1.69779e-09,-1.15227e-13,42469.9,-59.8753], Tmin=(952.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C(C)[C]=CC=C(28080)',
    structure = SMILES('C=[C]C(C)[C]=CC=C'),
    E0 = (582.636,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.62905,'amu*angstrom^2'), symmetry=1, barrier=(14.4631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629287,'amu*angstrom^2'), symmetry=1, barrier=(14.4685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629529,'amu*angstrom^2'), symmetry=1, barrier=(14.4741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629059,'amu*angstrom^2'), symmetry=1, barrier=(14.4633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.305564,0.0817683,-7.39902e-05,3.71916e-08,-7.79582e-12,70207.5,29.2571], Tmin=(100,'K'), Tmax=(1125.85,'K')), NASAPolynomial(coeffs=[12.4935,0.0384656,-1.62961e-05,3.02786e-09,-2.09497e-13,67463.1,-30.9871], Tmin=(1125.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.636,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CC(C)[C]=C(28081)',
    structure = SMILES('C=[C]C=CC(C)[C]=C'),
    E0 = (543.789,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.7818,'amu*angstrom^2'), symmetry=1, barrier=(17.9751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781839,'amu*angstrom^2'), symmetry=1, barrier=(17.976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781771,'amu*angstrom^2'), symmetry=1, barrier=(17.9745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781781,'amu*angstrom^2'), symmetry=1, barrier=(17.9747,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194054,0.0844721,-8.08726e-05,4.40912e-08,-1.00072e-11,65539.2,28.7822], Tmin=(100,'K'), Tmax=(1049.12,'K')), NASAPolynomial(coeffs=[12.0739,0.039178,-1.61131e-05,2.93984e-09,-2.0115e-13,63046.5,-29.1009], Tmin=(1049.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CC(C)C=[C]C=C(27252)',
    structure = SMILES('[CH]=CC(C)C=[C]C=C'),
    E0 = (553.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.857238,'amu*angstrom^2'), symmetry=1, barrier=(19.7096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856641,'amu*angstrom^2'), symmetry=1, barrier=(19.6959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.856446,'amu*angstrom^2'), symmetry=1, barrier=(19.6914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.857173,'amu*angstrom^2'), symmetry=1, barrier=(19.7081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0236462,0.0849807,-7.85072e-05,3.99054e-08,-8.30831e-12,66661,29.3187], Tmin=(100,'K'), Tmax=(1148.22,'K')), NASAPolynomial(coeffs=[14.165,0.0357162,-1.41484e-05,2.53748e-09,-1.72109e-13,63413.6,-40.8596], Tmin=(1148.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)C=[C]C=C(26570)',
    structure = SMILES('[CH2]C(C=C)C=[C]C=C'),
    E0 = (511.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.550407,'amu*angstrom^2'), symmetry=1, barrier=(12.6549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.548943,'amu*angstrom^2'), symmetry=1, barrier=(12.6213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.14273,'amu*angstrom^2'), symmetry=1, barrier=(72.2574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05201,'amu*angstrom^2'), symmetry=1, barrier=(24.1878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3606.27,'J/mol'), sigma=(6.17911,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.29 K, Pc=34.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00957748,0.0846348,-7.96538e-05,4.20261e-08,-9.04714e-12,61610.1,30.815], Tmin=(100,'K'), Tmax=(1118.3,'K')), NASAPolynomial(coeffs=[13.9114,0.034841,-1.28637e-05,2.20926e-09,-1.45848e-13,58496.5,-37.9018], Tmin=(1118.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(511.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]C(C)C=C(27248)',
    structure = SMILES('[CH2][CH]C#CC(C)C=C'),
    E0 = (483.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.124692,0.0786326,-6.65731e-05,3.1608e-08,-6.15107e-12,58298.3,32.5288], Tmin=(100,'K'), Tmax=(1229.18,'K')), NASAPolynomial(coeffs=[13.8689,0.0339058,-1.19913e-05,2.0043e-09,-1.2998e-13,54919.6,-36.6147], Tmin=(1229.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=[C][C](C)C=CC=C(28082)',
    structure = SMILES('[CH2]C=CC=C(C)[C]=C'),
    E0 = (389.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0777078,0.0748848,-3.83494e-05,-8.11842e-09,1.02399e-11,46952.8,28.2694], Tmin=(100,'K'), Tmax=(952.813,'K')), NASAPolynomial(coeffs=[17.1906,0.029961,-1.00027e-05,1.69781e-09,-1.15228e-13,42469.9,-59.875], Tmin=(952.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC=CC(C)[C]=C(28083)',
    structure = SMILES('[CH]=CC=CC(C)[C]=C'),
    E0 = (591.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.701889,'amu*angstrom^2'), symmetry=1, barrier=(16.1378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701229,'amu*angstrom^2'), symmetry=1, barrier=(16.1226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702952,'amu*angstrom^2'), symmetry=1, barrier=(16.1622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.70101,'amu*angstrom^2'), symmetry=1, barrier=(16.1176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0964679,0.0827097,-7.30219e-05,3.46488e-08,-6.72588e-12,71331,29.9338], Tmin=(100,'K'), Tmax=(1224.1,'K')), NASAPolynomial(coeffs=[14.8006,0.034661,-1.41439e-05,2.58292e-09,-1.77034e-13,67731.1,-43.9786], Tmin=(1224.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(591.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=C)C=CC=C(27251)',
    structure = SMILES('[CH2]C([C]=C)C=CC=C'),
    E0 = (549.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,230.952,230.965,230.966],'cm^-1')),
        HinderedRotor(inertia=(3.15952,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318958,'amu*angstrom^2'), symmetry=1, barrier=(12.0723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318993,'amu*angstrom^2'), symmetry=1, barrier=(12.0719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318871,'amu*angstrom^2'), symmetry=1, barrier=(12.0722,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0664527,0.0823168,-7.39762e-05,3.65044e-08,-7.35287e-12,66280,31.4195], Tmin=(100,'K'), Tmax=(1188.3,'K')), NASAPolynomial(coeffs=[14.4223,0.0339929,-1.29766e-05,2.28215e-09,-1.5303e-13,62868.2,-40.3158], Tmin=(1188.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)[CH]C=C[CH2](27701)',
    structure = SMILES('C#CC(C)[CH]C=C[CH2]'),
    E0 = (453.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376063,0.067864,-2.58547e-05,-1.51028e-08,1.09517e-11,54682.5,28.8353], Tmin=(100,'K'), Tmax=(988.294,'K')), NASAPolynomial(coeffs=[15.5295,0.0325768,-1.18266e-05,2.09959e-09,-1.45151e-13,50415.4,-50.528], Tmin=(988.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_P) + radical(Allyl_S)"""),
)

species(
    label = 'C=[C][C]=CC(C)C=C(27254)',
    structure = SMILES('[CH2]C#C[CH]C(C)C=C'),
    E0 = (458.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503659,0.0677288,-3.27883e-05,-5.05392e-09,7.18473e-12,55277.3,31.7469], Tmin=(100,'K'), Tmax=(972.674,'K')), NASAPolynomial(coeffs=[13.6255,0.0339965,-1.19652e-05,2.05622e-09,-1.38477e-13,51767.7,-36.114], Tmin=(972.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CC(C)C=C(27255)',
    structure = SMILES('[CH]C=C=CC(C)C=C'),
    E0 = (530.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459244,0.0783655,-5.5088e-05,2.0881e-08,-3.40392e-12,63921.1,28.8878], Tmin=(100,'K'), Tmax=(1365.73,'K')), NASAPolynomial(coeffs=[11.4314,0.0462294,-1.97921e-05,3.65157e-09,-2.49993e-13,60924.1,-27.4666], Tmin=(1365.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C(C)C=C1[CH]C1(28084)',
    structure = SMILES('C=[C]C(C)C=C1[CH]C1'),
    E0 = (578.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448659,0.0677431,-3.53563e-05,2.43806e-09,2.60855e-12,69765.6,27.3478], Tmin=(100,'K'), Tmax=(1137.55,'K')), NASAPolynomial(coeffs=[14.1838,0.0356269,-1.43438e-05,2.62811e-09,-1.81365e-13,65593.8,-45.2881], Tmin=(1137.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(578.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Allyl_S)"""),
)

species(
    label = 'C=C1C[CH][C]=CC1C(28085)',
    structure = SMILES('C=C1C[CH][C]=CC1C'),
    E0 = (408.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34643,0.039656,5.17881e-05,-8.40352e-08,3.11047e-11,49280.2,19.3496], Tmin=(100,'K'), Tmax=(1038.46,'K')), NASAPolynomial(coeffs=[12.6846,0.0420227,-1.8133e-05,3.54585e-09,-2.57793e-13,44442.9,-47.7309], Tmin=(1038.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C=C=C(C)C=CC=C(28086)',
    structure = SMILES('C=C=C(C)C=CC=C'),
    E0 = (248.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0944117,0.0769951,-4.19201e-05,-6.16184e-09,9.48088e-12,30068.5,26.7863], Tmin=(100,'K'), Tmax=(982.092,'K')), NASAPolynomial(coeffs=[19.2485,0.0270022,-9.53525e-06,1.70366e-09,-1.19696e-13,24880.8,-73.2512], Tmin=(982.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=C=C=CC(C)C=C(27261)',
    structure = SMILES('C=C=C=CC(C)C=C'),
    E0 = (336.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51276,0.0778454,-6.34769e-05,2.84516e-08,-5.40623e-12,40576.5,26.6472], Tmin=(100,'K'), Tmax=(1213.46,'K')), NASAPolynomial(coeffs=[11.7318,0.0408632,-1.77618e-05,3.33594e-09,-2.3182e-13,37853.7,-29.6488], Tmin=(1213.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = 'C=[C]CC=[C]C=C(27963)',
    structure = SMILES('C=[C]CC=[C]C=C'),
    E0 = (577.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1670,1700,300,440,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.900581,'amu*angstrom^2'), symmetry=1, barrier=(20.7061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901126,'amu*angstrom^2'), symmetry=1, barrier=(20.7187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.90036,'amu*angstrom^2'), symmetry=1, barrier=(20.7011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.65636,0.0656573,-5.43309e-05,2.3989e-08,-4.27542e-12,69567.2,27.5351], Tmin=(100,'K'), Tmax=(1344.25,'K')), NASAPolynomial(coeffs=[14.2856,0.0251019,-9.07671e-06,1.54569e-09,-1.01491e-13,65903,-42.2501], Tmin=(1344.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=C[C]=C[CH]C(=C)C(28087)',
    structure = SMILES('[CH2]C(C)=CC=[C]C=C'),
    E0 = (387.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0777217,0.0748846,-3.83488e-05,-8.11927e-09,1.02403e-11,46776.7,28.4254], Tmin=(100,'K'), Tmax=(952.81,'K')), NASAPolynomial(coeffs=[17.1905,0.0299611,-1.00028e-05,1.69782e-09,-1.15229e-13,42293.8,-59.7187], Tmin=(952.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CC(=C)[CH]C(26578)',
    structure = SMILES('[CH2]C(C=[C]C=C)=CC'),
    E0 = (422.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20964,0.0813771,-5.33489e-05,4.86344e-09,6.0616e-12,50985.3,28.0233], Tmin=(100,'K'), Tmax=(965.094,'K')), NASAPolynomial(coeffs=[18.7617,0.0279536,-9.49238e-06,1.631e-09,-1.11448e-13,46149.6,-68.9099], Tmin=(965.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC(C)C1=C(28088)',
    structure = SMILES('C=CC1=CC(C)C1=C'),
    E0 = (257.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865079,0.0514728,2.27911e-05,-6.54477e-08,2.84718e-11,31041.7,21.0026], Tmin=(100,'K'), Tmax=(975.992,'K')), NASAPolynomial(coeffs=[15.9895,0.0315931,-1.13688e-05,2.08898e-09,-1.503e-13,26084,-61.8704], Tmin=(975.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
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
    label = 'C=C[C]=C[CH]C(26604)',
    structure = SMILES('[CH2]C=[C]C=CC'),
    E0 = (338.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180],'cm^-1')),
        HinderedRotor(inertia=(2.72871,'amu*angstrom^2'), symmetry=1, barrier=(62.7384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968464,'amu*angstrom^2'), symmetry=1, barrier=(22.2669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7301,'amu*angstrom^2'), symmetry=1, barrier=(62.7705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53745,0.0463667,-1.1901e-05,-1.57978e-08,9.31935e-12,40858.4,21.6], Tmin=(100,'K'), Tmax=(968.725,'K')), NASAPolynomial(coeffs=[10.6019,0.0268685,-9.47277e-06,1.63757e-09,-1.11014e-13,38260.9,-26.1849], Tmin=(968.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]C1C(=C)C1C(28089)',
    structure = SMILES('[CH2]C=[C]C1C(=C)C1C'),
    E0 = (580.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281391,0.0704847,-3.57132e-05,-2.51723e-09,5.61998e-12,69926.7,26.3372], Tmin=(100,'K'), Tmax=(1048.18,'K')), NASAPolynomial(coeffs=[15.6153,0.0332581,-1.29068e-05,2.35498e-09,-1.63834e-13,65542.6,-53.9404], Tmin=(1048.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]C(C)C=C=[C]C(28090)',
    structure = SMILES('C=[C]C(C)[CH]C#CC'),
    E0 = (539.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.503788,0.0717872,-4.84835e-05,1.16407e-08,1.88608e-12,65068.2,31.1657], Tmin=(100,'K'), Tmax=(911.165,'K')), NASAPolynomial(coeffs=[11.9765,0.0355375,-1.20448e-05,1.98157e-09,-1.28563e-13,62391.6,-26.3315], Tmin=(911.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)[C]=C=CC(28091)',
    structure = SMILES('C=[C]C(C)C#C[CH]C'),
    E0 = (516.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247493,0.080223,-7.29872e-05,3.84611e-08,-8.41072e-12,62210.7,30.0829], Tmin=(100,'K'), Tmax=(1092.75,'K')), NASAPolynomial(coeffs=[12.014,0.0371507,-1.38609e-05,2.38825e-09,-1.57732e-13,59639.2,-27.7271], Tmin=(1092.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=[C][C](C)C=C=CC(28092)',
    structure = SMILES('C=[C]C(C)=C[C]=CC'),
    E0 = (470.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.253269,0.0900129,-8.95961e-05,4.92993e-08,-1.09542e-11,56692.3,27.4963], Tmin=(100,'K'), Tmax=(1090.65,'K')), NASAPolynomial(coeffs=[15.0491,0.0338909,-1.24098e-05,2.1186e-09,-1.39349e-13,53354.4,-47.6564], Tmin=(1090.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([C]=C)C=C=CC(27270)',
    structure = SMILES('[CH2]C([C]=C)C=C=CC'),
    E0 = (602.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.40376,'amu*angstrom^2'), symmetry=1, barrier=(9.28323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403483,'amu*angstrom^2'), symmetry=1, barrier=(9.27687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403044,'amu*angstrom^2'), symmetry=1, barrier=(9.26678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785712,'amu*angstrom^2'), symmetry=1, barrier=(18.0651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292918,0.0863241,-9.77794e-05,6.95025e-08,-2.10489e-11,72612,30.6472], Tmin=(100,'K'), Tmax=(791.829,'K')), NASAPolynomial(coeffs=[8.46935,0.0450194,-1.95325e-05,3.62267e-09,-2.48617e-13,71317.2,-6.89075], Tmin=(791.829,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(C)C=[C][CH]C(26416)',
    structure = SMILES('C#CC(C)C=[C][CH]C'),
    E0 = (517.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347132,0.0794049,-6.63539e-05,3.02794e-08,-5.77196e-12,62314.6,26.9399], Tmin=(100,'K'), Tmax=(1225.55,'K')), NASAPolynomial(coeffs=[12.9023,0.0384269,-1.61992e-05,2.99662e-09,-2.06539e-13,59237.2,-36.1851], Tmin=(1225.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C1[CH]C(C)C1=C(28093)',
    structure = SMILES('[CH2]C=C1[CH]C(C)C1=C'),
    E0 = (377.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19133,0.0406398,5.67745e-05,-1.01144e-07,4.1188e-11,45528.8,24.3069], Tmin=(100,'K'), Tmax=(964.532,'K')), NASAPolynomial(coeffs=[15.9418,0.0318532,-1.10276e-05,2.0281e-09,-1.48125e-13,40246.6,-58.9547], Tmin=(964.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=CC=CCJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=[C]C(C)C1[C]=CC1(28094)',
    structure = SMILES('C=[C]C(C)C1[C]=CC1'),
    E0 = (651.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.829767,0.0544658,8.12985e-06,-4.70351e-08,2.11218e-11,78486.5,29.8954], Tmin=(100,'K'), Tmax=(997.763,'K')), NASAPolynomial(coeffs=[15.1089,0.0326764,-1.24149e-05,2.30654e-09,-1.64833e-13,73872.2,-47.8052], Tmin=(997.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C=C=C(C)C=C=CC(28095)',
    structure = SMILES('C=C=C(C)C=C=CC'),
    E0 = (301.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.165,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0673712,0.081858,-6.89193e-05,3.07885e-08,-5.61863e-12,36402.9,26.2374], Tmin=(100,'K'), Tmax=(1296.74,'K')), NASAPolynomial(coeffs=[15.461,0.0343734,-1.39912e-05,2.54917e-09,-1.74274e-13,32410.6,-52.0279], Tmin=(1296.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH]=C=CC(C)[C]=C(28096)',
    structure = SMILES('C#C[CH]C(C)[C]=C'),
    E0 = (582.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2175,525,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,1685,370,330.332,330.333],'cm^-1')),
        HinderedRotor(inertia=(0.00154487,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159566,'amu*angstrom^2'), symmetry=1, barrier=(12.356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.948403,'amu*angstrom^2'), symmetry=1, barrier=(73.439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.948389,'amu*angstrom^2'), symmetry=1, barrier=(73.4393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.1384,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883625,0.0640722,-5.08247e-05,1.84394e-08,-1.17582e-12,70127.1,26.3806], Tmin=(100,'K'), Tmax=(929.697,'K')), NASAPolynomial(coeffs=[12.0634,0.0271316,-9.22947e-06,1.52396e-09,-9.91679e-14,67566,-29.3342], Tmin=(929.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(Cds_S)"""),
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
    E0 = (543.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (620.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (675.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (715.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (802.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (672.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (654.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (617.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (644.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (652.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (779.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (758.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (658.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (695.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (685.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (736.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (814.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (969.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (879.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (758.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (789.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (812.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (681.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (604.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (636.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (583.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (997.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (688.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (638.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (552.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (740.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (590.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (751.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (781.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (619.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (655.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (756.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (670.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (662.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (561.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (963.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['CH3CHCCH2(18175)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['[CH2]C1[C]=CC(C)C1=C(28074)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.70833e+09,'s^-1'), n=0.576561, Ea=(76.2219,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C=C(C)C=[C]C=C(28075)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(77.9769,'m^3/(mol*s)'), n=1.64, Ea=(16.5268,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeCs_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=[C]C(C)C#CC=C(28076)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=[C]C(C)C=C=C=C(28077)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C#CC(C)C=[C]C=C(28078)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', 'C=C=CC=[C]C=C(28079)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00970044,'m^3/(mol*s)'), n=2.41, Ea=(33.193,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Ca;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C(4699)', 'CH3CHCCH2(18175)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00472174,'m^3/(mol*s)'), n=2.41, Ea=(20.1294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C][CH]C(18176)', 'CH2CHCCH(26391)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C]=C[C](C)C=C(27245)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.03e+09,'s^-1'), n=1.31, Ea=(263.174,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 156 used for R2H_S;C_rad_out_OneDe/Cs;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C)[C]=CC=C(28080)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C=CC(C)[C]=C(28081)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CC(C)C=[C]C=C(27252)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['[CH2]C(C=C)C=[C]C=C(26570)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C[C]=[C]C(C)C=C(27248)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C][C](C)C=CC=C(28082)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=CC=CC(C)[C]=C(28083)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([C]=C)C=CC=C(27251)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C#CC(C)[CH]C=C[CH2](27701)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.41973e+10,'s^-1'), n=1.543, Ea=(335.233,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleDe_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C][C]=CC(C)C=C(27254)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_doubleC] for rate rule [R5HJ_3;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['[CH]=C[C]=CC(C)C=C(27255)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0.627, Ea=(245.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;Cd_H_out_singleH] for rate rule [R6HJ_3;Cd_rad_out_Cd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]C=C(4699)', 'C=[C][CH]C(18176)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C]C(C)C=C1[CH]C1(28084)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C1C[CH][C]=CC1C(28085)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C=C(C)C=CC=C(28086)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.4e+09,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C=C=CC(C)C=C(27261)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(S)(23)', 'C=[C]CC=[C]C=C(27963)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C[C]=C[CH]C(=C)C(28087)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.95888e+10,'s^-1'), n=0.7315, Ea=(144.62,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCs(-HC)CJ;CJ;CH3] + [cCs(-HC)CJ;CdsJ;C] for rate rule [cCs(-HC)CJ;CdsJ;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C[C]=CC(=C)[CH]C(26578)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 5 used for cCs(-HC)CJ;CdsJ;C
Exact match found for rate rule [cCs(-HC)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=CC1=CC(C)C1=C(28088)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H2CC(41)', 'C=C[C]=C[CH]C(26604)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['[CH2]C=[C]C1C(=C)C1C(28089)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(46.3011,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C]C(C)C=C=[C]C(28090)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C]C(C)[C]=C=CC(28091)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C][C](C)C=C=CC(28092)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.21176e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([C]=C)C=C=CC(27270)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C#CC(C)C=[C][CH]C(26416)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(291.349,'s^-1'), n=2.95077, Ea=(212.39,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_singleH] for rate rule [R7HJ_5;C_rad_out_2H;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['[CH2]C=C1[CH]C(C)C1=C(28093)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=[C]C(C)C1[C]=CC1(28094)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C]C(C)C=[C]C=C(26582)'],
    products = ['C=C=C(C)C=C=CC(28095)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CH2(19)', '[CH]=C=CC(C)[C]=C(28096)'],
    products = ['C=[C]C(C)C=[C]C=C(26582)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4554',
    isomers = [
        'C=[C]C(C)C=[C]C=C(26582)',
    ],
    reactants = [
        ('CH3CHCCH2(18175)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4554',
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

