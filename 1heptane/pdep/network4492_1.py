species(
    label = 'C=C[C]=CC[CH]O(26524)',
    structure = SMILES('C=C[C]=CC[CH]O'),
    E0 = (254.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,193.955,193.989,193.999],'cm^-1')),
        HinderedRotor(inertia=(0.522666,'amu*angstrom^2'), symmetry=1, barrier=(13.9536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522574,'amu*angstrom^2'), symmetry=1, barrier=(13.9535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.47971,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522428,'amu*angstrom^2'), symmetry=1, barrier=(13.9534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430935,0.077993,-8.62319e-05,5.20268e-08,-1.25918e-11,30764.2,27.3759], Tmin=(100,'K'), Tmax=(1006.21,'K')), NASAPolynomial(coeffs=[13.2122,0.0271836,-1.04881e-05,1.84252e-09,-1.23158e-13,28192.1,-34.3653], Tmin=(1006.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCsJOH) + radical(C=CJC=C)"""),
)

species(
    label = 'CH2CHOH(42)',
    structure = SMILES('C=CO'),
    E0 = (-138.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.72808,'amu*angstrom^2'), symmetry=1, barrier=(39.7321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.11,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-138.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = '[CH2]C1[C]=CCC1O(27792)',
    structure = SMILES('[CH2]C1[C]=CCC1O'),
    E0 = (280.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53571,0.0395906,2.39899e-05,-6.08028e-08,2.67678e-11,33821.9,24.7362], Tmin=(100,'K'), Tmax=(951.091,'K')), NASAPolynomial(coeffs=[14.0429,0.0232488,-7.42337e-06,1.30114e-09,-9.30772e-14,29802.9,-43.5982], Tmin=(951.091,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl)"""),
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
    label = 'C=C[C]=CCC=O(27818)',
    structure = SMILES('C=C[C]=CCC=O'),
    E0 = (157.629,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05087,'amu*angstrom^2'), symmetry=1, barrier=(24.1616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05082,'amu*angstrom^2'), symmetry=1, barrier=(24.1605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05052,'amu*angstrom^2'), symmetry=1, barrier=(24.1536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14085,0.052814,-2.28167e-05,-5.61669e-09,4.8368e-12,19069.7,24.168], Tmin=(100,'K'), Tmax=(1123.45,'K')), NASAPolynomial(coeffs=[13.9442,0.0257035,-1.12871e-05,2.17957e-09,-1.55476e-13,15027,-44.2801], Tmin=(1123.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CC=CO(27819)',
    structure = SMILES('C=C[C]=CC=CO'),
    E0 = (135.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,3615,1277.5,1000,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28045,'amu*angstrom^2'), symmetry=1, barrier=(29.44,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27915,'amu*angstrom^2'), symmetry=1, barrier=(29.4101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27854,'amu*angstrom^2'), symmetry=1, barrier=(29.3961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.272626,0.0636521,-1.35483e-05,-5.12935e-08,3.17303e-11,16496,23.4682], Tmin=(100,'K'), Tmax=(896.394,'K')), NASAPolynomial(coeffs=[25.3262,0.00261242,3.65767e-06,-9.21272e-10,6.43408e-14,9965.2,-106.035], Tmin=(896.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC#CC[CH]O(27820)',
    structure = SMILES('C=CC#CC[CH]O'),
    E0 = (226.504,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,376.71,376.71,376.711],'cm^-1')),
        HinderedRotor(inertia=(0.254645,'amu*angstrom^2'), symmetry=1, barrier=(25.6434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127179,'amu*angstrom^2'), symmetry=1, barrier=(12.8072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127178,'amu*angstrom^2'), symmetry=1, barrier=(12.8072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762983,'amu*angstrom^2'), symmetry=1, barrier=(76.8349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76391,0.0657047,-6.13489e-05,3.06452e-08,-6.13007e-12,27363.3,27.3161], Tmin=(100,'K'), Tmax=(1210.8,'K')), NASAPolynomial(coeffs=[13.741,0.0228341,-8.23935e-06,1.40348e-09,-9.24499e-14,24220.7,-37.7733], Tmin=(1210.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CCsJOH)"""),
)

species(
    label = 'C=C=C=CC[CH]O(27821)',
    structure = SMILES('C=C=C=CC[CH]O'),
    E0 = (285.106,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3025,407.5,1350,352.5,3615,1277.5,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.721843,'amu*angstrom^2'), symmetry=1, barrier=(16.5966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.721796,'amu*angstrom^2'), symmetry=1, barrier=(16.5955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72117,'amu*angstrom^2'), symmetry=1, barrier=(16.5811,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.739971,0.0755133,-9.13459e-05,6.1998e-08,-1.72016e-11,34404.4,25.7628], Tmin=(100,'K'), Tmax=(873.046,'K')), NASAPolynomial(coeffs=[10.6299,0.0302001,-1.34906e-05,2.54546e-09,-1.76779e-13,32677.5,-20.6073], Tmin=(873.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(CCsJOH)"""),
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
    label = '[CH2][CH]O(284)',
    structure = SMILES('[CH2][CH]O'),
    E0 = (143.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.217851,'amu*angstrom^2'), symmetry=1, barrier=(5.00882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197382,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84769,0.0229973,-1.85068e-05,7.77211e-09,-1.30725e-12,17300.5,13.0245], Tmin=(100,'K'), Tmax=(1418.34,'K')), NASAPolynomial(coeffs=[7.97636,0.00853337,-3.21015e-06,5.82141e-10,-3.99268e-14,15845.7,-13.5107], Tmin=(1418.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = 'C=C[C]=CCC[O](27822)',
    structure = SMILES('C=C[C]=CCC[O]'),
    E0 = (300.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21068,0.0646208,-6.0073e-05,3.44826e-08,-8.67336e-12,36194.1,25.8981], Tmin=(100,'K'), Tmax=(928.851,'K')), NASAPolynomial(coeffs=[7.51013,0.037493,-1.62645e-05,3.03985e-09,-2.10573e-13,35023.9,-4.02809], Tmin=(928.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=C[CH]CO(27823)',
    structure = SMILES('[CH2]C=[C]C=CCO'),
    E0 = (185.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.776434,0.0633631,-5.08699e-05,2.24993e-08,-4.06284e-12,22467.2,28.0893], Tmin=(100,'K'), Tmax=(1322,'K')), NASAPolynomial(coeffs=[12.8845,0.026728,-9.30259e-06,1.53768e-09,-9.88991e-14,19265.8,-33.705], Tmin=(1322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC=[C]C[CH]O(27824)',
    structure = SMILES('C=CC=[C]C[CH]O'),
    E0 = (293.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,275.148,275.157,275.158],'cm^-1')),
        HinderedRotor(inertia=(0.239966,'amu*angstrom^2'), symmetry=1, barrier=(12.8919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239952,'amu*angstrom^2'), symmetry=1, barrier=(12.8914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239947,'amu*angstrom^2'), symmetry=1, barrier=(12.8916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.239918,'amu*angstrom^2'), symmetry=1, barrier=(12.8911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.544967,0.0752313,-7.90402e-05,4.46221e-08,-1.01377e-11,35432.5,27.8439], Tmin=(100,'K'), Tmax=(1065.97,'K')), NASAPolynomial(coeffs=[13.4596,0.0267695,-1.08456e-05,1.97221e-09,-1.34991e-13,32679.2,-35.2864], Tmin=(1065.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CC[CH]O(27825)',
    structure = SMILES('C=[C]C=CC[CH]O'),
    E0 = (254.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,193.955,193.989,193.999],'cm^-1')),
        HinderedRotor(inertia=(0.522666,'amu*angstrom^2'), symmetry=1, barrier=(13.9536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522574,'amu*angstrom^2'), symmetry=1, barrier=(13.9535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.47971,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522428,'amu*angstrom^2'), symmetry=1, barrier=(13.9534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430935,0.077993,-8.62319e-05,5.20268e-08,-1.25918e-11,30764.2,27.3759], Tmin=(100,'K'), Tmax=(1006.21,'K')), NASAPolynomial(coeffs=[13.2122,0.0271836,-1.04881e-05,1.84252e-09,-1.23158e-13,28192.1,-34.3653], Tmin=(1006.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCsJOH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=[C]CCO(27826)',
    structure = SMILES('[CH2][CH]C#CCCO'),
    E0 = (276.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2100,2250,500,550,3615,1277.5,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856979,0.0659794,-6.51422e-05,3.80998e-08,-9.03939e-12,33375.6,30.7778], Tmin=(100,'K'), Tmax=(1089.1,'K')), NASAPolynomial(coeffs=[10.542,0.0277271,-8.76458e-06,1.32879e-09,-7.97531e-14,31425,-16.0432], Tmin=(1089.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH2]C=C[CH]C=CO(23616)',
    structure = SMILES('[CH2]C=C[CH]C=CO'),
    E0 = (97.1557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861405,0.0446295,4.46072e-05,-1.07119e-07,4.95686e-11,11820.8,25.7058], Tmin=(100,'K'), Tmax=(918.405,'K')), NASAPolynomial(coeffs=[23.9367,0.00763448,1.30669e-06,-3.95217e-10,2.17645e-14,4904.04,-98.2363], Tmin=(918.405,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.1557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=CC=CC[CH]O(27827)',
    structure = SMILES('[CH]=CC=CC[CH]O'),
    E0 = (302.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3615,1277.5,1000,3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.617472,'amu*angstrom^2'), symmetry=1, barrier=(14.1969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.616991,'amu*angstrom^2'), symmetry=1, barrier=(14.1858,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.617834,'amu*angstrom^2'), symmetry=1, barrier=(14.2052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.617282,'amu*angstrom^2'), symmetry=1, barrier=(14.1925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.356497,0.0759278,-7.72275e-05,4.10324e-08,-8.65161e-12,36555.2,28.447], Tmin=(100,'K'), Tmax=(1154.31,'K')), NASAPolynomial(coeffs=[15.5271,0.0233576,-8.91394e-06,1.57837e-09,-1.06705e-13,33052.9,-46.9197], Tmin=(1154.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C[CH]CC=O(20125)',
    structure = SMILES('[CH2]C=C[CH]CC=O'),
    E0 = (192.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,309.265,2254.03],'cm^-1')),
        HinderedRotor(inertia=(0.453707,'amu*angstrom^2'), symmetry=1, barrier=(30.8066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125951,'amu*angstrom^2'), symmetry=1, barrier=(8.5531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.453326,'amu*angstrom^2'), symmetry=1, barrier=(30.8054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0027132,'amu*angstrom^2'), symmetry=1, barrier=(30.8059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05788,0.0569622,-3.72418e-05,1.21239e-08,-1.59598e-12,23237.5,29.5197], Tmin=(100,'K'), Tmax=(1748.47,'K')), NASAPolynomial(coeffs=[14.876,0.0253504,-1.01224e-05,1.78373e-09,-1.17529e-13,18405.4,-44.8658], Tmin=(1748.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCJCC=O)"""),
)

species(
    label = 'C=[C][C]=CCCO(27828)',
    structure = SMILES('[CH2]C#C[CH]CCO'),
    E0 = (226.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924357,0.0664099,-6.41603e-05,3.63512e-08,-8.512e-12,27342.6,28.2236], Tmin=(100,'K'), Tmax=(1027.37,'K')), NASAPolynomial(coeffs=[10.2701,0.0300223,-1.10322e-05,1.87563e-09,-1.22598e-13,25422.3,-17.1165], Tmin=(1027.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(Propargyl)"""),
)

species(
    label = '[CH]=C[C]=CCCO(27829)',
    structure = SMILES('[CH]C=C=CCCO'),
    E0 = (298.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,540,610,2055,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.994934,0.0669244,-5.08825e-05,2.15322e-08,-3.93459e-12,36055.3,27.3412], Tmin=(100,'K'), Tmax=(1240.2,'K')), NASAPolynomial(coeffs=[9.80612,0.0385058,-1.65108e-05,3.05585e-09,-2.10116e-13,33869.8,-17.0645], Tmin=(1240.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O[CH]CC=C1[CH]C1(27830)',
    structure = SMILES('O[CH]CC=C1[CH]C1'),
    E0 = (289.856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910229,0.0586791,-3.20215e-05,-3.31457e-10,4.3186e-12,34980.9,25.1317], Tmin=(100,'K'), Tmax=(1038.85,'K')), NASAPolynomial(coeffs=[14.0437,0.025735,-9.90229e-06,1.80537e-09,-1.25839e-13,31301.1,-43.3075], Tmin=(1038.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(CCsJOH)"""),
)

species(
    label = 'OC1C[CH][C]=CC1(27831)',
    structure = SMILES('OC1C[CH][C]=CC1'),
    E0 = (182.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81546,0.030713,4.93896e-05,-8.29019e-08,3.27923e-11,22068.4,19.3898], Tmin=(100,'K'), Tmax=(981.202,'K')), NASAPolynomial(coeffs=[13.2749,0.0267867,-1.0022e-05,1.90936e-09,-1.40769e-13,17759.8,-46.174], Tmin=(981.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cds_S) + radical(cyclohexene-allyl)"""),
)

species(
    label = 'C=CC=CC=CO(23629)',
    structure = SMILES('C=CC=CC=CO'),
    E0 = (-63.0995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.342679,0.0572292,1.49765e-05,-8.31724e-08,4.29321e-11,-7435.74,23.6184], Tmin=(100,'K'), Tmax=(911.518,'K')), NASAPolynomial(coeffs=[26.7946,0.00264365,3.6095e-06,-8.48075e-10,5.4456e-14,-14812.6,-115.559], Tmin=(911.518,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.0995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC=CCC=O(20142)',
    structure = SMILES('C=CC=CCC=O'),
    E0 = (-41.3668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25991,0.045833,7.57032e-06,-3.97998e-08,1.69904e-11,-4864.2,24.1411], Tmin=(100,'K'), Tmax=(1049.03,'K')), NASAPolynomial(coeffs=[14.962,0.0264616,-1.17388e-05,2.34542e-09,-1.72883e-13,-9547.88,-51.2411], Tmin=(1049.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.3668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C=C=CCCO(27832)',
    structure = SMILES('C=C=C=CCCO'),
    E0 = (104.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02793,0.0666574,-6.01591e-05,3.01908e-08,-6.35638e-12,12711.4,25.1727], Tmin=(100,'K'), Tmax=(1115.82,'K')), NASAPolynomial(coeffs=[10.596,0.0323575,-1.40495e-05,2.64176e-09,-1.83984e-13,10576.2,-22.0364], Tmin=(1115.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH2]C(O)C=[C]C=C(26528)',
    structure = SMILES('[CH2]C(O)C=[C]C=C'),
    E0 = (273.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,237.521,237.832],'cm^-1')),
        HinderedRotor(inertia=(0.439806,'amu*angstrom^2'), symmetry=1, barrier=(17.5378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439616,'amu*angstrom^2'), symmetry=1, barrier=(17.538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438304,'amu*angstrom^2'), symmetry=1, barrier=(17.5428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438059,'amu*angstrom^2'), symmetry=1, barrier=(17.5398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3876.19,'J/mol'), sigma=(6.45323,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.45 K, Pc=32.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.404317,0.071352,-5.97882e-05,2.04165e-08,-7.17957e-13,32974,29.1159], Tmin=(100,'K'), Tmax=(975.943,'K')), NASAPolynomial(coeffs=[16.5018,0.0212801,-7.27465e-06,1.24335e-09,-8.41369e-14,29074.5,-52.0336], Tmin=(975.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CJCO) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CCC1O(27780)',
    structure = SMILES('C=CC1=CCC1O'),
    E0 = (21.8677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23519,0.0397063,4.1413e-05,-9.0466e-08,3.97191e-11,2748.91,23.2278], Tmin=(100,'K'), Tmax=(945.772,'K')), NASAPolynomial(coeffs=[19.6081,0.0150754,-3.69823e-06,6.67519e-10,-5.48568e-14,-3100.11,-76.9353], Tmin=(945.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.8677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = 'HCOH(T)(285)',
    structure = SMILES('[CH]O'),
    E0 = (205.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,403.876,3308.82],'cm^-1')),
        HinderedRotor(inertia=(0.0103144,'amu*angstrom^2'), symmetry=1, barrier=(22.7121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75938,0.0029613,8.90411e-06,-1.35016e-08,5.39816e-12,24775.6,6.76286], Tmin=(100,'K'), Tmax=(940.429,'K')), NASAPolynomial(coeffs=[5.09112,0.00321239,-9.31686e-07,1.59615e-10,-1.15729e-14,24263.5,-0.971], Tmin=(940.429,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HCOH(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C=[C]C=C(26613)',
    structure = SMILES('[CH2]C=[C]C=C'),
    E0 = (374.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(2.07984,'amu*angstrom^2'), symmetry=1, barrier=(47.8195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07746,'amu*angstrom^2'), symmetry=1, barrier=(47.7649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22509,0.0302876,1.05277e-05,-3.68274e-08,1.72727e-11,45167.7,17.3269], Tmin=(100,'K'), Tmax=(923.516,'K')), NASAPolynomial(coeffs=[10.3879,0.0174192,-5.09531e-06,8.16549e-10,-5.51179e-14,42701,-26.5965], Tmin=(923.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]C1CC1O(27813)',
    structure = SMILES('[CH2]C=[C]C1CC1O'),
    E0 = (295.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18441,0.0469305,9.42544e-06,-4.92688e-08,2.31605e-11,35615.3,25.9617], Tmin=(100,'K'), Tmax=(964.221,'K')), NASAPolynomial(coeffs=[16.1575,0.0211919,-7.12295e-06,1.29871e-09,-9.49352e-14,31036.9,-54.4978], Tmin=(964.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=C=CC[CH]O(27833)',
    structure = SMILES('CC#C[CH]C[CH]O'),
    E0 = (250.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.566692,0.0761077,-9.18267e-05,6.37077e-08,-1.74881e-11,30227.2,27.9971], Tmin=(100,'K'), Tmax=(998.178,'K')), NASAPolynomial(coeffs=[10.6605,0.0284645,-9.42073e-06,1.44967e-09,-8.67996e-14,28570.5,-18.8862], Tmin=(998.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(Sec_Propargyl)"""),
)

species(
    label = 'CC=C=[C]C[CH]O(27834)',
    structure = SMILES('C[CH]C#CC[CH]O'),
    E0 = (251.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986571,0.0682221,-6.22652e-05,2.10735e-08,3.74557e-12,30366.2,27.4194], Tmin=(100,'K'), Tmax=(707.696,'K')), NASAPolynomial(coeffs=[10.8538,0.0277505,-8.91132e-06,1.36065e-09,-8.19783e-14,28586.5,-19.4796], Tmin=(707.696,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.589,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
)

species(
    label = 'CC=[C][CH]C=CO(26217)',
    structure = SMILES('CC=[C][CH]C=CO'),
    E0 = (183.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83044,0.0517212,1.11494e-05,-6.40823e-08,3.22735e-11,22200.5,26.0093], Tmin=(100,'K'), Tmax=(922.293,'K')), NASAPolynomial(coeffs=[20.4771,0.0131485,-1.96362e-06,2.2116e-10,-1.79061e-14,16593.1,-77.9382], Tmin=(922.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = 'CC=[C][CH]CC=O(24868)',
    structure = SMILES('C[CH][C]=CCC=O'),
    E0 = (227.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40417,0.0481489,-1.20072e-05,-1.0811e-08,4.98103e-12,27419.8,25.0914], Tmin=(100,'K'), Tmax=(1248.56,'K')), NASAPolynomial(coeffs=[12.5252,0.0315227,-1.48615e-05,2.90254e-09,-2.05555e-13,23161.7,-36.9609], Tmin=(1248.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C1[CH]CC1O(27835)',
    structure = SMILES('[CH2]C=C1[CH]CC1O'),
    E0 = (184.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[10.4565,-0.0156326,0.000136662,-1.31131e-07,3.04373e-11,21803,-16.4712], Tmin=(100,'K'), Tmax=(1695.29,'K')), NASAPolynomial(coeffs=[72.2709,0.0296871,-7.25854e-05,1.76716e-08,-1.31532e-12,-26626.8,-428.341], Tmin=(1695.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'O[CH]CC1[C]=CC1(27836)',
    structure = SMILES('O[CH]CC1[C]=CC1'),
    E0 = (361.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.955633,0.05798,-3.20781e-05,4.77698e-10,3.92247e-12,43642.8,26.5051], Tmin=(100,'K'), Tmax=(1041.19,'K')), NASAPolynomial(coeffs=[13.7405,0.0257108,-9.86024e-06,1.79234e-09,-1.24617e-13,40067.3,-40.076], Tmin=(1041.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(CCsJOH)"""),
)

species(
    label = 'CC=C=CC=CO(26219)',
    structure = SMILES('CC=C=CC=CO'),
    E0 = (-10.3183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35892,0.0638707,-1.86754e-05,-3.67018e-08,2.32702e-11,-1094.9,23.5891], Tmin=(100,'K'), Tmax=(924.187,'K')), NASAPolynomial(coeffs=[22.218,0.0112934,-1.55934e-06,1.61829e-10,-1.35158e-14,-6930.27,-89.8562], Tmin=(924.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.3183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'CC=C=CCC=O(24870)',
    structure = SMILES('CC=C=CCC=O'),
    E0 = (11.4144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28232,0.0519897,-2.23989e-05,-1.05919e-09,1.94195e-12,1476.92,24.1167], Tmin=(100,'K'), Tmax=(1364.33,'K')), NASAPolynomial(coeffs=[13.8745,0.0297125,-1.40034e-05,2.7041e-09,-1.8895e-13,-3321.68,-45.5388], Tmin=(1364.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.4144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cdd-CdsCds)"""),
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
    label = 'C#C[CH]C[CH]O(26949)',
    structure = SMILES('C#C[CH]C[CH]O'),
    E0 = (292.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,3615,1277.5,1000,3000,3050,390,425,1340,1360,335,370,303.289,305.12],'cm^-1')),
        HinderedRotor(inertia=(0.00388476,'amu*angstrom^2'), symmetry=1, barrier=(10.6425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161691,'amu*angstrom^2'), symmetry=1, barrier=(10.6485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05523,'amu*angstrom^2'), symmetry=1, barrier=(68.8377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0553,'amu*angstrom^2'), symmetry=1, barrier=(68.8341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99162,0.0678312,-9.20601e-05,6.75543e-08,-1.91876e-11,35284.2,23.0522], Tmin=(100,'K'), Tmax=(981.921,'K')), NASAPolynomial(coeffs=[10.6377,0.020251,-6.71881e-06,1.01928e-09,-5.96874e-14,33789.3,-21.2748], Tmin=(981.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCsJOH)"""),
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
    E0 = (254.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (280.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (405.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (355.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (450.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (513.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (344.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (427.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (384.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (371.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (490.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (469.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (405.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (447.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (525.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (481.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (473.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (393.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (595.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (392.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (281.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (354.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (293.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (263.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (431.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (262.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (580.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (319.421,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (462.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (492.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (330.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (335.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (379.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (373.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (289.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (279.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (674.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['CH2CHOH(42)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['[CH2]C1[C]=CCC1O(27792)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(187000,'s^-1'), n=1.48, Ea=(25.6475,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 20.4 to 25.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C[C]=CCC=O(27818)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2818 used for Od_CO-CsH;HJ
Exact match found for rate rule [Od_CO-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C[C]=CC=CO(27819)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC#CC[CH]O(27820)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=C=C=CC[CH]O(27821)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2CHOH(42)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00668046,'m^3/(mol*s)'), n=2.5095, Ea=(31.5264,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-OsH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]O(284)', 'CH2CHCCH(26391)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=C[C]=CCC[O](27822)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4500,'s^-1'), n=2.62, Ea=(129.286,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 322 used for R2H_S;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=C[C]=C[CH]CO(27823)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.58526e-05,'s^-1'), n=4.94167, Ea=(116.966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_1H;Cs_H_out_H/OneDe] + [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_1H] for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=CC=[C]C[CH]O(27824)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C=CC[CH]O(27825)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=[C]CCO(27826)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00524938,'s^-1'), n=4.55, Ea=(129.286,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeO] + [R3H_SS_Cs;Cd_rad_out;Cs_H_out_1H] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['[CH2]C=C[CH]C=CO(23616)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.56742e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC=CC[CH]O(27827)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['[CH2]C=C[CH]CC=O(20125)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=[C][C]=CCCO(27828)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.68838e+09,'s^-1'), n=0.99528, Ea=(219.261,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_1H;Cd_H_out_doubleC] for rate rule [R5HJ_3;C_rad_out_H/NonDeO;Cd_H_out_doubleC]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C[C]=CCCO(27829)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.70099e+07,'s^-1'), n=1.485, Ea=(94.7468,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]O(284)', '[CH]=[C]C=C(4699)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['O[CH]CC=C1[CH]C1(27830)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['OC1C[CH][C]=CC1(27831)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHNd] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=CC=CC=CO(23629)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.94203e+09,'s^-1'), n=0.37, Ea=(99.6901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad_NDe] + [R3radExo;Y_rad;XH_Rrad_NDe] for rate rule [R3radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=CC=CCC=O(20142)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=C=C=CCCO(27832)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)C=[C]C=C(26528)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C=CC1=CCC1O(27780)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeO;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HCOH(T)(285)', '[CH2]C=[C]C=C(26613)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.13464e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['[CH2]C=[C]C1CC1O(27813)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(64.7047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['C[C]=C=CC[CH]O(27833)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['CC=C=[C]C[CH]O(27834)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['CC=[C][CH]C=CO(26217)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.42353e+06,'s^-1'), n=1.41298, Ea=(75.8094,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;C_rad_out_2H;XH_out] for rate rule [R5H_SMMS;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['CC=[C][CH]CC=O(24868)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(60975.7,'s^-1'), n=1.58648, Ea=(80.9836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7HJ_5;C_rad_out_2H;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['[CH2]C=C1[CH]CC1O(27835)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.25757e+07,'s^-1'), n=1.165, Ea=(125.102,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHNd] for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['O[CH]CC1[C]=CC1(27836)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['CC=C=CC=CO(26219)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C[C]=CC[CH]O(26524)'],
    products = ['CC=C=CCC=O(24870)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(19)', 'C#C[CH]C[CH]O(26949)'],
    products = ['C=C[C]=CC[CH]O(26524)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4492',
    isomers = [
        'C=C[C]=CC[CH]O(26524)',
    ],
    reactants = [
        ('CH2CHOH(42)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4492',
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

