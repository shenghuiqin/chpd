species(
    label = '[O]OOC1=CC=[C]C=CC1(6853)',
    structure = SMILES('[O]OOC1=CC=[C]C=CC1'),
    E0 = (462.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18526,0.095679,-9.75964e-05,4.86507e-08,-9.21629e-12,55818.2,34.1704], Tmin=(100,'K'), Tmax=(1422.27,'K')), NASAPolynomial(coeffs=[24.7041,0.0146237,-3.41674e-06,4.30079e-10,-2.39269e-14,49287.7,-96.9194], Tmin=(1422.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)(2)',
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
    label = 'C7H6O(490)(489)',
    structure = SMILES('C1=CC2=CC=C(C1)O2'),
    E0 = (41.2019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4154,'J/mol'), sigma=(6.51034,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=648.85 K, Pc=34.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20282,0.0266771,3.98414e-05,-5.77521e-08,1.86939e-11,5030.4,17.6585], Tmin=(100,'K'), Tmax=(1194.95,'K')), NASAPolynomial(coeffs=[12.8261,0.0324323,-1.92456e-05,4.14704e-09,-3.09659e-13,-458.206,-47.8272], Tmin=(1194.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.2019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s3_5_6_ane) - ring(Oxane) - ring(Tetrahydrofuran) + ring(36dihydro2hpyran) + ring(Furan)"""),
)

species(
    label = '[C]1=C[CH]C2(CC=C1)OOO2(7990)',
    structure = SMILES('[C]1=C[CH]C2(CC=C1)OOO2'),
    E0 = (290.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.49673,0.163389,-0.000297656,2.79424e-07,-1.00874e-10,35208.4,-15.6188], Tmin=(100,'K'), Tmax=(823.422,'K')), NASAPolynomial(coeffs=[11.4847,0.056378,-3.15048e-05,6.28422e-09,-4.4041e-13,34231.1,-72.3072], Tmin=(823.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cyclobutane) + ring(1,3-Cycloheptadiene) + ring(Cyclobutane) + radical(C=CJC=C) + radical(C=CCJCO)"""),
)

species(
    label = '[C]1=CC=C2C[CH]C1OOO2(7991)',
    structure = SMILES('[C]1=CC=C2C[CH]C1OOO2'),
    E0 = (463.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.63263,0.160285,-0.000345506,3.67357e-07,-1.42151e-10,55921.7,-8.32414], Tmin=(100,'K'), Tmax=(848.937,'K')), NASAPolynomial(coeffs=[-2.84083,0.0779603,-4.45265e-05,8.87058e-09,-6.17185e-13,59298.5,15.9871], Tmin=(848.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cycloheptane) + ring(1,3-Cycloheptadiene) + ring(Cycloheptane) + radical(CCJCOOH) + radical(Cds_S)"""),
)

species(
    label = 'H(3)(3)',
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
    label = '[O]OOC1=CC=C=C=CC1(7992)',
    structure = SMILES('[O]OOC1=CC=C=C=CC1'),
    E0 = (584.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.492068,0.0635485,-2.55254e-05,-1.8058e-08,1.28279e-11,70422.8,30.0004], Tmin=(100,'K'), Tmax=(986.278,'K')), NASAPolynomial(coeffs=[18.4975,0.0211878,-7.7352e-06,1.43931e-09,-1.04508e-13,65379.7,-64.1774], Tmin=(986.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + ring(1_2_3_cycloheptatriene) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1=CC#CC=CC1(7993)',
    structure = SMILES('[O]OOC1=CC#CC=CC1'),
    E0 = (530.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.826424,0.0518361,8.3335e-06,-4.75697e-08,2.05622e-11,63900.2,28.529], Tmin=(100,'K'), Tmax=(1051.23,'K')), NASAPolynomial(coeffs=[18.6822,0.0249865,-1.199e-05,2.50419e-09,-1.89416e-13,57875.5,-69.3064], Tmin=(1051.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-Ct(Cds-Cds)) + ring(cycloheptyne) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1=CC=C[C]=CC1(6903)',
    structure = SMILES('[O]OOC1=CC=C[C]=CC1'),
    E0 = (462.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18525,0.0956789,-9.7596e-05,4.86504e-08,-9.21619e-12,55818.2,34.1703], Tmin=(100,'K'), Tmax=(1422.28,'K')), NASAPolynomial(coeffs=[24.7039,0.014624,-3.4169e-06,4.30115e-10,-2.39297e-14,49287.8,-96.9182], Tmin=(1422.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]OOC1=C[C]=CC=CC1(7994)',
    structure = SMILES('[O]OOC1=C[C]=CC=CC1'),
    E0 = (462.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18525,0.0956789,-9.7596e-05,4.86504e-08,-9.21619e-12,55818.2,34.1703], Tmin=(100,'K'), Tmax=(1422.28,'K')), NASAPolynomial(coeffs=[24.7039,0.014624,-3.4169e-06,4.30115e-10,-2.39297e-14,49287.8,-96.9182], Tmin=(1422.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]OOC1=CC=CC=[C]C1(7995)',
    structure = SMILES('[O]OOC1=CC=CC=[C]C1'),
    E0 = (501.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.509227,0.0864336,-6.84032e-05,1.3682e-08,4.57448e-12,60461.9,32.6126], Tmin=(100,'K'), Tmax=(966.308,'K')), NASAPolynomial(coeffs=[23.1049,0.0173119,-5.54485e-06,9.74505e-10,-6.99252e-14,54561.6,-87.4182], Tmin=(966.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1=[C]C=CC=CC1(7996)',
    structure = SMILES('[O]OOC1=[C]C=CC=CC1'),
    E0 = (462.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18525,0.0956789,-9.7596e-05,4.86504e-08,-9.21619e-12,55818.2,34.1703], Tmin=(100,'K'), Tmax=(1422.28,'K')), NASAPolynomial(coeffs=[24.7039,0.014624,-3.4169e-06,4.30115e-10,-2.39297e-14,49287.8,-96.9182], Tmin=(1422.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]OOC1[CH]C=CC=CC=1(6847)',
    structure = SMILES('[O]OOC1C=C[CH]C=CC=1'),
    E0 = (365.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.197498,0.073427,-2.15188e-05,-4.10472e-08,2.57462e-11,44165.5,29.1833], Tmin=(100,'K'), Tmax=(934.108,'K')), NASAPolynomial(coeffs=[24.8968,0.0138447,-2.71972e-06,4.03936e-10,-3.22248e-14,37388.7,-101.352], Tmin=(934.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CCJC=C) + radical(ROOJ)"""),
)

species(
    label = 'OOOC1C=C[CH][C]=CC=1(6890)',
    structure = SMILES('OOOC1C=C[C]=C[CH]C=1'),
    E0 = (412.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.561521,0.0838166,-4.67596e-05,-1.64359e-08,1.73225e-11,49828.1,30.0211], Tmin=(100,'K'), Tmax=(933.624,'K')), NASAPolynomial(coeffs=[25.2079,0.0151405,-3.46679e-06,5.24581e-10,-3.85387e-14,43197.6,-102.271], Tmin=(933.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CCJC=C)"""),
)

species(
    label = 'OOOC1=[C]C=[C]C=CC1(7997)',
    structure = SMILES('OOOC1=[C]C=[C]C=CC1'),
    E0 = (509.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.917405,0.0987547,-9.7893e-05,4.17982e-08,-4.59059e-12,61453.2,32.0389], Tmin=(100,'K'), Tmax=(950.7,'K')), NASAPolynomial(coeffs=[23.0607,0.0192025,-6.03723e-06,9.89484e-10,-6.63816e-14,55929.8,-87.5], Tmin=(950.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'OOOC1=CC=[C]C=[C]C1(7998)',
    structure = SMILES('OOOC1=CC=[C]C=[C]C1'),
    E0 = (548.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23833,0.101054,-0.000108057,5.6371e-08,-1.12728e-11,66140.4,34.0714], Tmin=(100,'K'), Tmax=(1273.41,'K')), NASAPolynomial(coeffs=[25.0948,0.0158351,-4.72659e-06,7.31182e-10,-4.64108e-14,59636.7,-98.5386], Tmin=(1273.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'OOOC1=C[C]=[C]C=CC1(7999)',
    structure = SMILES('OOOC1=C[C]=[C]C=CC1'),
    E0 = (509.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.917405,0.0987547,-9.7893e-05,4.17982e-08,-4.59059e-12,61453.2,32.0389], Tmin=(100,'K'), Tmax=(950.7,'K')), NASAPolynomial(coeffs=[23.0607,0.0192025,-6.03723e-06,9.89484e-10,-6.63816e-14,55929.8,-87.5], Tmin=(950.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'OOOC1=CC=[C][C]=CC1(8000)',
    structure = SMILES('OOOC1=CC=[C][C]=CC1'),
    E0 = (509.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.917405,0.0987547,-9.7893e-05,4.17982e-08,-4.59059e-12,61453.2,32.0389], Tmin=(100,'K'), Tmax=(950.7,'K')), NASAPolynomial(coeffs=[23.0607,0.0192025,-6.03723e-06,9.89484e-10,-6.63816e-14,55929.8,-87.5], Tmin=(950.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C1=CC=[C]C=CC1(5766)',
    structure = SMILES('[O]C1=CC=[C]C=CC1'),
    E0 = (287.745,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269515,0.0682813,-3.44407e-05,-1.92435e-08,1.69793e-11,34754.8,21.5267], Tmin=(100,'K'), Tmax=(918.739,'K')), NASAPolynomial(coeffs=[21.4963,0.0121483,-2.03359e-06,2.27115e-10,-1.59173e-14,29323.1,-87.4149], Tmin=(918.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]O[O](1885)',
    structure = SMILES('[O]O[O]'),
    E0 = (192.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([767.866,3898.78,3898.78],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4747,0.0046895,-6.33891e-06,3.99272e-09,-7.67843e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.29,'K')), NASAPolynomial(coeffs=[2.72453,0.000671697,1.37806e-06,-3.54977e-10,2.60903e-14,24449.8,18.0441], Tmin=(1873.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[C]1=CC=[C]CC=C1(6891)',
    structure = SMILES('[C]1=CC=[C]CC=C1'),
    E0 = (602.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.1225,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.86783,0.0581809,-2.94151e-05,-1.13792e-08,1.09896e-11,72567.4,18.2935], Tmin=(100,'K'), Tmax=(943.72,'K')), NASAPolynomial(coeffs=[16.9001,0.0165006,-4.92624e-06,8.21747e-10,-5.74805e-14,68371.4,-64.3231], Tmin=(943.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[O]OOC1=CC=C2[CH]C2C1(8001)',
    structure = SMILES('[O]OOC1=CC=C2[CH]C2C1'),
    E0 = (513.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627627,0.0546444,1.33072e-05,-6.3915e-08,3.02155e-11,61867.1,33.0477], Tmin=(100,'K'), Tmax=(961.991,'K')), NASAPolynomial(coeffs=[20.5403,0.0195257,-6.27795e-06,1.17895e-09,-9.01727e-14,55829.7,-73.7147], Tmin=(961.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_0_2) + radical(ROOJ) + radical(Allyl_S)"""),
)

species(
    label = '[O]OOC12[CH]C=C1C=CC2(8002)',
    structure = SMILES('[O]OOC12C=C[C]1C=CC2'),
    E0 = (482.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482022,0.0528415,3.33849e-05,-9.49511e-08,4.3899e-11,58178.9,33.1351], Tmin=(100,'K'), Tmax=(940.754,'K')), NASAPolynomial(coeffs=[24.318,0.0134123,-2.47444e-06,4.24346e-10,-3.89824e-14,50954.2,-94.9665], Tmin=(940.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(CCJ(C)CO) + radical(ROOJ)"""),
)

species(
    label = '[C]1=CC2OOO[C]2CC=C1(8003)',
    structure = SMILES('[C]1=CC2OOO[C]2CC=C1'),
    E0 = (445.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0869447,0.122626,-0.000260197,2.96407e-07,-1.21579e-10,53747.8,15.4272], Tmin=(100,'K'), Tmax=(837.373,'K')), NASAPolynomial(coeffs=[-10.0585,0.0873613,-4.85307e-05,9.66823e-09,-6.76565e-13,58324.1,79.1185], Tmin=(837.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_5_7) - ring(Cycloheptane) - ring(123trioxolane) + ring(1,3-Cycloheptadiene) + ring(123trioxolane) + radical(C=CJC=C) + radical(C2CsJOO)"""),
)

species(
    label = '[O]OO[C]1CC=CC2=CC12(8004)',
    structure = SMILES('[O]OO[C]1CC=CC2=CC12'),
    E0 = (621.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141438,0.0810931,-7.99292e-05,4.19317e-08,-8.81836e-12,74898.1,29.1369], Tmin=(100,'K'), Tmax=(1151.08,'K')), NASAPolynomial(coeffs=[15.3923,0.0280954,-1.08653e-05,1.93141e-09,-1.30627e-13,71387.2,-46.585], Tmin=(1151.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(s2_3_6_ane) - ring(Cyclohexane) - ring(Cyclopropane) + ring(Cyclohexene) + ring(Cyclopropene) + radical(C2CsJOO) + radical(ROOJ)"""),
)

species(
    label = '[C]1[CH]C2CC(=CC=1)OOO2(8005)',
    structure = SMILES('[C]1[CH]C=C2CC(C=1)OOO2'),
    E0 = (472.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23711,0.00616441,0.000150492,-1.96348e-07,7.36784e-11,56891,22.628], Tmin=(100,'K'), Tmax=(970.403,'K')), NASAPolynomial(coeffs=[17.3954,0.0267793,-9.82063e-06,2.02994e-09,-1.62516e-13,50036.5,-70.206], Tmin=(970.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(472.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_6_7_ane) - ring(Cycloheptane) - ring(123trioxane) + ring(1,4-Cycloheptadiene) + ring(123trioxane) + radical(Cds_S) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH]=C(CC=CC#C)OO[O](8006)',
    structure = SMILES('[CH]=C(CC=CC#C)OO[O]'),
    E0 = (662.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,750,770,3400,2100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.424521,0.0902147,-9.30455e-05,4.80036e-08,-9.69721e-12,79824.2,39.6127], Tmin=(100,'K'), Tmax=(1209.85,'K')), NASAPolynomial(coeffs=[20.1141,0.0223102,-8.85633e-06,1.61284e-09,-1.11194e-13,74854.5,-63.3871], Tmin=(1209.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'OOOC1=CC=C=C=CC1(8007)',
    structure = SMILES('OOOC1=CC=C=C=CC1'),
    E0 = (432.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.205966,0.0673874,-2.16036e-05,-2.649e-08,1.62828e-11,52153.3,30.2691], Tmin=(100,'K'), Tmax=(991.588,'K')), NASAPolynomial(coeffs=[20.3778,0.0223498,-8.43776e-06,1.6117e-09,-1.18955e-13,46366.6,-75.8852], Tmin=(991.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + ring(1_2_3_cycloheptatriene)"""),
)

species(
    label = 'O(4)(4)',
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
    label = 'C1=CC2=CC=C(C1)OO2(6904)',
    structure = SMILES('C1=CC2=CC=C(C1)OO2'),
    E0 = (107.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02738,0.0776636,-6.30776e-05,2.35743e-08,-3.63966e-12,13026.1,-8.6986], Tmin=(100,'K'), Tmax=(1429.23,'K')), NASAPolynomial(coeffs=[14.77,0.039202,-2.27116e-05,4.74556e-09,-3.46156e-13,9097.86,-79.9068], Tmin=(1429.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxane) + ring(Cycloheptane) + ring(1,3-Cyclohexadiene)"""),
)

species(
    label = 'C1=CC2=CC=C(C1)OOO2(6906)',
    structure = SMILES('C1=CC2=CC=C(C1)OOO2'),
    E0 = (129.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48805,0.153762,-0.000322969,3.3931e-07,-1.30499e-10,15784.9,-1.56941], Tmin=(100,'K'), Tmax=(848.671,'K')), NASAPolynomial(coeffs=[-1.3241,0.0737655,-4.15506e-05,8.24733e-09,-5.73097e-13,18610.1,14.4754], Tmin=(848.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cycloheptane) + ring(Cycloheptane) + ring(1,3,5-Cycloheptatriene)"""),
)

species(
    label = '[O]OC1C=[C]C=CCC1=O(8008)',
    structure = SMILES('[O]OC1C=[C]C=CCC1=O'),
    E0 = (222.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620223,0.0182043,0.000192553,-3.06172e-07,1.31127e-10,26955.2,31.1678], Tmin=(100,'K'), Tmax=(913.032,'K')), NASAPolynomial(coeffs=[43.2492,-0.0222354,1.86076e-05,-3.64301e-09,2.31268e-13,13072.2,-204.012], Tmin=(913.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(378.308,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cycloheptane) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]OC1=CC=[C]C=CC1(8009)',
    structure = SMILES('[O]OC1=CC=[C]C=CC1'),
    E0 = (426.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.35274,0.0866808,-9.01709e-05,4.74667e-08,-9.72796e-12,51453.3,27.0604], Tmin=(100,'K'), Tmax=(1200.33,'K')), NASAPolynomial(coeffs=[19.6903,0.0198878,-6.70131e-06,1.10672e-09,-7.21236e-14,46641.7,-73.2954], Tmin=(1200.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]OOC1C=CC=C=CC=1(6880)',
    structure = SMILES('[O]OOC1C=CC=C=CC=1'),
    E0 = (494.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (137.113,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.578537,0.0628271,-2.85863e-05,-1.16148e-08,9.75043e-12,59600.4,27.9299], Tmin=(100,'K'), Tmax=(1005.12,'K')), NASAPolynomial(coeffs=[17.3442,0.0228061,-8.70676e-06,1.62845e-09,-1.17105e-13,54881.4,-59.7493], Tmin=(1005.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1,2,4,6-Cycloheptatetraene) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1[CH]CC=C=CC=1(8010)',
    structure = SMILES('[O]OOC1C=C[C]=CCC=1'),
    E0 = (463.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.92897,0.096267,-0.000101672,5.35423e-08,-1.08831e-11,55883.8,32.4295], Tmin=(100,'K'), Tmax=(1237.93,'K')), NASAPolynomial(coeffs=[22.7023,0.0184569,-5.62888e-06,8.71836e-10,-5.48124e-14,50144.3,-86.1716], Tmin=(1237.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1=CC=C=[C]CC1(8011)',
    structure = SMILES('[O]OO[C]1C=CC#CCC1'),
    E0 = (594.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.552527,0.0633375,-2.57301e-05,-7.84396e-09,5.76462e-12,71614.4,30.3291], Tmin=(100,'K'), Tmax=(1163.61,'K')), NASAPolynomial(coeffs=[16.4898,0.0316566,-1.46748e-05,2.88647e-09,-2.07047e-13,66341.2,-55.6954], Tmin=(1163.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(594.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + ring(cycloheptyne) + radical(ROOJ) + radical(C2CsJOO)"""),
)

species(
    label = '[O]OOC1=[C]C=C=CCC1(8012)',
    structure = SMILES('[O]OOC1=[C]C=C=CCC1'),
    E0 = (592.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.344204,0.0692909,-3.86614e-05,-1.75551e-09,6.32289e-12,71422.3,29.5673], Tmin=(100,'K'), Tmax=(1007.66,'K')), NASAPolynomial(coeffs=[16.6069,0.0274754,-1.02669e-05,1.85523e-09,-1.29513e-13,66990.3,-54.7439], Tmin=(1007.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=CJC=C) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1=C[C]=C=CCC1(8013)',
    structure = SMILES('[O]OOC1=CC#C[CH]CC1'),
    E0 = (565.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546218,0.0631661,-2.42803e-05,-1.23441e-08,8.42877e-12,68106.2,30.1447], Tmin=(100,'K'), Tmax=(1070.6,'K')), NASAPolynomial(coeffs=[16.3531,0.0294266,-1.24815e-05,2.3981e-09,-1.71884e-13,63270.7,-53.969], Tmin=(1070.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + ring(cycloheptyne) + radical(ROOJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[O]OO[C]1CC2C=C=CC12(8014)',
    structure = SMILES('[O]OO[C]1CC2C=C=CC12'),
    E0 = (824.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.530737,0.0660052,-3.88257e-05,7.96742e-09,1.60363e-13,99277.7,27.4275], Tmin=(100,'K'), Tmax=(1329.61,'K')), NASAPolynomial(coeffs=[16.307,0.031213,-1.38677e-05,2.61987e-09,-1.81636e-13,93962.6,-57.3893], Tmin=(1329.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s2_4_5_ane) - ring(Cyclopentane) - ring(Cyclobutane) + ring(1,2-Cyclopentadiene) + ring(Cyclobutane) + radical(C2CsJOO) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC1=CC2[C]=CC2C1(8015)',
    structure = SMILES('[O]OOC1=CC2[C]=CC2C1'),
    E0 = (610.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435284,0.0608292,-5.22661e-06,-4.52494e-08,2.38511e-11,73595,34.4306], Tmin=(100,'K'), Tmax=(962.511,'K')), NASAPolynomial(coeffs=[20.7142,0.0189666,-6.08326e-06,1.12442e-09,-8.47552e-14,67726.7,-72.8342], Tmin=(962.511,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_4_5_diene_1_5) + radical(cyclobutene-vinyl) + radical(ROOJ)"""),
)

species(
    label = '[O]OOC12[CH]C=C=CC1C2(8016)',
    structure = SMILES('[O]OOC12C=C[C]=CC1C2'),
    E0 = (464.664,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41717,0.0624133,-8.00233e-06,-4.15761e-08,2.24773e-11,56030.1,29.9396], Tmin=(100,'K'), Tmax=(956.215,'K')), NASAPolynomial(coeffs=[19.5758,0.0220372,-7.04749e-06,1.25098e-09,-9.08459e-14,50548,-71.1383], Tmin=(956.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]1C=C2[CH]C=C(C1)OOO2(8017)',
    structure = SMILES('[CH]1C=C2CC=C[C]1OOO2'),
    E0 = (299.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0878488,0.0682674,-3.68095e-06,-5.08968e-08,2.51101e-11,36196.6,5.11352], Tmin=(100,'K'), Tmax=(1007.45,'K')), NASAPolynomial(coeffs=[23.893,0.0226897,-9.72465e-06,2.00838e-09,-1.54491e-13,28845.8,-123.26], Tmin=(1007.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(Cycloheptane) + ring(Cycloheptane) + ring(1,4-Cycloheptadiene) + radical(C2CsJOO) + radical(C=CCJC(O)C=C)"""),
)

species(
    label = '[C]1=C[CH]CC2=CC1OOO2(8018)',
    structure = SMILES('[C]1=C[CH]CC2=CC1OOO2'),
    E0 = (502.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0459,0.0236935,7.58579e-05,-9.81618e-08,3.25731e-11,60473.2,23.4832], Tmin=(100,'K'), Tmax=(1117.17,'K')), NASAPolynomial(coeffs=[12.4262,0.0417311,-2.24817e-05,4.75309e-09,-3.55086e-13,54709,-43.164], Tmin=(1117.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + Estimated bicyclic component: polycyclic(s3_6_7_ane) - ring(Cycloheptane) - ring(123trioxane) + ring(1,4-Cycloheptadiene) + ring(123trioxene) + radical(Allyl_S) + radical(Cds_S)"""),
)

species(
    label = 'C=CC=C=CC=[C]OO[O](8019)',
    structure = SMILES('C=CC=C=CC=[C]OO[O]'),
    E0 = (629.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,610,2055,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.531726,0.0905781,-8.82951e-05,3.76218e-08,-4.38464e-12,75906.7,39.659], Tmin=(100,'K'), Tmax=(970.87,'K')), NASAPolynomial(coeffs=[21.4469,0.0187163,-6.14509e-06,1.04079e-09,-7.09399e-14,70758.2,-70.262], Tmin=(970.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'OOOC1C=CC=C=CC=1(6901)',
    structure = SMILES('OOOC1C=CC=C=CC=1'),
    E0 = (342.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294302,0.0666423,-2.45729e-05,-2.01788e-08,1.32677e-11,41330.8,28.192], Tmin=(100,'K'), Tmax=(1007.4,'K')), NASAPolynomial(coeffs=[19.224,0.0239693,-9.41031e-06,1.80111e-09,-1.31578e-13,35868.2,-71.4552], Tmin=(1007.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(1,2,4,6-Cycloheptatetraene)"""),
)

species(
    label = '[O]OOC1C=CC#C[CH]C1(8020)',
    structure = SMILES('[O]OOC1C=CC#C[CH]C1'),
    E0 = (553.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.5058,0.0601814,-5.28006e-06,-3.772e-08,1.84932e-11,66759.7,29.5227], Tmin=(100,'K'), Tmax=(1018.49,'K')), NASAPolynomial(coeffs=[18.5735,0.0266696,-1.10759e-05,2.17365e-09,-1.60317e-13,61137.2,-67.5092], Tmin=(1018.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + ring(cycloheptyne) + radical(ROOJ) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C1=CC=C2CC(C=1)O2(5926)',
    structure = SMILES('C1=CC=C2CC(C=1)O2'),
    E0 = (109.032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.304266,0.0800348,-6.6118e-05,2.04597e-08,-1.24757e-12,13280.7,-10.0654], Tmin=(100,'K'), Tmax=(1203.84,'K')), NASAPolynomial(coeffs=[25.103,0.014797,-8.73306e-06,1.91757e-09,-1.458e-13,5773.41,-143.127], Tmin=(1203.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(oxepane) - ring(Oxetane) + ring(1_2_cycloheptadiene) + ring(Oxetane)"""),
)

species(
    label = 'C1=CC=C2CC(C=1)OO2(8021)',
    structure = SMILES('C1=CC=C2CC(C=1)OO2'),
    E0 = (212.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.173611,0.123557,-0.000263822,2.96597e-07,-1.20566e-10,25645,-9.46177], Tmin=(100,'K'), Tmax=(836.351,'K')), NASAPolynomial(coeffs=[-8.07075,0.0817289,-4.60455e-05,9.20965e-09,-6.45509e-13,29749.8,43.8689], Tmin=(836.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(PolycyclicRing) - ring(Cycloheptane) - ring(12dioxolane) + ring(1_2_cycloheptadiene) + ring(12dioxolane)"""),
)

species(
    label = '[CH2]C1C=C=CC=C1OO[O](8022)',
    structure = SMILES('[CH2]C1C=C=CC=C1OO[O]'),
    E0 = (657.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0728419,0.0774418,-6.26735e-05,2.05299e-08,-6.01222e-13,79280.4,33.1371], Tmin=(100,'K'), Tmax=(997.879,'K')), NASAPolynomial(coeffs=[17.4775,0.0246229,-8.7522e-06,1.52573e-09,-1.04075e-13,74963,-55.0213], Tmin=(997.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + ring(124cyclohexatriene) + radical(ROOJ) + radical(Isobutyl)"""),
)

species(
    label = 'C1=CC=C2CC(C=1)OOO2(8023)',
    structure = SMILES('C1=CC=C2CC(C=1)OOO2'),
    E0 = (338.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (138.121,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.31186,0.0329049,1.55492e-05,-2.09832e-08,4.61745e-12,40615.3,5.87991], Tmin=(100,'K'), Tmax=(1891.01,'K')), NASAPolynomial(coeffs=[19.7656,0.0417293,-2.43798e-05,4.70255e-09,-3.13081e-13,27348.3,-98.1466], Tmin=(1891.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cdd-CdsCds) + Estimated bicyclic component: polycyclic(s3_6_7_ane) - ring(Cycloheptane) - ring(123trioxane) + ring(1_2_cycloheptadiene) + ring(123trioxane)"""),
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
    label = 'Ar(8)',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (517.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (584.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (483.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (812.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (754.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (676.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (674.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (574.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (574.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (506.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (511.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (542.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (608.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (596.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (591.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (462.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (794.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (707.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (593.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (608.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (621.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (505.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (793.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (480.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (550.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (468.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (605.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (669.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (716.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (570.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (774.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (636.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (609.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (824.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (610.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (472.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (516.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (590.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (670.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (480.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (772.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (545.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (554.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (817.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (471.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['O2(2)(2)', 'C7H6O(490)(489)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.63066e+10,'s^-1'), n=0, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO;Y_rad_intra;OO_intra] for rate rule [R4OO_DSD;Cd_rad_in/OneDe;OO_intra]
Euclidian distance = 3.16227766017
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[C]1=C[CH]C2(CC=C1)OOO2(7990)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra_O] for rate rule [R5_SS_D;doublebond_intra_HDe_secNd;radadd_intra_O]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[C]1=CC=C2C[CH]C1OOO2(7991)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(414000,'s^-1'), n=1.14, Ea=(20.92,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7plus;doublebond_intra_HNd_pri;radadd_intra] for rate rule [R9;doublebond_intra_HNd_pri;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)(3)', '[O]OOC1=CC=C=C=CC1(7992)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)(3)', '[O]OOC1=CC#CC=CC1(7993)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.9e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.7612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2707 used for Ct-Cd_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cd_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OOC1=CC=C[C]=CC1(6903)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]OOC1=C[C]=CC=CC1(7994)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.0945e+07,'s^-1'), n=1.951, Ea=(212.263,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OOC1=CC=CC=[C]C1(7995)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.48e+06,'s^-1'), n=1.85, Ea=(112.257,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 128 used for R3H_SD;Cd_rad_out_Cd;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;Cd_rad_out_Cd;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OOC1=[C]C=CC=CC1(7996)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.48e+06,'s^-1'), n=1.85, Ea=(112.257,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SD;Cd_rad_out_Cd;Cd_H_out_single] for rate rule [R3H_SD;Cd_rad_out_Cd;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OOC1[CH]C=CC=CC=1(6847)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out;Cs_H_out_H/(Cd-Cd-Cd)] for rate rule [R4H_SDS;Cd_rad_out_Cd;Cs_H_out_H/(Cd-Cd-Cd)]
Euclidian distance = 2.82842712475
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['OOOC1C=C[CH][C]=CC=1(6890)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(12044.4,'s^-1'), n=1.9, Ea=(49.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5H_SSSS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OOOC1=[C]C=[C]C=CC1(7997)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 3.74165738677
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OOOC1=CC=[C]C=[C]C1(7998)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(564.492,'s^-1'), n=2.19647, Ea=(60.7245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;Y_rad_out;XH_out] for rate rule [R6H_SSSSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.16227766017
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OOOC1=C[C]=[C]C=CC1(7999)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5470,'s^-1'), n=1.94, Ea=(87.4456,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SMSSR;Y_rad_out;XH_out] for rate rule [R6H_SMSSR;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OOOC1=CC=[C][C]=CC1(8000)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;Cd_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O2(2)(2)', '[O]C1=CC=[C]C=CC1(5766)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.27681e+06,'m^3/(mol*s)'), n=0.151, Ea=(183.294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/OneDe;Y_rad] for rate rule [O_rad/OneDe;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.7 to 183.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]O[O](1885)', '[C]1=CC=[C]CC=C1(6891)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.93887e+07,'m^3/(mol*s)'), n=-0.0900543, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Y_rad;O_rad/NonDe] + [Cd_rad/NonDe;Y_rad] for rate rule [Cd_rad/NonDe;O_rad/NonDe]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OOC1=CC=C2[CH]C2C1(8001)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.00063e+13,'s^-1'), n=-0.283562, Ea=(245.321,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_beta;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_cddouble]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OOC12[CH]C=C1C=CC2(8002)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra_pri;radadd_intra] for rate rule [Rn0c7_gamma_short;doublebond_intra_pri_NdNd;radadd_intra_cdsingleDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[C]1=CC2OOO[C]2CC=C1(8003)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.16118e+10,'s^-1'), n=-0.279759, Ea=(145.715,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn3cx_alpha;doublebond_intra_secNd;radadd_intra] for rate rule [Rn3c7_alpha_short;doublebond_intra_secNd_HCd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OO[C]1CC=CC2=CC12(8004)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.46e+06,'s^-1'), n=1.02, Ea=(159.139,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;doublebond_intra_secNd;radadd_intra] for rate rule [Rn0c7_beta_long;doublebond_intra_secNd_HCd;radadd_intra_cddouble]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 157.7 to 159.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[C]1[CH]C2CC(=CC=1)OOO2(8005)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.95208e+10,'s^-1'), n=0.228823, Ea=(43.0818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn3c7_beta_long;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(CC=CC#C)OO[O](8006)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 830 used for R7_linear;triplebond_intra_H;radadd_intra_cdsingleH
Exact match found for rate rule [R7_linear;triplebond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['OOOC1=CC=C=C=CC1(8007)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.07e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad_De] for rate rule [R7radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['O(4)(4)', 'C1=CC2=CC=C(C1)OO2(6904)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.23936e+08,'s^-1'), n=0.739184, Ea=(87.8828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R5OO;Y_rad_intra;OO] for rate rule [R5OO;Cd_rad_in/OneDe;OOJ]
Euclidian distance = 3.16227766017
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['C1=CC2=CC=C(C1)OOO2(6906)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.35773e+13,'s^-1'), n=0.0154583, Ea=(6.01101,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Ypri_rad_out] for rate rule [R7;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OC1C=[C]C=CCC1=O(8008)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CHR;R2_doublebond_CsC;R_O_R]
Euclidian distance = 1.41421356237
family: ketoenol"""),
)

reaction(
    label = 'reaction28',
    reactants = ['O(4)(4)', '[O]OC1=CC=[C]C=CC1(8009)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(54738.4,'m^3/(mol*s)'), n=0.884925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)(3)', '[O]OOC1C=CC=C=CC=1(6880)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.61283,'m^3/(mol*s)'), n=2.04274, Ea=(10.5229,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 101 used for Cds-CdH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OOC1[CH]CC=C=CC=1(8010)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.03429e+08,'s^-1'), n=1.36289, Ea=(107.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/OneDe;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_H/(Cd-Cd-Cd);Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OOC1=CC=C=[C]CC1(8011)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.62853e+10,'s^-1'), n=0.944167, Ea=(180.133,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]OOC1=[C]C=C=CCC1(8012)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;Cs_H_out_H/Cd] for rate rule [R4H_DSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]OOC1=C[C]=C=CCC1(8013)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_single;Cs_H_out_1H] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OO[C]1CC2C=C=CC12(8014)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(361.916,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_gamma_short;doublebond_intra;radadd_intra_cs] for rate rule [Rn0c7_gamma_short;doublebond_intra_secNd_HCd;radadd_intra_csHCd]
Euclidian distance = 3.60555127546
family: Intra_R_Add_Endocyclic
Ea raised from 358.3 to 361.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OOC1=CC2[C]=CC2C1(8015)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(148.288,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Rn0c7_gamma_short_SSDS_D;doublebond_intra;radadd_intra_csHCs] for rate rule [Rn0c7_gamma_short_SSDS_D;doublebond_intra_CdCdd;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 144.3 to 148.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[O]OOC12[CH]C=C=CC1C2(8016)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.01748e+11,'s^-1'), n=0.19, Ea=(10.2508,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6_cyclic;doublebond_intra_pri;radadd_intra_csHCs] + [Rn0c7_beta_long_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_NdNd;radadd_intra_csHCs]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[CH]1C=C2[CH]C=C(C1)OOO2(8017)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.79849e+13,'s^-1'), n=-0.00708316, Ea=(54.5141,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn3c7_gamma_short;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['[C]1=C[CH]CC2=CC1OOO2(8018)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.7039e+11,'s^-1'), n=-0.18575, Ea=(127.616,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_CdCdd;radadd_intra] for rate rule [Rn3c7_beta_long;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 3.16227766017
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=CC=C=CC=[C]OO[O](8019)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.625e+11,'s^-1'), n=0.16, Ea=(41.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_pri;radadd_intra_cdsingle] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleNd]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['OOOC1C=CC=C=CC=1(6901)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.14e+09,'s^-1'), n=0.137, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_De] for rate rule [R5radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]OOC1C=CC#C[CH]C1(8020)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(4.16475e+09,'s^-1'), n=0.737748, Ea=(218.723,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_pentyn_3_ene;CH_end;unsaturated_end]
Euclidian distance = 0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['O2(2)(2)', 'C1=CC=C2CC(C=1)O2(5926)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_sec_rad_intra;OO_intra] for rate rule [R3OO_SS;C_rad/H/OneDe_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['O(4)(4)', 'C1=CC=C2CC(C=1)OO2(8021)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.0895e+10,'s^-1'), n=0.53, Ea=(92.048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;C_sec_rad_intra;OOJ] + [R4OO;C_sec_rad_intra;OO] for rate rule [R4OO;C_rad/H/OneDe_intra;OOJ]
Euclidian distance = 1.41421356237
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C1C=C=CC=C1OO[O](8022)'],
    products = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O]OOC1=CC=[C]C=CC1(6853)'],
    products = ['C1=CC=C2CC(C=1)OOO2(8023)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSSS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

network(
    label = '408',
    isomers = [
        '[O]OOC1=CC=[C]C=CC1(6853)',
    ],
    reactants = [
        ('O2(2)(2)', 'C7H6O(490)(489)'),
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '408',
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

