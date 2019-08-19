species(
    label = 'CC[CH]C(=O)COO[O](3482)',
    structure = SMILES('CCC=C([O])COO[O]'),
    E0 = (-22.9502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,180,180,180,180,1600,1853.31,2655.91,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160696,'amu*angstrom^2'), symmetry=1, barrier=(3.69471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160696,'amu*angstrom^2'), symmetry=1, barrier=(3.69471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160696,'amu*angstrom^2'), symmetry=1, barrier=(3.69471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160696,'amu*angstrom^2'), symmetry=1, barrier=(3.69471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160696,'amu*angstrom^2'), symmetry=1, barrier=(3.69471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0404142,0.0830742,-7.87346e-05,3.93171e-08,-7.8521e-12,-2609.93,37.3351], Tmin=(100,'K'), Tmax=(1210.97,'K')), NASAPolynomial(coeffs=[16.5675,0.0282147,-1.07801e-05,1.90584e-09,-1.28527e-13,-6632.19,-45.9675], Tmin=(1210.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.9502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=C(C)OJ)"""),
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
    label = 'CCC1OCC1=O(3198)',
    structure = SMILES('CCC1OCC1=O'),
    E0 = (-248.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4038.67,'J/mol'), sigma=(6.4294,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=630.83 K, Pc=34.48 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71322,0.0296461,5.9882e-05,-9.84158e-08,3.90783e-11,-29768.4,23.7596], Tmin=(100,'K'), Tmax=(979.251,'K')), NASAPolynomial(coeffs=[15.843,0.0231785,-8.71379e-06,1.72785e-09,-1.32224e-13,-34992.9,-56.6588], Tmin=(979.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-248.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutane)"""),
)

species(
    label = 'CC[CH]C1([O])COOO1(3483)',
    structure = SMILES('CC[CH]C1([O])COOO1'),
    E0 = (86.9127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0259177,0.0747505,-2.98692e-05,-2.77723e-08,2.08544e-11,10611.1,31.4217], Tmin=(100,'K'), Tmax=(901.283,'K')), NASAPolynomial(coeffs=[20.5453,0.0215733,-4.81023e-06,6.20287e-10,-3.82913e-14,5354.71,-74.2738], Tmin=(901.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.9127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(123trioxolane) + radical(CCJCOOH) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH2]OO[O](153)',
    structure = SMILES('[CH2]OO[O]'),
    E0 = (234.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,723.211,723.221,723.277],'cm^-1')),
        HinderedRotor(inertia=(0.00984315,'amu*angstrom^2'), symmetry=1, barrier=(3.65401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00984722,'amu*angstrom^2'), symmetry=1, barrier=(3.6543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.0248,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0178059,-4.92945e-06,-9.05903e-09,5.69358e-12,28256.2,14.7563], Tmin=(100,'K'), Tmax=(922.628,'K')), NASAPolynomial(coeffs=[8.13509,0.00563938,-1.46153e-06,2.21976e-10,-1.50714e-14,26884.3,-11.8494], Tmin=(922.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-OsHHH) + radical(ROOJ) + radical(CsJOO)"""),
)

species(
    label = 'CCC=C=O(2396)',
    structure = SMILES('CCC=C=O'),
    E0 = (-101.683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,311.417],'cm^-1')),
        HinderedRotor(inertia=(0.144733,'amu*angstrom^2'), symmetry=1, barrier=(9.96061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143318,'amu*angstrom^2'), symmetry=1, barrier=(9.96905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24367,0.0334186,-1.32661e-05,-1.87751e-09,1.89793e-12,-12162.2,17.1242], Tmin=(100,'K'), Tmax=(1185.57,'K')), NASAPolynomial(coeffs=[8.36996,0.0212602,-8.65138e-06,1.58277e-09,-1.08606e-13,-14213,-15.9971], Tmin=(1185.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.683,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH)"""),
)

species(
    label = 'CCC=C(O)[CH]OO[O](3484)',
    structure = SMILES('CCC=C(O)[CH]OO[O]'),
    E0 = (-43.4599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.25755,0.0766989,-3.38085e-05,-2.47994e-08,1.88478e-11,-5058.28,37.2101], Tmin=(100,'K'), Tmax=(948.058,'K')), NASAPolynomial(coeffs=[23.7558,0.0168234,-4.64088e-06,7.95861e-10,-5.94746e-14,-11473.9,-87.1816], Tmin=(948.058,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-43.4599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(ROOJ)"""),
)

species(
    label = 'CC[C]=C(O)COO[O](3485)',
    structure = SMILES('CC[C]=C(O)COO[O]'),
    E0 = (77.0867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,1685,370,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.441672,0.0916805,-9.61744e-05,5.18161e-08,-1.09772e-11,9436.46,37.7475], Tmin=(100,'K'), Tmax=(1154.78,'K')), NASAPolynomial(coeffs=[18.8859,0.0247323,-9.21207e-06,1.61173e-09,-1.0839e-13,4972.65,-58.2782], Tmin=(1154.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.0867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'CCC=C([O])[CH]OOO(3486)',
    structure = SMILES('CCC=C([O])[CH]OOO'),
    E0 = (-57.6598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196173,0.0767964,-3.82811e-05,-1.48828e-08,1.36725e-11,-6769.89,37.732], Tmin=(100,'K'), Tmax=(973.881,'K')), NASAPolynomial(coeffs=[21.9793,0.0212051,-7.31955e-06,1.34113e-09,-9.77893e-14,-12772.1,-77.3054], Tmin=(973.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-57.6598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CCJO)"""),
)

species(
    label = 'C[CH]C=C(O)COO[O](3487)',
    structure = SMILES('C[CH]C=C(O)COO[O]'),
    E0 = (-19.6424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.232734,0.0821321,-6.34903e-05,1.56172e-08,2.19442e-12,-2200.29,35.4166], Tmin=(100,'K'), Tmax=(986.533,'K')), NASAPolynomial(coeffs=[19.9045,0.0230985,-8.11654e-06,1.43395e-09,-9.97804e-14,-7273.99,-67.0384], Tmin=(986.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(ROOJ)"""),
)

species(
    label = '[CH2]CC=C(O)COO[O](3488)',
    structure = SMILES('[CH2]CC=C(O)COO[O]'),
    E0 = (44.4913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3615,1277.5,1000,3000,3100,440,815,1455,1000,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.434073,0.0885634,-8.45222e-05,3.84012e-08,-6.04891e-12,5518.52,39.0324], Tmin=(100,'K'), Tmax=(1028.98,'K')), NASAPolynomial(coeffs=[19.8316,0.022994,-8.19548e-06,1.42662e-09,-9.69664e-14,648.586,-62.7149], Tmin=(1028.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.4913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'CC[C]=C([O])COOO(3489)',
    structure = SMILES('CC[C]=C([O])COOO'),
    E0 = (62.8869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,180,180,180,180,1070.99,1146.69,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.165158,'amu*angstrom^2'), symmetry=1, barrier=(3.79731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165158,'amu*angstrom^2'), symmetry=1, barrier=(3.79731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165158,'amu*angstrom^2'), symmetry=1, barrier=(3.79731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165158,'amu*angstrom^2'), symmetry=1, barrier=(3.79731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165158,'amu*angstrom^2'), symmetry=1, barrier=(3.79731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165158,'amu*angstrom^2'), symmetry=1, barrier=(3.79731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.13971,0.0889112,-9.04902e-05,4.83775e-08,-1.03726e-11,7714.46,37.4085], Tmin=(100,'K'), Tmax=(1127.35,'K')), NASAPolynomial(coeffs=[16.3548,0.0303864,-1.26199e-05,2.32839e-09,-1.60867e-13,3995.44,-44.1452], Tmin=(1127.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.8869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C[CH]C=C([O])COOO(3490)',
    structure = SMILES('C[CH]C=C([O])COOO'),
    E0 = (-33.8423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.395797,0.084704,-7.57597e-05,3.43308e-08,-6.13472e-12,-3901.98,36.7546], Tmin=(100,'K'), Tmax=(1358.23,'K')), NASAPolynomial(coeffs=[20.1251,0.0242702,-9.01835e-06,1.57202e-09,-1.05089e-13,-9476.43,-68.53], Tmin=(1358.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.8423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC=C([O])COOO(3491)',
    structure = SMILES('[CH2]CC=C([O])COOO'),
    E0 = (30.2914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,350,440,435,1725,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.269319,0.0874026,-8.43616e-05,4.18811e-08,-8.25526e-12,3802.39,39.1855], Tmin=(100,'K'), Tmax=(1229.62,'K')), NASAPolynomial(coeffs=[18.3015,0.0269907,-1.06654e-05,1.92474e-09,-1.31502e-13,-764.593,-54.2465], Tmin=(1229.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.2914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(RCCJ)"""),
)

species(
    label = 'CCC=C([O])C[O](3185)',
    structure = SMILES('CCC=C([O])C[O]'),
    E0 = (-56.7318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,306.647,306.815,310.866,3900.7],'cm^-1')),
        HinderedRotor(inertia=(0.00178854,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095552,'amu*angstrom^2'), symmetry=1, barrier=(6.46406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.095449,'amu*angstrom^2'), symmetry=1, barrier=(6.46141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5192,0.0582795,-4.76049e-05,2.39043e-08,-5.44871e-12,-6737.06,26.9742], Tmin=(100,'K'), Tmax=(991.402,'K')), NASAPolynomial(coeffs=[6.56758,0.03791,-1.67842e-05,3.17815e-09,-2.21998e-13,-7738.01,2.66249], Tmin=(991.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.7318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ)"""),
)

species(
    label = '[O]O[O](145)',
    structure = SMILES('[O]O[O]'),
    E0 = (192.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([767.866,3898.78,3898.78],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3633.48,'J/mol'), sigma=(5.77645,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.54 K, Pc=42.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4747,0.0046895,-6.33891e-06,3.99272e-09,-7.67843e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.29,'K')), NASAPolynomial(coeffs=[2.72453,0.000671697,1.37806e-06,-3.54977e-10,2.60903e-14,24449.8,18.0441], Tmin=(1873.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C([O])[CH]CC(2377)',
    structure = SMILES('[CH2]C([O])=CCC'),
    E0 = (29.6155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,331.592,331.647],'cm^-1')),
        HinderedRotor(inertia=(0.00153279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110666,'amu*angstrom^2'), symmetry=1, barrier=(8.63422,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11059,'amu*angstrom^2'), symmetry=1, barrier=(8.63483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2878,0.0517027,-2.50274e-05,-6.51759e-09,7.14957e-12,3666.72,23.8661], Tmin=(100,'K'), Tmax=(960.293,'K')), NASAPolynomial(coeffs=[12.8498,0.0219906,-7.43267e-06,1.27278e-09,-8.66509e-14,595.546,-35.8741], Tmin=(960.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.6155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC[CH][C]=O(2394)',
    structure = SMILES('CC[CH][C]=O'),
    E0 = (100.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,663.318],'cm^-1')),
        HinderedRotor(inertia=(0.0020755,'amu*angstrom^2'), symmetry=1, barrier=(15.6417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314446,'amu*angstrom^2'), symmetry=1, barrier=(7.22973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.63497,'amu*angstrom^2'), symmetry=1, barrier=(83.575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02045,0.0379613,-2.55997e-05,9.29985e-09,-1.39601e-12,12149.6,19.3286], Tmin=(100,'K'), Tmax=(1547.01,'K')), NASAPolynomial(coeffs=[9.61093,0.0183349,-6.56951e-06,1.09894e-09,-7.07096e-14,9801.14,-20.6029], Tmin=(1547.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC1O[C]1COO[O](3492)',
    structure = SMILES('CCC1O[C]1COO[O]'),
    E0 = (116.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0617082,0.0844423,-7.71335e-05,2.9142e-08,1.94723e-13,14111.4,33.8119], Tmin=(100,'K'), Tmax=(806.959,'K')), NASAPolynomial(coeffs=[15.0576,0.0288343,-8.57375e-06,1.25669e-09,-7.47391e-14,11081.5,-39.0958], Tmin=(806.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(ROOJ) + radical(C2CsJO)"""),
)

species(
    label = 'CCC1OOOC[C]1[O](3493)',
    structure = SMILES('CCC1OOOC[C]1[O]'),
    E0 = (78.6873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.941121,0.0925753,-9.05907e-05,4.61942e-08,-9.05278e-12,9655.55,29.6867], Tmin=(100,'K'), Tmax=(1399.07,'K')), NASAPolynomial(coeffs=[20.2125,0.0225031,-5.17803e-06,5.93465e-10,-2.86007e-14,4675.38,-76.1152], Tmin=(1399.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.6873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(123trioxane) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
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
    label = 'CC=C([O])COO[O](3494)',
    structure = SMILES('CC=C([O])COO[O]'),
    E0 = (-0.30475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (118.088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716296,0.0699335,-7.19259e-05,3.93702e-08,-8.62259e-12,83.5664,31.7072], Tmin=(100,'K'), Tmax=(1108.1,'K')), NASAPolynomial(coeffs=[13.5486,0.0236116,-9.22165e-06,1.64555e-09,-1.11504e-13,-2760.35,-31.5187], Tmin=(1108.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.30475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC=C1COO1(3194)',
    structure = SMILES('CCC=C1COO1'),
    E0 = (32.1945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2750,3150,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75612,0.0347556,3.23581e-05,-6.02328e-08,2.32253e-11,3965.58,24.9512], Tmin=(100,'K'), Tmax=(1026.74,'K')), NASAPolynomial(coeffs=[11.905,0.0304365,-1.2785e-05,2.48762e-09,-1.81144e-13,25.1527,-33.3192], Tmin=(1026.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.1945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutane)"""),
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
    label = 'CCC=C1COOO1(3495)',
    structure = SMILES('CCC=C1COOO1'),
    E0 = (-13.5951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74505,0.0271161,7.5635e-05,-1.16676e-07,4.60139e-11,-1533.67,29.2944], Tmin=(100,'K'), Tmax=(967.305,'K')), NASAPolynomial(coeffs=[15.7986,0.0252944,-8.83242e-06,1.70063e-09,-1.29474e-13,-6886.06,-51.6515], Tmin=(967.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.5951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentane)"""),
)

species(
    label = 'CCC=C1COOOO1(3496)',
    structure = SMILES('CCC=C1COOOO1'),
    E0 = (0.95187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16856,0.0330718,8.71539e-05,-1.37886e-07,5.40803e-11,242.456,31.0675], Tmin=(100,'K'), Tmax=(988.342,'K')), NASAPolynomial(coeffs=[21.2664,0.0244439,-1.01056e-05,2.15596e-09,-1.71912e-13,-7281.55,-83.623], Tmin=(988.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.95187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclohexanone)"""),
)

species(
    label = 'CC[CH]C(=O)CO[O](2927)',
    structure = SMILES('CCC=C([O])CO[O]'),
    E0 = (-58.9272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3010,987.5,1337.5,450,1655,492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,294.116,294.122,3815.49],'cm^-1')),
        HinderedRotor(inertia=(0.177638,'amu*angstrom^2'), symmetry=1, barrier=(10.9092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177695,'amu*angstrom^2'), symmetry=1, barrier=(10.9093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177674,'amu*angstrom^2'), symmetry=1, barrier=(10.9092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177684,'amu*angstrom^2'), symmetry=1, barrier=(10.9091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640894,0.0756591,-7.60292e-05,4.33127e-08,-1.02458e-11,-6967.81,30.7824], Tmin=(100,'K'), Tmax=(1009.5,'K')), NASAPolynomial(coeffs=[11.0892,0.0342594,-1.45147e-05,2.68918e-09,-1.85617e-13,-9077.35,-19.7235], Tmin=(1009.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-58.9272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCC=[C]COO[O](3497)',
    structure = SMILES('CCC=[C]COO[O]'),
    E0 = (291.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891437,0.0690029,-5.96404e-05,2.86507e-08,-5.79079e-12,35188.2,32.9079], Tmin=(100,'K'), Tmax=(1156.62,'K')), NASAPolynomial(coeffs=[10.9151,0.034337,-1.46825e-05,2.73705e-09,-1.8959e-13,32869.5,-16.9091], Tmin=(1156.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(ROOJ)"""),
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
    label = 'CC=CC(=O)COO[O](3498)',
    structure = SMILES('CC=CC(=O)COO[O]'),
    E0 = (-48.2863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (131.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808991,0.0754731,-7.14835e-05,3.78069e-08,-8.54791e-12,-5696.98,32.2862], Tmin=(100,'K'), Tmax=(1029.5,'K')), NASAPolynomial(coeffs=[10.1157,0.0393131,-1.87979e-05,3.68963e-09,-2.63002e-13,-7613.23,-12.8837], Tmin=(1029.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-48.2863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + radical(ROOJ)"""),
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
    label = 'C=CC(=O)COO[O](3499)',
    structure = SMILES('C=CC(=O)COO[O]'),
    E0 = (-12.2607,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32156,0.0616024,-5.73956e-05,2.83687e-08,-5.84138e-12,-1380.45,28.6308], Tmin=(100,'K'), Tmax=(1137.31,'K')), NASAPolynomial(coeffs=[10.7778,0.0283441,-1.35311e-05,2.65628e-09,-1.89351e-13,-3531.38,-18.2068], Tmin=(1137.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.2607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C[CH]CC(=O)COO[O](3500)',
    structure = SMILES('C[CH]CC(=O)COO[O]'),
    E0 = (32.1609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301001,0.0830515,-8.15141e-05,4.36284e-08,-9.63501e-12,3999.88,36.6099], Tmin=(100,'K'), Tmax=(1077.59,'K')), NASAPolynomial(coeffs=[13.0094,0.0358782,-1.58494e-05,3.00421e-09,-2.10271e-13,1260.97,-25.6507], Tmin=(1077.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.1609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(ROOJ) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2]CCC(=O)COO[O](3501)',
    structure = SMILES('[CH2]CCC(=O)COO[O]'),
    E0 = (37.5052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.198297,0.0882246,-9.73271e-05,6.04702e-08,-1.55879e-11,4643.82,36.2779], Tmin=(100,'K'), Tmax=(929.567,'K')), NASAPolynomial(coeffs=[11.6583,0.0389105,-1.77501e-05,3.3984e-09,-2.38649e-13,2513.28,-18.173], Tmin=(929.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.5052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(RCCJ) + radical(ROOJ)"""),
)

species(
    label = 'CCCC(=O)[CH]OO[O](3502)',
    structure = SMILES('CCCC(=O)[CH]OO[O]'),
    E0 = (-29.2564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.122024,0.0906367,-0.000103784,6.81879e-08,-1.86574e-11,-3383.73,34.5234], Tmin=(100,'K'), Tmax=(877.42,'K')), NASAPolynomial(coeffs=[11.0314,0.0409031,-1.87619e-05,3.5885e-09,-2.51432e-13,-5298.16,-16.6814], Tmin=(877.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.2564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(OCJC=O) + radical(ROOJ)"""),
)

species(
    label = 'CC[CH][C]1COOOO1(3503)',
    structure = SMILES('CC[CH][C]1COOOO1'),
    E0 = (249.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.820191,0.0508693,2.69233e-05,-6.5146e-08,2.58356e-11,30184.1,32.1548], Tmin=(100,'K'), Tmax=(1048.14,'K')), NASAPolynomial(coeffs=[16.4662,0.036113,-1.62915e-05,3.25914e-09,-2.39958e-13,24435,-55.843], Tmin=(1048.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclohexane) + radical(CCJCOOH) + radical(C2CsJOO)"""),
)

species(
    label = 'CC=CC(=O)COOO(3504)',
    structure = SMILES('CC=CC(=O)COOO'),
    E0 = (-200.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550978,0.0791643,-6.77443e-05,3.03059e-08,-5.67967e-12,-23968.2,32.4399], Tmin=(100,'K'), Tmax=(1230.17,'K')), NASAPolynomial(coeffs=[12.9681,0.0387891,-1.85131e-05,3.62602e-09,-2.57695e-13,-27023.3,-30.0379], Tmin=(1230.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H)"""),
)

species(
    label = 'CCC1OOCC1=O(3505)',
    structure = SMILES('CCC1OOCC1=O'),
    E0 = (-279.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37168,0.0342765,6.56613e-05,-1.10414e-07,4.42381e-11,-33491.9,26.8578], Tmin=(100,'K'), Tmax=(976.842,'K')), NASAPolynomial(coeffs=[17.9741,0.0246552,-9.18387e-06,1.82776e-09,-1.40833e-13,-39520,-67.1028], Tmin=(976.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane)"""),
)

species(
    label = '[CH2]C(C)C(=O)COO[O](3506)',
    structure = SMILES('[CH2]C(C)C(=O)COO[O]'),
    E0 = (37.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.182467,0.0962558,-0.000115396,7.56995e-08,-2.01521e-11,4598.1,36.0694], Tmin=(100,'K'), Tmax=(910.776,'K')), NASAPolynomial(coeffs=[13.5362,0.0360043,-1.61631e-05,3.06204e-09,-2.13439e-13,2099.22,-28.833], Tmin=(910.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + radical(ROOJ) + radical(CJC(C)C=O)"""),
)

species(
    label = 'CCC([C]=O)COO[O](3507)',
    structure = SMILES('CCC([C]=O)COO[O]'),
    E0 = (1.70365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263608,0.0845229,-8.67822e-05,4.93612e-08,-1.15855e-11,337.532,36.6587], Tmin=(100,'K'), Tmax=(1018.76,'K')), NASAPolynomial(coeffs=[12.5183,0.0364072,-1.59384e-05,3.00218e-09,-2.09255e-13,-2159.41,-22.691], Tmin=(1018.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.70365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CCC1OOOCC1=O(3508)',
    structure = SMILES('CCC1OOOCC1=O'),
    E0 = (-264.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797369,0.0401989,7.73363e-05,-1.31892e-07,5.24515e-11,-31715.9,28.6236], Tmin=(100,'K'), Tmax=(995.316,'K')), NASAPolynomial(coeffs=[23.4696,0.0237602,-1.04324e-05,2.27745e-09,-1.82815e-13,-39928,-99.2324], Tmin=(995.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclohexanone)"""),
)

species(
    label = 'C=C([CH]CC)OOO[O](3509)',
    structure = SMILES('[CH2]C(=CCC)OOO[O]'),
    E0 = (240.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.115,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0920014,0.0754549,-5.23803e-05,7.76204e-09,4.36302e-12,29047.2,39.688], Tmin=(100,'K'), Tmax=(977.082,'K')), NASAPolynomial(coeffs=[18.1809,0.0240788,-8.32162e-06,1.45369e-09,-1.0044e-13,24429.8,-52.7005], Tmin=(977.082,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(O)CJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]CC(144)',
    structure = SMILES('[CH]CC'),
    E0 = (328.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,328.86,330.153,2894.73],'cm^-1')),
        HinderedRotor(inertia=(0.00155435,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00154734,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00371,0.0161454,1.44268e-05,-2.62511e-08,1.01943e-11,39516.8,12.6846], Tmin=(100,'K'), Tmax=(1003.27,'K')), NASAPolynomial(coeffs=[6.63826,0.0156478,-5.75067e-06,1.05874e-09,-7.51289e-14,38083.3,-8.37163], Tmin=(1003.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]OOC[C]=O(3510)',
    structure = SMILES('[O]OOC[C]=O'),
    E0 = (90.9977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,180,180,180,662.235,666.021],'cm^-1')),
        HinderedRotor(inertia=(0.00181314,'amu*angstrom^2'), symmetry=1, barrier=(20.5865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00638396,'amu*angstrom^2'), symmetry=1, barrier=(2.00553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125726,'amu*angstrom^2'), symmetry=1, barrier=(39.5446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03199,0.0388353,-3.88225e-05,1.91823e-08,-3.69297e-12,11019.1,23.4073], Tmin=(100,'K'), Tmax=(1271.71,'K')), NASAPolynomial(coeffs=[11.5822,0.00879649,-3.39143e-06,6.08345e-10,-4.16118e-14,8590.07,-24.9624], Tmin=(1271.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.9977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(ROOJ)"""),
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
    E0 = (59.9479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (86.9127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (187.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (146.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (227.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (79.6336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (90.7708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (133.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (125.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (53.0731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (171.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-22.9502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (222.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (334.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (208.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (79.5451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (419.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (59.9479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (229.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (0.95187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (184.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (534.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (163.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (135.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (157.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (156.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (142.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (249.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (2.02302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (69.0978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (231.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (144.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (-14.0802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (554.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (419.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['O2(2)', 'CCC1OCC1=O(3198)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.04e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO;C_rad/H/NonDeC_intra;OO_intra]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CC[CH]C1([O])COOO1(3483)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.72e+06,'s^-1'), n=1.2, Ea=(109.863,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra_HNd;radadd_intra] for rate rule [R6_SSS_D;doublebond_intra_HNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 109.5 to 109.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]OO[O](153)', 'CCC=C=O(2396)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CsJ] for rate rule [Ck_O;CsJ-OsHH]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCC=C(O)[CH]OO[O](3484)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.65528e+06,'s^-1'), n=1.92506, Ea=(169.809,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC[C]=C(O)COO[O](3485)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCC=C([O])[CH]OOO(3486)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(657459,'s^-1'), n=1.87237, Ea=(102.584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['C[CH]C=C(O)COO[O](3487)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 289 used for R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]CC=C(O)COO[O](3488)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(484628,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;XH_out] for rate rule [R5H_SSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC[C]=C([O])COOO(3489)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(13085.7,'s^-1'), n=1.78967, Ea=(62.3074,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSSR;Y_rad_out;XH_out] for rate rule [R6H_DSSSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 4.58257569496
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['C[CH]C=C([O])COOO(3490)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.74e+06,'s^-1'), n=0.99, Ea=(76.0233,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 267 used for R7H_OOCCCC(Cs/Cs);O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R7H_OOCCCC(Cs/Cs);O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC=C([O])COOO(3491)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.27e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H;C_rad_out_2H;XH_out] for rate rule [R8H;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O2(2)', 'CCC=C([O])C[O](3185)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.75733e+07,'m^3/(mol*s)'), n=-0.180109, Ea=(42.4033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_rad/NonDe;Y_rad] for rate rule [O_rad/NonDe;O2_birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -0.2 to 42.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]O[O](145)', 'C=C([O])[CH]CC(2377)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.74095e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cs_rad;O_rad/NonDe] for rate rule [C_rad/H2/Cd;O_rad/NonDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]OO[O](153)', 'CC[CH][C]=O(2394)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.68901e+06,'m^3/(mol*s)'), n=0.124131, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;C_pri_rad] for rate rule [Y_rad;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Recombination
Ea raised from -1.2 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCC1O[C]1COO[O](3492)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_HNd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCC1OOOC[C]1[O](3493)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(684195,'s^-1'), n=1.09371, Ea=(102.495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS_D;doublebond_intra;radadd_intra] for rate rule [R6_SSS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(S)(23)', 'CC=C([O])COO[O](3494)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.87e+13,'cm^3/(mol*s)','*|/',0.25), n=-0.146, Ea=(0.0118826,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 6 used for carbene;C_pri/Cd
Exact match found for rate rule [carbene;C_pri/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['O2(2)', 'CCC=C1COO1(3194)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO_SS;Y_rad_intra;OO_intra]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['O(4)', 'CCC=C1COOO1(3495)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.14131e+09,'s^-1'), n=0.54, Ea=(252.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R4OO;Y_rad_intra;OO] for rate rule [R4OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation
Ea raised from 250.9 to 252.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCC=C1COOOO1(3496)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(23.9021,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination
Ea raised from 17.8 to 23.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['O(4)', 'CC[CH]C(=O)CO[O](2927)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', 'CCC=[C]COO[O](3497)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', 'CC=CC(=O)COO[O](3498)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.81571,'m^3/(mol*s)'), n=1.94461, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-OneDeH;HJ] for rate rule [Cds-CsH_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -2.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH3(17)', 'C=CC(=O)COO[O](3499)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.0263591,'m^3/(mol*s)'), n=2.28811, Ea=(11.3029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-OneDeH;CsJ-HHH] for rate rule [Cds-HH_Cds-COH;CsJ-HHH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[CH]CC(=O)COO[O](3500)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]CCC(=O)COO[O](3501)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(25000,'s^-1'), n=2.28, Ea=(119.244,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/OneDe] for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCCC(=O)[CH]OO[O](3502)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.03437e+07,'s^-1'), n=1.66522, Ea=(165.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out] for rate rule [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_OOH/H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CC[CH][C]1COOOO1(3503)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(463580,'s^-1'), n=1.14062, Ea=(272.827,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSS;multiplebond_intra;radadd_intra] for rate rule [R6_SSS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.73205080757
family: Intra_R_Add_Endocyclic
Ea raised from 269.1 to 272.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CC=CC(=O)COOO(3504)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['O(4)', 'CCC1OOCC1=O(3505)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.0895e+10,'s^-1'), n=0.53, Ea=(92.048,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;C_rad/H/NonDeC_intra;OOJ] + [R4OO;C_rad/H/NonDeC_intra;OO] for rate rule [R4OO;C_rad/H/NonDeC_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C)C(=O)COO[O](3506)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CCC([C]=O)COO[O](3507)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CC[CH]C(=O)COO[O](3482)'],
    products = ['CCC1OOOCC1=O(3508)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R6_SSSSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([CH]CC)OOO[O](3509)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]CC(144)', '[O]OOC[C]=O(3510)'],
    products = ['CC[CH]C(=O)COO[O](3482)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '1064',
    isomers = [
        'CC[CH]C(=O)COO[O](3482)',
    ],
    reactants = [
        ('O2(2)', 'CCC1OCC1=O(3198)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '1064',
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

