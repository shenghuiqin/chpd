species(
    label = '[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)',
    structure = SMILES('[CH2]C(O)(C=[C]C=C)C(=O)CC'),
    E0 = (79.4579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,370.055,781.76,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154794,'amu*angstrom^2'), symmetry=1, barrier=(3.55902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.47701,0.149744,-0.000209107,1.59927e-07,-4.88411e-11,9782.93,41.946], Tmin=(100,'K'), Tmax=(827.897,'K')), NASAPolynomial(coeffs=[17.5995,0.0494283,-2.13438e-05,3.89161e-09,-2.62174e-13,6572.33,-50.4332], Tmin=(827.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.4579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C=C(O)C(=O)CC(4626)',
    structure = SMILES('C=C(O)C(=O)CC'),
    E0 = (-358.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,332.629,332.633],'cm^-1')),
        HinderedRotor(inertia=(0.152178,'amu*angstrom^2'), symmetry=1, barrier=(11.9524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152152,'amu*angstrom^2'), symmetry=1, barrier=(11.9478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152157,'amu*angstrom^2'), symmetry=1, barrier=(11.9534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152209,'amu*angstrom^2'), symmetry=1, barrier=(11.954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4071.48,'J/mol'), sigma=(6.49965,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.96 K, Pc=33.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.785843,0.0702816,-6.54766e-05,3.23543e-08,-6.53762e-12,-42997.7,22.9463], Tmin=(100,'K'), Tmax=(1175.98,'K')), NASAPolynomial(coeffs=[12.9689,0.0288414,-1.26178e-05,2.38822e-09,-1.6711e-13,-45863,-37.8049], Tmin=(1175.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsHH)"""),
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
    label = '[CH2]C(O)=C([O])CC(4557)',
    structure = SMILES('[CH2]C(O)=C([O])CC'),
    E0 = (-184.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3615,1277.5,1000,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,402.631,402.631],'cm^-1')),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100034,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100035,'amu*angstrom^2'), symmetry=1, barrier=(11.5078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4303.36,'J/mol'), sigma=(7.05412,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.18 K, Pc=27.82 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.516209,0.0797578,-7.96062e-05,3.90546e-08,-7.19184e-12,-22064.1,30.4076], Tmin=(100,'K'), Tmax=(1528.49,'K')), NASAPolynomial(coeffs=[20.9225,0.0119894,-1.65417e-06,6.23996e-11,2.30461e-15,-27255.3,-77.6604], Tmin=(1528.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C[C]=CC1(O)CC1([O])CC(29783)',
    structure = SMILES('C=C[C]=CC1(O)CC1([O])CC'),
    E0 = (166.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8669,0.112881,-7.80936e-05,7.42939e-09,9.3575e-12,20205.1,40.6983], Tmin=(100,'K'), Tmax=(952.863,'K')), NASAPolynomial(coeffs=[26.5126,0.0317279,-1.01325e-05,1.71236e-09,-1.17797e-13,13072.5,-103.894], Tmin=(952.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=CJC=C) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1(O)C=C(C=C)C1([O])CC(29784)',
    structure = SMILES('[CH2]C1(O)C=C(C=C)C1([O])CC'),
    E0 = (170.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.68917,0.106934,-6.20846e-05,-7.45733e-09,1.38812e-11,20689.9,41.687], Tmin=(100,'K'), Tmax=(965.876,'K')), NASAPolynomial(coeffs=[26.4887,0.0320415,-1.06927e-05,1.87745e-09,-1.32321e-13,13296.8,-103.37], Tmin=(965.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CC(C)(O)CJ) + radical(C=CC(C)2OJ)"""),
)

species(
    label = '[CH2]C1[C]=CC(O)(C1)C(=O)CC(29785)',
    structure = SMILES('[CH2]C1[C]=CC(O)(C1)C(=O)CC'),
    E0 = (84.8227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23799,0.110523,-9.95958e-05,4.85505e-08,-9.65498e-12,10394.3,38.6636], Tmin=(100,'K'), Tmax=(1200.58,'K')), NASAPolynomial(coeffs=[18.3617,0.0452213,-1.80075e-05,3.24507e-09,-2.20798e-13,5688.17,-59.4761], Tmin=(1200.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.8227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(cyclopentene-vinyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C(O)(C#CC=C)C(=O)CC(29786)',
    structure = SMILES('[CH2]C(O)(C#CC=C)C(=O)CC'),
    E0 = (49.4809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2100,2250,500,550,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,246.18,893.193,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151916,'amu*angstrom^2'), symmetry=1, barrier=(3.49284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.182,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.05111,0.134458,-0.000161691,1.02775e-07,-2.59381e-11,6167.97,40.8964], Tmin=(100,'K'), Tmax=(969.095,'K')), NASAPolynomial(coeffs=[20.6512,0.0407589,-1.66695e-05,3.01745e-09,-2.0515e-13,1767.54,-67.9179], Tmin=(969.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.4809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCtH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C(O)(C=C=C=C)C(=O)CC(29787)',
    structure = SMILES('[CH2]C(O)(C=C=C=C)C(=O)CC'),
    E0 = (109.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1645.22,2868.82,3200],'cm^-1')),
        HinderedRotor(inertia=(0.152255,'amu*angstrom^2'), symmetry=1, barrier=(3.50065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152255,'amu*angstrom^2'), symmetry=1, barrier=(3.50065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152255,'amu*angstrom^2'), symmetry=1, barrier=(3.50065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152255,'amu*angstrom^2'), symmetry=1, barrier=(3.50065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152255,'amu*angstrom^2'), symmetry=1, barrier=(3.50065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152255,'amu*angstrom^2'), symmetry=1, barrier=(3.50065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.182,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.32257,0.149265,-0.000222097,1.81477e-07,-5.90558e-11,13429.6,40.876], Tmin=(100,'K'), Tmax=(823.397,'K')), NASAPolynomial(coeffs=[15.4725,0.051633,-2.3862e-05,4.47731e-09,-3.0589e-13,10878.3,-39.2144], Tmin=(823.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = 'C2H5CO(71)',
    structure = SMILES('CC[C]=O'),
    E0 = (-44.1874,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.192008,'amu*angstrom^2'), symmetry=1, barrier=(4.41465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191563,'amu*angstrom^2'), symmetry=1, barrier=(4.40441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3133.67,'J/mol'), sigma=(5.35118,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.47 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.91313,0.0225012,-8.59328e-06,4.80319e-10,2.48743e-13,-5274.24,14.655], Tmin=(100,'K'), Tmax=(1673.18,'K')), NASAPolynomial(coeffs=[7.46306,0.0158512,-6.42141e-06,1.12499e-09,-7.32055e-14,-7388.54,-11.406], Tmin=(1673.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-44.1874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2CO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C[C]=CC(=C)O(27798)',
    structure = SMILES('C=C[C]=CC(=C)O'),
    E0 = (130.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21632,'amu*angstrom^2'), symmetry=1, barrier=(27.9656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20741,'amu*angstrom^2'), symmetry=1, barrier=(27.7607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21858,'amu*angstrom^2'), symmetry=1, barrier=(28.0177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.642274,0.081296,-8.44264e-05,4.23515e-08,-7.89523e-12,15916.1,26.1388], Tmin=(100,'K'), Tmax=(1529.91,'K')), NASAPolynomial(coeffs=[21.645,0.00924721,-2.77545e-07,-2.03306e-10,2.05e-14,10709,-85.5914], Tmin=(1529.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
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
    label = 'OH(5)',
    structure = SMILES('[OH]'),
    E0 = (28.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3287.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0073,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4858,0.00133397,-4.70043e-06,5.64379e-09,-2.06318e-12,3411.96,1.99788], Tmin=(100,'K'), Tmax=(1005.25,'K')), NASAPolynomial(coeffs=[2.88225,0.00103869,-2.35652e-07,1.40229e-11,6.34581e-16,3669.56,5.59053], Tmin=(1005.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.372,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=C[C]=CC(=C)C(=O)CC(29788)',
    structure = SMILES('C=C[C]=CC(=C)C(=O)CC'),
    E0 = (160.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.183,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03904,0.107514,-0.000102175,5.15381e-08,-1.05245e-11,19447.4,33.8573], Tmin=(100,'K'), Tmax=(1174.45,'K')), NASAPolynomial(coeffs=[18.5782,0.0407014,-1.68424e-05,3.09975e-09,-2.13659e-13,14839.5,-63.9389], Tmin=(1174.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(O)([C]=CC=C)C(=O)CC(29789)',
    structure = SMILES('[CH2]C(O)([C]=CC=C)C(=O)CC'),
    E0 = (118.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,551.63,598.986,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158566,'amu*angstrom^2'), symmetry=1, barrier=(3.64575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.30448,0.146196,-0.000198717,1.47759e-07,-4.40894e-11,14448.8,42.2106], Tmin=(100,'K'), Tmax=(819.641,'K')), NASAPolynomial(coeffs=[17.4766,0.0496726,-2.20935e-05,4.11617e-09,-2.82019e-13,11205.7,-49.2896], Tmin=(819.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC(O)(C=O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C=C[C]=C)C(=O)CC(29790)',
    structure = SMILES('[CH2]C(O)([CH]C=C=C)C(=O)CC'),
    E0 = (61.9284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.46296,0.147029,-0.00018258,1.19892e-07,-3.14303e-11,7676.85,40.2103], Tmin=(100,'K'), Tmax=(930.749,'K')), NASAPolynomial(coeffs=[20.9227,0.0465271,-2.06116e-05,3.87892e-09,-2.69408e-13,3323.59,-70.9338], Tmin=(930.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.9284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJCO) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=C[C]=CC(C)([O])C(=O)CC(29791)',
    structure = SMILES('C=C[C]=CC(C)([O])C(=O)CC'),
    E0 = (110.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1050.74,1182.61,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162874,'amu*angstrom^2'), symmetry=1, barrier=(3.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162874,'amu*angstrom^2'), symmetry=1, barrier=(3.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162874,'amu*angstrom^2'), symmetry=1, barrier=(3.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162874,'amu*angstrom^2'), symmetry=1, barrier=(3.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162874,'amu*angstrom^2'), symmetry=1, barrier=(3.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162874,'amu*angstrom^2'), symmetry=1, barrier=(3.7448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.00016,0.142523,-0.000205648,1.70062e-07,-5.62173e-11,13498.6,41.1684], Tmin=(100,'K'), Tmax=(842.132,'K')), NASAPolynomial(coeffs=[12.3685,0.0579555,-2.59508e-05,4.79678e-09,-3.2478e-13,11657.1,-22.2478], Tmin=(842.132,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C[C]=[C]C(C)(O)C(=O)CC(29792)',
    structure = SMILES('[CH2][CH]C#CC(C)(O)C(=O)CC'),
    E0 = (68.2005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92466,0.134126,-0.000172813,1.25373e-07,-3.63191e-11,8412.3,44.1676], Tmin=(100,'K'), Tmax=(900.992,'K')), NASAPolynomial(coeffs=[16.0895,0.0477284,-1.82817e-05,3.11845e-09,-2.01468e-13,5426.89,-39.4153], Tmin=(900.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.2005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[CH]=CC=CC([CH2])(O)C(=O)CC(29793)',
    structure = SMILES('[CH]=CC=CC([CH2])(O)C(=O)CC'),
    E0 = (127.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3000,3100,440,815,1455,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1677.97,2832.06,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153773,'amu*angstrom^2'), symmetry=1, barrier=(3.53554,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.27944,0.144284,-0.000187354,1.31147e-07,-3.67448e-11,15562.4,42.054], Tmin=(100,'K'), Tmax=(872.193,'K')), NASAPolynomial(coeffs=[18.8288,0.0474774,-2.08636e-05,3.8881e-09,-2.67481e-13,11880.3,-56.8941], Tmin=(872.193,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([O])(C=CC=C)C(=O)CC(29794)',
    structure = SMILES('[CH2]C([O])(C=CC=C)C(=O)CC'),
    E0 = (124.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,180,1081.47,1148.63,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.165584,'amu*angstrom^2'), symmetry=1, barrier=(3.80709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165584,'amu*angstrom^2'), symmetry=1, barrier=(3.80709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165584,'amu*angstrom^2'), symmetry=1, barrier=(3.80709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165584,'amu*angstrom^2'), symmetry=1, barrier=(3.80709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165584,'amu*angstrom^2'), symmetry=1, barrier=(3.80709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165584,'amu*angstrom^2'), symmetry=1, barrier=(3.80709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.1585,0.145766,-0.000210989,1.73087e-07,-5.70377e-11,15215.7,42.1943], Tmin=(100,'K'), Tmax=(820.278,'K')), NASAPolynomial(coeffs=[13.6708,0.0565761,-2.59503e-05,4.86624e-09,-3.32847e-13,13022.6,-28.5761], Tmin=(820.278,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=C[C]=CC(C)(O)C(=O)[CH]C(29795)',
    structure = SMILES('C=C[C]=CC(C)(O)C([O])=CC'),
    E0 = (-2.36471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89781,0.120679,-0.000116543,5.87463e-08,-1.17787e-11,-64.4317,42.592], Tmin=(100,'K'), Tmax=(1210.12,'K')), NASAPolynomial(coeffs=[22.9912,0.0384103,-1.45692e-05,2.56885e-09,-1.73132e-13,-6088.27,-82.2304], Tmin=(1210.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-2.36471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(O)(C=CC=C)C(=O)[CH]C(29796)',
    structure = SMILES('[CH2]C(O)(C=CC=C)C([O])=CC'),
    E0 = (12.0833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99334,0.12065,-0.000114667,5.57808e-08,-1.07366e-11,1678.44,44.2163], Tmin=(100,'K'), Tmax=(1261.08,'K')), NASAPolynomial(coeffs=[24.7367,0.0358655,-1.38193e-05,2.46807e-09,-1.67751e-13,-5063.3,-90.941], Tmin=(1261.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.0833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]CC(=O)C(C)(O)C=[C]C=C(29797)',
    structure = SMILES('[CH2]CC(=O)C(C)(O)C=[C]C=C'),
    E0 = (77.8236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.27261,0.143638,-0.000188724,1.35204e-07,-3.87485e-11,9580.7,41.798], Tmin=(100,'K'), Tmax=(854.122,'K')), NASAPolynomial(coeffs=[18.2065,0.0477192,-2.02506e-05,3.6901e-09,-2.49833e-13,6082.79,-53.7704], Tmin=(854.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.8236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CJCC=O)"""),
)

species(
    label = 'C=[C][C]=CC(C)(O)C(=O)CC(29798)',
    structure = SMILES('C=[C][C]=CC(C)(O)C(=O)CC'),
    E0 = (65.2299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.33718,0.146745,-0.000204754,1.58421e-07,-4.87988e-11,8066.54,40.9849], Tmin=(100,'K'), Tmax=(860.636,'K')), NASAPolynomial(coeffs=[16.319,0.0507679,-2.13202e-05,3.81624e-09,-2.53604e-13,5198.58,-44.226], Tmin=(860.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.2299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C=CC=C(29799)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C=CC=C'),
    E0 = (92.0515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,511.167,637.401,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157676,'amu*angstrom^2'), symmetry=1, barrier=(3.62529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.37177,0.146124,-0.000191106,1.33854e-07,-3.7412e-11,11295.4,42.6153], Tmin=(100,'K'), Tmax=(875.253,'K')), NASAPolynomial(coeffs=[19.434,0.0464722,-2.03286e-05,3.77846e-09,-2.5949e-13,7478.13,-59.6803], Tmin=(875.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.0515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(CJCC=O)"""),
)

species(
    label = '[CH]=C[C]=CC(C)(O)C(=O)CC(29800)',
    structure = SMILES('[CH]C=C=CC(C)(O)C(=O)CC'),
    E0 = (90.7046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.81822,0.135856,-0.000163221,1.15083e-07,-3.37813e-11,11111.6,41.0853], Tmin=(100,'K'), Tmax=(821.681,'K')), NASAPolynomial(coeffs=[13.4139,0.0617039,-2.7852e-05,5.25035e-09,-3.6366e-13,8608.48,-29.4091], Tmin=(821.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.7046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(O)(C=C1[CH]C1)C(=O)CC(29801)',
    structure = SMILES('[CH2]C(O)(C=C1[CH]C1)C(=O)CC'),
    E0 = (114.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.49356,0.124331,-0.000132817,7.77628e-08,-1.86331e-11,13977.9,37.9035], Tmin=(100,'K'), Tmax=(1003.07,'K')), NASAPolynomial(coeffs=[17.0761,0.0502807,-2.2083e-05,4.167e-09,-2.90746e-13,10252.6,-51.7416], Tmin=(1003.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + ring(Methylene_cyclopropane) + radical(C=CC(O)(C=O)CJ) + radical(Allyl_S)"""),
)

species(
    label = 'C=C[C]=CC1(O)CO[C]1CC(29802)',
    structure = SMILES('C=C[C]=CC1(O)CO[C]1CC'),
    E0 = (157.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.09071,0.121145,-0.000122929,6.78661e-08,-1.46407e-11,19151.9,42.1586], Tmin=(100,'K'), Tmax=(1241.16,'K')), NASAPolynomial(coeffs=[21.4994,0.036263,-9.64349e-06,1.26807e-09,-6.8302e-14,13978.2,-73.9986], Tmin=(1241.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.323,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Oxetane) + radical(C2CsJOCs) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1(O)C=C(C=C)O[C]1CC(29803)',
    structure = SMILES('[CH2]C1(O)C=C(C=C)O[C]1CC'),
    E0 = (49.9768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.66427,0.121541,-0.000115471,5.56926e-08,-1.02929e-11,6272.11,41.2149], Tmin=(100,'K'), Tmax=(1473.01,'K')), NASAPolynomial(coeffs=[28.0111,0.0263095,-6.34373e-06,8.03601e-10,-4.37793e-14,-1470.45,-114.263], Tmin=(1473.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.9768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(C2CsJOC(O)) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'CCC(=O)C1(O)C=[C][CH]CC1(29804)',
    structure = SMILES('CCC(=O)C1(O)[CH][C]=CCC1'),
    E0 = (-33.3868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69209,0.111377,-9.0342e-05,3.73339e-08,-6.18185e-12,-3799.2,35.6322], Tmin=(100,'K'), Tmax=(1439.14,'K')), NASAPolynomial(coeffs=[23.9952,0.0399804,-1.59267e-05,2.86188e-09,-1.93551e-13,-11192.7,-97.6455], Tmin=(1439.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.3868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(Cds_S) + radical(C=CCJCO)"""),
)

species(
    label = 'CH3CH2OH(54)',
    structure = SMILES('CCO'),
    E0 = (-248.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00248522,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168417,'amu*angstrom^2'), symmetry=1, barrier=(8.09857,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0684,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3014.83,'J/mol'), sigma=(4.53,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.8587,-0.00374017,6.95554e-05,-8.86548e-08,3.51688e-11,-29996.1,4.80185], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.56244,0.0152042,-5.38968e-06,8.6225e-10,-5.12898e-14,-31525.6,-9.47302], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-248.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""CH3CH2OH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = '[CH2]C(=C=O)C=[C]C=C(29805)',
    structure = SMILES('C=C[C]=CC(=C)[C]=O'),
    E0 = (394.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21824,'amu*angstrom^2'), symmetry=1, barrier=(28.0097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21782,'amu*angstrom^2'), symmetry=1, barrier=(28.0002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21774,'amu*angstrom^2'), symmetry=1, barrier=(27.9982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.237791,0.0827695,-8.64105e-05,4.41231e-08,-8.69492e-12,47589.4,27.241], Tmin=(100,'K'), Tmax=(1251.16,'K')), NASAPolynomial(coeffs=[20.9168,0.015137,-5.32615e-06,9.17891e-10,-6.1819e-14,42295.9,-79.5577], Tmin=(1251.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=C(C)CJ=O)"""),
)

species(
    label = 'C=C=C=CC(C)(O)C(=O)CC(29806)',
    structure = SMILES('C=C=C=CC(C)(O)C(=O)CC'),
    E0 = (-103.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.82832,0.136199,-0.000175033,1.27477e-07,-3.79495e-11,-12230.6,39.0637], Tmin=(100,'K'), Tmax=(815.898,'K')), NASAPolynomial(coeffs=[14.8456,0.0544548,-2.47513e-05,4.68437e-09,-3.24872e-13,-14951.5,-37.9858], Tmin=(815.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-103.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
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
    label = '[CH2]C(O)(C=[C]C=C)C(C)=O(29807)',
    structure = SMILES('[CH2]C(O)(C=[C]C=C)C(C)=O'),
    E0 = (104.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,180,180,180,180,1575.27,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.1495,'amu*angstrom^2'), symmetry=1, barrier=(3.43731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1495,'amu*angstrom^2'), symmetry=1, barrier=(3.43731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1495,'amu*angstrom^2'), symmetry=1, barrier=(3.43731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1495,'amu*angstrom^2'), symmetry=1, barrier=(3.43731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1495,'amu*angstrom^2'), symmetry=1, barrier=(3.43731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1495,'amu*angstrom^2'), symmetry=1, barrier=(3.43731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55379,0.125786,-0.000163302,1.1327e-07,-3.11959e-11,12810.2,37.3056], Tmin=(100,'K'), Tmax=(890.299,'K')), NASAPolynomial(coeffs=[17.9265,0.0382658,-1.58487e-05,2.85804e-09,-1.92626e-13,9341.49,-54.4127], Tmin=(890.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC(O)(C=O)CJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=C[C](O)CC(=O)CC(29808)',
    structure = SMILES('[CH2]C=[C]C=C(O)CC(=O)CC'),
    E0 = (-55.1875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7686,0.116996,-0.000103075,4.64258e-08,-8.37405e-12,-6421.73,40.9513], Tmin=(100,'K'), Tmax=(1327.69,'K')), NASAPolynomial(coeffs=[23.5939,0.0405846,-1.67458e-05,3.0775e-09,-2.11632e-13,-13156.4,-88.5966], Tmin=(1327.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.1875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=C[C]=CC[C](O)C(=O)CC(29691)',
    structure = SMILES('C=C[C]=CCC(O)=C([O])CC'),
    E0 = (-10.4663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.80906,0.128661,-0.000128478,6.4508e-08,-1.24665e-11,-995.693,45.8349], Tmin=(100,'K'), Tmax=(1361.15,'K')), NASAPolynomial(coeffs=[29.5573,0.0260379,-7.11321e-06,1.01321e-09,-6.02576e-14,-9111.27,-117.739], Tmin=(1361.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.4663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=CC1=CC(O)(C1)C(=O)CC(29809)',
    structure = SMILES('C=CC1=CC(O)(C1)C(=O)CC'),
    E0 = (-173.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9416,0.117939,-0.000108408,5.12244e-08,-9.59266e-12,-20623,37.6979], Tmin=(100,'K'), Tmax=(1294.3,'K')), NASAPolynomial(coeffs=[24.5265,0.0361409,-1.36104e-05,2.39634e-09,-1.61353e-13,-27474.6,-96.8233], Tmin=(1294.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-173.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH2]C(C=[C]C=C)=C(CC)OO(29810)',
    structure = SMILES('[CH2]C(C=[C]C=C)=C(CC)OO'),
    E0 = (309.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.16653,0.127449,-0.00012776,6.63247e-08,-1.36508e-11,37496.8,42.0616], Tmin=(100,'K'), Tmax=(1181.59,'K')), NASAPolynomial(coeffs=[24.0926,0.0385561,-1.49149e-05,2.65682e-09,-1.8021e-13,31291.2,-89.0053], Tmin=(1181.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C[C]=CO[C](CC)C(=C)O(29687)',
    structure = SMILES('[CH2]C(O)=C(CC)OC=[C]C=C'),
    E0 = (64.4767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.99563,0.134082,-0.000138956,7.21034e-08,-1.43781e-11,8023.48,45.8519], Tmin=(100,'K'), Tmax=(1312.91,'K')), NASAPolynomial(coeffs=[30.3406,0.0254975,-6.8787e-06,9.65328e-10,-5.67758e-14,-125.027,-121.748], Tmin=(1312.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.4767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(O)(C=[C]C=C)C(=C)OC(29811)',
    structure = SMILES('[CH2]C(O)(C=[C]C=C)C(=C)OC'),
    E0 = (127.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,180,846.878,1359.77,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156362,'amu*angstrom^2'), symmetry=1, barrier=(3.59508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.36027,0.127958,-0.000127649,6.47323e-08,-1.28931e-11,15632,43.4566], Tmin=(100,'K'), Tmax=(1226.3,'K')), NASAPolynomial(coeffs=[26.457,0.0339605,-1.26717e-05,2.22548e-09,-1.50116e-13,8564.28,-101.449], Tmin=(1226.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]C(O)(C=[C]C=C)C(O)=CC(29812)',
    structure = SMILES('[CH2]C(O)(C=[C]C=C)C(O)=CC'),
    E0 = (73.2739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180,206.63,1532.11,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148893,'amu*angstrom^2'), symmetry=1, barrier=(3.42333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.43398,0.131178,-0.000136561,7.24302e-08,-1.50561e-11,9053.28,43.8903], Tmin=(100,'K'), Tmax=(1178.34,'K')), NASAPolynomial(coeffs=[26.3143,0.0335883,-1.23308e-05,2.14421e-09,-1.43868e-13,2278.27,-99.5213], Tmin=(1178.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.2739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC(C)(O)CJ)"""),
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
    label = 'C=C[C]=C[C](O)C(=O)CC(29813)',
    structure = SMILES('[CH2]C=[C]C=C(O)C(=O)CC'),
    E0 = (-25.6701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0292,0.108847,-0.000112703,6.24267e-08,-1.39171e-11,-2904.73,35.0585], Tmin=(100,'K'), Tmax=(1085.31,'K')), NASAPolynomial(coeffs=[18.0082,0.0386829,-1.57285e-05,2.85872e-09,-1.95557e-13,-7037,-58.3441], Tmin=(1085.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.6701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=[C]C1CC1(O)C(=O)CC(29814)',
    structure = SMILES('[CH2]C=[C]C1CC1(O)C(=O)CC'),
    E0 = (100.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69432,0.116107,-0.000106475,5.07641e-08,-9.68022e-12,12331.4,40.936], Tmin=(100,'K'), Tmax=(1263.61,'K')), NASAPolynomial(coeffs=[22.4983,0.0395235,-1.55645e-05,2.80014e-09,-1.90666e-13,6217.49,-81.4394], Tmin=(1263.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(100.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1(O)C=C=CCC1([O])CC(29815)',
    structure = SMILES('[CH2]C1(O)C=C=CCC1([O])CC'),
    E0 = (151.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-3.58418,0.175248,-0.000240739,1.75774e-07,-5.11097e-11,18493.4,28.9183], Tmin=(100,'K'), Tmax=(842.652,'K')), NASAPolynomial(coeffs=[22.0309,0.0536595,-2.43084e-05,4.55118e-09,-3.12678e-13,14176.3,-90.2753], Tmin=(842.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(C=CC(C)(O)CJ) + radical(CC(C)2OJ)"""),
)

species(
    label = '[CH2]C(O)(C=C=[C]C)C(=O)CC(29816)',
    structure = SMILES('[CH2]C(O)([CH]C#CC)C(=O)CC'),
    E0 = (128.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,3000,3100,440,815,1455,1000,2100,2250,500,550,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,182.517,963.202,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151371,'amu*angstrom^2'), symmetry=1, barrier=(3.48032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18786,0.138076,-0.000173855,1.18866e-07,-3.22731e-11,15646.8,44.9103], Tmin=(100,'K'), Tmax=(903.933,'K')), NASAPolynomial(coeffs=[19.238,0.0432679,-1.65335e-05,2.84292e-09,-1.85771e-13,11773.1,-56.294], Tmin=(903.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(128.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CJC(C)(C=O)O) + radical(CCJCO)"""),
)

species(
    label = '[CH2]C(O)([C]=C=CC)C(=O)CC(29817)',
    structure = SMILES('[CH2]C(O)(C#C[CH]C)C(=O)CC'),
    E0 = (74.5655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.97481,0.138466,-0.000166336,9.4943e-08,-1.48266e-11,9177.15,41.5369], Tmin=(100,'K'), Tmax=(696.447,'K')), NASAPolynomial(coeffs=[18.8319,0.0438907,-1.63255e-05,2.73699e-09,-1.75117e-13,5674.49,-55.6562], Tmin=(696.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.5655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(523.812,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Sec_Propargyl) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = '[CH2]C([O])(C=C=CC)C(=O)CC(29818)',
    structure = SMILES('[CH2]C([O])(C=C=CC)C(=O)CC'),
    E0 = (177.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,375,552.5,462.5,1710,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,180,180,180,493.254,666.86,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.162025,'amu*angstrom^2'), symmetry=1, barrier=(3.72528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162025,'amu*angstrom^2'), symmetry=1, barrier=(3.72528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162025,'amu*angstrom^2'), symmetry=1, barrier=(3.72528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162025,'amu*angstrom^2'), symmetry=1, barrier=(3.72528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162025,'amu*angstrom^2'), symmetry=1, barrier=(3.72528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162025,'amu*angstrom^2'), symmetry=1, barrier=(3.72528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.18293,0.152929,-0.000246672,2.2249e-07,-7.80837e-11,21558.3,42.3081], Tmin=(100,'K'), Tmax=(850.012,'K')), NASAPolynomial(coeffs=[9.20527,0.0650255,-3.09987e-05,5.84681e-09,-3.98332e-13,20861.9,-3.49115], Tmin=(850.012,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ) + radical(C=CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH2]C(O)(C=C=CC)C(=O)[CH]C(29819)',
    structure = SMILES('[CH2]C(O)(C=C=CC)C([O])=CC'),
    E0 = (64.8645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.33417,0.119613,-0.000121083,6.64879e-08,-1.4965e-11,7991.64,41.8881], Tmin=(100,'K'), Tmax=(1062.63,'K')), NASAPolynomial(coeffs=[17.3808,0.0491647,-2.16385e-05,4.09873e-09,-2.86882e-13,4014.23,-49.5377], Tmin=(1062.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.8645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = '[CH2]CC(=O)C([CH2])(O)C=C=CC(29820)',
    structure = SMILES('[CH2]CC(=O)C([CH2])(O)C=C=CC'),
    E0 = (144.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,180,180,180,180,1600,1681.84,2834.87,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154019,'amu*angstrom^2'), symmetry=1, barrier=(3.54121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.53571,0.155085,-0.000233864,1.93759e-07,-6.36345e-11,17643.8,43.2201], Tmin=(100,'K'), Tmax=(834.899,'K')), NASAPolynomial(coeffs=[15.2006,0.0545066,-2.51291e-05,4.699e-09,-3.19901e-13,15226.1,-35.89], Tmin=(834.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C=C1[CH]C(O)(C1)C(=O)CC(29821)',
    structure = SMILES('[CH2]C=C1[CH]C(O)(C1)C(=O)CC'),
    E0 = (-34.1324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[9.12131,0.0438529,6.21896e-05,-8.68599e-08,2.0735e-11,-4427.17,-9.75556], Tmin=(100,'K'), Tmax=(1810.59,'K')), NASAPolynomial(coeffs=[104.217,0.0164661,-6.64818e-05,1.62491e-08,-1.19721e-12,-68809.6,-607.688], Tmin=(1810.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-34.1324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(536.283,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + ring(methylenecyclobutane) + radical(Allyl_P) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C(O)(C(=O)CC)C1[C]=CC1(29743)',
    structure = SMILES('[CH2]C(O)(C(=O)CC)C1[C]=CC1'),
    E0 = (186.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.86349,0.124539,-0.000126482,6.71913e-08,-1.42546e-11,22604.3,40.9863], Tmin=(100,'K'), Tmax=(1142.02,'K')), NASAPolynomial(coeffs=[22.028,0.0408565,-1.6567e-05,3.02636e-09,-2.08095e-13,17147.5,-77.4488], Tmin=(1142.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CJC(C)(C=O)O) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C1(O)C=C=CCO[C]1CC(29822)',
    structure = SMILES('[CH2]C1(O)C=C=CCO[C]1CC'),
    E0 = (209.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.19,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.660899,0.108445,-0.00010093,5.49534e-08,-1.29085e-11,25300,34.0652], Tmin=(100,'K'), Tmax=(993.502,'K')), NASAPolynomial(coeffs=[11.5684,0.0592076,-2.65887e-05,5.06794e-09,-3.55385e-13,22870.1,-24.8542], Tmin=(993.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C2CsJOCs) + radical(C=CC(C)(O)CJ)"""),
)

species(
    label = 'C#C[CH]C([CH2])(O)C(=O)CC(29755)',
    structure = SMILES('[CH]=C=CC([CH2])(O)C(=O)CC'),
    E0 = (123.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.79541,0.136139,-0.00020414,1.64883e-07,-5.22823e-11,15083.6,38.7498], Tmin=(100,'K'), Tmax=(865.614,'K')), NASAPolynomial(coeffs=[15.2566,0.0432788,-1.8855e-05,3.41366e-09,-2.27319e-13,12658.4,-38.0118], Tmin=(865.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CC(O)(C=O)CJ) + radical(C=C=CJ)"""),
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
    E0 = (79.4579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (166.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (170.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (199.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (273.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (338.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (102.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (108.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (98.9891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (188.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (314.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (293.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (249.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (197.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (349.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (544.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (166.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (112.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (132.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (267.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (142.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (222.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (266.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (217.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (265.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (135.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (166.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (474.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (104.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (524.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (236.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (236.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (87.7422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (452.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (207.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (386.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (216.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (355.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (117.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (267.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (342.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (317.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (271.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (134.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (285.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (205.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (197.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (209.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (505.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C(O)C(=O)CC(4626)', 'CH2CHCCH(26391)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C[C]=CC1(O)CC1([O])CC(29783)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(86.6587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 81.9 to 86.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C1(O)C=C(C=C)C1([O])CC(29784)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(90.7287,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_Nd;radadd_intra_cdsingleDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 85.1 to 90.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C1[C]=CC(O)(C1)C(=O)CC(29785)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.40382e+09,'s^-1'), n=0.352, Ea=(120.148,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(O)(C#CC=C)C(=O)CC(29786)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.66e+08,'cm^3/(mol*s)'), n=1.64, Ea=(12.2591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2703 used for Ct-Cs_Ct-Cd;HJ
Exact match found for rate rule [Ct-Cs_Ct-Cd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]C(O)(C=C=C=C)C(=O)CC(29787)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 193 used for Ca_Ca;HJ
Exact match found for rate rule [Ca_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C2H5CO(71)', 'C=C[C]=CC(=C)O(27798)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-OneDeOs_Cds;CO_rad/NonDe]
Euclidian distance = 3.16227766017
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C=C(4699)', 'C=C(O)C(=O)CC(4626)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-COOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', 'CH2CHCCH(26391)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.0117197,'m^3/(mol*s)'), n=2.33233, Ea=(9.74423,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)', 'C=C[C]=CC(=C)C(=O)CC(29788)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(931.236,'m^3/(mol*s)'), n=1.015, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;OJ_pri] for rate rule [Cds-TwoDe_Cds;OJ_pri]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -7.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(O)([C]=CC=C)C(=O)CC(29789)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C(O)(C=C[C]=C)C(=O)CC(29790)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.806e+09,'s^-1'), n=1.172, Ea=(214.463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[C]=CC(C)([O])C(=O)CC(29791)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C[C]=[C]C(C)(O)C(=O)CC(29792)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=CC=CC([CH2])(O)C(=O)CC(29793)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])(C=CC=C)C(=O)CC(29794)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.71035e+12,'s^-1'), n=1.11009, Ea=(419.408,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Y_rad_out;Cd_H_out_single] for rate rule [R4H_SSD;O_rad_out;Cd_H_out_singleDe]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C[C]=CC(C)(O)C(=O)[CH]C(29795)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C(O)(C=CC=C)C(=O)[CH]C(29796)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_single;Cs_H_out] for rate rule [R5H_DSSS;Cd_rad_out_singleDe_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 4.12310562562
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]CC(=O)C(C)(O)C=[C]C=C(29797)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5H_CCC(O2d)C;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=[C][C]=CC(C)(O)C(=O)CC(29798)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.11562e+07,'s^-1'), n=1.49763, Ea=(188.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_2H;Cd_H_out_doubleC] for rate rule [R5HJ_3;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C=CC=C(29799)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(456.249,'s^-1'), n=2.07385, Ea=(50.4486,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSR;C_rad_out_2H;XH_out] for rate rule [R6H_SSSSD;C_rad_out_2H;Cd_H_out_singleDe]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C[C]=CC(C)(O)C(=O)CC(29800)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_2;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(O)=C([O])CC(4557)', '[CH]=[C]C=C(4699)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C(O)(C=C1[CH]C1)C(=O)CC(29801)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 133 used for R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble
Exact match found for rate rule [R3_D;doublebond_intra_pri_2H;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C[C]=CC1(O)CO[C]1CC(29802)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.01137e+09,'s^-1'), n=0.572544, Ea=(186.425,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S;multiplebond_intra;radadd_intra_cs2H] + [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C1(O)C=C(C=C)O[C]1CC(29803)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.72e+09,'s^-1'), n=0.186, Ea=(55.7727,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleDe] for rate rule [R5_DS_CO;carbonyl_intra_Nd;radadd_intra_cdsingleDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['CCC(=O)C1(O)C=[C][CH]CC1(29804)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.63396e+08,'s^-1'), n=0.483333, Ea=(87.4777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH3CH2OH(54)', '[CH2]C(=C=O)C=[C]C=C(29805)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;Cs_OH] for rate rule [Cdd_Cd;C_pri_OH]
Euclidian distance = 1.41421356237
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C=C=CC(C)(O)C(=O)CC(29806)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH2(S)(23)', '[CH2]C(O)(C=[C]C=C)C(C)=O(29807)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;C_pri/De]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C[C]=C[C](O)CC(=O)CC(29808)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=C[C]=CC[C](O)C(=O)CC(29691)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['C=CC1=CC(O)(C1)C(=O)CC(29809)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C=[C]C=C)=C(CC)OO(29810)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_R]
Euclidian distance = 1.0
family: ketoenol"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C[C]=CO[C](CC)C(=C)O(29687)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O] for rate rule [R_ROR;R1_doublebond;R2_doublebond_CH2CH3;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=C)OC(29811)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(37989.5,'s^-1'), n=2.515, Ea=(258.99,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C] + [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O] for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_C]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(O)=CC(29812)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4235.27,'s^-1'), n=2.8, Ea=(143.114,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_ROR;R1_doublebond;R2_doublebond_CsC;R_O_H] for rate rule [R_ROR;R1_doublebond_CHCH3;R2_doublebond_CsC;R_O_H]
Euclidian distance = 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CH2(19)', 'C=C[C]=C[C](O)C(=O)CC(29813)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/TDMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C=[C]C1CC1(O)C(=O)CC(29814)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(38.4463,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C1(O)C=C=CCC1([O])CC(29815)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0.21, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7_SMMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_SMMS;carbonylbond_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(O)(C=C=[C]C)C(=O)CC(29816)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(606700,'s^-1'), n=2.347, Ea=(214.468,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 119 used for R2H_S;Cd_rad_out_double;Cs_H_out_2H
Exact match found for rate rule [R2H_S;Cd_rad_out_double;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C(O)([C]=C=CC)C(=O)CC(29817)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(0.0029655,'s^-1'), n=4.271, Ea=(238.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SMM;C_rad_out_2H;Cd_H_out_single] for rate rule [R4H_SMM;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([O])(C=C=CC)C(=O)CC(29818)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6H;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C(O)(C=C=CC)C(=O)[CH]C(29819)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1062,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 116 used for R7H;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R7H;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]CC(=O)C([CH2])(O)C=C=CC(29820)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.81e+06,'s^-1'), n=1.5, Ea=(141.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R8H;C_rad_out_2H;Cs_H_out] for rate rule [R8H;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C=C1[CH]C(O)(C1)C(=O)CC(29821)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C(O)(C(=O)CC)C1[C]=CC1(29743)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.53664e+10,'s^-1'), n=0.43543, Ea=(118.482,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs2H] + [R4;doublebond_intra;radadd_intra_cs2H] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    products = ['[CH2]C1(O)C=C=CCO[C]1CC(29822)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.57149e+06,'s^-1'), n=1.0129, Ea=(129.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R7_linear;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 124.6 to 129.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction49',
    reactants = ['CH2(19)', 'C#C[CH]C([CH2])(O)C(=O)CC(29755)'],
    products = ['[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '4714',
    isomers = [
        '[CH2]C(O)(C=[C]C=C)C(=O)CC(29695)',
    ],
    reactants = [
        ('C=C(O)C(=O)CC(4626)', 'CH2CHCCH(26391)'),
        ('[CH2]C(O)=C([O])CC(4557)', 'CH2CHCCH(26391)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '4714',
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

