species(
    label = '[CH2]C(OO)O[C]=O(13097)',
    structure = SMILES('[CH2]C(OO)O[C]=O'),
    E0 = (-125.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,414.613,424.598],'cm^-1')),
        HinderedRotor(inertia=(0.00851514,'amu*angstrom^2'), symmetry=1, barrier=(41.4327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321769,'amu*angstrom^2'), symmetry=1, barrier=(41.3524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78004,'amu*angstrom^2'), symmetry=1, barrier=(41.3883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.64884,'amu*angstrom^2'), symmetry=1, barrier=(18.7019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219237,'amu*angstrom^2'), symmetry=1, barrier=(27.8834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32399,0.0630277,-7.11174e-05,4.24895e-08,-1.04248e-11,-15053.3,25.6899], Tmin=(100,'K'), Tmax=(975.121,'K')), NASAPolynomial(coeffs=[10.5838,0.0250435,-1.26873e-05,2.5421e-09,-1.8312e-13,-16859.2,-18.7498], Tmin=(975.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(CJCOOH)"""),
)

species(
    label = 'CH2CHOOH(64)',
    structure = SMILES('C=COO'),
    E0 = (-53.0705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1310,387.5,850,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.754187,'amu*angstrom^2'), symmetry=1, barrier=(17.3402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.754176,'amu*angstrom^2'), symmetry=1, barrier=(17.34,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3284.22,'J/mol'), sigma=(4.037,'angstroms'), dipoleMoment=(1.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79821,0.0274377,-2.03468e-05,7.62127e-09,-1.12671e-12,-6315.53,9.11829], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.82519,0.0274167,-2.04456e-05,7.77399e-09,-1.18661e-12,-6325.27,8.96641], Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(298,'K'), Tmax=(2000,'K'), E0=(-53.0705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CH2CHOOH""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
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
    label = 'C=C(OO)O[C]=O(15804)',
    structure = SMILES('C=C(OO)O[C]=O'),
    E0 = (-130.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1855,455,950,3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,222.798,222.809],'cm^-1')),
        HinderedRotor(inertia=(1.04934,'amu*angstrom^2'), symmetry=1, barrier=(36.961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04921,'amu*angstrom^2'), symmetry=1, barrier=(36.9615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256062,'amu*angstrom^2'), symmetry=1, barrier=(9.02141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349599,'amu*angstrom^2'), symmetry=1, barrier=(12.315,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922268,0.0757603,-0.000130508,1.17968e-07,-4.14674e-11,-15579.6,22.6287], Tmin=(100,'K'), Tmax=(822.419,'K')), NASAPolynomial(coeffs=[8.40412,0.0255046,-1.35568e-05,2.66498e-09,-1.85687e-13,-16341.3,-9.15321], Tmin=(822.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[CH2][CH]OO(104)',
    structure = SMILES('[CH2][CH]OO'),
    E0 = (224.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.00920734,'amu*angstrom^2'), symmetry=1, barrier=(3.53679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00921023,'amu*angstrom^2'), symmetry=1, barrier=(3.53685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43223,'amu*angstrom^2'), symmetry=1, barrier=(32.9297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52392,0.0279389,-1.79903e-05,3.58397e-09,4.58838e-13,27095.5,18.6054], Tmin=(100,'K'), Tmax=(1150.47,'K')), NASAPolynomial(coeffs=[9.41961,0.0107482,-4.42273e-06,8.47943e-10,-6.05179e-14,25059.9,-17.5802], Tmin=(1150.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CCsJOOH) + radical(CJCOOH)"""),
)

species(
    label = 'HO2(10)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26738e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.913,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11274e-07,6.34957e-11,-4.86385e-15,83.4208,3.09341], Tmin=(923.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=CO[C]=O(735)',
    structure = SMILES('C=CO[C]=O'),
    E0 = (-97.7509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43009,'amu*angstrom^2'), symmetry=1, barrier=(32.8805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42165,'amu*angstrom^2'), symmetry=1, barrier=(32.6866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06242,0.0351364,-1.90366e-05,-6.74125e-09,6.69243e-12,-11680.2,14.7084], Tmin=(100,'K'), Tmax=(965.053,'K')), NASAPolynomial(coeffs=[13.1388,0.00674642,-2.14129e-06,3.99152e-10,-3.05804e-14,-14633.9,-42.5616], Tmin=(965.053,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.7509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = 'C[C](OO)O[C]=O(15805)',
    structure = SMILES('C[C](OO)O[C]=O'),
    E0 = (-134.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46319,0.0605803,-6.66659e-05,3.97411e-08,-9.88641e-12,-16107.2,23.8851], Tmin=(100,'K'), Tmax=(954.492,'K')), NASAPolynomial(coeffs=[9.50409,0.026883,-1.37096e-05,2.75343e-09,-1.9854e-13,-17642.2,-14.533], Tmin=(954.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(Cs_P)"""),
)

species(
    label = '[CH2][C](OO)OC=O(15806)',
    structure = SMILES('[CH2][C](OO)OC=O'),
    E0 = (-117.135,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,360,370,350,499.802,500.022],'cm^-1')),
        HinderedRotor(inertia=(0.476021,'amu*angstrom^2'), symmetry=1, barrier=(32.2989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0423451,'amu*angstrom^2'), symmetry=1, barrier=(7.76638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0161822,'amu*angstrom^2'), symmetry=1, barrier=(32.2981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.38888,'amu*angstrom^2'), symmetry=1, barrier=(32.2984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0161732,'amu*angstrom^2'), symmetry=1, barrier=(32.2986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63968,0.0553172,-5.16139e-05,2.46121e-08,-4.87419e-12,-14006,27.3005], Tmin=(100,'K'), Tmax=(1175.18,'K')), NASAPolynomial(coeffs=[10.5818,0.0248804,-1.27644e-05,2.57313e-09,-1.8576e-13,-16107.7,-17.2836], Tmin=(1175.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(Cs_P) + radical(CJCOOH)"""),
)

species(
    label = 'CC(O[O])O[C]=O(15807)',
    structure = SMILES('CC(O[O])O[C]=O'),
    E0 = (-187.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,492.5,1135,1000,1380,1390,370,380,2900,435,276.68,278.263],'cm^-1')),
        HinderedRotor(inertia=(0.00223495,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172425,'amu*angstrom^2'), symmetry=1, barrier=(9.44106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764272,'amu*angstrom^2'), symmetry=1, barrier=(42.2275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51272,'amu*angstrom^2'), symmetry=1, barrier=(79.6272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59822,0.0573444,-6.47818e-05,4.19259e-08,-1.14714e-11,-22515.5,23.4708], Tmin=(100,'K'), Tmax=(869.836,'K')), NASAPolynomial(coeffs=[8.05629,0.0276469,-1.35703e-05,2.67655e-09,-1.90883e-13,-23639,-6.78517], Tmin=(869.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-187.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(ROOJ) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C(O[O])OC=O(15808)',
    structure = SMILES('[CH2]C(O[O])OC=O'),
    E0 = (-170.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83772,0.0513008,-4.68676e-05,2.29513e-08,-4.7759e-12,-20416.9,25.9701], Tmin=(100,'K'), Tmax=(1109.61,'K')), NASAPolynomial(coeffs=[8.79829,0.0262089,-1.29479e-05,2.57197e-09,-1.8435e-13,-21961.6,-8.33459], Tmin=(1109.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(ROOJ) + radical(CJCOOH)"""),
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
    label = 'C=C(OO)OC=O(13090)',
    structure = SMILES('C=C(OO)OC=O'),
    E0 = (-326.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60542,0.0609934,-6.88754e-05,2.61093e-08,9.60039e-12,-39231.9,21.6275], Tmin=(100,'K'), Tmax=(530.732,'K')), NASAPolynomial(coeffs=[7.37499,0.0298495,-1.57304e-05,3.16161e-09,-2.26324e-13,-40018.1,-4.18968], Tmin=(530.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-326.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-OdOsH)"""),
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
    label = 'O=[C]OC1CO1(15809)',
    structure = SMILES('O=[C]OC1CO1'),
    E0 = (-221.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46051,0.0413135,-3.51799e-05,1.61193e-08,-2.79238e-12,-26576.8,19.4791], Tmin=(100,'K'), Tmax=(1707.35,'K')), NASAPolynomial(coeffs=[9.858,0.012669,-2.1329e-06,1.38064e-10,-1.71159e-15,-28136.8,-21.697], Tmin=(1707.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-221.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsH) + ring(Ethylene_oxide) + radical((O)CJOC)"""),
)

species(
    label = '[CH2]C1OC(=O)O1(15810)',
    structure = SMILES('[CH2]C1OC(=O)O1'),
    E0 = (-287.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69577,0.00732323,8.11705e-05,-1.15292e-07,4.5114e-11,-34568.8,17.7697], Tmin=(100,'K'), Tmax=(963.683,'K')), NASAPolynomial(coeffs=[15.8998,0.00664039,-2.01192e-06,5.33332e-10,-5.26563e-14,-39626.9,-58.4834], Tmin=(963.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-287.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsOs) + ring(Cyclobutane) + radical(CJCO)"""),
)

species(
    label = 'O=C1CC(OO)O1(15811)',
    structure = SMILES('O=C1CC(OO)O1'),
    E0 = (-425.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0535,0.0306954,2.14794e-05,-5.2383e-08,2.32398e-11,-51057.5,21.7196], Tmin=(100,'K'), Tmax=(940.553,'K')), NASAPolynomial(coeffs=[12.454,0.0169099,-5.09105e-06,8.66632e-10,-6.19347e-14,-54360.6,-34.9785], Tmin=(940.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone)"""),
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
    label = '[CH2]C([O])OO(3363)',
    structure = SMILES('[CH2]C([O])OO'),
    E0 = (54.3707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,180],'cm^-1')),
        HinderedRotor(inertia=(0.0808977,'amu*angstrom^2'), symmetry=1, barrier=(1.86,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00660492,'amu*angstrom^2'), symmetry=1, barrier=(38.014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6563,'amu*angstrom^2'), symmetry=1, barrier=(38.0816,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0514,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06051,0.0476767,-7.19216e-05,6.60197e-08,-2.46593e-11,6604.25,21.5552], Tmin=(100,'K'), Tmax=(779.71,'K')), NASAPolynomial(coeffs=[5.02733,0.0246558,-1.2627e-05,2.49025e-09,-1.7554e-13,6378.72,9.50081], Tmin=(779.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.3707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCOOH) + radical(CCOJ)"""),
)

species(
    label = '[O]C(CO)O[C]=O(15812)',
    structure = SMILES('[O]C(CO)O[C]=O'),
    E0 = (-345.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40349,0.0641645,-9.91639e-05,9.15321e-08,-3.3923e-11,-41470.1,26.1943], Tmin=(100,'K'), Tmax=(797.514,'K')), NASAPolynomial(coeffs=[5.30008,0.0324144,-1.64884e-05,3.22937e-09,-2.26288e-13,-41703.5,10.7109], Tmin=(797.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-345.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-OdOsH) + radical(CCOJ) + radical((O)CJOC)"""),
)

species(
    label = '[C]=O(361)',
    structure = SMILES('[C]=O'),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.00200416,-1.61661e-05,2.55058e-08,-1.16424e-11,52802.7,4.52505], Tmin=(100,'K'), Tmax=(856.11,'K')), NASAPolynomial(coeffs=[0.961625,0.00569045,-3.48044e-06,7.19202e-10,-5.08041e-14,53738.7,21.4663], Tmin=(856.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: Klippenstein_Glarborg2016 + radical(CdCdJ2_triplet)"""),
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
    label = 'O=[C]O[CH]OO(8304)',
    structure = SMILES('O=[C]O[CH]OO'),
    E0 = (-105.113,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,1855,455,950,360.194,360.359],'cm^-1')),
        HinderedRotor(inertia=(0.607621,'amu*angstrom^2'), symmetry=1, barrier=(56.0038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.394604,'amu*angstrom^2'), symmetry=1, barrier=(36.3294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212764,'amu*angstrom^2'), symmetry=1, barrier=(19.6007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608009,'amu*angstrom^2'), symmetry=1, barrier=(56.0034,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5953,0.0558388,-7.72888e-05,5.38691e-08,-1.48203e-11,-12558.2,19.4262], Tmin=(100,'K'), Tmax=(890.718,'K')), NASAPolynomial(coeffs=[10.8592,0.0142361,-7.22697e-06,1.42974e-09,-1.01729e-13,-14208.4,-24.1946], Tmin=(890.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-105.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-OsOsHH) + group(Cds-OdOsH) + radical((O)CJOC) + radical(OCJO)"""),
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
    E0 = (-125.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (92.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-125.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-64.9496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (6.35172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (47.6468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-44.7558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-84.1468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (256.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-62.5328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-80.2896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-50.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-117.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-52.2962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-13.3829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (493.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (276.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['CH2CHOOH(64)', 'CO2(13)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(3)', 'C=C(OO)O[C]=O(15804)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsOs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]OO(104)', 'CO2(13)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(23.3993,'m^3/(mol*s)'), n=2.021, Ea=(52.3867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd-O2d;CJ]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 49.8 to 52.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['HO2(10)', 'C=CO[C]=O(735)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.91057e-05,'m^3/(mol*s)'), n=3.20111, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;OJ-O2s] for rate rule [Cds-OsH_Cds;OJ-O2s]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['C[C](OO)O[C]=O(15805)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(309968,'s^-1'), n=2.08546, Ea=(132.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][C](OO)OC=O(15806)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.9172e+08,'s^-1'), n=1.32036, Ea=(164.782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC(O[O])O[C]=O(15807)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.15, Ea=(143.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 246 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['[CH2]C(O[O])OC=O(15808)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;CO_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH]OO(104)', '[O][C]=O(669)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.49215e+07,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['C=C(OO)OC=O(13090)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['OH(5)', 'O=[C]OC1CO1(15809)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.98e+12,'s^-1','*|/',1.2), n=0, Ea=(45.6433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 1 used for R2OO_S;C_pri_rad_intra;OOH
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['OH(5)', '[CH2]C1OC(=O)O1(15810)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO;Y_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['O=C1CC(OO)O1(15811)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CO(12)', '[CH2]C([O])OO(3363)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.41e+07,'cm^3/(mol*s)'), n=0, Ea=(12.552,'kJ/mol'), T0=(1,'K'), Tmin=(250,'K'), Tmax=(2500,'K'), comment="""From training reaction 9 used for COm;O_rad/NonDe
Exact match found for rate rule [COm;O_rad/NonDe]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(OO)O[C]=O(13097)'],
    products = ['[O]C(CO)O[C]=O(15812)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.39e+11,'s^-1'), n=0, Ea=(112.55,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 1 used for R2OOH_S;C_rad_out_2H
Exact match found for rate rule [R2OOH_S;C_rad_out_2H]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]=O(361)', '[CH2]C([O])OO(3363)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1355.7,'m^3/(mol*s)'), n=1.40819, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -12.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(19)', 'O=[C]O[CH]OO(8304)'],
    products = ['[CH2]C(OO)O[C]=O(13097)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/O2;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '2411',
    isomers = [
        '[CH2]C(OO)O[C]=O(13097)',
    ],
    reactants = [
        ('CH2CHOOH(64)', 'CO2(13)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '2411',
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

