species(
    label = '[CH2]C(OO)C([O])=O(13095)',
    structure = SMILES('[CH2]C(OO)C([O])=O'),
    E0 = (-127.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92832,0.0488518,-3.73716e-05,1.40741e-08,-2.21429e-12,-15289.7,25.8672], Tmin=(100,'K'), Tmax=(1417.84,'K')), NASAPolynomial(coeffs=[10.1489,0.0256601,-1.28362e-05,2.53763e-09,-1.80141e-13,-17620.8,-16.6624], Tmin=(1417.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJCOOH) + radical(CCOJ)"""),
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
    label = '[O]C1([O])CC1OO(15813)',
    structure = SMILES('[O]C1([O])CC1OO'),
    E0 = (20.4502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93686,0.047813,-3.60231e-05,1.34719e-08,-2.10244e-12,2531.32,21.4641], Tmin=(100,'K'), Tmax=(1433.04,'K')), NASAPolynomial(coeffs=[10.093,0.025047,-1.21934e-05,2.38606e-09,-1.6847e-13,193.708,-20.8191], Tmin=(1433.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.4502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + ring(Cyclopropane) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'C=C(OO)C([O])=O(15814)',
    structure = SMILES('C=C(OO)C([O])=O'),
    E0 = (-110.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3615,1310,387.5,850,1000,2950,3100,1380,975,1025,1650,180,180,180,1437.64,1437.8,1437.98],'cm^-1')),
        HinderedRotor(inertia=(0.220865,'amu*angstrom^2'), symmetry=1, barrier=(5.07813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220854,'amu*angstrom^2'), symmetry=1, barrier=(5.07786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221049,'amu*angstrom^2'), symmetry=1, barrier=(5.08235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (103.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44362,0.0684653,-0.000130679,1.32816e-07,-5.0585e-11,-13174,24.0215], Tmin=(100,'K'), Tmax=(840.366,'K')), NASAPolynomial(coeffs=[2.8336,0.0336293,-1.81281e-05,3.57012e-09,-2.48202e-13,-12411.2,23.4861], Tmin=(840.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-110.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ)"""),
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
    label = 'tSPC_1553(806)',
    structure = SMILES('C=CC([O])=O'),
    E0 = (-95.7795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,519.883,519.884,519.886,519.889,519.893,519.897],'cm^-1')),
        HinderedRotor(inertia=(0.000623705,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66944,0.0212367,8.85353e-06,-2.77913e-08,1.22625e-11,-11464.5,15.0665], Tmin=(100,'K'), Tmax=(983.401,'K')), NASAPolynomial(coeffs=[10.0894,0.0103164,-3.86785e-06,7.48921e-10,-5.61021e-14,-13855.2,-25.3414], Tmin=(983.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.7795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""tSPC_1553""", comment="""Thermo library: CBS_QB3_1dHR"""),
)

species(
    label = 'C[C](OO)C([O])=O(15815)',
    structure = SMILES('CC(OO)=C([O])[O]'),
    E0 = (-166.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16537,0.0667239,-8.72113e-05,6.08965e-08,-1.71789e-11,-19980.5,23.9206], Tmin=(100,'K'), Tmax=(861.492,'K')), NASAPolynomial(coeffs=[10.5405,0.0231953,-1.14221e-05,2.24807e-09,-1.59766e-13,-21595.8,-19.9112], Tmin=(861.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C](OO)C(=O)O(15816)',
    structure = SMILES('[CH2][C](OO)C(=O)O'),
    E0 = (-166.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64611,0.0546502,-5.15659e-05,2.54125e-08,-5.21218e-12,-19945.9,27.7119], Tmin=(100,'K'), Tmax=(1140.09,'K')), NASAPolynomial(coeffs=[10.1558,0.0247936,-1.22836e-05,2.44196e-09,-1.75117e-13,-21886.3,-14.458], Tmin=(1140.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-166.524,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJCOOH) + radical(C2CsJOOH)"""),
)

species(
    label = 'CC(O[O])C([O])=O(15817)',
    structure = SMILES('CC(O[O])C([O])=O'),
    E0 = (-189.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,492.5,1135,1000,1380,1390,370,380,2900,435,180,180,919.874,921.772,922.176,923.489],'cm^-1')),
        HinderedRotor(inertia=(0.00513639,'amu*angstrom^2'), symmetry=1, barrier=(3.0973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133453,'amu*angstrom^2'), symmetry=1, barrier=(3.06835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.33966,'amu*angstrom^2'), symmetry=1, barrier=(53.7933,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37071,0.0414259,-2.57823e-05,7.4754e-09,-8.9296e-13,-22760.1,23.0247], Tmin=(100,'K'), Tmax=(1750.91,'K')), NASAPolynomial(coeffs=[9.22923,0.0257573,-1.2359e-05,2.36438e-09,-1.63189e-13,-25161.8,-13.9054], Tmin=(1750.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(O[O])C(=O)O(15818)',
    structure = SMILES('[CH2]C(O[O])C(=O)O'),
    E0 = (-201.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38297,0.0531685,-4.75301e-05,2.12559e-08,-3.78799e-12,-24127.3,28.4254], Tmin=(100,'K'), Tmax=(1344.43,'K')), NASAPolynomial(coeffs=[13.4366,0.0173061,-7.51796e-06,1.41502e-09,-9.85489e-14,-27368.4,-33.2941], Tmin=(1344.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJCOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C(OO)[C]1OO1(15819)',
    structure = SMILES('[CH2]C(OO)[C]1OO1'),
    E0 = (231.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332868,0.0662112,-6.98335e-05,3.5215e-08,-6.65724e-12,28032.2,29.7416], Tmin=(100,'K'), Tmax=(1465.47,'K')), NASAPolynomial(coeffs=[18.9552,0.00712898,-9.1199e-07,1.87463e-11,2.56382e-15,23460.3,-64.1937], Tmin=(1465.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(dioxirane) + radical(Cs_P) + radical(CJCOOH)"""),
)

species(
    label = '[O][C]1OCC1OO(15820)',
    structure = SMILES('[O][C]1OCC1OO'),
    E0 = (33.4356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65043,0.0502884,-4.36816e-05,1.88748e-08,-2.15354e-12,4107.36,23.7379], Tmin=(100,'K'), Tmax=(852.301,'K')), NASAPolynomial(coeffs=[9.27694,0.0221605,-7.66769e-06,1.25625e-09,-8.05853e-14,2528.95,-13.4699], Tmin=(852.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.4356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Oxetane) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'C=C(OO)C(=O)O(13092)',
    structure = SMILES('C=C(OO)C(=O)O'),
    E0 = (-335.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11947,0.0707328,-0.000116563,1.05149e-07,-3.72786e-11,-40303.7,25.0196], Tmin=(100,'K'), Tmax=(816.041,'K')), NASAPolynomial(coeffs=[7.42963,0.0273808,-1.40431e-05,2.74279e-09,-1.9109e-13,-40920,-1.60638], Tmin=(816.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
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
    label = '[O]C(=O)C1CO1(15821)',
    structure = SMILES('[O]C(=O)C1CO1'),
    E0 = (-223.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76084,0.0187749,2.84266e-05,-5.17062e-08,2.24743e-11,-26843,17.1715], Tmin=(100,'K'), Tmax=(893.748,'K')), NASAPolynomial(coeffs=[8.56842,0.0153479,-3.69332e-06,5.02031e-10,-3.13276e-14,-28782.3,-15.2362], Tmin=(893.748,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Ethylene_oxide) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C1OOC1=O(15822)',
    structure = SMILES('[CH2]C1OOC1=O'),
    E0 = (-99.1481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.76697,-0.0022377,0.000125492,-1.71216e-07,6.73607e-11,-11853.5,20.626], Tmin=(100,'K'), Tmax=(945.699,'K')), NASAPolynomial(coeffs=[19.7729,0.000791902,1.79134e-06,-1.99674e-10,-5.13069e-15,-18422,-78.1906], Tmin=(945.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.1481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclobutane) + radical(CJCOOH)"""),
)

species(
    label = '[O]C(=O)C[CH]OO(13094)',
    structure = SMILES('[O]C(=O)C[CH]OO'),
    E0 = (-143.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33243,0.0703082,-0.000123034,1.24925e-07,-4.85233e-11,-17119.5,24.0516], Tmin=(100,'K'), Tmax=(826.663,'K')), NASAPolynomial(coeffs=[1.96287,0.040774,-2.13883e-05,4.19794e-09,-2.92778e-13,-16318.9,26.6034], Tmin=(826.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CCsJOOH)"""),
)

species(
    label = 'O=C1OCC1OO(15802)',
    structure = SMILES('O=C1OCC1OO'),
    E0 = (-394.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22482,0.0193299,6.6156e-05,-1.0788e-07,4.53761e-11,-47406.3,22.7597], Tmin=(100,'K'), Tmax=(922.771,'K')), NASAPolynomial(coeffs=[16.1578,0.00969191,-6.86448e-07,2.08912e-11,-6.46877e-15,-52138.7,-55.0483], Tmin=(922.771,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-394.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone)"""),
)

species(
    label = '[O]C(=O)C([O])CO(15823)',
    structure = SMILES('[O]C(=O)C([O])CO'),
    E0 = (-329.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06545,0.0471182,-4.49458e-05,2.71914e-08,-7.57185e-12,-39537.6,26.1265], Tmin=(100,'K'), Tmax=(822.376,'K')), NASAPolynomial(coeffs=[5.43088,0.0307488,-1.50883e-05,2.9871e-09,-2.13791e-13,-40091.1,10.5484], Tmin=(822.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C=C([O])OOO(15824)',
    structure = SMILES('[CH2][CH]C(=O)OOO'),
    E0 = (-23.7023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1056.86,1056.86,1056.86,1056.86,1056.86,1056.86,1056.86,1056.86,1056.86,1056.86,2258.2],'cm^-1')),
        HinderedRotor(inertia=(0.0812829,'amu*angstrom^2'), symmetry=1, barrier=(1.86885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0812829,'amu*angstrom^2'), symmetry=1, barrier=(1.86885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0812829,'amu*angstrom^2'), symmetry=1, barrier=(1.86885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0812829,'amu*angstrom^2'), symmetry=1, barrier=(1.86885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0812829,'amu*angstrom^2'), symmetry=1, barrier=(1.86885,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27228,0.0464284,-9.1011e-06,-3.10281e-08,1.76011e-11,-2740.21,30.6582], Tmin=(100,'K'), Tmax=(956.997,'K')), NASAPolynomial(coeffs=[18.0198,0.00939133,-2.71582e-06,5.15979e-10,-4.1716e-14,-7455.14,-57.2892], Tmin=(956.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.7023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJCC=O) + radical(CCJCO)"""),
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
    label = '[O]C(=O)[CH]OO(8453)',
    structure = SMILES('[O]C(=O)[CH]OO'),
    E0 = (-153.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1310,387.5,850,1000,493.634,493.634,493.634,493.634,493.634,493.635],'cm^-1')),
        HinderedRotor(inertia=(0.000691815,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000691815,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.297774,'amu*angstrom^2'), symmetry=1, barrier=(51.4902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0349,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11716,0.0468317,-7.09665e-05,6.8055e-08,-2.68166e-11,-18357.3,19.4694], Tmin=(100,'K'), Tmax=(758.629,'K')), NASAPolynomial(coeffs=[4.10804,0.027342,-1.46501e-05,2.94037e-09,-2.09485e-13,-18400.6,12.12], Tmin=(758.629,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(OCJC=O)"""),
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
    label = '[CH2]C([C]=O)OO(15825)',
    structure = SMILES('[CH2]C([C]=O)OO'),
    E0 = (84.0428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3615,1310,387.5,850,1000,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435],'cm^-1')),
        HinderedRotor(inertia=(0.0683778,'amu*angstrom^2'), symmetry=1, barrier=(9.49386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.619144,'amu*angstrom^2'), symmetry=1, barrier=(14.2353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.619678,'amu*angstrom^2'), symmetry=1, barrier=(14.2476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618106,'amu*angstrom^2'), symmetry=1, barrier=(14.2115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (88.0621,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14167,0.0658139,-9.54084e-05,6.99681e-08,-2.00304e-11,10208.3,24.4411], Tmin=(100,'K'), Tmax=(861.374,'K')), NASAPolynomial(coeffs=[12.121,0.0148286,-6.62265e-06,1.25183e-09,-8.66531e-14,8316.82,-26.8896], Tmin=(861.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.0428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCCJ=O) + radical(CJCOOH)"""),
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
    E0 = (-127.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (20.4502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (111.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-12.4024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-124.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-62.9782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (4.56619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-16.6739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-46.5413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-68.9358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (256.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (232.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (33.4356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-64.3183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-82.0752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-51.8626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-4.96827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-119.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-15.1685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (290.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (228.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (327.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['CH2CHOOH(64)', 'CO2(13)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[O]C1([O])CC1OO(15813)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34238e+09,'s^-1'), n=0.889391, Ea=(148.169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic
Ea raised from 147.9 to 148.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(3)', 'C=C(OO)C([O])=O(15814)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-COOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2CHOOH(64)', '[O][C]=O(669)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.30847,'m^3/(mol*s)'), n=2.06429, Ea=(8.71744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds_Cds;CJ] + [Cds-OsH_Cds;YJ] for rate rule [Cds-OsH_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]OO(104)', 'CO2(13)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.04,'m^3/(mol*s)'), n=1.68, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;CJ] for rate rule [CO2;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['HO2(10)', 'tSPC_1553(806)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.91057e-05,'m^3/(mol*s)'), n=3.20111, Ea=(30.1248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;OJ-O2s] for rate rule [Cds-COH_Cds;OJ-O2s]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['C[C](OO)C([O])=O(15815)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(309968,'s^-1'), n=2.08546, Ea=(132.285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_2H;Cs_H_out] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[CH2][C](OO)C(=O)O(15816)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.99823e+07,'s^-1'), n=1.57622, Ea=(111.045,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC(O[O])C([O])=O(15817)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.15, Ea=(143.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 246 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[CH2]C(O[O])C(=O)O(15818)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.74454e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]OO(104)', '[O][C]=O(669)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[CH2]C(OO)[C]1OO1(15819)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(360.602,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[O][C]1OCC1OO(15820)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.54148e+08,'s^-1'), n=0.924088, Ea=(161.154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 160.7 to 161.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['C=C(OO)C(=O)O(13092)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['OH(5)', '[O]C(=O)C1CO1(15821)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.98e+12,'s^-1','*|/',1.2), n=0, Ea=(45.6433,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 1 used for R2OO_S;C_pri_rad_intra;OOH
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['OH(5)', '[CH2]C1OOC1=O(15822)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.22709e+11,'s^-1'), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3OO;Y_rad_intra;OOH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[O]C(=O)C[CH]OO(13094)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;CO] for rate rule [cCs(-HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['O=C1OCC1OO(15802)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(OO)C([O])=O(13095)'],
    products = ['[O]C(=O)C([O])CO(15823)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.39e+11,'s^-1'), n=0, Ea=(112.55,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 1 used for R2OOH_S;C_rad_out_2H
Exact match found for rate rule [R2OOH_S;C_rad_out_2H]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C=C([O])OOO(15824)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH2(19)', '[O]C(=O)[CH]OO(8453)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH2]C([C]=O)OO(15825)'],
    products = ['[CH2]C(OO)C([O])=O(13095)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

network(
    label = '2409',
    isomers = [
        '[CH2]C(OO)C([O])=O(13095)',
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
    label = '2409',
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

