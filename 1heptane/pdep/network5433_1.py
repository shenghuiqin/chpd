species(
    label = 'C=C(O)C(O)([CH]C)C[CH]C(6977)',
    structure = SMILES('C=C(O)C(O)([CH]C)C[CH]C'),
    E0 = (-129.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,538.322,600.994,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99999,0.117727,-0.000108969,5.28002e-08,-1.01255e-11,-15371.2,47.435], Tmin=(100,'K'), Tmax=(1270.04,'K')), NASAPolynomial(coeffs=[23.9978,0.0358442,-1.2258e-05,2.03341e-09,-1.321e-13,-21974.7,-84.2032], Tmin=(1270.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = 'C3H6(72)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36194e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35108e-06,1.1662e-09,-8.2762e-14,-487.138,-4.54468], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(O)C(O)=CC(5562)',
    structure = SMILES('C=C(O)C(O)=CC'),
    E0 = (-369.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.87493,'amu*angstrom^2'), symmetry=1, barrier=(20.1164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.877997,'amu*angstrom^2'), symmetry=1, barrier=(20.1869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.874595,'amu*angstrom^2'), symmetry=1, barrier=(20.1087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.875946,'amu*angstrom^2'), symmetry=1, barrier=(20.1397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4302.09,'J/mol'), sigma=(6.83849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.98 K, Pc=30.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.884099,0.0877395,-9.34163e-05,4.76671e-08,-9.0713e-12,-44294.8,25.8217], Tmin=(100,'K'), Tmax=(1473.72,'K')), NASAPolynomial(coeffs=[23.5689,0.00887539,-4.29798e-07,-1.49469e-10,1.60597e-14,-50145.5,-97.0299], Tmin=(1473.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-369.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH2]C1(O)C(C)C1(O)C[CH]C(32421)',
    structure = SMILES('[CH2]C1(O)C(C)C1(O)C[CH]C'),
    E0 = (-60.3158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.87874,0.121324,-0.000117332,6.01743e-08,-1.23289e-11,-7035.9,40.3947], Tmin=(100,'K'), Tmax=(1183.25,'K')), NASAPolynomial(coeffs=[21.9406,0.0408038,-1.52588e-05,2.66539e-09,-1.78521e-13,-12672.9,-78.5281], Tmin=(1183.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.3158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CJC(C)2O) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1(O)C(C)CC1(O)[CH]C(32265)',
    structure = SMILES('[CH2]C1(O)C(C)CC1(O)[CH]C'),
    E0 = (-59.2794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66383,0.107711,-6.60119e-05,5.05326e-10,9.93188e-12,-6910.73,39.3975], Tmin=(100,'K'), Tmax=(982.617,'K')), NASAPolynomial(coeffs=[24.8892,0.0362157,-1.27361e-05,2.26165e-09,-1.58058e-13,-13895.7,-97.2297], Tmin=(982.617,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.2794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CJC(C)2O) + radical(CCJCO)"""),
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
    label = 'C=C(O)C(O)([CH]C)C=CC(32422)',
    structure = SMILES('C=C(O)C(O)([CH]C)C=CC'),
    E0 = (-214.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,180,180,180,180,1600,1668.43,2825.49,3200],'cm^-1')),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153126,'amu*angstrom^2'), symmetry=1, barrier=(3.52067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.45543,0.104552,-7.2216e-05,1.31878e-08,3.9034e-12,-25626.4,42.0119], Tmin=(100,'K'), Tmax=(1026.4,'K')), NASAPolynomial(coeffs=[23.7306,0.0353442,-1.33747e-05,2.44399e-09,-1.71693e-13,-32321.3,-87.5799], Tmin=(1026.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCJCO)"""),
)

species(
    label = 'C=CC(O)(C[CH]C)C(=C)O(32423)',
    structure = SMILES('C=CC(O)(C[CH]C)C(=C)O'),
    E0 = (-208.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,433.478,704.668,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155782,'amu*angstrom^2'), symmetry=1, barrier=(3.58173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53479,0.110078,-9.62385e-05,4.35843e-08,-7.888e-12,-24810.7,42.9801], Tmin=(100,'K'), Tmax=(1329.21,'K')), NASAPolynomial(coeffs=[22.5478,0.0376061,-1.44546e-05,2.56551e-09,-1.73115e-13,-31212.9,-80.0583], Tmin=(1329.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-208.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(RCCJC)"""),
)

species(
    label = 'C=CCC(O)([CH]C)C(=C)O(32424)',
    structure = SMILES('C=CCC(O)([CH]C)C(=C)O'),
    E0 = (-196.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3580,3650,1210,1345,900,1100,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,546.606,587.881,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158051,'amu*angstrom^2'), symmetry=1, barrier=(3.63391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.58859,0.106256,-6.8238e-05,2.11116e-10,1.12569e-11,-23400.8,43.6923], Tmin=(100,'K'), Tmax=(958.782,'K')), NASAPolynomial(coeffs=[25.9603,0.0300807,-9.69786e-06,1.66796e-09,-1.16424e-13,-30464.9,-97.3457], Tmin=(958.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJCO)"""),
)

species(
    label = 'C3H6(T)(143)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238389,'amu*angstrom^2'), symmetry=1, barrier=(5.48103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00909639,'amu*angstrom^2'), symmetry=1, barrier=(22.1005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93778,0.0190991,4.26842e-06,-1.44873e-08,5.74941e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.81,'K')), NASAPolynomial(coeffs=[5.93909,0.0171892,-6.69152e-06,1.21546e-09,-8.39795e-14,33151.2,-4.14888], Tmin=(1046.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2COH(99)',
    structure = SMILES('C=[C]O'),
    E0 = (103.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3615,1277.5,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.989114,'amu*angstrom^2'), symmetry=1, barrier=(22.7417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1624,0.0134245,5.56346e-06,-1.95511e-08,9.36369e-12,12455.2,10.1544], Tmin=(100,'K'), Tmax=(925.618,'K')), NASAPolynomial(coeffs=[8.19875,0.00453462,-8.93448e-07,1.26083e-10,-9.46513e-15,10971.3,-16.733], Tmin=(925.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CH2COH""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH]CC(O)=CC(32425)',
    structure = SMILES('C[CH]CC(O)=CC'),
    E0 = (-96.4399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,280.483,280.491],'cm^-1')),
        HinderedRotor(inertia=(0.0021403,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173169,'amu*angstrom^2'), symmetry=1, barrier=(9.66566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173159,'amu*angstrom^2'), symmetry=1, barrier=(9.66481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173129,'amu*angstrom^2'), symmetry=1, barrier=(9.66408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.309929,0.0723519,-5.4739e-05,2.21173e-08,-3.64838e-12,-11458.9,30.1739], Tmin=(100,'K'), Tmax=(1430.67,'K')), NASAPolynomial(coeffs=[15.0456,0.0311527,-1.15433e-05,1.98897e-09,-1.31091e-13,-15675.3,-46.1947], Tmin=(1430.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.4399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(O)=C(O)[CH]C(4609)',
    structure = SMILES('[CH2]C(O)=C(O)[CH]C'),
    E0 = (-122.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,322.653],'cm^-1')),
        HinderedRotor(inertia=(0.00160507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208419,'amu*angstrom^2'), symmetry=1, barrier=(15.514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20977,'amu*angstrom^2'), symmetry=1, barrier=(15.5164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209392,'amu*angstrom^2'), symmetry=1, barrier=(15.5175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207973,'amu*angstrom^2'), symmetry=1, barrier=(15.5093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141743,0.0688338,-2.99459e-05,-3.072e-08,2.30897e-11,-14621,29.0351], Tmin=(100,'K'), Tmax=(905.863,'K')), NASAPolynomial(coeffs=[24.2058,0.00628117,1.26066e-06,-4.23747e-10,2.91719e-14,-20774,-94.5791], Tmin=(905.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(CCJCO)"""),
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
    label = 'C=C(O)C(=CC)C[CH]C(32426)',
    structure = SMILES('C=C(O)C(=CC)C[CH]C'),
    E0 = (-45.5673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2950,3100,1380,975,1025,1650,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.188,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.06733,0.103515,-9.25226e-05,4.42391e-08,-8.53579e-12,-5291.42,35.726], Tmin=(100,'K'), Tmax=(1245.46,'K')), NASAPolynomial(coeffs=[19.0546,0.0388894,-1.46881e-05,2.57552e-09,-1.72574e-13,-10303.6,-65.7671], Tmin=(1245.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.5673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC)"""),
)

species(
    label = 'C=C(O)C(O)([CH]C)[CH]CC(32301)',
    structure = SMILES('C=C(O)C(O)([CH]C)[CH]CC'),
    E0 = (-124.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,542.312,590.022,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158525,'amu*angstrom^2'), symmetry=1, barrier=(3.6448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.71335,0.108969,-6.97089e-05,9.90927e-10,1.09882e-11,-14722.3,46.1622], Tmin=(100,'K'), Tmax=(960.201,'K')), NASAPolynomial(coeffs=[25.8899,0.0325671,-1.06364e-05,1.82967e-09,-1.27032e-13,-21802.1,-95.1496], Tmin=(960.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CC(O)(C[CH]C)C(=C)O(6980)',
    structure = SMILES('[CH2]CC(O)(C[CH]C)C(=C)O'),
    E0 = (-124.355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,498.316,642.516,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8536,0.12016,-0.000116071,5.94474e-08,-1.21461e-11,-14738.5,46.1957], Tmin=(100,'K'), Tmax=(1187.86,'K')), NASAPolynomial(coeffs=[21.9915,0.0398639,-1.4675e-05,2.54033e-09,-1.69249e-13,-20403.4,-72.9481], Tmin=(1187.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC(O)([CH]C)C(=C)O(32427)',
    structure = SMILES('[CH2]CCC(O)([CH]C)C(=C)O'),
    E0 = (-118.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,544.683,589.321,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88504,0.115086,-8.93716e-05,2.34698e-08,2.40634e-12,-14075.6,46.0683], Tmin=(100,'K'), Tmax=(978.998,'K')), NASAPolynomial(coeffs=[25.2331,0.034376,-1.18128e-05,2.0494e-09,-1.40726e-13,-20827.3,-91.5499], Tmin=(978.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C([O])(CC)C[CH]C(6359)',
    structure = SMILES('C=C(O)C([O])(CC)C[CH]C'),
    E0 = (-100.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1025.12,1193.47,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160955,'amu*angstrom^2'), symmetry=1, barrier=(3.70068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4667.54,'J/mol'), sigma=(7.97399,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=729.06 K, Pc=20.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37847,0.111562,-9.67436e-05,4.49979e-08,-8.52632e-12,-11888.1,43.919], Tmin=(100,'K'), Tmax=(1257.13,'K')), NASAPolynomial(coeffs=[19.371,0.0455405,-1.79667e-05,3.2218e-09,-2.18469e-13,-17105,-60.9332], Tmin=(1257.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=CC(C)2OJ)"""),
)

species(
    label = 'C=C(O)C(O)([CH][CH]C)CC(6978)',
    structure = SMILES('C=C(O)C(O)([CH][CH]C)CC'),
    E0 = (-129.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,180,180,180,538.322,600.994,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157741,'amu*angstrom^2'), symmetry=1, barrier=(3.62678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.99999,0.117727,-0.000108969,5.28002e-08,-1.01255e-11,-15371.2,47.435], Tmin=(100,'K'), Tmax=(1270.04,'K')), NASAPolynomial(coeffs=[23.9978,0.0358442,-1.2258e-05,2.03341e-09,-1.321e-13,-21974.7,-84.2032], Tmin=(1270.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = 'C=C(O)C([O])([CH]C)CCC(6979)',
    structure = SMILES('C=C(O)C([O])([CH]C)CCC'),
    E0 = (-95.0432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1064.56,1147.69,1600,1800,3000,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161769,'amu*angstrom^2'), symmetry=1, barrier=(3.71939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.60905,0.108824,-7.81789e-05,1.95543e-08,1.52878e-12,-11216.5,44.5071], Tmin=(100,'K'), Tmax=(1033.85,'K')), NASAPolynomial(coeffs=[22.9346,0.0394958,-1.47805e-05,2.65393e-09,-1.83543e-13,-17661.2,-81.3439], Tmin=(1033.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-95.0432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C([O])C(O)(CC)C[CH]C(6981)',
    structure = SMILES('C=C([O])C(O)(CC)C[CH]C'),
    E0 = (-191.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.31921,0.113029,-0.000104673,5.33633e-08,-1.11128e-11,-22873,43.993], Tmin=(100,'K'), Tmax=(1150.43,'K')), NASAPolynomial(coeffs=[17.7632,0.0466791,-1.8161e-05,3.2294e-09,-2.18026e-13,-27263.5,-50.7425], Tmin=(1150.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-191.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)C(O)(CC)C[CH]C(6982)',
    structure = SMILES('[CH]=C(O)C(O)(CC)C[CH]C'),
    E0 = (-82.5055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1600,1686.15,2816.21,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154138,'amu*angstrom^2'), symmetry=1, barrier=(3.54393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.94568,0.122732,-0.000121557,6.3712e-08,-1.32909e-11,-9702.21,45.1425], Tmin=(100,'K'), Tmax=(1165.81,'K')), NASAPolynomial(coeffs=[22.2955,0.0395585,-1.45402e-05,2.51448e-09,-1.67467e-13,-15354.3,-75.526], Tmin=(1165.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.5055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]C(O)(CCC)C(=C)O(32428)',
    structure = SMILES('[CH2][CH]C(O)(CCC)C(=C)O'),
    E0 = (-118.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,544.683,589.321,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157746,'amu*angstrom^2'), symmetry=1, barrier=(3.62688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.88504,0.115086,-8.93716e-05,2.34698e-08,2.40634e-12,-14075.6,46.0683], Tmin=(100,'K'), Tmax=(978.998,'K')), NASAPolynomial(coeffs=[25.2331,0.034376,-1.18128e-05,2.0494e-09,-1.40726e-13,-20827.3,-91.5499], Tmin=(978.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C([O])C(O)([CH]C)CCC(32429)',
    structure = SMILES('C=C([O])C(O)([CH]C)CCC'),
    E0 = (-186.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.62131,0.111139,-8.90396e-05,3.15764e-08,-2.5237e-12,-22198.4,44.8367], Tmin=(100,'K'), Tmax=(1028.51,'K')), NASAPolynomial(coeffs=[22.1194,0.039329,-1.42387e-05,2.49044e-09,-1.6908e-13,-28167.3,-75.6425], Tmin=(1028.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-186.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCO)"""),
)

species(
    label = '[CH]=C(O)C(O)([CH]C)CCC(32430)',
    structure = SMILES('[CH]=C(O)C(O)([CH]C)CCC'),
    E0 = (-77.0497,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,180,180,180,180,1600,1709.07,2787.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154673,'amu*angstrom^2'), symmetry=1, barrier=(3.55623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.9976,0.117894,-9.56498e-05,2.87001e-08,8.81081e-13,-9038.5,45.0889], Tmin=(100,'K'), Tmax=(977.04,'K')), NASAPolynomial(coeffs=[25.6961,0.0338095,-1.15312e-05,1.98948e-09,-1.36157e-13,-15848.3,-95.0293], Tmin=(977.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.0497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = '[CH2][CH]CC(O)(CC)C(=C)O(6986)',
    structure = SMILES('[CH2][CH]CC(O)(CC)C(=C)O'),
    E0 = (-124.355,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,180,180,180,498.316,642.516,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15705,'amu*angstrom^2'), symmetry=1, barrier=(3.61089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.8536,0.12016,-0.000116071,5.94474e-08,-1.21461e-11,-14738.5,46.1957], Tmin=(100,'K'), Tmax=(1187.86,'K')), NASAPolynomial(coeffs=[21.9915,0.0398639,-1.4675e-05,2.54033e-09,-1.69249e-13,-20403.4,-72.9481], Tmin=(1187.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC1(O)[C](O)CC1C(32431)',
    structure = SMILES('C[CH]CC1(O)[C](O)CC1C'),
    E0 = (-82.1398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20266,0.105899,-8.48032e-05,3.62771e-08,-6.35229e-12,-9684.9,40.4967], Tmin=(100,'K'), Tmax=(1347.87,'K')), NASAPolynomial(coeffs=[19.056,0.0457789,-1.78976e-05,3.18529e-09,-2.14514e-13,-15146.1,-63.2875], Tmin=(1347.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-82.1398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(557.07,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(C2CsJOH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C1(O)CC(C)C[C]1O(32315)',
    structure = SMILES('C[CH]C1(O)CC(C)C[C]1O'),
    E0 = (-158.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.912574,0.0899723,-2.72154e-05,-3.20648e-08,1.98356e-11,-18864.2,39.4047], Tmin=(100,'K'), Tmax=(976.192,'K')), NASAPolynomial(coeffs=[21.4525,0.0398558,-1.40148e-05,2.4959e-09,-1.74939e-13,-25209.3,-78.0891], Tmin=(976.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(C2CsJOH) + radical(CCJCO)"""),
)

species(
    label = 'C=C(O)C(O)(C=CC)CC(6988)',
    structure = SMILES('C=C(O)C(O)(C=CC)CC'),
    E0 = (-414.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.69317,0.110236,-8.22759e-05,2.41786e-08,-4.03963e-13,-49661.2,40.3709], Tmin=(100,'K'), Tmax=(1068.41,'K')), NASAPolynomial(coeffs=[23.6251,0.0388845,-1.50046e-05,2.73361e-09,-1.90067e-13,-56408.9,-89.7104], Tmin=(1068.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-414.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CC(O)(CCC)C(=C)O(32432)',
    structure = SMILES('C=CC(O)(CCC)C(=C)O'),
    E0 = (-402.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7123,0.109712,-7.67187e-05,1.56904e-08,3.22141e-12,-48186.7,40.8764], Tmin=(100,'K'), Tmax=(1033.79,'K')), NASAPolynomial(coeffs=[24.394,0.0377155,-1.43543e-05,2.62244e-09,-1.83825e-13,-55134.9,-93.4383], Tmin=(1033.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=CCC(O)(CC)C(=C)O(6989)',
    structure = SMILES('C=CCC(O)(CC)C(=C)O'),
    E0 = (-396.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.79539,0.111614,-7.73566e-05,1.03135e-08,7.16615e-12,-47436.9,41.938], Tmin=(100,'K'), Tmax=(973.27,'K')), NASAPolynomial(coeffs=[25.5103,0.0341835,-1.16427e-05,2.03036e-09,-1.40733e-13,-54399.9,-97.5218], Tmin=(973.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = 'C=C(O)[C](O)C(C)C[CH]C(32069)',
    structure = SMILES('[CH2]C(O)=C(O)C(C)C[CH]C'),
    E0 = (-207.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.58566,0.12526,-0.000120971,5.97929e-08,-1.14815e-11,-24716.5,45.5555], Tmin=(100,'K'), Tmax=(1346.14,'K')), NASAPolynomial(coeffs=[27.7437,0.0298575,-8.78173e-06,1.31841e-09,-8.07295e-14,-32403.6,-108.004], Tmin=(1346.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-207.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(RCCJC)"""),
)

species(
    label = 'C=C(O)C(C)[C](O)C[CH]C(32433)',
    structure = SMILES('C=C(O)C(C)[C](O)C[CH]C'),
    E0 = (-135.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50616,0.120642,-0.000122298,6.79737e-08,-1.53768e-11,-16106.5,43.155], Tmin=(100,'K'), Tmax=(1064.36,'K')), NASAPolynomial(coeffs=[17.9186,0.047643,-1.94218e-05,3.53783e-09,-2.42184e-13,-20241.6,-51.7699], Tmin=(1064.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(RCCJC) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(C)C(O)([CH]C)C(=C)O(6787)',
    structure = SMILES('[CH2]C(C)C(O)([CH]C)C(=C)O'),
    E0 = (-124.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3580,3650,1210,1345,900,1100,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    molecularWeight = (142.196,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4683.75,'J/mol'), sigma=(8.00266,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=731.59 K, Pc=20.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.80199,0.124005,-0.000116902,5.59015e-08,-1.02595e-11,-14738.2,49.0896], Tmin=(100,'K'), Tmax=(1476.78,'K')), NASAPolynomial(coeffs=[28.6857,0.0271586,-6.79237e-06,8.94335e-10,-5.0253e-14,-22777.8,-110.829], Tmin=(1476.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = 'C=C(O)C1(O)CC(C)C1C(32073)',
    structure = SMILES('C=C(O)C1(O)CC(C)C1C'),
    E0 = (-383.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05364,0.0828225,1.85882e-05,-9.71294e-08,4.72252e-11,-45971.1,36.4239], Tmin=(100,'K'), Tmax=(942.306,'K')), NASAPolynomial(coeffs=[28.5436,0.0281883,-7.46857e-06,1.26922e-09,-9.55535e-14,-54701.4,-121.334], Tmin=(942.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-383.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(561.227,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'C[CH]CC(O)([CH]C)C(C)=O(32434)',
    structure = SMILES('C[CH]CC(O)([CH]C)C(C)=O'),
    E0 = (-138.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.196,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.885994,0.111205,-0.000107255,5.99795e-08,-1.41691e-11,-16542.1,42.5425], Tmin=(100,'K'), Tmax=(1001.74,'K')), NASAPolynomial(coeffs=[13.0713,0.0554715,-2.37983e-05,4.43731e-09,-3.07336e-13,-19338.3,-24.8173], Tmin=(1001.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(548.755,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(RCCJC) + radical(CCJCO)"""),
)

species(
    label = 'CHCH3(T)(95)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438699,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82363,-0.000909515,3.2138e-05,-3.7348e-08,1.3309e-11,41371.4,7.10948], Tmin=(100,'K'), Tmax=(960.812,'K')), NASAPolynomial(coeffs=[4.30487,0.00943069,-3.27559e-06,5.95121e-10,-4.27307e-14,40709.1,1.84202], Tmin=(960.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(O)[C](O)C[CH]C(32435)',
    structure = SMILES('[CH2]C(O)=C(O)C[CH]C'),
    E0 = (-152.082,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.14433,0.0966464,-9.86265e-05,5.05028e-08,-9.87759e-12,-18091.9,36.134], Tmin=(100,'K'), Tmax=(1378.16,'K')), NASAPolynomial(coeffs=[23.1209,0.0180387,-4.16648e-06,5.02428e-10,-2.62427e-14,-24003.4,-85.8967], Tmin=(1378.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.082,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(O)CJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(O)([CH]C)C(=C)O(32436)',
    structure = SMILES('[CH2]C(O)([CH]C)C(=C)O'),
    E0 = (-63.1419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3580,3650,1210,1345,900,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.2712,0.096693,-9.67083e-05,4.78631e-08,-9.00838e-12,-7388.02,38.4738], Tmin=(100,'K'), Tmax=(1443.67,'K')), NASAPolynomial(coeffs=[24.3371,0.0164422,-3.66614e-06,4.36706e-10,-2.3063e-14,-13813.1,-91.1192], Tmin=(1443.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-63.1419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CC(C)(O)CJ) + radical(CCJCO)"""),
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
    E0 = (-129.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-59.2683,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-59.2794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (2.70913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (6.94646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (17.9389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-69.4183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (14.6474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-104.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-1.28046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (34.7482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (23.3399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (28.7957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-24.1376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-6.66656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-41.1533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-45.6012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-38.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-72.4571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-100.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-38.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-77.9129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (162.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-5.43481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-57.7348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-66.2995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-121.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-121.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (27.6188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (27.6188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (35.1773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-121.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (74.4794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (191.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (280.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C3H6(72)', 'C=C(O)C(O)=CC(5562)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['[CH2]C1(O)C(C)C1(O)C[CH]C(32421)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.39513e+10,'s^-1'), n=0.560608, Ea=(70.4313,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd] + [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R4_S_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['[CH2]C1(O)C(C)CC1(O)[CH]C(32265)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.53628e+06,'s^-1'), n=1.28, Ea=(70.4202,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra_2H;radadd_intra_csHNd] + [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_2H_secNd;radadd_intra_csHNd]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 69.9 to 70.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', 'C=C(O)C(O)([CH]C)C=CC(32422)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C=CC(O)(C[CH]C)C(=C)O(32423)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.014e+08,'cm^3/(mol*s)'), n=1.733, Ea=(3.17984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2823 used for Cds-HH_Cds-Cs\O2s/H;HJ
Exact match found for rate rule [Cds-HH_Cds-Cs\O2s/H;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=CCC(O)([CH]C)C(=C)O(32424)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C3H6(T)(143)', 'C=C(O)C(O)=CC(5562)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00392164,'m^3/(mol*s)'), n=2.41519, Ea=(15.6067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;CJ] for rate rule [Cds-CdOs_Cds;CJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2COH(99)', 'C[CH]CC(O)=CC(32425)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00712612,'m^3/(mol*s)'), n=2.40979, Ea=(7.81798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ] for rate rule [Cds-OsCs_Cds;CdsJ-O2s]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C3H6(72)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(12.4666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(5)', 'C=C(O)C(=CC)C[CH]C(32426)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.567844,'m^3/(mol*s)'), n=2.025, Ea=(15.9149,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdCs_Cds-CsH;YJ] for rate rule [Cds-CdCs_Cds-CsH;OJ_pri]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C(O)C(O)([CH]C)[CH]CC(32301)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC(O)(C[CH]C)C(=C)O(6980)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CCC(O)([CH]C)C(=C)O(32427)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C(O)C([O])(CC)C[CH]C(6359)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.71,'s^-1'), n=3.021, Ea=(105.562,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 319 used for R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C(O)C(O)([CH][CH]C)CC(6978)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.3913e+06,'s^-1'), n=1.66106, Ea=(123.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C(O)C([O])([CH]C)CCC(6979)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C([O])C(O)(CC)C[CH]C(6981)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2960,'s^-1'), n=2.11, Ea=(84.0984,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;O_H_out] for rate rule [R4H_SS(Cd)S;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C(O)C(O)(CC)C[CH]C(6982)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C(O)(CCC)C(=C)O(32428)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C([O])C(O)([CH]C)CCC(32429)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0508,'s^-1'), n=3.24, Ea=(29.037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_H/NonDeC;XH_out] for rate rule [R5H_CCC(Cd);C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(O)C(O)([CH]C)CCC(32430)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(315594,'s^-1'), n=1.73223, Ea=(38.2227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;Cd_rad_out_singleH;Cs_H_out] + [R5H_RSSR;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5H_DSSS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 3.60555127546
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC(O)(CC)C(=C)O(6986)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C3H6(T)(143)', '[CH2]C(O)=C(O)[CH]C(4609)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.46075e+06,'m^3/(mol*s)'), n=0.027223, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -14.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C[CH]CC1(O)[C](O)CC1C(32431)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.71e+08,'s^-1'), n=0.99, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra_secNd_2H;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C[CH]C1(O)CC(C)C[C]1O(32315)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.51884e+07,'s^-1'), n=0.875, Ea=(71.9648,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_SS_D;doublebond_intra;radadd_intra_csHCs] + [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R5_SS_D;doublebond_intra_secNd_2H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C(O)C(O)(C=CC)CC(6988)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for R3radExo;Y_rad_NDe;XH_Rrad_NDe
Exact match found for rate rule [R3radExo;Y_rad_NDe;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=CC(O)(CCC)C(=C)O(32432)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=CCC(O)(CC)C(=C)O(6989)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.63e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C(O)[C](O)C(C)C[CH]C(32069)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C(O)C(C)[C](O)C[CH]C(32433)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C)C(O)([CH]C)C(=C)O(6787)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C=C(O)C1(O)CC(C)C1C(32073)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    products = ['C[CH]CC(O)([CH]C)C(C)=O(32434)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(204.179,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond_CsC;R_O_H]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction34',
    reactants = ['CHCH3(T)(95)', 'C=C(O)[C](O)C[CH]C(32435)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction35',
    reactants = ['CHCH3(T)(95)', '[CH2]C(O)([CH]C)C(=C)O(32436)'],
    products = ['C=C(O)C(O)([CH]C)C[CH]C(6977)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.06732e+06,'m^3/(mol*s)'), n=0.472793, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -3.5 to 0 kJ/mol."""),
)

network(
    label = '5433',
    isomers = [
        'C=C(O)C(O)([CH]C)C[CH]C(6977)',
    ],
    reactants = [
        ('C3H6(72)', 'C=C(O)C(O)=CC(5562)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = '5433',
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

