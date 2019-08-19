species(
    label = 'C7H7O(477)(476)',
    structure = SMILES('[O]C1C=CC2CC2C=1'),
    E0 = (119.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4298.03,'J/mol'), sigma=(6.88159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.34 K, Pc=29.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50151,0.0369302,3.97388e-05,-7.87999e-08,3.28524e-11,14535.1,19.764], Tmin=(100,'K'), Tmax=(968.373,'K')), NASAPolynomial(coeffs=[15.5218,0.0234322,-8.15075e-06,1.53248e-09,-1.14138e-13,9737.23,-58.1775], Tmin=(968.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ)"""),
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
    label = 'O=C1C=CC2CC2=C1(4373)',
    structure = SMILES('O=C1C=CC2CC2=C1'),
    E0 = (201.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49891,0.0469449,-1.55302e-05,-6.26513e-09,3.48984e-12,24353.9,21.7667], Tmin=(100,'K'), Tmax=(1279.82,'K')), NASAPolynomial(coeffs=[12.4105,0.0290741,-1.36099e-05,2.64479e-09,-1.86467e-13,20231.6,-38.761], Tmin=(1279.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + polycyclic(s2_3_6_diene_0_3)"""),
)

species(
    label = '[O]C1C=C[C]2CC2C=1(4362)',
    structure = SMILES('[O]C1C=C[C]2CC2C=1'),
    E0 = (216.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38878,0.0481031,-1.4027e-05,-1.06903e-08,5.78564e-12,26183.6,19.6637], Tmin=(100,'K'), Tmax=(1135.44,'K')), NASAPolynomial(coeffs=[11.9411,0.0289374,-1.24988e-05,2.38136e-09,-1.68121e-13,22626.4,-37.6977], Tmin=(1135.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = '[O]C1C=CC2C[C]2C=1(4363)',
    structure = SMILES('[O]C1C=CC2C[C]2C=1'),
    E0 = (216.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3888,0.0481029,-1.40264e-05,-1.06909e-08,5.78588e-12,26183.6,19.6636], Tmin=(100,'K'), Tmax=(1135.43,'K')), NASAPolynomial(coeffs=[11.9411,0.0289375,-1.24989e-05,2.38138e-09,-1.68122e-13,22626.5,-37.6972], Tmin=(1135.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Aromatic_pi_S_1_3) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1C=CC2[CH]C2C=1(4364)',
    structure = SMILES('[O]C1C=CC2[CH]C2C=1'),
    E0 = (345.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66427,0.0375281,2.14983e-05,-5.25522e-08,2.1968e-11,41696.3,21.5279], Tmin=(100,'K'), Tmax=(994.079,'K')), NASAPolynomial(coeffs=[12.9931,0.0250497,-9.62855e-06,1.82495e-09,-1.32711e-13,37808.1,-41.2874], Tmin=(994.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(cyclopropane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1=[C]C2CC2C=C1(4365)',
    structure = SMILES('[O]C1=[C]C2CC2C=C1'),
    E0 = (357.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46675,0.0415493,1.48879e-05,-4.94522e-08,2.18533e-11,43138.5,20.3626], Tmin=(100,'K'), Tmax=(981.37,'K')), NASAPolynomial(coeffs=[14.369,0.0228708,-8.39337e-06,1.57327e-09,-1.1479e-13,38973.2,-49.9604], Tmin=(981.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1C=[C]C2CC2C=1(4366)',
    structure = SMILES('[O]C1C=[C]C2CC2C=1'),
    E0 = (357.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46675,0.0415493,1.48879e-05,-4.94522e-08,2.18533e-11,43138.5,20.3626], Tmin=(100,'K'), Tmax=(981.37,'K')), NASAPolynomial(coeffs=[14.369,0.0228708,-8.39337e-06,1.57327e-09,-1.1479e-13,38973.2,-49.9604], Tmin=(981.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[C]=CC2CC2C=1(4367)',
    structure = SMILES('[O]C1[C]=CC2CC2C=1'),
    E0 = (318.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.122,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41721,0.0435529,1.03431e-05,-4.54917e-08,2.08816e-11,38467.4,19.663], Tmin=(100,'K'), Tmax=(966.975,'K')), NASAPolynomial(coeffs=[14.0027,0.0234822,-8.14739e-06,1.46951e-09,-1.0508e-13,34537.8,-48.3656], Tmin=(966.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=CJC=C) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C1=CCC2CC2=C1(4371)',
    structure = SMILES('[O]C1=CCC2CC2=C1'),
    E0 = (198.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16901,0.044305,2.50237e-05,-6.97978e-08,3.14777e-11,23946.2,23.9812], Tmin=(100,'K'), Tmax=(951.634,'K')), NASAPolynomial(coeffs=[17.7696,0.0191666,-5.71334e-06,1.02646e-09,-7.7377e-14,18765.4,-65.904], Tmin=(951.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_0_2) + radical(C=C(C)OJ)"""),
)

species(
    label = 'OC1=[C]C2CC2C=C1(4356)',
    structure = SMILES('OC1=[C]C2CC2C=C1'),
    E0 = (220.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11149,0.0453961,2.28444e-05,-6.7119e-08,3.01361e-11,26581,20.1364], Tmin=(100,'K'), Tmax=(962.411,'K')), NASAPolynomial(coeffs=[18.0041,0.0196844,-6.43492e-06,1.20431e-09,-9.123e-14,21268.7,-71.4195], Tmin=(962.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Cds_S)"""),
)

species(
    label = 'OC1[C]=CC2CC2C=1(4357)',
    structure = SMILES('OC1[C]=CC2CC2C=1'),
    E0 = (181.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06226,0.0473898,1.83676e-05,-6.33037e-08,2.92557e-11,21909.9,19.4361], Tmin=(100,'K'), Tmax=(951.028,'K')), NASAPolynomial(coeffs=[17.665,0.0202513,-6.16402e-06,1.09479e-09,-8.1051e-14,16821.3,-69.9791], Tmin=(951.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=CJC=C)"""),
)

species(
    label = 'OC1C=CC2C[C]2C=1(4358)',
    structure = SMILES('OC1C=CC2C[C]2C=1'),
    E0 = (79.0521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05204,0.0518219,-6.07857e-06,-2.76341e-08,1.34286e-11,9625.18,19.3654], Tmin=(100,'K'), Tmax=(1036.46,'K')), NASAPolynomial(coeffs=[14.9244,0.026793,-1.11147e-05,2.14356e-09,-1.55163e-13,5218.3,-55.4443], Tmin=(1036.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.0521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'OC1C=[C]C2CC2C=1(4359)',
    structure = SMILES('OC1C=[C]C2CC2C=1'),
    E0 = (220.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11149,0.0453961,2.28444e-05,-6.7119e-08,3.01361e-11,26581,20.1364], Tmin=(100,'K'), Tmax=(962.411,'K')), NASAPolynomial(coeffs=[18.0041,0.0196844,-6.43492e-06,1.20431e-09,-9.123e-14,21268.7,-71.4195], Tmin=(962.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Cds_S)"""),
)

species(
    label = 'OC1C=C[C]2CC2C=1(4360)',
    structure = SMILES('OC1C=C[C]2CC2C=1'),
    E0 = (79.0521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05204,0.0518219,-6.07857e-06,-2.76341e-08,1.34286e-11,9625.18,19.3654], Tmin=(100,'K'), Tmax=(1036.46,'K')), NASAPolynomial(coeffs=[14.9244,0.026793,-1.11147e-05,2.14356e-09,-1.55163e-13,5218.3,-55.4443], Tmin=(1036.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.0521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(Aromatic_pi_S_1_3)"""),
)

species(
    label = 'OC1C=CC2[CH]C2C=1(4361)',
    structure = SMILES('OC1C=CC2[CH]C2C=1'),
    E0 = (208.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,300,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3095,0.0413747,2.94246e-05,-7.0128e-08,3.01857e-11,25138.7,21.2997], Tmin=(100,'K'), Tmax=(971.275,'K')), NASAPolynomial(coeffs=[16.5975,0.0219133,-7.69816e-06,1.46247e-09,-1.09681e-13,20117.2,-62.5727], Tmin=(971.275,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(cyclopropane)"""),
)

species(
    label = 'C1=CC2CC2C2O[C]12(4368)',
    structure = SMILES('C1=CC2CC2C2O[C]12'),
    E0 = (299.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79955,0.034411,3.10313e-05,-5.64688e-08,2.11649e-11,36099.7,18.3436], Tmin=(100,'K'), Tmax=(1051.92,'K')), NASAPolynomial(coeffs=[11.3444,0.032059,-1.40159e-05,2.75522e-09,-2.00476e-13,32213.7,-37.1141], Tmin=(1051.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + polycyclic(s2_3_6_ene_1) - ring(Cyclohexene) + radical(C2CsJO)"""),
)

species(
    label = '[CH]1C2=CC3CC3C1O2(4369)',
    structure = SMILES('[CH]1C2=CC3CC3C1O2'),
    E0 = (398.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22415,0.0301506,9.0736e-05,-1.44044e-07,5.69285e-11,48091,12.7322], Tmin=(100,'K'), Tmax=(985.93,'K')), NASAPolynomial(coeffs=[23.2932,0.0176982,-7.59452e-06,1.74493e-09,-1.46611e-13,39992.8,-112.425], Tmin=(985.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s3_4_6_ene_0) + polycyclic(s2_3_6_ene_1) - ring(3,4-Dihydro-2H-pyran) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]1C2OC2=CC2CC12(4326)',
    structure = SMILES('[CH]1C2OC2=CC2CC12'),
    E0 = (307.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84735,0.0202455,9.5743e-05,-1.46864e-07,5.98309e-11,37051.2,20.9235], Tmin=(100,'K'), Tmax=(940.522,'K')), NASAPolynomial(coeffs=[19.9093,0.0124406,-1.87358e-06,3.45844e-10,-3.65275e-14,30601.3,-81.3345], Tmin=(940.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_3_6_ene_1) + Estimated bicyclic component: polycyclic(s2_3_6_ane) - ring(Cyclohexane) - ring(Ethylene_oxide) + ring(Cyclohexene) + ring(Ethylene_oxide) - ring(Cyclohexene) + radical(CCJCO)"""),
)

species(
    label = 'CH2(S)(21)(22)',
    structure = SMILES('[CH2]'),
    E0 = (418.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1358.21,2621.43,3089.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(418.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]c1ccccc1(4370)',
    structure = SMILES('[O]c1ccccc1'),
    E0 = (39.8872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22108,0.0226145,5.07101e-05,-8.26688e-08,3.32557e-11,4876.28,16.2889], Tmin=(100,'K'), Tmax=(962.489,'K')), NASAPolynomial(coeffs=[13.3453,0.0176701,-5.92751e-06,1.1285e-09,-8.61634e-14,822.531,-46.8882], Tmin=(962.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.8872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CbH) + group(Cb-O2s) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + group(Cb-H) + ring(Benzene) + radical(RC=COJ)"""),
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
    label = '[C]1=CC2CC2C=C1(4372)',
    structure = SMILES('[C]1=CC2CC2C=C1'),
    E0 = (395.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.1305,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12445,0.0244179,5.74202e-05,-8.86402e-08,3.46591e-11,47596.5,16.3811], Tmin=(100,'K'), Tmax=(970.105,'K')), NASAPolynomial(coeffs=[12.0678,0.0254255,-9.08952e-06,1.70148e-09,-1.25083e-13,43690.6,-41.4757], Tmin=(970.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + polycyclic(s2_3_6_diene_1_3) + radical(C=CJC=C)"""),
)

species(
    label = 'O=C1C=CC2C[C]2C1(4374)',
    structure = SMILES('O=C1C=CC2C[C]2C1'),
    E0 = (179.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21303,0.0397065,-4.01653e-06,-9.22743e-09,2.86111e-12,21685.3,19.0272], Tmin=(100,'K'), Tmax=(1570.7,'K')), NASAPolynomial(coeffs=[10.3319,0.0361185,-1.69084e-05,3.17052e-09,-2.14602e-13,17027,-30.5172], Tmin=(1570.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(Tertalkyl)"""),
)

species(
    label = 'O=C1C=C[C]2CC2C1(4375)',
    structure = SMILES('O=C1C=C[C]2CC2C1'),
    E0 = (126.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07228,0.0312266,3.1274e-05,-4.93262e-08,1.69952e-11,15272.7,16.6212], Tmin=(100,'K'), Tmax=(1104.6,'K')), NASAPolynomial(coeffs=[8.12472,0.0380089,-1.69092e-05,3.27591e-09,-2.33717e-13,12184.7,-21.1057], Tmin=(1104.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(Allyl_T)"""),
)

species(
    label = 'O=C1C=CC2[CH]C2C1(4376)',
    structure = SMILES('O=C1C=CC2[CH]C2C1'),
    E0 = (220.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00393,0.0364731,1.11255e-05,-2.62807e-08,8.57985e-12,26569.1,20.2815], Tmin=(100,'K'), Tmax=(1245.95,'K')), NASAPolynomial(coeffs=[8.30319,0.0383678,-1.77832e-05,3.43499e-09,-2.41397e-13,23282.6,-18.3834], Tmin=(1245.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(cyclopropane)"""),
)

species(
    label = 'O=C1[C]=CC2CC2C1(4377)',
    structure = SMILES('O=C1[C]=CC2CC2C1'),
    E0 = (238.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64809,0.0415136,6.45271e-06,-2.76648e-08,1.02391e-11,28755.1,19.9082], Tmin=(100,'K'), Tmax=(1164.89,'K')), NASAPolynomial(coeffs=[10.7693,0.0345614,-1.59733e-05,3.12726e-09,-2.2322e-13,24976.7,-32.5853], Tmin=(1164.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(C=CJC=O)"""),
)

species(
    label = 'O=C1C=[C]C2CC2C1(4378)',
    structure = SMILES('O=C1C=[C]C2CC2C1'),
    E0 = (232.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77344,0.0409159,2.88791e-06,-2.08533e-08,7.38064e-12,28012.7,19.2323], Tmin=(100,'K'), Tmax=(1234.66,'K')), NASAPolynomial(coeffs=[9.50672,0.0364411,-1.66775e-05,3.21122e-09,-2.25629e-13,24534.6,-26.0586], Tmin=(1234.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_3_6_ene_1) + radical(Cds_S)"""),
)

species(
    label = 'O=C1[CH]C2C3CC3C12(4379)',
    structure = SMILES('[O]C1=CC2C3CC3C12'),
    E0 = (161.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65806,0.0336628,4.48555e-05,-8.11329e-08,3.28661e-11,19555.7,17.7466], Tmin=(100,'K'), Tmax=(976.859,'K')), NASAPolynomial(coeffs=[14.55,0.0248277,-9.07105e-06,1.73106e-09,-1.28731e-13,14939.8,-54.8818], Tmin=(976.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + polycyclic(s2_4_4_ene_1) + polycyclic(s2_3_4_ane) - ring(Cyclobutane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=C1C2[CH]C3CC3C12(4380)',
    structure = SMILES('O=C1C2[CH]C3CC3C12'),
    E0 = (264.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84553,0.0357654,2.21674e-05,-4.30527e-08,1.54008e-11,31914.9,19.4453], Tmin=(100,'K'), Tmax=(1111.46,'K')), NASAPolynomial(coeffs=[9.9526,0.03534,-1.60599e-05,3.1497e-09,-2.26314e-13,28336.9,-28.5123], Tmin=(1111.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsCs) + polycyclic(s2_3_5_ane) + polycyclic(s2_3_5_ane) - ring(Cyclopentane) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]=CC1CC1C=C=O(4381)',
    structure = SMILES('[CH]=CC1CC1C=C=O'),
    E0 = (347.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,3120,650,792.5,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96196,0.0521027,-9.16749e-07,-3.9835e-08,1.98663e-11,41864.9,26.8422], Tmin=(100,'K'), Tmax=(975.323,'K')), NASAPolynomial(coeffs=[16.8985,0.0216115,-7.64757e-06,1.42014e-09,-1.03896e-13,37097.8,-58.1459], Tmin=(975.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(347.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_P)"""),
)

species(
    label = 'C7H7O(476)(475)',
    structure = SMILES('[O]C1C=CCC=CC=1'),
    E0 = (89.4153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4422.05,'J/mol'), sigma=(6.98257,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=690.71 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.277261,0.0660949,-2.22778e-05,-3.09166e-08,2.02327e-11,10902.7,21.0868], Tmin=(100,'K'), Tmax=(939.654,'K')), NASAPolynomial(coeffs=[21.5531,0.0149921,-3.7016e-06,6.01876e-10,-4.50583e-14,5161.97,-89.5036], Tmin=(939.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.4153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(1,3,5-Cycloheptatriene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C1C=CC(=O)C=C1(2919)',
    structure = SMILES('[CH2]C1C=CC(=O)C=C1'),
    E0 = (163.305,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51968,0.0527733,-2.97185e-05,8.20686e-09,-9.28549e-13,19731.1,22.5423], Tmin=(100,'K'), Tmax=(1920.95,'K')), NASAPolynomial(coeffs=[12.2463,0.0304374,-1.22773e-05,2.15393e-09,-1.40803e-13,15610.1,-36.2099], Tmin=(1920.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.305,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + ring(25cyclohexadienone) + radical(Isobutyl)"""),
)

species(
    label = 'C7H7O(481)(480)',
    structure = SMILES('O=C1C=CC2C[CH]C12'),
    E0 = (216.387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4297.64,'J/mol'), sigma=(6.68421,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=671.28 K, Pc=32.65 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04598,0.0293115,4.08646e-05,-6.24334e-08,2.21792e-11,26107.3,20.8423], Tmin=(100,'K'), Tmax=(1068.74,'K')), NASAPolynomial(coeffs=[9.73992,0.035063,-1.56965e-05,3.09528e-09,-2.24596e-13,22489.7,-26.0187], Tmin=(1068.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s2_4_5_ene_1) + radical(CCJCC=O)"""),
)

species(
    label = 'O=C1C=CC2[CH]C1C2(4382)',
    structure = SMILES('O=C1C=CC2[CH]C1C2'),
    E0 = (223.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15511,0.0242874,5.82629e-05,-8.27459e-08,2.98877e-11,26968.6,20.0686], Tmin=(100,'K'), Tmax=(1032.24,'K')), NASAPolynomial(coeffs=[10.517,0.0330884,-1.44019e-05,2.85489e-09,-2.09803e-13,23047.1,-31.1713], Tmin=(1032.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + polycyclic(s3_4_6_ene_1) + radical(CCJCC=O)"""),
)

species(
    label = 'O=[C]C1C=CC2CC12(4383)',
    structure = SMILES('O=[C]C1C=CC2CC12'),
    E0 = (190.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.13,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67008,0.03538,3.51862e-05,-6.69177e-08,2.66117e-11,23020.9,22.0522], Tmin=(100,'K'), Tmax=(1002.33,'K')), NASAPolynomial(coeffs=[13.1742,0.0279051,-1.11456e-05,2.1545e-09,-1.58046e-13,18784.1,-43.1062], Tmin=(1002.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + polycyclic(s2_3_5_ene_1) + radical(CC(C)CJ=O)"""),
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
    E0 = (420.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (429.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (429.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (557.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (570.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (569.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (534.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (383.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (370.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (323.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (220.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (264.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (208.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (267.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (333.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (398.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (307.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (458.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (638.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (301.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (253.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (347.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (379.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (276.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (250.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (375.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (363.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (119.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (199.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (376.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (383.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (333.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction108',
    reactants = ['H(3)(3)', 'O=C1C=CC2CC2=C1(4373)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(89.419,'m^3/(mol*s)'), n=1.598, Ea=(7.0333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsCs_Cds-OneDeH;HJ] for rate rule [Cds-CsCs_Cds-COH;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction96',
    reactants = ['H(3)(3)', '[O]C1C=C[C]2CC2C=1(4362)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction97',
    reactants = ['H(3)(3)', '[O]C1C=CC2C[C]2C=1(4363)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.92e+13,'cm^3/(mol*s)'), n=0.18, Ea=(0.518816,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""From training reaction 123 used for H_rad;C_rad/OneDeCs2
Exact match found for rate rule [C_rad/OneDeCs2;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction98',
    reactants = ['H(3)(3)', '[O]C1C=CC2[CH]C2C=1(4364)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+13,'cm^3/(mol*s)','*|/',3.16), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 59 used for H_rad;C_rad/H/NonDeC
Exact match found for rate rule [C_rad/H/NonDeC;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction99',
    reactants = ['H(3)(3)', '[O]C1=[C]C2CC2C=C1(4365)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.78711e+07,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction100',
    reactants = ['H(3)(3)', '[O]C1C=[C]C2CC2C=1(4366)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 40 used for Cd_rad/NonDe;H_rad
Exact match found for rate rule [Cd_rad/NonDe;H_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction101',
    reactants = ['H(3)(3)', '[O]C1[C]=CC2CC2C=1(4367)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.117e+14,'cm^3/(mol*s)'), n=-0.152, Ea=(4.19655,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 49 used for Cd_rad/Cd;H_rad
Exact match found for rate rule [H_rad;Cd_rad/Cd]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction106',
    reactants = ['[O]C1=CCC2CC2=C1(4371)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.66e+08,'s^-1'), n=1.162, Ea=(185.28,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [1_3_unsaturated_pentane_backbone;CH2(C)_1;unsaturated_end] for rate rule [1_3_unsaturated_pentane_backbone;CH2(C)_1;Cd(C)C_2]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction90',
    reactants = ['OC1=[C]C2CC2C=C1(4356)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_Cs;XH_out] for rate rule [R3H_DS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction91',
    reactants = ['OC1[C]=CC2CC2C=1(4357)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction92',
    reactants = ['C7H7O(477)(476)'],
    products = ['OC1C=CC2C[C]2C=1(4358)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+11,'s^-1'), n=0, Ea=(100.834,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H;O_rad_out;Cs_H_out_Cs2] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_Cs2_cy3]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction93',
    reactants = ['OC1C=[C]C2CC2C=1(4359)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction94',
    reactants = ['C7H7O(477)(476)'],
    products = ['OC1C=C[C]2CC2C=1(4360)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(126000,'s^-1'), n=1.85, Ea=(88.2824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SMSS;Y_rad_out;Cs_H_out] for rate rule [R5H_SMSS;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction95',
    reactants = ['OC1C=CC2[CH]C2C=1(4361)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(148.499,'s^-1'), n=2.47528, Ea=(59.2966,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H;C_rad_out_H/NonDeC;XH_out] + [R5H_SSMS;C_rad_out_single;XH_out] for rate rule [R5H_SSMS;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction102',
    reactants = ['C7H7O(477)(476)'],
    products = ['C1=CC2CC2C2O[C]12(4368)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.4423e+12,'s^-1'), n=-0.141781, Ea=(213.662,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_cyclic;doublebond_intra;radadd_intra] for rate rule [Rn1c6_alpha_short;doublebond_intra_secDe_HNd;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction103',
    reactants = ['C7H7O(477)(476)'],
    products = ['[CH]1C2=CC3CC3C1O2(4369)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(278.822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_cyclic;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c6_beta_short;doublebond_intra_pri_HNd_Cs;radadd_intra_O]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 278.6 to 278.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction104',
    reactants = ['C7H7O(477)(476)'],
    products = ['[CH]1C2OC2=CC2CC12(4326)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.48e+11,'s^-1'), n=0.26, Ea=(187.242,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra] for rate rule [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 185.0 to 187.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction105',
    reactants = ['CH2(S)(21)(22)', '[O]c1ccccc1(4370)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.54e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [elec_def;mb_db_twocdisub] for rate rule [carbene;mb_db_twocdisub_twoDe]
Euclidian distance = 1.41421356237
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction107',
    reactants = ['O(4)(4)', '[C]1=CC2CC2C=C1(4372)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2085.55,'m^3/(mol*s)'), n=1.09077, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction109',
    reactants = ['O=C1C=CC2C[C]2C1(4374)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_Cs2;Cs_H_out_H/OneDe] for rate rule [R2H_S;C_rad_out_Cs2_cy3;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction110',
    reactants = ['O=C1C=C[C]2CC2C1(4375)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.28e+07,'s^-1'), n=1.56, Ea=(126.775,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;Cs_H_out_H/OneDe] for rate rule [R3H_SS_12cy3;C_rad_out_OneDe/Cs;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction111',
    reactants = ['O=C1C=CC2[CH]C2C1(4376)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.82e+09,'s^-1'), n=0.73, Ea=(127.612,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeC;Cs_H_out_H/CO]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction112',
    reactants = ['O=C1[C]=CC2CC2C1(4377)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/NonDeC] for rate rule [R3H_SS;Cd_rad_out_Cd;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction113',
    reactants = ['O=C1C=[C]C2CC2C1(4378)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;Cs_H_out_H/OneDe] for rate rule [R4H_SSS;Cd_rad_out_Cd;Cs_H_out_H/CO]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction114',
    reactants = ['C7H7O(477)(476)'],
    products = ['O=C1[CH]C2C3CC3C12(4379)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.4227e+18,'s^-1'), n=-0.859165, Ea=(130.857,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0cx_gamma;doublebond_intra_pri;radadd_intra_csHCs] for rate rule [Rn0c6_gamma;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction115',
    reactants = ['C7H7O(477)(476)'],
    products = ['O=C1C2[CH]C3CC3C12(4380)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.98e+12,'s^-1'), n=0.5, Ea=(255.224,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c6_beta_long_SS_D;doublebond_intra_pri_HDe;radadd_intra_csHDe] for rate rule [Rn0c6_beta_long_SS_D;doublebond_intra_pri_HCO;radadd_intra_csHCO]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction116',
    reactants = ['[CH]=CC1CC1C=C=O(4381)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.03e+10,'s^-1'), n=0.19, Ea=(16.736,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6_DSS_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C7H7O(476)(475)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.1e+12,'s^-1'), n=0.14, Ea=(30.5535,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra] for rate rule [Rn0c7_beta_long_SS_D;doublebond_intra_pri_HNd_Cs;radadd_intra_csHCs]
Euclidian distance = 3.0
family: Intra_R_Add_Endocyclic
Ea raised from 26.3 to 30.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction117',
    reactants = ['[CH2]C1C=CC(=O)C=C1(2919)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.53542e+11,'s^-1'), n=0.21, Ea=(35.9824,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_cyclic;doublebond_intra_pri;radadd_intra_cs2H] + [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_cs] for rate rule [Rn1c6_alpha_long;doublebond_intra_pri;radadd_intra_cs2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction118',
    reactants = ['C7H7O(481)(480)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction119',
    reactants = ['O=C1C=CC2[CH]C1C2(4382)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction120',
    reactants = ['O=[C]C1C=CC2CC12(4383)'],
    products = ['C7H7O(477)(476)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

network(
    label = '362',
    isomers = [
        'C7H7O(477)(476)',
    ],
    reactants = [
    ],
    bathGas = {
        'Ne': 0.333333,
        'N2': 0.333333,
        'Ar(8)': 0.333333,
    },
)

pressureDependence(
    label = '362',
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

