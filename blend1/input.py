#Data sources
database(
    thermoLibraries =['BurkeH2O2','Klippenstein_Glarborg2016','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo','primaryThermoLibrary'], # 'FFCM1(-)','primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo','CBS_QB3_1dHR'
    reactionLibraries = [('2005_Senosiain_OH_C2H2',False),('Glarborg/C3', False)], # 
    seedMechanisms = ['BurkeH2O2inN2','Klippenstein_Glarborg2016'], #
    kineticsDepositories = ['training'], 
    kineticsFamilies = ['default'],
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=9,
    maximumOxygenAtoms=4,
    allowSingletO2 = True,
)

# List of species

species(
    label='CHPD',
    reactive=True,
    structure=SMILES("C1C=CCCCC=1"),
)
species(
    label='C7H16',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
        """),
)
species(
    label='O2',
    reactive=True,
        structure=adjacencyList(
        """
        1 O u1 p2 {2,S}
        2 O u1 p2 {1,S}
        """),
)
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N")
)
species(
    label='C7H9',
    reactive=True,
    structure=SMILES("[CH]1CC=CC=CC1"),
)
species(
    label='C7H10',
    reactive=True,
    structure=SMILES("C=CC=CCC=C"),
)
species(
    label='C7H11',
    reactive=True,
    structure=SMILES("[CH2]CCC=CC=C"),
)
species(
    label='C5H7O2',
    reactive=True,
    structure=SMILES("C=CC=CCOO"),
)
species(
    label='C6H12O',
    reactive=True,
    structure=SMILES("CCCCC=CO"),
)
species(
    label='C7H9O2',
    reactive=True,
    structure=SMILES("[O]OC1CC=CC=CC1"),
)
species(
    label='C7H9O',
    reactive=True,
    structure=SMILES("[O]C1CC=CC=CC1"),
)
species(
    label='C7H9O2(2)',
    reactive=True,
    structure=SMILES("OO[C]1CC=CC=CC1"),
)
species(
    label='C7H8',
    reactive=True,
    structure=SMILES("C1=CC=CCC=C1"),
)
species(
    label='C7H10O2',
    reactive=True,
    structure=SMILES("OOC1CC=CC=CC1"),
)
species(
    label='C7H8O(1)',
    reactive=True,
    structure=SMILES("C1=CC2CC(C=C1)O2"),
)
species(
    label='C7H8O(2)',
    reactive=True,
    structure=SMILES("C1=CC2OC2C=CC1"),
)
species(
    label='C7H9O4',
    reactive=True,
    structure=SMILES("[O]OC1C=CCC=CC1OO"),
)
species(
    label='C7H8O3',
    reactive=True,
    structure=SMILES("O=C1C=CCC=CC1OO"),
)
species(
    label='C7H11O4',
    reactive=True,
    structure=SMILES("[O]OC1(COO)C=CCCC1"),
)
species(
    label='C7H13O2',
    reactive=True,
    structure=SMILES("[O]C1CCCCCC1O"),
)

# Reaction systems
simpleReactor(
    temperature=[(600,'K'),(1000,'K')],
    pressure=[(20.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        #"C7H10": 1,
        "O2":   0.1010803, # phi=1 means 9.5 O2 per C7H10
        "N2":   0.36, # 8.1 times as much N2 as O2
        "CHPD":0.00989, # Cycloheptadiene C7H10  
        "C7H16":0.00946, 
        "C7H9":0,
        "C7H10":0,
        "C7H11":0,
        "C5H7O2":0,
        "C6H12O":0,
        "C7H9O2":0,
        "C7H9O":0,
        "C7H9O2(2)":0,
        "C7H8":0,
        "C7H10O2":0,
        "C7H8O(1)":0,
        "C7H8O(2)":0,
        "C7H9O4":0,
        "C7H8O3":0,
        "C7H11O4":0,
        "C7H13O2":0,     
    },
    terminationTime = (10.0, 's'),
    terminationRateRatio = 0.001,
    terminationConversion={
                'C7H16': 0.99,
        },
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0, # No pruning to start
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=1,
    maxNumObjsPerIter=2, 
    terminateAtMaxObjects=True,
    maxNumSpecies=100, # first stage is until core reaches 100 species
    filterReactions=True, # should speed up model generation
)

model(
    toleranceMoveToCore=0.4,
    toleranceInterruptSimulation=1e8,
    toleranceKeepInEdge=0.01, # Pruning enabled for stage 2
    maximumEdgeSpecies=100000,
    minCoreSizeForPrune=100,
    minSpeciesExistIterationsForPrune=2,
    filterReactions=True,
    )

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=False,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
    #generateSeedEachIteration=True,
)


"""
uncertainty(
    localAnalysis=False,
    globalAnalysis=False,
    uncorrelated=True,
    correlated=True,
    localNumber=10,
    globalNumber=5,
    terminationTime=None,
    pceRunTime=1800,
    logx=True
)
""" 

#quantumMechanics(
#    software='mopac',
#    method='pm3',
    # fileStore='QMfiles', # relative to where you run it from. Defaults to inside the output folder if not defined.
#    scratchDirectory = None, # not currently used
#    onlyCyclics = True,
#    maxRadicalNumber = 0,
#    )

#pressureDependence(
#    method='modified strong collision',
#    maximumGrainSize=(0.5,'kcal/mol'),
#    minimumNumberOfGrains=250,
#    temperatures=(300,2000,'K',8),
#    pressures=(0.01,100,'bar',5),
#    interpolation=('Chebyshev', 6, 4),
#)


