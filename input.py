#Data sources
database(
    thermoLibraries = ['BurkeH2O2','FFCM1(-)','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR'],
    reactionLibraries =[('FFCM1(-)',False),('2005_Senosiain_OH_C2H2',False)],
    seedMechanisms = ['BurkeH2O2inN2'],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=10,
    allowSingletO2 = False,
)

# List of species
species(
    label='C7H10',
    reactive=True,
    structure=SMILES("C1C=CCCCC=1"),
)
species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)
species(
    label='CO2',
    reactive=True,
    structure=SMILES("C(=O)=O"),
)
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)
species(
    label='C2H2',
    reactive=True,
    structure=SMILES("C#C"),
)
species(
    label='C4H4',
    reactive=True,
    structure=SMILES("C#CC=C"),
)
species(
    label='C5H6',
    reactive=True,
    structure=SMILES("C1C=CCC=1"),
)
species(
    label='C5H5',
    reactive=True,
    structure=SMILES("[CH]1C=CC=C1"),
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
    structure=SMILES("N#N"),
)

# Reaction systems
simpleReactor(
    temperature=[(700,'K'),(1500,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        "C7H10": 0.010,
        "O2":    0.095, # phi=1 means 9.5 O2 per C7H10
        "N2":    0.095 * 8.1, # 8.1 times as much N2 as O2
        "C2H2": 0,
        "C4H4": 0,
        "C5H6": 0,
        "C5H5": 0,
    },
    terminationRateRatio=0.01,
    terminationTime=(1e0,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0008,
    toleranceMoveToCore=0.008,
    toleranceInterruptSimulation=1E8,
    maximumEdgeSpecies=100000,
    filterReactions=True,
    maxNumObjsPerIter=2,
)

quantumMechanics(
    software='mopac',
    method='pm3',
    # fileStore='QMfiles', # relative to where it is run. Defaults within the output folder.
    scratchDirectory = None, # not currently used
    onlyCyclics = True,
    maxRadicalNumber = 0, 
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,3000,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=True,
    generatePlots=True,
    saveSimulationProfiles=False,
    saveEdgeSpecies=True,
)


