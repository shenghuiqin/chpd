# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)
species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N")
)
# Reaction systems

# Reaction systems
simpleReactor(
    temperature=[(700,'K'),(1000,'K')],
    pressure=[(20.0,'bar'),(40.0,'bar')],
    nSims=8,
    initialMoleFractions={
        #"C7H10": 1,
        "O2":  0.105263158, # phi=1 means 9.5 O2 per C7H10
        "N2":0.842105263 , # 8.1 times as much N2 as O2
        #"C7H16":1,  
        "CH4": 0.052631579,
    },
    terminationTime=(1e0,'s'), 
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.3,
    toleranceInterruptSimulation=0.3,
    maximumEdgeSpecies=1000
)

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=False,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
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


)