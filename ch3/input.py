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
simpleReactor(
    temperature=(1000,'K'),
    pressure=(20,'bar'),
    initialMoleFractions={
        "CH4": 1,
        "O2": 2,
        "N2": 16,
    },
    terminationConversion={
        'O2': 0.95,
    },
    terminationTime=(1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=0.5,
    maximumEdgeSpecies=1000
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
)