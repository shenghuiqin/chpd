#Data sources
database(
    thermoLibraries =['BurkeH2O2','FFCM1(-)','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo','primaryThermoLibrary'], # 'FFCM1(-)','primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo','CBS_QB3_1dHR'
    reactionLibraries = [('2005_Senosiain_OH_C2H2',False),('Glarborg/C3', False)], # 
    seedMechanisms = ['BurkeH2O2inN2','FFCM1(-)'], #
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
    label='H(3)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """),
)


species(
    label='O(4)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u2 p2 c0
        """),
)


species(
    label='OH(5)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """),
)


species(
    label='H2(6)',
    reactive=True,
    structure=adjacencyList(
        """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """),
)


species(
    label='H2O(7)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='Ar(8)',
    reactive=False,
    structure=adjacencyList(
        """
1 Ar u0 p4 c0
        """),
)


species(
    label='He(9)',
    reactive=False,
    structure=adjacencyList(
        """
1 He u0 p1 c0
        """),
)


species(
    label='HO2(10)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='H2O2(11)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CO(12)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
        """),
)


species(
    label='CO2(13)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
        """),
)


species(
    label='HCO(14)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C(T)(15)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u2 p1 c0
        """),
)


species(
    label='CH(16)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u1 p1 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH2(T)(17)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH3(18)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH2O(19)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='HCCO(20)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H(21)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2H2(22)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='H2CC(23)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p1 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH2(S)(24)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH3OH(25)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH3O(26)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH2CO(27)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H3(28)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H4(29)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH4(30)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2H6(31)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H5(32)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH2OH(33)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH3CO(34)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH2CHO(35)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
        """),
)


species(
    label='CH3CHO(36)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='N2',
    reactive=False,
    structure=adjacencyList(
        """
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
        """),
)


species(
    label='Ne',
    reactive=False,
    structure=adjacencyList(
        """
1 Ne u0 p4 c0
        """),
)


species(
    label='C7H16(1)',
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
    label='O2(2)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
        """),
)


species(
    label='C4H9(69)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C5H11(71)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {3,S} {15,S} {16,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C3H7(70)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H15(72)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {3,S} {4,S} {22,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H15(73)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {2,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {3,S} {4,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H15(75)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {4,S} {6,S} {22,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(150)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(152)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {5,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {22,S} {23,S} {24,S}
9  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(151)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C2H5O2(54)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(107)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(106)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(108)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C6H12(114)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,D} {16,S}
6  C u0 p0 c0 {5,D} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C3H6(59)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C4H8(81)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(214)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {6,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(666)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(459)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {15,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {5,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(472)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {12,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1097)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {17,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='CH3OO(41)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C5H10(88)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='HCCOH(37)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {5,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
        """),
)


species(
    label='HOC2H2(38)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u1 p0 c0 {2,D} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
        """),
)


species(
    label='HOCO(39)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CH3OOH(40)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {7,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H5O(42)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='cC2H4O(43)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H15(76)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {5,S} {21,S} {22,S}
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
        """),
)


species(
    label='C2H3O2(44)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(422)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {5,D}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,D} {4,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C2H4(78)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p1 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2(45)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u1 p0 c0 {2,T}
2 C u1 p0 c0 {1,T}
        """),
)


species(
    label='CHO3(486)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(412)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
        """),
)



species(
    label='CHO4(489)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {1,S}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C4H7(639)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u1 p0 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1246)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {8,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C2H6O(47)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {9,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2H4O(111)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1398)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(1399)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u0 p0 c0 {4,S} {5,D} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(1460)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {4,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C2H5O(48)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {8,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(179)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(699)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {11,S}
3  O u0 p2 c0 {1,S} {25,S}
4  O u0 p2 c0 {2,S} {26,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
10 C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
11 C u1 p0 c0 {2,S} {7,S} {8,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {3,S}
26 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(201)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H5O(49)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {8,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C4H6(1293)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,D} {5,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {7,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1896)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {24,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {10,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
8  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
10 C u0 p0 c0 {3,D} {6,S} {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2085)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
9  C u0 p0 c0 {2,D} {5,S} {6,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C2H4O(50)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {7,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u0 p0 c0 {2,D} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C4H6(1947)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2088)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {14,S} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='cC2H3O(51)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
        """),
)


species(
    label='OCHCHO(52)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,S} {5,S}
4 C u0 p0 c0 {2,D} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C2H6O2(53)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1285)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
        """),
)




species(
    label='S(519)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2156)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2836)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {2,D} {5,S} {7,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(452)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u0 p0 c0 {7,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {6,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3108)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C2H5O2(55)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {9,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H4O2(56)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {8,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C4H6(1295)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,D} {3,S} {6,S}
2  C u0 p0 c0 {1,D} {4,S} {5,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {7,S} {8,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C3H6O(57)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,D} {2,S} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C3H5(61)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C2H5CO(58)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4 C u1 p0 c0 {1,D} {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
        """),
)


species(
    label='CHO3(477)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C3H4(60)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C3H5(62)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,D} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C3H5(63)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u1 p0 c0 {2,D} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C3H4O(64)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3 C u0 p0 c0 {2,S} {4,D} {8,S}
4 C u0 p0 c0 {1,D} {3,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1454)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='H3CCCH(65)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C3H4O(66)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {2,D} {6,S} {7,S}
4 C u0 p0 c0 {1,D} {2,S} {8,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3313)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 C u0 p0 c0 {3,D} {8,S} {9,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='H2CCCH(67)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {3,D}
3 C u1 p0 c0 {2,D} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
        """),
)



species(
    label='S(3121)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u0 p2 c0 {2,S} {20,S}
4  O u0 p2 c0 {1,S} {19,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
6  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
9  C u1 p0 c0 {6,S} {7,S} {18,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3901)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {17,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {15,S} {16,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3912)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3296)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {14,S} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(4476)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C6H13(74)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u1 p0 c0 {4,S} {18,S} {19,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(1952)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u1 p0 c0 {3,S} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {3,S} {7,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4094)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u0 p0 c0 {4,T} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(3316)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u1 p0 c0 {4,S} {6,S} {10,S}
6  C u0 p0 c0 {1,S} {2,D} {5,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='CCCCCC(77)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C2H4(79)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4480)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
4 C u1 p0 c0 {2,S} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4646)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {3,D} {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C3H6(82)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 C u1 p0 c0 {1,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C4H8(83)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {10,S}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H8(84)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {4,S} {10,S}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4662)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,D} {5,S}
3 C u0 p0 c0 {2,D} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C6H13(124)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {3,S} {5,S} {19,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(196)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {5,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C4H8(85)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H9(86)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(544)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {3,D} {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(3136)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {2,S} {7,S} {8,S}
4 C u1 p0 c0 {1,D} {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C5H10(89)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {2,S} {13,S}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H10(90)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(3528)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,D} {7,S}
4 C u0 p0 c0 {3,D} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C5H10(91)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {5,S} {13,S}
5  C u1 p0 c0 {4,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(174)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(4636)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {2,S} {5,D} {8,S}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C5H10(92)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H11(93)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H11(94)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {10,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {2,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C3H7(96)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {10,S}
4  H u0 p0 c0 {2,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3085)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {17,S}
3  O u1 p2 c0 {8,S}
4  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {8,D} {16,S}
8  C u0 p0 c0 {3,S} {5,S} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C3H6(97)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,S} {7,S}
3 C u1 p0 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(723)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(5151)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4 C u0 p0 c0 {1,D} {2,S} {3,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(102)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H10(103)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
        """),
)


species(
    label='CCCCC(104)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(2947)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {3,D} {5,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C3H8(105)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(5291)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {11,S}
3  O u0 p2 c0 {1,S} {12,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(5293)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {5,S} {11,S}
3  O u0 p2 c0 {1,S} {12,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
7  C u1 p0 c0 {4,D} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(5328)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {3,D} {4,S} {9,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(5314)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {10,S}
2  O u1 p2 c0 {4,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2245)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {1,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u1 p0 c0 {2,S} {4,S} {8,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='CH2O(110)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1162)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {3,S} {13,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(1108)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {1,S} {16,S}
4  O u0 p2 c0 {2,S} {17,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
8  C u1 p0 c0 {2,S} {6,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3071)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {25,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
9  C u0 p0 c0 {5,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
11 C u1 p0 c0 {9,S} {23,S} {24,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1934)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(4489)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {15,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
7  C u0 p0 c0 {3,D} {5,S} {14,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4619)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {2,D} {4,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(3254)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
7  C u0 p0 c0 {2,D} {3,S} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H14(112)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {2,S} {6,D} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(218)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {25,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1525)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {8,S}
2  O u0 p2 c0 {1,S} {25,S}
3  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1313)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {3,S} {18,S} {19,S}
8  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(6658)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {19,S} {20,S}
8  C u0 p0 c0 {7,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
11 C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H14(113)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {5,S} {6,D} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1316)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {1,S} {6,S} {22,S} {23,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(204)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(5632)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {4,D} {7,S}
6  C u0 p0 c0 {3,D} {4,S} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(6402)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {7,S} {15,S}
3  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {2,S} {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(6835)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {9,S}
3  O u0 p2 c0 {2,S} {26,S}
4  O u0 p2 c0 {1,S} {25,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
9  C u0 p0 c0 {2,S} {11,S} {19,S} {20,S}
10 C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
11 C u1 p0 c0 {8,S} {9,S} {24,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(6834)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {9,S}
3  O u0 p2 c0 {2,S} {26,S}
4  O u0 p2 c0 {1,S} {25,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
9  C u0 p0 c0 {2,S} {7,S} {19,S} {20,S}
10 C u0 p0 c0 {11,S} {21,S} {22,S} {23,S}
11 C u1 p0 c0 {8,S} {10,S} {24,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(7425)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {15,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {1,S} {5,D} {7,S}
7  C u1 p0 c0 {6,S} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H14(115)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
2  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {1,S} {3,S} {20,S}
7  C u1 p0 c0 {2,S} {3,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(5641)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {3,D} {5,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(5640)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,D} {4,S} {8,S}
6  C u0 p0 c0 {3,D} {4,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C6H12(116)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {2,S} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(7780)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {24,S}
4  C u0 p0 c0 {2,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {7,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {5,S} {19,S} {20,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(7902)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {17,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {3,D} {5,S} {6,S}
9  C u1 p0 c0 {5,S} {15,S} {16,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='O2(1283)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,D}
2 O u0 p2 c0 {1,D}
        """),
)


species(
    label='S(7823)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {15,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {2,D} {3,S} {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(8143)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {4,S} {6,S} {11,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {3,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H14(117)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {1,S} {3,S} {20,S}
7  C u1 p0 c0 {2,S} {5,S} {21,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(5303)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,D} {3,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(4606)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H14(118)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  C u1 p0 c0 {3,S} {7,S} {20,S}
7  C u1 p0 c0 {2,S} {6,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(2492)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {10,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {2,S} {8,S} {21,S} {22,S}
11 C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C6H12(119)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {3,S} {6,S} {16,S}
6  C u1 p0 c0 {5,S} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H14(120)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {3,S} {7,S} {20,S}
7  C u1 p0 c0 {5,S} {6,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(121)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {2,S} {3,S} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(122)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {3,S} {4,S} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C6H13(123)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u1 p0 c0 {2,S} {3,S} {19,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(8464)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {15,S}
2  O u0 p2 c0 {4,S} {14,S}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {5,D} {6,S}
5  C u0 p0 c0 {3,S} {4,D} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(9208)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,S} {15,S}
2  O u0 p2 c0 {6,S} {14,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {1,S} {4,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H15(125)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {1,S} {21,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(173)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(7008)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {6,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u1 p0 c0 {2,S} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {4,D} {8,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(217)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {25,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H15(126)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {15,S} {16,S} {17,S}
7  C u1 p0 c0 {1,S} {21,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(8466)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {15,S}
2  O u0 p2 c0 {6,S} {16,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {7,S}
7  C u1 p0 c0 {6,S} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C6H12(128)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {3,S} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H14(129)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {2,S} {3,S} {20,S}
7  C u1 p0 c0 {3,S} {5,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(130)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {2,S} {3,S} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(132)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {7,D} {19,S}
7  C u0 p0 c0 {6,D} {20,S} {21,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C6H12(133)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H14(134)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {2,S} {4,S} {20,S}
7  C u1 p0 c0 {3,S} {5,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(2176)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {23,S}
2  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {2,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
8  C u1 p0 c0 {4,S} {5,S} {22,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(10607)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {25,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {7,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(9707)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {10,S}
3  O u0 p2 c0 {2,S} {26,S}
4  O u0 p2 c0 {1,S} {27,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {11,S} {19,S} {20,S}
10 C u0 p0 c0 {2,S} {8,S} {21,S} {22,S}
11 C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
27 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(10779)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,S} {24,S}
3  O u0 p2 c0 {1,S} {25,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {6,S} {18,S} {19,S}
9  C u0 p0 c0 {10,S} {20,S} {21,S} {22,S}
10 C u1 p0 c0 {7,S} {9,S} {23,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {2,S}
25 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4093)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {5,D} {6,S}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C7H14(135)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {3,S} {5,S} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(136)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {4,S} {7,S} {19,S}
7  C u1 p0 c0 {6,S} {20,S} {21,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(10905)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(8813)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {25,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {3,S} {7,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(139)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u0 p0 c0 {1,D} {2,S} {25,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(140)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u0 p0 c0 {1,D} {2,S} {25,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(10933)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,S} {24,S}
3  O u0 p2 c0 {1,S} {25,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {5,S} {18,S} {19,S}
9  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
10 C u1 p0 c0 {4,S} {7,S} {23,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {2,S}
25 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(141)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {2,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u0 p0 c0 {1,D} {2,S} {25,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(142)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u0 p0 c0 {1,D} {7,S} {25,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(143)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
9  C u0 p0 c0 {1,D} {6,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(10782)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {10,S} {24,S}
3  O u0 p2 c0 {1,S} {25,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u1 p0 c0 {2,S} {8,S} {23,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {2,S}
25 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(10958)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,S} {23,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(10965)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {24,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {3,D} {8,S} {23,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CO4(145)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {1,S}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
        """),
)


species(
    label='CH2O3(147)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2502)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {18,S}
7  C u1 p0 c0 {6,S} {19,S} {20,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(11215)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {16,S}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {9,D} {20,S}
9  C u0 p0 c0 {8,D} {21,S} {22,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(11214)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {9,D} {21,S}
9  C u0 p0 c0 {7,S} {8,D} {22,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3105)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {27,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {3,S}
27 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(11112)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {5,D} {15,S}
5  C u0 p0 c0 {4,D} {6,S} {17,S}
6  C u0 p0 c0 {5,S} {7,D} {16,S}
7  C u0 p0 c0 {6,D} {18,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C4H8(149)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C2H2O(153)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u1 p0 c0 {2,D} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C2H2O(154)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 C u0 p1 c0 {2,S} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(11094)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {2,D} {7,S} {22,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3165)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,D}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,D} {3,S} {12,S}
6  C u0 p0 c0 {2,D} {4,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C2H2(155)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 C u1 p0 c0 {2,D} {3,S}
2 C u1 p0 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(11713)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {1,S} {7,D} {19,S}
7  C u0 p0 c0 {4,S} {6,D} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(3075)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {25,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {6,S} {8,S} {16,S}
8  C u0 p0 c0 {7,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1938)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {15,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(725)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(10423)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {4,S} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C2H2O(158)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C2H3O(159)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {6,S}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(160)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {3,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(161)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {7,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {1,S} {5,S} {17,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(162)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {9,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
6  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
9  C u1 p0 c0 {1,S} {5,S} {6,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(163)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(164)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {6,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {3,S} {6,S} {23,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(165)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(5223)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {1,D} {4,S} {16,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(10299)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {1,D} {4,S} {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(166)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(167)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
8  C u1 p0 c0 {5,S} {19,S} {20,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(168)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {3,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(169)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
8  C u1 p0 c0 {6,S} {19,S} {20,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(170)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {6,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(171)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(172)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(11870)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
7  C u1 p0 c0 {4,S} {8,S} {18,S}
8  C u0 p0 c0 {7,S} {9,D} {19,S}
9  C u0 p0 c0 {8,D} {20,S} {21,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(175)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {1,S} {5,S} {6,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(5149)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {6,D} {9,S}
6  C u0 p0 c0 {5,D} {10,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(176)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {3,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1407)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
4  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {3,S} {18,S}
6  C u0 p0 c0 {2,S} {7,D} {20,S}
7  C u0 p0 c0 {4,S} {6,D} {19,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(177)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {3,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(178)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(180)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {6,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(5273)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(12790)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {3,S} {4,S} {6,D}
6 C u0 p0 c0 {2,S} {5,D} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(12838)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(13006)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {6,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(12842)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {4,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(5636)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u0 p0 c0 {2,D} {4,S} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(5624)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,S} {9,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {4,D} {7,S}
6 C u0 p0 c0 {3,D} {4,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(181)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1409)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {2,S} {6,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(183)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(184)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {8,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u1 p0 c0 {1,S} {6,S} {20,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(185)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {9,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {1,S} {6,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(186)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(187)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {4,S} {6,S} {23,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(188)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(13005)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
7  C u0 p0 c0 {9,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {9,D} {22,S}
9  C u0 p0 c0 {7,S} {8,D} {21,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(3074)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(189)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(190)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(191)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13690)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {11,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {25,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
8  C u0 p0 c0 {10,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
10 C u0 p0 c0 {8,S} {19,S} {20,S} {21,S}
11 C u1 p0 c0 {2,S} {7,S} {8,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(192)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
8  C u1 p0 c0 {6,S} {19,S} {20,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(193)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {6,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(194)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {6,S} {22,S} {23,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(4711)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {1,D} {5,S} {7,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(14115)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {24,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {10,S} {14,S} {15,S}
7  C u0 p0 c0 {9,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {3,D} {6,S} {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(195)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(197)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {1,S} {6,S} {8,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(198)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
9  C u1 p0 c0 {4,S} {6,S} {23,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(199)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
4  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {4,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {6,S} {22,S} {23,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(200)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
9  C u1 p0 c0 {5,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(202)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {6,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(701)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(203)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(205)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {1,S} {4,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(14326)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {9,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {8,S} {14,S} {15,S}
10 C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(206)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {9,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {1,S} {5,S} {6,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(207)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(208)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {4,S} {6,S} {23,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(209)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
8  C u1 p0 c0 {6,S} {19,S} {20,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(210)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u1 p0 c0 {6,S} {8,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(211)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(212)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {1,S} {5,S} {6,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(213)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
9  C u1 p0 c0 {3,S} {6,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(14146)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,S} {22,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {7,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {1,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {21,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(215)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {24,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {7,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(14155)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {2,D} {5,S} {6,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(14411)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {1,D} {4,S} {5,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(216)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(219)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {25,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH2O2(220)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
        """),
)


species(
    label='CHO2(221)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
        """),
)


species(
    label='CH3O3(223)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH3O3(224)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2711)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(3082)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='CH2O3(228)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {3,D} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C3H6O(229)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,D} {2,S} {3,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(230)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,D} {3,S} {9,S}
5 C u0 p0 c0 {2,D} {3,S} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(5226)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
7  C u0 p0 c0 {1,D} {5,S} {19,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(14750)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {27,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {9,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {8,S} {16,S} {17,S}
10 C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
27 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(231)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {4,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(702)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {25,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
10 C u0 p0 c0 {11,S} {21,S} {22,S} {23,S}
11 C u1 p0 c0 {8,S} {10,S} {24,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(232)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {4,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(5182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {14,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(13729)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,D} {5,S} {6,S}
5 C u0 p0 c0 {2,D} {4,S} {7,S}
6 C u0 p0 c0 {3,D} {4,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(14107)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(233)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {2,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13826)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,D}
5 C u0 p0 c0 {2,D} {4,S} {7,S}
6 C u0 p0 c0 {3,D} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(609)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {2,D} {3,S} {4,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(4634)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,D} {6,S}
4 C u0 p0 c0 {3,D} {5,D}
5 C u1 p0 c0 {4,D} {7,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C3H2O(503)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {5,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {6,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(16528)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {2,D} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(234)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
8  C u1 p0 c0 {2,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(235)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {7,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {4,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {6,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(236)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {4,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {6,S} {24,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13726)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {3,D} {4,S} {9,S}
6  C u1 p0 c0 {2,S} {10,S} {11,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(237)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {9,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {4,S} {5,S} {24,S}
9  C u0 p0 c0 {1,D} {3,S} {7,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(726)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(17462)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,D} {4,S} {6,S}
6 C u0 p0 c0 {3,D} {4,S} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(5335)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u1 p0 c0 {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(238)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {4,S} {5,S} {24,S}
9  C u0 p0 c0 {1,D} {5,S} {7,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(12794)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(2724)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(239)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {1,D} {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(240)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {4,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13166)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {5,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(13880)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
8  C u0 p0 c0 {9,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {5,S} {8,S} {22,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(696)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {11,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {25,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
9  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
10 C u0 p0 c0 {11,S} {22,S} {23,S} {24,S}
11 C u1 p0 c0 {2,S} {8,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(241)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {2,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(242)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {2,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {4,S} {5,S} {23,S}
9  C u0 p0 c0 {1,D} {6,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(243)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
4  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {3,S} {5,S} {24,S}
9  C u0 p0 c0 {1,D} {4,S} {7,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(244)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {1,D} {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(17880)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {5,S} {8,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u1 p2 c0 {2,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(18819)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {24,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {10,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {3,D} {7,S} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(245)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {7,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(246)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {7,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(247)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {7,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(248)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {2,S} {7,S} {23,S}
9  C u0 p0 c0 {1,D} {2,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(249)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {7,S} {23,S}
9  C u0 p0 c0 {1,D} {6,S} {24,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(250)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {8,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {6,S} {23,S}
9  C u0 p0 c0 {1,D} {6,S} {24,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(251)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {9,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {5,S} {7,S} {24,S}
9  C u0 p0 c0 {1,D} {4,S} {6,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(252)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,D}
2  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {4,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {1,D} {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(253)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {10,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(254)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(255)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(256)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(257)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {6,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(258)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {4,S} {10,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {11,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {9,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(259)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {11,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {9,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(260)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {11,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,D} {8,S} {10,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(261)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {4,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {11,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,D} {8,S} {10,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(262)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {11,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {10,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {7,S} {24,S} {25,S} {26,S}
11 C u1 p0 c0 {2,S} {3,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(18882)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {23,S}
3  O u0 p2 c0 {1,S} {24,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
8  C u0 p0 c0 {9,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {2,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {22,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {2,S}
24 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(263)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {7,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {6,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(264)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(265)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {4,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(266)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {4,S} {5,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(7422)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {3,S} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(7807)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {12,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {3,S} {4,D}
6  C u1 p0 c0 {4,S} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(18862)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
8  C u0 p0 c0 {9,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {2,D} {6,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(267)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {7,S} {16,S} {17,S}
9  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(268)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {4,S} {11,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {9,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(269)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {11,S} {21,S} {22,S}
10 C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {9,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(270)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {11,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,D} {4,S} {10,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(271)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {11,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,D} {8,S} {10,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(272)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {11,D}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {4,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u1 p0 c0 {2,S} {3,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(8159)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {14,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u1 p0 c0 {5,S} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(273)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(274)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(275)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {4,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(276)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,S} {21,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {2,D} {9,S} {26,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(277)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {11,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {11,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,D} {8,S} {10,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C3H5O(560)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {8,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(278)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {11,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
11 C u1 p0 c0 {2,S} {3,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2720)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u0 p0 c0 {1,S} {2,D} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(8460)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {12,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {4,D} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(738)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(279)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,D} {6,S}
4 C u0 p0 c0 {3,D} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,D} {9,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(280)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 C u0 p0 c0 {3,D} {7,S} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(281)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {8,S} {9,S}
5 C u0 p0 c0 {1,S} {2,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(282)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u1 p0 c0 {3,S} {8,S} {9,S}
5 C u1 p0 c0 {1,S} {2,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(283)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,D} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(284)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(285)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {26,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(286)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {26,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(4629)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {7,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 C u0 p0 c0 {2,D} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(287)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {26,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(14151)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {9,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(21667)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u1 p2 c0 {7,S}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {4,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {4,S} {9,S} {11,S}
8  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(6972)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
6  C u1 p0 c0 {4,S} {15,S} {16,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(288)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {1,S} {6,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {26,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(289)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {10,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {1,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {1,S} {2,D} {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(290)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {10,S} {26,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(291)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {10,S} {26,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(292)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {10,S} {26,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(293)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {10,S} {26,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {10,S} {21,S} {22,S}
9  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(294)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {5,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {10,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {1,S} {2,D} {9,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(295)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
10 C u1 p0 c0 {1,S} {2,D}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(18850)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(15977)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {5,S} {8,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(9205)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {15,S}
2  O u0 p2 c0 {5,S} {14,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(296)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(297)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
4  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(298)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(299)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {15,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {3,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(300)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {8,S} {9,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {4,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(301)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {8,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {4,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21970)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {15,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {3,S} {7,D} {12,S}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(302)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(303)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {10,S} {13,S} {14,S}
6  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {1,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {4,S} {6,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {5,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(304)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {9,S} {10,S} {17,S} {18,S}
7  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {1,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {5,S} {6,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {6,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(305)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(306)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(307)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {3,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(21969)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {14,S}
2  O u0 p2 c0 {6,S} {15,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {13,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(308)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
9  C u1 p0 c0 {3,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(15968)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {3,D} {4,S} {6,S}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(20975)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {11,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {9,S} {10,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(309)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {10,S} {19,S} {20,S}
7  C u0 p0 c0 {8,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {7,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {6,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(310)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {10,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {5,S} {6,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(311)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
5  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {10,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {4,S} {5,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {8,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(312)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {10,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {5,S} {6,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {8,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(313)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u1 p0 c0 {1,S} {2,D}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(314)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(12285)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S} {6,S} {7,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(315)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(3473)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {15,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,D} {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(316)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(317)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(318)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {17,S}
4  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {4,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {3,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(319)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {5,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(320)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {6,S} {7,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {25,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(321)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {10,S} {15,S} {16,S}
7  C u0 p0 c0 {9,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {1,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {5,S} {7,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {6,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(322)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(23622)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(323)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(324)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(325)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {3,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(326)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {10,S} {19,S} {20,S}
8  C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {8,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(327)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {25,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {9,S} {10,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {7,S} {24,S}
10 C u0 p0 c0 {1,S} {2,D} {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(328)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {10,D}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {4,S} {17,S} {18,S}
7  C u0 p0 c0 {9,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {10,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {5,S} {7,S} {25,S}
10 C u0 p0 c0 {1,S} {2,D} {8,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(330)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {12,S}
3  O u1 p2 c0 {12,S}
4  O u0 p2 c0 {12,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {11,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {9,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {8,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {2,S} {3,S} {4,D}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(331)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {10,S} {13,S} {14,S}
10 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(332)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {10,S} {15,S} {16,S}
10 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(333)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
10 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(334)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {19,S} {20,S}
8  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(335)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {5,S} {11,S} {19,S} {20,S}
10 C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(336)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {11,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {1,S} {9,S} {22,S} {23,S}
11 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(13887)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
8  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {5,S} {6,S} {22,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(337)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {10,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {1,S} {8,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(8779)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {13,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {5,D} {11,S}
5  C u0 p0 c0 {1,S} {4,D} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(338)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {11,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {12,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {1,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {9,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(339)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {11,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {12,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {1,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {9,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(340)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {10,S} {13,S} {14,S}
10 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3146)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {9,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(341)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {10,S} {15,S} {16,S}
10 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(342)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(343)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {19,S} {20,S}
8  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(344)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {10,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {7,S} {17,S} {18,S}
9  C u0 p0 c0 {6,S} {11,S} {19,S} {20,S}
10 C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(6996)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {5,S} {9,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u1 p0 c0 {1,S} {4,S} {7,S}
6 C u0 p0 c0 {3,D} {4,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(575)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5 C u0 p0 c0 {2,D} {3,S} {4,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(18983)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {4,S} {9,S}
3 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4 C u0 p0 c0 {2,S} {3,S} {5,D}
5 C u0 p0 c0 {1,S} {4,D} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4624)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {9,S}
2 O u1 p2 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,D}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(26109)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {9,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4 C u0 p0 c0 {1,S} {3,S} {5,D}
5 C u0 p0 c0 {2,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(345)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {11,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {12,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {10,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(14777)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {1,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(26141)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4 C u0 p0 c0 {1,D} {3,S} {5,S}
5 C u0 p0 c0 {2,D} {4,S} {9,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(22106)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {1,S} {4,D} {9,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(346)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {10,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {8,S} {12,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {10,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4621)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {10,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {9,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(347)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {1,S} {7,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {12,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {11,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(348)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {5,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {12,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {11,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(14409)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {25,S}
3  O u1 p2 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {18,S}
9  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4490)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(350)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {12,S}
3  O u1 p2 c0 {12,S}
4  O u0 p2 c0 {12,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {5,S} {25,S} {26,S} {27,S}
11 C u0 p0 c0 {9,S} {22,S} {23,S} {24,S}
12 C u0 p0 c0 {2,S} {3,S} {4,D}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(26265)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {12,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(351)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {7,S} {17,S} {18,S}
9  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(352)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {19,S} {20,S}
8  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(353)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
9  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(354)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(355)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {8,S} {17,S} {18,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(356)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {1,S} {5,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(357)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {11,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {1,S} {9,S} {22,S} {23,S}
11 C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(3123)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(358)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {11,S} {12,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {7,S} {8,S} {20,S} {21,S}
6  C u0 p0 c0 {7,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
9  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {9,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {1,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(359)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {11,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {8,S} {12,S} {20,S} {21,S}
10 C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {1,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {9,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(26264)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {12,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {3,D} {5,S} {11,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(27259)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,S} {11,S}
3  O u0 p2 c0 {1,S} {12,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(360)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {7,S} {17,S} {18,S}
9  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(361)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {19,S} {20,S}
8  C u0 p0 c0 {7,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(362)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {7,S} {19,S} {20,S}
9  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(363)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(13884)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {22,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {8,S} {16,S} {17,S} {18,S}
7  C u1 p0 c0 {4,S} {9,S} {19,S}
8  C u0 p0 c0 {6,S} {9,D} {20,S}
9  C u0 p0 c0 {7,S} {8,D} {21,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(364)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {10,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {11,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {9,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {8,S} {17,S} {18,S}
10 C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(365)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {11,S} {14,S} {15,S}
10 C u0 p0 c0 {5,S} {12,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {10,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(366)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {7,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {9,S} {12,S} {22,S} {23,S}
11 C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {10,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(367)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {12,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {12,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(368)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
10 C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {12,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {11,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(28250)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {8,D} {20,S}
8  C u0 p0 c0 {4,S} {7,D} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(370)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {12,S}
3  O u1 p2 c0 {12,S}
4  O u0 p2 c0 {12,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {9,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {2,S} {3,S} {4,D}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(371)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(28225)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(28359)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,D} {13,S}
5  C u0 p0 c0 {2,S} {4,D} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(372)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(373)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {5,S} {9,S} {19,S} {20,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(28453)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {14,S} {15,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(374)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {10,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {1,S} {8,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {1,S} {3,D} {27,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(375)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {11,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {12,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {1,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {9,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(376)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(377)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {14,S}
7  C u0 p0 c0 {6,S} {9,S} {19,S} {20,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(378)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {10,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {9,S} {19,S} {20,S}
9  C u0 p0 c0 {8,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(379)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {12,S} {27,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {8,S} {12,S} {22,S} {23,S}
11 C u0 p0 c0 {9,S} {24,S} {25,S} {26,S}
12 C u0 p0 c0 {2,S} {3,D} {10,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1322)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {5,D}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(11686)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {6,S}
3  C u0 p0 c0 {1,D} {4,S} {8,S}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  C u0 p0 c0 {2,D} {9,S} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(29601)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {7,D} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {13,S} {14,S}
7  C u0 p0 c0 {4,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(29494)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,D} {6,S}
2  C u0 p0 c0 {1,S} {5,D} {7,S}
3  C u0 p0 c0 {1,D} {8,S} {9,S}
4  C u0 p0 c0 {5,D} {10,S} {11,S}
5  C u0 p0 c0 {2,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(380)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {12,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {12,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {7,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
11 C u0 p0 c0 {12,S} {25,S} {26,S} {27,S}
12 C u0 p0 c0 {1,S} {3,D} {11,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(11918)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  C u1 p0 c0 {1,S} {11,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(30810)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,D}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {7,S} {8,S}
6  C u0 p0 c0 {3,D} {12,S} {13,S}
7  C u1 p0 c0 {5,S} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(381)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {2,S} {26,S} {27,S} {28,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(382)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {20,S} {21,S} {22,S}
10 C u0 p0 c0 {2,S} {26,S} {27,S} {28,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(383)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(30684)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,D} {4,S} {6,S}
2  C u0 p0 c0 {1,D} {7,S} {8,S}
3  C u1 p0 c0 {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {5,T}
5  C u0 p0 c0 {3,S} {4,T}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(31796)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {7,D}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u0 p0 c0 {4,S} {5,D} {11,S}
7  C u0 p0 c0 {4,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(32045)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(384)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(385)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u1 p0 c0 {1,S} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(386)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(387)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(31941)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,D}
4  C u0 p0 c0 {3,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u0 p0 c0 {3,D} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C3H6O(388)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C3H6O(389)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3  C u1 p0 c0 {2,S} {7,S} {8,S}
4  C u1 p0 c0 {1,S} {9,S} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C3H6O(390)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(32051)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {2,D} {4,S} {7,S}
6  C u0 p0 c0 {3,S} {7,D} {12,S}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(33188)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {4,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {7,D}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u0 p0 c0 {4,S} {5,D} {11,S}
7  C u0 p0 c0 {1,S} {4,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(6999)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,S} {8,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {4,D} {7,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {1,S}
9 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(33340)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {13,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {6,D}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,D} {12,S}
7  C u1 p0 c0 {5,S} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(34412)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {13,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {9,S}
6  C u0 p0 c0 {1,S} {4,D} {10,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(33346)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {4,S} {5,D} {10,S}
7  C u0 p0 c0 {2,D} {3,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(391)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(33345)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {13,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {3,S} {6,D} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(34446)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {3,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {7,D} {8,S}
6  C u0 p0 c0 {2,D} {3,S} {12,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(392)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {3,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(393)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u1 p2 c0 {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u1 p0 c0 {1,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(394)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u1 p0 c0 {1,S} {3,S} {9,S}
5  C u1 p0 c0 {2,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(34617)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,D} {3,S} {6,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {2,D} {4,S} {13,S}
7  C u0 p0 c0 {5,D} {11,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(34619)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {6,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {6,D}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,D} {10,S}
7  C u1 p0 c0 {5,S} {11,S} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(395)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(34629)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {15,S}
3  O u0 p2 c0 {1,S} {16,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {4,S} {6,D} {14,S}
9  C u0 p0 c0 {7,D} {12,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(34636)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {15,S}
3  O u0 p2 c0 {6,S} {16,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,D} {11,S}
9  C u0 p0 c0 {8,D} {13,S} {14,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(396)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {3,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(397)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {1,S} {25,S} {26,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(36707)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {15,S}
3  O u0 p2 c0 {8,S} {16,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {3,S} {6,D} {12,S}
9  C u0 p0 c0 {7,D} {13,S} {14,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(398)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {7,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {9,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {25,S} {26,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(36721)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {14,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u0 p0 c0 {3,D} {4,S} {11,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(399)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {1,S} {25,S} {26,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(400)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {2,S} {25,S} {26,S}
8  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(31940)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {5,D} {7,S} {10,S}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {7,T}
7  C u0 p0 c0 {4,S} {6,T}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(401)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {4,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
9  C u1 p0 c0 {1,S} {25,S} {26,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(402)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {4,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
7  C u0 p0 c0 {2,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {1,S} {2,S} {25,S} {26,S}
9  C u0 p0 c0 {6,S} {19,S} {20,S} {21,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(404)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u1 p2 c0 {11,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {4,S} {10,S} {19,S} {20,S}
8  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {7,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,S} {3,S} {27,S} {28,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {11,S}
28 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(406)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u1 p2 c0 {11,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {6,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {4,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {8,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {2,S} {3,S} {27,S} {28,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {11,S}
28 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(10412)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {5,S} {12,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {6,S}
6  C u1 p0 c0 {5,S} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(408)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u1 p2 c0 {11,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {8,S} {24,S} {25,S} {26,S}
11 C u0 p0 c0 {2,S} {3,S} {27,S} {28,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {11,S}
28 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(37340)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {4,S} {14,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u0 p0 c0 {6,S} {8,D} {9,S}
6  C u1 p0 c0 {1,S} {5,S} {10,S}
7  C u0 p0 c0 {3,S} {4,D} {11,S}
8  C u0 p0 c0 {5,D} {12,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(409)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u1 p0 c0 {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(411)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H8O(413)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u1 p0 c0 {3,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H8O(414)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {3,S} {11,S}
5  C u1 p0 c0 {2,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H8O(415)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {6,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C9H19(416)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {6,S} {27,S} {28,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(33317)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
4  C u0 p0 c0 {3,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {9,S} {10,S}
6  C u0 p0 c0 {7,D} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C9H19(417)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {15,S} {16,S}
3  C u0 p0 c0 {1,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {6,S} {27,S} {28,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H19(418)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {17,S} {18,S}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {1,S} {24,S} {25,S} {26,S}
8  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
9  C u1 p0 c0 {6,S} {27,S} {28,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(419)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {2,S} {11,S} {27,S} {28,S}
11 C u1 p0 c0 {10,S} {29,S} {30,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {11,S}
30 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(420)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
5  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {3,S} {24,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {2,S} {11,S} {27,S} {28,S}
11 C u1 p0 c0 {10,S} {29,S} {30,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {11,S}
30 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(37348)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {6,S} {14,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {8,D} {10,S}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {11,S}
8  C u0 p0 c0 {5,D} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(421)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {7,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {2,S} {11,S} {27,S} {28,S}
11 C u1 p0 c0 {10,S} {29,S} {30,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {11,S}
30 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(33324)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {5,S} {7,D}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  C u1 p0 c0 {3,S} {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(31942)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {7,D} {11,S} {12,S}
6  C u0 p0 c0 {4,D} {7,D}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(425)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(426)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
9  C u0 p0 c0 {10,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {9,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(37858)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,D} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {4,D} {10,S}
6  C u1 p0 c0 {1,S} {7,S} {11,S}
7  C u0 p0 c0 {2,S} {3,D} {6,S}
8  C u1 p0 c0 {4,S} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(5436)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {7,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(37864)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {2,S} {3,D} {5,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(427)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {6,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {1,S} {2,S} {10,S} {20,S}
7  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {5,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {6,S} {27,S} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(428)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
9  C u0 p0 c0 {10,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {9,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(429)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
3  C u0 p0 c0 {2,S} {6,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {7,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {2,S} {10,S} {20,S}
6  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
9  C u0 p0 c0 {7,S} {24,S} {25,S} {26,S}
10 C u0 p0 c0 {5,S} {27,S} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(38356)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
5  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {8,D} {12,S}
7  C u0 p0 c0 {2,S} {3,D} {4,S}
8  C u0 p0 c0 {1,S} {6,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(430)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {18,S} {19,S}
4  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {4,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
9  C u0 p0 c0 {10,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {9,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(431)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {6,S}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {18,S} {19,S}
4  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {4,S} {16,S} {17,S}
6  C u0 p0 c0 {1,S} {2,S} {10,S} {20,S}
7  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {21,S} {22,S} {23,S}
10 C u0 p0 c0 {6,S} {27,S} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(433)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {9,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {18,S} {19,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {11,S} {20,S} {21,S}
8  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {3,S} {12,S} {22,S}
10 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {7,S} {26,S} {27,S} {28,S}
12 C u0 p0 c0 {9,S} {29,S} {30,S} {31,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
28 H u0 p0 c0 {11,S}
29 H u0 p0 c0 {12,S}
30 H u0 p0 c0 {12,S}
31 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(33325)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,D} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {6,D}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(435)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {9,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
5  C u0 p0 c0 {4,S} {7,S} {20,S} {21,S}
6  C u0 p0 c0 {7,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {11,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {3,S} {12,S} {22,S}
10 C u0 p0 c0 {4,S} {26,S} {27,S} {28,S}
11 C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
12 C u0 p0 c0 {9,S} {29,S} {30,S} {31,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {12,S}
30 H u0 p0 c0 {12,S}
31 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(39096)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,D} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {8,S} {9,S}
6  C u1 p0 c0 {1,D} {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H5(1942)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {4,D} {8,S} {9,S}
4 C u1 p0 c0 {1,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(39292)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 C u0 p0 c0 {2,D} {8,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(39306)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {7,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(39307)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {6,D} {10,S} {11,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(437)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {9,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {13,S}
5  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {2,S} {3,S} {12,S} {22,S}
10 C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
11 C u0 p0 c0 {8,S} {26,S} {27,S} {28,S}
12 C u0 p0 c0 {9,S} {29,S} {30,S} {31,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {11,S}
27 H u0 p0 c0 {11,S}
28 H u0 p0 c0 {11,S}
29 H u0 p0 c0 {12,S}
30 H u0 p0 c0 {12,S}
31 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(39288)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,D} {3,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {3,T} {8,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(39419)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {6,S}
4  C u0 p0 c0 {2,D} {7,S} {8,S}
5  C u0 p0 c0 {3,D} {9,S} {10,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(39155)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,D}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u0 p0 c0 {4,D} {12,S} {13,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(39753)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {6,D} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(42071)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {3,D} {4,S} {5,S}
7  C u1 p0 c0 {4,S} {11,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(39914)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
5  C u0 p0 c0 {3,D} {4,S} {6,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u1 p0 c0 {4,S} {10,S} {11,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(42076)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {3,S} {4,S} {6,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(39669)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,D} {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u1 p0 c0 {4,D} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(39686)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,D}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,D} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {11,S} {12,S}
7  C u1 p0 c0 {2,D} {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(42564)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {4,S} {7,S}
7  C u1 p0 c0 {2,S} {6,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(39912)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {3,D} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u1 p0 c0 {4,S} {10,S} {11,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(441)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {9,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {3,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(16582)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {7,S}
7  C u0 p0 c0 {3,D} {6,S} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='CH2O2(442)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(443)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(444)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {3,S} {4,S} {17,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(445)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {6,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(446)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(447)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u0 p0 c0 {7,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {6,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(448)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(449)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {1,S} {5,S} {17,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(450)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {5,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(451)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {3,S} {4,S} {17,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(453)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(454)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {6,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {4,S} {16,S} {17,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(455)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
6  C u1 p0 c0 {3,S} {4,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(43250)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(16587)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,D}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {5,D} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(456)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
5  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {5,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(457)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {15,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {1,S} {4,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(458)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {15,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
6  C u1 p0 c0 {3,S} {4,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(43209)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,D}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,D} {4,S} {6,S}
6  C u0 p0 c0 {2,D} {5,S} {7,S}
7  C u0 p0 c0 {3,D} {6,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(460)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {15,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(714)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {27,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
11 C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
27 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H9O(461)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(462)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(463)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {16,S}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(43235)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {6,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(27195)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {12,S}
2  O u0 p2 c0 {5,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {2,D} {4,S} {6,S}
6  C u0 p0 c0 {3,D} {5,S} {7,S}
7  C u1 p0 c0 {1,S} {6,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2HO(467)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u1 p0 c0 {2,T}
4 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2HO(468)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {3,D} {4,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(469)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {4,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(470)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {12,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {1,S} {3,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(39920)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {8,D}
7  C u0 p0 c0 {3,S} {5,D} {6,S}
8  C u0 p0 c0 {6,D} {12,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(471)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {4,S} {11,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C3H7O(473)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
        """),
)


species(
    label='CO3t2(474)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {1,S} {3,D}
        """),
)


species(
    label='CO3t1(475)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
        """),
)


species(
    label='CHO3(476)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {1,S} {3,D}
5 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH3O2(478)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {6,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CH3O2(479)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {6,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(480)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {24,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {4,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2H2O(481)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(39051)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {7,S} {13,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
5  C u1 p0 c0 {4,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {2,S} {3,D} {4,S}
8  C u0 p0 c0 {1,S} {6,D} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(44093)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {8,S} {12,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,D} {8,S}
5  C u0 p0 c0 {4,D} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {1,S} {6,D} {11,S}
8  C u0 p0 c0 {2,S} {3,D} {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(482)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {3,S} {8,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C2H4O(483)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u0 p2 c0 {2,S} {7,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(8475)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {17,S}
3  O u0 p2 c0 {8,S} {18,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {8,D} {15,S} {16,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(14960)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {10,S} {19,S} {20,S}
10 C u0 p0 c0 {2,S} {9,S} {21,S} {22,S}
11 C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C2H4O(484)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(45761)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {1,S} {26,S}
4  O u0 p2 c0 {2,S} {27,S}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {10,S} {19,S} {20,S}
10 C u0 p0 c0 {2,S} {9,S} {21,S} {22,S}
11 C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {11,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {11,S}
26 H u0 p0 c0 {3,S}
27 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(45625)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {10,S}
3  O u0 p2 c0 {11,S} {22,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {14,S}
6  C u1 p0 c0 {2,S} {5,S} {11,S}
7  C u0 p0 c0 {1,S} {9,S} {12,D}
8  C u0 p0 c0 {5,S} {10,D} {15,S}
9  C u0 p0 c0 {7,S} {13,D} {16,S}
10 C u0 p0 c0 {2,S} {8,D} {17,S}
11 C u0 p0 c0 {3,S} {4,D} {6,S}
12 C u0 p0 c0 {7,D} {18,S} {19,S}
13 C u0 p0 c0 {9,D} {20,S} {21,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {12,S}
19 H u0 p0 c0 {12,S}
20 H u0 p0 c0 {13,S}
21 H u0 p0 c0 {13,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='CHPD(1)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,D} {14,S}
5  C u0 p0 c0 {3,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {7,S} {16,S}
7  C u0 p0 c0 {4,D} {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H9(3)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u1 p0 c0 {1,S} {2,S} {12,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u0 p0 c0 {2,S} {7,D} {14,S}
6  C u0 p0 c0 {4,D} {7,S} {15,S}
7  C u0 p0 c0 {5,D} {6,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H10(4)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,D} {10,S}
3  C u0 p0 c0 {1,S} {6,D} {11,S}
4  C u0 p0 c0 {2,D} {5,S} {13,S}
5  C u0 p0 c0 {4,S} {7,D} {12,S}
6  C u0 p0 c0 {3,D} {16,S} {17,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H11(5)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {5,S} {13,S}
5  C u0 p0 c0 {4,S} {7,D} {14,S}
6  C u1 p0 c0 {2,S} {15,S} {16,S}
7  C u0 p0 c0 {5,D} {17,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H7O2(6)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {15,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C6H12O(7)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,S} {19,S}
2  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {7,D} {17,S}
7  C u0 p0 c0 {1,S} {6,D} {18,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9O2(8)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {5,S} {9,D} {16,S}
8  C u0 p0 c0 {6,D} {9,S} {17,S}
9  C u0 p0 c0 {7,D} {8,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H9O(9)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u0 p0 c0 {6,D} {7,S} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H9O2(2)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  C u0 p0 c0 {3,S} {8,D} {14,S}
7  C u0 p0 c0 {4,S} {9,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H8(11)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,D} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {2,D} {7,S} {14,S}
5  C u0 p0 c0 {3,D} {6,S} {15,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {4,S} {6,D} {13,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(12)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {5,S} {9,D} {16,S}
8  C u0 p0 c0 {6,D} {9,S} {17,S}
9  C u0 p0 c0 {7,D} {8,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H8O(1)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
4  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {13,S}
6  C u0 p0 c0 {3,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {15,S}
8  C u0 p0 c0 {6,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H8O(2)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {7,D} {14,S}
6  C u0 p0 c0 {3,S} {8,D} {15,S}
7  C u0 p0 c0 {4,S} {5,D} {13,S}
8  C u0 p0 c0 {4,S} {6,D} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H9O4(15)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {10,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {10,D} {16,S}
9  C u0 p0 c0 {5,S} {11,D} {19,S}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H8O3(16)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {18,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {16,S}
7  C u0 p0 c0 {3,D} {4,S} {10,S}
8  C u0 p0 c0 {5,S} {6,D} {15,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(17)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {1,S} {5,S} {18,S} {19,S}
10 C u0 p0 c0 {5,S} {11,D} {20,S}
11 C u0 p0 c0 {8,S} {10,D} {21,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(18)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {22,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {9,S} {20,S} {21,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H10(950)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {10,S}
3  C u0 p0 c0 {2,D} {4,S} {12,S}
4  C u0 p0 c0 {3,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {7,S} {11,S}
6  C u1 p0 c0 {1,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H10(953)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {7,D} {16,S}
7  C u0 p0 c0 {1,S} {6,D} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H9(937)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {6,D} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {7,S} {15,S}
6  C u0 p0 c0 {3,D} {7,S} {16,S}
7  C u1 p0 c0 {5,S} {6,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H11(954)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
5  C u1 p0 c0 {4,S} {7,S} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {16,S}
7  C u0 p0 c0 {5,S} {6,D} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1027)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,D} {16,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {9,S} {18,S}
9  C u0 p0 c0 {6,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1107)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {9,D} {20,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1026)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
6  C u0 p0 c0 {4,S} {9,D} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {15,S}
8  C u0 p0 c0 {5,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {6,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(1182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {9,S} {17,S}
5  C u0 p0 c0 {3,S} {6,S} {11,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {9,S} {18,S} {19,S}
9  C u1 p0 c0 {4,S} {8,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(1136)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {21,S}
3  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {16,S}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {5,S} {9,D} {20,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {2,S}
        """),
)







species(
    label='S(1423)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,D} {2,S} {3,S}
5  C u0 p0 c0 {2,S} {7,D} {13,S}
6  C u0 p0 c0 {3,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {15,S}
8  C u0 p0 c0 {6,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H9(1021)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {7,S} {14,S}
6  C u0 p0 c0 {1,S} {7,D} {15,S}
7  C u0 p0 c0 {5,S} {6,D} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1395)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {6,D} {13,S} {14,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1147)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {16,S}
6  C u1 p0 c0 {3,S} {9,S} {15,S}
7  C u0 p0 c0 {4,S} {5,D} {14,S}
8  C u0 p0 c0 {4,S} {9,D} {13,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H8(1072)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,D} {14,S}
5  C u0 p0 c0 {2,S} {7,D} {15,S}
6  C u0 p0 c0 {2,S} {4,D} {13,S}
7  C u0 p0 c0 {3,S} {5,D} {12,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(1189)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {8,D} {18,S}
8  C u0 p0 c0 {6,S} {7,D} {19,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1660)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {8,S} {16,S}
7  C u1 p0 c0 {4,S} {17,S} {18,S}
8  C u0 p0 c0 {1,D} {6,S} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2171)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,D} {18,S}
8  C u0 p0 c0 {1,S} {7,D} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2308)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {8,S} {12,S} {17,S}
7  C u1 p0 c0 {4,S} {5,S} {18,S}
8  C u0 p0 c0 {1,D} {6,S} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2380)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {9,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {4,S} {10,S} {20,S}
10 C u0 p0 c0 {2,D} {9,S} {21,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2504)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {9,S} {11,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {4,S} {10,S} {19,S} {20,S}
10 C u0 p0 c0 {2,D} {9,S} {21,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2317)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,S} {20,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {2,S} {8,D} {18,S}
8  C u0 p0 c0 {1,S} {7,D} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2325)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {22,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {1,S} {4,S} {10,S} {20,S}
8  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
10 C u0 p0 c0 {3,D} {7,S} {21,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1122)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {16,S}
6  C u1 p0 c0 {3,S} {9,S} {15,S}
7  C u0 p0 c0 {4,S} {5,D} {14,S}
8  C u0 p0 c0 {4,S} {9,D} {13,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1558)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {8,D} {17,S}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u0 p0 c0 {4,S} {6,D} {15,S}
9  C u0 p0 c0 {5,S} {7,D} {14,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(2823)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {13,S}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {8,S} {12,S}
7  C u0 p0 c0 {1,D} {2,S} {15,S}
8  C u1 p0 c0 {6,S} {16,S} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C5H6(960)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,D} {3,S} {7,S}
2  C u0 p0 c0 {1,D} {4,S} {6,S}
3  C u0 p0 c0 {1,S} {5,D} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u1 p0 c0 {3,D} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(2292)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u0 p0 c0 {1,D} {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C5H6(2966)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {1,S} {2,D} {9,S}
4  C u0 p0 c0 {1,S} {5,D} {7,S}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H6(118)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {2,D} {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(2795)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u1 p0 c0 {3,S} {7,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {8,S} {14,S}
7  C u0 p0 c0 {4,S} {9,D} {13,S}
8  C u1 p0 c0 {2,S} {6,S} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3354)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u1 p0 c0 {3,S} {7,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {8,S} {14,S}
7  C u0 p0 c0 {4,S} {9,D} {13,S}
8  C u0 p0 c0 {2,D} {6,S} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3382)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,D} {3,S} {6,S}
5  C u0 p0 c0 {3,S} {8,D} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {9,S} {14,S}
8  C u0 p0 c0 {5,D} {16,S} {17,S}
9  C u0 p0 c0 {2,D} {7,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1055)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  C u0 p0 c0 {2,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {16,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1190)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {17,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,D} {13,S}
5  C u0 p0 c0 {3,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {8,S} {16,S}
7  C u0 p0 c0 {4,D} {8,S} {15,S}
8  C u1 p0 c0 {6,S} {7,S} {14,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2980)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {15,S}
7  C u0 p0 c0 {5,S} {8,D} {16,S}
8  C u0 p0 c0 {1,S} {7,D} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2999)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {20,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {18,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u0 p0 c0 {2,D} {7,S} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4179)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {6,S} {10,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {18,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4177)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {5,S} {9,D} {15,S}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 C u0 p0 c0 {6,S} {7,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4178)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2929)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {17,S}
2  C u0 p0 c0 {3,D} {5,S} {10,S}
3  C u0 p0 c0 {2,D} {4,S} {12,S}
4  C u0 p0 c0 {3,S} {6,D} {13,S}
5  C u0 p0 c0 {2,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {8,S} {11,S}
7  C u0 p0 c0 {1,S} {5,D} {14,S}
8  C u1 p0 c0 {6,S} {15,S} {16,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4584)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,D} {13,S}
6  C u0 p0 c0 {4,S} {10,D} {12,S}
7  C u0 p0 c0 {5,D} {8,S} {15,S}
8  C u0 p0 c0 {7,S} {9,D} {14,S}
9  C u0 p0 c0 {2,S} {8,D} {16,S}
10 C u0 p0 c0 {6,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4586)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {9,S} {16,S}
9  C u0 p0 c0 {8,S} {10,D} {14,S}
10 C u0 p0 c0 {9,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4585)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,D} {13,S}
6  C u0 p0 c0 {4,S} {9,D} {12,S}
7  C u0 p0 c0 {5,D} {8,S} {15,S}
8  C u0 p0 c0 {7,S} {10,D} {14,S}
9  C u0 p0 c0 {2,S} {6,D} {16,S}
10 C u0 p0 c0 {8,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2907)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,D} {5,S} {10,S}
3  C u0 p0 c0 {2,D} {4,S} {11,S}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {2,S} {7,D} {9,S}
6  C u0 p0 c0 {4,D} {8,S} {13,S}
7  C u0 p0 c0 {5,D} {14,S} {15,S}
8  C u0 p0 c0 {1,D} {6,S} {16,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(4183)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,S} {22,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
8  C u0 p0 c0 {2,S} {5,S} {11,S} {20,S}
9  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
10 C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
11 C u0 p0 c0 {3,D} {8,S} {21,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4181)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {17,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 C u0 p0 c0 {2,D} {7,S} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(4182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {2,D} {9,S} {19,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2408)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
8  C u1 p0 c0 {1,S} {7,S} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(5142)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {15,S}
5  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {8,S} {9,D}
8  C u1 p0 c0 {5,S} {7,S} {16,S}
9  C u0 p0 c0 {6,S} {7,D} {17,S}
10 C u0 p0 c0 {3,D} {4,S} {18,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(5164)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {8,D} {9,S}
8  C u0 p0 c0 {6,S} {7,D} {16,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {3,S} {9,D} {18,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C5H7(596)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u1 p0 c0 {2,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {5,D} {11,S}
5  C u0 p0 c0 {3,S} {4,D} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(3269)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H9(1521)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {1,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {6,S} {13,S}
6  C u0 p0 c0 {4,D} {5,S} {14,S}
7  C u1 p0 c0 {2,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H5(98)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u1 p0 c0 {2,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,D} {7,S}
3  C u0 p0 c0 {2,D} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {1,S} {4,D} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(5256)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {16,S}
9  C u0 p0 c0 {2,D} {5,S} {17,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(5475)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {4,D} {7,S} {11,S}
7  C u0 p0 c0 {5,D} {6,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(5116)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {10,S} {21,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {8,S} {14,S} {15,S}
10 C u0 p0 c0 {2,S} {5,S} {11,S} {20,S}
11 C u1 p0 c0 {4,D} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4588)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u1 p0 c0 {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {9,D} {12,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {10,S} {14,S}
9  C u0 p0 c0 {6,D} {17,S} {18,S}
10 C u0 p0 c0 {2,D} {8,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3423)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u1 p0 c0 {4,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {6,D} {9,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {14,S}
9  C u0 p0 c0 {2,D} {7,S} {16,S}
10 C u0 p0 c0 {8,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(5668)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {7,S} {16,S} {17,S}
9  C u1 p0 c0 {2,S} {4,S} {18,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3863)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,S} {17,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {1,S} {6,D} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(5102)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {14,S}
8  C u0 p0 c0 {7,D} {10,S} {13,S}
9  C u0 p0 c0 {2,D} {4,S} {16,S}
10 C u1 p0 c0 {8,S} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(5729)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {14,S}
8  C u0 p0 c0 {5,S} {6,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {3,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(5680)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {5,S} {6,D} {15,S}
8  C u0 p0 c0 {4,S} {10,D} {13,S}
9  C u0 p0 c0 {3,D} {5,S} {16,S}
10 C u0 p0 c0 {8,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(3425)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {14,S}
8  C u0 p0 c0 {5,S} {6,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {3,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(5535)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {4,D} {7,S} {11,S}
7  C u0 p0 c0 {5,D} {6,S} {12,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(6339)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u1 p0 c0 {3,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {5,S} {6,D} {12,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(6025)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {5,S} {6,D} {15,S}
8  C u0 p0 c0 {4,S} {10,D} {13,S}
9  C u0 p0 c0 {3,D} {5,S} {16,S}
10 C u0 p0 c0 {8,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(1590)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {8,D}
5  C u0 p0 c0 {2,S} {6,D} {13,S}
6  C u0 p0 c0 {3,S} {5,D} {12,S}
7  C u1 p0 c0 {3,S} {8,S} {14,S}
8  C u0 p0 c0 {4,D} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(6567)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {2,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(6576)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {8,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {6,S} {9,S} {14,S}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {15,S}
8  C u0 p0 c0 {1,S} {7,S} {10,D}
9  C u1 p0 c0 {4,S} {7,S} {16,S}
10 C u0 p0 c0 {5,S} {8,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(6568)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
6  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {10,D} {17,S}
8  C u0 p0 c0 {5,S} {9,D} {15,S}
9  C u0 p0 c0 {6,S} {8,D} {14,S}
10 C u0 p0 c0 {6,S} {7,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(6577)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u1 p0 c0 {5,S} {7,S} {16,S}
10 C u0 p0 c0 {6,S} {8,D} {17,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(6673)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u1 p0 c0 {6,S} {7,S} {15,S}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(6575)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {4,S} {10,S} {14,S} {15,S}
8  C u1 p0 c0 {1,S} {5,S} {6,S}
9  C u0 p0 c0 {6,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(6463)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,D} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {7,S} {14,S}
6  C u0 p0 c0 {4,D} {8,S} {15,S}
7  C u0 p0 c0 {5,S} {8,D} {13,S}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(6441)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {2,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {7,S} {11,S}
6  C u0 p0 c0 {4,D} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(7159)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,D} {15,S}
7  C u0 p0 c0 {5,S} {8,D} {14,S}
8  C u0 p0 c0 {7,D} {10,S} {17,S}
9  C u0 p0 c0 {6,D} {10,S} {16,S}
10 C u0 p0 c0 {2,D} {8,S} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(7158)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {16,S}
7  C u0 p0 c0 {5,S} {6,D} {14,S}
8  C u0 p0 c0 {2,D} {4,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {15,S}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(6769)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u1 p0 c0 {5,S} {8,S} {12,S}
3  C u0 p0 c0 {4,S} {5,D} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {2,S} {3,D} {10,S}
6  C u0 p0 c0 {4,D} {7,S} {11,S}
7  C u0 p0 c0 {6,S} {8,D} {13,S}
8  C u0 p0 c0 {2,S} {7,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(6460)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,D} {2,S} {3,S}
6  C u1 p0 c0 {2,S} {4,S} {13,S}
7  C u0 p0 c0 {3,S} {8,D} {15,S}
8  C u0 p0 c0 {4,S} {7,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(7425)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {2,S} {4,D} {13,S}
6  C u0 p0 c0 {1,D} {3,S} {8,S}
7  C u0 p0 c0 {2,S} {8,D} {11,S}
8  C u0 p0 c0 {6,S} {7,D} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(7426)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {7,S}
4  C u0 p0 c0 {2,S} {3,D} {10,S}
5  C u0 p0 c0 {2,S} {6,D} {11,S}
6  C u0 p0 c0 {5,D} {8,S} {14,S}
7  C u0 p0 c0 {3,S} {8,D} {12,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(6844)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {1,D} {2,S} {3,S}
5  C u0 p0 c0 {2,S} {7,D} {11,S}
6  C u0 p0 c0 {3,S} {8,D} {12,S}
7  C u0 p0 c0 {5,D} {8,S} {14,S}
8  C u0 p0 c0 {6,D} {7,S} {13,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(6794)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u0 p0 c0 {2,S} {7,D} {15,S}
6  C u0 p0 c0 {3,S} {4,D} {12,S}
7  C u0 p0 c0 {3,S} {5,D} {13,S}
8  C u1 p0 c0 {1,D} {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C6H7(363)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u1 p0 c0 {1,S} {4,S} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {6,S} {13,S}
6  C u0 p0 c0 {4,D} {5,S} {11,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H5O(143)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {2,S} {6,D} {9,S}
5  C u0 p0 c0 {3,D} {6,S} {10,S}
6  C u0 p0 c0 {4,D} {5,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(8088)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {12,S}
7  C u0 p0 c0 {2,D} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(8087)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {4,S} {8,D} {13,S}
8  C u0 p0 c0 {5,S} {7,D} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(7883)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {8,D} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {6,D} {8,S} {15,S}
8  C u0 p0 c0 {5,D} {7,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(7884)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
4  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u0 p0 c0 {3,S} {8,D} {15,S}
7  C u0 p0 c0 {4,S} {5,D} {12,S}
8  C u0 p0 c0 {4,S} {6,D} {13,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H7(7867)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  C u1 p0 c0 {2,S} {6,S} {12,S}
5  C u0 p0 c0 {1,S} {6,D} {11,S}
6  C u0 p0 c0 {4,S} {5,D} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H6(100)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,B} {6,B} {7,S}
2  C u0 p0 c0 {1,B} {3,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {11,S}
6  C u0 p0 c0 {1,B} {5,B} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(8095)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {4,S} {5,S} {13,S}
8  C u0 p0 c0 {3,D} {5,S} {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H4O(142)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,D}
2  C u0 p0 c0 {1,D} {3,S} {6,S}
3  C u0 p0 c0 {2,S} {4,D} {7,S}
4  C u0 p0 c0 {3,D} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,D} {9,S}
6  C u0 p0 c0 {2,S} {5,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(7443)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {3,D} {7,S} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {14,S}
7  C u1 p0 c0 {5,S} {6,S} {13,S}
8  C u1 p0 c0 {1,D} {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(3596)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {11,D} {16,S}
9  C u0 p0 c0 {6,S} {10,D} {15,S}
10 C u0 p0 c0 {9,D} {18,S} {19,S}
11 C u0 p0 c0 {3,S} {8,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C6H6(7863)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u1 p0 c0 {2,S} {5,S} {11,S}
2  C u1 p0 c0 {1,S} {6,S} {12,S}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {3,D} {9,S}
6  C u0 p0 c0 {2,S} {4,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(8171)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {13,S}
7  C u0 p0 c0 {5,D} {8,S} {14,S}
8  C u0 p0 c0 {6,D} {7,S} {12,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(2953)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {16,S}
7  C u0 p0 c0 {5,S} {6,D} {14,S}
8  C u0 p0 c0 {4,S} {10,D} {15,S}
9  C u0 p0 c0 {2,D} {5,S} {17,S}
10 C u0 p0 c0 {8,D} {18,S} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2954)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {8,S} {16,S}
8  C u0 p0 c0 {7,S} {10,D} {15,S}
9  C u0 p0 c0 {2,D} {5,S} {17,S}
10 C u0 p0 c0 {8,D} {18,S} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(4610)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u1 p0 c0 {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {9,D} {12,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {10,S} {14,S}
9  C u0 p0 c0 {6,D} {17,S} {18,S}
10 C u0 p0 c0 {3,D} {8,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3429)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u1 p0 c0 {4,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {6,D} {9,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {14,S}
9  C u0 p0 c0 {3,D} {7,S} {16,S}
10 C u0 p0 c0 {8,D} {17,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3385)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u0 p0 c0 {4,S} {8,D} {13,S}
7  C u0 p0 c0 {5,D} {9,S} {14,S}
8  C u0 p0 c0 {6,D} {16,S} {17,S}
9  C u0 p0 c0 {2,D} {7,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(8779)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {7,S}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {4,S} {7,D} {10,S}
4  C u0 p0 c0 {3,S} {9,D} {12,S}
5  C u0 p0 c0 {6,S} {8,D} {11,S}
6  C u1 p0 c0 {1,S} {5,S} {14,S}
7  C u0 p0 c0 {1,S} {3,D} {13,S}
8  C u0 p0 c0 {5,D} {16,S} {17,S}
9  C u0 p0 c0 {2,S} {4,D} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(8926)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {9,S} {15,S}
8  C u0 p0 c0 {1,S} {6,D} {16,S}
9  C u0 p0 c0 {2,D} {7,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(8925)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u0 p0 c0 {4,S} {9,D} {13,S}
7  C u0 p0 c0 {1,S} {5,D} {14,S}
8  C u0 p0 c0 {2,D} {3,S} {15,S}
9  C u0 p0 c0 {6,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3570)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {18,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {9,D} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {8,D} {14,S}
8  C u0 p0 c0 {1,S} {7,D} {15,S}
9  C u0 p0 c0 {5,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(8137)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u1 p0 c0 {4,S} {8,S} {10,S}
6  C u0 p0 c0 {2,D} {4,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {12,S}
8  C u0 p0 c0 {5,S} {7,D} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(5511)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {9,S} {17,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {7,S} {8,D}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {6,D} {15,S}
9  C u0 p0 c0 {2,S} {7,D} {16,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1604)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {8,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {3,S} {8,D} {13,S}
6  C u0 p0 c0 {4,D} {7,S} {14,S}
7  C u1 p0 c0 {1,S} {6,S} {16,S}
8  C u0 p0 c0 {1,S} {5,D} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(9457)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {9,S} {17,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {6,S} {7,D} {8,S}
6  C u1 p0 c0 {3,S} {5,S} {14,S}
7  C u0 p0 c0 {4,S} {5,D} {15,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {1,S} {8,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(9465)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {17,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {8,S} {9,D}
8  C u0 p0 c0 {2,D} {5,S} {7,S}
9  C u0 p0 c0 {6,S} {7,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1209)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {9,S} {17,S}
8  C u0 p0 c0 {6,S} {9,D} {16,S}
9  C u0 p0 c0 {7,S} {8,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(9656)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {6,D} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {5,D} {15,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {9,D} {14,S}
9  C u0 p0 c0 {5,S} {8,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(10145)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {6,S} {10,S} {15,S}
9  C u0 p0 c0 {3,S} {7,S} {11,S} {18,S}
10 C u0 p0 c0 {8,S} {11,D} {19,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(9964)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {8,D}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {3,S} {9,D} {14,S}
8  C u0 p0 c0 {5,D} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(10199)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {3,S} {8,D} {13,S}
7  C u0 p0 c0 {3,S} {9,D} {14,S}
8  C u0 p0 c0 {6,D} {9,S} {15,S}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(10144)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {16,S}
9  C u0 p0 c0 {6,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {8,S} {11,D} {20,S}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(7901)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u0 p0 c0 {2,S} {6,D} {10,S}
5  C u0 p0 c0 {1,S} {3,D} {6,S}
6  C u0 p0 c0 {4,D} {5,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(9472)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {8,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {8,D} {14,S}
6  C u0 p0 c0 {3,S} {7,D} {15,S}
7  C u0 p0 c0 {4,S} {6,D} {13,S}
8  C u0 p0 c0 {1,S} {5,D} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(8372)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {4,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(10568)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u1 p0 c0 {2,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {5,D} {14,S}
8  C u1 p0 c0 {6,S} {15,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(10726)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {3,S} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {13,S}
6  C u0 p0 c0 {3,S} {8,D} {12,S}
7  C u0 p0 c0 {1,S} {5,D} {14,S}
8  C u0 p0 c0 {6,D} {15,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1585)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {6,D} {13,S}
6  C u0 p0 c0 {3,S} {5,D} {14,S}
7  C u0 p0 c0 {3,S} {8,D} {15,S}
8  C u0 p0 c0 {4,S} {7,D} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1397)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {12,S} {13,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(10864)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u1 p0 c0 {2,S} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {6,D} {11,S} {12,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(11105)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {15,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {8,D}
7  C u0 p0 c0 {5,D} {6,S} {12,S}
8  C u0 p0 c0 {6,D} {13,S} {14,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(11104)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {15,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {12,S}
7  C u0 p0 c0 {8,D} {13,S} {14,S}
8  C u0 p0 c0 {6,D} {7,D}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(11116)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {4,S}
3  O u0 p2 c0 {5,S} {15,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u1 p0 c0 {4,S} {13,S} {14,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(11125)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {14,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {12,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(11364)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {16,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {7,D} {8,S} {15,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(11373)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {7,S} {16,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {2,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {5,S} {7,S} {15,S}
9  C u0 p0 c0 {4,D} {5,S} {6,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(10567)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {7,S} {12,S}
5  C u0 p0 c0 {2,S} {8,D} {14,S}
6  C u0 p0 c0 {3,S} {7,D} {13,S}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  C u0 p0 c0 {1,S} {5,D} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(11278)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {14,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {11,S}
7  C u0 p0 c0 {1,S} {6,D} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(11483)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,D} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u0 p0 c0 {2,D} {6,S} {13,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(10210)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {6,S} {18,S}
3  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,D} {8,S}
6  C u0 p0 c0 {2,S} {4,S} {5,D}
7  C u0 p0 c0 {3,S} {9,D} {14,S}
8  C u1 p0 c0 {5,S} {9,S} {15,S}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(11446)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {8,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {12,S}
5  C u0 p0 c0 {2,S} {6,D} {13,S}
6  C u0 p0 c0 {4,S} {5,D} {15,S}
7  C u0 p0 c0 {2,S} {8,D} {14,S}
8  C u0 p0 c0 {1,S} {7,D} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(10576)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {7,D} {15,S}
7  C u0 p0 c0 {5,S} {6,D} {14,S}
8  C u0 p0 c0 {1,D} {4,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(11847)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {11,S}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(6034)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {10,S} {19,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u1 p0 c0 {4,S} {9,S} {14,S}
8  C u0 p0 c0 {5,S} {6,D} {16,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {3,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3657)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {14,S}
7  C u0 p0 c0 {3,D} {5,S} {10,S}
8  C u0 p0 c0 {6,S} {10,D} {16,S}
9  C u0 p0 c0 {5,S} {11,D} {15,S}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {9,D} {18,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(10227)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
6  C u0 p0 c0 {8,S} {11,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {7,S} {8,D}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {6,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(11846)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
6  C u0 p0 c0 {8,S} {11,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {6,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(10230)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {7,S} {20,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {6,S} {10,D}
9  C u0 p0 c0 {6,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {18,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3569)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {8,D} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,D} {14,S}
8  C u0 p0 c0 {5,D} {16,S} {17,S}
9  C u0 p0 c0 {2,S} {7,D} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(12205)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
6  C u0 p0 c0 {4,S} {10,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,D} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {6,S} {9,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(9078)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {9,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {14,S}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {1,S} {6,D} {17,S}
9  C u0 p0 c0 {2,S} {7,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(9433)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {10,S} {20,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {15,S}
10 C u0 p0 c0 {2,S} {8,D} {17,S}
11 C u0 p0 c0 {9,D} {18,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3597)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u1 p2 c0 {8,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {8,D} {16,S}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u0 p0 c0 {6,S} {11,D} {15,S}
10 C u0 p0 c0 {3,D} {5,S} {17,S}
11 C u0 p0 c0 {9,D} {18,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(2819)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {19,S}
9  C u0 p0 c0 {6,S} {11,D} {18,S}
10 C u0 p0 c0 {6,S} {8,D} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(6674)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {6,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
8  C u1 p0 c0 {5,S} {6,S} {15,S}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(12794)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {20,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {15,S}
10 C u0 p0 c0 {3,D} {5,S} {17,S}
11 C u0 p0 c0 {9,D} {18,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(12545)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {7,D}
7  C u0 p0 c0 {3,S} {5,S} {6,D}
8  C u0 p0 c0 {4,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {16,S}
10 C u1 p0 c0 {1,S} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(9125)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {11,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {7,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {3,S} {9,S} {15,S}
8  C u0 p0 c0 {5,S} {9,D} {16,S}
9  C u0 p0 c0 {7,S} {8,D} {18,S}
10 C u0 p0 c0 {6,S} {11,D} {17,S}
11 C u0 p0 c0 {1,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(12558)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {17,S}
3  O u0 p2 c0 {10,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
6  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {4,S} {10,D}
8  C u0 p0 c0 {5,S} {9,D} {15,S}
9  C u0 p0 c0 {6,S} {8,D} {14,S}
10 C u0 p0 c0 {3,S} {7,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(6317)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {5,S} {15,S}
8  C u0 p0 c0 {5,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 C u0 p0 c0 {3,D} {4,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(12938)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {7,S} {10,S}
3  O u0 p2 c0 {6,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {1,S} {6,S} {7,D}
9  C u0 p0 c0 {4,S} {10,D} {16,S}
10 C u0 p0 c0 {2,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(12937)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {7,D} {8,S}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  C u0 p0 c0 {3,D} {5,S} {6,S}
9  C u0 p0 c0 {4,S} {10,D} {16,S}
10 C u0 p0 c0 {1,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(12555)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {6,S} {18,S}
4  C u0 p0 c0 {1,S} {2,S} {9,S} {11,S}
5  C u0 p0 c0 {7,S} {10,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {2,S} {6,S} {7,D}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {5,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(12939)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {11,S}
6  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,D} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {1,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3652)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {6,S} {10,D} {14,S}
9  C u0 p0 c0 {7,D} {11,S} {16,S}
10 C u0 p0 c0 {8,D} {18,S} {19,S}
11 C u0 p0 c0 {3,D} {9,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(9900)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {17,S}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(8981)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u1 p0 c0 {3,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {9,S} {13,S}
8  C u0 p0 c0 {1,S} {6,D} {15,S}
9  C u0 p0 c0 {2,D} {7,S} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(1496)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {16,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {3,D} {7,S} {12,S}
6  C u0 p0 c0 {4,D} {8,S} {15,S}
7  C u0 p0 c0 {5,S} {8,D} {13,S}
8  C u0 p0 c0 {6,S} {7,D} {14,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(13750)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {4,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u1 p0 c0 {3,S} {8,S} {12,S}
7  C u0 p0 c0 {4,S} {5,D} {14,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {2,S} {8,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13818)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {14,S}
9  C u0 p0 c0 {6,S} {7,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {3,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(13927)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {13,S}
7  C u0 p0 c0 {1,S} {3,S} {10,S} {14,S}
8  C u1 p0 c0 {5,S} {6,S} {15,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 C u0 p0 c0 {4,D} {6,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(13817)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {6,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {11,S} {16,S}
10 C u0 p0 c0 {1,S} {8,D} {17,S}
11 C u0 p0 c0 {3,D} {9,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(13868)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {17,S}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(13736)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {4,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
5  C u1 p0 c0 {3,S} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {8,D} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {1,S} {6,D} {15,S}
9  C u0 p0 c0 {2,S} {7,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3386)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u0 p0 c0 {3,S} {8,D} {13,S}
7  C u0 p0 c0 {4,S} {5,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u0 p0 c0 {2,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(4224)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {10,D}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u1 p0 c0 {4,S} {10,S} {14,S}
9  C u0 p0 c0 {5,S} {7,D} {15,S}
10 C u0 p0 c0 {6,D} {8,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(13821)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  C u1 p0 c0 {4,S} {10,S} {16,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(13723)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {7,D} {11,S}
5  C u0 p0 c0 {3,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {8,D} {13,S}
7  C u0 p0 c0 {4,D} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {6,D} {14,S}
9  C u0 p0 c0 {2,D} {7,S} {15,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13877)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u0 p0 c0 {6,S} {9,D} {16,S}
9  C u0 p0 c0 {4,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(14214)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {7,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {6,S} {13,S}
6  C u1 p0 c0 {5,S} {7,S} {14,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {2,D} {3,S} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(14206)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {13,S}
5  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u0 p0 c0 {2,D} {6,S} {9,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(14213)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {6,S} {14,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {9,S} {12,S}
8  C u0 p0 c0 {2,D} {3,S} {15,S}
9  C u1 p0 c0 {7,S} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(14837)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {1,D} {3,S} {5,S}
7  C u0 p0 c0 {3,S} {8,D} {16,S}
8  C u0 p0 c0 {4,S} {7,D} {15,S}
9  C u0 p0 c0 {2,D} {5,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(13344)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {7,S} {17,S}
3  O u0 p2 c0 {10,S} {18,S}
4  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {10,D}
8  C u0 p0 c0 {4,S} {9,D} {14,S}
9  C u0 p0 c0 {1,S} {8,D} {16,S}
10 C u0 p0 c0 {3,S} {7,D} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(11447)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {7,D} {14,S}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  C u0 p0 c0 {1,D} {3,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(14679)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {11,S} {14,S}
10 C u0 p0 c0 {1,S} {8,D} {16,S}
11 C u0 p0 c0 {3,D} {9,S} {17,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(13770)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {16,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {2,S} {3,D} {11,S}
5  C u0 p0 c0 {2,S} {7,D} {12,S}
6  C u0 p0 c0 {3,S} {8,D} {13,S}
7  C u0 p0 c0 {5,D} {8,S} {15,S}
8  C u0 p0 c0 {6,D} {7,S} {14,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(14182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {18,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {5,S} {10,D} {13,S}
8  C u0 p0 c0 {6,D} {9,S} {14,S}
9  C u0 p0 c0 {8,S} {11,D} {15,S}
10 C u0 p0 c0 {1,S} {7,D} {16,S}
11 C u0 p0 c0 {4,S} {9,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4240)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {1,S} {19,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {10,D}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u1 p0 c0 {4,S} {10,S} {15,S}
9  C u0 p0 c0 {5,S} {7,D} {14,S}
10 C u0 p0 c0 {6,D} {8,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(6855)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {16,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {7,D}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {2,S} {8,D} {12,S}
6  C u0 p0 c0 {3,S} {4,D} {14,S}
7  C u0 p0 c0 {3,D} {8,S} {13,S}
8  C u0 p0 c0 {5,D} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(14951)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {3,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u0 p0 c0 {6,D} {7,S} {15,S}
9  C u0 p0 c0 {2,D} {3,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3799)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u0 p0 c0 {2,S} {5,D} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {3,D} {7,S} {12,S}
6  C u0 p0 c0 {4,D} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {8,D} {14,S}
8  C u0 p0 c0 {6,S} {7,D} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(14950)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {9,D}
7  C u0 p0 c0 {5,D} {6,S} {14,S}
8  C u0 p0 c0 {2,D} {4,S} {15,S}
9  C u0 p0 c0 {6,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(14952)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,D}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {9,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {8,D} {15,S}
7  C u0 p0 c0 {1,D} {4,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {16,S}
9  C u0 p0 c0 {2,D} {5,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(12935)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {4,S} {17,S}
3  O u0 p2 c0 {5,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {1,S} {5,D} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {13,S}
8  C u0 p0 c0 {7,D} {9,S} {14,S}
9  C u0 p0 c0 {8,S} {10,D} {15,S}
10 C u0 p0 c0 {1,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(15164)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,S} {17,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {16,S}
8  C u0 p0 c0 {3,S} {9,D} {15,S}
9  C u0 p0 c0 {5,S} {8,D} {14,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(15738)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
6  C u1 p0 c0 {4,S} {7,S} {16,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {17,S}
9  C u0 p0 c0 {5,S} {7,D} {14,S}
10 C u0 p0 c0 {5,S} {8,D} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(15736)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {14,S}
8  C u1 p0 c0 {4,S} {10,S} {16,S}
9  C u0 p0 c0 {5,S} {10,D} {15,S}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(14682)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {8,D}
7  C u1 p0 c0 {5,S} {6,S} {13,S}
8  C u0 p0 c0 {6,D} {9,S} {14,S}
9  C u0 p0 c0 {8,S} {10,D} {15,S}
10 C u0 p0 c0 {1,S} {9,D} {16,S}
11 C u0 p0 c0 {3,D} {5,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(394)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(14342)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {18,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u1 p0 c0 {3,S} {8,S} {13,S}
7  C u0 p0 c0 {4,S} {5,D} {15,S}
8  C u0 p0 c0 {6,S} {9,D} {16,S}
9  C u0 p0 c0 {2,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1184)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {20,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u1 p0 c0 {6,S} {9,S} {18,S}
9  C u0 p0 c0 {7,D} {8,S} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(6893)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {4,S} {18,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
5  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,D} {16,S}
7  C u1 p0 c0 {4,S} {10,S} {15,S}
8  C u0 p0 c0 {5,S} {6,D} {14,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(14886)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {1,D} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,D} {12,S}
8  C u0 p0 c0 {2,D} {3,S} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(14373)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {18,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u0 p0 c0 {5,S} {11,D} {13,S}
8  C u0 p0 c0 {6,D} {9,S} {14,S}
9  C u0 p0 c0 {8,S} {10,D} {15,S}
10 C u0 p0 c0 {1,S} {9,D} {16,S}
11 C u0 p0 c0 {4,S} {7,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(14327)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {9,S} {16,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {6,S} {11,S}
5  C u0 p0 c0 {3,D} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {9,D} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {13,S}
8  C u0 p0 c0 {1,S} {7,D} {15,S}
9  C u0 p0 c0 {2,S} {6,D} {14,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(12924)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u1 p0 c0 {6,S} {7,S} {17,S}
10 C u0 p0 c0 {5,S} {8,D} {18,S}
11 C u0 p0 c0 {4,D} {5,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(13833)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {5,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {11,S} {16,S}
10 C u1 p0 c0 {7,S} {17,S} {18,S}
11 C u0 p0 c0 {4,D} {9,S} {19,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(16697)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {7,S} {25,S}
4  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
5  C u0 p0 c0 {10,S} {11,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
8  C u0 p0 c0 {4,S} {10,D} {23,S}
9  C u1 p0 c0 {4,S} {12,S} {22,S}
10 C u0 p0 c0 {5,S} {8,D} {21,S}
11 C u0 p0 c0 {5,S} {12,D} {20,S}
12 C u0 p0 c0 {9,S} {11,D} {24,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {12,S}
25 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1649)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,D} {17,S}
7  C u0 p0 c0 {1,D} {5,S} {8,S}
8  C u0 p0 c0 {6,D} {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(17376)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {11,S} {18,S}
9  C u0 p0 c0 {7,S} {10,S} {19,S} {20,S}
10 C u0 p0 c0 {1,S} {9,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {21,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(17377)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {8,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {5,S} {11,D} {21,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1653)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u1 p0 c0 {4,S} {8,S} {17,S}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(14370)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {6,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {10,S} {13,S}
7  C u0 p0 c0 {1,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {11,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {19,S}
11 C u0 p0 c0 {7,S} {8,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(16896)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {20,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {18,S}
8  C u1 p0 c0 {6,S} {7,S} {19,S}
9  C u0 p0 c0 {2,D} {5,S} {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(17375)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {11,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {10,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {18,S}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {3,S} {9,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(16286)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {9,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(16253)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {10,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {8,S} {13,S}
8  C u0 p0 c0 {7,S} {10,D} {14,S}
9  C u0 p0 c0 {3,D} {5,S} {15,S}
10 C u0 p0 c0 {2,S} {8,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(19452)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {20,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {3,S} {4,S}
7  C u0 p0 c0 {4,S} {9,D} {15,S}
8  C u1 p0 c0 {5,S} {16,S} {17,S}
9  C u0 p0 c0 {7,D} {18,S} {19,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(19451)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {20,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {8,D} {18,S}
8  C u0 p0 c0 {1,S} {7,D} {19,S}
9  C u1 p0 c0 {2,D} {6,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(19568)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {22,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {2,S} {5,S} {11,S} {19,S}
8  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {10,S} {15,S} {16,S}
10 C u0 p0 c0 {9,S} {11,S} {20,S} {21,S}
11 C u0 p0 c0 {3,D} {7,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(19454)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {20,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {2,D} {3,S} {7,S}
9  C u1 p0 c0 {1,S} {3,S} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(20155)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {18,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u1 p0 c0 {4,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {5,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(20203)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {9,S} {21,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {9,D} {18,S}
8  C u0 p0 c0 {2,D} {6,S} {19,S}
9  C u0 p0 c0 {1,S} {7,D} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(20445)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {11,S} {22,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {11,D} {20,S}
10 C u0 p0 c0 {1,S} {3,D} {8,S}
11 C u0 p0 c0 {2,S} {9,D} {21,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(20818)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {9,S} {22,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {11,S} {18,S} {19,S}
9  C u0 p0 c0 {1,S} {3,S} {10,S} {20,S}
10 C u1 p0 c0 {7,S} {9,S} {21,S}
11 C u0 p0 c0 {2,S} {4,D} {8,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(5811)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {18,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {5,S} {15,S} {16,S}
7  C u1 p0 c0 {1,S} {2,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(20947)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {9,S} {20,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {9,D} {18,S}
9  C u0 p0 c0 {2,S} {8,D} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(20771)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {12,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {6,S} {17,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(20811)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {8,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {10,S} {22,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {1,S} {6,S} {10,S} {18,S}
9  C u0 p0 c0 {7,S} {11,S} {19,S} {20,S}
10 C u1 p0 c0 {3,S} {8,S} {21,S}
11 C u0 p0 c0 {2,S} {4,D} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(16854)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {17,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {7,D}
5  C u0 p0 c0 {3,S} {6,D} {12,S}
6  C u0 p0 c0 {4,S} {5,D} {15,S}
7  C u0 p0 c0 {4,D} {8,S} {14,S}
8  C u0 p0 c0 {7,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {8,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(20232)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {11,S} {22,S}
3  O u0 p2 c0 {1,S} {23,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {11,D} {20,S}
10 C u0 p0 c0 {1,S} {4,D} {8,S}
11 C u0 p0 c0 {2,S} {9,D} {21,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {2,S}
23 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(20813)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {11,S} {21,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {8,S} {10,S} {16,S} {17,S}
8  C u1 p0 c0 {5,S} {7,S} {18,S}
9  C u0 p0 c0 {6,S} {11,D} {19,S}
10 C u0 p0 c0 {1,S} {4,D} {7,S}
11 C u0 p0 c0 {2,S} {9,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(20094)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {7,S} {9,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {5,D} {7,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u0 p0 c0 {5,S} {8,D} {12,S}
7  C u1 p0 c0 {2,S} {4,S} {14,S}
8  C u0 p0 c0 {1,S} {6,D} {15,S}
9  C u0 p0 c0 {2,S} {10,D} {13,S}
10 C u0 p0 c0 {3,S} {9,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
        """),
)



species(
    label='S(5810)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
7  C u0 p0 c0 {1,D} {2,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(21071)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {8,S} {20,S}
4  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {18,S}
9  C u1 p0 c0 {6,S} {8,S} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(21064)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {9,S} {20,S}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {16,S}
7  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {2,S} {7,S} {17,S} {18,S}
9  C u1 p0 c0 {3,S} {6,S} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(19487)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {9,S} {17,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {9,D} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {8,D} {14,S}
8  C u0 p0 c0 {1,S} {7,D} {16,S}
9  C u0 p0 c0 {2,S} {5,D} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(20106)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {5,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {7,S} {12,S}
7  C u0 p0 c0 {6,S} {10,D} {13,S}
8  C u0 p0 c0 {1,S} {9,D} {14,S}
9  C u0 p0 c0 {2,S} {8,D} {15,S}
10 C u0 p0 c0 {3,S} {7,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21215)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {8,S} {9,D}
7  C u0 p0 c0 {4,S} {8,D} {13,S}
8  C u0 p0 c0 {6,S} {7,D} {14,S}
9  C u0 p0 c0 {2,S} {6,D} {15,S}
10 C u0 p0 c0 {3,D} {5,S} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21213)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {8,S}
2  O u0 p2 c0 {6,S} {9,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {13,S}
8  C u0 p0 c0 {1,S} {6,D} {14,S}
9  C u0 p0 c0 {2,S} {7,D} {15,S}
10 C u0 p0 c0 {3,D} {5,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21190)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {11,S}
2  O u0 p2 c0 {10,S} {21,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {7,S} {16,S} {17,S}
9  C u0 p0 c0 {5,S} {11,S} {18,S} {19,S}
10 C u1 p0 c0 {2,S} {6,S} {20,S}
11 C u0 p0 c0 {1,S} {4,D} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C4H2(86)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,T}
2 C u0 p0 c0 {1,S} {4,T}
3 C u0 p0 c0 {1,T} {5,S}
4 C u0 p0 c0 {2,T} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
        """),
)


species(
    label='nC4H3(87)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,T}
3 C u1 p0 c0 {1,D} {6,S}
4 C u0 p0 c0 {2,T} {7,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='iC4H3(88)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,D} {5,S} {6,S}
2 C u1 p0 c0 {1,D} {3,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {7,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(22782)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u1 p2 c0 {10,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u1 p0 c0 {4,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {8,D} {9,S}
7  C u0 p0 c0 {4,S} {10,D} {12,S}
8  C u0 p0 c0 {5,S} {6,D} {14,S}
9  C u0 p0 c0 {3,D} {6,S} {15,S}
10 C u0 p0 c0 {2,S} {7,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='nC4H5(91)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,D} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 C u1 p0 c0 {2,D} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H7(94)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u1 p0 c0 {1,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(23384)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {2,S} {8,D} {16,S}
10 C u0 p0 c0 {3,D} {4,S} {15,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(20164)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {6,S} {22,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {7,S} {13,S}
7  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8  C u0 p0 c0 {5,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {7,S} {11,S} {20,S} {21,S}
10 C u0 p0 c0 {8,S} {11,S} {14,S} {15,S}
11 C u0 p0 c0 {9,S} {10,S} {16,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(23383)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {9,D}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {13,S}
7  C u0 p0 c0 {4,S} {8,D} {14,S}
8  C u0 p0 c0 {1,S} {7,D} {10,S}
9  C u0 p0 c0 {2,D} {5,S} {15,S}
10 C u0 p0 c0 {3,D} {8,S} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21631)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {20,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(23699)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
9  C u1 p0 c0 {7,S} {8,S} {18,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(22863)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {9,S}
2  O u0 p2 c0 {11,S} {21,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {11,D} {17,S}
9  C u0 p0 c0 {1,S} {4,D} {6,S}
10 C u1 p0 c0 {7,S} {18,S} {19,S}
11 C u0 p0 c0 {2,S} {8,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(23382)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {14,S}
7  C u0 p0 c0 {1,S} {6,D} {10,S}
8  C u0 p0 c0 {4,S} {9,D} {13,S}
9  C u0 p0 c0 {2,S} {8,D} {15,S}
10 C u0 p0 c0 {3,D} {7,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21148)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {10,S}
2  O u0 p2 c0 {1,S} {22,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {7,S} {11,D} {20,S}
10 C u0 p0 c0 {1,S} {3,D} {8,S}
11 C u0 p0 c0 {4,S} {9,D} {21,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C4H612(97)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C6H5(99)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {8,S}
2  C u0 p0 c0 {1,B} {4,B} {7,S}
3  C u0 p0 c0 {1,B} {5,B} {9,S}
4  C u0 p0 c0 {2,B} {6,B} {10,S}
5  C u0 p0 c0 {3,B} {6,B} {11,S}
6  C u1 p0 c0 {4,B} {5,B}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C7H7(101)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {7,S}
2  C u0 p0 c0 {1,B} {4,B} {8,S}
3  C u0 p0 c0 {1,B} {6,B} {12,S}
4  C u0 p0 c0 {2,B} {5,B} {9,S}
5  C u0 p0 c0 {4,B} {6,B} {10,S}
6  C u0 p0 c0 {3,B} {5,B} {11,S}
7  C u1 p0 c0 {1,S} {13,S} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H8(102)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,B} {4,B}
3  C u0 p0 c0 {2,B} {5,B} {11,S}
4  C u0 p0 c0 {2,B} {7,B} {15,S}
5  C u0 p0 c0 {3,B} {6,B} {12,S}
6  C u0 p0 c0 {5,B} {7,B} {13,S}
7  C u0 p0 c0 {4,B} {6,B} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C6H5O(138)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,B} {7,B}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {8,S}
6  C u0 p0 c0 {5,B} {7,B} {11,S}
7  C u0 p0 c0 {2,B} {6,B} {12,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(583)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,B} {5,B}
4  C u0 p0 c0 {3,B} {6,B} {9,S}
5  C u0 p0 c0 {3,B} {8,B} {13,S}
6  C u0 p0 c0 {4,B} {7,B} {10,S}
7  C u0 p0 c0 {6,B} {8,B} {11,S}
8  C u0 p0 c0 {5,B} {7,B} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(24214)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {4,S} {7,D} {11,S}
6  C u0 p0 c0 {4,S} {8,D} {12,S}
7  C u0 p0 c0 {5,D} {9,S} {13,S}
8  C u0 p0 c0 {6,D} {9,S} {14,S}
9  C u0 p0 c0 {2,D} {7,S} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(24003)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
4  C u1 p0 c0 {3,S} {7,S} {10,S}
5  C u0 p0 c0 {2,S} {7,D} {8,S}
6  C u0 p0 c0 {3,S} {8,D} {11,S}
7  C u0 p0 c0 {4,S} {5,D} {12,S}
8  C u0 p0 c0 {5,S} {6,D} {13,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(24213)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {2,D} {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {14,S}
9  C u0 p0 c0 {7,S} {8,D} {13,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(23997)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4  C u1 p0 c0 {3,S} {6,S} {9,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {5,D} {8,S} {13,S}
8  C u0 p0 c0 {6,D} {7,S} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(24545)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {5,S} {10,D} {15,S}
9  C u0 p0 c0 {6,S} {7,D} {12,S}
10 C u0 p0 c0 {6,S} {8,D} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(588)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {4,D} {8,S} {12,S}
7  C u0 p0 c0 {5,D} {8,S} {13,S}
8  C u0 p0 c0 {2,D} {6,S} {7,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(24544)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {7,S} {12,S}
7  C u0 p0 c0 {6,S} {8,D} {14,S}
8  C u0 p0 c0 {2,S} {7,D} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {15,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(24044)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {15,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {3,D} {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {9,D} {14,S}
9  C u0 p0 c0 {7,S} {8,D} {13,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(24609)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {14,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {4,D} {8,S} {12,S}
7  C u0 p0 c0 {5,D} {8,S} {13,S}
8  C u0 p0 c0 {2,D} {6,S} {7,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(24558)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {5,B} {8,B}
4  C u0 p0 c0 {5,B} {6,B} {9,S}
5  C u0 p0 c0 {3,B} {4,B} {10,S}
6  C u0 p0 c0 {2,S} {4,B} {7,B}
7  C u0 p0 c0 {6,B} {8,B} {11,S}
8  C u0 p0 c0 {3,B} {7,B} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(24566)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,D} {10,S}
5  C u0 p0 c0 {3,S} {8,D} {11,S}
6  C u0 p0 c0 {4,D} {7,S} {12,S}
7  C u0 p0 c0 {2,D} {6,S} {8,S}
8  C u1 p0 c0 {5,D} {7,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(24561)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {8,D} {12,S}
6  C u0 p0 c0 {3,S} {7,D} {11,S}
7  C u0 p0 c0 {6,D} {8,S} {13,S}
8  C u0 p0 c0 {2,S} {5,D} {7,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(24923)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {8,D}
6  C u1 p0 c0 {3,S} {4,S} {11,S}
7  C u0 p0 c0 {2,D} {3,S} {5,S}
8  C u0 p0 c0 {4,S} {5,D} {12,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(25021)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
8  C u0 p0 c0 {6,S} {10,D} {14,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {15,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(24565)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {6,S} {11,S}
6  C u0 p0 c0 {2,D} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {12,S}
8  C u1 p0 c0 {3,S} {7,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(25016)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {9,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {2,D} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {14,S}
9  C u0 p0 c0 {3,S} {7,D} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(25096)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {10,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u1 p0 c0 {6,S} {10,S} {14,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {4,S} {8,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(25092)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {3,D} {5,S} {9,S}
8  C u0 p0 c0 {4,D} {6,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {8,S} {9,D} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(141)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,D}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,D} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {6,S} {10,S}
6  C u0 p0 c0 {2,D} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {7,D} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(24928)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u1 p0 c0 {3,S} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {4,S} {7,S} {8,D}
7  C u0 p0 c0 {5,D} {6,S} {12,S}
8  C u0 p0 c0 {2,D} {6,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='iC4H7(103)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
        """),
)


species(
    label='cC3H4(104)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 C u0 p0 c0 {1,S} {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='iC4H8(106)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
        """),
)



species(
    label='S(24603)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {14,S}
3  O u0 p2 c0 {1,S} {15,S}
4  C u0 p0 c0 {1,S} {5,B} {6,B}
5  C u0 p0 c0 {2,S} {4,B} {7,B}
6  C u0 p0 c0 {4,B} {8,B} {10,S}
7  C u0 p0 c0 {5,B} {9,B} {13,S}
8  C u0 p0 c0 {6,B} {9,B} {11,S}
9  C u0 p0 c0 {7,B} {8,B} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(25709)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {15,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u1 p2 c0 {2,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {9,S} {10,D}
8  C u0 p0 c0 {6,S} {9,D} {13,S}
9  C u0 p0 c0 {7,S} {8,D} {14,S}
10 C u0 p0 c0 {4,D} {7,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(25576)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {13,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,B} {5,B}
4  C u0 p0 c0 {2,S} {3,B} {6,B}
5  C u0 p0 c0 {3,B} {8,B} {12,S}
6  C u0 p0 c0 {4,B} {7,B} {9,S}
7  C u0 p0 c0 {6,B} {8,B} {10,S}
8  C u0 p0 c0 {5,B} {7,B} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(25710)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {15,S}
2  O u0 p2 c0 {4,S} {10,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {7,D} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u0 p0 c0 {6,S} {8,D} {14,S}
10 C u0 p0 c0 {2,S} {3,D} {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C3H6O(108)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
4  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(25532)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {5,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {3,D} {5,S} {9,S}
8  C u0 p0 c0 {4,D} {6,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {8,S} {9,D} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(25550)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {13,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {7,D}
5  C u0 p0 c0 {6,D} {7,S} {8,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u0 p0 c0 {4,D} {5,S} {12,S}
8  C u1 p0 c0 {2,D} {5,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(25542)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {12,S}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {1,S} {5,D} {6,S}
4  C u0 p0 c0 {5,S} {7,S} {8,D}
5  C u0 p0 c0 {3,D} {4,S} {9,S}
6  C u0 p0 c0 {3,S} {7,D} {10,S}
7  C u0 p0 c0 {4,S} {6,D} {11,S}
8  C u0 p0 c0 {2,D} {4,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C3H3O(109)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {2,D} {6,S} {7,S}
4 C u1 p0 c0 {1,D} {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(25873)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {8,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {9,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,D} {7,S} {12,S}
6  C u0 p0 c0 {5,D} {8,S} {11,S}
7  C u0 p0 c0 {2,D} {5,S} {10,S}
8  C u0 p0 c0 {1,S} {6,S} {9,D}
9  C u0 p0 c0 {3,S} {8,D} {14,S}
10 C u0 p0 c0 {4,D} {7,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(24564)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {13,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,B} {5,B}
4  C u0 p0 c0 {3,B} {6,B} {9,S}
5  C u0 p0 c0 {3,B} {7,B} {10,S}
6  C u0 p0 c0 {4,B} {8,B} {11,S}
7  C u0 p0 c0 {5,B} {8,B} {12,S}
8  C u0 p0 c0 {2,S} {6,B} {7,B}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(25881)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {14,S}
2  O u1 p2 c0 {5,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {1,S} {6,D} {9,S}
8  C u0 p0 c0 {4,D} {5,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {12,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(26372)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {14,S}
3  O u1 p2 c0 {1,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {5,S} {9,S} {10,D}
7  C u1 p0 c0 {5,S} {8,S} {12,S}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {13,S}
10 C u0 p0 c0 {4,D} {6,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(25994)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {15,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {11,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {6,D} {10,S} {13,S}
9  C u0 p0 c0 {7,D} {10,S} {14,S}
10 C u0 p0 c0 {3,D} {8,S} {9,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(26654)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {8,D}
3  O u0 p2 c0 {9,D}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {1,S} {6,D} {9,S}
8  C u0 p0 c0 {2,D} {5,S} {10,S}
9  C u0 p0 c0 {3,D} {7,S} {14,S}
10 C u0 p0 c0 {4,D} {8,S} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(26715)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {9,S} {14,S}
3  O u0 p2 c0 {7,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u0 p0 c0 {4,D} {5,S} {9,S}
9  C u0 p0 c0 {2,S} {8,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(25129)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {15,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {2,S} {6,D} {9,S}
8  C u0 p0 c0 {3,D} {5,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(26369)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,S} {14,S}
3  O u1 p2 c0 {1,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u1 p0 c0 {5,S} {8,S} {11,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  C u0 p0 c0 {7,D} {8,S} {13,S}
10 C u0 p0 c0 {4,D} {8,D}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(26365)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {15,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {6,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,D} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {9,D} {11,S}
8  C u0 p0 c0 {6,S} {10,D} {12,S}
9  C u0 p0 c0 {7,D} {10,S} {14,S}
10 C u0 p0 c0 {8,D} {9,S} {13,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(26367)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {15,S}
3  O u0 p2 c0 {6,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {3,D} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {12,S}
9  C u0 p0 c0 {7,D} {10,S} {14,S}
10 C u0 p0 c0 {8,D} {9,S} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(26371)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u1 p2 c0 {1,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {9,D}
7  C u1 p0 c0 {5,S} {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {6,D} {8,S} {13,S}
10 C u0 p0 c0 {4,D} {8,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='H2C4O(111)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {5,D}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {2,D} {4,D}
4 C u0 p0 c0 {3,D} {5,D}
5 C u0 p0 c0 {1,D} {4,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C6H2(112)',
    reactive=True,
    structure=adjacencyList(
        """
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 C u0 p0 c0 {1,S} {5,T}
4 C u0 p0 c0 {2,S} {6,T}
5 C u0 p0 c0 {3,T} {7,S}
6 C u0 p0 c0 {4,T} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H3(113)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {4,D} {7,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {5,S}
4 C u1 p0 c0 {1,D} {8,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H4(114)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,D} {3,S} {7,S}
2  C u0 p0 c0 {1,D} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,T}
4  C u0 p0 c0 {3,T} {5,S}
5  C u0 p0 c0 {4,S} {6,T}
6  C u0 p0 c0 {5,T} {10,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H4(115)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,D} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u0 p0 c0 {1,D} {6,S} {10,S}
4  C u0 p0 c0 {2,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,T}
6  C u0 p0 c0 {3,S} {5,T}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C4H5O(116)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,D} {6,S}
3  C u0 p0 c0 {2,S} {5,D} {7,S}
4  C u0 p0 c0 {2,D} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,D} {10,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H5(117)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2 C u0 p0 c0 {4,D} {8,S} {9,S}
3 C u1 p0 c0 {1,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C4H5O(119)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u1 p0 c0 {1,D} {4,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(25788)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u1 p2 c0 {1,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {5,S} {9,D}
7  C u0 p0 c0 {5,S} {8,S} {10,D}
8  C u1 p0 c0 {7,S} {9,S} {13,S}
9  C u0 p0 c0 {6,D} {8,S} {12,S}
10 C u0 p0 c0 {4,D} {7,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(26366)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {15,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {2,S} {6,D} {9,S}
8  C u0 p0 c0 {5,S} {10,D} {12,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C4H6O(120)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {2,S} {4,D} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(27339)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {8,D} {9,S} {10,S}
7  C u0 p0 c0 {3,S} {5,S} {9,D}
8  C u0 p0 c0 {5,S} {6,D} {12,S}
9  C u0 p0 c0 {6,S} {7,D} {13,S}
10 C u0 p0 c0 {2,S} {4,D} {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(26686)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {9,D}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,D} {7,S} {12,S}
6  C u0 p0 c0 {5,D} {8,S} {11,S}
7  C u0 p0 c0 {2,D} {5,S} {10,S}
8  C u0 p0 c0 {1,D} {6,S} {9,S}
9  C u0 p0 c0 {3,D} {8,S} {14,S}
10 C u0 p0 c0 {4,D} {7,S} {13,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(26660)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {10,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {6,D}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
6  C u0 p0 c0 {3,D} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {13,S}
8  C u0 p0 c0 {7,D} {9,S} {12,S}
9  C u0 p0 c0 {4,S} {8,S} {10,D}
10 C u0 p0 c0 {1,S} {9,D} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(25852)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {6,S} {15,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {8,S} {9,S} {11,S}
7  C u0 p0 c0 {5,S} {9,D} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {13,S}
9  C u0 p0 c0 {6,S} {7,D} {14,S}
10 C u0 p0 c0 {2,S} {4,D} {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C4H6O(121)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {5,D} {9,S}
5  C u0 p0 c0 {4,D} {10,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H6(122)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,T}
4  C u0 p0 c0 {2,S} {3,T}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C4H6O(123)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {5,D} {10,S}
5  C u0 p0 c0 {1,S} {4,D} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H6O(124)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {1,D} {4,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H4O(125)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {5,S}
2 C u0 p0 c0 {3,S} {5,D} {6,S}
3 C u0 p0 c0 {2,S} {4,D} {7,S}
4 C u0 p0 c0 {1,S} {3,D} {8,S}
5 C u0 p0 c0 {1,S} {2,D} {9,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C4H82(128)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(27643)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {5,S} {14,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u0 p0 c0 {8,D} {9,S} {10,S}
8  C u0 p0 c0 {5,S} {7,D} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p0 c0 {2,S} {4,D} {7,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(27644)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {8,S} {14,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {9,D} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {13,S}
10 C u0 p0 c0 {2,S} {4,D} {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(582)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,D}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {1,D} {4,S} {8,S}
4  C u0 p0 c0 {2,D} {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {9,S}
6  C u0 p0 c0 {5,D} {7,S} {10,S}
7  C u0 p0 c0 {6,S} {8,D} {11,S}
8  C u0 p0 c0 {3,S} {7,D} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(25947)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {13,S}
2  O u0 p2 c0 {4,S} {14,S}
3  C u0 p0 c0 {1,S} {4,B} {5,B}
4  C u0 p0 c0 {2,S} {3,B} {6,B}
5  C u0 p0 c0 {3,B} {7,B} {9,S}
6  C u0 p0 c0 {4,B} {8,B} {11,S}
7  C u0 p0 c0 {5,B} {8,B} {10,S}
8  C u0 p0 c0 {6,B} {7,B} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='iC4H9(129)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C4H10(131)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
        """),
)


species(
    label='tC4H9(132)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H7(133)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,B} {4,B}
3  C u0 p0 c0 {2,B} {5,B} {11,S}
4  C u0 p0 c0 {2,B} {6,B} {12,S}
5  C u0 p0 c0 {3,B} {7,B} {13,S}
6  C u0 p0 c0 {4,B} {7,B} {14,S}
7  C u1 p0 c0 {5,B} {6,B}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H7O(134)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,B} {5,B}
4  C u0 p0 c0 {3,B} {6,B} {12,S}
5  C u0 p0 c0 {3,B} {8,B} {15,S}
6  C u0 p0 c0 {4,B} {7,B} {13,S}
7  C u0 p0 c0 {6,B} {8,B} {14,S}
8  C u0 p0 c0 {1,S} {5,B} {7,B}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(27807)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {9,S} {14,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {7,D} {12,S}
7  C u0 p0 c0 {1,S} {6,D} {9,S}
8  C u0 p0 c0 {5,S} {9,D} {13,S}
9  C u0 p0 c0 {3,S} {7,S} {8,D}
10 C u0 p0 c0 {2,S} {4,D} {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(27812)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {8,S} {14,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {10,S}
4  O u1 p2 c0 {10,S}
5  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {10,D}
7  C u0 p0 c0 {2,D} {5,S} {8,S}
8  C u0 p0 c0 {1,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {13,S}
10 C u0 p0 c0 {3,S} {4,S} {6,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(27789)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {14,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {10,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {7,S} {9,D} {10,S}
7  C u1 p0 c0 {5,S} {6,S} {12,S}
8  C u0 p0 c0 {2,D} {5,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {13,S}
10 C u0 p0 c0 {3,S} {4,D} {6,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(27314)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u1 p0 c0 {5,S} {9,S} {12,S}
7  C u0 p0 c0 {2,D} {5,S} {8,S}
8  C u0 p0 c0 {3,D} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(28121)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {11,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {2,D} {3,S} {6,S}
5  C u1 p0 c0 {3,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u1 p0 c0 {5,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28063)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {8,S} {14,S}
2  O u1 p2 c0 {9,S}
3  O u1 p2 c0 {10,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {8,D} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {1,S} {6,D} {9,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 C u0 p0 c0 {3,S} {4,D} {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28079)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {11,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {5,D} {6,S}
5  C u0 p0 c0 {2,S} {3,S} {4,D}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u1 p0 c0 {3,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28120)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {11,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u0 p0 c0 {2,D} {3,S} {6,S}
5  C u0 p0 c0 {3,S} {7,D} {9,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28324)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {11,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {7,D}
6  C u0 p0 c0 {2,D} {3,S} {4,S}
7  C u0 p0 c0 {3,S} {5,D} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28266)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {11,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
4  C u1 p0 c0 {3,S} {6,S} {9,S}
5  C u0 p0 c0 {2,D} {3,S} {7,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u1 p0 c0 {5,S} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28321)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {11,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {3,S} {5,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {3,S} {7,D} {10,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28658)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {7,S} {13,S}
3  O u1 p2 c0 {8,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {9,D}
7  C u0 p0 c0 {2,S} {8,D} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {7,D}
9  C u0 p0 c0 {6,D} {7,S} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(29234)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u0 p0 c0 {3,S} {7,D} {9,S}
7  C u0 p0 c0 {5,S} {6,D} {10,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(29272)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {8,S} {13,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {3,S} {7,D} {9,S}
9  C u0 p0 c0 {4,D} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(29273)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {8,S} {13,S}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u1 p0 c0 {5,S} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {6,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(29265)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {8,D}
6  C u0 p0 c0 {2,D} {4,S} {5,S}
7  C u0 p0 c0 {3,D} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(29377)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {8,S} {13,S}
2  O u1 p2 c0 {5,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {7,D} {9,S}
9  C u0 p0 c0 {4,D} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28268)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {11,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {7,S} {9,S}
5  C u0 p0 c0 {2,S} {3,D} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {10,S}
7  C u0 p0 c0 {4,S} {6,D} {8,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(29232)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {11,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(28260)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {11,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {2,D} {3,S} {6,S}
5  C u0 p0 c0 {3,D} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {7,D} {10,S}
7  C u0 p0 c0 {5,S} {6,D} {8,S}
8  H u0 p0 c0 {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H6O(135)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,B} {4,B} {8,S}
3  C u0 p0 c0 {2,B} {5,B} {9,S}
4  C u0 p0 c0 {2,B} {7,B} {13,S}
5  C u0 p0 c0 {3,B} {6,B} {10,S}
6  C u0 p0 c0 {5,B} {7,B} {11,S}
7  C u0 p0 c0 {4,B} {6,B} {12,S}
8  C u0 p0 c0 {1,D} {2,S} {14,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(24207)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {1,S} {4,D} {6,S}
6  C u1 p0 c0 {5,S} {11,S} {12,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C7H8O(136)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,B} {5,B}
4  C u0 p0 c0 {3,B} {6,B} {11,S}
5  C u0 p0 c0 {3,B} {8,B} {15,S}
6  C u0 p0 c0 {4,B} {7,B} {12,S}
7  C u0 p0 c0 {6,B} {8,B} {13,S}
8  C u0 p0 c0 {5,B} {7,B} {14,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C6H6O(137)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {13,S}
2  C u0 p0 c0 {1,S} {3,B} {4,B}
3  C u0 p0 c0 {2,B} {5,B} {8,S}
4  C u0 p0 c0 {2,B} {7,B} {12,S}
5  C u0 p0 c0 {3,B} {6,B} {9,S}
6  C u0 p0 c0 {5,B} {7,B} {10,S}
7  C u0 p0 c0 {4,B} {6,B} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H8O(139)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {16,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,B} {6,B}
4  C u0 p0 c0 {1,S} {5,B} {7,B}
5  C u0 p0 c0 {3,B} {4,B} {15,S}
6  C u0 p0 c0 {3,B} {8,B} {12,S}
7  C u0 p0 c0 {4,B} {8,B} {14,S}
8  C u0 p0 c0 {6,B} {7,B} {13,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H5O(140)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,B} {4,B} {8,S}
3  C u0 p0 c0 {2,B} {5,B} {9,S}
4  C u0 p0 c0 {2,B} {7,B} {13,S}
5  C u0 p0 c0 {3,B} {6,B} {10,S}
6  C u0 p0 c0 {5,B} {7,B} {11,S}
7  C u0 p0 c0 {4,B} {6,B} {12,S}
8  C u1 p0 c0 {1,D} {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C5H5O(144)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {2,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u0 p0 c0 {2,S} {6,D} {9,S}
5  C u0 p0 c0 {3,D} {6,S} {10,S}
6  C u0 p0 c0 {4,D} {5,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C5H5O(145)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {5,S} {7,S}
4  C u0 p0 c0 {2,D} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {6,D} {8,S}
6  C u0 p0 c0 {4,S} {5,D} {9,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(29864)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {13,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {2,D} {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u0 p0 c0 {3,S} {4,D} {5,S}
9  C u1 p0 c0 {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(29422)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {6,S} {12,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {8,D}
6  C u0 p0 c0 {2,S} {7,D} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {6,D}
8  C u0 p0 c0 {5,D} {6,S} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(29421)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {13,S}
2  O u0 p2 c0 {6,S} {12,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {11,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(30089)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {7,S} {13,S}
3  O u0 p2 c0 {6,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {3,D} {5,S} {7,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {1,S} {4,D} {5,S}
9  C u0 p0 c0 {1,S} {7,D} {12,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(17259)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {10,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {1,S} {4,S} {6,D}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {7,S} {8,S}
6  C u1 p0 c0 {3,D} {9,S}
7  H u0 p0 c0 {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(30129)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u0 p2 c0 {3,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H5O(592)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,D}
2  C u0 p0 c0 {3,D} {4,S} {7,S}
3  C u0 p0 c0 {2,D} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {6,D} {9,S}
5  C u0 p0 c0 {1,D} {3,S} {10,S}
6  C u1 p0 c0 {4,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(28329)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {2,S} {4,D} {7,S}
7  C u0 p0 c0 {5,D} {6,S} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(30866)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {6,S} {11,D}
8  C u1 p0 c0 {4,S} {10,S} {18,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {19,S}
11 C u0 p0 c0 {3,S} {7,D} {20,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='C5H6O(146)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {5,D} {8,S}
4  C u0 p0 c0 {2,S} {6,D} {9,S}
5  C u0 p0 c0 {3,D} {6,S} {10,S}
6  C u0 p0 c0 {4,D} {5,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C4H5(147)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {4,S} {7,S}
3 C u0 p0 c0 {1,S} {4,D} {8,S}
4 C u0 p0 c0 {2,S} {3,D} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C5H9(148)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u1 p0 c0 {2,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(27152)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {7,S} {14,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {12,S}
9  C u0 p0 c0 {1,S} {8,D} {13,S}
10 C u1 p0 c0 {4,D} {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='cC5H9(151)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {4,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='cC5H8(152)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {2,S} {4,D} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(1097)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,D} {15,S}
5  C u0 p0 c0 {1,S} {6,D} {14,S}
6  C u0 p0 c0 {3,S} {5,D} {16,S}
7  C u0 p0 c0 {3,S} {4,D} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H9(153)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {2,S} {11,S}
4  C u0 p0 c0 {1,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H8(154)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,S} {5,D} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H9(155)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H9(156)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  C u1 p0 c0 {4,S} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)




species(
    label='S(30985)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  C u0 p0 c0 {3,D} {4,S} {7,S}
3  C u0 p0 c0 {2,D} {5,S} {8,S}
4  C u0 p0 c0 {2,S} {6,D} {9,S}
5  C u1 p0 c0 {1,S} {3,S} {10,S}
6  C u0 p0 c0 {1,S} {4,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(16923)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,S} {5,S}
4 C u1 p0 c0 {2,D} {3,S}
5 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H9(1019)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {13,S}
4  C u0 p0 c0 {2,S} {3,D} {14,S}
5  C u1 p0 c0 {1,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {7,D} {15,S}
7  C u0 p0 c0 {5,S} {6,D} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
        """),
)



species(
    label='S(7919)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {12,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {2,S} {6,D} {9,S}
5  C u0 p0 c0 {3,D} {6,S} {10,S}
6  C u0 p0 c0 {4,D} {5,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(31426)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
4  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {6,D} {11,S}
6  C u0 p0 c0 {4,S} {5,D} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C6H11(161)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {13,S}
5  C u1 p0 c0 {3,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H10(162)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {5,S} {14,S}
5  C u0 p0 c0 {4,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(164)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {14,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {16,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(165)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {14,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(166)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {6,S} {15,S}
6  C u1 p0 c0 {5,S} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)



species(
    label='S(23905)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {21,S}
3  O u0 p2 c0 {9,D}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {10,D} {18,S}
9  C u0 p0 c0 {1,S} {3,D} {7,S}
10 C u0 p0 c0 {8,D} {11,S} {19,S}
11 C u0 p0 c0 {4,D} {10,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(23820)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {18,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(23795)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(23807)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {1,S} {19,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {11,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {9,D} {16,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(31869)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {4,S} {14,S}
6  C u0 p0 c0 {2,S} {7,D} {16,S}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(32442)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {8,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u0 p0 c0 {6,S} {8,D} {12,S}
8  C u0 p0 c0 {1,S} {7,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1123)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,D} {18,S}
7  C u0 p0 c0 {4,S} {8,D} {17,S}
8  C u0 p0 c0 {5,S} {7,D} {15,S}
9  C u0 p0 c0 {5,S} {6,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(32443)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {8,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,D} {10,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {1,S} {5,D} {12,S}
8  C u0 p0 c0 {1,S} {6,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H13(171)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {16,S}
6  C u1 p0 c0 {4,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)




species(
    label='S(30079)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u1 p0 c0 {3,S} {6,S} {10,S}
5  C u0 p0 c0 {2,S} {6,D} {9,S}
6  C u0 p0 c0 {4,S} {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(9973)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {12,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {2,S} {3,D} {9,S}
5  C u0 p0 c0 {2,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(31888)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
6  C u1 p0 c0 {3,S} {4,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(9915)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {13,S}
6  C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
7  C u1 p0 c0 {3,S} {6,S} {16,S}
8  C u0 p0 c0 {5,S} {9,D} {17,S}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(33011)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
7  C u1 p0 c0 {3,S} {4,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H9(1525)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
5  C u1 p0 c0 {1,S} {3,S} {14,S}
6  C u0 p0 c0 {2,S} {7,D} {15,S}
7  C u0 p0 c0 {4,S} {6,D} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H9(1939)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u0 p0 c0 {2,S} {3,D} {12,S}
6  C u0 p0 c0 {2,S} {4,D} {14,S}
7  C u1 p0 c0 {1,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C8H15(178)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,D} {19,S}
7  C u1 p0 c0 {5,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H16(179)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {8,D} {22,S}
8  C u0 p0 c0 {7,D} {23,S} {24,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H17(180)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {6,S} {24,S} {25,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H17(181)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
8  C u1 p0 c0 {5,S} {7,S} {25,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H17(182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {22,S} {23,S} {24,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u1 p0 c0 {4,S} {5,S} {25,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H17(183)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
3  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {2,S} {22,S} {23,S} {24,S}
8  C u1 p0 c0 {4,S} {5,S} {25,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(33085)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
4  C u0 p0 c0 {3,S} {5,D} {8,S}
5  C u0 p0 c0 {2,S} {4,D} {9,S}
6  C u1 p0 c0 {3,S} {10,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(33029)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
6  C u1 p0 c0 {5,S} {8,S} {12,S}
7  C u0 p0 c0 {4,S} {8,D} {11,S}
8  C u0 p0 c0 {6,S} {7,D} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(33030)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {8,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {5,S}
4  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  C u1 p0 c0 {1,S} {4,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2835)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
7  C u1 p0 c0 {3,S} {5,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {18,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(5985)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u1 p0 c0 {2,S} {6,S} {12,S}
8  C u0 p0 c0 {3,D} {4,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(33038)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
7  C u1 p0 c0 {4,S} {6,S} {16,S}
8  C u0 p0 c0 {5,S} {9,D} {18,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(33549)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {8,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {8,D} {15,S}
7  C u0 p0 c0 {3,S} {9,D} {17,S}
8  C u0 p0 c0 {4,S} {6,D} {16,S}
9  C u0 p0 c0 {4,S} {7,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(33087)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
6  C u1 p0 c0 {2,S} {3,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C8H18(184)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(33550)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {8,D} {12,S}
8  C u0 p0 c0 {2,S} {7,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C9H17(185)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {20,S} {21,S}
6  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,D} {22,S}
8  C u1 p0 c0 {6,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H18(186)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {9,D} {25,S}
9  C u0 p0 c0 {8,D} {26,S} {27,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H9(8388)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u1 p0 c0 {1,S} {5,S} {12,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u0 p0 c0 {4,D} {7,S} {16,S}
7  C u0 p0 c0 {5,D} {6,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(33322)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {2,S} {4,D} {15,S}
6  C u0 p0 c0 {1,S} {7,D} {16,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(34440)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
3  C u0 p0 c0 {4,S} {5,D} {7,S}
4  C u0 p0 c0 {3,S} {6,D} {8,S}
5  C u0 p0 c0 {1,S} {3,D} {15,S}
6  C u0 p0 c0 {2,S} {4,D} {16,S}
7  C u1 p0 c0 {3,S} {17,S} {18,S}
8  C u1 p0 c0 {4,S} {19,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C9H19(187)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {7,S} {27,S} {28,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H19(188)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {9,S} {25,S} {26,S} {27,S}
9  C u1 p0 c0 {6,S} {8,S} {28,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(2570)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
5  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {8,D} {16,S}
7  C u0 p0 c0 {3,S} {9,D} {17,S}
8  C u0 p0 c0 {5,S} {6,D} {15,S}
9  C u0 p0 c0 {5,S} {7,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(33910)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {7,S} {14,S}
6  C u0 p0 c0 {4,S} {8,D} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {2,S} {6,D} {16,S}
9  C u0 p0 c0 {7,D} {17,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(6403)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {8,S} {11,S}
7  C u0 p0 c0 {2,D} {4,S} {12,S}
8  C u0 p0 c0 {3,D} {6,S} {13,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(33965)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {3,S} {14,S} {15,S}
6  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
7  C u1 p0 c0 {4,S} {6,S} {16,S}
8  C u0 p0 c0 {3,S} {9,D} {17,S}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3379)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,D}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,D} {5,S} {7,S}
4  C u0 p0 c0 {3,D} {6,S} {8,S}
5  C u0 p0 c0 {1,D} {3,S} {9,S}
6  C u0 p0 c0 {2,D} {4,S} {10,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(34749)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {8,D} {11,S}
7  C u0 p0 c0 {2,D} {5,S} {12,S}
8  C u0 p0 c0 {3,S} {6,D} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(34760)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {13,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,D} {9,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,D} {11,S}
8  C u0 p0 c0 {3,D} {6,S} {12,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(34744)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,D} {6,S} {10,S}
5  C u0 p0 c0 {4,D} {7,S} {9,S}
6  C u0 p0 c0 {1,D} {4,S} {8,S}
7  C u0 p0 c0 {2,D} {5,S} {11,S}
8  C u0 p0 c0 {3,D} {6,S} {12,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(8001)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {7,D} {16,S}
6  C u0 p0 c0 {2,S} {8,D} {17,S}
7  C u0 p0 c0 {3,S} {5,D} {14,S}
8  C u0 p0 c0 {3,S} {6,D} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(35144)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {11,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {10,D} {18,S}
9  C u0 p0 c0 {6,S} {11,D} {17,S}
10 C u0 p0 c0 {7,S} {8,D} {19,S}
11 C u0 p0 c0 {2,S} {9,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(35143)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {10,D} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {17,S}
10 C u0 p0 c0 {2,S} {8,D} {18,S}
11 C u0 p0 c0 {9,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(34771)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {8,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {9,S} {12,S}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {13,S}
7  C u0 p0 c0 {3,S} {4,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
9  C u1 p0 c0 {4,S} {5,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(33332)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u1 p0 c0 {1,S} {3,S} {10,S}
6  C u0 p0 c0 {2,S} {4,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(321)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u1 p0 c0 {1,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(34706)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
6  C u0 p0 c0 {2,D} {3,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(33975)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {1,S} {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
8  C u1 p0 c0 {4,S} {6,S} {13,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C9H19(189)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {9,S} {20,S} {21,S}
6  C u0 p0 c0 {8,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {25,S} {26,S} {27,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {5,S} {6,S} {28,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H19(190)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {4,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {3,S} {25,S} {26,S} {27,S}
9  C u1 p0 c0 {5,S} {6,S} {28,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H19(191)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {3,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {25,S} {26,S} {27,S}
9  C u1 p0 c0 {5,S} {6,S} {28,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(28322)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,D} {3,S} {6,S}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {2,D} {4,S} {7,S}
7  C u0 p0 c0 {5,D} {6,S} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(34185)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u0 p0 c0 {2,D} {4,S} {9,S}
8  C u0 p0 c0 {6,D} {10,S} {15,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(27153)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {2,S} {7,D} {8,S}
7  C u0 p0 c0 {5,S} {6,D} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {1,S} {8,D} {13,S}
10 C u1 p0 c0 {4,D} {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(34187)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {12,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {9,S} {10,D}
7  C u0 p0 c0 {4,S} {9,D} {15,S}
8  C u0 p0 c0 {2,D} {4,S} {10,S}
9  C u0 p0 c0 {6,S} {7,D} {16,S}
10 C u0 p0 c0 {6,D} {8,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C9H20(192)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {5,S} {18,S} {19,S}
5  C u0 p0 c0 {4,S} {7,S} {20,S} {21,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(36048)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u0 p2 c0 {3,S} {8,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {3,S} {9,S} {15,S}
8  C u0 p0 c0 {4,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {6,S} {7,S} {18,S}
10 C u0 p0 c0 {5,S} {11,D} {19,S}
11 C u0 p0 c0 {8,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(36061)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u0 p2 c0 {2,S} {8,S}
4  O u0 p2 c0 {1,S} {10,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {17,S}
10 C u1 p0 c0 {4,S} {8,S} {18,S}
11 C u0 p0 c0 {9,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(36040)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {8,S}
4  O u0 p2 c0 {2,S} {11,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {6,S} {11,S} {15,S}
8  C u0 p0 c0 {3,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {18,S}
10 C u0 p0 c0 {8,S} {9,D} {19,S}
11 C u1 p0 c0 {4,S} {7,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(36062)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {8,S}
4  O u0 p2 c0 {2,S} {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {16,S}
9  C u1 p0 c0 {7,S} {8,S} {17,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {10,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(36256)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,D}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4  C u1 p0 c0 {3,S} {9,S} {10,S}
5  C u0 p0 c0 {1,D} {3,S} {8,S}
6  C u0 p0 c0 {2,D} {3,S} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(35475)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,D}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,D} {9,S}
5  C u0 p0 c0 {1,D} {3,S} {10,S}
6  C u0 p0 c0 {2,S} {4,D} {11,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(193)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
6  C u0 p0 c0 {5,S} {8,S} {23,S} {24,S}
7  C u0 p0 c0 {4,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {10,D} {25,S}
9  C u1 p0 c0 {7,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(194)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {9,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {10,D} {28,S}
10 C u0 p0 c0 {9,D} {29,S} {30,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(3333)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {8,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {8,D} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,D} {14,S}
8  C u0 p0 c0 {2,S} {5,D} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(34723)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
7  C u1 p0 c0 {5,S} {6,S} {15,S}
8  C u0 p0 c0 {3,S} {9,D} {16,S}
9  C u0 p0 c0 {8,D} {17,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(34003)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {5,D}
2 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {5,S}
5 C u1 p0 c0 {1,D} {4,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(7726)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
3  C u1 p0 c0 {2,S} {4,S} {6,S}
4  C u0 p0 c0 {1,D} {2,S} {3,S}
5  C u0 p0 c0 {2,S} {7,D} {10,S}
6  C u0 p0 c0 {3,S} {8,D} {12,S}
7  C u0 p0 c0 {5,D} {8,S} {13,S}
8  C u0 p0 c0 {6,D} {7,S} {11,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(3535)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {8,D} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,D} {14,S}
8  C u0 p0 c0 {2,S} {5,D} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(195)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {19,S} {20,S}
5  C u0 p0 c0 {4,S} {7,S} {21,S} {22,S}
6  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 C u1 p0 c0 {8,S} {30,S} {31,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(196)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {10,S} {28,S} {29,S} {30,S}
10 C u1 p0 c0 {7,S} {9,S} {31,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(197)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {10,S} {23,S} {24,S}
7  C u0 p0 c0 {9,S} {10,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
9  C u0 p0 c0 {7,S} {25,S} {26,S} {27,S}
10 C u1 p0 c0 {6,S} {7,S} {31,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(28348)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {12,S}
2  O u1 p2 c0 {4,S}
3  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u0 p0 c0 {3,S} {7,D} {10,S}
7  C u0 p0 c0 {5,S} {6,D} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(27170)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {6,S} {9,S}
3  O u0 p2 c0 {7,S} {14,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
6  C u0 p0 c0 {2,S} {7,S} {10,D}
7  C u0 p0 c0 {3,S} {6,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {12,S}
9  C u0 p0 c0 {2,S} {4,D} {5,S}
10 C u0 p0 c0 {1,S} {6,D} {13,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(36998)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {2,D} {7,S}
6  C u0 p0 c0 {4,S} {7,T}
7  C u0 p0 c0 {5,S} {6,T}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(29930)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {12,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {10,S}
6  C u0 p0 c0 {2,S} {4,D} {7,S}
7  C u0 p0 c0 {5,D} {6,S} {11,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(198)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {4,S} {17,S} {18,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {9,S} {19,S} {20,S}
5  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {5,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {4,S} {28,S} {29,S} {30,S}
10 C u1 p0 c0 {6,S} {7,S} {31,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(199)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
2  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
10 C u1 p0 c0 {6,S} {7,S} {31,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(200)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {19,S} {20,S}
5  C u0 p0 c0 {4,S} {6,S} {21,S} {22,S}
6  C u0 p0 c0 {5,S} {8,S} {23,S} {24,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 C u0 p0 c0 {8,S} {30,S} {31,S} {32,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(201)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
5  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {22,S} {23,S}
7  C u0 p0 c0 {6,S} {9,S} {26,S} {27,S}
8  C u0 p0 c0 {5,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {11,D} {28,S}
10 C u1 p0 c0 {8,S} {29,S} {30,S}
11 C u0 p0 c0 {9,D} {31,S} {32,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {11,S}
32 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(202)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {7,S} {22,S} {23,S}
6  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {24,S} {25,S}
8  C u0 p0 c0 {6,S} {10,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {11,D} {31,S}
11 C u0 p0 c0 {10,D} {32,S} {33,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(30452)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,B} {6,B} {10,S}
5  C u0 p0 c0 {4,B} {7,B} {11,S}
6  C u0 p0 c0 {4,B} {9,B} {15,S}
7  C u0 p0 c0 {5,B} {8,B} {12,S}
8  C u0 p0 c0 {7,B} {9,B} {13,S}
9  C u0 p0 c0 {6,B} {8,B} {14,S}
10 C u0 p0 c0 {1,S} {2,D} {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(30103)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {12,S}
2  O u0 p2 c0 {6,D}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,D} {4,S} {8,S}
7  C u0 p0 c0 {5,D} {8,S} {11,S}
8  C u0 p0 c0 {3,D} {6,S} {7,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(203)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {6,S} {22,S} {23,S}
6  C u0 p0 c0 {5,S} {8,S} {24,S} {25,S}
7  C u0 p0 c0 {1,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {11,S} {28,S} {29,S}
10 C u0 p0 c0 {8,S} {30,S} {31,S} {32,S}
11 C u1 p0 c0 {9,S} {33,S} {34,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(204)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {7,S} {22,S} {23,S}
6  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {24,S} {25,S}
8  C u0 p0 c0 {6,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {11,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {8,S} {10,S} {34,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(205)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
5  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {22,S} {23,S}
7  C u0 p0 c0 {5,S} {11,S} {26,S} {27,S}
8  C u0 p0 c0 {10,S} {11,S} {24,S} {25,S}
9  C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
10 C u0 p0 c0 {8,S} {28,S} {29,S} {30,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {9,S}
32 H u0 p0 c0 {9,S}
33 H u0 p0 c0 {9,S}
34 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(206)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
2  C u0 p0 c0 {1,S} {3,S} {18,S} {19,S}
3  C u0 p0 c0 {2,S} {5,S} {20,S} {21,S}
4  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {10,S} {22,S} {23,S}
6  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {11,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {6,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {5,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(207)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
2  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
3  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {10,S} {22,S} {23,S}
7  C u0 p0 c0 {3,S} {11,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(208)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
2  C u0 p0 c0 {4,S} {6,S} {20,S} {21,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {8,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {10,S} {22,S} {23,S}
7  C u0 p0 c0 {3,S} {11,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(209)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {6,S} {22,S} {23,S}
6  C u0 p0 c0 {5,S} {7,S} {24,S} {25,S}
7  C u0 p0 c0 {6,S} {9,S} {26,S} {27,S}
8  C u0 p0 c0 {1,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {11,S} {28,S} {29,S}
10 C u0 p0 c0 {8,S} {30,S} {31,S} {32,S}
11 C u0 p0 c0 {9,S} {33,S} {34,S} {35,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(210)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {7,S} {23,S} {24,S}
6  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {25,S} {26,S}
8  C u0 p0 c0 {7,S} {10,S} {29,S} {30,S}
9  C u0 p0 c0 {6,S} {11,S} {27,S} {28,S}
10 C u0 p0 c0 {8,S} {12,D} {31,S}
11 C u1 p0 c0 {9,S} {32,S} {33,S}
12 C u0 p0 c0 {10,D} {34,S} {35,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {12,S}
35 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(211)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {8,S} {25,S} {26,S}
7  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {27,S} {28,S}
9  C u0 p0 c0 {7,S} {11,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {9,S} {12,D} {34,S}
12 C u0 p0 c0 {11,D} {35,S} {36,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {12,S}
36 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(212)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {7,S} {25,S} {26,S}
7  C u0 p0 c0 {6,S} {9,S} {27,S} {28,S}
8  C u0 p0 c0 {1,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {11,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {12,S} {31,S} {32,S}
11 C u0 p0 c0 {9,S} {33,S} {34,S} {35,S}
12 C u1 p0 c0 {10,S} {36,S} {37,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {7,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {12,S}
37 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(213)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {8,S} {25,S} {26,S}
7  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {27,S} {28,S}
9  C u0 p0 c0 {7,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {12,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {9,S} {11,S} {37,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(214)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {7,S} {23,S} {24,S}
6  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {10,S} {25,S} {26,S}
8  C u0 p0 c0 {6,S} {12,S} {29,S} {30,S}
9  C u0 p0 c0 {11,S} {12,S} {27,S} {28,S}
10 C u0 p0 c0 {7,S} {34,S} {35,S} {36,S}
11 C u0 p0 c0 {9,S} {31,S} {32,S} {33,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {11,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {10,S}
35 H u0 p0 c0 {10,S}
36 H u0 p0 c0 {10,S}
37 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(215)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
2  C u0 p0 c0 {1,S} {3,S} {19,S} {20,S}
3  C u0 p0 c0 {2,S} {4,S} {21,S} {22,S}
4  C u0 p0 c0 {3,S} {6,S} {23,S} {24,S}
5  C u0 p0 c0 {1,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {11,S} {25,S} {26,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {12,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {7,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {6,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(216)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {19,S} {20,S}
2  C u0 p0 c0 {1,S} {3,S} {21,S} {22,S}
3  C u0 p0 c0 {2,S} {7,S} {23,S} {24,S}
4  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {9,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {11,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {12,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {7,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {2,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {3,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(217)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
2  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
3  C u0 p0 c0 {2,S} {7,S} {23,S} {24,S}
4  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {1,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {11,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {12,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {7,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {2,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {3,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(218)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {9,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {7,S} {25,S} {26,S}
7  C u0 p0 c0 {6,S} {8,S} {27,S} {28,S}
8  C u0 p0 c0 {7,S} {10,S} {29,S} {30,S}
9  C u0 p0 c0 {1,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {12,S} {31,S} {32,S}
11 C u0 p0 c0 {9,S} {33,S} {34,S} {35,S}
12 C u0 p0 c0 {10,S} {36,S} {37,S} {38,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {7,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {12,S}
37 H u0 p0 c0 {12,S}
38 H u0 p0 c0 {12,S}
        """),
)


species(
    label='S(219)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {13,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {12,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {19,S} {20,S}
5  C u0 p0 c0 {4,S} {6,S} {21,S} {22,S}
6  C u0 p0 c0 {5,S} {7,S} {23,S} {24,S}
7  C u0 p0 c0 {6,S} {8,S} {25,S} {26,S}
8  C u0 p0 c0 {7,S} {9,S} {27,S} {28,S}
9  C u0 p0 c0 {8,S} {10,S} {29,S} {30,S}
10 C u0 p0 c0 {9,S} {11,S} {31,S} {32,S}
11 C u0 p0 c0 {10,S} {13,S} {33,S} {34,S}
12 C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
13 C u0 p0 c0 {1,S} {11,S} {35,S} {36,S}
14 C u0 p0 c0 {12,S} {37,S} {38,S} {39,S}
15 H u0 p0 c0 {12,S}
16 H u0 p0 c0 {12,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {13,S}
36 H u0 p0 c0 {13,S}
37 H u0 p0 c0 {14,S}
38 H u0 p0 c0 {14,S}
39 H u0 p0 c0 {14,S}
        """),
)


species(
    label='S(220)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {12,S}
2  O u0 p2 c0 {1,S} {39,S}
3  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {19,S} {20,S}
5  C u0 p0 c0 {4,S} {6,S} {21,S} {22,S}
6  C u0 p0 c0 {5,S} {7,S} {23,S} {24,S}
7  C u0 p0 c0 {6,S} {8,S} {25,S} {26,S}
8  C u0 p0 c0 {7,S} {10,S} {27,S} {28,S}
9  C u0 p0 c0 {3,S} {13,S} {15,S} {16,S}
10 C u0 p0 c0 {8,S} {14,S} {29,S} {30,S}
11 C u0 p0 c0 {12,S} {14,S} {31,S} {32,S}
12 C u0 p0 c0 {1,S} {11,S} {33,S} {34,S}
13 C u0 p0 c0 {9,S} {35,S} {36,S} {37,S}
14 C u1 p0 c0 {10,S} {11,S} {38,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {11,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {12,S}
34 H u0 p0 c0 {12,S}
35 H u0 p0 c0 {13,S}
36 H u0 p0 c0 {13,S}
37 H u0 p0 c0 {13,S}
38 H u0 p0 c0 {14,S}
39 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(221)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {13,S}
2  O u0 p2 c0 {1,S} {39,S}
3  O u0 p2 c0 {15,D}
4  C u0 p0 c0 {5,S} {11,S} {18,S} {19,S}
5  C u0 p0 c0 {4,S} {6,S} {20,S} {21,S}
6  C u0 p0 c0 {5,S} {7,S} {22,S} {23,S}
7  C u0 p0 c0 {6,S} {8,S} {24,S} {25,S}
8  C u0 p0 c0 {7,S} {9,S} {26,S} {27,S}
9  C u0 p0 c0 {8,S} {10,S} {28,S} {29,S}
10 C u0 p0 c0 {9,S} {12,S} {30,S} {31,S}
11 C u0 p0 c0 {4,S} {14,S} {16,S} {17,S}
12 C u0 p0 c0 {10,S} {15,S} {32,S} {33,S}
13 C u0 p0 c0 {1,S} {15,S} {37,S} {38,S}
14 C u0 p0 c0 {11,S} {34,S} {35,S} {36,S}
15 C u0 p0 c0 {3,D} {12,S} {13,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {12,S}
33 H u0 p0 c0 {12,S}
34 H u0 p0 c0 {14,S}
35 H u0 p0 c0 {14,S}
36 H u0 p0 c0 {14,S}
37 H u0 p0 c0 {13,S}
38 H u0 p0 c0 {13,S}
39 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(20127)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {5,D}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u0 p0 c0 {2,S} {5,D} {7,S}
4 C u0 p0 c0 {2,D} {8,S} {9,S}
5 C u0 p0 c0 {1,D} {3,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(37603)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {1,D} {4,S} {5,S}
7  C u0 p0 c0 {4,S} {9,D} {14,S}
8  C u0 p0 c0 {2,D} {3,S} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(37605)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {1,D} {4,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {6,S} {7,D} {16,S}
9  C u0 p0 c0 {2,D} {3,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(36500)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {8,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {2,D} {3,S} {4,S}
6  C u0 p0 c0 {4,S} {8,D} {13,S}
7  C u0 p0 c0 {3,S} {9,D} {14,S}
8  C u0 p0 c0 {1,S} {6,D} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(222)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {21,S} {22,S}
4  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {9,S} {23,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(223)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {1,S} {8,S} {26,S} {27,S}
5  C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {7,S} {20,S} {21,S}
7  C u0 p0 c0 {6,S} {8,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {7,S} {24,S} {25,S}
9  C u0 p0 c0 {5,S} {10,S} {12,S} {13,S}
10 C u0 p0 c0 {9,S} {28,S} {29,S} {30,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {4,S}
27 H u0 p0 c0 {4,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C8H15(224)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H13(225)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
7  C u1 p0 c0 {1,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C6H11(226)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  C u1 p0 c0 {4,S} {5,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(37726)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {15,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {5,B} {9,S} {10,B}
5  C u0 p0 c0 {4,B} {6,B} {13,S}
6  C u0 p0 c0 {5,B} {7,B} {12,S}
7  C u0 p0 c0 {6,B} {8,B} {11,S}
8  C u0 p0 c0 {7,B} {10,B} {14,S}
9  C u0 p0 c0 {1,S} {3,D} {4,S}
10 C u1 p0 c0 {4,B} {8,B}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(37730)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,D} {10,S}
6  C u1 p0 c0 {4,S} {8,S} {12,S}
7  C u0 p0 c0 {5,D} {9,S} {15,S}
8  C u0 p0 c0 {6,S} {9,D} {13,S}
9  C u0 p0 c0 {7,S} {8,D} {14,S}
10 C u0 p0 c0 {2,S} {3,D} {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(37725)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u1 p0 c0 {4,S} {8,S} {11,S}
6  C u0 p0 c0 {4,S} {9,D} {12,S}
7  C u0 p0 c0 {2,S} {3,D} {4,S}
8  C u0 p0 c0 {5,S} {10,D} {13,S}
9  C u0 p0 c0 {6,D} {10,S} {15,S}
10 C u0 p0 c0 {8,D} {9,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(37775)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {10,S}
2  O u0 p2 c0 {1,S} {16,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,B} {6,B} {10,S}
5  C u0 p0 c0 {4,B} {7,B} {11,S}
6  C u0 p0 c0 {4,B} {8,B} {12,S}
7  C u0 p0 c0 {5,B} {9,B} {13,S}
8  C u0 p0 c0 {6,B} {9,B} {14,S}
9  C u0 p0 c0 {7,B} {8,B} {15,S}
10 C u0 p0 c0 {1,S} {3,D} {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(29306)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,S} {14,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,D} {13,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(38213)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
3 C u0 p0 c0 {2,S} {5,D} {8,S}
4 C u0 p0 c0 {1,D} {2,S} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {9,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(37708)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {3,D} {7,S}
7  C u1 p0 c0 {5,D} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(227)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {26,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {10,D} {29,S}
10 C u0 p0 c0 {6,S} {9,D} {30,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(228)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
2  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {9,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {26,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u0 p0 c0 {9,D} {29,S} {30,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
        """),
)



species(
    label='C7H13(230)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u1 p0 c0 {5,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(231)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(232)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {22,S} {23,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {6,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(39999)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {7,D}
3  C u1 p0 c0 {4,S} {7,S} {11,S}
4  C u0 p0 c0 {3,S} {8,D} {10,S}
5  C u0 p0 c0 {6,S} {9,D} {12,S}
6  C u1 p0 c0 {1,S} {5,S} {13,S}
7  C u0 p0 c0 {1,S} {2,D} {3,S}
8  C u0 p0 c0 {4,D} {14,S} {15,S}
9  C u0 p0 c0 {5,D} {16,S} {17,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(12629)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {3,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(28649)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {5,D}
3 O u1 p2 c0 {1,S}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u0 p0 c0 {4,D} {8,S} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(39153)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {1,S}
4  C u1 p0 c0 {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,D} {8,S}
6  C u0 p0 c0 {1,S} {2,D} {4,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(233)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
4  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {5,S} {18,S} {19,S}
7  C u0 p0 c0 {1,S} {10,S} {24,S} {25,S}
8  C u0 p0 c0 {9,S} {10,S} {22,S} {23,S}
9  C u0 p0 c0 {8,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {7,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(234)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {7,S} {22,S} {23,S}
4  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
8  C u0 p0 c0 {4,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {10,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {8,S} {9,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(235)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {7,S} {24,S} {25,S}
5  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
6  C u0 p0 c0 {5,S} {7,S} {20,S} {21,S}
7  C u0 p0 c0 {4,S} {6,S} {22,S} {23,S}
8  C u0 p0 c0 {2,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {8,S} {10,S} {26,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(236)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {17,S} {18,S}
2  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {2,S} {9,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {5,S} {27,S} {28,S} {29,S}
10 C u1 p0 c0 {6,S} {7,S} {8,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(237)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {5,S} {22,S} {23,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(41373)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {5,D} {7,S} {11,S}
7  C u0 p0 c0 {2,S} {3,D} {6,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(41371)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u0 p0 c0 {2,S} {3,D} {4,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(238)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {5,S} {20,S} {21,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {10,S} {22,S} {23,S}
8  C u0 p0 c0 {5,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {6,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {7,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(21432)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,D} {9,S}
6  C u0 p0 c0 {2,S} {3,D} {4,S}
7  C u0 p0 c0 {5,D} {10,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(2982)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u0 p0 c0 {1,S} {8,D} {17,S}
8  C u0 p0 c0 {6,D} {7,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(41366)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u0 p2 c0 {2,S} {4,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5 C u1 p0 c0 {4,S} {6,S} {9,S}
6 C u0 p0 c0 {2,S} {3,D} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(41449)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {5,D}
4  C u0 p0 c0 {5,S} {6,D} {7,S}
5  C u0 p0 c0 {1,S} {3,D} {4,S}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(43582)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,D} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {4,D} {8,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(239)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {7,S} {18,S} {19,S}
4  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
5  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {10,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {6,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {7,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C9H16(240)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
8  C u0 p0 c0 {7,S} {9,D} {23,S}
9  C u0 p0 c0 {8,D} {24,S} {25,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C8H14(241)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,D} {20,S}
8  C u0 p0 c0 {7,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(242)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {9,S} {21,S} {22,S}
6  C u0 p0 c0 {1,S} {10,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {9,D} {26,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 C u1 p0 c0 {6,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C7H12(243)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {18,S} {19,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(244)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {5,S} {10,D}
9  C u1 p0 c0 {6,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C6H10(245)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,D} {15,S}
6  C u0 p0 c0 {3,S} {5,D} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(246)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {9,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {10,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {9,D} {26,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 C u1 p0 c0 {6,S} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(3101)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {10,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
9  C u0 p0 c0 {4,S} {10,D} {19,S}
10 C u0 p0 c0 {2,S} {9,D} {20,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(247)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {1,S} {10,D} {25,S}
9  C u1 p0 c0 {6,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(44328)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {6,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {3,D} {4,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(44330)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u0 p0 c0 {3,D} {4,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(248)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u1 p0 c0 {4,S} {5,S} {26,S}
9  C u0 p0 c0 {6,S} {10,D} {27,S}
10 C u0 p0 c0 {9,D} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(249)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
6  C u0 p0 c0 {1,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {5,S} {10,D} {25,S}
9  C u1 p0 c0 {6,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(250)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {10,D} {25,S}
9  C u1 p0 c0 {1,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(251)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
8  C u1 p0 c0 {5,S} {10,S} {27,S}
9  C u0 p0 c0 {4,S} {10,D} {28,S}
10 C u0 p0 c0 {8,S} {9,D} {29,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C8H15(252)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {2,S} {6,D} {21,S}
8  C u1 p0 c0 {4,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(253)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {9,D} {10,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 C u1 p0 c0 {8,S} {28,S} {29,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C8H15(254)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(255)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {3,S} {26,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {10,D} {29,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C8H15(256)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {21,S}
7  C u0 p0 c0 {4,S} {6,D} {20,S}
8  C u1 p0 c0 {4,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(257)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
8  C u1 p0 c0 {4,S} {5,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {27,S}
10 C u0 p0 c0 {9,D} {28,S} {29,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C8H15(258)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {1,S} {8,D} {19,S}
7  C u1 p0 c0 {1,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(259)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {8,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {6,S} {9,D} {26,S}
9  C u0 p0 c0 {8,D} {10,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C7H12(260)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u0 p0 c0 {3,S} {7,D} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  C u0 p0 c0 {5,D} {18,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(261)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {1,S} {9,D} {26,S}
9  C u0 p0 c0 {8,D} {10,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C8H15(262)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {3,S} {4,S} {20,S}
7  C u0 p0 c0 {4,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(263)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {1,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {5,S} {9,D} {26,S}
9  C u0 p0 c0 {8,D} {10,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
        """),
)


species(
    label='C6H10(264)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  C u0 p0 c0 {4,D} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C8H14(266)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {6,D} {18,S}
6  C u0 p0 c0 {5,D} {7,S} {20,S}
7  C u0 p0 c0 {6,S} {8,D} {19,S}
8  C u0 p0 c0 {7,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H14(267)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {18,S}
7  C u0 p0 c0 {5,D} {21,S} {22,S}
8  C u0 p0 c0 {6,D} {19,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1086)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u1 p0 c0 {1,S} {6,S} {13,S}
4  C u0 p0 c0 {2,S} {6,D} {15,S}
5  C u0 p0 c0 {1,S} {7,D} {14,S}
6  C u0 p0 c0 {3,S} {4,D} {16,S}
7  C u0 p0 c0 {5,D} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H12(268)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {15,S}
6  C u0 p0 c0 {4,D} {18,S} {19,S}
7  C u0 p0 c0 {5,D} {16,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(1085)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
3  C u1 p0 c0 {1,S} {2,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {6,S} {16,S}
6  C u0 p0 c0 {5,S} {7,D} {15,S}
7  C u0 p0 c0 {6,D} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(44543)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {7,S} {16,S}
7  C u1 p0 c0 {6,S} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1117)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {4,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {17,S}
9  C u0 p0 c0 {8,D} {19,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(44547)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {5,S} {14,S}
5  C u0 p0 c0 {4,S} {7,D} {13,S}
6  C u1 p0 c0 {1,S} {15,S} {16,S}
7  C u0 p0 c0 {5,D} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(44550)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {5,S} {6,D} {17,S}
8  C u0 p0 c0 {4,S} {9,D} {16,S}
9  C u0 p0 c0 {8,D} {19,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(1115)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {16,S}
7  C u0 p0 c0 {6,D} {8,S} {18,S}
8  C u0 p0 c0 {7,S} {9,D} {17,S}
9  C u0 p0 c0 {8,D} {19,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(44690)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
7  C u1 p0 c0 {5,S} {8,S} {17,S}
8  C u0 p0 c0 {7,S} {9,D} {18,S}
9  C u0 p0 c0 {8,D} {19,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(44884)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {5,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {16,S}
9  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
10 C u0 p0 c0 {8,S} {11,D} {20,S}
11 C u0 p0 c0 {10,D} {21,S} {22,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(44885)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {10,S} {15,S}
8  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
9  C u0 p0 c0 {3,S} {11,S} {19,S} {20,S}
10 C u0 p0 c0 {7,S} {11,D} {21,S}
11 C u0 p0 c0 {9,S} {10,D} {22,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(44670)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
6  C u0 p0 c0 {8,S} {14,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {18,S}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u1 p0 c0 {5,S} {19,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(44698)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {4,S} {9,D} {18,S}
8  C u1 p0 c0 {6,S} {9,S} {19,S}
9  C u0 p0 c0 {7,D} {8,S} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H18(269)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {7,S} {23,S} {24,S}
5  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {7,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {6,S} {21,S} {22,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {8,S} {25,S} {26,S} {27,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H18(270)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {2,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {9,D} {26,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H18(271)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {2,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {8,D} {26,S} {27,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C6H11(272)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H12(273)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C9H17(274)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {8,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {1,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(44641)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {5,S}
3  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
4  C u0 p0 c0 {7,S} {9,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
7  C u1 p0 c0 {4,S} {5,S} {18,S}
8  C u0 p0 c0 {3,S} {9,D} {20,S}
9  C u0 p0 c0 {4,S} {8,D} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C9H17(275)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {7,S} {8,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(44655)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {4,D} {7,S} {16,S}
6  C u0 p0 c0 {7,D} {17,S} {18,S}
7  C u1 p0 c0 {5,S} {6,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(44654)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,D} {15,S}
5  C u0 p0 c0 {6,D} {7,S} {16,S}
6  C u0 p0 c0 {5,D} {17,S} {18,S}
7  C u1 p0 c0 {4,D} {5,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C9H17(276)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {3,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C9H17(277)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {1,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(278)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {6,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(279)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {7,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {6,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(280)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {9,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {8,D} {23,S}
8  C u0 p0 c0 {5,S} {7,D} {24,S}
9  C u1 p0 c0 {4,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(281)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {8,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {8,D} {23,S}
8  C u0 p0 c0 {4,S} {7,D} {24,S}
9  C u1 p0 c0 {5,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(282)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {1,S} {9,D} {22,S}
8  C u1 p0 c0 {5,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(283)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {7,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {5,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(284)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {7,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u1 p0 c0 {3,S} {4,S} {23,S}
8  C u0 p0 c0 {5,S} {9,D} {24,S}
9  C u0 p0 c0 {8,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(285)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {4,S} {9,D} {22,S}
8  C u1 p0 c0 {5,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(286)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {9,D} {22,S}
8  C u1 p0 c0 {1,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H13(287)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {18,S}
6  C u0 p0 c0 {3,S} {5,D} {17,S}
7  C u1 p0 c0 {3,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C9H17(288)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {3,S} {9,S} {24,S}
8  C u0 p0 c0 {4,S} {9,D} {25,S}
9  C u0 p0 c0 {7,S} {8,D} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(289)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
7  C u1 p0 c0 {3,S} {9,S} {24,S}
8  C u0 p0 c0 {4,S} {9,D} {25,S}
9  C u0 p0 c0 {7,S} {8,D} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H13(290)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {7,D} {16,S}
6  C u1 p0 c0 {1,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(1082)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,D} {13,S}
4  C u0 p0 c0 {1,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {7,S} {14,S}
6  C u0 p0 c0 {4,D} {17,S} {18,S}
7  C u1 p0 c0 {5,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C9H17(291)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
7  C u1 p0 c0 {3,S} {4,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {24,S}
9  C u0 p0 c0 {8,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C7H13(292)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C9H17(293)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {3,S} {8,D} {9,S}
8  C u0 p0 c0 {4,S} {7,D} {24,S}
9  C u1 p0 c0 {7,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(294)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {8,D} {23,S}
8  C u0 p0 c0 {7,D} {9,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(36444)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {8,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {3,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {17,S}
10 C u0 p0 c0 {4,D} {8,S} {18,S}
11 C u0 p0 c0 {9,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='C7H13(295)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {3,S} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H9(958)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,D} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {10,S}
3  C u0 p0 c0 {2,D} {5,S} {11,S}
4  C u0 p0 c0 {1,D} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u1 p0 c0 {4,S} {13,S} {14,S}
7  C u0 p0 c0 {5,D} {15,S} {16,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C9H17(296)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {1,S} {8,D} {23,S}
8  C u0 p0 c0 {7,D} {9,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C9H17(297)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {8,D} {23,S}
8  C u0 p0 c0 {7,D} {9,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C6H10(298)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u0 p0 c0 {6,D} {15,S} {16,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C8H16(299)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {20,S} {21,S}
4  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {4,S} {22,S} {23,S} {24,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H16(300)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {8,S} {15,S} {16,S}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {8,D} {23,S}
8  C u0 p0 c0 {3,S} {7,D} {24,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H16(301)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {23,S} {24,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C5H9(302)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C5H10(303)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C8H15(304)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {1,S} {7,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(305)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(45160)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,S} {8,D} {11,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,D} {14,S}
8  C u0 p0 c0 {5,D} {15,S} {16,S}
9  C u0 p0 c0 {7,D} {17,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(45036)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {9,D}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {11,D} {14,S}
8  C u1 p0 c0 {6,S} {16,S} {17,S}
9  C u0 p0 c0 {3,D} {5,S} {15,S}
10 C u0 p0 c0 {4,D} {6,S} {18,S}
11 C u0 p0 c0 {7,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(45159)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {6,S} {13,S}
6  C u0 p0 c0 {5,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {8,S} {15,S}
8  C u0 p0 c0 {7,S} {9,D} {16,S}
9  C u0 p0 c0 {8,D} {17,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(45037)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {9,S}
3  O u0 p2 c0 {10,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {9,S} {11,D} {15,S}
8  C u0 p0 c0 {4,D} {5,S} {18,S}
9  C u1 p0 c0 {2,S} {7,S} {17,S}
10 C u0 p0 c0 {3,D} {6,S} {16,S}
11 C u0 p0 c0 {7,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='C8H15(306)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {1,S} {6,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(45280)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {2,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {15,S}
8  C u0 p0 c0 {7,S} {11,D} {16,S}
9  C u1 p0 c0 {1,S} {6,S} {17,S}
10 C u0 p0 c0 {4,D} {5,S} {18,S}
11 C u0 p0 c0 {8,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(45262)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {6,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
5  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {3,D} {4,S} {12,S}
7  C u0 p0 c0 {2,D} {5,S} {11,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(45279)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {7,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {3,S} {8,S} {15,S}
8  C u0 p0 c0 {7,S} {11,D} {16,S}
9  C u1 p0 c0 {1,S} {5,S} {17,S}
10 C u0 p0 c0 {4,D} {6,S} {18,S}
11 C u0 p0 c0 {8,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='C8H15(307)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {6,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(308)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {6,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(309)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {5,S} {6,D} {21,S}
8  C u1 p0 c0 {4,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(45205)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {8,S} {14,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {3,S} {9,D} {13,S}
8  C u0 p0 c0 {5,S} {6,D} {16,S}
9  C u0 p0 c0 {7,D} {17,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(45554)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {5,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {9,D} {18,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 C u0 p0 c0 {6,S} {11,D} {16,S}
11 C u0 p0 c0 {10,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(45553)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 C u0 p0 c0 {5,S} {11,D} {16,S}
11 C u0 p0 c0 {10,D} {19,S} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(45501)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {6,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {8,D} {16,S}
8  C u0 p0 c0 {5,S} {7,D} {15,S}
9  C u1 p0 c0 {4,S} {17,S} {18,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(45480)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {8,S}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u0 p0 c0 {4,D} {7,S} {14,S}
7  C u0 p0 c0 {5,D} {6,S} {15,S}
8  C u0 p0 c0 {4,S} {9,D} {13,S}
9  C u0 p0 c0 {8,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(45483)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,S} {9,D} {11,S}
6  C u0 p0 c0 {4,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {8,D} {14,S}
8  C u0 p0 c0 {2,S} {7,D} {15,S}
9  C u0 p0 c0 {5,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(45260)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {4,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u1 p0 c0 {4,S} {8,S} {15,S}
6  C u0 p0 c0 {3,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {9,S} {16,S}
8  C u0 p0 c0 {5,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(310)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(46144)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  C u1 p0 c0 {3,S} {5,S} {12,S}
5  C u0 p0 c0 {4,S} {7,D} {15,S}
6  C u0 p0 c0 {2,D} {7,S} {8,S}
7  C u0 p0 c0 {5,D} {6,S} {14,S}
8  C u0 p0 c0 {6,S} {9,D} {13,S}
9  C u0 p0 c0 {8,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(5707)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {3,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u1 p0 c0 {3,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {8,D} {11,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {9,S} {14,S}
8  C u0 p0 c0 {5,D} {15,S} {16,S}
9  C u0 p0 c0 {2,D} {7,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(23031)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 C u1 p0 c0 {1,S} {3,S} {8,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C8H15(311)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {7,D} {20,S}
7  C u0 p0 c0 {4,S} {6,D} {21,S}
8  C u1 p0 c0 {3,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(46722)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  O u0 p2 c0 {8,D}
3  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {6,S} {14,S}
6  C u0 p0 c0 {1,D} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {9,D} {13,S}
8  C u0 p0 c0 {2,D} {3,S} {15,S}
9  C u0 p0 c0 {7,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(46727)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {2,D} {6,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {14,S}
9  C u0 p0 c0 {8,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(46152)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u0 p0 c0 {2,D} {3,S} {8,S}
7  C u0 p0 c0 {4,S} {5,D} {13,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {16,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(46721)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {8,S} {17,S}
2  O u0 p2 c0 {4,D}
3  C u0 p0 c0 {5,D} {6,S} {13,S}
4  C u0 p0 c0 {2,D} {5,S} {7,S}
5  C u0 p0 c0 {3,D} {4,S} {12,S}
6  C u0 p0 c0 {3,S} {8,D} {10,S}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {6,D} {14,S}
9  C u0 p0 c0 {7,D} {15,S} {16,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C8H15(312)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {1,S} {8,D} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(313)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {2,S} {4,S} {20,S}
7  C u0 p0 c0 {3,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(314)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {8,D} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(315)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,D} {19,S}
7  C u1 p0 c0 {1,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(956)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u1 p0 c0 {5,S} {9,S} {18,S}
8  C u0 p0 c0 {6,S} {9,D} {17,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(48836)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {7,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {3,S} {8,S} {14,S}
6  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {17,S}
8  C u0 p0 c0 {5,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(48886)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {6,S}
2  O u1 p2 c0 {7,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {16,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {17,S}
8  C u0 p0 c0 {6,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C8H15(316)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {7,S} {18,S} {19,S} {20,S}
6  C u1 p0 c0 {3,S} {8,S} {21,S}
7  C u0 p0 c0 {5,S} {8,D} {22,S}
8  C u0 p0 c0 {6,S} {7,D} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(48904)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {9,S} {17,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {9,S} {18,S}
9  C u1 p0 c0 {4,S} {8,S} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(48913)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {19,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {15,S}
7  C u1 p0 c0 {5,S} {6,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {17,S}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(48907)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {15,S}
7  C u0 p0 c0 {6,D} {9,S} {16,S}
8  C u1 p0 c0 {5,S} {17,S} {18,S}
9  C u0 p0 c0 {2,D} {7,S} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C6H11(317)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  C u1 p0 c0 {2,S} {16,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(318)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(319)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {2,S} {4,D} {14,S}
6  C u1 p0 c0 {2,S} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C8H15(320)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u1 p0 c0 {2,S} {8,S} {21,S}
7  C u0 p0 c0 {3,S} {8,D} {22,S}
8  C u0 p0 c0 {6,S} {7,D} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(322)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {7,D} {20,S}
7  C u0 p0 c0 {6,D} {8,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1084)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {2,S} {3,D} {14,S}
5  C u0 p0 c0 {1,S} {7,D} {13,S}
6  C u1 p0 c0 {2,S} {15,S} {16,S}
7  C u0 p0 c0 {5,D} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C8H15(323)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {1,S} {7,D} {20,S}
7  C u0 p0 c0 {6,D} {8,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C8H15(324)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {6,D} {8,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2983)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u1 p0 c0 {1,D} {7,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(3571)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,D}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {1,D} {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {9,S} {15,S}
8  C u1 p0 c0 {4,S} {16,S} {17,S}
9  C u0 p0 c0 {2,D} {7,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(48983)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {9,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {14,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {17,S}
8  C u0 p0 c0 {7,D} {9,S} {18,S}
9  C u0 p0 c0 {2,D} {8,S} {19,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(48984)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u1 p0 c0 {6,S} {17,S} {18,S}
9  C u0 p0 c0 {2,S} {7,D} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='C8H15(325)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {3,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {6,D} {21,S}
8  C u1 p0 c0 {6,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C6H10(326)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {15,S} {16,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(5118)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u1 p2 c0 {1,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
9  C u0 p0 c0 {4,S} {10,D} {19,S}
10 C u0 p0 c0 {3,D} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(15233)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {8,D} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u0 p0 c0 {2,S} {6,D} {9,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(49070)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {10,S}
2  O u0 p2 c0 {10,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {8,D} {10,S} {19,S}
10 C u0 p0 c0 {1,S} {2,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(49013)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {4,S} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {16,S}
7  C u0 p0 c0 {6,D} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(49198)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {8,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {19,S}
9  C u1 p0 c0 {6,S} {7,S} {8,S}
10 C u0 p0 c0 {2,S} {3,D} {8,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(49279)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {1,S} {8,D} {10,S}
10 C u0 p0 c0 {2,S} {3,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(49278)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {6,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
7  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {18,S}
10 C u0 p0 c0 {2,S} {3,D} {6,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(49284)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {9,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {1,S} {8,D} {19,S}
10 C u1 p0 c0 {2,S} {3,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(49295)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {7,S}
2  O u0 p2 c0 {1,S} {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {4,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
9  C u0 p0 c0 {6,S} {8,S} {19,S} {20,S}
10 C u0 p0 c0 {2,S} {3,D} {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(5797)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {1,S} {6,D} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(49552)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {9,D}
2  O u1 p2 c0 {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u1 p0 c0 {6,S} {7,S} {9,S}
9  C u0 p0 c0 {1,D} {8,S} {10,S}
10 C u0 p0 c0 {2,S} {3,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(49184)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='CO2(12610)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
        """),
)


species(
    label='C8H15(327)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u1 p0 c0 {2,S} {3,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
        """),
)



species(
    label='C7H15(329)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {3,S} {6,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H14(330)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {19,S} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(333)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {7,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(49575)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {7,S}
2  O u1 p2 c0 {10,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {17,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 C u0 p0 c0 {2,S} {3,D} {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(49950)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {10,S} {18,S}
2  O u0 p2 c0 {9,D}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {8,D} {9,S}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u0 p0 c0 {2,D} {7,S} {10,S}
10 C u0 p0 c0 {1,S} {3,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(49942)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {1,D} {6,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(49557)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {9,D}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
9  C u0 p0 c0 {2,D} {4,S} {10,S}
10 C u0 p0 c0 {1,S} {3,D} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(49945)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  C u1 p0 c0 {4,S} {5,S} {7,S}
7  C u1 p0 c0 {1,D} {6,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C7H13(334)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {1,S} {5,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(335)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {4,S} {5,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(32085)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {7,D}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {7,S}
6  C u0 p0 c0 {4,S} {5,D} {14,S}
7  C u0 p0 c0 {1,D} {5,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(50151)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
5  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
9  C u0 p0 c0 {2,D} {4,S} {18,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(49944)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {7,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u1 p0 c0 {4,S} {5,S} {14,S}
7  C u0 p0 c0 {1,S} {5,D} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(336)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {4,S} {5,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(337)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(338)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {17,S}
6  C u0 p0 c0 {4,S} {5,D} {18,S}
7  C u1 p0 c0 {3,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(339)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {7,D} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(341)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {7,D} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(342)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,D} {16,S}
6  C u1 p0 c0 {1,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(343)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {18,S}
7  C u1 p0 c0 {6,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(344)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {18,S}
7  C u1 p0 c0 {6,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(345)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {6,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {7,S} {18,S}
6  C u0 p0 c0 {4,S} {7,D} {19,S}
7  C u0 p0 c0 {5,S} {6,D} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H9(346)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C7H13(347)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {18,S}
7  C u1 p0 c0 {5,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H8(348)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u0 p0 c0 {2,D} {12,S} {13,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C7H13(349)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H13(350)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {18,S}
7  C u1 p0 c0 {5,S} {19,S} {20,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C5H8(351)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {5,D} {12,S} {13,S}
5  C u0 p0 c0 {3,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C6H12(352)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {5,S} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(353)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {4,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u1 p0 c0 {3,S} {7,S} {18,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(354)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(355)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {10,S} {19,S} {20,S}
9  C u0 p0 c0 {7,S} {10,S} {15,S} {16,S}
10 C u0 p0 c0 {8,S} {9,S} {17,S} {18,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(356)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {6,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {9,S} {19,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {8,S} {17,S} {18,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(357)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
6  C u1 p0 c0 {1,S} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(358)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  C u1 p0 c0 {3,S} {4,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H11(359)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,D} {13,S}
5  C u1 p0 c0 {1,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
        """),
)


species(
    label='C6H9(360)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,D} {6,S} {11,S}
5  C u0 p0 c0 {3,D} {14,S} {15,S}
6  C u1 p0 c0 {4,S} {12,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C6H8(361)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,D} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {6,S} {13,S}
6  C u0 p0 c0 {3,D} {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
        """),
)


# Reaction systems
simpleReactor(
    temperature=[(600,'K'),(1500,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=6,
    initialMoleFractions={
        #"C7H10": 1,
        "O2(2)": 0.00622, # phi=1 means 9.5 O2 per C7H10
        "N2":   0.36, # 8.1 times as much N2 as O2
        "CHPD(1)":0.00989, # Cycloheptadiene C7H10    
        "C7H16(1)": 0.00463,
"H(3)": 0,
"O(4)": 0,
"OH(5)": 0,
"H2(6)": 0,
"H2O(7)": 0,
"Ar(8)": 0,
"He(9)": 0,
"HO2(10)": 0,
"H2O2(11)": 0,
"CO(12)": 0,
"CO2(13)": 0,
"HCO(14)": 0,
"C(T)(15)": 0,
"CH(16)": 0,
"CH2(T)(17)": 0,
"CH3(18)": 0,
"CH2O(19)": 0,
"HCCO(20)": 0,
"C2H(21)": 0,
"C2H2(22)": 0,
"H2CC(23)": 0,
"CH2(S)(24)": 0,
"CH3OH(25)": 0,
"CH3O(26)": 0,
"CH2CO(27)": 0,
"C2H3(28)": 0,
"C2H4(29)": 0,
"CH4(30)": 0,
"C2H6(31)": 0,
"C2H5(32)": 0,
"CH2OH(33)": 0,
"CH3CO(34)": 0,
"CH2CHO(35)": 0,
"CH3CHO(36)": 0,
"Ne": 0,
"C4H9(69)": 0,
"C5H11(71)": 0,
"C3H7(70)": 0,
"C7H15(72)": 0,
"C7H15(73)": 0,
"C7H15(75)": 0,
"S(150)": 0,
"S(152)": 0,
"S(151)": 0,
"C2H5O2(54)": 0,
"S(107)": 0,
"S(106)": 0,
"S(108)": 0,
"C6H12(114)": 0,
"C3H6(59)": 0,
"C4H8(81)": 0,
"S(214)": 0,
"S(666)": 0,
"S(459)": 0,
"S(472)": 0,
"S(1097)": 0,
"CH3OO(41)": 0,
"C5H10(88)": 0,
"HCCOH(37)": 0,
"HOC2H2(38)": 0,
"HOCO(39)": 0,
"CH3OOH(40)": 0,
"C2H5O(42)": 0,
"cC2H4O(43)": 0,
"C7H15(76)": 0,
"C2H3O2(44)": 0,
"S(422)": 0,
"C2H4(78)": 0,
"C2(45)": 0,
"CHO3(486)": 0,
"S(412)": 0,
"CHO4(489)": 0,
"C4H7(639)": 0,
"S(1246)": 0,
"C2H6O(47)": 0,
"C2H4O(111)": 0,
"S(1398)": 0,
"S(1399)": 0,
"S(1460)": 0,
"C2H5O(48)": 0,
"S(179)": 0,
"S(699)": 0,
"S(201)": 0,
"C2H5O(49)": 0,
"C4H6(1293)": 0,
"S(1896)": 0,
"S(2085)": 0,
"C2H4O(50)": 0,
"C4H6(1947)": 0,
"S(2088)": 0,
"cC2H3O(51)": 0,
"OCHCHO(52)": 0,
"C2H6O2(53)": 0,
"S(1285)": 0,
"S(519)": 0,
"S(2156)": 0,
"S(2836)": 0,
"S(452)": 0,
"S(3108)": 0,
"C2H5O2(55)": 0,
"C2H4O2(56)": 0,
"C4H6(1295)": 0,
"C3H6O(57)": 0,
"C3H5(61)": 0,
"C2H5CO(58)": 0,
"CHO3(477)": 0,
"C3H4(60)": 0,
"C3H5(62)": 0,
"C3H5(63)": 0,
"C3H4O(64)": 0,
"S(1454)": 0,
"H3CCCH(65)": 0,
"C3H4O(66)": 0,
"S(3313)": 0,
"H2CCCH(67)": 0,
"S(3121)": 0,
"S(3901)": 0,
"S(3912)": 0,
"S(3296)": 0,
"S(4476)": 0,
"C6H13(74)": 0,
"S(1952)": 0,
"S(4094)": 0,
"S(3316)": 0,
"CCCCCC(77)": 0,
"C2H4(79)": 0,
"S(4480)": 0,
"S(4646)": 0,
"C3H6(82)": 0,
"C4H8(83)": 0,
"C4H8(84)": 0,
"S(4662)": 0,
"C6H13(124)": 0,
"S(196)": 0,
"C4H8(85)": 0,
"C4H9(86)": 0,
"S(544)": 0,
"S(3136)": 0,
"C5H10(89)": 0,
"C5H10(90)": 0,
"S(3528)": 0,
"C5H10(91)": 0,
"S(174)": 0,
"S(4636)": 0,
"C5H10(92)": 0,
"C5H11(93)": 0,
"C5H11(94)": 0,
"C3H7(96)": 0,
"S(3085)": 0,
"C3H6(97)": 0,
"S(723)": 0,
"S(5151)": 0,
"S(102)": 0,
"C4H10(103)": 0,
"CCCCC(104)": 0,
"S(2947)": 0,
"C3H8(105)": 0,
"S(5291)": 0,
"S(5293)": 0,
"S(5328)": 0,
"S(5314)": 0,
"S(2245)": 0,
"CH2O(110)": 0,
"S(1162)": 0,
"S(1108)": 0,
"S(3071)": 0,
"S(1934)": 0,
"S(4489)": 0,
"S(4619)": 0,
"S(3254)": 0,
"C7H14(112)": 0,
"S(218)": 0,
"S(1525)": 0,
"S(1313)": 0,
"S(6658)": 0,
"C7H14(113)": 0,
"S(1316)": 0,
"S(204)": 0,
"S(5632)": 0,
"S(6402)": 0,
"S(6835)": 0,
"S(6834)": 0,
"S(7425)": 0,
"C7H14(115)": 0,
"S(5641)": 0,
"S(5640)": 0,
"C6H12(116)": 0,
"S(7780)": 0,
"S(7902)": 0,
"O2(1283)": 0,
"S(7823)": 0,
"S(8143)": 0,
"C7H14(117)": 0,
"S(5303)": 0,
"S(4606)": 0,
"C7H14(118)": 0,
"S(2492)": 0,
"C6H12(119)": 0,
"C7H14(120)": 0,
"C7H14(121)": 0,
"C7H14(122)": 0,
"C6H13(123)": 0,
"S(8464)": 0,
"S(9208)": 0,
"C7H15(125)": 0,
"S(173)": 0,
"S(7008)": 0,
"S(217)": 0,
"S(182)": 0,
"C7H15(126)": 0,
"S(8466)": 0,
"C6H12(128)": 0,
"C7H14(129)": 0,
"C7H14(130)": 0,
"C7H14(132)": 0,
"C6H12(133)": 0,
"C7H14(134)": 0,
"S(2176)": 0,
"S(10607)": 0,
"S(9707)": 0,
"S(10779)": 0,
"S(4093)": 0,
"C7H14(135)": 0,
"C7H14(136)": 0,
"S(10905)": 0,
"S(8813)": 0,
"S(139)": 0,
"S(140)": 0,
"S(10933)": 0,
"S(141)": 0,
"S(142)": 0,
"S(143)": 0,
"S(10782)": 0,
"S(10958)": 0,
"S(10965)": 0,
"CO4(145)": 0,
"CH2O3(147)": 0,
"S(2502)": 0,
"S(11215)": 0,
"S(11214)": 0,
"S(3105)": 0,
"S(11112)": 0,
"C4H8(149)": 0,
"C2H2O(153)": 0,
"C2H2O(154)": 0,
"S(11094)": 0,
"S(3165)": 0,
"C2H2(155)": 0,
"S(11713)": 0,
"S(3075)": 0,
"S(1938)": 0,
"S(725)": 0,
"S(10423)": 0,
"C2H2O(158)": 0,
"C2H3O(159)": 0,
"S(160)": 0,
"S(161)": 0,
"S(162)": 0,
"S(163)": 0,
"S(164)": 0,
"S(165)": 0,
"S(5223)": 0,
"S(10299)": 0,
"S(166)": 0,
"S(167)": 0,
"S(168)": 0,
"S(169)": 0,
"S(170)": 0,
"S(171)": 0,
"S(172)": 0,
"S(11870)": 0,
"S(175)": 0,
"S(5149)": 0,
"S(176)": 0,
"S(1407)": 0,
"S(177)": 0,
"S(178)": 0,
"S(180)": 0,
"S(5273)": 0,
"S(12790)": 0,
"S(12838)": 0,
"S(13006)": 0,
"S(12842)": 0,
"S(5636)": 0,
"S(5624)": 0,
"S(181)": 0,
"S(1409)": 0,
"S(183)": 0,
"S(184)": 0,
"S(185)": 0,
"S(186)": 0,
"S(187)": 0,
"S(188)": 0,
"S(13005)": 0,
"S(3074)": 0,
"S(189)": 0,
"S(190)": 0,
"S(191)": 0,
"S(13690)": 0,
"S(192)": 0,
"S(193)": 0,
"S(194)": 0,
"S(4711)": 0,
"S(14115)": 0,
"S(195)": 0,
"S(197)": 0,
"S(198)": 0,
"S(199)": 0,
"S(200)": 0,
"S(202)": 0,
"S(701)": 0,
"S(203)": 0,
"S(205)": 0,
"S(14326)": 0,
"S(206)": 0,
"S(207)": 0,
"S(208)": 0,
"S(209)": 0,
"S(210)": 0,
"S(211)": 0,
"S(212)": 0,
"S(213)": 0,
"S(14146)": 0,
"S(215)": 0,
"S(14155)": 0,
"S(14411)": 0,
"S(216)": 0,
"S(219)": 0,
"CH2O2(220)": 0,
"CHO2(221)": 0,
"CH3O3(223)": 0,
"CH3O3(224)": 0,
"S(2711)": 0,
"S(3082)": 0,
"CH2O3(228)": 0,
"C3H6O(229)": 0,
"S(230)": 0,
"S(5226)": 0,
"S(14750)": 0,
"S(231)": 0,
"S(702)": 0,
"S(232)": 0,
"S(5182)": 0,
"S(13729)": 0,
"S(14107)": 0,
"S(233)": 0,
"S(13826)": 0,
"S(609)": 0,
"S(4634)": 0,
"C3H2O(503)": 0,
"S(16528)": 0,
"S(234)": 0,
"S(235)": 0,
"S(236)": 0,
"S(13726)": 0,
"S(237)": 0,
"S(726)": 0,
"S(17462)": 0,
"S(5335)": 0,
"S(238)": 0,
"S(12794)": 0,
"S(2724)": 0,
"S(239)": 0,
"S(240)": 0,
"S(13166)": 0,
"S(13880)": 0,
"S(696)": 0,
"S(241)": 0,
"S(242)": 0,
"S(243)": 0,
"S(244)": 0,
"S(17880)": 0,
"S(18819)": 0,
"S(245)": 0,
"S(246)": 0,
"S(247)": 0,
"S(248)": 0,
"S(249)": 0,
"S(250)": 0,
"S(251)": 0,
"S(252)": 0,
"S(253)": 0,
"S(254)": 0,
"S(255)": 0,
"S(256)": 0,
"S(257)": 0,
"S(258)": 0,
"S(259)": 0,
"S(260)": 0,
"S(261)": 0,
"S(262)": 0,
"S(18882)": 0,
"S(263)": 0,
"S(264)": 0,
"S(265)": 0,
"S(266)": 0,
"S(7422)": 0,
"S(7807)": 0,
"S(18862)": 0,
"S(267)": 0,
"S(268)": 0,
"S(269)": 0,
"S(270)": 0,
"S(271)": 0,
"S(272)": 0,
"S(8159)": 0,
"S(273)": 0,
"S(274)": 0,
"S(275)": 0,
"S(276)": 0,
"S(277)": 0,
"C3H5O(560)": 0,
"S(278)": 0,
"S(2720)": 0,
"S(8460)": 0,
"S(738)": 0,
"S(279)": 0,
"S(280)": 0,
"S(281)": 0,
"S(282)": 0,
"S(283)": 0,
"S(284)": 0,
"S(285)": 0,
"S(286)": 0,
"S(4629)": 0,
"S(287)": 0,
"S(14151)": 0,
"S(21667)": 0,
"S(6972)": 0,
"S(288)": 0,
"S(289)": 0,
"S(290)": 0,
"S(291)": 0,
"S(292)": 0,
"S(293)": 0,
"S(294)": 0,
"S(295)": 0,
"S(18850)": 0,
"S(15977)": 0,
"S(9205)": 0,
"S(296)": 0,
"S(297)": 0,
"S(298)": 0,
"S(299)": 0,
"S(300)": 0,
"S(301)": 0,
"S(21970)": 0,
"S(302)": 0,
"S(303)": 0,
"S(304)": 0,
"S(305)": 0,
"S(306)": 0,
"S(307)": 0,
"S(21969)": 0,
"S(308)": 0,
"S(15968)": 0,
"S(20975)": 0,
"S(309)": 0,
"S(310)": 0,
"S(311)": 0,
"S(312)": 0,
"S(313)": 0,
"S(314)": 0,
"S(12285)": 0,
"S(315)": 0,
"S(3473)": 0,
"S(316)": 0,
"S(317)": 0,
"S(318)": 0,
"S(319)": 0,
"S(320)": 0,
"S(321)": 0,
"S(322)": 0,
"S(23622)": 0,
"S(323)": 0,
"S(324)": 0,
"S(325)": 0,
"S(326)": 0,
"S(327)": 0,
"S(328)": 0,
"S(330)": 0,
"S(331)": 0,
"S(332)": 0,
"S(333)": 0,
"S(334)": 0,
"S(335)": 0,
"S(336)": 0,
"S(13887)": 0,
"S(337)": 0,
"S(8779)": 0,
"S(338)": 0,
"S(339)": 0,
"S(340)": 0,
"S(3146)": 0,
"S(341)": 0,
"S(342)": 0,
"S(343)": 0,
"S(344)": 0,
"S(6996)": 0,
"S(575)": 0,
"S(18983)": 0,
"S(4624)": 0,
"S(26109)": 0,
"S(345)": 0,
"S(14777)": 0,
"S(26141)": 0,
"S(22106)": 0,
"S(346)": 0,
"S(4621)": 0,
"S(347)": 0,
"S(348)": 0,
"S(14409)": 0,
"S(4490)": 0,
"S(350)": 0,
"S(26265)": 0,
"S(351)": 0,
"S(352)": 0,
"S(353)": 0,
"S(354)": 0,
"S(355)": 0,
"S(356)": 0,
"S(357)": 0,
"S(3123)": 0,
"S(358)": 0,
"S(359)": 0,
"S(26264)": 0,
"S(27259)": 0,
"S(360)": 0,
"S(361)": 0,
"S(362)": 0,
"S(363)": 0,
"S(13884)": 0,
"S(364)": 0,
"S(365)": 0,
"S(366)": 0,
"S(367)": 0,
"S(368)": 0,
"S(28250)": 0,
"S(370)": 0,
"S(371)": 0,
"S(28225)": 0,
"S(28359)": 0,
"S(372)": 0,
"S(373)": 0,
"S(28453)": 0,
"S(374)": 0,
"S(375)": 0,
"S(376)": 0,
"S(377)": 0,
"S(378)": 0,
"S(379)": 0,
"S(1322)": 0,
"S(11686)": 0,
"S(29601)": 0,
"S(29494)": 0,
"S(380)": 0,
"S(11918)": 0,
"S(30810)": 0,
"S(381)": 0,
"S(382)": 0,
"S(383)": 0,
"S(30684)": 0,
"S(31796)": 0,
"S(32045)": 0,
"S(384)": 0,
"S(385)": 0,
"S(386)": 0,
"S(387)": 0,
"S(31941)": 0,
"C3H6O(388)": 0,
"C3H6O(389)": 0,
"C3H6O(390)": 0,
"S(32051)": 0,
"S(33188)": 0,
"S(6999)": 0,
"S(33340)": 0,
"S(34412)": 0,
"S(33346)": 0,
"S(391)": 0,
"S(33345)": 0,
"S(34446)": 0,
"S(392)": 0,
"S(393)": 0,
"S(394)": 0,
"S(34617)": 0,
"S(34619)": 0,
"S(395)": 0,
"S(34629)": 0,
"S(34636)": 0,
"S(396)": 0,
"S(397)": 0,
"S(36707)": 0,
"S(398)": 0,
"S(36721)": 0,
"S(399)": 0,
"S(400)": 0,
"S(31940)": 0,
"S(401)": 0,
"S(402)": 0,
"S(404)": 0,
"S(406)": 0,
"S(10412)": 0,
"S(408)": 0,
"S(37340)": 0,
"S(409)": 0,
"S(411)": 0,
"C4H8O(413)": 0,
"C4H8O(414)": 0,
"C4H8O(415)": 0,
"C9H19(416)": 0,
"S(33317)": 0,
"C9H19(417)": 0,
"C9H19(418)": 0,
"S(419)": 0,
"S(420)": 0,
"S(37348)": 0,
"S(421)": 0,
"S(33324)": 0,
"S(31942)": 0,
"S(425)": 0,
"S(426)": 0,
"S(37858)": 0,
"S(5436)": 0,
"S(37864)": 0,
"S(427)": 0,
"S(428)": 0,
"S(429)": 0,
"S(38356)": 0,
"S(430)": 0,
"S(431)": 0,
"S(433)": 0,
"S(33325)": 0,
"S(435)": 0,
"S(39096)": 0,
"C4H5(1942)": 0,
"S(39292)": 0,
"S(39306)": 0,
"S(39307)": 0,
"S(437)": 0,
"S(39288)": 0,
"S(39419)": 0,
"S(39155)": 0,
"S(39753)": 0,
"S(42071)": 0,
"S(39914)": 0,
"S(42076)": 0,
"S(39669)": 0,
"S(39686)": 0,
"S(42564)": 0,
"S(39912)": 0,
"S(441)": 0,
"S(16582)": 0,
"CH2O2(442)": 0,
"S(443)": 0,
"S(444)": 0,
"S(445)": 0,
"S(446)": 0,
"S(447)": 0,
"S(448)": 0,
"S(449)": 0,
"S(450)": 0,
"S(451)": 0,
"S(453)": 0,
"S(454)": 0,
"S(455)": 0,
"S(43250)": 0,
"S(16587)": 0,
"S(456)": 0,
"S(457)": 0,
"S(458)": 0,
"S(43209)": 0,
"S(460)": 0,
"S(714)": 0,
"C4H9O(461)": 0,
"S(462)": 0,
"S(463)": 0,
"S(43235)": 0,
"S(27195)": 0,
"C2HO(467)": 0,
"C2HO(468)": 0,
"S(469)": 0,
"S(470)": 0,
"S(39920)": 0,
"S(471)": 0,
"C3H7O(473)": 0,
"CO3t2(474)": 0,
"CO3t1(475)": 0,
"CHO3(476)": 0,
"CH3O2(478)": 0,
"CH3O2(479)": 0,
"S(480)": 0,
"C2H2O(481)": 0,
"S(39051)": 0,
"S(44093)": 0,
"S(482)": 0,
"C2H4O(483)": 0,
"S(8475)": 0,
"S(14960)": 0,
"C2H4O(484)": 0,
"S(45761)": 0,
"S(45625)": 0,
"C7H9(3)": 0,
"C7H10(4)": 0,
"C7H11(5)": 0,
"C5H7O2(6)": 0,
"C6H12O(7)": 0,
"C7H9O2(8)": 0,
"C7H9O(9)": 0,
"C7H9O2(2)": 0,
"C7H8(11)": 0,
"S(12)": 0,
"C7H8O(1)": 0,
"C7H8O(2)": 0,
"C7H9O4(15)": 0,
"C7H8O3(16)": 0,
"S(17)": 0,
"S(18)": 0,
"C7H10(950)": 0,
"C7H10(953)": 0,
"C7H9(937)": 0,
"C7H11(954)": 0,
"S(1027)": 0,
"S(1107)": 0,
"S(1026)": 0,
"S(1182)": 0,
"S(1136)": 0,
"S(1423)": 0,
"C7H9(1021)": 0,
"S(1395)": 0,
"S(1147)": 0,
"C7H8(1072)": 0,
"S(1189)": 0,
"S(1660)": 0,
"S(2171)": 0,
"S(2308)": 0,
"S(2380)": 0,
"S(2504)": 0,
"S(2317)": 0,
"S(2325)": 0,
"S(1122)": 0,
"S(1558)": 0,
"S(2823)": 0,
"C5H6(960)": 0,
"S(2292)": 0,
"C5H6(2966)": 0,
"C5H6(118)": 0,
"S(2795)": 0,
"S(3354)": 0,
"S(3382)": 0,
"S(1055)": 0,
"S(1190)": 0,
"S(2980)": 0,
"S(2999)": 0,
"S(4179)": 0,
"S(4177)": 0,
"S(4178)": 0,
"S(2929)": 0,
"S(4584)": 0,
"S(4586)": 0,
"S(4585)": 0,
"S(2907)": 0,
"S(4183)": 0,
"S(4181)": 0,
"S(4182)": 0,
"S(2408)": 0,
"S(5142)": 0,
"S(5164)": 0,
"C5H7(596)": 0,
"S(3269)": 0,
"C7H9(1521)": 0,
"C5H5(98)": 0,
"S(5256)": 0,
"S(5475)": 0,
"S(5116)": 0,
"S(4588)": 0,
"S(3423)": 0,
"S(5668)": 0,
"S(3863)": 0,
"S(5102)": 0,
"S(5729)": 0,
"S(5680)": 0,
"S(3425)": 0,
"S(5535)": 0,
"S(6339)": 0,
"S(6025)": 0,
"S(1590)": 0,
"S(6567)": 0,
"S(6576)": 0,
"S(6568)": 0,
"S(6577)": 0,
"S(6673)": 0,
"S(6575)": 0,
"S(6463)": 0,
"S(6441)": 0,
"S(7159)": 0,
"S(7158)": 0,
"S(6769)": 0,
"S(6460)": 0,
"S(7425)": 0,
"S(7426)": 0,
"S(6844)": 0,
"S(6794)": 0,
"C6H7(363)": 0,
"C5H5O(143)": 0,
"S(8088)": 0,
"S(8087)": 0,
"S(7883)": 0,
"S(7884)": 0,
"C6H7(7867)": 0,
"C6H6(100)": 0,
"S(8095)": 0,
"C5H4O(142)": 0,
"S(7443)": 0,
"S(3596)": 0,
"C6H6(7863)": 0,
"S(8171)": 0,
"S(2953)": 0,
"S(2954)": 0,
"S(4610)": 0,
"S(3429)": 0,
"S(3385)": 0,
"S(8779)": 0,
"S(8926)": 0,
"S(8925)": 0,
"S(3570)": 0,
"S(8137)": 0,
"S(5511)": 0,
"S(1604)": 0,
"S(9457)": 0,
"S(9465)": 0,
"S(1209)": 0,
"S(9656)": 0,
"S(10145)": 0,
"S(9964)": 0,
"S(10199)": 0,
"S(10144)": 0,
"S(7901)": 0,
"S(9472)": 0,
"S(8372)": 0,
"S(10568)": 0,
"S(10726)": 0,
"S(1585)": 0,
"S(1397)": 0,
"S(10864)": 0,
"S(11105)": 0,
"S(11104)": 0,
"S(11116)": 0,
"S(11125)": 0,
"S(11364)": 0,
"S(11373)": 0,
"S(10567)": 0,
"S(11278)": 0,
"S(11483)": 0,
"S(10210)": 0,
"S(11446)": 0,
"S(10576)": 0,
"S(11847)": 0,
"S(6034)": 0,
"S(3657)": 0,
"S(10227)": 0,
"S(11846)": 0,
"S(10230)": 0,
"S(3569)": 0,
"S(12205)": 0,
"S(9078)": 0,
"S(9433)": 0,
"S(3597)": 0,
"S(2819)": 0,
"S(6674)": 0,
"S(12794)": 0,
"S(12545)": 0,
"S(9125)": 0,
"S(12558)": 0,
"S(6317)": 0,
"S(12938)": 0,
"S(12937)": 0,
"S(12555)": 0,
"S(12939)": 0,
"S(3652)": 0,
"S(9900)": 0,
"S(8981)": 0,
"S(1496)": 0,
"S(13750)": 0,
"S(13818)": 0,
"S(13927)": 0,
"S(13817)": 0,
"S(13868)": 0,
"S(13736)": 0,
"S(3386)": 0,
"S(4224)": 0,
"S(13821)": 0,
"S(13723)": 0,
"S(13877)": 0,
"S(14214)": 0,
"S(14206)": 0,
"S(14213)": 0,
"S(14837)": 0,
"S(13344)": 0,
"S(11447)": 0,
"S(14679)": 0,
"S(13770)": 0,
"S(14182)": 0,
"S(4240)": 0,
"S(6855)": 0,
"S(14951)": 0,
"S(3799)": 0,
"S(14950)": 0,
"S(14952)": 0,
"S(12935)": 0,
"S(15164)": 0,
"S(15738)": 0,
"S(15736)": 0,
"S(14682)": 0,
"S(394)": 0,
"S(14342)": 0,
"S(1184)": 0,
"S(6893)": 0,
"S(14886)": 0,
"S(14373)": 0,
"S(14327)": 0,
"S(12924)": 0,
"S(13833)": 0,
"S(16697)": 0,
"S(1649)": 0,
"S(17376)": 0,
"S(17377)": 0,
"S(1653)": 0,
"S(14370)": 0,
"S(16896)": 0,
"S(17375)": 0,
"S(16286)": 0,
"S(16253)": 0,
"S(19452)": 0,
"S(19451)": 0,
"S(19568)": 0,
"S(19454)": 0,
"S(20155)": 0,
"S(20203)": 0,
"S(20445)": 0,
"S(20818)": 0,
"S(5811)": 0,
"S(20947)": 0,
"S(20771)": 0,
"S(20811)": 0,
"S(16854)": 0,
"S(20232)": 0,
"S(20813)": 0,
"S(20094)": 0,
"S(5810)": 0,
"S(21071)": 0,
"S(21064)": 0,
"S(19487)": 0,
"S(20106)": 0,
"S(21215)": 0,
"S(21213)": 0,
"S(21190)": 0,
"C4H2(86)": 0,
"nC4H3(87)": 0,
"iC4H3(88)": 0,
"S(22782)": 0,
"nC4H5(91)": 0,
"C4H7(94)": 0,
"S(23384)": 0,
"S(20164)": 0,
"S(23383)": 0,
"S(21631)": 0,
"S(23699)": 0,
"S(22863)": 0,
"S(23382)": 0,
"S(21148)": 0,
"C4H612(97)": 0,
"C6H5(99)": 0,
"C7H7(101)": 0,
"C7H8(102)": 0,
"C6H5O(138)": 0,
"S(583)": 0,
"S(24214)": 0,
"S(24003)": 0,
"S(24213)": 0,
"S(23997)": 0,
"S(24545)": 0,
"S(588)": 0,
"S(24544)": 0,
"S(24044)": 0,
"S(24609)": 0,
"S(24558)": 0,
"S(24566)": 0,
"S(24561)": 0,
"S(24923)": 0,
"S(25021)": 0,
"S(24565)": 0,
"S(25016)": 0,
"S(25096)": 0,
"S(25092)": 0,
"S(141)": 0,
"S(24928)": 0,
"iC4H7(103)": 0,
"cC3H4(104)": 0,
"iC4H8(106)": 0,
"S(24603)": 0,
"S(25709)": 0,
"S(25576)": 0,
"S(25710)": 0,
"C3H6O(108)": 0,
"S(25532)": 0,
"S(25550)": 0,
"S(25542)": 0,
"C3H3O(109)": 0,
"S(25873)": 0,
"S(24564)": 0,
"S(25881)": 0,
"S(26372)": 0,
"S(25994)": 0,
"S(26654)": 0,
"S(26715)": 0,
"S(25129)": 0,
"S(26369)": 0,
"S(26365)": 0,
"S(26367)": 0,
"S(26371)": 0,
"H2C4O(111)": 0,
"C6H2(112)": 0,
"C6H3(113)": 0,
"C6H4(114)": 0,
"C6H4(115)": 0,
"C4H5O(116)": 0,
"C4H5(117)": 0,
"C4H5O(119)": 0,
"S(25788)": 0,
"S(26366)": 0,
"C4H6O(120)": 0,
"S(27339)": 0,
"S(26686)": 0,
"S(26660)": 0,
"S(25852)": 0,
"C4H6O(121)": 0,
"C4H6(122)": 0,
"C4H6O(123)": 0,
"C4H6O(124)": 0,
"C4H4O(125)": 0,
"C4H82(128)": 0,
"S(27643)": 0,
"S(27644)": 0,
"S(582)": 0,
"S(25947)": 0,
"iC4H9(129)": 0,
"C4H10(131)": 0,
"tC4H9(132)": 0,
"C7H7(133)": 0,
"C7H7O(134)": 0,
"S(27807)": 0,
"S(27812)": 0,
"S(27789)": 0,
"S(27314)": 0,
"S(28121)": 0,
"S(28063)": 0,
"S(28079)": 0,
"S(28120)": 0,
"S(28324)": 0,
"S(28266)": 0,
"S(28321)": 0,
"S(28658)": 0,
"S(29234)": 0,
"S(29272)": 0,
"S(29273)": 0,
"S(29265)": 0,
"S(29377)": 0,
"S(28268)": 0,
"S(29232)": 0,
"S(28260)": 0,
"C7H6O(135)": 0,
"S(24207)": 0,
"C7H8O(136)": 0,
"C6H6O(137)": 0,
"C7H8O(139)": 0,
"C7H5O(140)": 0,
"C5H5O(144)": 0,
"C5H5O(145)": 0,
"S(29864)": 0,
"S(29422)": 0,
"S(29421)": 0,
"S(30089)": 0,
"S(17259)": 0,
"S(30129)": 0,
"C5H5O(592)": 0,
"S(28329)": 0,
"S(30866)": 0,
"C5H6O(146)": 0,
"C4H5(147)": 0,
"C5H9(148)": 0,
"S(27152)": 0,
"cC5H9(151)": 0,
"cC5H8(152)": 0,
"S(1097)": 0,
"C5H9(153)": 0,
"C5H8(154)": 0,
"C5H9(155)": 0,
"C5H9(156)": 0,
"S(30985)": 0,
"S(16923)": 0,
"C7H9(1019)": 0,
"S(7919)": 0,
"S(31426)": 0,
"C6H11(161)": 0,
"C6H10(162)": 0,
"C6H11(164)": 0,
"C6H11(165)": 0,
"C6H11(166)": 0,
"S(23905)": 0,
"S(23820)": 0,
"S(23795)": 0,
"S(23807)": 0,
"S(31869)": 0,
"S(32442)": 0,
"S(1123)": 0,
"S(32443)": 0,
"C7H13(171)": 0,
"S(30079)": 0,
"S(9973)": 0,
"S(31888)": 0,
"S(9915)": 0,
"S(33011)": 0,
"C7H9(1525)": 0,
"C7H9(1939)": 0,
"C8H15(178)": 0,
"C8H16(179)": 0,
"C8H17(180)": 0,
"C8H17(181)": 0,
"C8H17(182)": 0,
"C8H17(183)": 0,
"S(33085)": 0,
"S(33029)": 0,
"S(33030)": 0,
"S(2835)": 0,
"S(5985)": 0,
"S(33038)": 0,
"S(33549)": 0,
"S(33087)": 0,
"C8H18(184)": 0,
"S(33550)": 0,
"C9H17(185)": 0,
"C9H18(186)": 0,
"C7H9(8388)": 0,
"S(33322)": 0,
"S(34440)": 0,
"C9H19(187)": 0,
"C9H19(188)": 0,
"S(2570)": 0,
"S(33910)": 0,
"S(6403)": 0,
"S(33965)": 0,
"S(3379)": 0,
"S(34749)": 0,
"S(34760)": 0,
"S(34744)": 0,
"S(8001)": 0,
"S(35144)": 0,
"S(35143)": 0,
"S(34771)": 0,
"S(33332)": 0,
"C6H11(321)": 0,
"S(34706)": 0,
"S(33975)": 0,
"C9H19(189)": 0,
"C9H19(190)": 0,
"C9H19(191)": 0,
"S(28322)": 0,
"S(34185)": 0,
"S(27153)": 0,
"S(34187)": 0,
"C9H20(192)": 0,
"S(36048)": 0,
"S(36061)": 0,
"S(36040)": 0,
"S(36062)": 0,
"S(36256)": 0,
"S(35475)": 0,
"S(193)": 0,
"S(194)": 0,
"S(3333)": 0,
"S(34723)": 0,
"S(34003)": 0,
"S(7726)": 0,
"S(3535)": 0,
"S(195)": 0,
"S(196)": 0,
"S(197)": 0,
"S(28348)": 0,
"S(27170)": 0,
"S(36998)": 0,
"S(29930)": 0,
"S(198)": 0,
"S(199)": 0,
"S(200)": 0,
"S(201)": 0,
"S(202)": 0,
"S(30452)": 0,
"S(30103)": 0,
"S(203)": 0,
"S(204)": 0,
"S(205)": 0,
"S(206)": 0,
"S(207)": 0,
"S(208)": 0,
"S(209)": 0,
"S(210)": 0,
"S(211)": 0,
"S(212)": 0,
"S(213)": 0,
"S(214)": 0,
"S(215)": 0,
"S(216)": 0,
"S(217)": 0,
"S(218)": 0,
"S(219)": 0,
"S(220)": 0,
"S(221)": 0,
"S(20127)": 0,
"S(37603)": 0,
"S(37605)": 0,
"S(36500)": 0,
"C9H17(222)": 0,
"S(223)": 0,
"C8H15(224)": 0,
"C7H13(225)": 0,
"C6H11(226)": 0,
"S(37726)": 0,
"S(37730)": 0,
"S(37725)": 0,
"S(37775)": 0,
"S(29306)": 0,
"S(38213)": 0,
"S(37708)": 0,
"S(227)": 0,
"S(228)": 0,
"C7H13(230)": 0,
"C7H14(231)": 0,
"S(232)": 0,
"S(39999)": 0,
"S(12629)": 0,
"S(28649)": 0,
"S(39153)": 0,
"S(233)": 0,
"S(234)": 0,
"S(235)": 0,
"S(236)": 0,
"S(237)": 0,
"S(41373)": 0,
"S(41371)": 0,
"S(238)": 0,
"S(21432)": 0,
"S(2982)": 0,
"S(41366)": 0,
"S(41449)": 0,
"S(43582)": 0,
"S(239)": 0,
"C9H16(240)": 0,
"C8H14(241)": 0,
"S(242)": 0,
"C7H12(243)": 0,
"S(244)": 0,
"C6H10(245)": 0,
"S(246)": 0,
"S(3101)": 0,
"S(247)": 0,
"S(44328)": 0,
"S(44330)": 0,
"S(248)": 0,
"S(249)": 0,
"S(250)": 0,
"S(251)": 0,
"C8H15(252)": 0,
"S(253)": 0,
"C8H15(254)": 0,
"S(255)": 0,
"C8H15(256)": 0,
"S(257)": 0,
"C8H15(258)": 0,
"S(259)": 0,
"C7H12(260)": 0,
"S(261)": 0,
"C8H15(262)": 0,
"S(263)": 0,
"C6H10(264)": 0,
"C8H14(266)": 0,
"C8H14(267)": 0,
"S(1086)": 0,
"C7H12(268)": 0,
"S(1085)": 0,
"S(44543)": 0,
"S(1117)": 0,
"S(44547)": 0,
"S(44550)": 0,
"S(1115)": 0,
"S(44690)": 0,
"S(44884)": 0,
"S(44885)": 0,
"S(44670)": 0,
"S(44698)": 0,
"C9H18(269)": 0,
"C9H18(270)": 0,
"C9H18(271)": 0,
"C6H11(272)": 0,
"C6H12(273)": 0,
"C9H17(274)": 0,
"S(44641)": 0,
"C9H17(275)": 0,
"S(44655)": 0,
"S(44654)": 0,
"C9H17(276)": 0,
"C9H17(277)": 0,
"C9H17(278)": 0,
"C9H17(279)": 0,
"C9H17(280)": 0,
"C9H17(281)": 0,
"C9H17(282)": 0,
"C9H17(283)": 0,
"C9H17(284)": 0,
"C9H17(285)": 0,
"C9H17(286)": 0,
"C7H13(287)": 0,
"C9H17(288)": 0,
"C9H17(289)": 0,
"C7H13(290)": 0,
"S(1082)": 0,
"C9H17(291)": 0,
"C7H13(292)": 0,
"C9H17(293)": 0,
"C9H17(294)": 0,
"S(36444)": 0,
"C7H13(295)": 0,
"C7H9(958)": 0,
"C9H17(296)": 0,
"C9H17(297)": 0,
"C6H10(298)": 0,
"C8H16(299)": 0,
"C8H16(300)": 0,
"C8H16(301)": 0,
"C5H9(302)": 0,
"C5H10(303)": 0,
"C8H15(304)": 0,
"C8H15(305)": 0,
"S(45160)": 0,
"S(45036)": 0,
"S(45159)": 0,
"S(45037)": 0,
"C8H15(306)": 0,
"S(45280)": 0,
"S(45262)": 0,
"S(45279)": 0,
"C8H15(307)": 0,
"C8H15(308)": 0,
"C8H15(309)": 0,
"S(45205)": 0,
"S(45554)": 0,
"S(45553)": 0,
"S(45501)": 0,
"S(45480)": 0,
"S(45483)": 0,
"S(45260)": 0,
"C8H15(310)": 0,
"S(46144)": 0,
"S(5707)": 0,
"S(23031)": 0,
"C8H15(311)": 0,
"S(46722)": 0,
"S(46727)": 0,
"S(46152)": 0,
"S(46721)": 0,
"C8H15(312)": 0,
"C8H15(313)": 0,
"C8H15(314)": 0,
"C8H15(315)": 0,
"S(956)": 0,
"S(48836)": 0,
"S(48886)": 0,
"C8H15(316)": 0,
"S(48904)": 0,
"S(48913)": 0,
"S(48907)": 0,
"C6H11(317)": 0,
"C6H11(318)": 0,
"C6H11(319)": 0,
"C8H15(320)": 0,
"C8H15(322)": 0,
"S(1084)": 0,
"C8H15(323)": 0,
"C8H15(324)": 0,
"S(2983)": 0,
"S(3571)": 0,
"S(48983)": 0,
"S(48984)": 0,
"C8H15(325)": 0,
"C6H10(326)": 0,
"S(5118)": 0,
"S(15233)": 0,
"S(49070)": 0,
"S(49013)": 0,
"S(49198)": 0,
"S(49279)": 0,
"S(49278)": 0,
"S(49284)": 0,
"S(49295)": 0,
"S(5797)": 0,
"S(49552)": 0,
"S(49184)": 0,
"CO2(12610)": 0,
"C8H15(327)": 0,
"C7H15(329)": 0,
"C7H14(330)": 0,
"C7H13(333)": 0,
"S(49575)": 0,
"S(49950)": 0,
"S(49942)": 0,
"S(49557)": 0,
"S(49945)": 0,
"C7H13(334)": 0,
"C7H13(335)": 0,
"S(32085)": 0,
"S(50151)": 0,
"S(49944)": 0,
"C7H13(336)": 0,
"C7H13(337)": 0,
"C7H13(338)": 0,
"C7H13(339)": 0,
"C7H13(341)": 0,
"C7H13(342)": 0,
"C7H13(343)": 0,
"C7H13(344)": 0,
"C7H13(345)": 0,
"C5H9(346)": 0,
"C7H13(347)": 0,
"C5H8(348)": 0,
"C7H13(349)": 0,
"C7H13(350)": 0,
"C5H8(351)": 0,
"C6H12(352)": 0,
"S(353)": 0,
"S(354)": 0,
"S(355)": 0,
"S(356)": 0,
"C6H11(357)": 0,
"C6H11(358)": 0,
"C6H11(359)": 0,
"C6H9(360)": 0,
"C6H8(361)": 0,

    },
    terminationTime = (5.0, 's'),
    terminationRateRatio = 0.01,
    terminationConversion={
                'C7H16(1)': 0.99,
    },
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0.01, #  pruning to start
    toleranceMoveToCore=0.8,
    toleranceInterruptSimulation=1,
    maxNumObjsPerIter=2,      #
    terminateAtMaxObjects=True,
    maxNumSpecies=1550, # first stage is until core reaches 100 species
    filterReactions=True, # should speed up model generation
    filterThreshold=2e8,
)
model(
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=1e8,
    toleranceKeepInEdge=0.01, # Pruning enabled for stage 2
    maximumEdgeSpecies=100000,
    minCoreSizeForPrune=100,
    minSpeciesExistIterationsForPrune=2,
    filterReactions=True,
    filterThreshold=2e8,
    maxNumObjsPerIter=4,
    terminateAtMaxObjects=True,
    )

options(
    units='si',
    # saveRestartPeriod=(1,'hour'),
    generateOutputHTML=False,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
    generateSeedEachIteration=False,
)