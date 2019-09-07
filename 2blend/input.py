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
    reactive=True,
    structure=adjacencyList(
        """
1 Ar u0 p4 c0
        """),
)


species(
    label='He(9)',
    reactive=True,
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
    reactive=True,
    structure=adjacencyList(
        """
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
        """),
)


species(
    label='Ne',
    reactive=True,
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
    structure=adjac