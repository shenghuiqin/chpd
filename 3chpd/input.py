#Data sources
database(
    thermoLibraries =['primaryThermoLibrary'], # 'FFCM1(-)','primaryThermoLibrary', 'BurkeH2O2','DFT_QCI_thermo','CBS_QB3_1dHR'
    reactionLibraries = [], # 
    seedMechanisms = [], #
    kineticsDepositories = ['training'], 
    kineticsFamilies = ['Ketoenol','Intra_R_Add_Exocyclic','Cyclic_Ether_Formation'],#,'Cyclic_Ether_Formation','HO2_Elimination_from_PeroxyRadical','Intra_Disproportionation','intra_H_migration','ketoenol'
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumRadicalElectrons = 2,
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=9,
    maximumOxygenAtoms=4,
)

# List of species

species(
    label='N2(1)',
    reactive=True,
    structure=adjacencyList(
        """
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
        """),
)


species(
    label='Ar(2)',
    reactive=True,
    structure=adjacencyList(
        """
1 Ar u0 p4 c0
        """),
)


species(
    label='He(3)',
    reactive=True,
    structure=adjacencyList(
        """
1 He u0 p1 c0
        """),
)


species(
    label='Ne(4)',
    reactive=True,
    structure=adjacencyList(
        """
1 Ne u0 p4 c0
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
    label='HO2(3)',
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
    label='C7H9(5)',
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
    label='C7H10(9)',
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
    label='OO(8)',
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
    label='H(12)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """),
)


species(
    label='C7H9O2(20)',
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
    label='C7H9O2(19)',
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
    label='S(36)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,D} {16,S}
7  C u0 p0 c0 {5,S} {8,D} {15,S}
8  C u0 p0 c0 {7,D} {9,S} {18,S}
9  C u0 p0 c0 {6,D} {8,S} {17,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(35)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {19,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,D} {16,S}
7  C u0 p0 c0 {4,S} {8,D} {15,S}
8  C u0 p0 c0 {3,S} {7,D} {17,S}
9  C u0 p0 c0 {3,S} {6,D} {18,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O(53)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {13,S}
5  C u0 p0 c0 {3,S} {8,D} {15,S}
6  C u0 p0 c0 {2,S} {7,D} {14,S}
7  C u0 p0 c0 {4,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {5,D} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='OH(D)(44)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H10O(67)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {18,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {8,D} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u0 p0 c0 {5,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9O(69)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {8,S} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u1 p0 c0 {6,S} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9O3(82)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O3(83)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {8,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O3(84)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
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
    label='S(94)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(93)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {16,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {8,D} {9,S} {18,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(92)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {17,S}
9  C u0 p0 c0 {6,S} {10,D} {16,S}
10 C u0 p0 c0 {5,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(133)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {18,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {14,S}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {15,S}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(98)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {6,S} {19,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(97)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {6,S} {19,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {14,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {15,S}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(99)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {18,S}
2  O u0 p2 c0 {3,S} {19,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {9,D} {15,S}
7  C u0 p0 c0 {5,S} {8,D} {14,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {6,D} {8,S} {16,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H9O(45)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {3,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {8,D} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u0 p0 c0 {5,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H8O(75)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {4,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {8,S} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u1 p0 c0 {6,S} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='C7H8O3(85)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {6,S}
2  O u1 p2 c0 {7,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(101)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {10,S} {16,S}
8  C u0 p0 c0 {2,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(179)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {15,S}
6  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u1 p0 c0 {5,S} {10,S} {17,S}
9  C u0 p0 c0 {6,S} {10,D} {16,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(180)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {16,S}
8  C u0 p0 c0 {2,D} {6,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(108)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {4,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {7,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(183)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
5  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
7  C u1 p0 c0 {4,S} {5,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(126)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u1 p0 c0 {4,S} {9,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(134)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {19,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {9,D} {16,S}
8  C u0 p0 c0 {4,S} {10,D} {17,S}
9  C u0 p0 c0 {3,S} {5,S} {7,D}
10 C u0 p0 c0 {6,S} {8,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(221)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {15,S}
8  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {7,S} {11,D} {18,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(222)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {17,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(234)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {19,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
8  C u0 p0 c0 {3,D} {5,S} {7,S}
9  C u0 p0 c0 {6,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(226)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {19,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {10,D} {17,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(225)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {4,D} {6,S} {7,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(255)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u1 p2 c0 {5,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {3,D} {4,S} {6,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(245)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u1 p0 c0 {4,S} {10,S} {16,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(273)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {10,D}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {9,S} {11,S} {17,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
9  C u1 p0 c0 {5,S} {6,S} {18,S}
10 C u0 p0 c0 {3,D} {5,S} {7,S}
11 C u0 p0 c0 {4,D} {6,S} {8,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(184)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {6,D}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
5  C u0 p0 c0 {3,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {2,D} {3,S} {4,S}
7  C u1 p0 c0 {4,S} {9,S} {16,S}
8  C u0 p0 c0 {5,S} {9,D} {15,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(191)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {17,S}
9  C u0 p0 c0 {3,D} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(188)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
7  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {3,D} {5,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {19,S}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(151)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {17,S}
2  O u0 p2 c0 {5,S} {18,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {4,S} {8,D}
6  C u0 p0 c0 {1,S} {3,S} {7,D}
7  C u0 p0 c0 {6,D} {9,S} {14,S}
8  C u0 p0 c0 {5,D} {9,S} {16,S}
9  C u1 p0 c0 {7,S} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(279)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {20,S}
2  O u0 p2 c0 {6,S} {21,S}
3  O u0 p2 c0 {10,D}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {9,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {3,D} {5,S} {8,S}
11 C u0 p0 c0 {4,D} {6,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(312)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {1,S} {21,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(310)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {18,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(311)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {1,S} {21,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {11,D}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 C u0 p0 c0 {7,S} {8,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(309)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {11,D}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 C u0 p0 c0 {7,S} {8,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(330)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,S} {10,S} {15,S}
7  C u0 p0 c0 {2,S} {5,S} {10,D}
8  C u0 p0 c0 {1,S} {4,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {16,S}
10 C u0 p0 c0 {6,S} {7,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(314)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {7,S} {20,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(315)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(313)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {7,S} {20,S}
4  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {15,S}
7  C u0 p0 c0 {3,S} {5,S} {10,D}
8  C u0 p0 c0 {2,S} {4,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {16,S}
10 C u0 p0 c0 {6,S} {7,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(301)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {17,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {9,S} {14,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u1 p0 c0 {7,S} {8,S} {15,S}
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
    label='S(362)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {11,D}
10 C u1 p0 c0 {7,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(363)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {3,D} {7,S} {11,S}
10 C u0 p0 c0 {6,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(335)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {10,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {16,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {7,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(324)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {5,S} {11,S} {16,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {9,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(317)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {5,S} {10,S} {16,S}
9  C u0 p0 c0 {4,D} {7,S} {11,S}
10 C u0 p0 c0 {8,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(366)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
6  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {4,D} {7,S} {11,S}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(323)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {6,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {6,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(331)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {11,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {4,S} {7,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(365)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
6  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {4,D} {5,S} {7,S}
9  C u0 p0 c0 {2,S} {6,S} {11,D}
10 C u1 p0 c0 {5,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(355)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {18,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {10,S} {15,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {16,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(316)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
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
    label='S(339)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {18,S}
2  O u1 p2 c0 {6,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {15,S}
7  C u0 p0 c0 {1,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {16,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(404)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {10,S} {20,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {7,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {5,S} {11,S} {17,S}
9  C u0 p0 c0 {4,D} {5,S} {7,S}
10 C u0 p0 c0 {2,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(401)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {7,S} {20,S}
3  O u1 p2 c0 {5,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,D} {17,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(274)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {20,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {14,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(388)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {8,S} {18,S}
2  O u1 p2 c0 {6,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {5,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {15,S}
7  C u0 p0 c0 {3,D} {5,S} {6,S}
8  C u0 p0 c0 {1,S} {4,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(427)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {16,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(426)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {10,S} {20,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {17,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {4,D} {6,S} {8,S}
10 C u0 p0 c0 {2,S} {7,S} {11,D}
11 C u0 p0 c0 {5,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(436)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u1 p2 c0 {10,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {16,S}
8  C u0 p0 c0 {2,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {4,S} {6,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H11(85)',
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
    label='C7H11(86)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  C u0 p0 c0 {3,S} {7,D} {17,S}
7  C u0 p0 c0 {4,S} {6,D} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(111)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {6,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {9,D} {20,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {8,S}
        """),
)


species(
    label='C7H8(110)',
    reactive=True,
    structure=adjacencyList(
        """
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,D} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {2,D} {6,S} {12,S}
5  C u0 p0 c0 {3,D} {7,S} {15,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(150)',
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
    label='S(113)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {7,S} {9,D} {19,S}
9  C u0 p0 c0 {6,S} {8,D} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='H2O(204)',
    reactive=True,
    structure=adjacencyList(
        """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9O(175)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {17,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
4  C u0 p0 c0 {2,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {6,D} {14,S}
6  C u0 p0 c0 {3,S} {5,D} {15,S}
7  C u0 p0 c0 {3,S} {8,D} {16,S}
8  C u1 p0 c0 {4,S} {7,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(256)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {19,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {9,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {18,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(334)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {1,S} {19,S}
4  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {9,D}
7  C u1 p0 c0 {5,S} {10,S} {15,S}
8  C u0 p0 c0 {2,S} {9,S} {10,D}
9  C u0 p0 c0 {6,D} {8,S} {17,S}
10 C u0 p0 c0 {7,S} {8,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(349)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {17,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u1 p0 c0 {4,S} {9,S} {14,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u0 p0 c0 {1,S} {8,S} {9,D}
8  C u0 p0 c0 {6,D} {7,S} {16,S}
9  C u0 p0 c0 {5,S} {7,D} {15,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(425)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,D}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {8,D} {14,S}
6  C u0 p0 c0 {1,D} {3,S} {7,S}
7  C u0 p0 c0 {6,S} {9,D} {16,S}
8  C u0 p0 c0 {5,D} {9,S} {15,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(418)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {10,S} {19,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u1 p0 c0 {6,S} {10,S} {17,S}
9  C u0 p0 c0 {3,D} {7,S} {11,S}
10 C u0 p0 c0 {2,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(420)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
6  C u0 p0 c0 {7,S} {11,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
8  C u0 p0 c0 {5,S} {11,D} {18,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 C u0 p0 c0 {3,S} {6,S} {8,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(493)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,D}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {2,D} {5,S} {6,S}
9  C u0 p0 c0 {7,S} {11,D} {17,S}
10 C u0 p0 c0 {3,D} {6,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(480)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {8,S} {10,S} {14,S} {15,S}
8  C u0 p0 c0 {3,D} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {11,S} {17,S}
10 C u0 p0 c0 {7,S} {11,D} {16,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(164)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {5,S}
2  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {3,S} {7,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {8,D} {19,S}
8  C u0 p0 c0 {6,S} {7,D} {18,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(524)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {19,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
6  C u1 p0 c0 {4,S} {5,S} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {18,S}
8  C u0 p0 c0 {5,S} {7,D} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(553)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,S} {21,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {14,S}
8  C u0 p0 c0 {4,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {7,S} {10,D} {20,S}
10 C u0 p0 c0 {8,S} {9,D} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='C7H8O(215)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,D}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,D} {2,S} {6,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u0 p0 c0 {6,D} {7,S} {15,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(612)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {6,S} {10,D}
9  C u1 p0 c0 {7,S} {10,S} {18,S}
10 C u0 p0 c0 {8,D} {9,S} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(555)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,S} {21,S}
3  O u0 p2 c0 {1,S} {22,S}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {14,S}
8  C u0 p0 c0 {4,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {7,S} {10,D} {20,S}
10 C u0 p0 c0 {8,S} {9,D} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(635)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {20,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,D} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {18,S}
10 C u0 p0 c0 {8,S} {9,D} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(661)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {21,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {10,S} {18,S} {19,S}
9  C u0 p0 c0 {2,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {20,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(577)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {1,S} {19,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {6,D} {10,S} {15,S}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u1 p0 c0 {8,S} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(519)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {19,S}
2  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u1 p0 c0 {5,S} {8,S} {17,S}
8  C u0 p0 c0 {6,D} {7,S} {18,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(766)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u1 p2 c0 {6,S}
2  O u0 p2 c0 {5,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {2,D} {3,S} {7,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {6,D} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(532)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,D}
2  C u0 p0 c0 {3,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,D} {4,S} {8,S}
7  C u0 p0 c0 {5,S} {8,D} {17,S}
8  C u0 p0 c0 {6,S} {7,D} {18,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(770)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {9,S} {21,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {17,S}
8  C u0 p0 c0 {6,S} {9,S} {18,S} {19,S}
9  C u0 p0 c0 {2,S} {8,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {20,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(771)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {21,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {7,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {4,S} {10,D} {20,S}
10 C u0 p0 c0 {8,S} {9,D} {19,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(869)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,D}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {7,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {2,D} {5,S} {6,S}
9  C u0 p0 c0 {3,D} {7,S} {11,S}
10 C u0 p0 c0 {6,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(495)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {10,D}
8  C u0 p0 c0 {3,D} {6,S} {11,S}
9  C u0 p0 c0 {2,S} {10,S} {11,D}
10 C u0 p0 c0 {7,D} {9,S} {16,S}
11 C u0 p0 c0 {8,S} {9,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(890)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {9,S} {17,S}
3  O u0 p2 c0 {1,S} {18,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {7,D} {10,S}
7  C u0 p0 c0 {5,S} {6,D} {14,S}
8  C u0 p0 c0 {4,S} {5,S} {11,D}
9  C u0 p0 c0 {2,S} {10,D} {11,S}
10 C u0 p0 c0 {6,S} {9,D} {15,S}
11 C u0 p0 c0 {8,D} {9,S} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(899)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {7,S} {16,S}
2  O u1 p2 c0 {6,S}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,D} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {1,S} {8,S} {10,D}
8  C u0 p0 c0 {6,D} {7,S} {15,S}
9  C u0 p0 c0 {3,S} {5,D} {10,S}
10 C u0 p0 c0 {7,D} {9,S} {14,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(931)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {16,S}
3  O u0 p2 c0 {1,S} {17,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,D} {8,S}
6  C u0 p0 c0 {2,S} {7,S} {9,D}
7  C u0 p0 c0 {5,D} {6,S} {14,S}
8  C u0 p0 c0 {5,S} {10,D} {13,S}
9  C u0 p0 c0 {6,D} {11,S} {15,S}
10 C u0 p0 c0 {8,D} {11,S} {12,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(950)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {15,S}
2  O u1 p2 c0 {7,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {13,S}
6  C u0 p0 c0 {4,D} {10,S} {14,S}
7  C u0 p0 c0 {2,S} {5,D} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {12,S}
9  C u0 p0 c0 {8,D} {10,S} {11,S}
10 C u0 p0 c0 {3,D} {6,S} {9,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='C7H9(93)',
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
    label='S(985)',
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
7  C u0 p0 c0 {4,S} {8,D} {15,S}
8  C u0 p0 c0 {5,S} {7,D} {16,S}
9  C u0 p0 c0 {5,S} {6,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(151)',
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
    label='S(140)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,D} {14,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {9,S} {15,S}
8  C u0 p0 c0 {5,D} {9,S} {16,S}
9  C u1 p0 c0 {7,S} {8,S} {17,S}
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
    label='S(230)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {10,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {18,S}
9  C u0 p0 c0 {6,S} {11,D} {19,S}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H8O(212)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {7,D} {15,S}
6  C u0 p0 c0 {3,S} {8,D} {16,S}
7  C u0 p0 c0 {4,S} {5,D} {14,S}
8  C u0 p0 c0 {4,S} {6,D} {13,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
        """),
)


species(
    label='S(231)',
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
    label='C7H8O(222)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,D} {14,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {8,S} {15,S}
8  C u0 p0 c0 {5,D} {7,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(191)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
7  C u1 p0 c0 {6,S} {9,S} {17,S}
8  C u0 p0 c0 {5,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(177)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
7  C u1 p0 c0 {5,S} {6,S} {17,S}
8  C u0 p0 c0 {3,S} {9,D} {19,S}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(189)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u0 p2 c0 {1,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {6,S} {17,S}
8  C u0 p0 c0 {5,S} {9,D} {19,S}
9  C u0 p0 c0 {6,S} {8,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(235)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {11,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {6,S} {9,D} {19,S}
9  C u0 p0 c0 {7,S} {8,D} {20,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(419)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {17,S}
9  C u0 p0 c0 {5,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {8,S} {11,D} {21,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(239)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {10,S} {17,S}
9  C u0 p0 c0 {7,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {8,S} {11,D} {21,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(178)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {20,S}
3  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u1 p0 c0 {4,S} {9,S} {18,S}
8  C u0 p0 c0 {6,S} {9,D} {17,S}
9  C u0 p0 c0 {7,S} {8,D} {19,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(238)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {9,S} {12,S} {18,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {19,S}
10 C u0 p0 c0 {5,S} {11,D} {20,S}
11 C u0 p0 c0 {9,S} {10,D} {21,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(683)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {7,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {18,S}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(237)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {15,S}
8  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
9  C u0 p0 c0 {5,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {7,S} {11,D} {20,S}
11 C u0 p0 c0 {9,S} {10,D} {21,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {11,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(686)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {2,S} {5,S} {10,S} {17,S}
8  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {7,S} {11,D} {21,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(518)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
5  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {2,D} {3,S} {5,S}
8  C u0 p0 c0 {4,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(401)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {8,D} {18,S}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(420)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {6,S} {11,S} {18,S} {19,S}
10 C u0 p0 c0 {8,S} {11,D} {21,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='C7H11(292)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {7,S} {10,S} {15,S}
5  C u1 p0 c0 {1,S} {2,S} {16,S}
6  C u0 p0 c0 {3,S} {7,D} {17,S}
7  C u0 p0 c0 {4,S} {6,D} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(556)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {4,S}
2  O u0 p2 c0 {1,S} {18,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {4,S} {10,S} {16,S}
8  C u0 p0 c0 {5,S} {6,D} {14,S}
9  C u0 p0 c0 {5,S} {10,D} {15,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
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
    label='C7H9O(286)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {17,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {14,S}
6  C u0 p0 c0 {4,S} {8,D} {15,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u1 p0 c0 {6,D} {7,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(656)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {2,S} {21,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {9,S} {10,S} {15,S}
8  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {6,S} {7,S} {18,S}
10 C u0 p0 c0 {7,S} {11,D} {20,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1074)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {2,S} {21,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {9,S} {11,S} {17,S} {18,S}
9  C u1 p0 c0 {2,S} {6,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {20,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(661)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {1,S} {22,S}
4  O u0 p2 c0 {2,S} {21,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {17,S}
9  C u1 p0 c0 {5,S} {11,S} {18,S}
10 C u0 p0 c0 {8,S} {11,D} {19,S}
11 C u0 p0 c0 {9,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1299)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {3,S} {6,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {10,D} {18,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {8,D} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1287)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {20,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {15,S}
7  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {6,S} {10,D} {19,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1298)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {3,S} {9,S}
3  O u1 p2 c0 {2,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {5,S} {10,D} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {17,S}
9  C u0 p0 c0 {2,S} {8,D} {10,S}
10 C u0 p0 c0 {7,D} {9,S} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1347)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {3,S} {4,S}
3  O u0 p2 c0 {2,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {17,S}
7  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {5,S} {10,D} {18,S}
10 C u0 p0 c0 {6,S} {9,D} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1346)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {2,S} {20,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {17,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {9,S} {16,S}
9  C u0 p0 c0 {8,S} {10,D} {18,S}
10 C u0 p0 c0 {5,S} {9,D} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1382)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {3,S} {8,S}
3  O u0 p2 c0 {2,S} {19,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u1 p0 c0 {4,S} {10,S} {15,S}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {2,S} {9,S} {10,D}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 C u0 p0 c0 {6,S} {8,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1364)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1493)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {4,S} {18,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {3,S} {9,S} {14,S}
7  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {7,S} {9,D} {17,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1457)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {17,S}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {8,S} {14,S}
6  C u0 p0 c0 {4,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {9,S} {15,S}
8  C u0 p0 c0 {5,S} {9,D} {16,S}
9  C u0 p0 c0 {2,S} {7,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1492)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {3,S} {18,S}
3  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u1 p0 c0 {6,S} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1586)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {18,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1584)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {19,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {15,S}
8  C u0 p0 c0 {6,S} {10,D} {16,S}
9  C u0 p0 c0 {7,S} {11,D} {17,S}
10 C u0 p0 c0 {8,D} {11,S} {18,S}
11 C u0 p0 c0 {3,S} {9,D} {10,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1587)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {20,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {3,S} {5,S} {10,S}
7  C u0 p0 c0 {2,S} {8,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {6,S} {11,D} {19,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1577)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {3,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {8,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u1 p0 c0 {3,S} {5,S} {16,S}
9  C u0 p0 c0 {6,S} {7,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1581)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {11,D} {17,S}
9  C u0 p0 c0 {7,S} {10,D} {16,S}
10 C u0 p0 c0 {9,D} {11,S} {18,S}
11 C u0 p0 c0 {3,S} {8,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1630)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {3,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {11,S} {15,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {16,S}
9  C u0 p0 c0 {1,S} {8,S} {11,D}
10 C u1 p0 c0 {5,S} {8,S} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1628)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {1,S} {5,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1678)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {7,S} {13,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {16,S}
8  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {8,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {1,S} {7,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1749)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {7,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {14,S}
7  C u0 p0 c0 {3,S} {5,S} {10,S} {13,S}
8  C u0 p0 c0 {10,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {7,S} {8,S} {17,S}
11 C u0 p0 c0 {8,S} {9,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='O(T)(511)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1 O u2 p2 c0
        """),
)


species(
    label='S(1734)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {4,S} {18,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {11,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
7  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
8  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1185)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {6,S} {10,D}
9  C u0 p0 c0 {3,D} {6,S} {7,S}
10 C u0 p0 c0 {4,S} {8,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(500)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {6,S} {18,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
5  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  C u0 p0 c0 {4,S} {9,D} {16,S}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1052)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {19,S}
2  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
6  C u1 p0 c0 {3,S} {5,S} {17,S}
7  C u0 p0 c0 {1,S} {5,S} {8,D}
8  C u0 p0 c0 {4,S} {7,D} {18,S}
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
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(994)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {17,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u1 p0 c0 {3,S} {4,S} {13,S}
6  C u1 p0 c0 {3,S} {9,S} {14,S}
7  C u0 p0 c0 {2,D} {4,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1057)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {11,S}
5  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {9,D} {15,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {16,S}
9  C u0 p0 c0 {5,S} {6,D} {14,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(1122)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {6,S} {7,S} {16,S}
9  C u0 p0 c0 {3,D} {5,S} {7,S}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {6,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1854)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {8,S} {10,S} {14,S} {15,S}
8  C u1 p0 c0 {5,S} {7,S} {16,S}
9  C u0 p0 c0 {6,S} {11,D} {17,S}
10 C u0 p0 c0 {3,D} {7,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1717)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {6,S} {10,S} {16,S}
9  C u0 p0 c0 {3,D} {7,S} {11,S}
10 C u0 p0 c0 {8,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1108)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,D} {4,S} {7,S}
9  C u0 p0 c0 {2,S} {6,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1930)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {8,S} {18,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u0 p0 c0 {7,D} {11,S} {15,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u1 p0 c0 {9,S} {10,S} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1143)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {11,S}
6  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1908)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {6,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {4,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {10,D} {16,S}
9  C u0 p0 c0 {3,D} {7,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1852)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {19,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {14,S}
8  C u0 p0 c0 {2,S} {5,S} {10,D}
9  C u0 p0 c0 {3,S} {6,S} {11,D}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1962)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {1,S} {19,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {5,S} {10,D}
9  C u0 p0 c0 {7,D} {11,S} {17,S}
10 C u0 p0 c0 {8,D} {11,S} {15,S}
11 C u1 p0 c0 {9,S} {10,S} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1953)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {17,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {4,S} {9,D}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u0 p0 c0 {5,S} {7,D} {15,S}
10 C u0 p0 c0 {5,S} {8,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1965)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {17,S}
3  O u0 p2 c0 {7,S} {18,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {14,S}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {15,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1997)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {7,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
6  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u0 p0 c0 {2,S} {10,S} {11,D}
9  C u0 p0 c0 {5,S} {7,D} {16,S}
10 C u1 p0 c0 {5,S} {8,S} {17,S}
11 C u0 p0 c0 {6,S} {8,D} {15,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1402)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {15,S}
8  C u1 p0 c0 {4,S} {10,S} {16,S}
9  C u0 p0 c0 {5,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {8,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1579)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {4,S} {10,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {16,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {6,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2064)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {5,S} {17,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {9,S} {11,S}
6  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {4,S} {10,D}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {16,S}
10 C u0 p0 c0 {6,S} {7,D} {15,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2072)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {10,D}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {15,S}
10 C u0 p0 c0 {6,S} {7,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2094)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u0 p2 c0 {3,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {3,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2095)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {4,S} {10,S}
4  O u0 p2 c0 {3,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
7  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
8  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
9  C u1 p0 c0 {5,S} {7,S} {17,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2093)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {4,S} {10,S}
4  O u0 p2 c0 {3,S} {20,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {15,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {1,S} {6,S} {7,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1624)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {3,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {3,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {9,D}
8  C u0 p0 c0 {6,S} {9,D} {16,S}
9  C u0 p0 c0 {7,D} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2143)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {8,S}
2  O u0 p2 c0 {4,S} {18,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {2,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {14,S}
6  C u0 p0 c0 {4,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
8  C u1 p0 c0 {1,S} {5,S} {6,S}
9  C u1 p0 c0 {7,S} {10,S} {17,S}
10 C u0 p0 c0 {3,D} {5,S} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2124)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {8,S} {9,S}
3  O u0 p2 c0 {5,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {15,S}
7  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
8  C u0 p0 c0 {2,S} {5,S} {10,S} {16,S}
9  C u0 p0 c0 {2,S} {6,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2135)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {6,S} {10,S} {14,S} {15,S}
8  C u1 p0 c0 {4,S} {6,S} {16,S}
9  C u0 p0 c0 {3,D} {5,S} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2117)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {15,S}
6  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
8  C u1 p0 c0 {7,S} {10,S} {16,S}
9  C u0 p0 c0 {3,D} {5,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2186)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {13,S}
7  C u0 p0 c0 {9,S} {10,S} {14,S} {15,S}
8  C u1 p0 c0 {4,S} {6,S} {17,S}
9  C u0 p0 c0 {3,D} {5,S} {7,S}
10 C u1 p0 c0 {6,S} {7,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2170)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {4,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {15,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {16,S}
8  C u1 p0 c0 {1,S} {6,S} {10,S}
9  C u0 p0 c0 {3,D} {6,S} {7,S}
10 C u1 p0 c0 {4,S} {8,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2200)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {7,S} {17,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {5,S} {10,D}
8  C u1 p0 c0 {2,S} {6,S} {9,S}
9  C u1 p0 c0 {8,S} {10,S} {15,S}
10 C u0 p0 c0 {7,D} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2192)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {13,S}
7  C u1 p0 c0 {4,S} {6,S} {14,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {15,S}
10 C u0 p0 c0 {8,D} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1036)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {19,S}
2  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
6  C u1 p0 c0 {3,S} {4,S} {17,S}
7  C u0 p0 c0 {1,S} {2,S} {8,D}
8  C u0 p0 c0 {5,S} {7,D} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2233)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {9,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {13,S}
7  C u1 p0 c0 {5,S} {9,S} {14,S}
8  C u1 p0 c0 {6,S} {10,S} {15,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1853)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {9,S} {21,S}
3  O u1 p2 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
8  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {2,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {20,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1306)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {16,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,D} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {9,D}
7  C u0 p0 c0 {5,D} {6,S} {14,S}
8  C u0 p0 c0 {2,S} {4,D} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2299)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {3,D} {5,S} {11,S}
8  C u0 p0 c0 {6,S} {10,D} {15,S}
9  C u0 p0 c0 {2,S} {10,S} {11,D}
10 C u0 p0 c0 {8,D} {9,S} {16,S}
11 C u0 p0 c0 {7,S} {9,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1172)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {16,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {7,S} {8,D}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u0 p0 c0 {3,S} {9,D} {12,S}
7  C u0 p0 c0 {4,S} {5,D} {14,S}
8  C u0 p0 c0 {4,D} {9,S} {13,S}
9  C u0 p0 c0 {6,D} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2308)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {3,S} {9,S} {10,S} {15,S}
9  C u1 p0 c0 {7,S} {8,S} {16,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2309)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {8,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {6,S} {8,S} {15,S}
10 C u0 p0 c0 {7,S} {11,D} {16,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2377)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {11,D} {15,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {2,S} {10,D} {11,S}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 C u0 p0 c0 {7,D} {9,S} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2375)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {8,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {10,D}
8  C u0 p0 c0 {3,D} {5,S} {6,S}
9  C u0 p0 c0 {6,S} {11,D} {15,S}
10 C u0 p0 c0 {7,D} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2251)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {9,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {4,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2376)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {9,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {11,D}
8  C u0 p0 c0 {5,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {15,S}
10 C u0 p0 c0 {3,D} {6,S} {11,S}
11 C u0 p0 c0 {7,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1448)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(2359)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {5,S} {11,S} {17,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2358)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {15,S}
10 C u1 p0 c0 {7,S} {11,S} {18,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2357)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {15,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {4,S} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2457)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {6,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {15,S}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u1 p0 c0 {6,S} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2454)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {5,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u1 p0 c0 {3,S} {6,S} {15,S}
8  C u0 p0 c0 {4,S} {9,D} {17,S}
9  C u0 p0 c0 {5,S} {8,D} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2455)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {6,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {15,S}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {17,S}
9  C u1 p0 c0 {7,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2522)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
9  C u0 p0 c0 {3,S} {5,S} {11,S} {18,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2493)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {9,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {11,S} {13,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u1 p0 c0 {4,S} {7,S} {8,S}
10 C u0 p0 c0 {1,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2524)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {15,S}
10 C u0 p0 c0 {7,S} {11,D} {19,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2523)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {17,S}
10 C u0 p0 c0 {6,S} {11,D} {19,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2355)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {13,S}
6  C u0 p0 c0 {3,S} {7,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {12,S} {15,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {16,S}
9  C u1 p0 c0 {1,S} {5,S} {8,S}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2578)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {7,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2577)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2575)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {7,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {13,S}
8  C u0 p0 c0 {1,S} {9,S} {11,S} {16,S}
9  C u1 p0 c0 {1,S} {5,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {17,S}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2630)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {6,S} {9,S}
3  O u0 p2 c0 {5,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {16,S}
8  C u0 p0 c0 {4,S} {5,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {6,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2613)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {9,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {12,S}
7  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2621)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {6,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {4,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {6,S} {10,D} {17,S}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2638)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {4,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2494)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {7,S} {18,S}
4  O u0 p2 c0 {8,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {13,S}
8  C u0 p0 c0 {4,S} {5,S} {11,S} {14,S}
9  C u1 p0 c0 {1,S} {7,S} {10,S}
10 C u0 p0 c0 {1,S} {9,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2479)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {3,S} {10,S} {11,S} {15,S}
10 C u1 p0 c0 {7,S} {9,S} {18,S}
11 C u0 p0 c0 {4,D} {8,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2353)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {7,S} {11,S} {17,S}
11 C u0 p0 c0 {4,S} {9,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2352)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
7  C u0 p0 c0 {2,S} {8,S} {10,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u1 p0 c0 {3,S} {5,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1857)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {7,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
7  C u0 p0 c0 {2,S} {8,S} {9,S} {14,S}
8  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {3,D} {5,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2567)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1401)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {9,S} {10,S} {12,S} {13,S}
6  C u1 p0 c0 {4,S} {7,S} {15,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {4,S} {10,D} {16,S}
9  C u0 p0 c0 {5,S} {7,D} {14,S}
10 C u0 p0 c0 {3,S} {5,S} {8,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2698)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2451)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {7,S} {18,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {13,S}
5  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u1 p0 c0 {4,S} {9,S} {16,S}
9  C u0 p0 c0 {7,D} {8,S} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2509)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {3,S} {18,S}
3  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {9,S} {12,S} {15,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u1 p0 c0 {4,S} {5,S} {16,S}
9  C u0 p0 c0 {6,S} {7,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2596)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {6,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2552)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {3,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {11,S} {13,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u1 p0 c0 {7,S} {8,S} {17,S}
11 C u0 p0 c0 {6,S} {9,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2829)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {15,S}
8  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {6,S} {8,S} {16,S} {17,S}
10 C u0 p0 c0 {6,S} {11,D} {19,S}
11 C u0 p0 c0 {7,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2812)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {18,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2828)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {11,S} {16,S}
9  C u0 p0 c0 {7,S} {10,S} {17,S} {18,S}
10 C u0 p0 c0 {3,S} {9,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2819)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {6,S} {18,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
9  C u0 p0 c0 {4,S} {10,D} {16,S}
10 C u0 p0 c0 {6,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(835)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {4,S} {19,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
5  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {18,S}
9  C u0 p0 c0 {7,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(977)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {10,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {14,S}
8  C u1 p0 c0 {4,S} {9,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2428)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {11,S} {12,S} {16,S}
9  C u1 p0 c0 {1,S} {6,S} {7,S}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2729)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {11,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
8  C u0 p0 c0 {9,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {2,S} {6,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {3,S} {8,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2525)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {9,S} {20,S}
3  O u0 p2 c0 {4,S} {10,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {18,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {2,S} {8,S} {11,S} {17,S}
10 C u0 p0 c0 {3,S} {6,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2476)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
8  C u0 p0 c0 {7,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {6,S} {11,S} {15,S} {18,S}
10 C u1 p0 c0 {3,S} {8,S} {11,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2859)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {10,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {7,S} {10,S} {16,S} {17,S}
9  C u1 p0 c0 {2,S} {6,S} {11,S}
10 C u0 p0 c0 {3,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2915)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {4,S} {10,S}
4  O u0 p2 c0 {3,S} {20,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {11,S} {14,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {17,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u1 p0 c0 {1,S} {7,S} {8,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2546)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {6,S} {17,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {15,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {4,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u0 p0 c0 {6,S} {9,D} {16,S}
9  C u0 p0 c0 {7,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2913)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {4,S} {10,S}
4  O u0 p2 c0 {3,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {17,S}
9  C u0 p0 c0 {1,S} {5,S} {10,D}
10 C u0 p0 c0 {3,S} {9,D} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2936)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {4,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {16,S}
7  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
8  C u1 p0 c0 {1,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {10,S} {17,S}
10 C u0 p0 c0 {3,D} {6,S} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2955)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {10,S}
3  O u0 p2 c0 {4,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {17,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {4,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {10,D}
10 C u0 p0 c0 {2,S} {5,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2962)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {16,S}
9  C u0 p0 c0 {2,S} {4,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2925)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u1 p0 c0 {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,D} {5,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3019)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {17,S}
3  O u0 p2 c0 {9,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u1 p0 c0 {2,S} {5,S} {10,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {3,S} {8,D} {10,S}
10 C u1 p0 c0 {7,S} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3017)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {8,S} {17,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {6,S} {15,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {2,S} {6,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2995)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {4,S} {17,S}
3  O u0 p2 c0 {8,S} {18,S}
4  C u0 p0 c0 {2,S} {6,S} {9,S} {11,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {14,S}
6  C u0 p0 c0 {4,S} {10,S} {12,S} {13,S}
7  C u1 p0 c0 {1,S} {5,S} {10,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {16,S}
10 C u1 p0 c0 {6,S} {7,S} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3014)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {10,D}
4  C u0 p0 c0 {1,S} {5,S} {9,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {10,S} {12,S}
6  C u0 p0 c0 {8,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {16,S}
8  C u1 p0 c0 {2,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {17,S}
10 C u0 p0 c0 {3,D} {5,S} {7,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2295)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {16,S}
2  O u0 p2 c0 {4,S} {15,S}
3  C u0 p0 c0 {1,S} {5,D} {7,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {3,D} {4,S} {14,S}
6  C u0 p0 c0 {4,D} {8,S} {11,S}
7  C u0 p0 c0 {3,S} {9,D} {13,S}
8  C u1 p0 c0 {6,S} {9,S} {10,S}
9  C u0 p0 c0 {7,D} {8,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(933)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {8,S} {10,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(961)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u1 p2 c0 {7,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {9,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {6,S} {7,D} {15,S}
10 C u0 p0 c0 {6,S} {8,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1313)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {16,S}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {3,S} {7,D} {12,S}
6  C u0 p0 c0 {4,D} {9,S} {15,S}
7  C u0 p0 c0 {5,D} {8,S} {13,S}
8  C u0 p0 c0 {7,S} {9,D} {14,S}
9  C u0 p0 c0 {2,S} {6,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(3112)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {17,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {18,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {13,S}
7  C u0 p0 c0 {5,S} {10,D} {12,S}
8  C u0 p0 c0 {3,S} {6,D} {9,S}
9  C u0 p0 c0 {8,S} {11,D} {16,S}
10 C u0 p0 c0 {7,D} {11,S} {14,S}
11 C u0 p0 c0 {9,D} {10,S} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3111)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {9,S} {18,S}
3  O u0 p2 c0 {7,S} {17,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {3,S} {6,D} {10,S}
8  C u0 p0 c0 {5,S} {11,D} {14,S}
9  C u0 p0 c0 {2,S} {10,D} {11,S}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 C u0 p0 c0 {8,D} {9,S} {15,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3113)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {8,S} {17,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,D}
7  C u0 p0 c0 {5,S} {11,D} {13,S}
8  C u0 p0 c0 {3,S} {9,S} {10,D}
9  C u0 p0 c0 {6,D} {8,S} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {15,S}
11 C u0 p0 c0 {7,D} {10,S} {14,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3114)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {7,S} {17,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,D}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {10,S} {14,S}
9  C u0 p0 c0 {6,D} {11,S} {16,S}
10 C u0 p0 c0 {8,S} {11,D} {13,S}
11 C u0 p0 c0 {9,S} {10,D} {15,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3117)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {11,S} {17,S}
10 C u0 p0 c0 {8,D} {11,S} {16,S}
11 C u0 p0 c0 {3,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3116)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {11,D}
8  C u0 p0 c0 {5,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {15,S}
10 C u0 p0 c0 {3,D} {5,S} {11,S}
11 C u0 p0 c0 {7,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3035)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {9,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {7,S} {13,S}
7  C u1 p0 c0 {6,S} {9,S} {15,S}
8  C u1 p0 c0 {5,S} {10,S} {14,S}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3044)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {9,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u1 p0 c0 {5,S} {6,S} {15,S}
9  C u0 p0 c0 {3,S} {7,D} {10,S}
10 C u1 p0 c0 {6,S} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2475)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {14,S}
9  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {1,S} {5,S} {6,S}
11 C u0 p0 c0 {4,D} {8,S} {9,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3082)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {16,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {11,S}
6  C u0 p0 c0 {2,S} {4,D} {8,S}
7  C u0 p0 c0 {5,D} {9,S} {13,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {7,S} {8,D} {14,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1311)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {16,S}
2  O u1 p2 c0 {9,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u0 p0 c0 {3,S} {8,D} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {15,S}
8  C u0 p0 c0 {6,D} {9,S} {14,S}
9  C u0 p0 c0 {2,S} {7,D} {8,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(930)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {4,S} {18,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {2,S} {7,S} {9,S} {11,S}
5  C u0 p0 c0 {1,S} {8,S} {10,S} {12,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {6,S} {15,S}
8  C u0 p0 c0 {3,D} {5,S} {6,S}
9  C u0 p0 c0 {4,S} {10,D} {16,S}
10 C u0 p0 c0 {5,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(3269)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {15,S}
9  C u0 p0 c0 {3,D} {5,S} {11,S}
10 C u0 p0 c0 {6,S} {11,D} {16,S}
11 C u0 p0 c0 {9,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2415)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {1,S} {18,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {4,D} {5,S} {9,S}
8  C u0 p0 c0 {6,D} {10,S} {15,S}
9  C u0 p0 c0 {7,S} {11,D} {13,S}
10 C u1 p0 c0 {8,S} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {14,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {11,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3267)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {8,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {14,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {7,D} {11,S} {15,S}
10 C u0 p0 c0 {8,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(3270)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u0 p0 c0 {6,S} {10,D} {15,S}
9  C u0 p0 c0 {7,D} {11,S} {17,S}
10 C u0 p0 c0 {8,D} {11,S} {16,S}
11 C u0 p0 c0 {3,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3312)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {16,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
6  C u0 p0 c0 {4,S} {9,D} {12,S}
7  C u0 p0 c0 {3,D} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {10,D} {13,S}
9  C u0 p0 c0 {6,D} {10,S} {15,S}
10 C u0 p0 c0 {8,D} {9,S} {14,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3237)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {3,S} {6,S} {10,S} {15,S}
9  C u1 p0 c0 {7,S} {11,S} {16,S}
10 C u0 p0 c0 {8,S} {11,D} {17,S}
11 C u0 p0 c0 {4,S} {9,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3238)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {11,S} {14,S}
8  C u0 p0 c0 {3,S} {9,S} {10,S} {15,S}
9  C u1 p0 c0 {6,S} {8,S} {16,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {7,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3239)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {8,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {3,S} {4,S} {9,S} {11,S}
9  C u1 p0 c0 {7,S} {8,S} {15,S}
10 C u0 p0 c0 {6,S} {11,D} {16,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2409)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {3,D} {5,S} {9,S}
8  C u0 p0 c0 {6,D} {10,S} {13,S}
9  C u0 p0 c0 {7,S} {11,D} {15,S}
10 C u1 p0 c0 {8,S} {11,S} {14,S}
11 C u0 p0 c0 {9,D} {10,S} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3061)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {16,S}
9  C u1 p0 c0 {3,S} {6,S} {11,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3384)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {6,S} {20,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3370)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {14,S}
9  C u0 p0 c0 {4,S} {6,S} {11,S} {17,S}
10 C u0 p0 c0 {3,S} {8,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3396)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {11,S} {13,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u1 p0 c0 {1,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {6,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3305)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {17,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {4,D} {5,S} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {9,D}
8  C u0 p0 c0 {6,S} {11,D} {13,S}
9  C u0 p0 c0 {7,D} {10,S} {14,S}
10 C u1 p0 c0 {9,S} {11,S} {16,S}
11 C u0 p0 c0 {8,D} {10,S} {15,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2288)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {15,S}
2  O u1 p2 c0 {7,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {7,S} {10,S}
5  C u0 p0 c0 {3,S} {8,D} {11,S}
6  C u1 p0 c0 {8,S} {9,S} {14,S}
7  C u0 p0 c0 {2,S} {4,S} {9,D}
8  C u0 p0 c0 {5,D} {6,S} {12,S}
9  C u0 p0 c0 {6,S} {7,D} {13,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(3311)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {16,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {5,S} {6,D} {14,S}
8  C u0 p0 c0 {3,D} {4,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {13,S}
10 C u0 p0 c0 {8,S} {9,D} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2315)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,S} {17,S}
3  O u0 p2 c0 {7,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u1 p0 c0 {5,S} {11,S} {13,S}
7  C u0 p0 c0 {3,D} {5,S} {9,S}
8  C u0 p0 c0 {2,S} {9,D} {10,S}
9  C u0 p0 c0 {7,S} {8,D} {14,S}
10 C u0 p0 c0 {8,S} {11,D} {15,S}
11 C u0 p0 c0 {6,S} {10,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3163)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {17,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u1 p0 c0 {5,S} {7,S} {13,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {14,S}
9  C u0 p0 c0 {7,D} {11,S} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {15,S}
11 C u0 p0 c0 {3,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2396)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,S} {17,S}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {11,D} {14,S}
7  C u0 p0 c0 {5,S} {9,D} {13,S}
8  C u0 p0 c0 {2,S} {9,S} {10,D}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {15,S}
11 C u0 p0 c0 {3,S} {6,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3147)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {17,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u1 p2 c0 {8,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {13,S}
7  C u0 p0 c0 {5,S} {9,D} {12,S}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  C u0 p0 c0 {7,D} {11,S} {14,S}
10 C u0 p0 c0 {8,S} {11,D} {16,S}
11 C u0 p0 c0 {9,S} {10,D} {15,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2434)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {11,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u1 p0 c0 {5,S} {9,S} {13,S}
8  C u0 p0 c0 {6,D} {11,S} {14,S}
9  C u0 p0 c0 {7,S} {10,D} {16,S}
10 C u0 p0 c0 {9,D} {11,S} {15,S}
11 C u0 p0 c0 {3,D} {8,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2220)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {16,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {3,S} {6,S} {10,S} {17,S}
10 C u1 p0 c0 {9,S} {11,S} {18,S}
11 C u0 p0 c0 {4,D} {7,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1056)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {6,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {14,S}
4  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
6  C u0 p0 c0 {1,S} {3,S} {9,D}
7  C u0 p0 c0 {2,D} {3,S} {4,S}
8  C u1 p0 c0 {5,S} {9,S} {15,S}
9  C u0 p0 c0 {6,D} {8,S} {16,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(1848)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {6,S} {16,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {9,D}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u0 p0 c0 {4,S} {6,D} {14,S}
8  C u1 p0 c0 {4,S} {9,S} {13,S}
9  C u0 p0 c0 {5,D} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2252)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {18,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {10,S} {19,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {16,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3552)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {6,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {11,S} {14,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {17,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u0 p0 c0 {3,D} {7,S} {8,S}
11 C u0 p0 c0 {6,S} {9,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(2248)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {1,S} {6,S} {11,S} {14,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {15,S}
9  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
10 C u1 p0 c0 {6,S} {8,S} {18,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3553)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
7  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {7,S} {11,S} {13,S} {14,S}
9  C u0 p0 c0 {3,D} {6,S} {7,S}
10 C u0 p0 c0 {5,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(2078)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {16,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {17,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u1 p0 c0 {4,S} {8,S} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {12,S}
7  C u0 p0 c0 {2,S} {6,D} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u1 p0 c0 {7,S} {10,S} {14,S}
10 C u0 p0 c0 {8,D} {9,S} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3240)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {8,S} {18,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {4,S} {9,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1300)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {9,S}
2  O u1 p2 c0 {6,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {2,S} {4,D} {8,S}
7  C u0 p0 c0 {5,D} {9,S} {15,S}
8  C u0 p0 c0 {6,S} {9,D} {14,S}
9  C u0 p0 c0 {1,S} {7,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(3554)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {9,S} {18,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {1,S} {6,S} {10,D}
9  C u0 p0 c0 {3,S} {6,S} {11,D}
10 C u0 p0 c0 {5,S} {8,D} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3630)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {7,S} {17,S}
11 C u0 p0 c0 {4,S} {8,S} {9,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3651)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {6,S} {17,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {12,S}
7  C u1 p0 c0 {5,S} {6,S} {13,S}
8  C u0 p0 c0 {6,S} {10,D} {14,S}
9  C u0 p0 c0 {5,S} {11,D} {15,S}
10 C u0 p0 c0 {4,S} {8,D} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {16,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3650)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {17,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {7,S} {20,S}
4  O u0 p2 c0 {10,S} {19,S}
5  C u0 p0 c0 {2,S} {7,S} {11,S} {13,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {12,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {14,S}
9  C u1 p0 c0 {6,S} {10,S} {15,S}
10 C u0 p0 c0 {4,S} {9,S} {11,D}
11 C u0 p0 c0 {5,S} {10,D} {16,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3609)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {16,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {7,S} {18,S}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {8,D}
6  C u0 p0 c0 {2,S} {4,S} {9,D}
7  C u1 p0 c0 {3,S} {8,S} {10,S}
8  C u0 p0 c0 {5,D} {7,S} {13,S}
9  C u0 p0 c0 {6,D} {10,S} {15,S}
10 C u1 p0 c0 {7,S} {9,S} {14,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2730)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
10 C u1 p0 c0 {3,S} {8,S} {11,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3531)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {6,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {8,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {1,S} {7,S} {11,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2727)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {10,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {10,S} {14,S} {15,S}
8  C u0 p0 c0 {9,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {1,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2728)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {10,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {11,S} {14,S} {17,S}
9  C u1 p0 c0 {1,S} {5,S} {6,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3730)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {14,S}
6  C u0 p0 c0 {2,S} {9,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {6,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {1,S} {5,S} {6,S}
11 C u0 p0 c0 {4,D} {8,S} {9,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2731)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,S} {19,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {7,S} {11,S} {17,S}
11 C u0 p0 c0 {4,S} {9,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2978)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {8,S} {21,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {10,S} {11,S} {13,S}
8  C u0 p0 c0 {4,S} {5,S} {11,D}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 C u0 p0 c0 {7,S} {8,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3366)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {8,S} {18,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {12,S} {15,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {16,S}
9  C u1 p0 c0 {1,S} {5,S} {8,S}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3555)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {18,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {10,S} {11,S} {13,S} {14,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {17,S}
10 C u0 p0 c0 {7,S} {8,D} {16,S}
11 C u0 p0 c0 {7,S} {9,D} {15,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2977)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {18,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {8,S} {21,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {13,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
7  C u0 p0 c0 {9,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {4,S} {5,S} {11,D}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {6,S} {9,D} {16,S}
11 C u0 p0 c0 {7,S} {8,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(934)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u1 p0 c0 {5,S} {10,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {2,S} {8,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(1301)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u1 p2 c0 {5,S}
2  O u1 p2 c0 {8,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {7,D} {12,S}
5  C u0 p0 c0 {1,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {8,S} {13,S}
7  C u0 p0 c0 {4,D} {9,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {14,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
        """),
)


species(
    label='S(2313)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {5,S} {10,S}
8  C u0 p0 c0 {6,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {11,S} {16,S}
10 C u0 p0 c0 {7,S} {11,D} {17,S}
11 C u0 p0 c0 {3,S} {9,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(1417)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u1 p2 c0 {7,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {16,S}
9  C u0 p0 c0 {3,D} {5,S} {11,S}
10 C u0 p0 c0 {6,S} {11,D} {15,S}
11 C u0 p0 c0 {9,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(1419)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {8,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {2,D} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u0 p0 c0 {6,S} {11,D} {15,S}
10 C u0 p0 c0 {8,D} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(2432)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u1 p2 c0 {8,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {9,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {9,D} {16,S}
8  C u0 p0 c0 {2,S} {5,S} {11,D}
9  C u0 p0 c0 {6,S} {7,D} {15,S}
10 C u0 p0 c0 {3,D} {6,S} {11,S}
11 C u0 p0 c0 {8,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
        """),
)


species(
    label='S(3805)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {18,S}
2  O u0 p2 c0 {6,S} {17,S}
3  O u0 p2 c0 {8,S} {19,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u0 p0 c0 {1,S} {4,S} {9,D}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  C u0 p0 c0 {7,D} {10,S} {15,S}
10 C u1 p0 c0 {8,S} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2394)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {11,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {10,D} {15,S}
8  C u0 p0 c0 {2,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {11,D} {17,S}
10 C u0 p0 c0 {7,D} {11,S} {16,S}
11 C u0 p0 c0 {3,S} {9,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(2425)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u1 p0 c0 {9,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(637)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {20,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {4,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {9,D} {19,S}
9  C u0 p0 c0 {7,S} {8,D} {18,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2713)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {11,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {6,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {9,S} {14,S} {17,S}
9  C u1 p0 c0 {1,S} {8,S} {11,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {1,S} {9,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3626)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
6  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {9,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {3,S} {6,S} {8,S}
10 C u0 p0 c0 {1,S} {5,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3901)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {19,S}
2  O u0 p2 c0 {6,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {15,S}
7  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
8  C u0 p0 c0 {7,S} {9,D} {18,S}
9  C u1 p0 c0 {6,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3900)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {19,S}
2  O u0 p2 c0 {7,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {2,S} {5,S} {9,D}
8  C u1 p0 c0 {6,S} {9,S} {17,S}
9  C u0 p0 c0 {7,D} {8,S} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3898)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {19,S}
2  O u0 p2 c0 {6,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {7,S} {8,S} {11,S}
7  C u1 p0 c0 {4,S} {6,S} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {18,S}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3899)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {19,S}
2  O u0 p2 c0 {6,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {6,S} {9,D} {18,S}
9  C u1 p0 c0 {7,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3937)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {21,S}
2  O u0 p2 c0 {6,S} {22,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
9  C u0 p0 c0 {5,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {6,S} {11,D} {20,S}
11 C u0 p0 c0 {9,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {10,S}
21 H u0 p0 c0 {1,S}
22 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(905)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {21,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {10,S} {22,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {16,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {11,S} {17,S}
9  C u0 p0 c0 {7,S} {10,S} {18,S} {19,S}
10 C u0 p0 c0 {3,S} {9,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {1,S}
22 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2429)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {5,S} {20,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {4,D} {7,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3628)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {10,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {11,S} {14,S} {17,S}
9  C u1 p0 c0 {3,S} {6,S} {7,S}
10 C u0 p0 c0 {1,S} {5,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2427)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
9  C u1 p0 c0 {5,S} {6,S} {16,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3992)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {5,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
9  C u0 p0 c0 {8,S} {11,S} {16,S} {17,S}
10 C u1 p0 c0 {5,S} {6,S} {18,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2569)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {7,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3993)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {5,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {15,S}
8  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {6,S} {11,S} {14,S} {18,S}
10 C u1 p0 c0 {1,S} {7,S} {8,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3989)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {6,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {10,S} {11,S} {16,S} {17,S}
10 C u1 p0 c0 {6,S} {9,S} {18,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3355)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {10,S} {14,S} {16,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {15,S}
9  C u1 p0 c0 {1,S} {6,S} {8,S}
10 C u0 p0 c0 {4,S} {7,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1438)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u1 p2 c0 {5,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {16,S}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {9,D} {17,S}
9  C u0 p0 c0 {4,S} {8,D} {18,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(3990)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {5,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {8,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {1,S} {6,S} {11,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3368)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {5,S} {10,S} {14,S}
9  C u0 p0 c0 {3,S} {7,S} {11,S} {15,S}
10 C u1 p0 c0 {8,S} {11,S} {18,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4040)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u1 p0 c0 {5,S} {7,S} {17,S}
11 C u0 p0 c0 {4,S} {8,S} {9,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4033)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {14,S}
6  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
7  C u1 p0 c0 {5,S} {6,S} {15,S}
8  C u0 p0 c0 {4,S} {9,D} {16,S}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4038)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
8  C u0 p0 c0 {10,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {5,S} {8,S} {17,S}
11 C u0 p0 c0 {4,S} {8,S} {9,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4035)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {3,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {3,S} {8,S} {11,S}
5  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {17,S}
9  C u1 p0 c0 {7,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2430)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {6,S} {10,D}
10 C u0 p0 c0 {4,S} {9,D} {11,S}
11 C u1 p0 c0 {8,S} {10,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4039)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {11,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
9  C u1 p0 c0 {1,S} {7,S} {11,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {1,S} {9,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4032)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {7,S}
2  O u0 p2 c0 {3,S} {18,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {10,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
7  C u0 p0 c0 {1,S} {3,S} {9,S} {16,S}
8  C u0 p0 c0 {4,S} {9,D} {17,S}
9  C u1 p0 c0 {7,S} {8,D}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4074)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {8,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {4,S} {6,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {3,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {1,S} {6,S} {11,S} {17,S}
9  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
10 C u0 p0 c0 {7,S} {11,D} {18,S}
11 C u0 p0 c0 {8,S} {10,D} {19,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {11,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3927)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {19,S}
2  O u0 p2 c0 {8,S} {20,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u1 p0 c0 {3,S} {6,S} {17,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {18,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4116)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {14,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {7,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(1807)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3693)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {17,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {9,S} {20,S}
5  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {9,D}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u0 p0 c0 {4,S} {7,D} {11,S}
10 C u0 p0 c0 {5,S} {8,D} {16,S}
11 C u1 p0 c0 {5,S} {9,S} {15,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4129)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {4,S} {18,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
6  C u0 p0 c0 {2,S} {4,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {12,S}
8  C u0 p0 c0 {4,S} {5,S} {14,S} {15,S}
9  C u0 p0 c0 {7,S} {10,D} {17,S}
10 C u0 p0 c0 {6,S} {9,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3353)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {3,S} {6,S} {10,S} {15,S}
9  C u0 p0 c0 {7,S} {11,S} {16,S} {17,S}
10 C u1 p0 c0 {8,S} {11,S} {18,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3991)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u1 p2 c0 {6,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {6,S} {10,S} {16,S} {17,S}
10 C u0 p0 c0 {9,S} {11,S} {18,S} {19,S}
11 C u0 p0 c0 {4,D} {7,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4117)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {9,S} {13,S}
9  C u1 p0 c0 {1,S} {8,S} {11,S}
10 C u0 p0 c0 {6,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4206)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u0 p2 c0 {7,S} {19,S}
4  O u0 p2 c0 {2,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {3,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {7,S} {11,D}
10 C u1 p0 c0 {6,S} {7,S} {17,S}
11 C u0 p0 c0 {8,S} {9,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4219)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {9,S}
3  O u0 p2 c0 {7,S} {18,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {13,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {11,S}
7  C u0 p0 c0 {3,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {2,S} {7,S} {10,D}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1568)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {7,S}
2  O u0 p2 c0 {6,S} {17,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {6,S} {8,S} {13,S} {14,S}
6  C u1 p0 c0 {2,S} {4,S} {5,S}
7  C u0 p0 c0 {1,S} {3,S} {9,D}
8  C u1 p0 c0 {5,S} {9,S} {15,S}
9  C u0 p0 c0 {7,D} {8,S} {16,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4165)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {4,S} {18,S}
4  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
7  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
8  C u0 p0 c0 {4,S} {10,S} {14,S} {15,S}
9  C u0 p0 c0 {5,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4041)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {11,S} {20,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
8  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {8,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {3,S} {9,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2426)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {6,S} {11,S} {16,S} {17,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(956)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {7,D}
4  C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
5  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {3,D} {4,S} {5,S}
8  C u0 p0 c0 {2,S} {4,S} {10,D}
9  C u1 p0 c0 {6,S} {10,S} {16,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3101)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {6,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u0 p0 c0 {2,S} {4,S} {10,D}
8  C u0 p0 c0 {5,S} {6,D} {14,S}
9  C u1 p0 c0 {5,S} {10,S} {15,S}
10 C u0 p0 c0 {7,D} {9,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3604)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {8,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {18,S}
4  O u1 p2 c0 {9,S}
5  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
6  C u0 p0 c0 {7,S} {9,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
8  C u0 p0 c0 {1,S} {5,S} {9,D}
9  C u0 p0 c0 {4,S} {6,S} {8,D}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {7,S} {10,D} {16,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3592)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {3,S} {5,S}
3  O u0 p2 c0 {2,S} {18,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {2,S} {7,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {11,S} {15,S}
7  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {1,S} {6,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {17,S}
10 C u0 p0 c0 {7,S} {11,D} {16,S}
11 C u0 p0 c0 {4,S} {6,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3350)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {2,S} {5,S} {11,S} {14,S}
8  C u0 p0 c0 {9,S} {10,S} {15,S} {16,S}
9  C u1 p0 c0 {3,S} {6,S} {8,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {7,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4303)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {8,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {14,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {8,S} {9,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u0 p0 c0 {3,D} {4,S} {7,S}
10 C u0 p0 c0 {5,S} {8,D} {16,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(4295)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {4,S} {5,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
5  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {3,D} {4,S} {6,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {7,S} {9,D} {15,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(2474)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {14,S}
9  C u0 p0 c0 {10,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {1,S} {6,S} {9,S}
11 C u0 p0 c0 {4,D} {8,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3063)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {9,S} {19,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {4,S} {7,S} {11,D}
10 C u1 p0 c0 {3,S} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3067)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {8,S} {9,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {16,S}
9  C u1 p0 c0 {1,S} {5,S} {8,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {7,S} {10,D} {17,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4338)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {11,S} {14,S}
8  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {8,S} {10,S} {17,S}
10 C u1 p0 c0 {9,S} {11,S} {18,S}
11 C u0 p0 c0 {4,D} {7,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4339)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {9,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {7,S} {10,S} {14,S}
9  C u0 p0 c0 {1,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {4,S} {9,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3530)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {13,S}
7  C u0 p0 c0 {2,S} {9,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {7,S} {11,S} {15,S} {18,S}
10 C u1 p0 c0 {3,S} {7,S} {8,S}
11 C u0 p0 c0 {4,D} {6,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3528)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {16,S}
7  C u0 p0 c0 {2,S} {8,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
9  C u0 p0 c0 {10,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {3,S} {7,S} {9,S}
11 C u0 p0 c0 {4,D} {6,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3731)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {11,S} {14,S}
9  C u0 p0 c0 {10,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {3,S} {6,S} {9,S}
11 C u0 p0 c0 {4,D} {8,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(144)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {7,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u1 p0 c0 {3,S} {9,S} {15,S}
7  C u0 p0 c0 {4,S} {5,D} {14,S}
8  C u0 p0 c0 {4,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(2909)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {4,S} {11,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {17,S}
6  C u0 p0 c0 {7,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {6,S} {10,S} {12,S}
9  C u0 p0 c0 {1,S} {5,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {18,S}
11 C u0 p0 c0 {3,S} {9,D} {10,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3533)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {6,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {10,S} {11,S} {17,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {11,S} {14,S} {18,S}
10 C u1 p0 c0 {1,S} {7,S} {8,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4378)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {5,S} {11,S} {13,S} {17,S}
9  C u1 p0 c0 {3,S} {6,S} {7,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='C7H8(122)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u1 p0 c0 {1,S} {5,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {6,S} {12,S}
5  C u0 p0 c0 {2,S} {7,D} {15,S}
6  C u1 p0 c0 {4,S} {7,S} {13,S}
7  C u0 p0 c0 {5,D} {6,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
        """),
)


species(
    label='C7H8O(213)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {2,S} {3,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
4  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {7,D} {16,S}
6  C u0 p0 c0 {3,S} {8,D} {15,S}
7  C u0 p0 c0 {3,S} {5,D} {14,S}
8  C u0 p0 c0 {4,S} {6,D} {13,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {5,S}
        """),
)


species(
    label='S(318)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {3,S} {7,S}
3  O u1 p2 c0 {2,S}
4  O u1 p2 c0 {1,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {10,D} {19,S}
9  C u0 p0 c0 {7,S} {11,D} {18,S}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {6,S} {9,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(298)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {3,S} {5,S}
2  O u0 p2 c0 {4,S} {6,S}
3  O u1 p2 c0 {1,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {13,S}
7  C u0 p0 c0 {10,S} {11,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {18,S}
9  C u0 p0 c0 {6,S} {11,D} {19,S}
10 C u0 p0 c0 {7,S} {8,D} {17,S}
11 C u0 p0 c0 {7,S} {9,D} {16,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {11,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(4358)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {6,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {3,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
9  C u0 p0 c0 {6,S} {11,S} {17,S} {18,S}
10 C u1 p0 c0 {1,S} {7,S} {11,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3728)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {2,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {14,S}
8  C u0 p0 c0 {6,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {3,S} {10,S} {11,S} {17,S}
10 C u1 p0 c0 {7,S} {9,S} {18,S}
11 C u0 p0 c0 {4,D} {8,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3407)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {5,S} {19,S}
4  O u0 p2 c0 {9,S} {20,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {9,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {11,S} {12,S} {16,S}
9  C u1 p0 c0 {4,S} {6,S} {7,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4269)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {15,S}
7  C u0 p0 c0 {2,S} {4,S} {8,D}
8  C u0 p0 c0 {6,S} {7,D} {17,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4377)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {6,S} {18,S}
3  O u0 p2 c0 {7,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {13,S}
6  C u0 p0 c0 {2,S} {5,S} {10,S} {12,S}
7  C u0 p0 c0 {3,S} {8,S} {10,S} {14,S}
8  C u0 p0 c0 {7,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {1,S} {5,S} {11,D}
10 C u1 p0 c0 {6,S} {7,S} {17,S}
11 C u0 p0 c0 {4,S} {8,S} {9,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3732)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {9,S} {20,S}
3  O u1 p2 c0 {8,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {14,S}
9  C u0 p0 c0 {2,S} {7,S} {11,S} {17,S}
10 C u0 p0 c0 {8,S} {11,S} {18,S} {19,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4460)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {20,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {10,S} {21,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {14,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {2,S} {9,S} {10,S} {15,S}
8  C u0 p0 c0 {6,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {4,D} {5,S} {7,S}
10 C u0 p0 c0 {3,S} {7,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4270)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {19,S}
2  O u0 p2 c0 {8,D}
3  O u1 p2 c0 {9,S}
4  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {4,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {2,D} {5,S} {6,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {18,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2973)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {11,S} {14,S} {15,S}
7  C u0 p0 c0 {2,S} {9,S} {10,S} {13,S}
8  C u0 p0 c0 {3,S} {5,S} {9,D}
9  C u0 p0 c0 {7,S} {8,D} {16,S}
10 C u0 p0 c0 {7,S} {11,D} {17,S}
11 C u0 p0 c0 {4,S} {6,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(458)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {7,D}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
6  C u1 p0 c0 {4,S} {9,S} {15,S}
7  C u0 p0 c0 {2,D} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {16,S}
9  C u0 p0 c0 {6,S} {8,D} {17,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2089)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {4,S} {9,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {13,S}
9  C u0 p0 c0 {3,S} {6,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4357)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {11,S} {14,S}
7  C u0 p0 c0 {2,S} {5,S} {10,S} {13,S}
8  C u0 p0 c0 {3,S} {9,S} {10,S} {15,S}
9  C u0 p0 c0 {8,S} {11,S} {16,S} {17,S}
10 C u1 p0 c0 {7,S} {8,S} {18,S}
11 C u0 p0 c0 {4,D} {6,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1871)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {7,D}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {15,S}
6  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u0 p0 c0 {3,D} {6,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {16,S}
10 C u0 p0 c0 {8,S} {9,D} {17,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(131)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {3,S}
2  O u1 p2 c0 {1,S}
3  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,S} {8,D} {14,S}
7  C u0 p0 c0 {5,D} {9,S} {16,S}
8  C u0 p0 c0 {6,D} {9,S} {17,S}
9  C u1 p0 c0 {7,S} {8,S} {15,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
        """),
)


species(
    label='S(4489)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {11,S} {20,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {11,D}
9  C u0 p0 c0 {4,D} {6,S} {7,S}
10 C u1 p0 c0 {5,S} {11,S} {17,S}
11 C u0 p0 c0 {3,S} {8,D} {10,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2121)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {9,S} {18,S}
4  O u0 p2 c0 {2,S} {19,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u1 p0 c0 {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {11,D}
10 C u1 p0 c0 {8,S} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2118)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {5,S} {18,S}
4  O u0 p2 c0 {2,S} {19,S}
5  C u0 p0 c0 {3,S} {7,S} {10,S} {12,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {13,S}
7  C u0 p0 c0 {5,S} {8,S} {14,S} {15,S}
8  C u1 p0 c0 {1,S} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {6,S} {11,D}
10 C u1 p0 c0 {5,S} {11,S} {16,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3470)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {2,S} {5,S}
2  O u0 p2 c0 {1,S} {17,S}
3  O u1 p2 c0 {7,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u1 p0 c0 {5,S} {9,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {10,S} {14,S}
9  C u0 p0 c0 {6,S} {11,D} {16,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {15,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {11,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4565)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,D}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {12,S}
6  C u0 p0 c0 {2,D} {4,S} {8,S}
7  C u0 p0 c0 {5,S} {9,D} {13,S}
8  C u0 p0 c0 {6,S} {10,D} {14,S}
9  C u0 p0 c0 {7,D} {10,S} {15,S}
10 C u0 p0 c0 {3,S} {8,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
        """),
)


species(
    label='S(4581)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {15,S}
3  O u1 p2 c0 {8,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,D}
6  C u0 p0 c0 {4,S} {9,D} {12,S}
7  C u0 p0 c0 {2,S} {9,S} {10,D}
8  C u0 p0 c0 {3,S} {5,D} {10,S}
9  C u0 p0 c0 {6,D} {7,S} {14,S}
10 C u0 p0 c0 {7,D} {8,S} {13,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4582)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {15,S}
3  O u1 p2 c0 {6,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,D}
6  C u0 p0 c0 {3,S} {4,S} {9,D}
7  C u0 p0 c0 {2,S} {9,S} {10,D}
8  C u0 p0 c0 {5,D} {10,S} {12,S}
9  C u0 p0 c0 {6,D} {7,S} {13,S}
10 C u0 p0 c0 {7,D} {8,S} {14,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4589)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {15,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {9,S} {12,S}
8  C u0 p0 c0 {6,D} {10,S} {13,S}
9  C u0 p0 c0 {7,S} {10,D} {14,S}
10 C u0 p0 c0 {3,S} {8,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4585)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {7,D}
3  O u1 p2 c0 {8,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {7,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {9,D}
7  C u0 p0 c0 {2,D} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {5,S} {10,D}
9  C u0 p0 c0 {6,D} {10,S} {14,S}
10 C u0 p0 c0 {8,D} {9,S} {15,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
        """),
)


species(
    label='S(4607)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {14,S}
3  O u0 p2 c0 {8,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,D}
6  C u0 p0 c0 {2,S} {4,S} {7,D}
7  C u0 p0 c0 {6,D} {8,S} {12,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,D} {10,S} {13,S}
10 C u1 p0 c0 {8,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4636)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {5,S} {16,S}
3  O u0 p2 c0 {7,S} {17,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {12,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u0 p0 c0 {5,S} {10,D} {13,S}
9  C u0 p0 c0 {7,D} {11,S} {14,S}
10 C u0 p0 c0 {8,D} {11,S} {15,S}
11 C u0 p0 c0 {4,D} {9,S} {10,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4587)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {15,S}
3  O u1 p2 c0 {10,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {6,D}
6  C u0 p0 c0 {2,S} {5,D} {8,S}
7  C u0 p0 c0 {4,S} {9,D} {12,S}
8  C u0 p0 c0 {6,S} {10,D} {13,S}
9  C u0 p0 c0 {7,D} {10,S} {14,S}
10 C u0 p0 c0 {3,S} {8,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4606)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {14,S}
3  O u0 p2 c0 {7,S} {15,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {8,D}
6  C u0 p0 c0 {1,S} {4,S} {10,D}
7  C u0 p0 c0 {3,S} {8,S} {9,D}
8  C u0 p0 c0 {5,D} {7,S} {12,S}
9  C u0 p0 c0 {7,D} {10,S} {13,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3882)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u0 p2 c0 {10,S} {21,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {4,S} {8,D} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2181)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {7,S} {9,S}
3  O u0 p2 c0 {8,S} {17,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {2,S} {8,S} {10,S} {15,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u0 p0 c0 {2,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4555)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {6,S}
2  O u0 p2 c0 {5,S} {8,S}
3  O u0 p2 c0 {5,S} {17,S}
4  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
7  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
8  C u1 p0 c0 {2,S} {6,S} {10,S}
9  C u0 p0 c0 {5,S} {10,D} {15,S}
10 C u0 p0 c0 {8,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4696)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {5,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {7,S} {17,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {1,S} {4,S} {8,S} {12,S}
6  C u0 p0 c0 {2,S} {4,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
8  C u1 p0 c0 {5,S} {7,S} {14,S}
9  C u0 p0 c0 {7,S} {10,D} {15,S}
10 C u0 p0 c0 {6,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4703)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {7,S}
2  O u0 p2 c0 {9,S} {10,S}
3  O u0 p2 c0 {5,S} {18,S}
4  O u0 p2 c0 {5,S} {19,S}
5  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
7  C u0 p0 c0 {1,S} {6,S} {10,S} {13,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
9  C u0 p0 c0 {2,S} {5,S} {11,S} {16,S}
10 C u0 p0 c0 {2,S} {7,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {17,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2182)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {8,S}
2  O u0 p2 c0 {5,S} {9,S}
3  O u0 p2 c0 {4,S} {17,S}
4  C u0 p0 c0 {3,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {2,S} {4,S} {10,S} {14,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {15,S}
7  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
8  C u1 p0 c0 {1,S} {6,S} {9,S}
9  C u0 p0 c0 {2,S} {8,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4701)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {9,S}
2  O u0 p2 c0 {5,S} {8,S}
3  O u0 p2 c0 {4,S} {17,S}
4  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,S} {11,S} {15,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {1,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {16,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4145)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {9,S} {11,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {2,S} {7,S} {8,D}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {4,S} {7,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4733)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {9,S}
2  O u0 p2 c0 {6,S} {8,S}
3  O u0 p2 c0 {5,S} {17,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {10,S} {12,S}
6  C u0 p0 c0 {2,S} {8,S} {10,S} {14,S}
7  C u0 p0 c0 {4,S} {9,S} {13,S} {15,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {1,S} {7,S} {8,D}
10 C u1 p0 c0 {5,S} {6,S} {16,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {10,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4698)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {4,S} {9,S}
3  O u0 p2 c0 {4,S} {17,S}
4  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
6  C u0 p0 c0 {1,S} {5,S} {9,S} {12,S}
7  C u0 p0 c0 {4,S} {5,S} {13,S} {14,S}
8  C u0 p0 c0 {4,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {2,S} {6,S} {10,D}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4740)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {7,S}
2  O u0 p2 c0 {6,S} {10,S}
3  O u0 p2 c0 {6,S} {18,S}
4  O u0 p2 c0 {5,S} {19,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
9  C u0 p0 c0 {6,S} {11,S} {15,S} {16,S}
10 C u0 p0 c0 {2,S} {7,S} {11,D}
11 C u0 p0 c0 {9,S} {10,D} {17,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4143)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {9,S} {19,S}
3  O u0 p2 c0 {8,S} {20,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {9,S} {11,S} {15,S} {16,S}
8  C u0 p0 c0 {3,S} {6,S} {11,D}
9  C u0 p0 c0 {2,S} {7,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {17,S}
11 C u0 p0 c0 {4,S} {7,S} {8,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4540)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {9,S} {20,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
8  C u0 p0 c0 {9,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {2,S} {8,S} {11,D}
10 C u0 p0 c0 {3,D} {6,S} {8,S}
11 C u0 p0 c0 {4,S} {7,S} {9,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4541)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {9,S} {20,S}
3  O u0 p2 c0 {10,D}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {10,S} {11,S} {17,S} {18,S}
9  C u0 p0 c0 {2,S} {7,S} {11,D}
10 C u0 p0 c0 {3,D} {6,S} {8,S}
11 C u0 p0 c0 {4,S} {8,S} {9,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(974)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {17,S}
2  O u0 p2 c0 {4,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {9,D} {15,S}
8  C u0 p0 c0 {3,S} {6,S} {10,D}
9  C u0 p0 c0 {7,D} {10,S} {16,S}
10 C u1 p0 c0 {8,D} {9,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4790)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {13,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {12,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {5,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {3,D} {5,S} {7,S}
10 C u0 p0 c0 {6,S} {11,D} {18,S}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4780)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {7,S} {19,S}
2  O u0 p2 c0 {6,S} {20,S}
3  O u0 p2 c0 {9,D}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {5,S} {9,S} {15,S}
7  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
8  C u0 p0 c0 {9,S} {11,S} {16,S} {17,S}
9  C u0 p0 c0 {3,D} {6,S} {8,S}
10 C u0 p0 c0 {7,S} {11,D} {18,S}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4375)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {9,S} {19,S}
4  O u0 p2 c0 {11,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {9,S} {11,S} {16,S} {17,S}
9  C u1 p0 c0 {3,S} {7,S} {8,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {4,S} {8,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4811)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {3,D} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {17,S}
10 C u1 p0 c0 {8,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4379)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {9,S} {11,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {7,S} {18,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {13,S}
8  C u0 p0 c0 {5,S} {10,S} {16,S} {17,S}
9  C u1 p0 c0 {1,S} {7,S} {11,S}
10 C u0 p0 c0 {4,S} {8,S} {11,D}
11 C u0 p0 c0 {1,S} {9,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4142)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {18,S}
2  O u0 p2 c0 {8,S} {19,S}
3  O u0 p2 c0 {9,S} {20,S}
4  O u1 p2 c0 {11,S}
5  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
6  C u0 p0 c0 {5,S} {11,S} {13,S} {14,S}
7  C u0 p0 c0 {8,S} {9,S} {15,S} {16,S}
8  C u0 p0 c0 {2,S} {7,S} {11,D}
9  C u0 p0 c0 {3,S} {7,S} {10,D}
10 C u0 p0 c0 {5,S} {9,D} {17,S}
11 C u0 p0 c0 {4,S} {6,S} {8,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4823)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {19,S}
2  O u1 p2 c0 {4,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {9,S} {15,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
8  C u0 p0 c0 {3,D} {7,S} {10,S}
9  C u0 p0 c0 {6,S} {10,D} {17,S}
10 C u0 p0 c0 {8,S} {9,D} {18,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(4822)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {9,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {7,S} {9,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {6,S} {15,S}
8  C u0 p0 c0 {5,S} {10,D} {16,S}
9  C u0 p0 c0 {3,D} {6,S} {10,S}
10 C u0 p0 c0 {8,D} {9,S} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4261)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 3
1  O u0 p2 c0 {5,S} {9,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {4,S} {7,S}
4  O u1 p2 c0 {3,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {5,S} {10,S} {15,S} {16,S}
9  C u1 p0 c0 {1,S} {7,S} {11,S}
10 C u0 p0 c0 {8,S} {11,D} {17,S}
11 C u0 p0 c0 {9,S} {10,D} {18,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(970)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u1 p0 c0 {4,S} {10,S} {15,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {16,S}
10 C u0 p0 c0 {7,S} {9,D} {17,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4820)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {4,S} {18,S}
2  O u0 p2 c0 {4,S} {19,S}
3  O u0 p2 c0 {8,D}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {10,S} {15,S} {16,S}
8  C u0 p0 c0 {3,D} {6,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {17,S}
10 C u1 p0 c0 {7,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(4858)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {6,S} {21,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {8,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,D} {17,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4894)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {8,S} {21,S}
4  O u0 p2 c0 {9,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {16,S}
9  C u0 p0 c0 {4,D} {7,S} {8,S}
10 C u0 p0 c0 {5,S} {11,D} {17,S}
11 C u0 p0 c0 {8,S} {10,D} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(4895)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {4,S} {8,S}
4  O u0 p2 c0 {3,S} {21,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u0 p0 c0 {5,S} {7,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {8,S} {14,S} {15,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u0 p0 c0 {5,S} {11,D} {16,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(4893)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {5,S} {19,S}
2  O u0 p2 c0 {5,S} {20,S}
3  O u0 p2 c0 {7,S} {21,S}
4  O u0 p2 c0 {10,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {5,S} {9,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {11,D} {17,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 C u0 p0 c0 {9,D} {10,S} {18,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {11,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(3627)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {11,S} {20,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {2,S} {7,S} {8,S} {12,S}
6  C u0 p0 c0 {1,S} {7,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
8  C u0 p0 c0 {4,S} {5,S} {9,S} {14,S}
9  C u0 p0 c0 {8,S} {11,S} {17,S} {18,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {3,S} {9,S} {10,D}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {9,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(1799)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {3,S} {20,S}
2  O u0 p2 c0 {3,S} {21,S}
3  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {6,S} {9,S} {16,S} {17,S}
8  C u0 p0 c0 {3,S} {9,D} {19,S}
9  C u0 p0 c0 {7,S} {8,D} {18,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(3897)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {20,S}
2  O u0 p2 c0 {5,S} {19,S}
3  C u0 p0 c0 {4,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {13,S}
7  C u1 p0 c0 {5,S} {9,S} {16,S}
8  C u0 p0 c0 {6,S} {9,D} {17,S}
9  C u0 p0 c0 {7,S} {8,D} {18,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(1798)',
    reactive=True,
    structure=adjacencyList(
        """
1  O u0 p2 c0 {6,S} {20,S}
2  O u0 p2 c0 {8,S} {21,S}
3  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {1,S} {4,S} {9,S} {16,S}
7  C u0 p0 c0 {5,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {2,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {19,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
        """),
)


species(
    label='S(2283)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {3,S} {4,S}
2  O u0 p2 c0 {8,S} {20,S}
3  O u0 p2 c0 {1,S} {21,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {8,S} {14,S} {15,S}
7  C u0 p0 c0 {9,S} {10,S} {16,S} {17,S}
8  C u0 p0 c0 {2,S} {6,S} {10,D}
9  C u1 p0 c0 {5,S} {7,S} {18,S}
10 C u0 p0 c0 {7,S} {8,D} {19,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(904)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {22,S}
2  O u0 p2 c0 {4,S} {8,S}
3  O u0 p2 c0 {9,S} {21,S}
4  O u1 p2 c0 {2,S}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {2,S} {5,S} {11,S} {18,S}
9  C u0 p0 c0 {3,S} {7,S} {10,S} {15,S}
10 C u0 p0 c0 {9,S} {11,D} {19,S}
11 C u0 p0 c0 {8,S} {10,D} {20,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {1,S}
        """),
)


species(
    label='S(2489)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {9,S} {20,S}
5  C u0 p0 c0 {2,S} {6,S} {7,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {16,S}
8  C u0 p0 c0 {3,S} {9,S} {11,S} {15,S}
9  C u1 p0 c0 {4,S} {7,S} {8,S}
10 C u0 p0 c0 {1,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(3729)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {6,S} {10,S}
2  O u0 p2 c0 {5,S} {19,S}
3  O u0 p2 c0 {7,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {2,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {3,S} {8,S} {11,S} {14,S}
8  C u0 p0 c0 {7,S} {10,S} {17,S} {18,S}
9  C u0 p0 c0 {5,S} {11,S} {15,S} {16,S}
10 C u1 p0 c0 {1,S} {6,S} {8,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)


species(
    label='S(2256)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {18,S}
3  O u0 p2 c0 {8,S} {19,S}
4  O u0 p2 c0 {10,S} {20,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {8,S} {9,S} {14,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {15,S}
9  C u1 p0 c0 {5,S} {7,S} {16,S}
10 C u0 p0 c0 {4,S} {6,S} {11,D}
11 C u0 p0 c0 {8,S} {10,D} {17,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {11,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
        """),
)


species(
    label='S(2477)',
    reactive=True,
    structure=adjacencyList(
        """
multiplicity 2
1  O u0 p2 c0 {5,S} {6,S}
2  O u0 p2 c0 {7,S} {19,S}
3  O u0 p2 c0 {10,S} {20,S}
4  O u0 p2 c0 {11,D}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {10,S} {13,S}
7  C u0 p0 c0 {2,S} {8,S} {11,S} {14,S}
8  C u0 p0 c0 {7,S} {10,S} {16,S} {17,S}
9  C u0 p0 c0 {5,S} {11,S} {15,S} {18,S}
10 C u1 p0 c0 {3,S} {6,S} {8,S}
11 C u0 p0 c0 {4,D} {7,S} {9,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
        """),
)

# Reaction systems
simpleReactor(
    temperature=[(600,'K'),(1000,'K')],
    pressure=[(10.0,'bar'),(40.0,'bar')],
    nSims=6,
    initialMoleFractions={
        #"C7H10": 1,
        "O2(2)":   0.108803, # phi=1 means 9.5 O2 per C7H10
        "N2(1)":   0.36, # 8.1 times as much N2 as O2
        "CHPD(1)":0.00989, # Cycloheptadiene C7H10    
"Ar(2)": 0,
"He(3)": 0,
"Ne(4)": 0,
"HO2(3)": 0,
"C7H9(5)": 0,
"C7H10(9)": 0,
"OO(8)": 0,
"H(12)": 0,
"C7H9O2(20)": 0,
"C7H9O2(19)": 0,
"S(36)": 0,
"S(35)": 0,
"C7H9O(53)": 0,
"OH(D)(44)": 0,
"C7H10O(67)": 0,
"C7H9O(69)": 0,
"C7H9O3(82)": 0,
"C7H9O3(83)": 0,
"C7H9O3(84)": 0,
"S(94)": 0,
"S(93)": 0,
"S(92)": 0,
"S(133)": 0,
"S(98)": 0,
"S(97)": 0,
"S(99)": 0,
"C7H9O(45)": 0,
"C7H8O(75)": 0,
"C7H8O3(85)": 0,
"S(101)": 0,
"S(179)": 0,
"S(180)": 0,
"S(108)": 0,
"S(183)": 0,
"S(126)": 0,
"S(134)": 0,
"S(221)": 0,
"S(222)": 0,
"S(234)": 0,
"S(226)": 0,
"S(225)": 0,
"S(255)": 0,
"S(245)": 0,
"S(273)": 0,
"S(184)": 0,
"S(191)": 0,
"S(188)": 0,
"S(151)": 0,
"S(279)": 0,
"S(312)": 0,
"S(310)": 0,
"S(311)": 0,
"S(309)": 0,
"S(330)": 0,
"S(314)": 0,
"S(315)": 0,
"S(313)": 0,
"S(301)": 0,
"S(362)": 0,
"S(363)": 0,
"S(335)": 0,
"S(324)": 0,
"S(317)": 0,
"S(366)": 0,
"S(323)": 0,
"S(331)": 0,
"S(365)": 0,
"S(355)": 0,
"S(316)": 0,
"S(339)": 0,
"S(404)": 0,
"S(401)": 0,
"S(274)": 0,
"S(388)": 0,
"S(427)": 0,
"S(426)": 0,
"S(436)": 0,
"C7H11(85)": 0,
"C7H11(86)": 0,
"S(111)": 0,
"C7H8(110)": 0,
"S(150)": 0,
"S(113)": 0,
"H2O(204)": 0,
"C7H9O(175)": 0,
"S(256)": 0,
"S(334)": 0,
"S(349)": 0,
"S(425)": 0,
"S(418)": 0,
"S(420)": 0,
"S(493)": 0,
"S(480)": 0,
"S(164)": 0,
"S(524)": 0,
"S(553)": 0,
"C7H8O(215)": 0,
"S(612)": 0,
"S(555)": 0,
"S(635)": 0,
"S(661)": 0,
"S(577)": 0,
"S(519)": 0,
"S(766)": 0,
"S(532)": 0,
"S(770)": 0,
"S(771)": 0,
"S(869)": 0,
"S(495)": 0,
"S(890)": 0,
"S(899)": 0,
"S(931)": 0,
"S(950)": 0,
"C7H9(93)": 0,
"S(985)": 0,
"S(151)": 0,
"S(140)": 0,
"S(230)": 0,
"C7H8O(212)": 0,
"S(231)": 0,
"C7H8O(222)": 0,
"S(191)": 0,
"S(177)": 0,
"S(189)": 0,
"S(235)": 0,
"S(419)": 0,
"S(239)": 0,
"S(178)": 0,
"S(238)": 0,
"S(683)": 0,
"S(237)": 0,
"S(686)": 0,
"S(518)": 0,
"S(401)": 0,
"S(420)": 0,
"C7H11(292)": 0,
"S(556)": 0,
"C7H9O(286)": 0,
"S(656)": 0,
"S(1074)": 0,
"S(661)": 0,
"S(1299)": 0,
"S(1287)": 0,
"S(1298)": 0,
"S(1347)": 0,
"S(1346)": 0,
"S(1382)": 0,
"S(1364)": 0,
"S(1493)": 0,
"S(1457)": 0,
"S(1492)": 0,
"S(1586)": 0,
"S(1584)": 0,
"S(1587)": 0,
"S(1577)": 0,
"S(1581)": 0,
"S(1630)": 0,
"S(1628)": 0,
"S(1678)": 0,
"S(1749)": 0,
"O(T)(511)": 0,
"S(1734)": 0,
"S(1185)": 0,
"S(500)": 0,
"S(1052)": 0,
"S(994)": 0,
"S(1057)": 0,
"S(1122)": 0,
"S(1854)": 0,
"S(1717)": 0,
"S(1108)": 0,
"S(1930)": 0,
"S(1143)": 0,
"S(1908)": 0,
"S(1852)": 0,
"S(1962)": 0,
"S(1953)": 0,
"S(1965)": 0,
"S(1997)": 0,
"S(1402)": 0,
"S(1579)": 0,
"S(2064)": 0,
"S(2072)": 0,
"S(2094)": 0,
"S(2095)": 0,
"S(2093)": 0,
"S(1624)": 0,
"S(2143)": 0,
"S(2124)": 0,
"S(2135)": 0,
"S(2117)": 0,
"S(2186)": 0,
"S(2170)": 0,
"S(2200)": 0,
"S(2192)": 0,
"S(1036)": 0,
"S(2233)": 0,
"S(1853)": 0,
"S(1306)": 0,
"S(2299)": 0,
"S(1172)": 0,
"S(2308)": 0,
"S(2309)": 0,
"S(2377)": 0,
"S(2375)": 0,
"S(2251)": 0,
"S(2376)": 0,
"S(1448)": 0,
"S(2359)": 0,
"S(2358)": 0,
"S(2357)": 0,
"S(2457)": 0,
"S(2454)": 0,
"S(2455)": 0,
"S(2522)": 0,
"S(2493)": 0,
"S(2524)": 0,
"S(2523)": 0,
"S(2355)": 0,
"S(2578)": 0,
"S(2577)": 0,
"S(2575)": 0,
"S(2630)": 0,
"S(2613)": 0,
"S(2621)": 0,
"S(2638)": 0,
"S(2494)": 0,
"S(2479)": 0,
"S(2353)": 0,
"S(2352)": 0,
"S(1857)": 0,
"S(2567)": 0,
"S(1401)": 0,
"S(2698)": 0,
"S(2451)": 0,
"S(2509)": 0,
"S(2596)": 0,
"S(2552)": 0,
"S(2829)": 0,
"S(2812)": 0,
"S(2828)": 0,
"S(2819)": 0,
"S(835)": 0,
"S(977)": 0,
"S(2428)": 0,
"S(2729)": 0,
"S(2525)": 0,
"S(2476)": 0,
"S(2859)": 0,
"S(2915)": 0,
"S(2546)": 0,
"S(2913)": 0,
"S(2936)": 0,
"S(2955)": 0,
"S(2962)": 0,
"S(2925)": 0,
"S(3019)": 0,
"S(3017)": 0,
"S(2995)": 0,
"S(3014)": 0,
"S(2295)": 0,
"S(933)": 0,
"S(961)": 0,
"S(1313)": 0,
"S(3112)": 0,
"S(3111)": 0,
"S(3113)": 0,
"S(3114)": 0,
"S(3117)": 0,
"S(3116)": 0,
"S(3035)": 0,
"S(3044)": 0,
"S(2475)": 0,
"S(3082)": 0,
"S(1311)": 0,
"S(930)": 0,
"S(3269)": 0,
"S(2415)": 0,
"S(3267)": 0,
"S(3270)": 0,
"S(3312)": 0,
"S(3237)": 0,
"S(3238)": 0,
"S(3239)": 0,
"S(2409)": 0,
"S(3061)": 0,
"S(3384)": 0,
"S(3370)": 0,
"S(3396)": 0,
"S(3305)": 0,
"S(2288)": 0,
"S(3311)": 0,
"S(2315)": 0,
"S(3163)": 0,
"S(2396)": 0,
"S(3147)": 0,
"S(2434)": 0,
"S(2220)": 0,
"S(1056)": 0,
"S(1848)": 0,
"S(2252)": 0,
"S(3552)": 0,
"S(2248)": 0,
"S(3553)": 0,
"S(2078)": 0,
"S(3240)": 0,
"S(1300)": 0,
"S(3554)": 0,
"S(3630)": 0,
"S(3651)": 0,
"S(3650)": 0,
"S(3609)": 0,
"S(2730)": 0,
"S(3531)": 0,
"S(2727)": 0,
"S(2728)": 0,
"S(3730)": 0,
"S(2731)": 0,
"S(2978)": 0,
"S(3366)": 0,
"S(3555)": 0,
"S(2977)": 0,
"S(934)": 0,
"S(1301)": 0,
"S(2313)": 0,
"S(1417)": 0,
"S(1419)": 0,
"S(2432)": 0,
"S(3805)": 0,
"S(2394)": 0,
"S(2425)": 0,
"S(637)": 0,
"S(2713)": 0,
"S(3626)": 0,
"S(3901)": 0,
"S(3900)": 0,
"S(3898)": 0,
"S(3899)": 0,
"S(3937)": 0,
"S(905)": 0,
"S(2429)": 0,
"S(3628)": 0,
"S(2427)": 0,
"S(3992)": 0,
"S(2569)": 0,
"S(3993)": 0,
"S(3989)": 0,
"S(3355)": 0,
"S(1438)": 0,
"S(3990)": 0,
"S(3368)": 0,
"S(4040)": 0,
"S(4033)": 0,
"S(4038)": 0,
"S(4035)": 0,
"S(2430)": 0,
"S(4039)": 0,
"S(4032)": 0,
"S(4074)": 0,
"S(3927)": 0,
"S(4116)": 0,
"S(1807)": 0,
"S(3693)": 0,
"S(4129)": 0,
"S(3353)": 0,
"S(3991)": 0,
"S(4117)": 0,
"S(4206)": 0,
"S(4219)": 0,
"S(1568)": 0,
"S(4165)": 0,
"S(4041)": 0,
"S(2426)": 0,
"S(956)": 0,
"S(3101)": 0,
"S(3604)": 0,
"S(3592)": 0,
"S(3350)": 0,
"S(4303)": 0,
"S(4295)": 0,
"S(2474)": 0,
"S(3063)": 0,
"S(3067)": 0,
"S(4338)": 0,
"S(4339)": 0,
"S(3530)": 0,
"S(3528)": 0,
"S(3731)": 0,
"S(144)": 0,
"S(2909)": 0,
"S(3533)": 0,
"S(4378)": 0,
"C7H8(122)": 0,
"C7H8O(213)": 0,
"S(318)": 0,
"S(298)": 0,
"S(4358)": 0,
"S(3728)": 0,
"S(3407)": 0,
"S(4269)": 0,
"S(4377)": 0,
"S(3732)": 0,
"S(4460)": 0,
"S(4270)": 0,
"S(2973)": 0,
"S(458)": 0,
"S(2089)": 0,
"S(4357)": 0,
"S(1871)": 0,
"S(131)": 0,
"S(4489)": 0,
"S(2121)": 0,
"S(2118)": 0,
"S(3470)": 0,
"S(4565)": 0,
"S(4581)": 0,
"S(4582)": 0,
"S(4589)": 0,
"S(4585)": 0,
"S(4607)": 0,
"S(4636)": 0,
"S(4587)": 0,
"S(4606)": 0,
"S(3882)": 0,
"S(2181)": 0,
"S(4555)": 0,
"S(4696)": 0,
"S(4703)": 0,
"S(2182)": 0,
"S(4701)": 0,
"S(4145)": 0,
"S(4733)": 0,
"S(4698)": 0,
"S(4740)": 0,
"S(4143)": 0,
"S(4540)": 0,
"S(4541)": 0,
"S(974)": 0,
"S(4790)": 0,
"S(4780)": 0,
"S(4375)": 0,
"S(4811)": 0,
"S(4379)": 0,
"S(4142)": 0,
"S(4823)": 0,
"S(4822)": 0,
"S(4261)": 0,
"S(970)": 0,
"S(4820)": 0,
"S(4858)": 0,
"S(4894)": 0,
"S(4895)": 0,
"S(4893)": 0,
"S(3627)": 0,
"S(1799)": 0,
"S(3897)": 0,
"S(1798)": 0,
"S(2283)": 0,
"S(904)": 0,
"S(2489)": 0,
"S(3729)": 0,
"S(2256)": 0,
"S(2477)": 0,

        
    },
    terminationTime = (5.0, 's'),
    terminationRateRatio = 0.01,
    terminationConversion={
                'CHPD(1)': 0.99,
    },
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)


model(
    toleranceKeepInEdge=0, # No pruning to start
    toleranceMoveToCore=0.4,
    toleranceInterruptSimulation=1,
    maxNumObjsPerIter=3,      #
    terminateAtMaxObjects=True,
    maxNumSpecies=100, # first stage is until core reaches 100 species
    filterReactions=True, # should speed up model generation
    filterThreshold=2e8,
)

model(
    toleranceMoveToCore=0.8,
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
    generateOutputHTML=True,
    generatePlots=True,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
    generateSeedEachIteration=False,
)



