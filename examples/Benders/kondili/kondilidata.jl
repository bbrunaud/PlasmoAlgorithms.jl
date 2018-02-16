H = 9
INF1 = 1000

tasks = 1:5
units = 1:4
periods = 1:10


#TABLE K(I,J) unit j that can perform task i
K = [
      1  0  0  0
      0  1  1  0
      0  1  1  0
      0  1  1  0
      0  0  0  1
    ]

UnitsForTask = Dict(
                    1 => [1],
                    2 => [2,3],
                    3 => [2,3],
                    4 => [2,3],
                    5 => [4]
)

TasksForUnit = Dict(
    1 => [1],
    2 => [2,3,4],
    3 => [2,3,4],
    4 => [5]
  )

#PARAMETER P(I) processing times /
ProcessingTime = [
            1
            2
            2
            1
            2
]


#PARAMETER V(J) capacity of process units /
UnitCapacity= [
            100
             80
             50
            200
]


## STN Data
states = 1:9

# TABLE TSET(I,S) task i receiving material in state s
TasksIn = Dict(
      1 => [1],
      2 => [2],
      3 => [],
      4 => [3],
      5 => [2],
      6 => [4],
      7 => [3,4],
      8 => [5],
      9 => []
)



#TABLE TSETB(I,S) task i producing material in state s
TasksOut = Dict(
      1 => [],
      2 => [1],
      3 => [2],
      4 => [],
      5 => [3],
      6 => [2,5],
      7 => [],
      8 => [4],
      9 => [5]
)

#TABLE RHO(I,S) proportion of input to task i of state s
ρin = [
     1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.4  0.0  0.0  0.6  0.0  0.0  0.0  0.0
     0.0  0.0  0.0  0.5  0.0  0.0  0.5  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  0.8  0.2  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0
]

#TABLE RHOB(I,S) proportion of output to task i of state s
ρout = [
     0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  0.4  0.0  0.0  0.6  0.0  0.0  0.0
     0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0
     0.0  0.0  0.0  0.0  0.0  0.1  0.0  0.0  0.9
     ]

#PARAMETER START(S)  initial amount  /
InitialInventory = [
            2000
             0
             0
            2000
             0
             0
            2000
             0
             0
]


#PARAMETER CD(S) sale price /
Price = [
              0
              -1
              10
              0
              -1
              -1
              0
              -1
              10
]

#PARAMETER CAP(S)  available storage capacity  /
StorageCapacity = [
            2000
            100
            2000
            2000
            150
            200
            2000
            100
            2000
]

ϵ =  0.01
