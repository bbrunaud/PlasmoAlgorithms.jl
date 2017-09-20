using MySQL
using JuMP
using JLD

numPlants = 14
numPWarehouses = 8
numSWarehouses = 10
numCustomers = 105
numProducts = 78
numPeriods = 60

# Customer, Products and Periods ranges
rCu = 1:numCustomers
rPr = 1:numProducts
pers = 1:numPeriods

firstP = rPr[1]
lastP = rPr[end]

t1 = pers[1]
tN = pers[end]

Ccode = 3000
firstCcode = rCu[1]+Ccode
lastCcode = rCu[end]+Ccode

function populateParam(ct, df)
	ar = convert(Array,df)
	for j in 1:size(ar,1)
		ct[ar[j,1:(end-1)]...] = ar[j,end]
	end
end

con = mysql_connect("localhost","bbrunaud","Chile1.","DFLdata")

plants = ["MFG$i" for i in 1:numPlants]
pwarehouses = ["WHP$i" for i in 1:numPWarehouses]
swarehouses = ["WHS$i" for i in 1:numSWarehouses]
warehouses = vcat(pwarehouses,swarehouses)
customers = ["CUS$i" for i in 1:numCustomers]
sites = vcat(plants,warehouses,customers)
products = ["PROD$i" for i in 1:numProducts]
periods = 1:numPeriods
levels = ["low", "high"]

#s0def = [(j,p) => 0 for j in warehouses for p in products]

modes = [0.001,
		3.5,
		4.5,
		8,
		10,
		15,
		17,
		20,
		21,
		24,
		25,
		27,
		30,
		32,
		33,
		35,
]

Capacity = Dict(
"MFG1"	=>	2827.875,
"MFG2"	=>	100,
"MFG3"	=>	2137.244,
"MFG4"	=>	3649.901,
"MFG5"	=>	5674.399,
"MFG6"	=>	2177.908,
"MFG7"	=>	100,
"MFG8"	=>	100,
"MFG9"	=>	100,
"MFG10"	=>	100,
"MFG11"	=>	3234.959,
"MFG12"	=>	5519.028,
"MFG13"	=>	307,
"MFG14"	=>	7444.689
)

HoldingCost = Dict(j=>0.02 for j in warehouses)
SizeCost = Dict(j=>0.1 for j in warehouses)

FixedCost = Dict(
"WHP1"	=>	5.906,
"WHP2"	=>	0.001,
"WHP3"	=>	0.241,
"WHP4"	=>	6.937,
"WHP5"	=>	0.351,
"WHP6"	=>	9.456,
"WHP7"	=>	10.956,
"WHP8"	=>	11.868,
"WHS1"	=>	3.76,
"WHS2"	=>	0.001,
"WHS3"	=>	8.567,
"WHS4"	=>	5.854,
"WHS5"	=>	0.001,
"WHS6"	=>	0.001,
"WHS7"	=>	10.781,
"WHS8"	=>	0.001,
"WHS9"	=>	4.2,
"WHS10"	=>	0.001
)

Threshold = Dict(
"WHP1"	=>	0	,
"WHP2"	=>	0	,
"WHP3"	=>	35	,
"WHP4"	=>	490	,
"WHP5"	=>	10.5,
"WHP6"	=>	840	,
"WHP7"	=>	840	,
"WHP8"	=>	840	,
"WHS1"	=>	0	,
"WHS2"	=>	0	,
"WHS3"	=>	0	,
"WHS4"	=>	105	,
"WHS5"	=>	175	,
"WHS6"	=>	0	,
"WHS7"	=>	0	,
"WHS8"	=>	0	,
"WHS9"	=>	140	,
"WHS10"	=>	0
)

FlowCost = Dict(
(	"WHP1"	,	"low"	)	=>	0,
(	"WHP2"	,	"low"	)	=>	0.007785714,
(	"WHP3"	,	"low"	)	=>	0,
(	"WHP4"	,	"low"	)	=>	0,
(	"WHP5"	,	"low"	)	=>	0,
(	"WHP6"	,	"low"	)	=>	0,
(	"WHP7"	,	"low"	)	=>	0.0101214286,
(	"WHP8"	,	"low"	)	=>	0.01695,
(	"WHS1"	,	"low"	)	=>	0,
(	"WHS2"	,	"low"	)	=>	0.0075,
(	"WHS3"	,	"low"	)	=>	0,
(	"WHS4"	,	"low"	)	=>	0,
(	"WHS5"	,	"low"	)	=>	0,
(	"WHS6"	,	"low"	)	=>	0.0107142857,
(	"WHS7"	,	"low"	)	=>	0,
(	"WHS8"	,	"low"	)	=>	0,
(	"WHS9"	,	"low"	)	=>	0,
(	"WHS10"	,	"low"	)	=>	0.008085714,
(	"WHP1"	,	"high"	)	=>	0	,
(	"WHP2"	,	"high"	)	=>	0.007785714	,
(	"WHP3"	,	"high"	)	=>	0.003442857	,
(	"WHP4"	,	"high"	)	=>	0.007078571	,
(	"WHP5"	,	"high"	)	=>	0.016657143	,
(	"WHP6"	,	"high"	)	=>	0.016457143	,
(	"WHP7"	,	"high"	)	=>	0.016642857	,
(	"WHP8"	,	"high"	)	=>	0.024014286	,
(	"WHS1"	,	"high"	)	=>	0	,
(	"WHS2"	,	"high"	)	=>	0.0075	,
(	"WHS3"	,	"high"	)	=>	0	,
(	"WHS4"	,	"high"	)	=>	0.037414964	,
(	"WHS5"	,	"high"	)	=>	0.009271429	,
(	"WHS6"	,	"high"	)	=>	0.0107142857	,
(	"WHS7"	,	"high"	)	=>	0	,
(	"WHS8"	,	"high"	)	=>	0	,
(	"WHS9"	,	"high"	)	=>	0.038014286	,
(	"WHS10"	,	"high"	)	=>	0.008085714
)

MaxYearlyDemand = Dict(
"PROD34"	=> 	126.001	,
"PROD16"	=> 	473.651	,
"PROD33"	=> 	121.563	,
"PROD31"	=> 	147.005	,
"PROD73"	=> 	14.997	,
"PROD50"	=> 	70	,
"PROD40"	=> 	119.212	,
"PROD26"	=> 	215.003	,
"PROD37"	=> 	159.994	,
"PROD17"	=> 	378.004	,
"PROD2"	=> 	7801.875	,
"PROD62"	=> 	35.001	,
"PROD60"	=> 	84.587	,
"PROD78"	=> 	0.03	,
"PROD15"	=> 	654.504	,
"PROD24"	=> 	264.664	,
"PROD22"	=> 	267.011	,
"PROD12"	=> 	587.499	,
"PROD35"	=> 	104.501	,
"PROD32"	=> 	118.192	,
"PROD76"	=> 	0.301	,
"PROD39"	=> 	109.9	,
"PROD64"	=> 	16.006	,
"PROD77"	=> 	0.198	,
"PROD18"	=> 	340.973	,
"PROD1"	=> 	15034.315	,
"PROD7"	=> 	1556.728	,
"PROD23"	=> 	258.507	,
"PROD59"	=> 	21.2	,
"PROD8"	=> 	1249.607	,
"PROD14"	=> 	537.068	,
"PROD52"	=> 	45.005	,
"PROD44"	=> 	78.738	,
"PROD58"	=> 	21.002	,
"PROD5"	=> 	3111.659	,
"PROD25"	=> 	187.397	,
"PROD36"	=> 	112.336	,
"PROD53"	=> 	38.003	,
"PROD66"	=> 	12.006	,
"PROD46"	=> 	64.395	,
"PROD11"	=> 	609.914	,
"PROD19"	=> 	305.004	,
"PROD4"	=> 	3454.559	,
"PROD47"	=> 	59.997	,
"PROD68"	=> 	9.001	,
"PROD30"	=> 	320.801	,
"PROD57"	=> 	27.898	,
"PROD43"	=> 	74.436	,
"PROD51"	=> 	38.999	,
"PROD38"	=> 	101.907	,
"PROD28"	=> 	181.787	,
"PROD45"	=> 	315.316	,
"PROD6"	=> 	2632.352	,
"PROD20"	=> 	357.61	,
"PROD61"	=> 	19.399	,
"PROD48"	=> 	58.504	,
"PROD54"	=> 	37.8	,
"PROD56"	=> 	84.997	,
"PROD63"	=> 	20.001	,
"PROD74"	=> 	2	,
"PROD21"	=> 	319.387	,
"PROD10"	=> 	737.667	,
"PROD75"	=> 	0.298	,
"PROD41"	=> 	107.906	,
"PROD70"	=> 	13.194	,
"PROD65"	=> 	11	,
"PROD9"	=> 	1407.993	,
"PROD13"	=> 	561.005	,
"PROD3"	=> 	4294.824	,
"PROD29"	=> 	130	,
"PROD71"	=> 	4.403	,
"PROD27"	=> 	150.009	,
"PROD67"	=> 	39.998	,
"PROD42"	=> 	67.594	,
"PROD55"	=> 	38.501	,
"PROD49"	=> 	98.509	,
"PROD69"	=> 	7.479	,
"PROD72"	=> 	18.694
)

# Get Demands
query = """
		SELECT Customer, Product, Period, Demand
		FROM Demands
		WHERE
			(ProductNumber BETWEEN $firstP AND $lastP)
			AND
			(Period BETWEEN $t1 AND $tN)
			AND
			(SiteCode BETWEEN $firstCcode AND $lastCcode)
		"""

demands = mysql_execute(con,query)

Demand = Dict((k,p,t) => 0.0 for k in customers for p in products for t in periods)
populateParam(Demand,demands)

# Get Routes
query = """
		SELECT Origin, Destination, Product, Active
		FROM Routes
		WHERE
			(ProductNumber BETWEEN $firstP AND $lastP)
			AND
			((DestinationType='Customer' AND (DestinationCode BETWEEN $firstCcode AND $lastCcode) AND OriginType<>'Plant')
			 OR
			 DestinationType='Warehouse')
		"""
Route = Dict((i,j,p) => 0 for i in sites for j in sites for p in products)

populateParam(Route,mysql_execute(con,query))

# Get transportation costs
query = """
		SELECT Origin, Destination, Mode, CostPerUnit
		FROM TransportationCosts
		WHERE
			Origin IN $(tuple(sites...))
			AND
			Destination IN $(tuple(sites...))
		"""

transpcosts = mysql_execute(con,query)
modes = unique(transpcosts[:Mode])
TransportationCost = Dict((i,j,f) => 0.0 for i in sites for j in sites for f in modes)
populateParam(TransportationCost,transpcosts)

origins = vcat(plants,warehouses)
destinations = vcat(warehouses,customers)

TCSingleMode = Dict((i,j) => 0.0 for i in origins for j in destinations)
for i in origins
  for j in destinations
    nzcosts = [TransportationCost[i,j,f]/f for f in modes if TransportationCost[i,j,f]>0]
    nzmodes = [f for f in modes if TransportationCost[i,j,f]>0]
    if length(nzmodes) == 0
      TCSingleMode[i,j] = 0.0
    else
      val = nzcosts'*nzmodes/sum(nzmodes)
      TCSingleMode[i,j] = val[1]
      end
  end
end
