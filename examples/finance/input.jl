using JuMP
#sets
Investments = 1:2 

InitialWealth = 55000 
TargetCapital = 80000 
ExcessPerUnitBenefit = 1 
ShortagePerUnitCost = 4 

#the indices of return are Investments(STOCKs, bonds) * (high, low)
Return = [1.25 1.06;
			1.14  1.12 ]






















