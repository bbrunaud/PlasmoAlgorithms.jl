abstract type CutData end

struct BendersCutData <: CutData
  θk
  λk
  xk
end

struct LLIntegerCutData <: CutData
  θlb
  yk
end

struct IntegerCutData <: CutData
  yk
end

struct LagrangeCutData <: CutData
  zk
  λc
  coeffs
  xk
end

struct LagrangeCrossCutData <: CutData
    zk
    λk
end

(==)(cd1::BendersCutData,cd2::BendersCutData) = (cd1.θk == cd2.θk) && (cd1.λk == cd2.λk) && (cd1.xk == cd2.xk)
(==)(cd1::LLIntegerCutData,cd2::LLIntegerCutData) = (cd1.θlb == cd2.θlb) &&  (cd1.yk == cd2.yk)
(==)(cd1::IntegerCutData,cd2::IntegerCutData) = (cd1.yk == cd2.yk)
(==)(cd1::LagrangeCutData,cd2::LagrangeCutData) = (cd1.zk == cd2.zk) && (cd1.λc == cd2.λc) && (cd1.coeffs == cd2.coeffs) && (cd1.xk == cd2.xk)
(==)(cd1::LagrangeCrossCutData,cd2::LagrangeCrossCutData) = (cd1.zk == cd2.zk) && (cd1.λk == cd2.λk)
