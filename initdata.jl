#===variable descriptions===#

# parameters and counters
# alpha       # capital-share parameter in production function
# delta       # depreciation rate
# sigL        # elasticity of hours-disutility
# fag         # first economically active age group
# rho         # discount rate households
# sigma       # elasticity of intertemporal substitution
# tt          # time index

# age-dependent variables
# parlv0      # multiplicative shift parameter hours-disutillity
# parlv1      # additive shift parameter hours-disutillity

# time-dependent variables
# A          # aggregate assets
# Cons       # aggregate consumption
# CG         # aggregate public consumption
# DG         # public debt
# Div        # dividends
# Exp        # government expenditure
# Inv        # investment
# K          # capital stock
# LD         # aggregate labor demand
# LS         # aggregate labor supply
# N          # aggregate population
# NB         # number of newborns
# Nc         # number of children
# Nr         # number of retirees
# Nw         # number of workers
# P          # pension expenditure
# PB         # primary balance
# Rev        # government revenue
# TaxF       # payroll tax revenue
# TFP        # TFP stock
# V          # value of representative firm
# Y          # output (GDP)
# ab         # aggregate accidental bequests
# eda        # excess demand asset market
# edab       # (excess demand) accidental bequests
# edg        # (excess demand) government budget
# ediv       # (excess demand) intervivo transfers
# edl        # excess demand labor market
# edw        # (excess demand) Walras' Law
# edy        # excess demand goods market
# iv         # aggregate received intervivo transfers
# pc         # price of consumption
# qTob       # Tobin's q
# r          # real interest rate
# rag        # retirement age group
# tauC       # consumption tax rate
# tauF       # payroll tax rate
# tauW       # average wage tax rate
# taul       # lump-sum tax rate
# tauprof    # corporate tax rate
# uck        # user cost of capital
# w          # wage rate

# age-dependent and time-dependent variables (with v and z suffix)
# A          # assets
# Cons       # consumption
# HH_nonconv # counter for non-converging households
# N          # population mass
# Sav        # end-of-period savings
# ab         # accidental bequest
# dis_tot    # disutility of labor
# cG         # public consumption
# ell        # hours worked
# gam        # conditional survival probability
# iv         # intervivo transfers
# lambda     # shadow price of assets
# notret     # not retired indicator
# p          # pension income
# pc         # price of consumption
# r          # real interest rate
# tauC       # consumption tax rate
# tauW       # wage tax rate
# taul       # lump-sum tax rate
# theta      # productivity (age-specific)
# w          # wage rate
# y          # per-period labor and pension income

#=====initialization=====#

const ncoh = tend + (nag-1)

# parameters and counters
alpha                    = 0.0
delta                    = 0.0
sigL                     = 0.0
fag                      = 0.0
rho                      = 0.0
sigma                    = 0.0
tt                       = 0.0

# age-dependent variables
parlv0                   = zeroscol(nag)
parlv1                   = zeroscol(nag)

# time-dependent variables
A                        = zerosrow(tend)
A0                       = 0.0
Cons                     = zerosrow(tend)
Cons0                    = 0.0
CG                       = zerosrow(tend)
CG0                      = 0.0
DG                       = zerosrow(tend)
DG0                      = 0.0
Div                      = zerosrow(tend)
Div0                     = 0.0
Exp                      = zerosrow(tend)
Exp0                     = 0.0
Inv                      = zerosrow(tend)
Inv0                     = 0.0
K                        = zerosrow(tend)
K0                       = 0.0
LD                       = zerosrow(tend)
LD0                      = 0.0
LS                       = zerosrow(tend)
LS0                      = 0.0
N                        = zerosrow(tend)
N0                       = 0.0
NB                       = zerosrow(tend)
NB0                      = 0.0
Nc                       = zerosrow(tend)
Nc0                      = 0.0
Nr                       = zerosrow(tend)
Nr0                      = 0.0
Nw                       = zerosrow(tend)
Nw0                      = 0.0
P                        = zerosrow(tend)
P0                       = 0.0
PB                       = zerosrow(tend)
PB0                      = 0.0
Rev                      = zerosrow(tend)
Rev0                     = 0.0
TaxF                     = zerosrow(tend)
TaxF0                    = 0.0
TFP                      = zerosrow(tend)
TFP0                     = 0.0
V                        = zerosrow(tend)
V0                       = 0.0
Y                        = zerosrow(tend)
Y0                       = 0.0
ab                       = zerosrow(tend)
ab0                      = 0.0
eda                      = zerosrow(tend)
eda0                     = 0.0
edab                     = zerosrow(tend)
edab0                    = 0.0
edg                      = zerosrow(tend)
edg0                     = 0.0
ediv                     = zerosrow(tend)
ediv0                    = 0.0
edl                      = zerosrow(tend)
edl0                     = 0.0
edw                      = zerosrow(tend)
edw0                     = 0.0
edy                      = zerosrow(tend)
edy0                     = 0.0
iv                       = zerosrow(tend)
iv0                      = 0.0
pc                       = zerosrow(tend)
pc0                      = 0.0
qTob                     = zerosrow(tend)
qTob0                    = 0.0
r                        = zerosrow(tend)
r0                       = 0.0
rag                      = zerosrow(tend)
rag0                     = 0.0
tauC                     = zerosrow(tend)
tauC0                    = 0.0
tauF                     = zerosrow(tend)
tauF0                    = 0.0
tauW                     = zerosrow(tend)
tauW0                    = 0.0
taul                     = zerosrow(tend)
taul0                    = 0.0
tauprof                  = zerosrow(tend)
tauprof0                 = 0.0
uck                      = zerosrow(tend)
uck0                     = 0.0
w                        = zerosrow(tend)
w0                       = 0.0

# age-dependent and time-dependent variables
Av                       = zerosmat(nag,tend)
Az                       = zerosmat(nag,ncoh)
Av0                      = zeroscol(nag)
Consv                    = zerosmat(nag,tend)
Consz                    = zerosmat(nag,ncoh)
Consv0                   = zeroscol(nag)
HH_nonconvv              = zerosmat(nag,tend)
HH_nonconvz              = zerosmat(nag,ncoh)
HH_nonconvv0             = zeroscol(nag)
Nv                       = zerosmat(nag,tend)
Nz                       = zerosmat(nag,ncoh)
Nv0                      = zeroscol(nag)
Savv                     = zerosmat(nag,tend)
Savz                     = zerosmat(nag,ncoh)
Savv0                    = zeroscol(nag)
abv                      = zerosmat(nag,tend)
abz                      = zerosmat(nag,ncoh)
abv0                     = zeroscol(nag)
dis_totv                 = zerosmat(nag,tend)
dis_totz                 = zerosmat(nag,ncoh)
dis_totv0                = zeroscol(nag)
cGv                      = zerosmat(nag,tend)
cGz                      = zerosmat(nag,ncoh)
cGv0                     = zeroscol(nag)
ellv                     = zerosmat(nag,tend)
ellz                     = zerosmat(nag,ncoh)
ellv0                    = zeroscol(nag)
gamv                     = zerosmat(nag,tend)
gamz                     = zerosmat(nag,ncoh)
gamv0                    = zeroscol(nag)
ivv                      = zerosmat(nag,tend)
ivz                      = zerosmat(nag,ncoh)
ivv0                     = zeroscol(nag)
lambdav                  = zerosmat(nag,tend)
lambdaz                  = zerosmat(nag,ncoh)
lambdav0                 = zeroscol(nag)
notretv                  = zerosmat(nag,tend)
notretz                  = zerosmat(nag,ncoh)
notretv0                 = zeroscol(nag)
pv                       = zerosmat(nag,tend)
pz                       = zerosmat(nag,ncoh)
pv0                      = zeroscol(nag)
pcv                      = zerosmat(nag,tend)
pcz                      = zerosmat(nag,ncoh)
pcv0                     = zeroscol(nag)
rv                       = zerosmat(nag,tend)
rz                       = zerosmat(nag,ncoh)
rv0                      = zeroscol(nag)
tauCv                    = zerosmat(nag,tend)
tauCz                    = zerosmat(nag,ncoh)
tauCv0                   = zeroscol(nag)
tauWv                    = zerosmat(nag,tend)
tauWz                    = zerosmat(nag,ncoh)
tauWv0                   = zeroscol(nag)
taulv                    = zerosmat(nag,tend)
taulz                    = zerosmat(nag,ncoh)
taulv0                   = zeroscol(nag)
thetav                   = zerosmat(nag,tend)
thetaz                   = zerosmat(nag,ncoh)
thetav0                  = zeroscol(nag)
wv                       = zerosmat(nag,tend)
wz                       = zerosmat(nag,ncoh)
wv0                      = zeroscol(nag)
yv                       = zerosmat(nag,tend)
yz                       = zerosmat(nag,ncoh)
yv0                      = zeroscol(nag)

@with_kw mutable struct ModelData

  parlv0::Matrix{Float64}       = parlv0
  parlv1::Matrix{Float64}       = parlv1

  A::Matrix{Float64}            = A
  Cons::Matrix{Float64}         = Cons
  CG::Matrix{Float64}           = CG
  DG::Matrix{Float64}           = DG
  Div::Matrix{Float64}          = Div
  Exp::Matrix{Float64}          = Exp
  Inv::Matrix{Float64}          = Inv
  K::Matrix{Float64}            = K
  LD::Matrix{Float64}           = LD
  LS::Matrix{Float64}           = LS
  N::Matrix{Float64}            = N
  NB::Matrix{Float64}           = NB
  Nc::Matrix{Float64}           = Nc
  Nr::Matrix{Float64}           = Nr
  Nw::Matrix{Float64}           = Nw
  P::Matrix{Float64}            = P
  PB::Matrix{Float64}           = PB
  Rev::Matrix{Float64}          = Rev
  TaxF::Matrix{Float64}         = TaxF
  TFP::Matrix{Float64}          = TFP
  V::Matrix{Float64}            = V
  Y::Matrix{Float64}            = Y
  ab::Matrix{Float64}           = ab
  eda::Matrix{Float64}          = eda
  edab::Matrix{Float64}         = edab
  edg::Matrix{Float64}          = edg
  ediv::Matrix{Float64}         = ediv
  edl::Matrix{Float64}          = edl
  edw::Matrix{Float64}          = edw
  edy::Matrix{Float64}          = edy
  iv::Matrix{Float64}           = iv
  pc::Matrix{Float64}           = pc
  qTob::Matrix{Float64}         = qTob
  r::Matrix{Float64}            = r
  rag::Matrix{Float64}          = rag
  tauC::Matrix{Float64}         = tauC
  tauF::Matrix{Float64}         = tauF
  tauW::Matrix{Float64}         = tauW
  taul::Matrix{Float64}         = taul
  tauprof::Matrix{Float64}      = tauprof
  uck::Matrix{Float64}          = uck
  w::Matrix{Float64}            = w

  Av::Matrix{Float64}           = Av
  Az::Matrix{Float64}           = Az
  Av0::Matrix{Float64}          = Av0
  Consv::Matrix{Float64}        = Consv
  Consz::Matrix{Float64}        = Consz
  Consv0::Matrix{Float64}       = Consv0
  HH_nonconvv::Matrix{Float64}  = HH_nonconvv
  HH_nonconvz::Matrix{Float64}  = HH_nonconvz
  HH_nonconvv0::Matrix{Float64} = HH_nonconvv0
  Nv::Matrix{Float64}           = Nv
  Nz::Matrix{Float64}           = Nz
  Nv0::Matrix{Float64}          = Nv0
  Savv::Matrix{Float64}         = Savv
  Savz::Matrix{Float64}         = Savz
  Savv0::Matrix{Float64}        = Savv0
  abv::Matrix{Float64}          = abv
  abz::Matrix{Float64}          = abz
  abv0::Matrix{Float64}         = abv0
  dis_totv::Matrix{Float64}     = dis_totv
  dis_totz::Matrix{Float64}     = dis_totz
  dis_totv0::Matrix{Float64}    = dis_totv0
  cGv::Matrix{Float64}          = cGv
  cGz::Matrix{Float64}          = cGz
  cGv0::Matrix{Float64}         = cGv0
  ellv::Matrix{Float64}         = ellv
  ellz::Matrix{Float64}         = ellz
  ellv0::Matrix{Float64}        = ellv0
  gamv::Matrix{Float64}         = gamv
  gamz::Matrix{Float64}         = gamz
  gamv0::Matrix{Float64}        = gamv0
  ivv::Matrix{Float64}          = ivv
  ivz::Matrix{Float64}          = ivz
  ivv0::Matrix{Float64}         = ivv0
  lambdav::Matrix{Float64}      = lambdav
  lambdaz::Matrix{Float64}      = lambdaz
  lambdav0::Matrix{Float64}     = lambdav0
  notretv::Matrix{Float64}      = notretv
  notretz::Matrix{Float64}      = notretz
  notretv0::Matrix{Float64}     = notretv0
  pv::Matrix{Float64}           = pv
  pz::Matrix{Float64}           = pz
  pv0::Matrix{Float64}          = pv0
  pcv::Matrix{Float64}          = pcv
  pcz::Matrix{Float64}          = pcz
  pcv0::Matrix{Float64}         = pcv0
  rv::Matrix{Float64}           = rv
  rz::Matrix{Float64}           = rz
  rv0::Matrix{Float64}          = rv0
  tauCv::Matrix{Float64}        = tauCv
  tauCz::Matrix{Float64}        = tauCz
  tauCv0::Matrix{Float64}       = tauCv0
  tauWv::Matrix{Float64}        = tauWv
  tauWz::Matrix{Float64}        = tauWz
  tauWv0::Matrix{Float64}       = tauWv0
  taulv::Matrix{Float64}        = taulv
  taulz::Matrix{Float64}        = taulz
  taulv0::Matrix{Float64}       = taulv0
  thetav::Matrix{Float64}       = thetav
  thetaz::Matrix{Float64}       = thetaz
  thetav0::Matrix{Float64}      = thetav0
  wv::Matrix{Float64}           = wv
  wz::Matrix{Float64}           = wz
  wv0::Matrix{Float64}          = wv0
  yv::Matrix{Float64}           = yv
  yz::Matrix{Float64}           = yz
  yv0::Matrix{Float64}          = yv0
end

;