using Pkg, DataFrames
using JuMP, Gurobi
using Distributions, LinearAlgebra
#################################################################################
rho = 0.7e-3 #obj scale 1e-3
timeset = 1:24
regset = 1:6 # 6states
PdownMax = 2400
UserGapIni = 1e-2
UserGap = 1e-2
Niter = 1000
K = 10
η = 0.03
RPsplit = 0
σ0_pv = 0.2
σ0_wind = 0.1
υ0 = 0.0
#################################################################################
include("NetworkDataType.jl")
include("NetworkLoad.jl")
include("PHiniKKT.jl")
include("PHiterKKT.jl")
include("Demand24Load.jl")
include("GenCurveLoad.jl")

testsystem = "ISONE8busTN"

stateset = 1:6
loadzone = Dict()
loadzone[1] = 1
loadzone[2] = 2
loadzone[3] = 3
loadzone[4] = [4, 5, 8] # WCMA, NEMA, SEMA
loadzone[5] = 6
loadzone[6] = 7


filename_Load24 = string(pwd(),"/data/load24data.csv")
loadcurve = Demand24Load(filename_Load24)

buses, lines, generators, datamat = NetworkLoad(testsystem)
lineset = 1:length(lines)
busset = 1:length(buses)
genset = 1:length(generators)
genset_hat = genset[end]+1 : genset[end]+3*length(busset)
Dp = 2*datamat["bus"][:,2]

B_g = []
for g in genset
  push!(B_g, generators[g].location)
end
for i in busset
  for j in 1:3
     push!(B_g, i)
  end
end

B_gn = Array{Array{Int64},1}(undef,length(buses))
for ii in 1:length(buses)
  B_gn[ii] = Int64[]
end
for g in union(genset,genset_hat)
  push!(B_gn[B_g[g]], g)
end


S_g = Dict()
for s in stateset
   for z in loadzone[s]
      for g in B_gn[z]
         S_g[g] = s
      end
   end
end

S_gs = Array{Array{Int64},1}(undef,length(stateset))
for ii in 1:length(stateset)
   S_gs[ii] = Int64[]
end
for s in stateset
   for n in loadzone[s]
      for g in B_gn[n]
         push!(S_gs[s], g)
      end
   end
end
#################################################################################
fp_ini = Array{Float64}(undef, length(lineset), length(timeset), length(busset))
gpbar_ini = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset))
gp_ini = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset))
lambda_ini = Array{Float64}(undef, length(busset), length(timeset), length(busset))
theta_ini = Array{Float64}(undef, length(busset), length(timeset), length(busset))
gpmax_ini = Array{Float64}(undef, length(genset_hat), length(busset))
status_ini = Vector{}(undef, length(busset))
solvetime_ini = Vector{Float64}(undef, length(busset))

for n in regset
   global temp = PHiniKKT(testsystem, n, K, timeset, PdownMax, UserGapIni, η, RPsplit, σ0_pv, σ0_wind, υ0)
   global status_ini[n] = termination_status(temp)
   global solvetime_ini[n] = MOI.get(temp, MOI.SolveTime())
   #
   global fp_ini[:,:,n] = value.(temp[:fp]).data
   global gpbar_ini[:,:,n] = value.(temp[:gpbar]).data
   global gp_ini[:,:,n] = value.(temp[:gp]).data
   global lambda_ini[:,:,n] = value.(temp[:lambda]).data
   global theta_ini[:,:,n] = value.(temp[:theta]).data
   global gpmax_ini[:,n] = value.(temp[:gpmax]).data
end
#################################################################################
fp = Array{Float64}(undef, length(lineset), length(timeset), length(busset))
gpbar = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset))
gp = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset))
lambda = Array{Float64}(undef, length(busset), length(timeset), length(busset))
theta = Array{Float64}(undef, length(busset), length(timeset), length(busset))
gpmax = Array{Float64}(undef, length(genset_hat), length(busset))
status = Vector{Any}(undef, length(busset))
solvetime = Vector{Float64}(undef, length(busset))
#################################################################################
fp_iter = Array{Float64}(undef, length(lineset), length(timeset), length(busset), Niter)
pdown_iter = Array{Float64}(undef, length(busset), length(timeset), length(busset), Niter)
pdown_settled_iter = Array{Float64}(undef, length(busset), length(timeset), length(busset), Niter)
gpbar_iter = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset), Niter)
gp_iter = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset), Niter)
lambda_iter = Array{Float64}(undef, length(busset), length(timeset), length(busset), Niter)
theta_iter = Array{Float64}(undef, length(busset), length(timeset), length(busset), Niter)
gpmax_iter = Array{Float64}(undef, length(genset_hat), length(busset), Niter)
status_iter = Array{Any}(undef, length(busset), Niter)
solvetime_iter = Array{Float64}(undef, length(busset), Niter)
objvalue_iter = Array{Float64}(undef, length(busset), Niter)
#
conv = fill(77.7, Niter)
penalty1_iter = Array{Float64}(undef, length(busset), Niter)
penalty2_iter = Array{Float64}(undef, length(busset), Niter)
E_fp_iter = Array{Float64}(undef, length(lineset), length(timeset), length(busset), Niter)
E_gp_iter = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), Niter)
E_lambda_iter = Array{Float64}(undef, length(busset), length(timeset), length(busset), Niter)
w_iter = zeros(length(lineset), length(timeset), length(busset), Niter)
m_iter = zeros(length(union(genset, genset_hat)), length(timeset), length(busset), Niter)
mm_iter = zeros(length(busset), length(timeset), length(busset), Niter)
w_old = zeros(length(lineset), length(timeset), length(busset))
m_old = zeros(length(union(genset, genset_hat)), length(timeset), length(busset))
mm_old = zeros(length(busset), length(timeset), length(busset))
w_new = zeros(length(lineset), length(timeset), length(busset))
m_new = zeros(length(union(genset, genset_hat)), length(timeset), length(busset))
mm_new = zeros(length(busset), length(timeset), length(busset))
#
fp_old = Array{Float64}(undef, length(lineset), length(timeset), length(busset))
gpbar_old = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset))
gp_old = Array{Float64}(undef, length(union(genset, genset_hat)), length(timeset), length(busset))
lambda_old = Array{Float64}(undef, length(busset), length(timeset), length(busset))
theta_old = Array{Float64}(undef, length(busset), length(timeset), length(busset))
gpmax_old = Array{Float64}(undef, length(genset_hat))
status_old = Vector{}(undef, length(busset))
solvetime_old = Vector{Float64}(undef, length(busset))
#################################################################################
for s in regset
   for i in intersect(S_gs[s], genset_hat)
      global gpmax_old[i-length(genset)] = gpmax_ini[i-length(genset),s]
   end
end
gpmax_old[findall(x->abs(x)<=1e-8, gpmax_old)].=0
fp_old[:,:,:] = fp_ini[:,:,:]
gp_old[:,:,:] = gp_ini[:,:,:]
lambda_old[:,:,:] = lambda_ini[:,:,:]
E_fp = sum(fp_old[:,:,n] for n in regset) / length(regset)
E_gp = sum(gp_old[:,:,n] for n in regset) / length(regset)
E_lambda = sum(lambda_old[:,:,n] for n in regset) / length(regset)
for kk in 1:Niter
   println("###iteration begin###")
   for n in busset
      global m_new[:,:,n] = m_old[:,:,n] + rho * (gp_old[:,:,n] - E_gp[:,:])
      global w_new[:,:,n] = w_old[:,:,n] + rho * (fp_old[:,:,n] - E_fp[:,:])
      global mm_new[:,:,n] = mm_old[:,:,n] + rho * (lambda_old[:,:,n] - E_lambda[:,:])
   end
   for n in regset
      global temp = PHiterKKT(testsystem, K, n, rho, w_new[:,:,n], m_new[:,:,n], mm_new[:,:,n], E_fp, E_gp, E_lambda, gpmax_old, timeset, PdownMax, UserGap, η, RPsplit, σ0_pv, σ0_wind, υ0);
      global fp_iter[:,:,n,kk] = value.(temp[:fp]).data[:,:]
      global pdown_iter[:,:,n,kk] = value.(temp[:pdown]).data[:,:]
      for n1 in busset, t1 in 1:length(timeset)
         global pdown_settled_iter[n1,t1,n,kk] = Dp[n1] * loadcurve[n1,timeset[t1]] - sum(value.(temp[:gp])[i,timeset[t1]] for i in B_gn[n1])
      end
      global gp_iter[:,:,n,kk] = value.(temp[:gp]).data[:,:]
      global gpbar_iter[:,:,n,kk] = value.(temp[:gpbar]).data[:,:]
      global lambda_iter[:,:,n,kk] = value.(temp[:lambda]).data[:,:]
      global theta_iter[:,:,n,kk] = value.(temp[:theta]).data[:,:]
      global gpmax_iter[:,n,kk] = value.(temp[:gpmax]).data[:]
      global status_iter[n,kk] = termination_status(temp)
      global solvetime_iter[n,kk] = MOI.get(temp, MOI.SolveTime())
      global objvalue_iter[n,kk] = MOI.get(temp, MOI.ObjectiveValue())
      for i in intersect(S_gs[n], genset_hat)
         global gpmax_old[i-length(genset)] = gpmax_iter[i-length(genset),n,kk]
      end
      global gpmax_old[findall(x->abs(x)<=1e-8, gpmax_old)].=0
   end
   global fp_old[:,:,:] = fp_iter[:,:,:,kk]
   global gp_old[:,:,:] = gp_iter[:,:,:,kk]
   global lambda_old[:,:,:] = lambda_iter[:,:,:,kk]
   global E_fp = sum(fp_old[:,:,n] for n in regset) / length(regset)
   global E_gp = sum(gp_old[:,:,n] for n in regset) / length(regset)
   global E_lambda = sum(lambda_old[:,:,n] for n in regset) / length(regset)

   global w_old = w_new
   global m_old = m_new
   global mm_old = mm_new
   global conv[kk] = sum(norm(fp_old[:,:,n] - E_fp[:,:]) for n in regset) + sum(norm(gp_old[:,:,n] - E_gp[:,:]) for n in regset) + sum(norm(lambda_old[:,:,n] - E_lambda[:,:]) for n in regset)

   println("iteration #:")
   display(kk)
   println()
   println("convergence:")
   display(conv[kk])
   println()
   println("Simulation time:")
   display(solvetime_iter[regset,kk])
   println()
   println("###iteration end###")
   println("--------------------")

end
