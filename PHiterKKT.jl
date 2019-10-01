using Pkg, DataFrames
using JuMP, Gurobi
using Distributions

function PHiterKKT(testsystem, K, s, rho, w, mm, mmm, E_fp, E_gp, E_lambda, gpmax_old, timeset, PdownMax, UserGap, η, RPsplit, σ0_pv, σ0_wind, υ0)
   ScaleParam = 1e-3
   ScaleObj = 1e-5

   include("NetworkDataType.jl")
   include("NetworkLoad.jl")
   include("Demand24Load.jl")
   filename_Load24 = string(pwd(),"/data/load24data.csv")
   loadcurve = Demand24Load(filename_Load24)

   include("GenCurveLoad.jl")
   filename_windcurve = string(pwd(),"/data/hourly_wind_gen_2018.csv")
   filename_pvcurve = string(pwd(),"/data/hourly_solar_gen_2018.csv")
   pvcurve = GenCurveLoad(filename_pvcurve)[2:31:end,:]/2883.8  # pick one day a month over the year / nameplate capacity = 2883MW
   windcurve = GenCurveLoad(filename_windcurve)[2:31:end,:]/1300 # pick one day a month over the year / nameplate capacity = 1300MW

   stateset = 1:6
   loadzone = Dict()
   loadzone[1] = 1
   loadzone[2] = 2
   loadzone[3] = 3
   loadzone[4] = [4, 5, 8] # WCMA, NEMA, SEMA
   loadzone[5] = 6
   loadzone[6] = 7

   buses, lines, generators, datamat = NetworkLoad(testsystem)
   lineset = 1:length(lines)
   busset = 1:length(buses)
   genset = 1:length(generators)
   Dp = 2*datamat["bus"][:,2]

   windset = Int64[]
   pvset = Int64[]
   hydroset = Int64[]
   nucset = Int64[]
   coalset = Int64[]
   oilset = Int64[]
   gasset = Int64[]

   for g in genset
      if datamat["gen"][g,2] == "Wind"
         push!(windset, generators[g].gindex)
      elseif datamat["gen"][g,2] == "PV"
         push!(pvset, generators[g].gindex)
      elseif datamat["gen"][g,2] == "Hydro"
         push!(hydroset, generators[g].gindex)
      elseif datamat["gen"][g,2] == "Nuclear"
         push!(nucset, generators[g].gindex)
      elseif datamat["gen"][g,2] == "Coal"
         push!(coalset, generators[g].gindex)
      elseif datamat["gen"][g,2] == "Oil"
         push!(oilset, generators[g].gindex)
      elseif datamat["gen"][g,2] == "Gas"
         push!(gasset, generators[g].gindex)
      end
   end

   for g in genset
      generators[g].Pmin = 0
   end

   genset_hat = genset[end]+1 : genset[end]+3*length(busset)
   for i in 1:length(busset)
      push!(windset, genset[end]+3*(i-1)+1)
      push!(pvset, genset[end]+3*(i-1)+2)
      push!(gasset, genset[end]+3*(i-1)+3)
   end

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

   # Investment cost source:  https://www.eia.gov/electricity/generatorcosts/
   Cinv = Dict()
   DCRF = 3.5481e-4 # h=10years / 5%

   Cinv["Wind"] = DCRF * 1630 * 1e3
   Cinv["PV"] = DCRF * 2434 * 1e3
   Cinv["Hydro"] = DCRF * 5312 * 1e3
   Cinv["Nuclear"] = DCRF * 1000 * 1e3
   Cinv["Coal"] = DCRF * 1000 * 1e3
   Cinv["Oil"] = DCRF * 1672 * 1e3
   Cinv["Gas"] = DCRF * 895 * 1e3

   Coper = Dict()
   Coper["Wind"] = [0, 2, 0]
   Coper["PV"] = [0, 3, 0]
   Coper["Hydro"] = [0, 4, 0]
   Coper["Nuclear"] = [0, 7, 1250]
   Coper["Coal"] = [0.0, 19.98, 1535.0]
   Coper["Oil"] = [0.0, 192.06, 5000.0]
   Coper["Gas"] = [0.0, 20.0, 612.84]

   Cinvg = Vector{Float64}(undef, length(union(genset,genset_hat)))
   Coperg = Vector{Vector{Float64}}(undef, length(union(genset,genset_hat)))
   for g in genset_hat
      if g in windset
         Cinvg[g] = Cinv["Wind"]
         Coperg[g] = Coper["Wind"]
      elseif g in pvset
         Cinvg[g] = Cinv["PV"]
         Coperg[g] = Coper["PV"]
      elseif g in hydroset
         Cinvg[g] = Cinv["Hydro"]
         Coperg[g] = Coper["Hydro"]
      elseif g in nucset
         Cinvg[g] = Cinv["Nuclear"]
         Coperg[g] = Coper["Nuclear"]
      elseif g in coalset
         Cinvg[g] = Cinv["Coal"]
         Coperg[g] = Coper["Coal"]
      elseif g in oilset
         Cinvg[g] = Cinv["Oil"]
         Coperg[g] = Coper["Oil"]
      elseif g in gasset
         Cinvg[g] = Cinv["Gas"]
         Coperg[g] = Coper["Gas"]
      end
   end

   for g in genset
      if g in windset
         Cinvg[g] = Cinv["Wind"]
         Coperg[g] = Coper["Wind"]
      elseif g in pvset
         Cinvg[g] = Cinv["PV"]
         Coperg[g] = Coper["PV"]
      elseif g in hydroset
         Cinvg[g] = Cinv["Hydro"]
         Coperg[g] = Coper["Hydro"]
      elseif g in nucset
         Cinvg[g] = Cinv["Nuclear"]
         Coperg[g] = generators[g].cost
      elseif g in coalset
         Cinvg[g] = Cinv["Coal"]
         Coperg[g] = generators[g].cost
      elseif g in oilset
         Cinvg[g] = Cinv["Oil"]
         Coperg[g] = generators[g].cost
      elseif g in gasset
         Cinvg[g] = Cinv["Gas"]
         Coperg[g] = generators[g].cost
      end
   end

   ρ = 1.0*ones(length(union(genset, genset_hat)),timeset[end])
   for g in pvset
      for t in timeset
         ρ[g,t] = pvcurve[1,t]
      end
   end
   for g in windset
      for t in timeset
         ρ[g,t] = windcurve[1,t]
      end
   end

   σ = zeros(length(union(genset, genset_hat)))
   for g in pvset
      σ[g] = σ0_pv
   end
   for g in windset
      σ[g] = σ0_wind
   end


   υ = υ0*ones(length(union(genset, genset_hat)), timeset[end])
   inv_ϕ = quantile.(Normal(), 1-η)

   #---------------------------------------------------------------------------------#
   #|RPS|  Maine;  NH; Vermont;WCMA;  NEMA;  CT;    RI;   SEMA
   # kappa = [1; 0.252; 0.75; 0.411; 0.411; 0.44; 0.385; 0.411]
   #|RPS|  Maine; NH;  VT;   MA;    CT;   RI
   kappa = [1; 0.252; 0.75; 0.411; 0.44; 0.385]
   # kappa = zeros(8)
   InvB = 1e13  * ones(length(stateset)) * DCRF# Investment Budget
   PolicyB = 1e8  * ones(length(stateset)) # Policy Budget
   PET = 3*ones(length(stateset), timeset[end])
   PCT = 3*ones(length(stateset))
   Dπ = ones(length(busset), timeset[end])

   Γ = 0
   I_R = intersect(genset, union(windset,pvset))
   I_R_hat = intersect(genset_hat, union(windset,pvset))
   I_C = setdiff(genset, union(windset,pvset))
   I_C_hat = setdiff(genset_hat,union(windset,pvset))

   α_Newg = zeros(length(union(genset, genset_hat)))
   for nnn in loadzone[s]
      for i in intersect(S_gs[s], I_C)
         α_Newg[i] = generators[i].Pmax/sum(generators[j].Pmax for j in intersect(S_gs[s], I_C))
      end
   end

   #---------------------------------------------------------------------------------#
   m = Model(with_optimizer(Gurobi.Optimizer, MIPGap=UserGap, NumericFocus = 3))
   #---------------------------------------------------------------------------------#
   @variables(m, begin
      gpbar[g in union(genset, genset_hat), t in timeset] >= 0
      gpmax[g in genset_hat] >= 0
      pdown[n in busset, t in timeset]
      x[n in busset, t in timeset, j in intersect(B_gn[n], union(I_R, I_R_hat))] >= 0
      yLB[g in union(I_C, I_C_hat), t in timeset] >= 0
      yUB[g in union(I_C, I_C_hat), t in timeset] >= 0
      Cinv_total[n in busset] >= 0
   end)
   #---------------------------------------------------------------------------------#
   @variables(m, begin
      fp[l in lineset, t in timeset]
      gp[g in union(genset, genset_hat), t in timeset] >= 0
      theta[n in busset, t in timeset]
      xi[l in lineset, t in timeset]
      lambda[n in busset, t in timeset]
      gammaLB[g in union(genset, genset_hat), t in timeset] >= 0
      gammaUB[g in union(genset, genset_hat), t in timeset] >= 0
      deltaLB[l in lineset, t in timeset] >= 0
      deltaUB[l in lineset, t in timeset] >= 0
      OCgen[g in union(genset, genset_hat), t in timeset]
      Obj_KKT
   end)
   #---------------------------------------------------------------------------------#
   lambdaLB = 0
   lambdaUB = 700
   M = 2^K - 1
   delPdown = (PdownMax-(-PdownMax))/M
   @variable(m, s1[n in loadzone[s], t in timeset, k in 1:K])
   @variable(m, z[n in loadzone[s], t in timeset, k in 1:K], Bin)
   #---------------------------------------------------------------------------------#
   @constraint(m, Linearization[n in loadzone[s], t in timeset], pdown[n,t] == - PdownMax + delPdown * sum(2^(k-1) * z[n,t,k] for k in 1:K))
   @constraint(m, BigM1[n in loadzone[s], t in timeset, k in 1:K], lambdaLB * z[n,t,k] <= s1[n,t,k])
   @constraint(m, BigM2[n in loadzone[s], t in timeset, k in 1:K], lambdaUB * z[n,t,k] >= s1[n,t,k])
   @constraint(m, BigM3[n in loadzone[s], t in timeset, k in 1:K], lambda[n,t] - lambdaUB * (1-z[n,t,k]) <= s1[n,t,k])
   @constraint(m, BigM4[n in loadzone[s], t in timeset, k in 1:K], lambda[n,t] - lambdaLB * (1-z[n,t,k]) >= s1[n,t,k])
   @constraint(m, MPlimitLB[n in busset, t in timeset], lambda[n,t] >= 0.0)

   @objective(m, Max,
      ScaleObj*
      (sum(
         sum(Dπ[n,t] * Dp[n] * loadcurve[n,t] - (lambda[n,t]*(-PdownMax)+delPdown*sum(2^(k-1) * s1[n,t,k] for k in 1:K)) for n in loadzone[s])
         - PET[s,t] * sum(gpbar[i,t] for i in intersect(union(I_R, I_R_hat), S_gs[s]))
         - sum(Coperg[i][2] * gpbar[i,t] for i in intersect(union(I_R, I_R_hat), S_gs[s]))
         - sum(Coperg[i][2] * (gpbar[i,t] - α_Newg[i]*υ[i,t]) for i in intersect(union(I_C, I_C_hat), S_gs[s]))
      for t in timeset)
      - PCT[s] * sum(gpmax[i] for i in intersect(I_R_hat, S_gs[s]))
      - sum(Cinv_total[n] for n in loadzone[s])
      - sum(mm[i,t] * gp[i,timeset[t]] for i in union(genset, genset_hat), t in 1:length(timeset))
      - sum(mmm[n,t] * lambda[n,timeset[t]] for n in busset, t in 1:length(timeset))
      - rho/2 * sum((gp[i,timeset[t]] - E_gp[i,t])^2 for t in 1:length(timeset), i in union(genset, genset_hat))
      - rho/2 * sum((lambda[n,timeset[t]] - E_lambda[n,t])^2 for t in 1:length(timeset), n in busset)
      )
      )

   @constraint(m, Cinv_def[n in busset], Cinv_total[n] ==
      sum(Cinvg[g] * gpmax[g] for g in intersect(genset_hat, B_gn[n])))

   @constraint(m, gpbarDef1[g in intersect(I_R, S_gs[s]), t in timeset], gpbar[g,t] == ρ[g,t] * generators[g].Pmax + υ[g,t])
   @constraint(m, gpbarDef2[g in intersect(I_R_hat, S_gs[s]), t in timeset], gpbar[g,t] == ρ[g,t] * gpmax[g] + υ[g,t])

   @constraint(m, xdef1[t in timeset, j in intersect(S_gs[s], I_R)], x[B_g[j],t,j] == σ[j] * generators[j].Pmax)
   @constraint(m, xdef2[t in timeset, j in intersect(S_gs[s], I_R_hat)], x[B_g[j],t,j] == σ[j] * gpmax[j])

   @constraint(m, yLBdef1[i in intersect(I_C, S_gs[s]),     t in timeset], yLB[i,t] == 1/(inv_ϕ*α_Newg[i])*(gpbar[i,t] - generators[i].Pmin))
   @constraint(m, yUBdef1[i in intersect(I_C, S_gs[s]),     t in timeset], yUB[i,t] == 1/(inv_ϕ*α_Newg[i])*(generators[i].Pmax - gpbar[i,t]))

   if RPsplit == 1
      @constraint(m, CCreform1[i in intersect(I_C, S_gs[s]), t in timeset], sum(x[B_g[i],t,j]^2 for j in intersect(B_gn[B_g[i]], union(I_R, I_R_hat))) <= yLB[i,t]^2)
      @constraint(m, CCreform2[i in intersect(I_C, S_gs[s]), t in timeset], sum(x[B_g[i],t,j]^2 for j in intersect(B_gn[B_g[i]], union(I_R, I_R_hat))) <= yUB[i,t]^2)
   else RPsplit == 0
      @constraint(m, CCreform1[i in intersect(I_C, S_gs[s]), t in timeset], sum(x[B_g[j],t,j]^2 for j in intersect(S_gs[S_g[i]], union(I_R, I_R_hat))) <= yLB[i,t]^2)
      @constraint(m, CCreform2[i in intersect(I_C, S_gs[s]), t in timeset], sum(x[B_g[j],t,j]^2 for j in intersect(S_gs[S_g[i]], union(I_R, I_R_hat))) <= yUB[i,t]^2)
   end
   @constraint(m, SupplyDemandBalance[t in timeset, n in loadzone[s]], sum(gpbar[g,t] for g in B_gn[n]) + pdown[n,t] == Dp[n] * loadcurve[n,t])
   @constraint(m, PdownLB[n in loadzone[s], t in timeset], pdown[n,t] + PdownMax >= 0) #zetaLB
   @constraint(m, PdownUB[n in loadzone[s], t in timeset], PdownMax - pdown[n,t] >= 0) #zetaUB
   @constraint(m, RPS, sum(sum(gpbar[g,t] for g in intersect(S_gs[s], union(I_R, I_R_hat))) for t in timeset) >= sum(sum(kappa[s] * Dp[n] * loadcurve[n,t] for t in timeset) for n in loadzone[s]))

   @constraint(m, InvBudget, ScaleParam*sum(Cinv_total[n] for n in loadzone[s]) <= ScaleParam*InvB[s])
   @constraint(m, PolicyBudget,
      ScaleParam*(sum(PET[s,t]*sum(gpbar[g,t] for g in intersect(S_gs[s], union(I_R, I_R_hat))) for t in timeset)
      + PCT[s]*sum(gpmax[g] for g in intersect(S_gs[s], I_R_hat)))
      <= ScaleParam*PolicyB[s])
   #---------------------------------------------------------------------------------#
   @constraint(m, DCPF[l in lineset, t in timeset], fp[l,t] == (1/lines[l].x) * (theta[lines[l].fbus,t] - theta[lines[l].tbus,t]))
   @constraint(m, slack[t in timeset], theta[1,t] == 0)
   @constraint(m, NodeBalance[n in busset, t in timeset], sum(gp[g,t] for g in B_gn[n]) + sum(fp[l,t] for l in buses[n].inline)
      - sum(fp[l,t] for l in buses[n].outline) == Dp[n] * loadcurve[n,t])
   @constraint(m, dual_gp[g in union(genset, genset_hat), t in timeset], - Coperg[g][2] + lambda[B_g[g],t] + gammaLB[g,t] - gammaUB[g,t] == 0)
   @constraint(m, dual_fp[l in lineset, t in timeset], xi[l,t] + lambda[lines[l].tbus,t] - lambda[lines[l].fbus,t] + deltaLB[l,t] - deltaUB[l,t] == 0)
   @constraint(m, dual_theta[n in busset, t in timeset], - sum(xi[l,t]/lines[l].x for l in buses[n].outline) + sum(xi[l,t]/lines[l].x for l in buses[n].inline) == 0)
   #---------------------------------------------------------------------------------#
   @constraint(m, gpmax_old_put[g in setdiff(genset_hat, S_gs[s])], gpmax[g] == gpmax_old[g-length(genset)])
   #---------------------------------------------------------------------------------#
   @constraint(m, gpLB2[g in genset, t in timeset], gp[g,t] >= generators[g].Pmin)
   @constraint(m, gpUB2_C[g in I_C, t in timeset], gp[g,t] <= generators[g].Pmax)
   @constraint(m, gpUB2_R[g in I_R, t in timeset], gp[g,t] <= ρ[g,t] * generators[g].Pmax + υ[g,t])
   @constraint(m, gpLB3[g in genset_hat, t in timeset], gp[g,t] >= Γ*gpmax[g])
   @constraint(m, gpUB3_C[g in I_C_hat, t in timeset], gp[g,t] <= gpmax[g])
   @constraint(m, gpUB3_R[g in I_R_hat, t in timeset], gp[g,t] <= ρ[g,t] * gpmax[g] + υ[g,t])
   @constraint(m, LineCapLB[l in lineset, t in timeset], fp[l,t] >= - lines[l].u)
   @constraint(m, LineCapUB[l in lineset, t in timeset], fp[l,t] <= + lines[l].u)
   #----------CS-----------#
   @constraint(m, CS_gpLB2[g in genset, t in timeset], [ρ[g,t] * gp[g,t] + υ[g,t] - generators[g].Pmin, gammaLB[g,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_gpUB2_C[g in I_C, t in timeset], [generators[g].Pmax - gp[g,t], gammaUB[g,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_gpUB2_R[g in I_R, t in timeset], [ρ[g,t] * generators[g].Pmax + υ[g,t] - gp[g,t], gammaUB[g,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_gpLB3[g in genset_hat, t in timeset], [ρ[g,t] * gp[g,t] + υ[g,t] - Γ*gpmax[g], gammaLB[g,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_gpUB3_C[g in I_C_hat, t in timeset], [gpmax[g] - gp[g,t], gammaUB[g,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_gpUB3_R[g in I_R_hat, t in timeset], [ρ[g,t] * gpmax[g] + υ[g,t] - gp[g,t], gammaUB[g,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_LineCapLB[l in lineset, t in timeset], [fp[l,t] + lines[l].u, deltaLB[l,t]] in MOI.SOS1([1.0,2.0]))
   @constraint(m, CS_LineCapUB[l in lineset, t in timeset], [lines[l].u - fp[l,t], deltaUB[l,t]] in MOI.SOS1([1.0,2.0]))
   #---------------------------#
   @constraint(m, GenCost[g in union(genset, genset_hat), t in timeset], OCgen[g,t] == Coperg[g][2]*gp[g,t])
   @constraint(m, Obj_KKT == - sum(sum(OCgen[g,t] for g in union(genset, genset_hat)) for t in timeset))
   ##########################################################################################
   optimize!(m)
   return m
end
