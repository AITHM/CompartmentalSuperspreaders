using CSV, DataFrames, Random, CompartmentalSuperspreaders

Random.seed!(123)

offspring_datasets =  [chop(dat, tail=4) for dat in readdir(".\\data\\offspring\\") if last(dat, 3) == "csv"]
clinical_datasets =  [chop(dat, tail=4) for dat in readdir(".\\data\\offspring\\clinical\\") if last(dat, 3) == "csv"]


# Run baseline parameters

print("Running baseline parameters ...")

include("params\\priors_baseline.jl")
include("params\\k_baseline.jl")
include("params\\c_baseline.jl")

offspring_out = map(x -> fit_ensemble(x, α_table, c_table, models=[:NegBin, :TwoType, :Clinical, :SEIR2, :ZeroInf, :SingleType, :SEIR1]), offspring_datasets)
CSV.write("outputs\\offspring\\baseline\\parm_summary.csv", vcat(getfield.(offspring_out, 1)...))
CSV.write("outputs\\offspring\\baseline\\score_summary.csv", vcat(getfield.(offspring_out, 2)...))
CSV.write("outputs\\offspring\\baseline\\full_chain.csv", vcat(getfield.(offspring_out, 3)...))
CSV.write("outputs\\offspring\\baseline\\model_fit.csv", vcat(getfield.(offspring_out, 4)...))

println(" done")

# Run clinically-segregated data

print("Running clincally-segregated data ...")

include("params\\priors_baseline.jl")
include("params\\k_baseline.jl")
include("params\\c_baseline.jl")

clinical_out = map(x -> fit_ensemble(x, α_table, c_table, dir=".\\data\\offspring\\clinical\\", models=[:NegBin, :SingleType, :SEIR1]), clinical_datasets)
CSV.write("outputs\\offspring\\clinical\\parm_summary.csv", vcat(getfield.(clinical_out, 1)...))
CSV.write("outputs\\offspring\\clinical\\score_summary.csv", vcat(getfield.(clinical_out, 2)...))
CSV.write("outputs\\offspring\\clinical\\full_chain.csv", vcat(getfield.(clinical_out, 3)...))
CSV.write("outputs\\offspring\\clinical\\model_fit.csv", vcat(getfield.(clinical_out, 4)...))

println(" done")


# Run k = 2 parameters

print("Running sensitivity: k = 2 ...")

include("params\\priors_baseline.jl")
include("params\\k2.jl")
include("params\\c_baseline.jl")

offspring_out = map(x -> fit_ensemble(x, α_table, c_table, models=[:NegBin, :Clinical, :ZeroInf, :SEIR2, :SEIR1]), offspring_datasets)
CSV.write("outputs\\offspring\\k2\\parm_summary.csv", vcat(getfield.(offspring_out, 1)...))
CSV.write("outputs\\offspring\\k2\\score_summary.csv", vcat(getfield.(offspring_out, 2)...))
CSV.write("outputs\\offspring\\k2\\full_chain.csv", vcat(getfield.(offspring_out, 3)...))
CSV.write("outputs\\offspring\\k2\\model_fit.csv", vcat(getfield.(offspring_out, 4)...))

println(" done")



# Run k = 3 parameters

print("Running sensitivity: k = 3 ...")

include("params\\priors_baseline.jl")
include("params\\k3.jl")
include("params\\c_baseline.jl")

offspring_out = map(x -> fit_ensemble(x, α_table, c_table, models=[:NegBin, :Clinical, :ZeroInf, :SEIR2, :SEIR1]), offspring_datasets)
CSV.write("outputs\\offspring\\k3\\parm_summary.csv", vcat(getfield.(offspring_out, 1)...))
CSV.write("outputs\\offspring\\k3\\score_summary.csv", vcat(getfield.(offspring_out, 2)...))
CSV.write("outputs\\offspring\\k3\\full_chain.csv", vcat(getfield.(offspring_out, 3)...))
CSV.write("outputs\\offspring\\k3\\model_fit.csv", vcat(getfield.(offspring_out, 4)...))

println(" done")



# Run lower clinical estimates parameters

print("Running sensitivity: c = c_lower ...")

include("params\\priors_baseline.jl")
include("params\\k_baseline.jl")
include("params\\c_lower.jl")

offspring_out = map(x -> fit_ensemble(x, α_table, c_table, models=[:NegBin, :Clinical]), offspring_datasets)
CSV.write("outputs\\offspring\\clinical_lower\\parm_summary.csv", vcat(getfield.(offspring_out, 1)...))
CSV.write("outputs\\offspring\\clinical_lower\\score_summary.csv", vcat(getfield.(offspring_out, 2)...))
CSV.write("outputs\\offspring\\clinical_lower\\full_chain.csv", vcat(getfield.(offspring_out, 3)...))
CSV.write("outputs\\offspring\\clinical_lower\\model_fit.csv", vcat(getfield.(offspring_out, 4)...))

println(" done")

# Run upper clinical estimates parameters

print("Running sensitivity: c = c_upper ...")

include("params\\priors_baseline.jl")
include("params\\k_baseline.jl")
include("params\\c_upper.jl")

offspring_out = map(x -> fit_ensemble(x, α_table, c_table, models=[:NegBin, :Clinical]), offspring_datasets)
CSV.write("outputs\\offspring\\clinical_upper\\parm_summary.csv", vcat(getfield.(offspring_out, 1)...))
CSV.write("outputs\\offspring\\clinical_upper\\score_summary.csv", vcat(getfield.(offspring_out, 2)...))
CSV.write("outputs\\offspring\\clinical_upper\\full_chain.csv", vcat(getfield.(offspring_out, 3)...))
CSV.write("outputs\\offspring\\clinical_upper\\model_fit.csv", vcat(getfield.(offspring_out, 4)...))

println(" done")
