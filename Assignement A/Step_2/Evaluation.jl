# Include the evaluate policy file
include("evaluate_policy.jl")
using JLD

# Different method 
nb_methods = 2
Methods = ["kmeans" "kmedoids"]
# "FastForward"

# Different number of day ahead strategy 
nb_day_aheads = 4
Day_aheads = [1 2 3 4 ; 120 80 60 48]

# Different number of generated scenarios
nb_scenarios = 3
Scenarios_generation = [500 1000 5000]

# Different granularity
nb_granularity = 3
Granularities = [1 0.5 0.1]

# Final_cost = Dict()
# Final_time = Dict()

# ##############################
# # Comparison of the parameters
# ##############################

# # 4 loops together
# for m in Methods 
#     for d in 1:nb_day_aheads
#         look_ahead_days = Day_aheads[1,d]
#         nb_reduced_scenarios = Day_aheads[2,d]
#         for s in Scenarios_generation
#             for g in Granularities
#                 println("Methods $m with $look_ahead_days ahead days, $s scenarios generated and a granularity of $g")
#                 Final_cost[(m,look_ahead_days,nb_reduced_scenarios,s,g)], Final_time[(m,look_ahead_days,nb_reduced_scenarios,s,g)] = evaluate_policy("ST", m, look_ahead_days, nb_reduced_scenarios, s, g)
#             end
#         end
#     end
# end

# # Output the values in a file
# script_directory = @__DIR__
# file_path = joinpath(script_directory, "results_comparison.jld")
# save(file_path, "Final_cost", Final_cost, "Final_time", Final_time)

# Load back the result file 
script_directory = @__DIR__
file_path = joinpath(script_directory, "results_comparison.jld")
dict = load(file_path)
Final_cost, Final_time = dict["Final_cost"], dict["Final_time"]

# Put everything in a 4 dimensions table
Table_cost = zeros(nb_methods, nb_day_aheads, nb_scenarios, nb_granularity)
Table_time = zeros(nb_methods, nb_day_aheads, nb_scenarios, nb_granularity)
for m in 1:nb_methods 
    for d in 1:nb_day_aheads
        look_ahead_days = Day_aheads[1,d]
        nb_reduced_scenarios = Day_aheads[2,d]
        for s in 1:nb_scenarios
            for g in 1:nb_granularity
                Table_cost[m,d,s,g] = Final_cost[(Methods[m], look_ahead_days, nb_reduced_scenarios, Scenarios_generation[s], Granularities[g])]
                Table_time[m,d,s,g] = Final_time[(Methods[m], look_ahead_days, nb_reduced_scenarios, Scenarios_generation[s], Granularities[g])]
            end
        end
    end
end

# Find the minimum and the index
Min_cost = minimum(Table_cost)
Argmin_cost = argmin(Table_cost)
Min_time = minimum(Table_time)
Argmin_time = argmin(Table_time)

##############################
# Comparison of the policy
##############################

# evaluate_policy("ST", "kmeans", 3, 60, 1000, 0.5)
# evaluate_policy("EV", "kmeans", 3, 60, 1000, 0.5)
# evaluate_policy("OiH", "kmeans", 3, 60, 1000, 0.5)