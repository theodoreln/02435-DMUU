using Random
using JuMP
using Gurobi
using Printf
using Clustering
using Distances
using Plots

include("two_stage_problem_data.jl")
include("price_process.jl")


# All of the scenario selection functions take the matrix of all scenarios (next_prices -> 3xN) 
# and the number the scenarios shall be reduced to (k) as arguments and returns the 
# reduced scenario matrix (3xk) and the sceanrios' new probabilities (Probs -> 1xk)

number_of_warehouses = 3
number_of_scenarios = 1000

#To try the function here
# prices=round.(10 * rand(3), digits=2)
# next_prices = Array{Float64}(undef, number_of_warehouses, number_of_scenarios)
# for w in 1:number_of_warehouses
#     for n in 1:number_of_scenarios
#         next_prices[w,n] = sample_next(prices[w])
#     end
# end    
# N=50


function kmeans_selection(next_prices, no_of_selected_scnearios)

    N = no_of_selected_scnearios
    clusters = kmeans(next_prices, N; maxiter=2000, display=:iter)
    clustered_prices = clusters.centers # get the cluster centers
        
    # Probabilities of data point belonging to cluster 
    scenario_assignments = assignments(clusters) #Assigning which scenario belongs to which cluster 
    
    Probs = zeros(N)
    for i in scenario_assignments
        Probs[i] = Probs[i] + 1/number_of_scenarios
    end
    reduced_next_prices= clustered_prices # Matrix of scenarios with only the reduced scenarios
    
    return reduced_next_prices, Probs
end

function kmedoids_selection(next_prices, no_of_selected_scnearios)
    N = no_of_selected_scnearios
    number_of_warehouses = 3
    # Calculate Euclidean Distance matrix (size: number_of_scenarios x number_of_scenarios)
    Distance_matrix = pairwise(Euclidean(), next_prices; dims=2)
    # Find the clusters 
    clusters = kmedoids(Distance_matrix, N; maxiter=200, display=:iter)
    medoids_indices = clusters.medoids
    medoids_values = zeros(number_of_warehouses,N)
         
    for i in 1:N
        medoids_values[:,i] = next_prices[:,medoids_indices[i]]
    end
        
    scenario_assignments = assignments(clusters) #Assigning which scenario belongs to which cluster 
        
    # Probabilities of the medoids
    Probs = zeros(N)
    for i in scenario_assignments
        Probs[i] = Probs[i] + 1/number_of_scenarios
    end
    
    reduced_next_prices = medoids_values
    return reduced_next_prices, Probs
end


#Performs fast forward selection for the given parameters
#D = Symmetric distance matrix
#p = vector of probabilities
#n = target number of scenarios
#Returns Array with 2 element, [1] = list of probabilities, [2] = list of selected scenario indices
function FastForwardSelection(next_prices, no_of_selected_scnearios)
    D = pairwise(Euclidean(), next_prices; dims=2)  # Create Euclidean Distance Distance_matrix
    p = fill(1/number_of_scenarios, number_of_scenarios)
    n = no_of_selected_scnearios
    init_d = copy(D)
    not_selected_scenarios = collect(range(1,length(D[:,1]);step=1))
    selected_scenarios = []
    while length(selected_scenarios) < n
        selected = select_scenario(D, p, not_selected_scenarios)
        deleteat!(not_selected_scenarios, findfirst(isequal(selected), not_selected_scenarios))
        push!(selected_scenarios, selected)
        D = UpdateDistanceMatrix(D, selected, not_selected_scenarios)
    end
    result_prob = RedistributeProbabilities(init_d, p, selected_scenarios, not_selected_scenarios)
    reduced_next_prices = Array{Float64}(undef, number_of_warehouses, n)
    for i in 1:n
        reduced_next_prices[:,i] = next_prices[:,selected_scenarios[i]]
    end
    return reduced_next_prices, result_prob
end

#Redistributes probabilities at the end of the fast forward selection
#D = original distance matrix
#p = probabilities
#selected_scenarios = indices of selected scenarios
#not_selected_scenarios = indices of non selected scenarios
function RedistributeProbabilities(D, p, selected_scenarios, not_selected_scenarios)
    probabilities = p
    for s in not_selected_scenarios
        min_idx = -1
        min_dist = Inf
        for i in selected_scenarios
            if D[s,i] < min_dist
                min_idx = i
                min_dist = D[s,i]
            end
        end
        probabilities[min_idx] = probabilities[min_idx] + p[s]
        probabilities[s] = 0.0
    end
    new_probabilities = [probabilities[i] for i in selected_scenarios]
    return new_probabilities
end

#Updates the distance matrix in the fast forward selection
#D = current distance matrix
#selected = index of scenario selected in this iteration
#scenarios = index list of not selected scenarios
function UpdateDistanceMatrix(D, selected, not_selected_scenarios)
    for s in not_selected_scenarios
        if s!=selected
            for s2 in not_selected_scenarios
                if s2!=selected
                    D[s,s2] = min(D[s,s2], D[s,selected])
                end
            end
        end
    end
    return D
end

#Selects the scenario idx with minimum Kantorovic distance
#D = Distance matrix
#p = probabilities
#scenarios = not selected scenarios
function select_scenario(D, p, not_selected_scenarios)
    min_dist = Inf
    min_idx = -1
    for s in not_selected_scenarios
        dist = sum(p[s2]*D[s2,s] for s2 in not_selected_scenarios if s!=s2)
        if dist < min_dist
            min_dist = dist
            min_idx = s
        end
    end
    return min_idx
end



# This function I am really not sure if it's doing the right thing honestly : 
# function FastForward(next_prices, no_of_selected_scnearios)
#     N = no_of_selected_scnearios
#     original = copy(next_prices)
#     subset=[]

#     Distance_matrix = pairwise(Euclidean(), next_prices; dims=2)  # Create Euclidean Distance Distance_matrix
#     Probs = fill(1/number_of_scenarios, number_of_scenarios) # Probability vector with same probability for each scenario
#     Probs_og = deepcopy(Probs)
#     Distance_matrix_og = deepcopy(Distance_matrix)
#     n_scenarios_updated = size(Distance_matrix,1)
    
#     # Selecting N scenarios with the shortest Kantorovich Distance
#     for iteration in 1:N
#         d_k = 0
#         d_k_min = 0
#         selected_s = 1
#         start = true
#         for line in 1:n_scenarios_updated
#             for column in 1:n_scenarios_updated
#                 d_k += Distance_matrix[line,column] * Probs[column]
#             end
#             if start == true
#                 start = false
#                 d_k_min = copy(d_k)
#             elseif d_k < d_k_min
#                 d_k_min = copy(d_k)
#                 selected_s = line #saves index of selected scenario
#             end
#             d_k = 0
#         end
#         # println(selected_s)
#         # Update Distance matrix 
#         Distance_matrix_old = deepcopy(Distance_matrix)
#         Distance_matrix = Array{Float64}(undef, n_scenarios_updated-1, n_scenarios_updated-1)
#         for line in 1:(selected_s-1)
#             for column in 1:(selected_s-1)
#                 Distance_matrix[line,column] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
#             end
#             for column in (selected_s+1):n_scenarios_updated
#                 Distance_matrix[line,column-1] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
#             end   
#         end 
#         for line in (selected_s+1):n_scenarios_updated
#             for column in 1:(selected_s-1)
#                 Distance_matrix[line-1,column] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
#             end
#             for column in (selected_s+1):n_scenarios_updated
#                 Distance_matrix[line-1,column-1] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
#             end
#         end                

#         # Slice Prob_matrix
#         to_keep = setdiff(1:n_scenarios_updated, selected_s)  # Get the indices of kept scenarios
#         Probs = Probs[to_keep]

#         n_scenarios_updated = n_scenarios_updated-1

#         # Scenario sets (original: scenarios not to keep, subset: scenarios to keep)
#         # println(original)
#         push!(subset,original[:,selected_s])
#         original = original[:, to_keep]
#     end

#     # Update Probability matrix
#     subset = hcat(subset...)
#     selected_indices=[] # column indices in next_prices of the selected scenarios
#     for col in eachcol(subset)
#         index_in_og = findfirst(x -> x == col, eachcol(next_prices))
#         push!(selected_indices, index_in_og)
#     end
#     removed_indices = setdiff(1:number_of_scenarios, selected_indices) # column indices in next_prices of the non-selected scenarios
    
#     start = true
#     closest_scen = selected_indices[1]
#     d_k_min=0
#     for r in removed_indices
#         for s in selected_indices
#             d_k = Distance_matrix_og[r,s]
#             if start == true
#                 start = false
#                 d_k_min = copy(d_k)
#             elseif d_k<d_k_min
#                 d_k_min = copy(d_k)
#                 closest_scen = s
#             end
#         end
#         Probs_og[closest_scen]+= Probs_og[r]
#         start = true
#         closest_scen = selected_indices[1]
#     end
#     Probs = Probs_og[selected_indices]
#     reduced_next_prices = subset
    
#     return reduced_next_prices, Probs
# end

# reduced_next_prices, Probs = FastForward(next_prices, N)
