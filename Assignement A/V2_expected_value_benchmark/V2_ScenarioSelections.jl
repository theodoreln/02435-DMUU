using Random
using JuMP
using Gurobi
using Printf
using Clustering
using Distances
using Plots

include("V2_02435_two_stage_problem_data.jl")
include("V2_price_process.jl")


# All of the scenario selection functions take the matrix of all scenarios (next_prices -> 3xN) 
# and the number the scenarios shall be reduced to (k) as arguments and returns the 
# reduced scenario matrix (3xk) and the sceanrios' new probabilities (Probs -> 1xk)
# prices=round.(10 * rand(3), digits=2)
number_of_warehouses = 3
number_of_scenarios = 1000
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


#print("kmeans: ",kmeans_selection(next_prices, N))

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
    
#print("kmedoids: ", kmedoids_selection(next_prices, N))

function FastForward(next_prices, no_of_selected_scnearios)
    N = no_of_selected_scnearios
    original = next_prices
    subset=[]

    Distance_matrix = pairwise(Euclidean(), next_prices; dims=2)  # Create Euclidean Distance Distance_matrix
    Probs = fill(1/number_of_scenarios, number_of_scenarios) # Probability vector with same probability for each scenario
    Probs_og = Probs
    Distance_matrix_og = Distance_matrix
    n_scenarios_updated = size(Distance_matrix,1)
    
    # Selecting N scenarios with the shortest Kantorovich Distance
    for iteration in 1:N
        d_k = 0
        d_k_min = 0
        selected_s = 1
        start = true
        for line in 1:n_scenarios_updated
            for column in 1:n_scenarios_updated
                d_k += Distance_matrix[line,column] * Probs[column]
            end
            if start == true
                start = false
                d_k_min = d_k
            elseif d_k < d_k_min
                d_k_min = d_k
                selected_s = line #saves index of selected scenario
            end
            d_k = 0
        end
        # Update Distance matrix 
        Distance_matrix_old = Distance_matrix
        Distance_matrix = Array{Float64}(undef, n_scenarios_updated-1, n_scenarios_updated-1)
        for line in 1:(selected_s-1)
            for column in 1:(selected_s-1)
                Distance_matrix[line,column] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
            end
            for column in (selected_s+1):n_scenarios_updated
                Distance_matrix[line,column-1] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
            end   
        end 
        for line in (selected_s+1):n_scenarios_updated
            for column in 1:(selected_s-1)
                Distance_matrix[line-1,column] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
            end
            for column in (selected_s+1):n_scenarios_updated
                Distance_matrix[line-1,column-1] = min(Distance_matrix_old[line,column], Distance_matrix_old[line,selected_s])
            end
        end                

        # Slice Prob_matrix
        to_keep = setdiff(1:n_scenarios_updated, selected_s)  # Get the indices of kept scenarios
        Probs = Probs[to_keep]

        n_scenarios_updated = n_scenarios_updated-1

        # Scenario sets (original: scenarios not to keep, subset: scenarios to keep)
        push!(subset,original[:,selected_s])
        original = original[:, to_keep]
    end

    # Update Probability matrix
    subset = hcat(subset...)
    selected_indices=[] # column indices in next_prices of the selected scenarios
    for col in eachcol(subset)
        index_in_og = findfirst(x -> x == col, eachcol(next_prices))
        push!(selected_indices, index_in_og)
    end
    removed_indices = setdiff(1:number_of_scenarios, selected_indices) # column indices in next_prices of the non-selected scenarios
    
    start = true
    closest_scen = selected_indices[1]
    d_k_min=0
    for r in removed_indices
        for s in selected_indices
            d_k = Distance_matrix_og[r,s]
            if start == true
                start = false
                d_k_min = d_k
            elseif d_k<d_k_min
                d_k_min = d_k
                closest_scen = s
            end
        end
        Probs_og[closest_scen]+= Probs_og[r]
        start = true
        closest_scen = selected_indices[1]
    end
    Probs = Probs_og[selected_indices]
    reduced_next_prices = subset
    
    return reduced_next_prices, Probs
end
