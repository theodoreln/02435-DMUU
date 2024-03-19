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

