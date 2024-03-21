# Policy for the multi-stage stochastic optimization

# Useful files
include("price_process.jl")
# include("ScenarioSelections.jl")

# Loading the problem's parameters
include("multistage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_sim_periods, sim_T, demand_trajectory = load_the_data()

#Import packages
using Random
using JuMP
using Gurobi
using Printf
using Clustering
using Distances

function make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices, look_ahead_days, nb_initial_scenarios, granularity, nb_reduced_scenarios, type_red)
    
    # Control the number of ahead day to not overpass the limit 
    # The variable "actual_look_ahead_days" will decide how much day ahead we are going to look
    # Its value can vary between "look_ahead_days" and 1 !!!
    if look_ahead_days > number_of_sim_periods - tau
        actual_look_ahead_days = number_of_sim_periods-tau+1
    else 
        actual_look_ahead_days = look_ahead_days+1
    end
    # Here the current day is included

    # Generate the initial scenarios
    prices_trajectory_scenarios = Generate_scenarios(number_of_warehouses, W, current_prices, nb_initial_scenarios, actual_look_ahead_days)

    # Discretize prices
    prices_trajectory_scenarios .= round.(prices_trajectory_scenarios ./ granularity) .* granularity

    #Reduce the number of scenarios
    if type_red == "kmeans" 
        prices_trajectory_reduced, Probs = kmeans_selection(prices_trajectory_scenarios, nb_initial_scenarios, nb_reduced_scenarios, granularity)
    elseif type_red == "kmedoids"
        prices_trajectory_reduced, Probs = kmedoids_selection(number_of_warehouses, prices_trajectory_scenarios, nb_initial_scenarios, nb_reduced_scenarios, actual_look_ahead_days)
    else 
        prices_trajectory_reduced, Probs = FastForwardSelection(number_of_warehouses, prices_trajectory_scenarios, nb_initial_scenarios, nb_reduced_scenarios, actual_look_ahead_days)
    end

    #Populate
    Sets = Populating_sets(actual_look_ahead_days, prices_trajectory_reduced, nb_reduced_scenarios)

    ######################
    # Optimization problem
    # Declare model with Gurobi solver
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

    #Declare the variables to optimize
    # Quantities of coffee ordered
    @variable(model, quantities_ordered[1:number_of_warehouses, 1:actual_look_ahead_days, 1:nb_reduced_scenarios]>=0)
    # Quantities send from w to q
    @variable(model, quantities_send[1:number_of_warehouses, 1:number_of_warehouses, 1:actual_look_ahead_days, 1:nb_reduced_scenarios]>=0)
    # Quantities recieved by w from q
    @variable(model, quantities_recieved[1:number_of_warehouses, 1:number_of_warehouses, 1:actual_look_ahead_days, 1:nb_reduced_scenarios]>=0)
    # Quantities in the warehouse stockage
    @variable(model, quantities_stocked[1:number_of_warehouses, 1:actual_look_ahead_days, 1:nb_reduced_scenarios]>=0)
    # Quantities mising to complete the demand
    @variable(model, quantities_missed[1:number_of_warehouses, 1:actual_look_ahead_days, 1:nb_reduced_scenarios]>=0)

    #Objective function
    @objective(model, Min, sum(Probs[s]*(sum(quantities_ordered[w,t,s]*prices_trajectory_reduced[w,t,s] for w in 1:number_of_warehouses, t in 1:actual_look_ahead_days) 
    + sum(quantities_send[w,q,t,s]*cost_tr[w,q] for w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:actual_look_ahead_days)
    + sum(quantities_missed[w,t,s]*cost_miss[w] for w in 1:number_of_warehouses, t in 1:actual_look_ahead_days)) for s in 1:nb_reduced_scenarios))

    #Constraints of the problem
    # Constraint on stockage capacities limited to the maximum capacities
    @constraint(model, Stockage_limit[w in 1:number_of_warehouses, t in 1:actual_look_ahead_days, s in 1:nb_reduced_scenarios], quantities_stocked[w,t,s] <= warehouse_capacities[w])
    # Constraint on transport capacities limited to the maximum capacities
    @constraint(model, Transport_limit[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:actual_look_ahead_days, s in 1:nb_reduced_scenarios], quantities_send[w,q,t,s] <= transport_capacities[w,q])
    # Constraint on quantity send equal quantity recieved
    @constraint(model, Send_recieved[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:actual_look_ahead_days, s in 1:nb_reduced_scenarios], quantities_send[w,q,t,s] == quantities_recieved[q,w,t,s])
    # Constraint on a warehouse can only send to others warehouse
    @constraint(model, Self_transport[w in 1:number_of_warehouses, t in 1:actual_look_ahead_days, s in 1:nb_reduced_scenarios], quantities_send[w,w,t,s] == 0)
    # Constraint on quantity send limited to previous stock
    @constraint(model, Transport_stock[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 2:actual_look_ahead_days, s in 1:nb_reduced_scenarios], sum(quantities_send[w,q,t,s] for q in 1:number_of_warehouses) <= quantities_stocked[w,t-1,s])
    @constraint(model, Transport_stock_start[w in 1:number_of_warehouses, q in 1:number_of_warehouses, s in 1:nb_reduced_scenarios], sum(quantities_send[w,q,1,s] for q in 1:number_of_warehouses) <= current_stock[w])
    # Constraint on quantity stock at time t with input and output
    @constraint(model, Stockage[w in 1:number_of_warehouses, t in 2:actual_look_ahead_days, s in 1:nb_reduced_scenarios], quantities_stocked[w,t,s] == quantities_stocked[w,t-1,s]+quantities_ordered[w,t,s]
    +sum(quantities_recieved[w,q,t,s] - quantities_send[w,q,t,s] for q in 1:number_of_warehouses)- demand_trajectory[w,t] + quantities_missed[w,t,s])
    @constraint(model, Stockage_start[w in 1:number_of_warehouses, s in nb_reduced_scenarios], quantities_stocked[w,1,s] == current_stock[w]+quantities_ordered[w,1,s]
    +sum(quantities_recieved[w,q,1,s] - quantities_send[w,q,1,s] for q in 1:number_of_warehouses)- demand_trajectory[w,1] + quantities_missed[w,1,s])
    # Constraints of non-anticipativity
    Keys_list = collect(keys(Sets))
    for key in Keys_list
        set = key[1]
        time = key[2]
        Others_sets = Sets[key]
        for setp in Others_sets 
            @constraint(model, [w in 1:number_of_warehouses], quantities_ordered[w,time,set] == quantities_ordered[w,time,setp])
            @constraint(model, [w in 1:number_of_warehouses, q in 1:number_of_warehouses], quantities_send[w,q,time,set] == quantities_send[w,q,time,setp])
            @constraint(model, [w in 1:number_of_warehouses, q in 1:number_of_warehouses], quantities_recieved[w,q,time,set] == quantities_recieved[w,q,time,setp])
            @constraint(model, [w in 1:number_of_warehouses], quantities_stocked[w,time,set] == quantities_stocked[w,time,setp])
            @constraint(model, [w in 1:number_of_warehouses], quantities_missed[w,time,set] == quantities_missed[w,time,setp])
        end
    end
    
    # Optimize the model
    optimize!(model)

    #Check if optimal solution was found
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found")
        
        #return interesting values
        return value.(quantities_ordered)[:,1,1],value.(quantities_send)[:,:,1,1],value.(quantities_recieved)[:,:,1,1],value.(quantities_stocked)[:,1,1],value.(quantities_missed)[:,1,1]
    else
        return error("No solution.")
    end
end

# Generate the number of scenarios we want
function Generate_scenarios(number_of_warehouses, W, current_prices, nb_initial_scenarios, actual_look_ahead_days)

    Scen = collect(1:nb_initial_scenarios)

    # Output vector with all scenarios shape : (number_of_warehouses[3], actual_look_ahead_days, nb_initial_scenarios)
    prices_trajectory_scenarios = zeros(number_of_warehouses, actual_look_ahead_days, nb_initial_scenarios)
    for s in Scen
        for w in W
            prices_trajectory_scenarios[w,1,s] = current_prices[w]
            for t in 2:actual_look_ahead_days
                prices_trajectory_scenarios[w,t,s] = sample_next(prices_trajectory_scenarios[w,t-1,s])
            end
        end
    end
    return prices_trajectory_scenarios
end

# Kmeans selections
function kmeans_selection(prices_trajectory_scenarios, nb_initial_scenarios, nb_reduced_scenarios, granularity)

    reshaped_array = reshape(prices_trajectory_scenarios, :, size(prices_trajectory_scenarios, 3))
    clusters = kmeans(reshaped_array, nb_reduced_scenarios; maxiter=2000)
    clustered_prices = clusters.centers # get the cluster centers
        
    # Probabilities of data point belonging to cluster 
    scenario_assignments = assignments(clusters) #Assigning which scenario belongs to which cluster 
    
    Probs = zeros(nb_reduced_scenarios)
    for i in scenario_assignments
        Probs[i] = Probs[i] + 1/nb_initial_scenarios
    end

    prices_trajectory_reduced = reshape(clustered_prices, size(prices_trajectory_scenarios)[1], size(prices_trajectory_scenarios)[2], :)
    prices_trajectory_reduced .= round.(prices_trajectory_reduced ./ granularity) .* granularity
    
    return prices_trajectory_reduced, Probs
end

# Kmedoids selections
function kmedoids_selection(number_of_warehouses, prices_trajectory_scenarios, nb_initial_scenarios, nb_reduced_scenarios, actual_look_ahead_days)

    if actual_look_ahead_days != 1
        reshaped_array = reshape(prices_trajectory_scenarios, :, size(prices_trajectory_scenarios, 3))
        # Calculate Euclidean Distance matrix (size: number_of_scenarios x number_of_scenarios)
        Distance_matrix = pairwise(Euclidean(), reshaped_array; dims=2)
        # Find the clusters 
        clusters = kmedoids(Distance_matrix, nb_reduced_scenarios; maxiter=200)
        medoids_indices = clusters.medoids
        medoids_values = zeros(number_of_warehouses*actual_look_ahead_days,nb_reduced_scenarios)
            
        for i in 1:nb_reduced_scenarios
            medoids_values[:,i] = reshaped_array[:,medoids_indices[i]]
        end
            
        scenario_assignments = assignments(clusters) #Assigning which scenario belongs to which cluster 
            
        # Probabilities of the medoids
        Probs = zeros(nb_reduced_scenarios)
        for i in scenario_assignments
            Probs[i] = Probs[i] + 1/nb_initial_scenarios
        end
        
        prices_trajectory_reduced = reshape(medoids_values, size(prices_trajectory_scenarios)[1], size(prices_trajectory_scenarios)[2], :)
    else
        prices_trajectory_reduced = prices_trajectory_scenarios[:,:,1:nb_reduced_scenarios]
        Probs = zeros(nb_reduced_scenarios) 
        for i in 1:nb_reduced_scenarios
            Probs[i] = 1/nb_reduced_scenarios
        end
    end

    return prices_trajectory_reduced, Probs
end


#Performs fast forward selection for the given parameters
#D = Symmetric distance matrix
#p = vector of probabilities
#n = target number of scenarios
#Returns Array with 2 element, [1] = list of probabilities, [2] = list of selected scenario indices
function FastForwardSelection(number_of_warehouses, prices_trajectory_scenarios, nb_initial_scenarios, nb_reduced_scenarios, actual_look_ahead_days)
    
    if actual_look_ahead_days != 1
        reshaped_array = reshape(prices_trajectory_scenarios, :, size(prices_trajectory_scenarios, 3))
        # Calculate Euclidean Distance matrix (size: number_of_scenarios x number_of_scenarios)
        Distance_matrix = pairwise(Euclidean(), reshaped_array; dims=2)
        p = fill(1/nb_initial_scenarios, nb_initial_scenarios)
        init_d = copy(Distance_matrix)
        n= nb_reduced_scenarios
        not_selected_scenarios = collect(range(1,length(Distance_matrix[:,1]);step=1))
        selected_scenarios = []
        while length(selected_scenarios) < n
            selected = select_scenario(Distance_matrix, p, not_selected_scenarios)
            deleteat!(not_selected_scenarios, findfirst(isequal(selected), not_selected_scenarios))
            push!(selected_scenarios, selected)
            Distance_matrix = UpdateDistanceMatrix(Distance_matrix, selected, not_selected_scenarios)
        end
        result_prob = RedistributeProbabilities(init_d, p, selected_scenarios, not_selected_scenarios)
        reduced_next_prices = Array{Float64}(undef, number_of_warehouses*actual_look_ahead_days, n)
        for i in 1:n
            reduced_next_prices[:,i] = reshaped_array[:,selected_scenarios[i]]
        end
        prices_trajectory_reduced = reshape(reduced_next_prices, size(prices_trajectory_scenarios)[1], size(prices_trajectory_scenarios)[2], :)
    else 
        prices_trajectory_reduced = prices_trajectory_scenarios[:,:,1:nb_reduced_scenarios]
        result_prob = zeros(nb_reduced_scenarios) 
        for i in 1:nb_reduced_scenarios
            result_prob[i] = 1/nb_reduced_scenarios
        end
    end 
    return prices_trajectory_reduced, result_prob
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


# Populating the non-anticipativity sets
function Populating_sets(actual_look_ahead_days, prices_trajectory_reduced, nb_reduced_scenarios)
    
    Sets = Dict()
    Scen = collect(1:nb_reduced_scenarios)
    Day = collect(1:actual_look_ahead_days)

    # Check every pair of scenarios and look at which point they are similar
    for s in Scen 
        # To save some memory and time we only look in one direction, it will reduce the number of constraints also
        for sp in s+1:nb_reduced_scenarios
            sim = 1
            for d in Day 
                if sim == 1
                    if prices_trajectory_reduced[:,d,s] == prices_trajectory_reduced[:,d,sp]
                        key = (s,d)
                        if haskey(Sets, key)
                            push!(Sets[key], sp)
                        else 
                            Sets[key] = [sp]
                        end
                    else 
                        sim = 0
                    end
                end
            end
        end
    end
    return(Sets)
end


##########
# To try #
##########

include("simulation_experiments.jl")
# Creating the random experiments on which the policy will be evaluated
number_of_experiments, Expers, Price_experiments = simulation_experiments_creation(number_of_warehouses, W, number_of_sim_periods)
# Price_experiments shape (number_of_experiments[40], number of warehouses[3], number of periods [5])

# current_prices = Price_experiments[1,:,1]
# nb_initial_scenarios = 1000
# nb_reduced_scenarios = 50
# actual_look_ahead_days = 1
# granularity = 0.5
# prices_trajectory = Generate_scenarios(number_of_warehouses, W, current_prices, nb_initial_scenarios, actual_look_ahead_days)
# reshaped_array = reshape(prices_trajectory, :, size(prices_trajectory, 3))
# Distance_matrix = pairwise(Euclidean(), reshaped_array; dims=2)
# clusters = kmedoids(Distance_matrix, nb_reduced_scenarios; maxiter=200, display=:iter)
# medoids_indices = clusters.medoids
# medoids_values = zeros(number_of_warehouses*actual_look_ahead_days,nb_reduced_scenarios)
# for i in 1:nb_reduced_scenarios
#     medoids_values[:,i] = reshaped_array[:,medoids_indices[i]]
# end
# scenario_assignments = assignments(clusters) #Assigning which scenario belongs to which cluster 
# # Probabilities of the medoids
# Probs = zeros(nb_reduced_scenarios)
# for i in scenario_assignments
#     Probs[i] = Probs[i] + 1/nb_initial_scenarios
# end
# prices_trajectory_reduced = reshape(medoids_values, size(prices_trajectory)[1], size(prices_trajectory)[2], :)

# prices_trajectory .= round.(prices_trajectory ./ granularity) .* granularity
# reshaped_data = reshape(prices_trajectory, :, size(prices_trajectory, 3))
# clusters = kmeans(reshaped_data, nb_reduced_scenarios; maxiter=2000)
# clustered_prices = clusters.centers 
# prices_trajectory_reduced = reshape(clustered_prices, size(prices_trajectory)[1], size(prices_trajectory)[2], :)
# prices_trajectory_reduced .= round.(prices_trajectory_reduced ./ granularity) .* granularity
# Sets = Populating_sets(actual_look_ahead_days, prices_trajectory_reduced, nb_reduced_scenarios)
# Keys_list = collect(keys(Sets))

# Testing FastForward
# current_prices = Price_experiments[1,:,1]
# nb_initial_scenarios = 100
# nb_reduced_scenarios = 20
# granularity = 0.5
# actual_look_ahead_days = 3
# prices_trajectory= Generate_scenarios(number_of_warehouses, W, current_prices, nb_initial_scenarios, actual_look_ahead_days)
# prices_trajectory .= round.(prices_trajectory ./ granularity) .* granularity
# reshaped_data = reshape(prices_trajectory, :, size(prices_trajectory, 3))

# print(kmeans_selection(prices_trajectory, nb_initial_scenarios, nb_reduced_scenarios, granularity))
# print(FastForwardSelection(number_of_warehouses, prices_trajectory, nb_initial_scenarios, nb_reduced_scenarios, actual_look_ahead_days))