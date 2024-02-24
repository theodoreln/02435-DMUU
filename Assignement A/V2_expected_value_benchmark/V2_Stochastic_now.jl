#Import useful files
include("V2_02435_two_stage_problem_data.jl")
include("V2_price_process.jl")

#Import packages
using Random
using JuMP
using Gurobi
using Printf
using Clustering
using Distances
using Plots

prices=round.(10 * rand(3), digits=2)

function Make_Stochastic_here_and_now_decision(prices, method)
# Inputs of function are the prices for today and the method to use for scenario selection. Those are 
#"K-means", "K-medoids" and "FastForward".
    
    #Importation of the inputs from the two stage problem
    number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data()

    #Number of scenarios
    number_of_scenarios = 100
    #number_list = collect(1:number_of_scenarios)
        
    #Computation of the prices of different scenarios at t=2
    next_prices = Array{Float64}(undef, number_of_warehouses, number_of_scenarios)
    for w in 1:number_of_warehouses
        for n in 1:number_of_scenarios
            next_prices[w,n] = sample_next(prices[w])
        end
    end    



    #Sceario Selection

    N = 10 # Reduced number of scenarios

    # K-means
    if method == "K-means" 
        clusters = kmeans(next_prices, N; maxiter=200, display=:iter)
        clustered_prices = clusters.centers # get the cluster centers
        
         # Probabilities of data point belonging to cluster 
        scenario_assignments = assignments(clusters) #Assigning which scenario belongs to which cluster 
        
        Probs = zeros(N)
        for i in scenario_assignments
            Probs[i] = Probs[i] + 1/number_of_scenarios
        end

        next_prices= clustered_prices # Update the matrix of scenarios to matrix with reduced scenarios
    end

    # K-medoids
    if method == "K-medoids"
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

        next_prices = medoids_values


#= 
        #Prepare medoid data for plotting
        medoid_data = zeros(Float64, 3, number_of_scenarios)
        for i = 1:N
            medoid_data[:,i] = next_prices[:,medoids_indices[i]]
        end
        
        #Plot medoids on sampled_scenarios
        plotdata = hcat(range(1,4;step=1), medoid_data)
        plot!(plotdata[:,1],plotdata[:,2:N+1], legend=false, color=:auto)
        savefig(figure,"kmedoids-result.pdf")
         =#

    end

    if method =="FastForward"

        original = next_prices
        subset=[]

        Distance_matrix = pairwise(Euclidean(), next_prices; dims=2)  # Create Euclidean Distance Distance_matrix
        Probs = fill(1/number_of_scenarios, number_of_scenarios) # Probability vector with same probability for each scenario
        Probs_og = Probs
        Distance_matrix_og = Distance_matrix
        n_scenarios_updated = number_of_scenarios
        
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
        next_prices = subset
    end



    #Creation of the parameters of the problem
    # Creation of the demand, W rows and T columns, 10 everywhere
    demand_coffee = demand_trajectory

    number_of_scenarios = N
    #Declare model with Gurobi solver
    model_ST = Model(Gurobi.Optimizer)

    #Declare the variables to optimize
    # Quantities of coffee ordered, W rows
    @variable(model_ST, quantities_ordered_now[1:number_of_warehouses]>=0)
    @variable(model_ST, quantities_ordered_sec[1:number_of_warehouses, 1:number_of_scenarios]>=0)
    # Quantities send from w to q, W rows W columns
    @variable(model_ST, quantities_send_now[1:number_of_warehouses, 1:number_of_warehouses]>=0)
    @variable(model_ST, quantities_send_sec[1:number_of_warehouses, 1:number_of_warehouses, 1:number_of_scenarios]>=0)
    # Quantities recieved by w from q, W rows W columns
    @variable(model_ST, quantities_recieved_now[1:number_of_warehouses, 1:number_of_warehouses]>=0)
    @variable(model_ST, quantities_recieved_sec[1:number_of_warehouses, 1:number_of_warehouses, 1:number_of_scenarios]>=0)
    # Quantities in the warehouse stockage, W rows
    @variable(model_ST, quantities_stocked_now[1:number_of_warehouses]>=0)
    @variable(model_ST, quantities_stocked_sec[1:number_of_warehouses, 1:number_of_scenarios]>=0)
    # Quantities mising to complete the demand, W rows
    @variable(model_ST, quantities_missed_now[1:number_of_warehouses]>=0)
    @variable(model_ST, quantities_missed_sec[1:number_of_warehouses, 1:number_of_scenarios]>=0)
    
    #Objective function
    @objective(model_ST, Min, sum(quantities_ordered_now[w]*prices[w] for w in 1:number_of_warehouses) + sum(quantities_send_now[w,q]*cost_tr[w,q] for w in 1:number_of_warehouses, q in 1:number_of_warehouses)
    + sum(quantities_missed_now[w]*cost_miss[w] for w in 1:number_of_warehouses) + sum(1/number_of_scenarios*(sum(quantities_ordered_sec[w,n]*next_prices[w,n] for w in 1:number_of_warehouses)
    + sum(quantities_send_sec[w,q,n]*cost_tr[w,q] for w in 1:number_of_warehouses, q in 1:number_of_warehouses)
    + sum(quantities_missed_sec[w,n]*cost_miss[w] for w in 1:number_of_warehouses)) for n in 1:number_of_scenarios))

    #Constraints of the problem
    # Constraint on stockage capacities limited to the maximum capacities
    @constraint(model_ST, Stockage_limit_now[w in 1:number_of_warehouses], quantities_stocked_now[w] <= warehouse_capacities[w])
    @constraint(model_ST, Stockage_limit_sec[w in 1:number_of_warehouses, n in 1:number_of_scenarios], quantities_stocked_sec[w,n] <= warehouse_capacities[w])
    # Constraint on transport capacities limited to the maximum capacities
    @constraint(model_ST, Transport_limit_now[w in 1:number_of_warehouses, q in 1:number_of_warehouses], quantities_send_now[w,q] <= transport_capacities[w,q])
    @constraint(model_ST, Transport_limit_sec[w in 1:number_of_warehouses, q in 1:number_of_warehouses, n in 1:number_of_scenarios], quantities_send_sec[w,q,n] <= transport_capacities[w,q])
    # Constraint on quantity send equal quantity recieved
    @constraint(model_ST, Send_recieved_now[w in 1:number_of_warehouses, q in 1:number_of_warehouses], quantities_send_now[w,q] == quantities_recieved_now[q,w])
    @constraint(model_ST, Send_recieved_sec[w in 1:number_of_warehouses, q in 1:number_of_warehouses, n in 1:number_of_scenarios], quantities_send_sec[w,q,n] == quantities_recieved_sec[q,w,n])
    # Constraint on a warehouse can only send to others warehouse
    # Useless cause the self-transport capacity is equal to 0
    @constraint(model_ST, Self_transport_now[w in 1:number_of_warehouses], quantities_send_now[w,w] == 0)
    @constraint(model_ST, Self_transport_sec[w in 1:number_of_warehouses, n in 1:number_of_scenarios], quantities_send_sec[w,w,n] == 0)
    # Constraint on quantity send limited to previous stock
    @constraint(model_ST, Transport_stock_now[w in 1:number_of_warehouses, q in 1:number_of_warehouses], sum(quantities_send_now[w,q] for q in 1:number_of_warehouses) <= initial_stock[w])
    @constraint(model_ST, Transport_stock_sec[w in 1:number_of_warehouses, q in 1:number_of_warehouses, n in 1:number_of_scenarios], sum(quantities_send_sec[w,q,n] for q in 1:number_of_warehouses) <= quantities_stocked_now[w])
    # Constraint on quantity stock at time t with input and output
    @constraint(model_ST, Stockage_now[w in 1:number_of_warehouses], quantities_stocked_now[w] == initial_stock[w]+quantities_ordered_now[w]
    +sum(quantities_recieved_now[w,q] - quantities_send_now[w,q] for q in 1:number_of_warehouses)- demand_coffee[w,1] + quantities_missed_now[w])
    @constraint(model_ST, Stockage_sec[w in 1:number_of_warehouses, n in 1:number_of_scenarios], quantities_stocked_sec[w,n] == quantities_stocked_now[w]+quantities_ordered_sec[w,n]
    +sum(quantities_recieved_sec[w,q,n] - quantities_send_sec[w,q,n] for q in 1:number_of_warehouses)- demand_coffee[w,2] + quantities_missed_sec[w,n])


    optimize!(model_ST)


    #Check if optimal solution was found
    if termination_status(model_ST) == MOI.OPTIMAL
        println("Optimal solution found")

        #return interesting values
        return value.(quantities_ordered_now),value.(quantities_send_now),value.(quantities_recieved_now),value.(quantities_stocked_now),value.(quantities_missed_now),objective_value(model_ST)
    else
        return error("No solution.")
    end

end

#qo_ST,qs_ST,qr_ST,qst_ST,qm_ST,ov_ST=Make_Stochastic_here_and_now_decision(prices,"K-means")
#qo_ST2,qs_ST2,qr_ST2,qst_ST2,qm_ST2,ov_ST2=Make_Stochastic_here_and_now_decision(prices,"K-medoids")
#qo_ST3,qs_ST3,qr_ST3,qst_ST3,qm_ST3,ov_ST3=Make_Stochastic_here_and_now_decision(prices,"FastForward")