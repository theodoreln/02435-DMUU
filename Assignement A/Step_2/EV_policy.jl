# Useful files
include("price_process.jl")

# Loading the problem's parameters
include("multistage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data()

#Import packages
using Random
using JuMP
using Gurobi
using Printf


function Make_EV_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices, look_ahead_days, nb_initial_scenarios)

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

    # Mean on those scenarios to have only one path of prices
    cost_coffee = reshape(mean(prices_trajectory_scenarios, dims=3), size(prices_trajectory_scenarios, 1), size(prices_trajectory_scenarios, 2))

    #Declare model with Gurobi solver
    model_EB = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

    #Declare the variables to optimize
    # Quantities of coffee ordered, W rows and T columns
    @variable(model_EB, quantities_ordered[1:number_of_warehouses, 1:actual_look_ahead_days]>=0)
    # Quantities send from w to q, W rows W columns and T layers
    @variable(model_EB, quantities_send[1:number_of_warehouses, 1:number_of_warehouses, 1:actual_look_ahead_days]>=0)
    # Quantities recieved by w from q, W rows W columns and T layers
    @variable(model_EB, quantities_recieved[1:number_of_warehouses, 1:number_of_warehouses, 1:actual_look_ahead_days]>=0)
    # Quantities in the warehouse stockage, W rows and T columns
    @variable(model_EB, quantities_stocked[1:number_of_warehouses, 1:actual_look_ahead_days]>=0)
    # Quantities mising to complete the demand, W rows and T columns
    @variable(model_EB, quantities_missed[1:number_of_warehouses, 1:actual_look_ahead_days]>=0)
    
    #Objective function
    @objective(model_EB, Min, sum(quantities_ordered[w,t]*cost_coffee[w,t] for w in 1:number_of_warehouses, t in 1:actual_look_ahead_days) 
    + sum(quantities_send[w,q,t]*cost_tr[w,q] for w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:actual_look_ahead_days)
    + sum(quantities_missed[w,t]*cost_miss[w] for w in 1:number_of_warehouses, t in 1:actual_look_ahead_days))

    #Constraints of the problem
    # Constraint on stockage capacities limited to the maximum capacities
    @constraint(model_EB, Stockage_limit[w in 1:number_of_warehouses, t in 1:actual_look_ahead_days], quantities_stocked[w,t] <= warehouse_capacities[w])
    # Constraint on transport capacities limited to the maximum capacities
    @constraint(model_EB, Transport_limit[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:actual_look_ahead_days], quantities_send[w,q,t] <= transport_capacities[w,q])
    # Constraint on quantity send equal quantity recieved
    @constraint(model_EB, Send_recieved[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:actual_look_ahead_days], quantities_send[w,q,t] == quantities_recieved[q,w,t])
    # Constraint on a warehouse can only send to others warehouse
    # Useless cause the self-transport capacity is equal to 0
    @constraint(model_EB, Self_transport[w in 1:number_of_warehouses, t in 1:actual_look_ahead_days], quantities_send[w,w,t] == 0)
    # Constraint on quantity send limited to previous stock
    @constraint(model_EB, Transport_stock[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 2:actual_look_ahead_days], sum(quantities_send[w,q,t] for q in 1:number_of_warehouses) <= quantities_stocked[w,t-1])
    @constraint(model_EB, Transport_stock_start[w in 1:number_of_warehouses, q in 1:number_of_warehouses], sum(quantities_send[w,q,1] for q in 1:number_of_warehouses) <= current_stock[w])
    # Constraint on quantity stock at time t with input and output
    @constraint(model_EB, Stockage[w in 1:number_of_warehouses, t in 2:actual_look_ahead_days], quantities_stocked[w,t] == quantities_stocked[w,t-1]+quantities_ordered[w,t]
    +sum(quantities_recieved[w,q,t] - quantities_send[w,q,t] for q in 1:number_of_warehouses)- demand_trajectory[w,t] + quantities_missed[w,t])
    @constraint(model_EB, Stockage_start[w in 1:number_of_warehouses], quantities_stocked[w,1] == current_stock[w]+quantities_ordered[w,1]
    +sum(quantities_recieved[w,q,1] - quantities_send[w,q,1] for q in 1:number_of_warehouses)- demand_trajectory[w,1] + quantities_missed[w,1])

    optimize!(model_EB)

    #Check if optimal solution was found
    if termination_status(model_EB) == MOI.OPTIMAL
        println("Optimal solution found")

        # Return interesting values
        return value.(quantities_ordered[:,1]),value.(quantities_send[:,:,1]),value.(quantities_recieved[:,:,1]),value.(quantities_stocked[:,1]),value.(quantities_missed[:,1])
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
