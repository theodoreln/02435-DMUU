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


function OiH_policy(cost_coffee)

    #Declare model with Gurobi solver
    model_OiH = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

    #Declare the variables to optimize
    # Quantities of coffee ordered, W rows and T columns
    @variable(model_OiH, quantities_ordered[1:number_of_warehouses, 1:number_of_simulation_periods]>=0)
    # Quantities send from w to q, W rows W columns and T layers
    @variable(model_OiH, quantities_send[1:number_of_warehouses, 1:number_of_warehouses, 1:number_of_simulation_periods]>=0)
    # Quantities recieved by w from q, W rows W columns and T layers
    @variable(model_OiH, quantities_recieved[1:number_of_warehouses, 1:number_of_warehouses, 1:number_of_simulation_periods]>=0)
    # Quantities in the warehouse stockage, W rows and T columns
    @variable(model_OiH, quantities_stocked[1:number_of_warehouses, 1:number_of_simulation_periods]>=0)
    # Quantities mising to complete the demand, W rows and T columns
    @variable(model_OiH, quantities_missed[1:number_of_warehouses, 1:number_of_simulation_periods]>=0)
    
    #Objective function
    @objective(model_OiH, Min, sum(quantities_ordered[w,t]*cost_coffee[w,t] for w in 1:number_of_warehouses, t in 1:number_of_simulation_periods) 
    + sum(quantities_send[w,q,t]*cost_tr[w,q] for w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:number_of_simulation_periods)
    + sum(quantities_missed[w,t]*cost_miss[w] for w in 1:number_of_warehouses, t in 1:number_of_simulation_periods))

    #Constraints of the problem
    # Constraint on stockage capacities limited to the maximum capacities
    @constraint(model_OiH, Stockage_limit[w in 1:number_of_warehouses, t in 1:number_of_simulation_periods], quantities_stocked[w,t] <= warehouse_capacities[w])
    # Constraint on transport capacities limited to the maximum capacities
    @constraint(model_OiH, Transport_limit[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:number_of_simulation_periods], quantities_send[w,q,t] <= transport_capacities[w,q])
    # Constraint on quantity send equal quantity recieved
    @constraint(model_OiH, Send_recieved[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 1:number_of_simulation_periods], quantities_send[w,q,t] == quantities_recieved[q,w,t])
    # Constraint on a warehouse can only send to others warehouse
    # Useless cause the self-transport capacity is equal to 0
    @constraint(model_OiH, Self_transport[w in 1:number_of_warehouses, t in 1:number_of_simulation_periods], quantities_send[w,w,t] == 0)
    # Constraint on quantity send limited to previous stock
    @constraint(model_OiH, Transport_stock[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 2:number_of_simulation_periods], sum(quantities_send[w,q,t] for q in 1:number_of_warehouses) <= quantities_stocked[w,t-1])
    @constraint(model_OiH, Transport_stock_start[w in 1:number_of_warehouses, q in 1:number_of_warehouses], sum(quantities_send[w,q,1] for q in 1:number_of_warehouses) <= initial_stock[w])
    # Constraint on quantity stock at time t with input and output
    @constraint(model_OiH, Stockage[w in 1:number_of_warehouses, t in 2:number_of_simulation_periods], quantities_stocked[w,t] == quantities_stocked[w,t-1]+quantities_ordered[w,t]
    +sum(quantities_recieved[w,q,t] - quantities_send[w,q,t] for q in 1:number_of_warehouses)- demand_trajectory[w,t] + quantities_missed[w,t])
    @constraint(model_OiH, Stockage_start[w in 1:number_of_warehouses], quantities_stocked[w,1] == initial_stock[w]+quantities_ordered[w,1]
    +sum(quantities_recieved[w,q,1] - quantities_send[w,q,1] for q in 1:number_of_warehouses)- demand_trajectory[w,1] + quantities_missed[w,1])

    optimize!(model_OiH)

    #Check if optimal solution was found
    if termination_status(model_OiH) == MOI.OPTIMAL
        println("Optimal solution found")

        # Return interesting values
        return objective_value(model_OiH)
    else
        return error("No solution.")
    end
end