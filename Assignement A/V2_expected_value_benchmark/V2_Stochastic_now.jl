#Import useful files
include("V2_02435_two_stage_problem_data.jl")
include("V2_price_process.jl")

#Import packages
using Random
using JuMP
using Gurobi
using Printf


prices=round.(10 * rand(3), digits=2)

function Make_Stochastic_here_and_now_decision(prices).
    
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


    #Creation of the parameters of the problem
    # Creation of the demand, W rows and T columns, 10 everywhere
    demand_coffee = demand_trajectory

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