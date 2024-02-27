#Import useful files
include("V2_02435_two_stage_problem_data.jl")
include("V2_price_process.jl")

#Import packages
using Random
using JuMP
using Gurobi
using Printf

# prices=round.(10 * rand(3), digits=2)

# new_prices = round.(map(Float64,map(sample_next, prices)),digits=2)

function Calculate_OiH_solution(prices, new_prices)
    
    # Importation of the inputs from the two stage problem
    number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, 
    initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data(2)

    # Transpose matrix to have columns = t and rows = w
    #prices = transpose(hcat(prices...))

    #Creation of the parameters of the problem
    # Coffee prices
    cost_coffee = [prices,new_prices]

    # Creation of the demand, W rows and T columns, 10 everywhere
    demand_coffee = demand_trajectory

    #Declare model with Gurobi solver
    model_OiH = Model(Gurobi.Optimizer)

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
    @objective(model_OiH, Min, sum(quantities_ordered[w,t]*cost_coffee[t][w] for w in 1:number_of_warehouses, t in 1:number_of_simulation_periods) 
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
    @constraint(model_OiH, Self_transport[w in 1:number_of_warehouses, t in 1:number_of_simulation_periods], quantities_send[w,w,t] == 0)
    # Constraint on quantity send limited to previous stock
    @constraint(model_OiH, Transport_stock[w in 1:number_of_warehouses, q in 1:number_of_warehouses, t in 2:number_of_simulation_periods], sum(quantities_send[w,q,t] for q in 1:number_of_warehouses) <= quantities_stocked[w,t-1])
    @constraint(model_OiH, Transport_stock_start[w in 1:number_of_warehouses, q in 1:number_of_warehouses], sum(quantities_send[w,q,1] for q in 1:number_of_warehouses) <= initial_stock[w])
    # Constraint on quantity stock at time t with input and output
    @constraint(model_OiH, Stockage[w in 1:number_of_warehouses, t in 2:number_of_simulation_periods], quantities_stocked[w,t] == quantities_stocked[w,t-1]+quantities_ordered[w,t]
    +sum(quantities_recieved[w,q,t] - quantities_send[w,q,t] for q in 1:number_of_warehouses)- demand_coffee[w,t] + quantities_missed[w,t])
    @constraint(model_OiH, Stockage_start[w in 1:number_of_warehouses], quantities_stocked[w,1] == initial_stock[w]+quantities_ordered[w,1]
    +sum(quantities_recieved[w,q,1] - quantities_send[w,q,1] for q in 1:number_of_warehouses)- demand_coffee[w,1] + quantities_missed[w,1])


    optimize!(model_OiH)


    #Check if optimal solution was found
    if termination_status(model_OiH) == MOI.OPTIMAL
        println("Optimal solution found")
        
        # Display of the results in a text file
        
        # # Get the directory of the current script
        # script_directory = @__DIR__
        # # Construct the full file path
        # file_path = joinpath(script_directory, "output.txt")
        # # Open or create a text file
        # file = open(file_path, "w")

        # # Write inside the text file
        # println(file,"Cost of the solution : $(round.(objective_value(model_OiH), digits=2))")
        # println(file,"-----------------")
        # for t in 1:number_of_simulation_periods
        #     println(file, "Time Step: $t")
        #     println(file,"----") 
        #     # Display other information for the current time step
        #     for w in 1:number_of_warehouses 
        #         println(file,"Warehouse $w : Demand $(demand_coffee[w,t]) / Ordered $(round.(value.(quantities_ordered)[w,t], digits=2)) / Price $(cost_coffee[t][w])")
        #         if t != 1 
        #             println(file,"Previous Stock $(round.(value.(quantities_stocked)[w,t-1], digits=2)) / Sent $(sum(round.(value.(quantities_send)[w,q,t], digits=2) for q in 1:number_of_warehouses)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q,t], digits=2) for q in 1:number_of_warehouses))")
        #         else 
        #             println(file,"Previous Stock 2.00 / Sent $(sum(round.(value.(quantities_send)[w,q,t], digits=2) for q in 1:number_of_warehouses)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q,t], digits=2) for q in 1:number_of_warehouses))")
        #         end
        #         println(file,"Missed $(round.(value.(quantities_missed)[w,t], digits=2)) / Stock $(round.(value.(quantities_stocked)[w,t], digits=2))")
        #         println(file,"----")  
        #     end
        #     println(file,"-----------------")  # Separator between time steps
        # end

        # # Flush the file to ensure all data is written
        # flush(file)
        # # Close the file
        # close(file)
        # # Open the file 
        # run(`cmd /c start notepad $file_path`)
        

        #return interesting values
        return value.(quantities_ordered),value.(quantities_send),value.(quantities_recieved),value.(quantities_stocked),value.(quantities_missed),objective_value(model_OiH)
    else
        return error("No solution.")
    end

end

# qo_OiH,qs_OiH,qr_OiH,qst_OiH,qm_OiH,ov_OiH=Calculate_OiH_solution(prices, new_prices)