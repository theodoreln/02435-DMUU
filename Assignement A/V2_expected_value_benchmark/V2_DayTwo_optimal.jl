#Import useful files
include("V2_02435_two_stage_problem_data.jl")

#Import packages
using Random
using JuMP
using Gurobi
using Printf


function optimal_stage2_decision(day_two_stock, day_two_prices)
    #Importation of the inputs from the two stage problem
    number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data(1)

    initial_stock = day_two_stock
    cost_coffee = day_two_prices
    # Creation of the demand, W rows and T columns, 10 everywhere
    demand_coffee = demand_trajectory

    #Declare model with Gurobi solver
    model_DTWO = Model(Gurobi.Optimizer)

    #Declare the variables to optimize
    # Quantities of coffee ordered, W rows 
    @variable(model_DTWO, quantities_ordered[1:number_of_warehouses]>=0)
    # Quantities send from w to q, W rows and W columns
    @variable(model_DTWO, quantities_send[1:number_of_warehouses, 1:number_of_warehouses]>=0)
    # Quantities recieved by w from q, W rows and W columns 
    @variable(model_DTWO, quantities_recieved[1:number_of_warehouses, 1:number_of_warehouses]>=0)
    # Quantities in the warehouse stockage, W rows
    @variable(model_DTWO, quantities_stocked[1:number_of_warehouses]>=0)
    # Quantities mising to complete the demand, W rows
    @variable(model_DTWO, quantities_missed[1:number_of_warehouses]>=0)
    
    #Objective function
    @objective(model_DTWO, Min, sum(quantities_ordered[w]*cost_coffee[w] for w in 1:number_of_warehouses) 
    + sum(quantities_send[w,q]*cost_tr[w,q] for w in 1:number_of_warehouses, q in 1:number_of_warehouses)
    + sum(quantities_missed[w]*cost_miss[w] for w in 1:number_of_warehouses))

    #Constraints of the problem
    # Constraint on stockage capacities limited to the maximum capacities
    @constraint(model_DTWO, Stockage_limit[w in 1:number_of_warehouses], quantities_stocked[w] <= warehouse_capacities[w])
    # Constraint on transport capacities limited to the maximum capacities
    @constraint(model_DTWO, Transport_limit[w in 1:number_of_warehouses, q in 1:number_of_warehouses], quantities_send[w,q] <= transport_capacities[w,q])
    # Constraint on quantity send equal quantity recieved
    @constraint(model_DTWO, Send_recieved[w in 1:number_of_warehouses, q in 1:number_of_warehouses], quantities_send[w,q] == quantities_recieved[q,w])
    # Constraint on a warehouse can only send to others warehouse
    # Useless cause the self-transport capacity is equal to 0
    @constraint(model_DTWO, Self_transport[w in 1:number_of_warehouses], quantities_send[w,w] == 0)
    # Constraint on quantity send limited to previous stock
    @constraint(model_DTWO, Transport_stock_start[w in 1:number_of_warehouses, q in 1:number_of_warehouses], sum(quantities_send[w,q] for q in 1:number_of_warehouses) <= initial_stock[w])
    # Constraint on quantity stock at time t with input and output
    @constraint(model_DTWO, Stockage_start[w in 1:number_of_warehouses], quantities_stocked[w] == initial_stock[w]+quantities_ordered[w]
    +sum(quantities_recieved[w,q] - quantities_send[w,q] for q in 1:number_of_warehouses)- demand_coffee[w] + quantities_missed[w])


    optimize!(model_DTWO)


    #Check if optimal solution was found
    if termination_status(model_DTWO) == MOI.OPTIMAL
        println("Optimal solution found")
        
        # Display of the results in a text file
        

        # Get the directory of the current script
        script_directory = @__DIR__
        # Construct the full file path
        file_path = joinpath(script_directory, "output.txt")
        # Open or create a text file
        file = open(file_path, "w")

        # Write inside the text file
        println(file,"Cost of the solution : $(round.(objective_value(model_DTWO), digits=2))")
        println(file,"-----------------")
        for w in 1:number_of_warehouses 
            println(file,"Warehouse $w : Demand $(demand_coffee[w]) / Ordered $(round.(value.(quantities_ordered)[w], digits=2)) / Price $(cost_coffee[w])")
            println(file,"Previous Stock 2.00 / Sent $(sum(round.(value.(quantities_send)[w,q], digits=2) for q in 1:number_of_warehouses)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q], digits=2) for q in 1:number_of_warehouses))")
            println(file,"Missed $(round.(value.(quantities_missed)[w,], digits=2)) / Stock $(round.(value.(quantities_stocked)[w], digits=2))")
            println(file,"----")  
        end

        # Flush the file to ensure all data is written
        flush(file)
        # Close the file
        close(file)
        # Open the file 
        run(`cmd /c start notepad $file_path`)
        

        #return interesting values
        return value.(quantities_ordered[:,1]),value.(quantities_send[:,:,1]),value.(quantities_recieved[:,:,1]),value.(quantities_stocked[:,1]),value.(quantities_missed[:,1]),objective_value(model_DTWO)
    else
        return error("No solution.")
    end

end

prices = round.(10 * rand(3), digits=2)
qo_ST,qs_ST,qr_ST,qst_ST,qm_ST,ov_ST=Make_Stochastic_here_and_now_decision(prices,50)
day_two_prices = round.(10 * rand(3), digits=2)
optimal_stage2_decision(qst_ST, day_two_prices)