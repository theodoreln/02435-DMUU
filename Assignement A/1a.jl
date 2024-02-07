#Import packages
using JuMP
using Gurobi
using Printf

#Declare the size of the problem
# Number of district 
W = 3
# Number of time steps
T = 10

#Creation of the parameters of the problem
# Creation of the price list, W rows and T columns, here we use random number between 0 and 10
cost_coffee = round.(10 * rand(W, T), digits=2)
# Creation of the penalty price, number between 10 and 20
cost_miss = round.(10 * rand(W) .+ 10, digits=2)
# Creation of the transportation cost, W rows and W columns, 5 everywhere
cost_tr = 3*ones(W, W)
# Creation of the storage capacities, 1 row and W columns, 10 everywhere
warehouse_capacities = 10*ones(W)
# Creation of the transport capacities, W rows and W columns, 3 everywhere
transport_capacities = 3*ones(W,W)
# Creation of the demand, W rows and T columns, 10 everywhere
demand_coffee = round.(10 * ones(W, T), digits=2)


#Declare model with Gurobi solver
model_one = Model(Gurobi.Optimizer)

#Declare the variables to optimize
# Quantities of coffee ordered, W rows and T columns
@variable(model_one, quantities_ordered[1:W, 1:T]>=0)
# Quantities send from w to q, W rows W columns and T layers
@variable(model_one, quantities_send[1:W, 1:W, 1:T]>=0)
# Quantities recieved by w from q, W rows W columns and T layers
@variable(model_one, quantities_recieved[1:W, 1:W, 1:T]>=0)
# Quantities in the warehouse stockage, W rows and T columns
@variable(model_one, quantities_stocked[1:W, 1:T]>=0)
# Quantities mising to complete the demand, W rows and T columns
@variable(model_one, quantities_missed[1:W, 1:T]>=0)

#Objective function
@objective(model_one, Min, sum(quantities_ordered[w,t]*cost_coffee[w,t] for w in 1:W, t in 1:T) 
+ sum(quantities_send[w,q,t]*cost_tr[w,q] for w in 1:W, q in 1:W, t in 1:T)
+ sum(quantities_missed[w,t]*cost_miss[w] for w in 1:W, t in 1:T))

#Constraints of the problem
# Constraint on stockage capacities limited to the maximum capacities
@constraint(model_one, Stockage_limit[w in 1:W, t in 1:T], quantities_stocked[w,t] <= warehouse_capacities[w])
# Constraint on transport capacities limited to the maximum capacities
@constraint(model_one, Transport_limit[w in 1:W, q in 1:W, t in 1:T], quantities_send[w,q,t] <= transport_capacities[w,q])
# Constraint on quantity send equal quantity recieved
@constraint(model_one, Send_recieved[w in 1:W, q in 1:W, t in 1:T], quantities_send[w,q,t] == quantities_recieved[q,w,t])
# Constraint on a warehouse can only send to others warehouse
@constraint(model_one, Self_transport[w in 1:W, t in 1:T], quantities_send[w,w,t] == 0)
# Constraint on quantity send limited to previous stock
@constraint(model_one, Transport_stock[w in 1:W, q in 1:W, t in 2:T], sum(quantities_send[w,q,t] for q in 1:W) <= quantities_stocked[w,t-1])
@constraint(model_one, Transport_stock_start[w in 1:W, q in 1:W], sum(quantities_send[w,q,1] for q in 1:W) <= 2)
# Constraint on quantity stock at time t with input and output
@constraint(model_one, Stockage[w in 1:W, t in 2:T], quantities_stocked[w,t] == quantities_stocked[w,t-1]+quantities_ordered[w,t]
+sum(quantities_recieved[w,q,t] - quantities_send[w,q,t] for q in 1:W)- demand_coffee[w,t] + quantities_missed[w,t])
@constraint(model_one, Stockage_start[w in 1:W], quantities_stocked[w,1] == 2+quantities_ordered[w,1]
+sum(quantities_recieved[w,q,1] - quantities_send[w,q,1] for q in 1:W)- demand_coffee[w,1] + quantities_missed[w,1])


optimize!(model_one)


#Check if optimal solution was found
if termination_status(model_one) == MOI.OPTIMAL
    println("Optimal solution found")

    # Display of the results
    # Get the directory of the current script
    script_directory = @__DIR__
    # Construct the full file path
    file_path = joinpath(script_directory, "output.txt")
    # Open or create a text file
    file = open(file_path, "w")

    # Write inside the text file
    println(file,"Cost of the solution : $(round.(objective_value(model_one), digits=2)) $ ")
    println(file,"-----------------")
    for t in 1:T
        println(file, "Time Step: $t")
        println(file,"----") 
        # Display other information for the current time step
        for w in 1:W 
            println(file,"Warehouse $w : Demand $(demand_coffee[w,t]) / Ordered $(round.(value.(quantities_ordered)[w,t], digits=2)) / Price $(cost_coffee[w,t])")
            if t != 1 
                println(file,"Previous Stock $(round.(value.(quantities_stocked)[w,t-1], digits=2)) / Sent $(sum(round.(value.(quantities_send)[w,q,t], digits=2) for q in 1:W)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q,t], digits=2) for q in 1:W))")
            else 
                println(file,"Previous Stock 2.00 / Sent $(sum(round.(value.(quantities_send)[w,q,t], digits=2) for q in 1:W)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q,t], digits=2) for q in 1:W))")
            end
            println(file,"Missed $(round.(value.(quantities_missed)[w,t], digits=2)) / Stock $(round.(value.(quantities_stocked)[w,t], digits=2))")
            println(file,"----")  
        end
        println(file,"-----------------")  # Separator between time steps
    end

    # Flush the file to ensure all data is written
    flush(file)
    # Close the file
    close(file)
    # Open the file 
    run(`cmd /c start notepad $file_path`)

else
    error("No solution.")
end