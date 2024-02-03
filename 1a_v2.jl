#Import packages
using JuMP
using Gurobi
using Printf


#Declare model with Gurobi solver
model_one = Model(Gurobi.Optimizer)


#Declare the size of the problem
# Number of district 
W = 3
# Number of time steps
T = 10


#Creation of the parameters of the problem
# Creation of the price list, W rows and T columns, here we use random number between 0 and 10
cost_coffee = round.(10 * rand(W, T), digits=2)
# Creation of the penalty price, constructed for 3 here
cost_miss = [10,15,20]
# Creation of the transportation cost, W rows and W columns, 5 everywhere
cost_tr = 3*ones(W, W)
# Creation of the storage capacities, 1 row and W columns, 10 everywhere
warehouse_capacities = 10*ones(W)
# Creation of the transport capacities, W rows and W columns, 3 everywhere
transport_capacities = 3*ones(W,W)
# Creation of the demand, W rows and T columns, randon between 5 and 15
demand_coffee = round.(12 * rand(W, T), digits=2)


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
# Constraint on stockage capacities limited to the maximum capacities and should be above the provided quantity
for w in 1:W
    for t in 1:T
        @constraint(model_one, quantities_stocked[w,t] <= warehouse_capacities[w])
        @constraint(model_one, quantities_stocked[w,t] >= demand_coffee[w,t] - quantities_missed[w,t])
    end
end

# Constraint on the quantity stocked at day t function of what's left from the day before and what is coming
for w in 1:W
    for t in 1:T
        if t == 1
            @constraint(model_one, quantities_stocked[w,t] == quantities_ordered[w,t]
            +sum(quantities_recieved[w,q,t] - quantities_send[w,q,t] for q in 1:W))
        else 
            @constraint(model_one, quantities_stocked[w,t] == quantities_stocked[w,t-1] +quantities_ordered[w,t]
            +sum(quantities_recieved[w,q,t] - quantities_send[w,q,t] for q in 1:W)
            - demand_coffee[w,t-1] + quantities_missed[w,t-1])
        end
    end
end

# Constraints on quantities send and recieved
for w in 1:W
    for q in 1:W
        for t in 1:T
            # Quantities send and recieved limited by the transport capacity
            @constraint(model_one, quantities_send[w,q,t] <= transport_capacities[w,q])
            @constraint(model_one, quantities_recieved[w,q,t] <= transport_capacities[q,w])
            # Quantity recieved by q from w equal the quantity send by w to q
            @constraint(model_one, quantities_send[w,q,t] == quantities_recieved[q,w,t])
            if t == 1
                # Quantities send and recieved at t = 1 are nul because there is no stock at the beginning
                @constraint(model_one, quantities_send[w,q,t] == 0)
                @constraint(model_one, quantities_recieved[w,q,t] == 0)
            end
            if t != 1
                # Quantities send and recieved are limited to the quantities left from the day before
                @constraint(model_one, quantities_send[w,q,t] <= quantities_stocked[w,t-1] - demand_coffee[w,t-1] + quantities_missed[w,t-1])
                @constraint(model_one, quantities_recieved[w,q,t] <= quantities_stocked[q,t-1] - demand_coffee[q,t-1] + quantities_missed[q,t-1])
            end
            if w == q 
                # Can't exchange quantitites between the same warehouse
                @constraint(model_one, quantities_send[w,q,t] == 0)
                @constraint(model_one, quantities_recieved[w,q,t] == 0)
            end
        end
    end
end

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
    for t in 1:T
        println(file, "Time Step: $t")
        println(file,"----") 
        # Display other information for the current time step
        for w in 1:W 
            println(file,"Warehouse $w : Ordered $(round.(value.(quantities_ordered)[w,t], digits=2)) / Price $(cost_coffee[w,t])")
            if t != 1 
                println(file,"Remaining Stock $(round.(value.(quantities_stocked)[w,t-1] - demand_coffee[w,t-1] + value.(quantities_missed)[w,t-1], digits=2)) / Sent $(sum(round.(value.(quantities_send)[w,q,t], digits=2) for q in 1:W)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q,t], digits=2) for q in 1:W))")
            else 
                println(file,"Remaining Stock 0 / Sent $(sum(round.(value.(quantities_send)[w,q,t], digits=2) for q in 1:W)) / Recieved $(sum(round.(value.(quantities_recieved)[w,q,t], digits=2) for q in 1:W))")
            end
            println(file,"Stock $(round.(value.(quantities_stocked)[w,t], digits=2)) / Demand $(demand_coffee[w,t])")
            println(file,"Provided $(round.(demand_coffee[w,t] - value.(quantities_missed)[w,t], digits=2)) / Missed $(round.(value.(quantities_missed)[w,t], digits=2))")
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
