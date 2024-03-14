# This code applies and evaluates a policy over 40 different experiments and gives you the average policy's cost
# Name you file "02435_multistage_policy" and put it in the same folder
# Then simply execute this file

# You should not change anything in this file (except for the file and function names in lines 9 and 48, if you want to evaluate different policies).


# Including your policy and the dummy policy
include("multistage_policy.jl")
include("dummy_policy.jl")

include("multistage_problem_data.jl")
# Loading the problem's parameters
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_sim_periods, sim_T, demand_trajectory = load_the_data()

include("simulation_experiments.jl")
# Creating the random experiments on which the policy will be evaluated
number_of_experiments, Expers, Price_experiments = simulation_experiments_creation(number_of_warehouses, W, number_of_sim_periods)

# Including a function that checks if your policy's decisions are feasible
include("feasibility_check.jl")

# Initialization of the decision variables and policy cost
x = Dict()
send = Dict()
receive = Dict()
z = Dict()
m = Dict()
policy_cost = 99999999*ones(number_of_experiments, number_of_sim_periods)
policy_cost_at_experiment = 99999999*ones(number_of_experiments)

# for each experiment
for e in Expers
    println("EXPERIMENT N°$e LAUNCHED")
    # and for each timeslot of the horizon
    for tau in sim_T
        println("    CURRENTLY OPTIMIZING DAY $tau")
        # Set each warehouse's stock level 
        if tau == 1
            current_stock = initial_stock
        else
            current_stock = z[(e,tau-1)]
        end
        # Observe current demands and prices
        current_demands = demand_trajectory[:,tau]
        current_prices = Price_experiments[e,:,tau]

        # Call policy to make a decision for here and now
        # This number of look ahead days does not include the current day
        look_ahead_days = 3
        nb_initial_scenarios = 1000 
        granularity = 0.5
        nb_reduced_scenarios = 60
        elapsed_time = @elapsed begin
            x[(e,tau)], send[(e,tau)], receive[(e,tau)], z[(e,tau)], m[(e,tau)] = make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices, look_ahead_days, nb_initial_scenarios, granularity, nb_reduced_scenarios)
        end

        # Print the elapsed time
        println("Elapsed time: $elapsed_time seconds")
        
        #Check whether the policy's here and now decisions are feasible/meaningful
        successful = check_feasibility(x[(e,tau)], send[(e,tau)], receive[(e,tau)], z[(e,tau)], m[(e,tau)], current_stock, current_demands,  warehouse_capacities, transport_capacities)
        # If not, then the policy's decisions are discarded for this timeslot, and the dummy policy is used instead
        if successful == 0
            println("DECISION DOES NOT MEET THE CONSTRAINTS FOR THIS TIMESLOT. THE DUMMY POLICY WILL BE USED INSTEAD")
            println(e, number_of_sim_periods, tau, current_stock, current_demands, x[(e,tau)], send[(e,tau)], receive[(e,tau)], z[(e,tau)], m[(e,tau)])
            global ep = e
            global taup = tau
            global current_stockp = current_stock
            global current_demandsp = current_demands
            global xp = x[(e,tau)]
            global sp = send[(e,tau)]
            global rp = receive[(e,tau)] 
            global zp = z[(e,tau)]
            global mp = m[(e,tau)]
            x[(e,tau)], send[(e,tau)], receive[(e,tau)], z[(e,tau)], m[(e,tau)] = make_dummy_decision(number_of_sim_periods, tau, current_stock, current_demands, current_prices)
        end

        policy_cost[e,tau] = sum(current_prices[w]*x[(e,tau)][w] + cost_miss[w]*m[(e,tau)][w] + sum(cost_tr[w,q]*receive[(e,tau)][w,q] for q in W) for w in W)  
    end
    policy_cost_at_experiment[e] = sum(policy_cost[e,tau] for tau in sim_T)
end

FINAL_POLICY_COST = sum(policy_cost_at_experiment[e] for e in Expers) / number_of_experiments
println("THE FINAL POLICY EXPECTED COST IS ", FINAL_POLICY_COST)