using JuMP
using Gurobi
using Random


function make_dummy_decision(number_of_sim_periods, tau, current_stock, current_prices)

    include("V2_02435_multistage_problem_data.jl")
    number_of_warehouses, W, c_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data()

    current_demands = demand_trajectory[:,tau]

    x = current_demands
    send = zeros(number_of_warehouses, number_of_warehouses)
    receive = zeros(number_of_warehouses, number_of_warehouses)
    z = current_stock
    m = zeros(number_of_warehouses)

    return x, send, receive, z, m
end