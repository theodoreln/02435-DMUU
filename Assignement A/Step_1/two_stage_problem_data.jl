using Random

function load_the_data(number_of_days)

    number_of_warehouses = 3
    W = collect(1:number_of_warehouses)

    number_of_simulation_periods = number_of_days
    sim_T = collect(1:number_of_simulation_periods)

    #Cost of missing demand at w
    #Call by cost_miss[w]
    cost_miss = [10,15,20]

    #Distance-based transportation cost for each pair of warehouses w1 and w2
    #Call by cost_tr[w1,w2]
    cost_tr = ones(number_of_warehouses, number_of_warehouses)*5

    #Capacity of warehouse w
    warehouse_capacities = 10*ones(number_of_warehouses)

    #Capacity of the transportation link for each pair of warehouses
    #Call by transport_capacities[w1,w2]
    transport_capacities = 4*ones(number_of_warehouses, number_of_warehouses)
    transport_capacities[3,1] = 0
    transport_capacities[1,3] = 0
    for w in W
        for q in W
            if w == q
                transport_capacities[w,q] = 0
            end
        end
    end

    #Initial stock of at w
    initial_stock = 2*ones(number_of_warehouses)

    demand_trajectory = 4*ones(number_of_warehouses, number_of_simulation_periods)

    return number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory
end

#end  # End of module
