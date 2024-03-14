using Random
include("price_process.jl")

function simulation_experiments_creation(number_of_warehouses, W, number_of_periods)

    number_of_experiments = 40
    Expers = collect(1:number_of_experiments)

    initial_prices = zeros(number_of_experiments, number_of_warehouses)
    price_trajectory = zeros(number_of_experiments, number_of_warehouses, number_of_periods)
    for e in Expers
        for w in W
            initial_prices[e,w] = rand()*10
            price_trajectory[e,w,1] = initial_prices[e,w]
            for t in 2:number_of_periods
                price_trajectory[e,w,t] = sample_next(price_trajectory[e,w,t-1])
            end
        end
    end


    return number_of_experiments, Expers, price_trajectory
end
#end