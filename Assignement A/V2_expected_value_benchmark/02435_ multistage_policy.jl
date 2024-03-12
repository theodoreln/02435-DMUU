#Import useful files
include("V2_02435_two_stage_problem_data.jl")
include("V2_price_process.jl")
include("V2_ScenarioSelections.jl")

#Import packages
using Random
using JuMP
using Gurobi
using Printf

# prices=round.(10 * rand(3), digits=2)

function Make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock,current_prices)
    Nb_look_ahead_day = 3
    Nb_of_reduced_scenarios = 50
    Nb_scenarios = 1000
    
    #Control the number of ahead day
    if Nb_look_ahead_day>=number_of_sim_periods-tau
        Actual_Nb_look_ahead_day = number_of_sim_periods-tau

    else 
        Actual_Nb_look_ahead_day = Nb_look_ahead_day
    
    end
    
    #Generate prices
    prices = zeros((length(W),Actual_Nb_look_ahead_day,Nb_scenarios))
    for w in  W
        for t in 1:Actual_Nb_look_ahead_day
            for s in 1:Nb_scenarios
                if t==1
                    prices[w,t,s] = sample_next(current_prices[w])
                else
                    prices[w,t,s] = sample_next(prices[w,t-1,s])
                end
            end
        end
    end


    #Discretize prices
    prices = round(prices,digits=2)


    #Reduce the number of scenarios
    reduced_next_prices, Probs = FastForward(prices, Nb_of_reduced_scenarios)

    #Populate
    Sets = fill([], (length(W),Actual_Nb_look_ahead_day,Nb_scenarios))

    #Can be optimized because we are writing several times the same thing (if two scenario share the same history at a stage)
    for w in  W
        for t in 1:Actual_Nb_look_ahead_day
            for s in 1:Nb_scenarios
                for d in 1:Nb_scenarios
                    if reduced_next_prices[w,:t-1,s]==reduced_next_prices[w,:t-1,d]
                        Sets[w,t,s].append!(d)
                    end
                end
            end
        end
    end

    

    return x[(e,tau)], send[(e,tau)], receive[(e,tau)], z[(e,tau)], m[(e,tau)]
end
