# Import useful files
include("two_stage_problem_data.jl")
include("price_process.jl")
include("ScenarioSelections.jl")
include("Expected_Value.jl")
include("Optimal_Hindsight.jl")
include("Stochastic_now.jl")
include("DayTwo_optimal.jl")

using JLD
using Plots

####################################################
# Evaluating the different functions for 100 samples 
####################################################

# Generating 100 values for the initial prices 
number_experiments = 100

function Evaluation_comparison(number_experiments)

    # Initialize the vector to keep the values we want to follow for each experiment 
    Prices_day1 = zeros(number_experiments,3)
    Prices_day2 = zeros(number_experiments,3)
    Qu_ordered_day1 = zeros(number_experiments,3, 5)
    Qu_send_day1 = zeros(number_experiments, 3, 3, 5)
    Qu_rec_day1 = zeros(number_experiments, 3, 3, 5)
    Qu_stock_day1 = zeros(number_experiments, 3, 5)
    Qu_miss_day1 = zeros(number_experiments, 3, 5)
    Qu_ordered_day2 = zeros(number_experiments, 3, 5)
    Qu_send_day2 = zeros(number_experiments, 3, 3, 5)
    Qu_rec_day2 = zeros(number_experiments, 3, 3, 5)
    Qu_stock_day2 = zeros(number_experiments, 3, 5)
    Qu_miss_day2 = zeros(number_experiments, 3, 5)
    Cost = zeros(number_experiments, 5)

    # Iterate for 100 experiments
    for e in range(1,number_experiments) 

        # Generate prices of day 1 for this experiment
        prices_exp = round.(10 * rand(3), digits=2)
        Prices_day1[e,:] = prices_exp

        # Here and now decisions for EV 
        qo_EB,qs_EB,qr_EB,qst_EB,qm_EB,cost_EB=Make_EV_here_and_now_decision(prices_exp)
        Qu_ordered_day1[e,:,1], Qu_stock_day1[e,:,1], Qu_miss_day1[e,:,1]  = qo_EB, qst_EB, qm_EB
        Qu_send_day1[e,:,:,1], Qu_rec_day1[e,:,:,1], Cost[e,1] = qs_EB, qr_EB, cost_EB
        # Here and now decisions for Stoch, N = 5 
        qo_ST_5,qs_ST_5,qr_ST_5,qst_ST_5,qm_ST_5,cost_ST_5=Make_Stochastic_here_and_now_decision(prices_exp,5)
        Qu_ordered_day1[e,:,2], Qu_stock_day1[e,:,2], Qu_miss_day1[e,:,2]  = qo_ST_5, qst_ST_5, qm_ST_5
        Qu_send_day1[e,:,:,2], Qu_rec_day1[e,:,:,2], Cost[e,2] = qs_ST_5, qr_ST_5, cost_ST_5
        # Here and now decisions for Stoch, N = 20
        qo_ST_20,qs_ST_20,qr_ST_20,qst_ST_20,qm_ST_20,cost_ST_20=Make_Stochastic_here_and_now_decision(prices_exp,20)
        Qu_ordered_day1[e,:,3], Qu_stock_day1[e,:,3], Qu_miss_day1[e,:,3]  = qo_ST_20, qst_ST_20, qm_ST_20
        Qu_send_day1[e,:,:,3], Qu_rec_day1[e,:,:,3], Cost[e,3] = qs_ST_20, qr_ST_20, cost_ST_20
        # Here and now decisions for Stoch, N = 50 
        qo_ST_50,qs_ST_50,qr_ST_50,qst_ST_50,qm_ST_50,cost_ST_50=Make_Stochastic_here_and_now_decision(prices_exp,50)
        Qu_ordered_day1[e,:,4], Qu_stock_day1[e,:,4], Qu_miss_day1[e,:,4]  = qo_ST_50, qst_ST_50, qm_ST_50
        Qu_send_day1[e,:,:,4], Qu_rec_day1[e,:,:,4], Cost[e,4] = qs_ST_50, qr_ST_50, cost_ST_50

        # Generate prices of day 2 for this experiment
        new_prices_exp = round.(map(Float64,map(sample_next, prices_exp)),digits=2)
        Prices_day2[e,:] = new_prices_exp

        # Optimizing in day two with decision from day one 
        qo_D2,qs_D2,qr_D2,qst_D2,qm_D2,ov_D2 = optimal_stage2_decision(qst_EB, new_prices_exp)
        Qu_ordered_day2[e,:,1], Qu_stock_day2[e,:,1], Qu_miss_day2[e,:,1]  = qo_D2, qst_D2, qm_D2
        Qu_send_day2[e,:,:,1], Qu_rec_day2[e,:,:,1] = qs_D2, qr_D2
        Cost[e,1] = Cost[e,1] + ov_D2
        qo_D2,qs_D2,qr_D2,qst_D2,qm_D2,ov_D2 = optimal_stage2_decision(qst_ST_5, new_prices_exp)
        Qu_ordered_day2[e,:,2], Qu_stock_day2[e,:,2], Qu_miss_day2[e,:,2]  = qo_D2, qst_D2, qm_D2
        Qu_send_day2[e,:,:,2], Qu_rec_day2[e,:,:,2] = qs_D2, qr_D2
        Cost[e,2] = Cost[e,2] + ov_D2
        qo_D2,qs_D2,qr_D2,qst_D2,qm_D2,ov_D2 = optimal_stage2_decision(qst_ST_20, new_prices_exp)
        Qu_ordered_day2[e,:,3], Qu_stock_day2[e,:,3], Qu_miss_day2[e,:,3]  = qo_D2, qst_D2, qm_D2
        Qu_send_day2[e,:,:,3], Qu_rec_day2[e,:,:,3] = qs_D2, qr_D2
        Cost[e,3] = Cost[e,3] + ov_D2
        qo_D2,qs_D2,qr_D2,qst_D2,qm_D2,ov_D2 = optimal_stage2_decision(qst_ST_50, new_prices_exp)
        Qu_ordered_day2[e,:,4], Qu_stock_day2[e,:,4], Qu_miss_day2[e,:,4]  = qo_D2, qst_D2, qm_D2
        Qu_send_day2[e,:,:,4], Qu_rec_day2[e,:,:,4] = qs_D2, qr_D2
        Cost[e,4] = Cost[e,4] + ov_D2

        # Optimal Hindshight optimization
        qo_OiH,qs_OiH,qr_OiH,qst_OiH,qm_OiH,ov_OiH=Calculate_OiH_solution(prices_exp, new_prices_exp)
        Qu_ordered_day1[e,:,5], Qu_stock_day1[e,:,5], Qu_miss_day1[e,:,5]  = qo_OiH[:,1], qst_OiH[:,1], qm_OiH[:,1]
        Qu_send_day1[e,:,:,5], Qu_rec_day1[e,:,:,5] = qs_OiH[:,:,1], qr_OiH[:,:,1]
        Qu_ordered_day2[e,:,5], Qu_stock_day2[e,:,5], Qu_miss_day2[e,:,5]  = qo_OiH[:,2], qst_OiH[:,2], qm_OiH[:,2]
        Qu_send_day2[e,:,:,5], Qu_rec_day2[e,:,:,5], Cost[e,5] = qs_OiH[:,:,2], qr_OiH[:,:,2], ov_OiH

    end

    # Output the values in a file
    script_directory = @__DIR__
    file_path = joinpath(script_directory, "results_ff.jld")
    save(file_path, "Prices_day1", Prices_day1, "Prices_day2", Prices_day2, "Qu_ordered_day1", Qu_ordered_day1, "Qu_send_day1", Qu_send_day1, "Qu_rec_day1", Qu_rec_day1, 
    "Qu_stock_day1", Qu_stock_day1, "Qu_miss_day1", Qu_miss_day1, "Qu_ordered_day2", Qu_ordered_day2, "Qu_send_day2", Qu_send_day2, "Qu_rec_day2", Qu_rec_day2, 
    "Qu_stock_day2", Qu_stock_day2, "Qu_miss_day2", Qu_miss_day2, "Cost", Cost)
end

# To launch the computation, watch out for the number of experiments, it can take a lot of time
Evaluation_comparison(number_experiments)

# Load back the result file 
script_directory = @__DIR__
file_path = joinpath(script_directory, "results_ff.jld")
dict = load(file_path)

Prices_day1, Prices_day2, Qu_ordered_day1, Qu_send_day1, Qu_rec_day1, Qu_stock_day1, Qu_miss_day1, Qu_ordered_day2, Qu_send_day2, Qu_rec_day2, Qu_stock_day2, Qu_miss_day2, Cost =
 dict["Prices_day1"], dict["Prices_day2"], dict["Qu_ordered_day1"], dict["Qu_send_day1"], dict["Qu_rec_day1"], dict["Qu_stock_day1"], dict["Qu_miss_day1"], dict["Qu_ordered_day2"], dict["Qu_send_day2"], dict["Qu_rec_day2"], dict["Qu_stock_day2"], dict["Qu_miss_day2"], dict["Cost"]


function Display_results_one_experiment(Dict, no_experiment, number_experiments)
    if no_experiment > number_experiments
        error("Index experiment exceeds the maximum number of the experiment")
    end

    # Load the data 
    Prices_day1, Prices_day2, Qu_ordered_day1, Qu_send_day1, Qu_rec_day1, Qu_stock_day1, Qu_miss_day1, Qu_ordered_day2, Qu_send_day2, Qu_rec_day2, Qu_stock_day2, Qu_miss_day2, Cost =
    dict["Prices_day1"], dict["Prices_day2"], dict["Qu_ordered_day1"], dict["Qu_send_day1"], dict["Qu_rec_day1"], dict["Qu_stock_day1"], dict["Qu_miss_day1"], dict["Qu_ordered_day2"], dict["Qu_send_day2"], dict["Qu_rec_day2"], dict["Qu_stock_day2"], dict["Qu_miss_day2"], dict["Cost"]

    # Dict to deal with the name of each method
    Dict_meth = Dict(1 => "EV", 2 => "ST_5", 3 => "ST_20", 4 => "ST_50", 5 => "OiH")

    # Write a note
    # Get the directory of the current script
    script_directory = @__DIR__
    # Construct the full file path
    file_path = joinpath(script_directory, "Evaluation", "evaluation "*string(no_experiment)*".txt")
    # Open or create a text file
    file = open(file_path, "w")

    println(file, "Evaluation of the different method for the optimization of the warehouse problem")
    println(file, " ")
    for m in 1:5 
        println(file, "Method : $(Dict_meth[m]) / Cost of the method : $(Cost[no_experiment,m])")
        println(file,"----")
        println(file, "Day 1")
        println(file,"----") 
        for w in 1:3 
            println(file,"Warehouse $w : Demand 4.00 / Ordered $(round.(Qu_ordered_day1[no_experiment,w,m], digits=2)) / Price $(Prices_day1[no_experiment,w])")
            println(file,"Previous Stock 2.00 / Sent $(sum(Qu_send_day1[no_experiment,w,:,m])) / Recieved $(sum(Qu_rec_day1[no_experiment,w,:,m]))")
            println(file,"Missed $(round.(Qu_miss_day1[no_experiment,w,m], digits=2)) / Stock $(round.(Qu_stock_day1[no_experiment,w,m], digits=2))")
            println(file,"----")  
        end
        println(file,"----")
        println(file, "Day 2")
        println(file,"----") 
        for w in 1:3 
            println(file,"Warehouse $w : Demand 4.00 / Ordered $(round.(Qu_ordered_day2[no_experiment,w,m], digits=2)) / Price $(Prices_day2[no_experiment,w])")
            println(file,"Previous Stock $(round.(Qu_stock_day1[no_experiment,w,m], digits=2)) / Sent $(sum(Qu_send_day2[no_experiment,w,:,m])) / Recieved $(sum(Qu_rec_day2[no_experiment,w,:,m]))")
            println(file,"Missed $(round.(Qu_miss_day2[no_experiment,w,m], digits=2)) / Stock $(round.(Qu_stock_day2[no_experiment,w,m], digits=2))")
            println(file,"----")  
        end
        println(file, "----------------------------------------")
        println(file, " ")
    end

    # Flush the file to ensure all data is written
    flush(file)
    # Close the file
    close(file)
    # Open the file 
    # run(`cmd /c start notepad $file_path`)

end

for no_experiment in 1:number_experiments
    Display_results_one_experiment(Dict, no_experiment, number_experiments)
end

# Print mean value of each method
println("Mean value of EV method "*string(mean(Cost[:,1])))
println("Mean value of ST-5 method "*string(mean(Cost[:,2])))
println("Mean value of ST_20 method "*string(mean(Cost[:,3])))
println("Mean value of ST_50 method "*string(mean(Cost[:,4])))
println("Mean value of OiH method "*string(mean(Cost[:,5])))

# Plot the values of all method together
# Extracting columns for plotting
x_values = 1:size(Cost,1)
y_values = [Cost[:, i] for i in 1:size(Cost, 2)]

# Plotting the curves
plot(x_values, y_values, label=["EV" "ST_5" "ST_20" "ST_50" "OiH"], xlabel="Experiments", ylabel="Cost", title="Cost of all methods on all experiments")
