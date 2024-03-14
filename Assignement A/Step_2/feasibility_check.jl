function check_feasibility(x, send, receive, z, m, current_stock, current_demands, warehouse_capacities, transport_capacities)
    

    function any_negative_element(matrices...)
        for matrix in matrices
            if any(matrix .< 0 - 0.01)
                return true  # Returns true if any element is negative in any matrix
            end
        end
        return false  # Returns false if no element is negative in any matrix
    end
    
    function any_element_outside_epsilon(matrices...; epsilon=0.1)
        for matrix in matrices
            if any(matrix .> epsilon) || any(matrix .< -epsilon)
                return true  # Returns true if any element satisfies the condition
            end
        end
        return false  # Returns false if no element satisfies the condition
    end

    balance_residual = zeros(number_of_warehouses)
    capacity_residual = ones(number_of_warehouses)
    transport_residual = ones(number_of_warehouses, number_of_warehouses)
    consistency_residual = zeros(number_of_warehouses, number_of_warehouses)
    stock_positivity = ones(number_of_warehouses)
    delay_residual = ones(number_of_warehouses)

    #Check decisions
    for w in W
            stock_positivity[w] = z[w]
            delay_residual[w] = current_stock[w] - sum(send[w,q] for q in W)
            balance_residual[w] = current_stock[w] + x[w] + sum(receive[w,q] - send[w,q] for q in W) - current_demands[w] + m[w] - z[w]
            capacity_residual[w] = warehouse_capacities[w] - current_stock[w]
        for q in W
            transport_residual[w,q] = transport_capacities[w,q] - send[w,q]
            consistency_residual[w,q] = send[w,q] - receive[q,w]
            
        end    
    end

    if any_negative_element(stock_positivity, delay_residual, transport_residual, capacity_residual) || any_element_outside_epsilon(balance_residual, consistency_residual)
        println("THE DECISIONS DO NOT RESPECT THE CONSTRAINTS. CHECK THAT THE CONSTRAINTS FOR THE CURRENT TIMESLOT ARE CORRECTLY IMPLEMENTED")
        success = 0
    else
        success = 1
    end

    return success

end