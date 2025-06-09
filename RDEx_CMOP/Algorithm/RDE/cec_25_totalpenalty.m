function total = cec_25_totalpenalty(func_num, PopulG, PopulH, PopulAllConstr)
    global tempG tempH ng_B nh_B epsilon0001

    for i = 1:3
        PopulG(i) = tempG(i);
    end

    for i = 1:6
        PopulH(i) = tempH(i);
        if PopulH(i) < epsilon0001 && PopulH(i) > -epsilon0001
            PopulH(i) = 0;
        end
    end

    total_constr = ng_B(func_num) + nh_B(func_num);
    for i = 1:total_constr
        if i <= ng_B(func_num)
            PopulAllConstr(i) = tempG(i);
        else
            h_idx = i - ng_B(func_num);
            if tempH(h_idx) >= epsilon0001
                PopulAllConstr(i) = tempH(h_idx);
            elseif tempH(h_idx) <= -epsilon0001
                PopulAllConstr(i) = -tempH(h_idx);
            else
                PopulAllConstr(i) = 0;
            end
        end
    end

    total = 0;
    for i = 1:ng_B(func_num)
        if tempG(i) > 0
            total = total + tempG(i);
        end
    end

    for i = 1:nh_B(func_num)
        if tempH(i) >= epsilon0001
            total = total + tempH(i);
        elseif tempH(i) <= -epsilon0001
            total = total + abs(tempH(i));
        end
    end

    total = total / (ng_B(func_num) + nh_B(func_num));
end