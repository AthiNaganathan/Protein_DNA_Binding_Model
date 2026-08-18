classdef binding_model
    methods
        % Count and all possible binding permutations 
        function result = generate_perms(~, v, n)
            for i = 1:n
                if i == 1
                    perms = {v};
                else
                    cur = zeros(repelem(length(v), i));
                    l_p = perms(i-1);
                    l_p = l_p{:};
                    for j = 1:length(v)
                        e = v(j);
        
                        idx_p = repmat({':'}, 1, i-1);
                        idx_v = repmat({':'}, 1, i-1);
                        if i <= 3
                            idx_p{1} = j:length(v);
                            idx_v{1} = 1:length(v)-j+1;
                        else
                            idx_p{end} = j:length(v);
                            idx_v{end} = 1:length(v)-j+1;
                        end
                        l_v = zeros(size(l_p));
                        l_v(idx_v{:}) = l_p(idx_p{:});
                        
                        idx_p = repmat({':'}, 1, i);
                        if i == 2
                            idx_p{1} = j;
                        else
                            idx_p{end} = j;
                        end
                        cur(idx_p{:}) = (e*(10^(i - 1))).*(l_v ~= 0) + l_v;
                    end
        
                    perms{1, end+1} = cur;
                end
            end
        
            result = [];
            for i = 1:length(perms)
                temp = perms{i};
                flat = reshape(temp, 1, []);
                result = [result flat];
            end
            result = result(result ~= 0);
            result = [0 result];
        end

        % Fill binding permutations with protein binding
        function f_results = fill_perms(~, perms, n_DNA, no_of_modes, p_sites)
            f_results = [];
        
            modes = zeros(no_of_modes, p_sites);
            for i = 1:p_sites:p_sites*no_of_modes
                modes(floor(i/p_sites)+1, :) = i:i+p_sites-1;
            end
        
            for perm = perms
                dgts = num2str(perm)-'0';
                results = zeros(1, n_DNA);
                if perm ~= 0
                    for dgt = dgts
                        new_results = [];
                        for i = 1:size(results, 1)
                            result = results(i, :);
                            for j = 1 : n_DNA - p_sites + 1
                                if all(result(j:j+p_sites-1) == 0)
                                    new_result = result(1, :);
                                    new_result(j:j+p_sites-1) = modes(dgt, :);
                                    if sum(dgts(dgts == dgt)) > 1 
                                        if isempty(new_results)
                                            new_results = [new_results; new_result];
                                        else
                                            if ~ismember(new_result, new_results, 'rows')
                                                new_results = [new_results; new_result];
                                            end
                                        end
                                    else
                                        new_results = [new_results; new_result];
                                    end
                                end
                            end
                        end
                        results = new_results;
                    end
                end
                f_results = [f_results; results];
            end
        end

        % Refer parameters
        function [ws, os] = refer_params(~, states, conc, e_l, c_l)
            global T;
            ws = zeros(1, size(states, 1));
            os = zeros(1, size(states, 1));
            for i = 1:size(states, 1)
                state = states(i, :);
                wi = zeros(1, length(state));
                oi = zeros(1, length(state));
                for j = 1:length(state)
                    if state(j) == 0
                        oi(j) = 0;
                        wi(j) = 1;
                    else
                        oi(j) = c_l(state(j));
                        wi(j) = exp(-e_l(state(j))/(0.001*8.315*T))*conc^(oi(j));
                    end
                end
                ws(i) = prod(wi);
                os(i) = sum(oi);
            end
        end

        % Get binding distribution
        function [ws, b_dis] = get_b_dis(~, fval, concs, states, no_of_modes, p_sites)
            global T;
            b_dis = zeros(length(concs), size(states, 2));
            ws = zeros(size(states, 1), size(states, 2), length(concs));
        
            if no_of_modes == 1
                be1 = fval(1);
            
                e_l = repelem(be1/p_sites, p_sites);
                c_l = repelem(1/p_sites, p_sites);
            else 
                be1 = fval(1);
                be2 = fval(2);
            
                e_l = [repelem(be1/p_sites, p_sites) repelem(be2/p_sites, p_sites)];
                c_l = [repelem(1/p_sites, p_sites) repelem(1/p_sites, p_sites)];
            end
        
            for c = 1:length(concs)
                conc = concs(c);
                for i = 1:size(states, 1)
                    state = states(i, :);
                    for j = 1:length(state)
                        if state(j) == 0
                            ws(i, j, c) = 1;
                        else
                            ws(i, j, c) = exp(-e_l(state(j))/(0.001*8.315*T))*conc^(c_l(state(j)));
                        end
                    end
                end
                temp = zeros(1, size(states, 2));
                for i = 1:size(states, 2)
                    temp(i) = sum(prod(ws(ws(:, i, c) ~= 1, :, c), 2));
                end
                Z = sum(prod(ws(:, :, c), 2));
                b_dis(c, :) = temp./Z;
            end
        end
    
        function [no_p, calc_S] = get_calc_sigs(bm, fval, concs, l_states, states, no_of_modes, p_sites)
            global s_min;
            calc_S = zeros(1, length(concs));
            no_p = zeros(1, length(concs));

            if no_of_modes == 1
                be1 = fval(1);
                s_max = fval(2);
            
                e_l = [be1 repelem(0, p_sites-1)];
                c_l = [1 repelem(0, p_sites-1)];
            else
                    be1 = fval(1);
                    be2 = fval(2);
                    s_max = fval(3);
                
                    e_l = [be1, repelem(0, p_sites-1),...
                           be2, repelem(0, p_sites-1)];
                    c_l = [1, repelem(0, p_sites-1),...
                           1, repelem(0, p_sites-1)];
            end
            
            for c = 1:length(concs)
                conc = concs(c);
                [wl, ol] = bm.refer_params(l_states, conc, e_l, c_l);
                [ws, ~] = bm.refer_params(states, conc, e_l, c_l);
        
                calc_S(c) = s_min + (s_max - s_min)*(sum(wl.*ol)/sum(ws));
                no_p(c) = sum(wl.*ol)/sum(ws);
            end
        end

        % Initial calculated signals
        function calc_S = init_calc_sigs(bm, vals, concs, no_of_modes, p_sites)
            if no_of_modes == 1
                fval = vals(1:2);
                n = vals(3);
                modes = vals(4);
                sites = vals(5);
            
                perms = bm.generate_perms(1:modes, floor(n/sites));
            else
                fval = vals(1:3);
                n = vals(4);
                modes = vals(5);
                sites = vals(6);
            
                perms = bm.generate_perms(1:modes, floor(n/sites));
                perms(perms == 11) = [];
            end
            DNA_states = bm.fill_perms(perms, n, modes, sites);
            [~, calc_S] = bm.get_calc_sigs(fval, concs, DNA_states, DNA_states, no_of_modes, p_sites);
        end
    end
end