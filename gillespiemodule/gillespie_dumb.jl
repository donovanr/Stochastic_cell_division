module gillespie_dumb
export gillespie_even_spacing_dumb


using Distributions


## return the index of a weighted draw, given a vector of weights
function weightedchoice(weights)
    normweights = weights/sum(weights)
    rnd = rand()
    for i = 1:length(normweights)
        rnd -= normweights[i]
        if rnd < 0.0
            return i
        end
    end
end


## return a random sample from a lognormal distribution
## with parameters mu and sigma and true mean and stdv m and s
function empirical_lognormal(m, s)
    mu = log( m^2/sqrt(m^2 + s^2) )
    sigma = sqrt( log(1 + s^2/m^2) )
    LogNormal(mu,sigma)
end


## gillespie for one dt
function gillespie_ssa_onestep(dt,state,stoich_M,react_propens,reaction_rates)
	t = 0.0
	while t < dt
		a = react_propens(state,reaction_rates)                              # update list of reaction propensities
		t += rand(Exponential(1/sum(a))) # take an exp distributed step
		if t < dt                                                    # if not too late,
            state += stoich_M[weightedchoice(a),:]                   # fire a reaction and update species vector
        end
    end
    return state
end


## this function randomly divides all the species counts in half binomially

function cell_divide(species_list)
    new_list = zeros(Int64,1,length(species_list))
    for i = 1:length(species_list)
        if species_list[i] > 0
            new_list[i] = rand(Binomial(species_list[i]))
        end
    end
    return new_list
end


## provide some feedback on divergent growth/death of species in cell
function cell_regulate(spec_list,init_list,reg_factor)
    thresh_factor = float(reg_factor)
    heal_factor = float(1/reg_factor)
    for j = 1:length(spec_list) # don't skip length variable
        if spec_list[j] < init_list[j]/thresh_factor
            spec_list[j] = iround(spec_list[j]*(1+heal_factor)) + 1
        end
        if spec_list[j] > init_list[j]*thresh_factor
            spec_list[j] = iround(spec_list[j]/(1+heal_factor))
        end
    end
    return spec_list
end

# function cell_regulate_feedback(spec_list,init_list,reg_factor,thresh_factor,reaction_rates,orignial_rates)
#     for j = 1:length(spec_list) # don't skip length variable
#         if (spec_list[j] > init_list[j]/thresh_factor) & (spec_list[j] < init_list[j]*thresh_factor)
#             reaction_rates[j] = orignial_rates[j]
#         elseif spec_list[j] <= init_list[j]/thresh_factor
#             reaction_rates[j] = reaction_rates[j]*(1.0 + reg_factor)
#         elseif spec_list[j] >= init_list[j]*thresh_factor
#             reaction_rates[j] = reaction_rates[j]*(1.0 - reg_factor)
#         end
#     end
#     return reaction_rates
# end## feedback on growth/death as a modification of the growth rate


function noisey_weights(reaction_rates, noisefactor)
    for j = 1:length(reaction_rates) # don't skip length variable
        reaction_rates[j] = rand(empirical_lognormal(reaction_rates[j],noisefactor*reaction_rates[j]))
        return reaction_rates
    end
end


## if a species drops to zero, set it to one
function cell_resuscitate(spec_list)
    for j = 1:length(spec_list) # check if length variable
        if spec_list[j] < 1
            spec_list[j] = 1
        end
    end
    return spec_list
end


## do n steps of gillespie_one_step, updating and recording

function gillespie_even_spacing_dumb(num_steps, delta_t, initial_state, stoich_Mat, react_props, reaction_rates, mean_divsize, sigma_divsize, regfac, ratenoisefactor)
    orignial_rates = copy(reaction_rates)
    z_list = zeros(num_steps,length(initial_state))
    t_list = zeros(num_steps,1)
    r_list = zeros(num_steps,length(reaction_rates))
    divsize = rand(empirical_lognormal(mean_divsize,sigma_divsize))
    z = copy(initial_state)
    for i = 1:num_steps
        z_list[i,:] = z
        t_list[i] = (i-1)*delta_t
        r_list[i,:] = reaction_rates
        z = gillespie_ssa_onestep(delta_t,z,stoich_Mat,react_props,reaction_rates)
        if z[1] > divsize
            z = cell_divide(z)
            z = cell_regulate(z, initial_state, regfac)
            #reaction_rates = cell_regulate_feedback(z,initial_state,regfac,thresh_factor,reaction_rates,orignial_rates)
            reaction_rates = noisey_weights(orignial_rates, ratenoisefactor)
            z = cell_resuscitate(z)
            divsize = rand(empirical_lognormal(mean_divsize,sigma_divsize))
        end
    end
    data = hcat(t_list,z_list,r_list)
    return data
end

end # module