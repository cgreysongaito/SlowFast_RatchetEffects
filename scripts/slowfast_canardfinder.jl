# include("packages.jl")
#
# include("slowfast_commoncode.jl")
#
#
# #Maximum value of resource isocline
# 0.9318181818181819
# # Time series and phase plots of stochastic canards
# let
#     test = figure()
#     pert_timeseries_plot(0.01, 0.6, 0.0, 1, 1234, 5000.0, 2000.0:1.0:5000.0)
#     return test
# end
#
# let
#     test = figure()
#     pert_phase_plot(0.01, 0.6, 0.0, 1, 1234, 5000.0, 2000.0:1.0:5000.0)
#     vlines(0.9318181818181819, ymin = 2.1, ymax = 2.4, linestyles = `dashed`)
#     return test
# end
#
#
# #Examination of what sto solver spits out
# test = RozMac_pert(0.01, 0.6, 0.0, 1, 1234, 5000.0, 2000.0:1.0:5000.0)



function orientation(p1, p2, p3)
    val = (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p3[2] - p2[2]) * (p2[1] - p1[1])
    if (val == 0)
        return 0 #colinear
    elseif (val > 0)
        return 1 #clockwise
    else
        return -1 #anticlockwise
    end
end

function dointersect(p1, p2, q1, q2) # https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    o1 = orientation(p1, p2, q1)
    o2 = orientation(p1, p2, q2)
    o3 = orientation(q1, q2, p1)
    o4 = orientation(q1, q2, p2)

    if (o1 != o2 && o3 != o4)
        return true
    else
        return false
    end
end

function cf_returnmap_check(sol, pass_points, rm_point1, rm_point2)
    rm_pass_points = []
    for j in 1:length(pass_points)
        for l in Int64(pass_points[j][1]):length(sol)-1
            if sol.u[l+1][1] < rm_point1[1] && sol.u[l][1] > rm_point1[1]
                if dointersect(sol.u[l],sol.u[l+1],rm_point1,rm_point2) == true #Stricter condition for intersection of return map
                    append!(rm_pass_points, [[l , sol.u[l][1], sol.u[l][2]]])
                end
            end
        end
    end
    if length(rm_pass_points) < 1
        return false
    else
        return [true, rm_pass_points]
    end
end

function cf_resaxis_check(sol, pass_points, res_Hopf_point, res_lims, con_lims)
    resaxis_pass_points = []
    for j in 1:length(pass_points)
        for l in Int64(pass_points[j][1]):length(sol)-1
            if sol.u[l+1][1] < res_Hopf_point
                if res_lims[1] < sol.u[l][1] < res_lims[2] && con_lims[1] < sol.u[l][2] < con_lims[2]
                    append!(resaxis_pass_points, [[l , sol.u[l][1], sol.u[l][2]]])
                    break
                end
            else
                break
            end
        end
    end
    if length(resaxis_pass_points) < 1
        return false
    else
        return [true, resaxis_pass_points]
    end
end

function cf_box_check(sol, pass_points, res_lims, con_lims)
    new_pass_points = []
    for j in 1:length(pass_points)
        for l in Int64(pass_points[j][1]):length(sol)-1
            if res_lims[1] < sol.u[l][1] < res_lims[2] && con_lims[1] < sol.u[l][2] < con_lims[2]
                append!(new_pass_points, [[l , sol.u[l][1], sol.u[l][2]]])
                break
            end
        end
    end
    if length(new_pass_points) < 1
        return false
    else
        return [true, new_pass_points]
    end
end

function cf_ressiocline_check(sol, pass_points, res_Hopf_point, par) #TODO maybe do it so checks whether stays within bubble for certain period of time?
    new_pass_points = []
    for j in 1:length(pass_points)
        for l in Int64(pass_points[j][1]):length(sol)-5
            if sol.u[l][1] > res_Hopf_point && isapprox(sol.u[l][2], res_iso(sol.u[l][1], par); atol = 1e-1)
                append!(new_pass_points, [[l , sol.u[l][1], sol.u[l][2]]])
                break
            end
        end
    end

    if length(new_pass_points) < 1
        return false
    else
        return [true, new_pass_points]
    end
end

function cf_returnmap(ep, eff, freq, r, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    par = RozMacPar(e = eff, ε = ep)
    res_Hopf_point = 0.9318181818181819 # TODO make this general
    rm_point1 = [0.9318181818181819, 2.1] #NOTE THIS ONLY WORKS IF DON"T CHANGE a or k #TODO need to code in more general method - ie calculating max resisocline then adding error
    rm_point2 = [0.9318181818181819, 2.4]#NOTE THIS ONLY WORKS IF DON"T CHANGE a or k
    rm_pass_points1 = cf_returnmap_check(sol, [[1 , sol.u[1][1], sol.u[1][2]]], rm_point1, rm_point2)
    if rm_pass_points1 == false
        return false
    else
        resaxis_pass_points = cf_resaxis_check(sol, rm_pass_points1[2], res_Hopf_point, [0.0, 0.1], [2.1, 2.3])
    end
    if resaxis_pass_points == false
        return false
    else
        resaxisbox_pass_points = cf_box_check(sol, resaxis_pass_points[2], [0.0,0.1], [0.0,1.8])
    end
    if resaxisbox_pass_points == false
        return false
    else
        resiso_pass_points = cf_ressiocline_check(sol, resaxisbox_pass_points[2], res_Hopf_point, par)
    end
    if resiso_pass_points == false
        return false
    else
        rm_pass_points2 = cf_returnmap_check(sol, resiso_pass_points[2],  rm_point1, rm_point2)
    end
    if rm_pass_points2 == false
        return false
    else
        return true
    end
end
