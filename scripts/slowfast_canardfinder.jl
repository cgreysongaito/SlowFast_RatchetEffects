#### Functions for quasi-canard finder
#### "Slow organisms exhibit sudden population disappearances in a reddened world" by Greyson-Gaito, Gellner, & McCann.

function orientation(p1, p2, p3)
    val = (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p3[2] - p2[2]) * (p2[1] - p1[1])
    if (val == 0)
        return 0 #colinear
    elseif (val > 0)
        return 1 #clockwise
    else
        return -1 #anticlockwise
    end
end # Part 1 function from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

function dointersect(p1, p2, q1, q2) 
    o1 = orientation(p1, p2, q1)
    o2 = orientation(p1, p2, q2)
    o3 = orientation(q1, q2, p1)
    o4 = orientation(q1, q2, p2)

    if (o1 != o2 && o3 != o4)
        return true
    else
        return false
    end
end # Part 2 function from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

function cf_returnmap_check(sol, pass_points, rm_point1, rm_point2, which_pass)
    if which_pass ∈ ["first", "second"] == false
        error("which_pass must be either first or second")
    end
    rm_pass_points = []
    for j in 1:length(pass_points)
        for l in Int64(pass_points[j][1]):length(sol)-1
            if sol.u[l+1][1] < rm_point1[1] && sol.u[l][1] > rm_point1[1]
                if dointersect(sol.u[l],sol.u[l+1],rm_point1,rm_point2) == true #Stricter condition for intersection of return map
                    if which_pass == "first"
                        append!(rm_pass_points, [[l , sol.u[l][1], sol.u[l][2]]])
                    else
                        append!(rm_pass_points, [[l , sol.u[l][1], sol.u[l][2]]])
                        break
                    end
                end
            end
        end
    end
    if length(rm_pass_points) < 1
        return false
    else
        return [true, rm_pass_points]
    end
end #First (and last) check in canard finder (Steps 1 and 5 in SI). Whether trajectory intersects small line at the maximum of the resource isocline.

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
end #Second check in canard finder (Step 2 in SI). Whether trajectory enters thin box on the consumer axis from just above the intersection of the resource isocline and the consumer axis up to the maximum of the resource isocline

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
end #Third check in canard finder (Step 3 in SI). Whether trajectory enters thin box on the consumer axis from 0.0 up to specified point on consumer axis.

function res_iso(model, R, p)
    if model == "RozMac"
        @unpack a, k, r, h = p
        return r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
    else
        @unpack x, y, δ, R₀, k  = p
        return (R*(R*δ - R + R₀*δ - R₀) + k*(-R*δ + R - R₀*δ + R₀))/(k*x*y)
    end
end #function to help calculate the resource isocline for both the Rosenzweig-MacArthur and Yodzis-Innes models

function cf_ressiocline_check(model, sol, pass_points, res_Hopf_point, par) #TODO maybe do it so checks whether stays within bubble for certain period of time?
    new_pass_points = []
    for j in 1:length(pass_points)
        for l in Int64(pass_points[j][1]):length(sol)-5
            if sol.u[l][1] > res_Hopf_point && isapprox(sol.u[l][2], res_iso(model, sol.u[l][1], par); atol = 1e-1)
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
end #Fourth check in canard finder (Step 4 in SI). Whether trajectory crosses the resource isocline.

function axial_checker(sol)
    if isapprox(sol.u[end][1], 3.0; atol = 1e-2)  && isapprox(sol.u[end][2], 0.00; atol = 1e-2)
            return "axial"
        else
            return "nothing"
    end
end #Function to check whether trajectory lands on the axial solution or not

function hopfpointsfinder(model, effR0)
    if model == "RozMac"
        @unpack k, a, h, r, m = RozMacPar(e = effR0)
        hopf_R = k / 2 - 1 / (2 * a * h)
        hopf_C = r * (a * h * k * hopf_R - a * h * hopf_R^2 + k - hopf_R) / (a * k)
        zero_C = r * k  / (a * k)
    elseif model == "YodInn"
        @unpack k, R₀, δ, x, y = YodInnScalePar(R₀ = effR0)
        hopf_R = -R₀/2 + k/2
        hopf_C = (hopf_R*(hopf_R*δ - hopf_R + R₀*δ - R₀) + k*(-hopf_R*δ + hopf_R - R₀*δ + R₀))/(k*x*y)
        zero_C = (k*(- R₀*δ + R₀))/(k*x*y)
    else
        error("model must be either RozMac or YodInn")
    end
    return vcat(hopf_R, hopf_C, zero_C)
end #function to find the resource and consumer values of the maximum point on the resource isocline. Also returns the consumer value of the intersection of the resource isocline with the consumer axis.

function cf_returnmap(model, ep, effR0, freq, r, seed, tsend, tvals)
    if model == "RozMac"
        sol = RozMac_pert(ep, effR0, freq, r, seed, tsend, tvals)
        par = RozMacPar(e = effR0, ε = ep)
    elseif model == "YodInn"
        sol = YodInn_pert(ep, effR0, freq, r, seed, tsend, tvals)
        par = YodInnScalePar(R₀ = effR0, ε = ep)
    else
        error("model must be either RozMac or YodInn")
    end
    hopfpoints = hopfpointsfinder(model, effR0)
    res_Hopf_point = hopfpoints[1]
    rm_point1 = [hopfpoints[1], hopfpoints[2]-(hopfpoints[2]*0.1)] 
    rm_point2 = [hopfpoints[1], hopfpoints[2]+(hopfpoints[2]*0.1)]
    resaxisboxcheckpoints = [0.0, hopfpoints[3]*0.8]
    rm_pass_points1 = cf_returnmap_check(sol, [[1 , sol.u[1][1], sol.u[1][2]]], rm_point1, rm_point2, "first")
    if rm_pass_points1 == false
        return axial_checker(sol)
    else
        resaxis_pass_points = cf_resaxis_check(sol, rm_pass_points1[2], res_Hopf_point, [0.0, 0.1], [hopfpoints[3], hopfpoints[2]])
    end
    if resaxis_pass_points == false
        return axial_checker(sol)
    else
        resaxisbox_pass_points = cf_box_check(sol, resaxis_pass_points[2], [0.0,0.1], resaxisboxcheckpoints)
    end
    if resaxisbox_pass_points == false
        return axial_checker(sol)
    else
        resiso_pass_points = cf_ressiocline_check(model, sol, resaxisbox_pass_points[2], res_Hopf_point, par)
    end
    if resiso_pass_points == false
        return axial_checker(sol)
    else
        rm_pass_points2 = cf_returnmap_check(sol, resiso_pass_points[2],  rm_point1, rm_point2, "second")
    end
    if rm_pass_points2 == false
        return axial_checker(sol)
    else
        return "canard"
    end
end #Function to run through each checker and returns either "canard", "axial", or "nothing" for a single simulation of the Rosenzweig-MacArthur or Yodzis-Innes models.
