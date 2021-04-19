function mycross(a::Vector{T}, b::Vector{T}) where T<:Real
    r1 = a[2]*b[3] - a[3]*b[2]
    r2 = a[3]*b[1] - a[1]*b[3]
    r3 = a[1]*b[2] - a[2]*b[1]
    return [r1,r2,r3]
end

function mynorm(a::Vector{T}) where T<:Real
    n = sqrt(a[1]*a[1] + a[2]*a[2] + a[3]*a[3])
    return n
end

function get_circumcircle_info_triangle(nodeA, nodeB, nodeC, center, vecCircle1, vecCircle2, vecTemp)
    # center of circumcircle of triangle (with 3 points)
    vecCircle1 .= nodeB .- nodeA
    vecCircle2 .= nodeC .- nodeA
    vecTemp .= mycross(vecCircle1, vecCircle2) # norm vector
    l_ab = mynorm(vecCircle1)
    l_ac = mynorm(vecCircle2)
    l_n = mynorm(vecTemp)
    center .= ((mycross(vecTemp, vecCircle1) .* l_ac^2) .+ (mycross(vecCircle2,vecTemp) * l_ab^2)) ./ (2 * l_n^2)
    center .+= nodeA
    vecTemp ./= l_n
    # get two orthogonal vectors of circle
    vecCircle1 .= nodeA .- center
    # radius = mynorm(vecCircle1)
    vecCircle2 .= mycross(vecCircle1, vecTemp)
end

function get_circlepoints(mym::Mesh3D{T1,T2}) where {T1<:AbstractFloat, T2<:Integer}
    # get points on circumcircle of elements
    n_elements = size(mym.elements,1)
    step = 30 # point angle in deg
    n_points = 0
    n_points_per_element = convert(T2,360/step)
    # preallocations
    circleElements = zeros(T2, n_elements, n_points_per_element + 2)
    circleNodes = zeros(n_elements * (n_points_per_element + 2), 3)
    nodeA = Vector{T1}(undef,3)
    nodeB = Vector{T1}(undef,3)
    nodeC = Vector{T1}(undef,3)
    center = Vector{T1}(undef,3)
    vecCircle1 = Vector{T1}(undef,3)
    vecCircle2 = Vector{T1}(undef,3)
    vecTemp = Vector{T1}(undef,3)
    point = Vector{T1}(undef,3)
    for i = 1:n_elements
        p = mym.elementstatus[i,3]
        n1 = mym.nodes2parts[p,3]
        n2 = mym.nodes2parts[p,4]
        nodeset = @view mym.nodes[n1:n2,1:3]
        nodeA[:] .= nodeset[mym.elements[i,1],:]
        nodeB[:] .= nodeset[mym.elements[i,2],:]
        nodeC[:] .= nodeset[mym.elements[i,3],:]
        get_circumcircle_info_triangle(nodeA, nodeB, nodeC, center, vecCircle1, vecCircle2, vecTemp)
        for i_deg = step:step:360
            i_rad = i_deg * pi / 180
            point .= center .+ (vecCircle1 .* sin(i_rad)) .+ (vecCircle2 .* cos(i_rad))
            n_points += 1
            pointnumber = convert(T2, i_deg / step)
            circleElements[i, pointnumber] = n_points
            circleNodes[n_points,:] = point
        end
        # circumcircle is constructed based on nodeA
        # nodeA is therefore one of the 12 points
        # adding nodeB and nodeC manually to the collection
        n_points += 1
        circleElements[i, n_points_per_element + 1] = n_points
        circleNodes[n_points,:] = nodeB
        n_points += 1
        circleElements[i, n_points_per_element + 2] = n_points
        circleNodes[n_points,:] = nodeC
    end
    return Circles(circleNodes, circleElements)
end
