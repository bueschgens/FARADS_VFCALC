function buckets_main(mym::Mesh3D, delta_XYZ::Vector{T1}; blockparts = (1:size(mym.elements2parts,1))) where T1<:AbstractFloat
    println("starting bucket creation and sorting")
    # create empty buckets
    b_empty = create_buckets(mym, delta_XYZ)
    println("    created ", size(b_empty.volumes,1), " buckets")
    # circle around elements
    circles = get_circlepoints(mym)
    # occupation of buckets
    # blockparts determines which parts participate in blocking
    tol = 0 # tolerance for bucket coords
    b2p = check_bucket_and_parts(b_empty, mym, tol, blockparts)
    b2e, bstatus = check_bucket_and_elements(b_empty, mym, b2p, circles, tol)
    println("bucket creation and sorting done")
    return OccBuckets(b_empty.nodes, b_empty.volumes, bstatus, b2p, b2e, b_empty.delta, b_empty.n_buckets_dir)
end

function buckets_main_parallel(mym::Mesh3D, delta_XYZ::Vector{T1}; blockparts = (1:size(mym.elements2parts,1))) where T1<:AbstractFloat
    println("starting bucket creation and sorting")
    # create empty buckets
    b_empty = create_buckets(mym, delta_XYZ)
    println("    created ", size(b_empty.volumes,1), " buckets")
    # circle around elements
    circles = get_circlepoints(mym)
    # occupation of buckets
    # blockparts determines which parts participate in blocking
    tol = 0 # tolerance for bucket coords
    b2p = check_bucket_and_parts(b_empty, mym, tol, blockparts)
    b2e, bstatus = check_bucket_and_elements_parallel(b_empty, mym, b2p, circles, tol)
    println("bucket creation and sorting done")
    return OccBuckets(b_empty.nodes, b_empty.volumes, bstatus, b2p, b2e, b_empty.delta, b_empty.n_buckets_dir)
end

function create_buckets(mym::Mesh3D{T1,T2}, delta_XYZ::Vector{T1}) where {T1<:AbstractFloat, T2<:Integer}
    # get all necessary data based on Mesh3D to prepare bucket creation
    # get min and max dimensions of Mesh3D
    min_XYZ = vec(minimum(mym.nodes, dims = 1))
    max_XYZ = vec(maximum(mym.nodes, dims = 1))
    # get number of buckets per dimension
    length_XYZ = [(max_XYZ[i] - min_XYZ[i]) for i in [1,2,3]]
    n_buckets_XYZ = [ceil(T2, length_XYZ[i] / delta_XYZ[i]) for i in [1,2,3]]
    # directly set at right start position
    x = 0 + min_XYZ[1] : delta_XYZ[1] : (n_buckets_XYZ[1] * delta_XYZ[1]) + min_XYZ[1]
    y = 0 + min_XYZ[2] : delta_XYZ[2] : (n_buckets_XYZ[2] * delta_XYZ[2]) + min_XYZ[2]
    z = 0 + min_XYZ[3] : delta_XYZ[3] : (n_buckets_XYZ[3] * delta_XYZ[3]) + min_XYZ[3]
    n_nodes_X = size(x,1)
    n_nodes_Y = size(y,1)
    n_nodes_Z = size(z,1)
    # create nodes in plane at z = 0
    n_nodes_XY = 0
    nodes_XY = zeros(n_nodes_X * n_nodes_Y, 2)
    nodes_quad = zeros(T2, n_nodes_X, n_nodes_Y)
    for iy = 1:n_nodes_Y
        for ix = 1:n_nodes_X
            n_nodes_XY += 1
            nodes_XY[n_nodes_XY,1] = x[ix]
            nodes_XY[n_nodes_XY,2] = y[iy]
            nodes_quad[ix,iy] = n_nodes_XY
        end
    end
    # create elements in the plane with nodes
    n_elements_XY = 0 # number of elements per plane
    elements_XY = zeros(T2, (n_nodes_X-1) * (n_nodes_Y-1), 4)
    elements_quad = zeros(T2, (n_nodes_X-1), (n_nodes_Y-1))
    for iy = 1:(n_nodes_Y - 1)
        for ix = 1:(n_nodes_X - 1)
            n_elements_XY += 1
            elements_XY[n_elements_XY,1] = nodes_quad[ix,iy]
            elements_XY[n_elements_XY,2] = nodes_quad[ix,iy+1]
            elements_XY[n_elements_XY,3] = nodes_quad[ix+1,iy+1]
            elements_XY[n_elements_XY,4] = nodes_quad[ix+1,iy]
            elements_quad[ix,iy] = n_elements_XY
        end
    end
    # extrude XY node plane in z direction
    nodes_XYZ = zeros(n_nodes_XY * n_nodes_Z, 3)
    elements_XY_Z = zeros(T2, n_elements_XY * n_nodes_Z, 4)
    # all elements of xy plane in all z positions
    plane = 0
    for iz = 1:n_nodes_Z
        n1 = (plane * n_nodes_XY) + 1
        n2 = (plane + 1) * n_nodes_XY
        nodes_XYZ[n1:n2,1:2] = nodes_XY[:,:]
        nodes_XYZ[n1:n2,3] .= z[iz]
        e1 = (plane * n_elements_XY) + 1
        e2 = (plane + 1) * n_elements_XY
        elements_XY_Z[e1:e2,:] = elements_XY[:,:] .+ (plane * n_nodes_XY)
        plane += 1
    end
    # create buckets / volumes
    buckets = zeros(T2, size(elements_XY,1) * (n_nodes_Z - 1),8)
    plane = 0
    for i = 1:(n_nodes_Z - 1) # number planes - 1
        b1 = (plane * n_elements_XY) + 1
        b2 = (plane + 1) * n_elements_XY
        unten1 = (plane * n_elements_XY) + 1
        unten2 = (plane + 1) * n_elements_XY
        oben1 = ((plane + 1) * n_elements_XY) + 1
        oben2 = (plane + 2) * n_elements_XY
        buckets[b1:b2,1:4] = elements_XY_Z[unten1:unten2,:]
        buckets[b1:b2,5:8] = elements_XY_Z[oben1:oben2,:]
        plane += 1
    end
    return EmptyBuckets(nodes_XYZ, buckets, delta_XYZ, n_buckets_XYZ)
end

function check_bucket_and_parts(myb::EmptyBuckets, mym::Mesh3D{T1,T2}, tolerance, parts) where {T1<:AbstractFloat, T2<:Integer}
    # check buckets for intersection with all selected parts
    # parts is defaulted one function higher
    n_buckets = size(myb.volumes,1)
    n_parts = size(mym.elements2parts,1)
    buckets2parts = Vector{Vector{T2}}(undef,n_buckets)
    for i_b = 1:n_buckets
        # get min and max coords of bucket
        min_bucket_XYZ = minimum(myb.nodes[myb.volumes[i_b,:],1:3], dims = 1)
        max_bucket_XYZ = maximum(myb.nodes[myb.volumes[i_b,:],1:3], dims = 1)
        # small tolerance with bucket coords
        min_bucket_XYZ .-= tolerance
        max_bucket_XYZ .+= tolerance
        temp_b2p = []
        for i_p in parts
            # get min and max coords of part
            n1 = mym.nodes2parts[i_p,3]
            n2 = mym.nodes2parts[i_p,4]
            min_part_XYZ = minimum(mym.nodes[n1:n2,1:3], dims = 1)
            max_part_XYZ = maximum(mym.nodes[n1:n2,1:3], dims = 1)
            hit = 1
            i_xyz = 1
            while hit == 1 && i_xyz <= 3
                if ((min_bucket_XYZ[i_xyz] >= min_part_XYZ[i_xyz] &&
                     min_bucket_XYZ[i_xyz] <= max_part_XYZ[i_xyz]) || 
                    (max_bucket_XYZ[i_xyz] >= min_part_XYZ[i_xyz] && 
                     max_bucket_XYZ[i_xyz] <= max_part_XYZ[i_xyz]) || 
                    (min_bucket_XYZ[i_xyz] <= min_part_XYZ[i_xyz] && 
                     max_bucket_XYZ[i_xyz] >= max_part_XYZ[i_xyz]))
                    # for all 3 dims (x,y,z)
                    # check if min of bucket is between min and max of part
                    # check if max of bucket is between min and max of part
                    # check if bucket min and max are both outside of part min and max
                    hit = 1
                else 
                    hit = 0
                end
                i_xyz += 1
            end
            if hit == 1
                push!(temp_b2p, i_p)
            end
        end
        buckets2parts[i_b] = temp_b2p
    end
    return buckets2parts
end

function check_bucket_and_elements(myb::EmptyBuckets, mym::Mesh3D{T1,T2}, buckets2parts, circles, tolerance) where {T1<:AbstractFloat, T2<:Integer}
    # check if elements of penetrating parts penetrate bucket
    n_buckets = size(myb.volumes,1)
    n_parts = size(mym.elements2parts,1)
    n_points_per_element = size(circles.elements,2)
    buckets2elements = Vector{Vector{T2}}(undef,n_buckets)
    bucketstatus = zeros(T2, n_buckets) # 0(empty) or 1(occupied)
    temp_b2e = Vector{T2}(undef,500)
    checks_tot = 0
    checks_com = 0
    checks_not_com = 0
    checks_near = 0
    checks_circle = 0
    progress = Progress(n_buckets, dt=1, barglyphs=BarGlyphs("[|| ]"), barlen=50)
    for i_b = 1:n_buckets
        next!(progress)
        # get min and max coords of bucket
        min_bucket_XYZ = minimum(myb.nodes[myb.volumes[i_b,:],1:3], dims = 1)
        max_bucket_XYZ = maximum(myb.nodes[myb.volumes[i_b,:],1:3], dims = 1)
        # small tolerance with bucket coords
        min_bucket_XYZ .-= tolerance
        max_bucket_XYZ .+= tolerance
        # prepare for circle check
        b_com = ((max_bucket_XYZ[:] .+ min_bucket_XYZ[:]) ./ 2)
        b_n = @view myb.nodes[myb.volumes[i_b,1],:]
        # go through all parts penetrating the bucket
        hit = 0
        for i_p in buckets2parts[i_b]
            n1 = mym.nodes2parts[i_p,3]
            n2 = mym.nodes2parts[i_p,4]
            nodes_of_part = @view mym.nodes[n1:n2,1:3]
            e1 = mym.elements2parts[i_p,3]
            e2 = mym.elements2parts[i_p,4]
            # go through all elements of part
            for i_e = e1:e2
                checks_tot += 1
                # first check com
                e_com = @view mym.com[i_e,:]
                if is_point_inside_quader(e_com, min_bucket_XYZ, max_bucket_XYZ)
                    checks_com += 1
                    # com inside bucket
                    hit += 1
                    temp_b2e[hit] = i_e
                else # check if element is near bucket
                    checks_not_com += 1
                    # presorting of the parts elements
                    e_nA = @view nodes_of_part[mym.elements[i_e,1],:]
                    e_nB = @view nodes_of_part[mym.elements[i_e,2],:]
                    e_nC = @view nodes_of_part[mym.elements[i_e,3],:]
                    if is_element_near_bucket(e_com, e_nA, e_nB, e_nC, b_com, b_n)
                        # if near check all points of circle
                        checks_near += 1
                        # TO DO: create cirle around element in situ
                        for i_points = 1:n_points_per_element
                            point = @view circles.nodes[circles.elements[i_e,i_points],:]
                            if is_point_inside_quader(point, min_bucket_XYZ, max_bucket_XYZ)
                                # circle point inside bucket
                                checks_circle += 1
                                hit += 1
                                temp_b2e[hit] = i_e
                                break
                            end
                        end
                    end
                end
            end
        end
        if hit > 0
            bucketstatus[i_b,1] = 1
            buckets2elements[i_b] = temp_b2e[1:hit]
        else
            buckets2elements[i_b] = []
        end
    end
    println("Result checking buckets and elements")
    println("    checks total: ",checks_tot)
    println("    com inside bucket: ",checks_com, " / others: ",checks_not_com)
    println("    others --> near bucket: ", checks_near, " --> circle inside: ", checks_circle)
    return buckets2elements, bucketstatus
end

function check_bucket_and_elements_parallel(myb::EmptyBuckets, mym::Mesh3D{T1,T2}, buckets2parts, circles, tolerance) where {T1<:AbstractFloat, T2<:Integer}
    # check if elements of penetrating parts penetrate bucket
    n_buckets = size(myb.volumes,1)
    n_parts = size(mym.elements2parts,1)
    n_points_per_element = size(circles.elements,2)
    buckets2elements = Vector{Vector{T2}}(undef,n_buckets)
    bucketstatus = zeros(T2, n_buckets) # 0(empty) or 1(occupied)
    n_threads = Threads.nthreads()
    println("starting blocking check with ",n_threads," threads")
    temp_b2e = Array{T2,2}(undef,500,n_threads)
    checks_tot = zeros(T2,n_threads)
    checks_com = zeros(T2,n_threads)
    checks_not_com = zeros(T2,n_threads)
    checks_near = zeros(T2,n_threads)
    checks_circle = zeros(T2,n_threads)
    progress = Progress(n_buckets, dt=1, barglyphs=BarGlyphs("[|| ]"), barlen=50)
    for i_b = 1:n_buckets
        next!(progress)
        # get min and max coords of bucket
        min_bucket_XYZ = minimum(myb.nodes[myb.volumes[i_b,:],1:3], dims = 1)
        max_bucket_XYZ = maximum(myb.nodes[myb.volumes[i_b,:],1:3], dims = 1)
        # small tolerance with bucket coords
        min_bucket_XYZ .-= tolerance
        max_bucket_XYZ .+= tolerance
        # prepare for circle check
        b_com = ((max_bucket_XYZ[:] .+ min_bucket_XYZ[:]) ./ 2)
        b_n = @view myb.nodes[myb.volumes[i_b,1],:]
        # go through all parts penetrating the bucket
        hit = zeros(T2,n_threads)
        for i_p in buckets2parts[i_b]
            n1 = mym.nodes2parts[i_p,3]
            n2 = mym.nodes2parts[i_p,4]
            nodes_of_part = @view mym.nodes[n1:n2,1:3]
            e1 = mym.elements2parts[i_p,3]
            e2 = mym.elements2parts[i_p,4]
            # go through all elements of part
            # Threads.@threads for i_e = e1:e2
            Threads.@threads for i_e = shuffle!(collect(e1:e2))
                checks_tot[Threads.threadid()] += 1
                # first check com
                e_com = @view mym.com[i_e,:]
                if is_point_inside_quader(e_com, min_bucket_XYZ, max_bucket_XYZ)
                    checks_com[Threads.threadid()] += 1
                    # com inside bucket
                    hit[Threads.threadid()] += 1
                    temp_b2e[hit[Threads.threadid()],Threads.threadid()] = i_e
                else # check if element is near bucket
                    checks_not_com[Threads.threadid()] += 1
                    # presorting of the parts elements
                    e_nA = @view nodes_of_part[mym.elements[i_e,1],:]
                    e_nB = @view nodes_of_part[mym.elements[i_e,2],:]
                    e_nC = @view nodes_of_part[mym.elements[i_e,3],:]
                    if is_element_near_bucket(e_com, e_nA, e_nB, e_nC, b_com, b_n)
                        # if near check all points of circle
                        checks_near[Threads.threadid()] += 1
                        # TO DO: create cirle around element in situ
                        for i_points = 1:n_points_per_element
                            point = @view circles.nodes[circles.elements[i_e,i_points],:]
                            if is_point_inside_quader(point, min_bucket_XYZ, max_bucket_XYZ)
                                # circle point inside bucket
                                checks_circle[Threads.threadid()] += 1
                                hit[Threads.threadid()] += 1
                                temp_b2e[hit[Threads.threadid()],Threads.threadid()] = i_e
                                break
                            end
                        end
                    end
                end
            end
        end
        if sum(hit) > 0
            bucketstatus[i_b,1] = 1
            temp_b2e_all = Vector{T2}(undef,sum(hit))
            i1 = 0
            i2 = 0
            for i = 1:n_threads
                i1 = i2 + 1
                i2 = i2 + hit[i]
                temp_b2e_all[i1:i2] = temp_b2e[1:hit[i],i]
            end
            buckets2elements[i_b] = temp_b2e_all
        else
            buckets2elements[i_b] = []
        end
    end
    println("Result checking buckets and elements")
    println("    checks total: ",sum(checks_tot))
    println("    com inside bucket: ",sum(checks_com), " / others: ",sum(checks_not_com))
    println("    others --> near bucket: ", sum(checks_near), " --> circle inside: ", sum(checks_circle))
    return buckets2elements, bucketstatus
end

function is_point_inside_quader(point, quader_min, quader_max)
    # quader_min und quader_max are min und max coords of bucket(X,Y,Z)
    isinside = false
    if (point[1] >= quader_min[1] && 
        point[1] <= quader_max[1] && 
        point[2] >= quader_min[2] && 
        point[2] <= quader_max[2] && 
        point[3] >= quader_min[3] && 
        point[3] <= quader_max[3])
        isinside = true
    end
    return isinside
end

function is_element_near_bucket(e_com, e_nA, e_nB, e_nC, b_com, b_n)
    # check if element is near bucket 
    # with distances center to corners
    vec_l(a,b) = mynorm(a-b)
    b_l = vec_l(b_com, b_n)
    e_l = [vec_l(e_nA, e_com), vec_l(e_nB, e_com), vec_l(e_nC, e_com)]
    con_l = vec_l(e_com,b_com)
    isnear = false
    if (con_l - maximum(e_l)) <= (b_l - 1E-8)
        isnear = true
    end
    return isnear
end

function get_elements_inside_bucket(mym::Mesh3D, myb::OccBuckets, bucket)
    # print all elements in the specific bucket
    println("All Elements inside Bucket ",bucket)
    b_elems = Vector{Int64}()
    for i_e in myb.buckets2elements[bucket]
        e_p = mym.elementstatus[i_e,3] # returns part number
        n1 = mym.nodes2parts[e_p,3]
        n2 = mym.nodes2parts[e_p,4]
        nodes = mym.nodes[n1:n2,1:3]
        elements = zeros(Int,1,3)
        elements[1,:] = mym.elements[i_e,:]
        push!(b_elems, i_e)
    end
    println(b_elems)
end

function get_number_of_elements_in_buckets(myb::OccBuckets)
    # get number of elements in the specific bucket
    n_buckets = size(myb.volumes,1)
    blist = zeros(Int64, n_buckets)
    for i_b = 1:n_buckets
        n_e_in_b = size(myb.buckets2elements[i_b],1)
        blist[i_b] = n_e_in_b
        # println(n_e_in_b)
    end
    return blist
end

function info_bucket_occupation(myb::OccBuckets)
    # get information on bucket occupation
    # max min and mean
    belemlist = get_number_of_elements_in_buckets(myb)
    bmax = maximum(belemlist)
    bmin = minimum(belemlist)
    belemlist_occ = belemlist[myb.occupied[:] .== 1]
    bmin_occ = minimum(belemlist_occ)
    n_buckets = size(myb.volumes,1)
    n_buckets_occ = sum(myb.occupied)
    bmean = sum(belemlist) / n_buckets
    bmean = round(bmean, digits=1)
    bmean_occ = sum(belemlist) / n_buckets_occ
    bmean_occ = round(bmean_occ, digits=1)
    println("Bucket occupation with elements:")
    println("    Buckets occupied: ", n_buckets_occ, " / ", n_buckets)
    println("    Max: ", bmax, " // Min: ", bmin, " // Min_occ: ", bmin_occ)
    println("    Mean: ", bmean, " // Mean_occ: ", bmean_occ)
end

function bucketcellarray2cellvec(cellArr::Vector{Vector{T2}}) where T2<:Integer
	# conversion from cell array to vector+assignment with empty cells
	n_cells = size(cellArr,1)
	cellvec = zeros(T2, 0) # append!(cellvec,current_cell)
	cellvec_assignment = zeros(T2, n_cells, 3) # first is empty or not
	for i = 1:n_cells
		size_old = size(cellvec,1)
		current_cell = cellArr[i]
		if size(current_cell,1) == 0 # isempty
			cellvec_assignment[i,1] = 0
		else
			append!(cellvec,current_cell)
			size_new = size(cellvec,1)
			cellvec_assignment[i,1] = 1
			cellvec_assignment[i,2] = size_old + 1
			cellvec_assignment[i,3] = size_new
		end
	end
	return cellvec, cellvec_assignment
end

function make_buckets_vector(myb::OccBuckets{T1, T2}) where {T1<:AbstractFloat, T2<:Integer}
	a = vec(myb.nodes)
	b = vec(myb.volumes)
	c = vec(myb.occupied)
	bucketcells, bucketcellassign = bucketcellarray2cellvec(myb.buckets2elements)
	d = vec(bucketcells)
	e = vec(bucketcellassign)
	f = vec(myb.delta)
	g = vec(myb.n_XYZ)
	h = size(myb.occupied,1)
	return VecOccBuckets(a, b, c, d, e, f, g, h)
end