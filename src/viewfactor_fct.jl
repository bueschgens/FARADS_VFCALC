function calculate_vf!(vfmat, refinement, mym::VecMesh3D{T1,T2}) where {T1<:AbstractFloat, T2<:Integer}
	# run the viewfactor calculation
	n_elements = mym.nelements
	if refinement <= -4 || refinement >= 6
		error("Refinement modus not supported")
	end
    if refinement < 0
        autorefine = 1
        autoscale = refinement
        ref_power = 5
    else
        autorefine = 0
        ref_power = refinement
    end
    println("refinement set to ",refinement)
    n_threads = Threads.nthreads()
    println("starting viewfactor calculation with ",n_threads," threads")
    # prealloc ref vf mat for each thread
    n_vf = 4 ^ ref_power
    vftemp = Array{Array{T1,2},2}(undef,n_threads,2)
    for i = 1:n_threads
        vftemp[i,1] = zeros(n_vf,n_vf)
        vftemp[i,2] = zeros(n_vf,n_vf)
    end
    progress = Progress(n_elements, dt=1, barglyphs=BarGlyphs("[|| ]"), barlen=50)
    for i1 = 1:1:n_elements
        next!(progress)
        #for i2 = (i1+1):1:n_elements
        Threads.@threads for i2 = (i1+1):1:n_elements
        #Threads.@threads for i2 = shuffle!(collect((i1+1):n_elements))
            if vfmat[i1,i2] > 0 # --> vf existing
                phi1, phi2, length_s = get_vfcalc_data(i1, i2, mym)
                if autorefine == 1
                    pair_refine = checkRefinementModus(i1, i2, mym, length_s, autoscale)
                else
                    pair_refine = refinement
                end
                # vfmat[i1, i2] = pair_refine
                # vfmat[i2, i1] = pair_refine
                if pair_refine == 0 # --> without refinement
                    vfmat[i1, i2] = phi1 * phi2 * mym.area[i2] / (pi * length_s * length_s)
                    vfmat[i2, i1] = phi1 * phi2 * mym.area[i1] / (pi * length_s * length_s)
                else # --> refinementModus
                    vf_refined1, vf_refined2 = doRefinement(i1, i2, mym, pair_refine, vftemp)
                    vfmat[i1, i2] = vf_refined1
                    vfmat[i2, i1] = vf_refined2
                end
            end
        end
    end 
    println("viewfactor calculation done")
end

@inline function get_vfcalc_data(i1, i2, mym)
    # calc vf data between two elements
    nvec1 = @SVector [mym.nvec[i1], mym.nvec[i1+1*mym.nelements], mym.nvec[i1+2*mym.nelements]]
    nvec2 = @SVector [mym.nvec[i2], mym.nvec[i2+1*mym.nelements], mym.nvec[i2+2*mym.nelements]]
    com1 = @SVector [mym.com[i1], mym.com[i1+1*mym.nelements], mym.com[i1+2*mym.nelements]]
    com2 = @SVector [mym.com[i2], mym.com[i2+1*mym.nelements], mym.com[i2+2*mym.nelements]]
    phi1, phi2, length_s = calcVFmath(nvec1, nvec2, com1, com2)
    return phi1, phi2, length_s
end

function checkRefinementModus(i1, i2, mym, length_s, autoscale)
	# get size of element
	size1 = getElemSize(i1, mym)
	size2 = getElemSize(i2, mym)
    size_max = max(size1, size2)

	check = length_s / size_max

	if autoscale == -3
		check = min(check, 7.4) # beyond 7.5 negativ
		modus_raw = (-1.2) * check + 9
    elseif autoscale == -2
		check = min(check, 7.4) # beyond 7.5 negativ
		modus_raw = (-0.8) * check + 6
	elseif autoscale == -1
		check = min(check, 6.2) # beyond 6.25 negativ
        modus_raw = (-0.8) * check + 5
    else
        error("Refinement: autoscale out of range")
    end

    modus = floor(Int64, modus_raw)

	if (modus > 5)
        return 5 # max refinement modus allowed
    else
        return modus
    end
end

function getElemSize(i_elem, mym)
    # get information on triangle
    part = mym.elementstatus[i_elem+2*mym.nelements] # part number
    nodes1 = mym.elements[i_elem]
	nodes2 = mym.elements[i_elem+1*mym.nelements]
	nodes3 = mym.elements[i_elem+2*mym.nelements]
    nodeoffset = mym.nodes2parts[part+2*mym.nparts]
    
    a1 = mym.nodes[nodeoffset - 1 + nodes1]
    a2 = mym.nodes[nodeoffset - 1 + nodes1+1*mym.nnodes]
    a3 = mym.nodes[nodeoffset - 1 + nodes1+2*mym.nnodes]
    b1 = mym.nodes[nodeoffset - 1 + nodes2]
    b2 = mym.nodes[nodeoffset - 1 + nodes2+1*mym.nnodes]
    b3 = mym.nodes[nodeoffset - 1 + nodes2+2*mym.nnodes]
    c1 = mym.nodes[nodeoffset - 1 + nodes3]
    c2 = mym.nodes[nodeoffset - 1 + nodes3+1*mym.nnodes]
    c3 = mym.nodes[nodeoffset - 1 + nodes3+2*mym.nnodes]

    l_ab = sqrt((b1-a1) * (b1-a1) + (b2-a2) * (b2-a2) + (b3-a3) * (b3-a3))
	l_bc = sqrt((c1-b1) * (c1-b1) + (c2-b2) * (c2-b2) + (c3-b3) * (c3-b3))
    l_ca = sqrt((a1-c1) * (a1-c1) + (a2-c2) * (a2-c2) + (a3-c3) * (a3-c3))

	l_max = max(l_ab, l_bc)
	l_max = max(l_max, l_ca)
	l_max = l_max / 2.0

	return l_max
end


function doRefinement(i1, i2, mym, modus, vftemp)
	# original element o1 o2
	# first refined i1 i2

	coords_o1, nvec_o1 = getSpecificTriangleElem(i1, mym)
    coords_o2, nvec_o2 = getSpecificTriangleElem(i2, mym)

    n_vf = 4 ^ modus
    
    if modus == 1
        coords_i1, com_i1, nvec_i1, area_i1 = getRefinedTriangles(coords_o1, nvec_o1)
        coords_i2, com_i2, nvec_i2, area_i2 = getRefinedTriangles(coords_o2, nvec_o2)
    elseif modus == 2
        com_i1, nvec_i1, area_i1 = getRefinedTriangles2(coords_o1, nvec_o1, n_vf)
        com_i2, nvec_i2, area_i2 = getRefinedTriangles2(coords_o2, nvec_o2, n_vf)
    elseif modus == 3
        com_i1, nvec_i1, area_i1 = getRefinedTriangles3(coords_o1, nvec_o1, n_vf)
        com_i2, nvec_i2, area_i2 = getRefinedTriangles3(coords_o2, nvec_o2, n_vf)
    elseif modus == 4
        com_i1, nvec_i1, area_i1 = getRefinedTriangles4(coords_o1, nvec_o1, n_vf)
        com_i2, nvec_i2, area_i2 = getRefinedTriangles4(coords_o2, nvec_o2, n_vf)
    elseif modus == 5
        com_i1, nvec_i1, area_i1 = getRefinedTriangles5(coords_o1, nvec_o1, n_vf)
        com_i2, nvec_i2, area_i2 = getRefinedTriangles5(coords_o2, nvec_o2, n_vf)
    else
        error("refinement modus not implemented: ", modus, " between ", i1, " and ", i2)
    end

    vf_refined1, vf_refined2 = getVFfromRefinement(n_vf, com_i1, nvec_i1, area_i1, com_i2, nvec_i2, area_i2, vftemp)

    return vf_refined1, vf_refined2
end

function getSpecificTriangleElem(i_e, mym)
    # get information on original triangle
    part = mym.elementstatus[i_e+2*mym.nelements] # part number
    nodes1 = mym.elements[i_e]
	nodes2 = mym.elements[i_e+1*mym.nelements]
	nodes3 = mym.elements[i_e+2*mym.nelements]
    nodeoffset = mym.nodes2parts[part+2*mym.nparts]
    
    coords = SMatrix{3,3}(mym.nodes[nodeoffset - 1 + nodes1],
                          mym.nodes[nodeoffset - 1 + nodes2],
                          mym.nodes[nodeoffset - 1 + nodes3],
                          mym.nodes[nodeoffset - 1 + nodes1+1*mym.nnodes],
                          mym.nodes[nodeoffset - 1 + nodes2+1*mym.nnodes],
                          mym.nodes[nodeoffset - 1 + nodes3+1*mym.nnodes],
                          mym.nodes[nodeoffset - 1 + nodes1+2*mym.nnodes],
                          mym.nodes[nodeoffset - 1 + nodes2+2*mym.nnodes],
                          mym.nodes[nodeoffset - 1 + nodes3+2*mym.nnodes])

    nvec = @SVector [mym.nvec[i_e], mym.nvec[i_e+1*mym.nelements], mym.nvec[i_e+2*mym.nelements]]

    return coords, nvec
end

function getRefinedTriangles(coords_o, nvec_o)
	# new coords
    middle = SMatrix{3,3}((coords_o[2,1] + coords_o[1,1])/2, 
                          (coords_o[3,1] + coords_o[2,1])/2,
                          (coords_o[1,1] + coords_o[3,1])/2,
                          (coords_o[2,2] + coords_o[1,2])/2, 
                          (coords_o[3,2] + coords_o[2,2])/2,
                          (coords_o[1,2] + coords_o[3,2])/2,
                          (coords_o[2,3] + coords_o[1,3])/2, 
                          (coords_o[3,3] + coords_o[2,3])/2,
                          (coords_o[1,3] + coords_o[3,3])/2)

    coords_n = SMatrix{6,3}(coords_o[1,1], coords_o[2,1], coords_o[3,1],
                            middle[1,1], middle[2,1], middle[3,1],
                            coords_o[1,2], coords_o[2,2], coords_o[3,2],
                            middle[1,2], middle[2,2], middle[3,2],
                            coords_o[1,3], coords_o[2,3], coords_o[3,3],
                            middle[1,3], middle[2,3], middle[3,3])

    # new nodes
    nodes_n = SMatrix{4,3}([1 4 6 3 4 2 4 6 6 5 5 5])

    # new com and area
    com_nA, area_nA = get_com_svec_wrapper(1 ,coords_n, nodes_n)
    com_nB, area_nB = get_com_svec_wrapper(2 ,coords_n, nodes_n)
    com_nC, area_nC = get_com_svec_wrapper(3 ,coords_n, nodes_n)
    com_nD, area_nD = get_com_svec_wrapper(4 ,coords_n, nodes_n)

    com_n = SMatrix{4,3}(com_nA[1], com_nB[1], com_nC[1], com_nD[1],
                         com_nA[2], com_nB[2], com_nC[2], com_nD[2],
                         com_nA[3], com_nB[3], com_nC[3], com_nD[3])

    area_n = @SVector [area_nA, area_nB, area_nC, area_nD]

    # copying old nvec to new nvec -> no need to calc it
    nvec_n = SMatrix{4,3}(nvec_o[1], nvec_o[1], nvec_o[1], nvec_o[1],
                          nvec_o[2], nvec_o[2], nvec_o[2], nvec_o[2],
                          nvec_o[3], nvec_o[3], nvec_o[3], nvec_o[3])

    return coords_n, com_n, nvec_n, area_n
end

function get_com_svec_wrapper(j, coords, nodes)
    # wrapper for calculating area and com in Svector form
    x = @SVector [coords[nodes[j,1],1], 
                    coords[nodes[j,2],1], 
                    coords[nodes[j,3],1]]
    y = @SVector [coords[nodes[j,1],2], 
                    coords[nodes[j,2],2], 
                    coords[nodes[j,3],2]]
    z = @SVector [coords[nodes[j,1],3], 
                    coords[nodes[j,2],3], 
                    coords[nodes[j,3],3]]
    dx1 = x[3] - x[1]
    dy1 = y[3] - y[1]
    dz1 = z[3] - z[1]
    dx2 = x[3] - x[2]
    dy2 = y[3] - y[2]
    dz2 = z[3] - z[2]
    cpx = dy1 * dz2 - dz1 * dy2
    cpy = dz1 * dx2 - dx1 * dz2
    cpz = dx1 * dy2 - dy1 * dx2
    area = sqrt(cpx * cpx + cpy * cpy + cpz * cpz) / 2
    cx = (x[1] + x[2] + x[3]) / 3
    cy = (y[1] + y[2] + y[3]) / 3
    cz = (z[1] + z[2] + z[3]) / 3
    com = @SVector [cx, cy, cz]
    return com, area
end

function getVFfromRefinement(n_vf, com_i1, nvec_i1, area_i1, com_i2, nvec_i2, area_i2, vftemp)
    # calculate small vf matrix

	for i1 = 1:1:n_vf

        nvec1 = @SVector [nvec_i1[i1,1], nvec_i1[i1,2], nvec_i1[i1,3]]
        com1 = @SVector [com_i1[i1,1], com_i1[i1,2], com_i1[i1,3]]

		for i2 = 1:1:n_vf

            nvec2 = @SVector [nvec_i2[i2,1], nvec_i2[i2,2], nvec_i2[i2,3]]
            com2 = @SVector [com_i2[i2,1], com_i2[i2,2], com_i2[i2,3]]
            
            phi1, phi2, length_s = calcVFmath(nvec1, nvec2, com1, com2)

            EPSILON = 1E-12
            if phi1 <= EPSILON || phi2 <= EPSILON
                # no viewfactor existing -> vf = 0
                vftemp[Threads.threadid(),1][i1,i2] = 0
                vftemp[Threads.threadid(),2][i2,i1] = 0
			else
                # viewfactor existing -> calc necessary
                vf = phi1 * phi2 / (pi * length_s * length_s)
                vftemp[Threads.threadid(),1][i1,i2] = vf * area_i2[i2]
                vftemp[Threads.threadid(),2][i2,i1] = vf * area_i1[i1]
            end
            
        end

    end

    vf_refined1, vf_refined2 = sumUpVF(n_vf, vftemp, area_i1, area_i2)

    return vf_refined1, vf_refined2
end

function sumUpVF(n_vf, vftemp, area_i1, area_i2)
    # sum up small vf mat
	sum12 = 0.0
	sum21 = 0.0

	tot_area_i1 = 0.0
	tot_area_i2 = 0.0

	for i = 1:1:n_vf
		for j = 1:1:n_vf
			sum12 = sum12 + vftemp[Threads.threadid(),1][i, j] * area_i1[i]
			sum21 = sum21 + vftemp[Threads.threadid(),2][i, j] * area_i2[i]
        end
		tot_area_i1 = tot_area_i1 + area_i1[i]
		tot_area_i2 = tot_area_i2 + area_i2[i]
    end

	vfsum1 = sum12 / tot_area_i1
    vfsum2 = sum21 / tot_area_i2

    return vfsum1, vfsum2
end

function getRefinedTriangles2(coords_o, nvec_o, n_vf)

    com_nn = zeros(n_vf,3)
    nvec_nn = zeros(n_vf,3)
    area_nn = zeros(n_vf)

    coords_n, com_n, nvec_n, area_n = getRefinedTriangles(coords_o, nvec_o)

    counter = 0

    for i = 1:4
        coords_ni, nvec_ni = getSpecificTriangleRefined(i, coords_n, nvec_n)
        coords_nnt, com_nnt, nvec_nnt, area_nnt = getRefinedTriangles(coords_ni, nvec_ni)
        com_nn[(1+counter*4):((counter+1)*4),:] = com_nnt
        nvec_nn[(1+counter*4):((counter+1)*4),:] = nvec_nnt
        area_nn[(1+counter*4):((counter+1)*4)] = area_nnt
        counter = counter + 1
    end

    return com_nn, nvec_nn, area_nn
end

function getSpecificTriangleRefined(number, coords_o, nvec_o)

	# node array same as above -> therefore copied
    nodes_o = SMatrix{4,3}([1 4 6 3 4 2 4 6 6 5 5 5])

	# split up refined triangels into single triangles
    
    nodes_n = @SVector [nodes_o[number,1], nodes_o[number,2], nodes_o[number,3]]

    coords_n = SMatrix{3,3}(coords_o[nodes_n[1], 1], coords_o[nodes_n[2], 1],
                            coords_o[nodes_n[3], 1],
                            coords_o[nodes_n[1], 2], coords_o[nodes_n[2], 2], coords_o[nodes_n[3], 2],
                            coords_o[nodes_n[1], 3], coords_o[nodes_n[2], 3], coords_o[nodes_n[3], 3])

    nvec_n = @SVector [nvec_o[number,1], nvec_o[number,2], nvec_o[number,3]]
    
    return coords_n, nvec_n
end

function getRefinedTriangles3(coords_o, nvec_o, n_vf)

    com_f = zeros(n_vf,3)
    nvec_f = zeros(n_vf,3)
    area_f = zeros(n_vf)

    coords_n1, com_n1, nvec_n1, area_n1 = getRefinedTriangles(coords_o, nvec_o)

    counter = 0

    for i = 1:4
        coords_n1i, nvec_n1i = getSpecificTriangleRefined(i, coords_n1, nvec_n1)
        coords_n2t, com_n2t, nvec_n2t, area_n2t = getRefinedTriangles(coords_n1i, nvec_n1i)

        for j = 1:4
            coords_n2tj, nvec_n2tj = getSpecificTriangleRefined(j, coords_n2t, nvec_n2t)
            coords_n3t, com_n3t, nvec_n3t, area_n3t = getRefinedTriangles(coords_n2tj, nvec_n2tj)
            com_f[(1+counter*4):((counter+1)*4),:] = com_n3t
            nvec_f[(1+counter*4):((counter+1)*4),:] = nvec_n3t
            area_f[(1+counter*4):((counter+1)*4)] = area_n3t
            counter = counter + 1
        end
        
    end

    return com_f, nvec_f, area_f
end

function getRefinedTriangles4(coords_o, nvec_o, n_vf)

    com_f = zeros(n_vf,3)
    nvec_f = zeros(n_vf,3)
    area_f = zeros(n_vf)

    coords_n1, com_n1, nvec_n1, area_n1 = getRefinedTriangles(coords_o, nvec_o)

    counter = 0

    for i = 1:4
        coords_n1i, nvec_n1i = getSpecificTriangleRefined(i, coords_n1, nvec_n1)
        coords_n2t, com_n2t, nvec_n2t, area_n2t = getRefinedTriangles(coords_n1i, nvec_n1i)

        for j = 1:4
            coords_n2tj, nvec_n2tj = getSpecificTriangleRefined(j, coords_n2t, nvec_n2t)
            coords_n3t, com_n3t, nvec_n3t, area_n3t = getRefinedTriangles(coords_n2tj, nvec_n2tj)

            for k = 1:4
                coords_n3tk, nvec_n3tk = getSpecificTriangleRefined(k, coords_n3t, nvec_n3t)
                coords_n4t, com_n4t, nvec_n4t, area_n4t = getRefinedTriangles(coords_n3tk, nvec_n3tk)
                com_f[(1+counter*4):((counter+1)*4),:] = com_n4t
                nvec_f[(1+counter*4):((counter+1)*4),:] = nvec_n4t
                area_f[(1+counter*4):((counter+1)*4)] = area_n4t
                counter = counter + 1
            end

        end
        
    end

    return com_f, nvec_f, area_f
end

function getRefinedTriangles5(coords_o, nvec_o, n_vf)

    com_f = zeros(n_vf,3)
    nvec_f = zeros(n_vf,3)
    area_f = zeros(n_vf)

    coords_n1, com_n1, nvec_n1, area_n1 = getRefinedTriangles(coords_o, nvec_o)

    counter = 0

    for i = 1:4
        coords_n1i, nvec_n1i = getSpecificTriangleRefined(i, coords_n1, nvec_n1)
        coords_n2t, com_n2t, nvec_n2t, area_n2t = getRefinedTriangles(coords_n1i, nvec_n1i)

        for j = 1:4
            coords_n2tj, nvec_n2tj = getSpecificTriangleRefined(j, coords_n2t, nvec_n2t)
            coords_n3t, com_n3t, nvec_n3t, area_n3t = getRefinedTriangles(coords_n2tj, nvec_n2tj)

            for k = 1:4
                coords_n3tk, nvec_n3tk = getSpecificTriangleRefined(k, coords_n3t, nvec_n3t)
                coords_n4t, com_n4t, nvec_n4t, area_n4t = getRefinedTriangles(coords_n3tk, nvec_n3tk)

                for l = 1:4
                    coords_n4tl, nvec_n4tl = getSpecificTriangleRefined(l, coords_n4t, nvec_n4t)
                    coords_n5t, com_n5t, nvec_n5t, area_n5t = getRefinedTriangles(coords_n4tl, nvec_n4tl)
                    com_f[(1+counter*4):((counter+1)*4),:] = com_n5t
                    nvec_f[(1+counter*4):((counter+1)*4),:] = nvec_n5t
                    area_f[(1+counter*4):((counter+1)*4)] = area_n5t
                    counter = counter + 1
                end
            end

        end
        
    end

    return com_f, nvec_f, area_f
end