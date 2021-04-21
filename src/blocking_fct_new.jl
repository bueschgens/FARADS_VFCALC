function blocking_vf_new!(vfmat, myb::VecOccBuckets, mym::VecMesh3D)
	# run the blocking calculation
	n_elements = mym.nelements
	println("starting blocking check with ",Threads.nthreads()," threads")
	offset = zeros(3) # get node offset -> for all coords positiv
	offset[1] = minimum(mym.nodes[1:mym.nnodes])
	offset[2] = minimum(mym.nodes[mym.nnodes+1:2*mym.nnodes])
	offset[3] = minimum(mym.nodes[2*mym.nnodes+1:3*mym.nnodes])
	pass = zeros(1) # not usable in parallel
	progress = Progress(n_elements, dt=1, barglyphs=BarGlyphs("[|| ]"), barlen=50)
	for i1 = 1:1:n_elements
		next!(progress)
		#for i2 = (i1+1):1:n_elements
		#Threads.@threads for i2 = (i1+1):1:n_elements
		Threads.@threads for i2 = shuffle!(collect((i1+1):n_elements))
            if vfmat[i1,i2] > 0 # --> vf existing
                hasShadowing = bucketwalk_with_return!(i1, i2, mym, myb, offset, pass)
                if hasShadowing == 1
                    # default 0, optional pass[1] -> hitten element
                    vfmat[i1, i2] = 0
                    vfmat[i2, i1] = 0
                end
            end
		end
	end
    println("blocking check done")
end

function blocking_vf_2elem_new(myb::VecOccBuckets, mym::VecMesh3D, i1, i2)
	# run the blocking calculation
    # necessary here: offset on p1 and p2
    # --> translate all nodes into positive coords
	# --> for bucket numbering
	# nodes @Ray triangle also with offset (pointA, pointB, pointC)
	offset = zeros(3) # get node offset -> for all coords positiv
	offset[1] = minimum(mym.nodes[1:mym.nnodes])
	offset[2] = minimum(mym.nodes[mym.nnodes+1:2*mym.nnodes])
	offset[3] = minimum(mym.nodes[2*mym.nnodes+1:3*mym.nnodes])
    max_steps = get_max_steps(myb)
    buckets = Vector{Int64}(undef,max_steps)
    nbuckets = 0
	hitelem = 0
	if existing_vf_2elem(mym, i1, i2)
        p1 = @SVector [mym.com[i1]-offset[1], 
                        mym.com[i1+1*mym.nelements]-offset[2], 
                        mym.com[i1+2*mym.nelements]-offset[3]]
        p2 = @SVector [mym.com[i2]-offset[1], 
                        mym.com[i2+1*mym.nelements]-offset[2], 
                        mym.com[i2+2*mym.nelements]-offset[3]]
        p1p2 = @SVector[p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3]]
        dir = SetDir(p1p2)
		nbuckets = bucketwalk_with_return!(buckets, i1, i2, myb, p1, p2, dir)
        for i = 1:nbuckets
            bNr = buckets[i]
            hit, hitelem = RayTriangleMain(myb, bNr, mym, p1, p1p2, i1, i2, offset)
			if hit == 1
				println("--> Element hitted: ", hitelem, " in bucket: ", bNr)
				break
			end
        end
	end
    println(buckets[1:nbuckets])
    println("blocking check done for ", i1, " and ", i2)
    return buckets[1:nbuckets], hitelem
end

function blocking_vf_2elem_bucket(myb::VecOccBuckets, mym::VecMesh3D, i1, i2, bNr)
	# run the blocking calculation
    # necessary here: offset on p1 and p2
    # --> translate all nodes into positive coords
	# --> for bucket numbering
	# nodes @Ray triangle also with offset (pointA, pointB, pointC)
	offset = zeros(3) # get node offset -> for all coords positiv
	offset[1] = minimum(mym.nodes[1:mym.nnodes])
	offset[2] = minimum(mym.nodes[mym.nnodes+1:2*mym.nnodes])
	offset[3] = minimum(mym.nodes[2*mym.nnodes+1:3*mym.nnodes])
	hitelem = 0
	if existing_vf_2elem(mym, i1, i2)
        p1 = @SVector [mym.com[i1]-offset[1], 
                        mym.com[i1+1*mym.nelements]-offset[2], 
                        mym.com[i1+2*mym.nelements]-offset[3]]
        p2 = @SVector [mym.com[i2]-offset[1], 
                        mym.com[i2+1*mym.nelements]-offset[2], 
                        mym.com[i2+2*mym.nelements]-offset[3]]
        p1p2 = @SVector[p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3]]
        dir = SetDir(p1p2)
		hit, hitelem = RayTriangleMain(myb, bNr, mym, p1, p1p2, i1, i2, offset)
		if hit == 1
			println("--> Element hitted: ", hitelem, " in bucket: ", bNr)
		end
	end
    println("blocking check done for ", i1, " and ", i2)
    return hitelem
end

function blocking_vf_2elem_elem(mym::VecMesh3D, i1, i2, ib)
	# run the blocking calculation
    # necessary here: offset on p1 and p2
    # --> translate all nodes into positive coords
	# --> for bucket numbering
	# nodes @Ray triangle also with offset (pointA, pointB, pointC)
	offset = zeros(3) # get node offset -> for all coords positiv
	offset[1] = minimum(mym.nodes[1:mym.nnodes])
	offset[2] = minimum(mym.nodes[mym.nnodes+1:2*mym.nnodes])
	offset[3] = minimum(mym.nodes[2*mym.nnodes+1:3*mym.nnodes])
	if existing_vf_2elem(mym, i1, i2)
        p1 = @SVector [mym.com[i1]-offset[1], 
                        mym.com[i1+1*mym.nelements]-offset[2], 
                        mym.com[i1+2*mym.nelements]-offset[3]]
        p2 = @SVector [mym.com[i2]-offset[1], 
                        mym.com[i2+1*mym.nelements]-offset[2], 
                        mym.com[i2+2*mym.nelements]-offset[3]]
        p1p2 = @SVector[p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3]]
        dir = SetDir(p1p2)
		# get information on potential blocking element
		part = mym.elementstatus[ib+2*mym.nelements] # returns part number
		nodes1 = mym.elements[ib]
		nodes2 = mym.elements[ib+1*mym.nelements]
		nodes3 = mym.elements[ib+2*mym.nelements]
		nodeoffset = mym.nodes2parts[part+2*mym.nparts]
		# offset muss hier auf die points drauf
		pointA = @SVector [mym.nodes[nodeoffset - 1 + nodes1]-offset[1], 
							mym.nodes[nodeoffset - 1 + nodes1+1*mym.nnodes]-offset[2],
							mym.nodes[nodeoffset - 1 + nodes1+2*mym.nnodes]-offset[3]]
		pointB = @SVector [mym.nodes[nodeoffset - 1 + nodes2]-offset[1],
							mym.nodes[nodeoffset - 1 + nodes2+1*mym.nnodes]-offset[2],
							mym.nodes[nodeoffset - 1 + nodes2+2*mym.nnodes]-offset[3]]
		pointC = @SVector [mym.nodes[nodeoffset - 1 + nodes3]-offset[1],
							mym.nodes[nodeoffset - 1 + nodes3+1*mym.nnodes]-offset[2],
							mym.nodes[nodeoffset - 1 + nodes3+2*mym.nnodes]-offset[3]]
		println(ib, " ", pointA, " ", pointB, " ", pointC)
		hit = MoellerTrumbore(pointA, pointC, pointB, p1, p1p2)
		if hit == 1
			println("--> Element hitted")
		end
	end
    println("blocking check done for ", i1, " and ", i2)
end

function get_max_steps(myb::VecOccBuckets)
    len = sqrt(myb.n_XYZ[1]^2 + myb.n_XYZ[2]^2 + myb.n_XYZ[3]^2)
    n_max = round(Int64, 2 * len)
    return n_max
end

function bucketwalk_with_return!(buckets, i1, i2, myb, p1, p2, dir)
	
	bNr_p1, bIdx_p1 = getBucketNumberFromPoint(p1, dir, myb, 1)
	bNr_p2, bIdx_p2 = getBucketNumberFromPoint(p2, dir, myb, 2)
	# bucket number here +1 (comparable to c)
	# println("Start bucket Nr. ", bNr_p1, "  / End bucket Nr. ", bNr_p2)

    buckets[1] = bNr_p1
    nbuckets = 1

	if bNr_p1 == bNr_p2 # Start Bucket == End Bucket?
		# println("Start Bucket equals End Bucket")
    else
        # Start Bucket run to End Bucket
        bIdx_delta = @SVector [bIdx_p2[1] - bIdx_p1[1], 
                                bIdx_p2[2] - bIdx_p1[2], 
                                bIdx_p2[3] - bIdx_p1[3]]
        sum_dir = abs(dir[1]) + abs(dir[2]) + abs(dir[3])
        println("bucket walk case: ", sum_dir)
        if sum_dir == 3
            nb = bucketWalk_3D(buckets, p1, i1, p2, i2, dir, bIdx_p1, bNr_p2, myb, bIdx_delta)
        else 
            if sum_dir == 2
                nb = bucketWalk_2D(buckets, p1, i1, p2, i2, dir, bIdx_p1, bNr_p2, myb, bIdx_delta)
            else 
                nb = bucketWalk_1D(buckets, p1, i1, p2, i2, dir, bIdx_p1, bNr_p2, myb, bIdx_delta)
            end
        end
	end
    nbuckets = nbuckets + nb
	return nbuckets
end

function SetDir(p1p2)
    # get 3D direction of bucketwalk
	EPSILON = 1E-8
    # dim 1
	if p1p2[1] > EPSILON	
		dir1 = 1
	elseif p1p2[1] < (-EPSILON)
			dir1 = (-1)
	else
		dir1 = 0
	end
    # dim 2
	if p1p2[2] > EPSILON
		dir2 = 1
	elseif p1p2[2] < (-EPSILON)
		dir2 = (-1)
	else
		dir2 = 0
	end
    # dim 3
	if p1p2[3] > EPSILON
		dir3 = 1
	elseif p1p2[3] < (-EPSILON)
		dir3 = (-1)
	else
		dir3 = 0
	end
	dir = @SVector [dir1, dir2, dir3]
	return dir
end

function getBucketNumberFromPoint(p, dir, myb, point1or2)
	# get bucket number from point
	tolerance = 1E-8 # default point 1
	if point1or2 == 2 
		tolerance = (-tolerance)
	end
	bCoord = @SVector [floor((p[1] + (dir[1] * tolerance)) / 
						myb.delta[1]) * myb.delta[1],
					   floor((p[2] + (dir[2] * tolerance)) / 
					   myb.delta[2]) * myb.delta[2],
					   floor((p[3] + (dir[3] * tolerance)) / 
					   myb.delta[3]) * myb.delta[3] ]
	bIdx = @SVector [round(Int64, bCoord[1] / myb.delta[1]),
						round(Int64, bCoord[2] / myb.delta[2]),
						round(Int64, bCoord[3] / myb.delta[3])]
	bNr = Int64(bIdx[1] + (bIdx[2] * myb.n_XYZ[1]) + (bIdx[3] * myb.n_XYZ[1] * myb.n_XYZ[2]))
	bNr = bNr + 1 # here: bNr + 1
	return bNr, bIdx
end

function bucketWalk_3D(buckets, p1, i1, p2, i2, dir, bIdx_p1, bNr_p2, myb, bIdx_delta)
	# bucket walk case for 3 dimensional connection between start- and endpoint

	# XY-Plane
    line_1_slope = (p2[2] - p1[2]) / (p2[1] - p1[1])
    line_1_axisIc = p1[2] - (line_1_slope * p1[1])
	# YZ-Plane
	line_2_slope = (p2[3] - p1[3]) / (p2[2] - p1[2])
    line_2_axisIc = p1[3] - (line_2_slope * p1[2])
	# ZX-Plane
	line_3_slope = (p2[1] - p1[1]) / (p2[3] - p1[3])
    line_3_axisIc = p1[1] - (line_3_slope * p1[3])

    bEdges = getBucketEdges(p1, dir, myb)
	bIdx_step = @SVector [bIdx_p1[1], bIdx_p1[2], bIdx_p1[3]]

    count = 0 # number of walked buckets
    while count < 10000 # max operations
		count = count + 1

        pointX1 = bEdges[1]
        pointX2 = line_1_slope * bEdges[1] + line_1_axisIc
        pointX3 = (bEdges[1] - line_3_axisIc) / line_3_slope

        pointY1 = (bEdges[2] - line_1_axisIc) / line_1_slope
        pointY2 = bEdges[2]
        pointY3 = line_2_slope * bEdges[2] + line_2_axisIc

        pointZ1 = line_3_slope * bEdges[3] + line_3_axisIc
        pointZ2 = (bEdges[3] - line_2_axisIc) / line_2_slope
        pointZ3 = bEdges[3]

		vecLenX = sqrt((pointX1-p1[1]) * (pointX1-p1[1]) + 
						(pointX2-p1[2]) * (pointX2-p1[2]) + 
						(pointX3-p1[3]) * (pointX3-p1[3]))

		vecLenY = sqrt((pointY1-p1[1]) * (pointY1-p1[1]) + 
						(pointY2-p1[2]) * (pointY2-p1[2]) + 
						(pointY3-p1[3]) * (pointY3-p1[3]))

		vecLenZ = sqrt((pointZ1-p1[1]) * (pointZ1-p1[1]) + 
						(pointZ2-p1[2]) * (pointZ2-p1[2]) + 
						(pointZ3-p1[3]) * (pointZ3-p1[3]))

		vecLen = @SVector [vecLenX, vecLenY, vecLenZ]

		# println("X: ", pointX1, ", ", pointX2, ", ", pointX3, " --> vecLen = ", vecLenX)
		# println("Y: ", pointY1, ", ", pointY2, ", ", pointY3, " --> vecLen = ", vecLenY)
		# println("Z: ", pointZ1, ", ", pointZ2, ", ", pointZ3, " --> vecLen = ", vecLenZ)

        bIdx_step, bEdges = getNextBucket(bEdges, dir, myb, vecLen, bIdx_step)
		bNr_step = bIdx_step[1] + (bIdx_step[2] * myb.n_XYZ[1]) + 
					(bIdx_step[3] * myb.n_XYZ[1] * myb.n_XYZ[2])
		bNr_step = bNr_step + 1 # here: bNr + 1
        # println("Next Bucket: ", bNr_step)
        buckets[count+1] = bNr_step

        if bNr_step == bNr_p2
            # println("End Bucket reached: ", bNr_p2)
            break
		end

        # checkup if more steps in direction as necessary between start- and endbucket
		for i = 1:1:3
			# println("dim ", i, " --> ", abs(bIdx_step[i] - bIdx_p1[i]), " / ", abs(bIdx_delta[i]))
			if abs(bIdx_step[i] - bIdx_p1[i]) > abs(bIdx_delta[i])
				error("Bucketwalk 3 did not reach endbucket: i1 = ", i1, ", i2 = ", i2)
			end
        end
	
	end

    return count
end


function getBucketEdges(p1, dir, myb)
	# done at the beginning of bucketwalk
	# takes startpoint and returns the next bucket planes 
	# which get cut by the connection between start- and endpoint
	# if less than 3 dimensions keeps p1 coordinate in unchanged dimension
    tol = 1E-8 # step tolernace
    # X-Richtung
    if dir[1] > 0
        bcoord1 = ceil((p1[1] + tol) / myb.delta[1]) * myb.delta[1]
    elseif dir[1] < 0
        bcoord1 = floor((p1[1] - tol) / myb.delta[1]) * myb.delta[1]
    else
        bcoord1 = p1[1]
    end
    # Y-Richtung
    if dir[2] > 0
        bcoord2 = ceil((p1[2] + tol) / myb.delta[2]) * myb.delta[2]
    elseif dir[2] < 0
        bcoord2 = floor((p1[2] - tol) / myb.delta[2]) * myb.delta[2]
    else
        bcoord2 = p1[2]
    end
    # Z-Richtung
    if dir[3] > 0
        bcoord3 = ceil((p1[3] + tol) / myb.delta[3]) * myb.delta[3]
    elseif dir[3] < 0
        bcoord3 = floor((p1[3] - tol) / myb.delta[3]) * myb.delta[3]
    else 
        bcoord3 = p1[3]
    end
	bCuttingPlaneCoord = @SVector [bcoord1, bcoord2, bcoord3]
	return bCuttingPlaneCoord
end


function getNextBucket(bEdge_old, dir, myb, vecLen, bIdx_step_old)
	# check which direction to go with shortest vector
	# update bIdx_step with adding/substracting 1 in the direction
	# update bCuttingPlaneCoord with adding/substracting delta

    tol = 1E-8 # length tolerance

	vecLen_min = argmin(vecLen) # returns index of minimum
	
	bIdx_step_new1 = bIdx_step_old[1]
	bIdx_step_new2 = bIdx_step_old[2]
	bIdx_step_new3 = bIdx_step_old[3]

	bEdge_new1 = bEdge_old[1]
	bEdge_new2 = bEdge_old[2]
	bEdge_new3 = bEdge_old[3]

    if vecLen_min == 1
        # println("Next Bucket in X-Direction")
        bIdx_step_new1 = bIdx_step_new1 + dir[1]
        bEdge_new1 = bEdge_new1 + dir[1] * myb.delta[1]
        if abs(vecLen[1] - vecLen[2]) < tol
            # println("Next Bucket in XY-Direction")
            bIdx_step_new2 = bIdx_step_new2 + dir[2]
			bEdge_new2 = bEdge_new2 + dir[2] * myb.delta[2]
		end
        if abs(vecLen[1] - vecLen[3]) < tol
            # println("Next Bucket in XZ-Direction")
            bIdx_step_new3 = bIdx_step_new3 + dir[3]
            bEdge_new3 = bEdge_new3 + dir[3] * myb.delta[3]
		end
    else
        if vecLen_min == 2
            # println("Next Bucket in Y-Direction")
            bIdx_step_new2 = bIdx_step_new2 + dir[2]
            bEdge_new2 = bEdge_new2 + dir[2] * myb.delta[2]
            if abs(vecLen[2] - vecLen[1]) < tol
                # println("Next Bucket in YX-Direction")
                bIdx_step_new1 = bIdx_step_new1 + dir[1]
                bEdge_new1 = bEdge_new1 + dir[1] * myb.delta[1]
			end
            if abs(vecLen[2] - vecLen[3]) < tol
                # println("Next Bucket in YZ-Direction")
                bIdx_step_new3 = bIdx_step_new3 + dir[3]
                bEdge_new3 = bEdge_new3 + dir[3] * myb.delta[3]
			end
        else
            if vecLen_min == 3
                # println("Next Bucket in Z-Direction")
                bIdx_step_new3 = bIdx_step_new3 + dir[3]
                bEdge_new3 = bEdge_new3 + dir[3] * myb.delta[3]
                if abs(vecLen[3] - vecLen[1]) < tol
                    # println("Next Bucket in ZX-Direction")
                    bIdx_step_new1 = bIdx_step_new1 + dir[1]
                    bEdge_new1 = bEdge_new1 + dir[1] * myb.delta[1]
				end
                if abs(vecLen[3] - vecLen[2]) < tol
                    # println("Next Bucket in ZY-Direction")
                    bIdx_step_new2 = bIdx_step_new2 + dir[2]
                    bEdge_new2 = bEdge_new2 + dir[2] * myb.delta[2]
				end
            end
        end
	end

	bIdx_step_new = @SVector [bIdx_step_new1, bIdx_step_new2, bIdx_step_new3]
	bEdge_new = @SVector [bEdge_new1, bEdge_new2, bEdge_new3]

	return bIdx_step_new, bEdge_new
end


function bucketWalk_2D(buckets, p1, i1, p2, i2, dir, bIdx_p1, bNr_p2, myb, bIdx_delta)
	# bucket walk case for 2 dimensional connection between start- and endpoint

	absDir = @SVector [abs(dir[1]), abs(dir[2]), abs(dir[3])]
	zero = argmin(absDir)
	first = findfirst(isequal(1), absDir)
	last = findlast(isequal(1), absDir)
	line_slope = (p2[last] - p1[last]) / (p2[first] - p1[first])
	line_axisIc = p1[last] - (line_slope * p1[first])
	bEdges = getBucketEdges(p1, dir, myb)
	bIdx_step = @SVector [bIdx_p1[1],bIdx_p1[2],bIdx_p1[3]]

    count = 0
    while count < 10000
		count = count + 1

		points_zero1 = bEdges[zero]
		points_zero2 = bEdges[zero]
		points_first1 = bEdges[first]
		points_first2 = (bEdges[last] - line_axisIc) / line_slope
		points_last1 = line_slope * bEdges[first] + line_axisIc
		points_last2 = bEdges[last]
		
		if zero == 1
			point1 = @SVector [points_zero1, points_first1, points_last1] 
			point2 = @SVector [points_zero2, points_first2, points_last2]
		elseif zero == 2
			point1 = @SVector [points_first1, points_zero1, points_last1] 
			point2 = @SVector [points_first2, points_zero2, points_last2]
		else #zero == 3
			point1 = @SVector [points_first1, points_last1, points_zero1] 
			point2 = @SVector [points_first2, points_last2, points_zero2]
		end

		vecLen_first = sqrt((point1[1]-p1[1]) * (point1[1]-p1[1]) + 
							(point1[2]-p1[2]) * (point1[2]-p1[2]) + 
							(point1[3]-p1[3]) * (point1[3]-p1[3]))

		vecLen_last = sqrt((point2[1]-p1[1]) * (point2[1]-p1[1]) + 
							(point2[2]-p1[2]) * (point2[2]-p1[2]) + 
							(point2[3]-p1[3]) * (point2[3]-p1[3]))

		if zero == 1
			vecLen = @SVector [10000, vecLen_first, vecLen_last]
		elseif zero == 2
			vecLen = @SVector [vecLen_first, 10000, vecLen_last]
		else #zero == 3
			vecLen = @SVector [vecLen_first, vecLen_last, 10000]
		end

		bIdx_step, bEdges = getNextBucket(bEdges, dir, myb, vecLen, bIdx_step)
		bNr_step = bIdx_step[1] + (bIdx_step[2] * myb.n_XYZ[1]) + 
					(bIdx_step[3] * myb.n_XYZ[1] * myb.n_XYZ[2])
		bNr_step = bNr_step + 1 # here: bNr + 1
        # println("Next Bucket: ", bNr_step)
        buckets[count+1] = bNr_step

        if bNr_step == bNr_p2
            # println("End Bucket reached: ", bNr_p2)
            break
		end

		# checkup if more steps in direction as necessary between start- and endbucket
		for i = 1:1:3
			# println("dim ", i, " --> ", abs(bIdx_step[i] - bIdx_p1[i]), " / ", abs(bIdx_delta[i]))
			if abs(bIdx_step[i] - bIdx_p1[i]) > abs(bIdx_delta[i])
				error("Bucketwalk 2 did not reach endbucket: i1 = ", i1, ", i2 = ", i2)
			end
        end
		
    end

    return count
end


function bucketWalk_1D(buckets, p1, i1, p2, i2, dir, bIdx_p1, bNr_p2, myb, bIdx_delta)
	# bucket walk case for 1 dimensional connection between start- and endpoint
	# only cares about the dimension which changes -> called one
	
	absDir = @SVector [abs(dir[1]), abs(dir[2]), abs(dir[3])]
	one = findfirst(isequal(1), absDir)
	bEdges = getBucketEdges(p1, dir, myb)
	bIdx_step = @SVector [bIdx_p1[1],bIdx_p1[2],bIdx_p1[3]] # init bIdx with startbucket

	count = 0
    while count < 10000
		count = count + 1

		point = @SVector [bEdges[1],bEdges[2],bEdges[3]]
		vecLen_one = sqrt((point[1]-p1[1]) * (point[1]-p1[1]) + 
						(point[2]-p1[2]) * (point[2]-p1[2]) + 
						(point[3]-p1[3]) * (point[3]-p1[3]))

		if one == 1
			vecLen = @SVector [vecLen_one, 10000, 10000]
		elseif one == 2
			vecLen = @SVector [10000, vecLen_one, 10000]
		else # one == 3
			vecLen = @SVector [10000, 10000, vecLen_one]
		end

        bIdx_step, bEdges = getNextBucket(bEdges, dir, myb, vecLen, bIdx_step)
		bNr_step = bIdx_step[1] + (bIdx_step[2] * myb.n_XYZ[1]) + 
					(bIdx_step[3] * myb.n_XYZ[1] * myb.n_XYZ[2])
		bNr_step = bNr_step + 1 # muss weil aus c so
        # println("Next Bucket: ", bNr_step)
        buckets[count+1] = bNr_step

        if bNr_step == bNr_p2
            # println("End Bucket reached: ", bNr_p2)
            break
		end

		# checkup if more steps in direction as necessary between start- and endbucket
		for i = 1:1:3
			# println("dim ", i, " --> ", abs(bIdx_step[i] - bIdx_p1[i]), " / ", abs(bIdx_delta[i]))
			if abs(bIdx_step[i] - bIdx_p1[i]) > abs(bIdx_delta[i])
				error("Bucketwalk 1 did not reach endbucket: i1 = ", i1, ", i2 = ", i2)
			end
        end

    end

	return count
end

function RayTriangleMain(myb, bNr, mym, orig, dir, i1, i2, offset)
	# run through all elements of occupied bucket

	hit = 0
    hit_element = 0

	#nbuckets = size(myb.occupied,1)
	nbuckets = myb.nbuckets
	e1 = myb.buckets2elements_assign[bNr+1*nbuckets]
	e2 = myb.buckets2elements_assign[bNr+2*nbuckets]
	# println("bucket ", bNr, " from ", e1, " to ", e2)
	
	#for i_c = 1:size(myb.buckets2elements[bNr],1)
	for i_c = e1:e2

		#i_elem = myb.buckets2elements[bNr][i_c]
		i_elem = myb.buckets2elements_cellvec[i_c]
		
		if i_elem != i1 && i_elem != i2

			# get information on triangle
			part = mym.elementstatus[i_elem+2*mym.nelements] # returns part number
			nodes1 = mym.elements[i_elem]
			nodes2 = mym.elements[i_elem+1*mym.nelements]
			nodes3 = mym.elements[i_elem+2*mym.nelements]
			nodeoffset = mym.nodes2parts[part+2*mym.nparts]
			# offset muss hier auf die points drauf
			pointA = @SVector [mym.nodes[nodeoffset - 1 + nodes1]-offset[1], 
								mym.nodes[nodeoffset - 1 + nodes1+1*mym.nnodes]-offset[2],
								mym.nodes[nodeoffset - 1 + nodes1+2*mym.nnodes]-offset[3]]
			pointB = @SVector [mym.nodes[nodeoffset - 1 + nodes2]-offset[1],
								mym.nodes[nodeoffset - 1 + nodes2+1*mym.nnodes]-offset[2],
								mym.nodes[nodeoffset - 1 + nodes2+2*mym.nnodes]-offset[3]]
			pointC = @SVector [mym.nodes[nodeoffset - 1 + nodes3]-offset[1],
								mym.nodes[nodeoffset - 1 + nodes3+1*mym.nnodes]-offset[2],
								mym.nodes[nodeoffset - 1 + nodes3+2*mym.nnodes]-offset[3]]

			# println(i_elem, " ", pointA, " ", pointB, " ", pointC)
			hit = MoellerTrumbore(pointA, pointC, pointB, orig, dir)

			if hit == 1
				hit_element = i_elem
				# println("Hitten element: ", hit_element)
				break
			end

		end

	end

	return hit, hit_element
end


