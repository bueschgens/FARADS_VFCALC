function vf_elements_to_faces(mym::Mesh3D, vfmat)
	# Elements to Face viewfactors
	n_faces = size(mym.elements2faces,1)

	facearea = zeros(n_faces)
	facevf = zeros(n_faces,n_faces)

	for i = 1:n_faces
		e11 = mym.elements2faces[i,3]
		e12 = mym.elements2faces[i,4]
		
		facearea[i] = sum(mym.area[e11:e12,1])
		
		for j = 1:n_faces
			e21 = mym.elements2faces[j,3]
			e22 = mym.elements2faces[j,4]
			
			vfsplit = vfmat[e11:e12,e21:e22]
			areasplit = zeros(mym.elements2faces[i,2],mym.elements2faces[j,2])
			areasplit[:,1] = mym.area[e11:e12,1]
			for k = 2:mym.elements2faces[j,2]
				areasplit[:,k] = areasplit[:,1] # repeat probieren
			end
			vfsplit .*= areasplit
			facevf[i,j] = sum(vfsplit) ./ facearea[i]
			
		end
	end

	controlvf = sum(facevf,dims = 2)
	print_array(controlvf)

	return facevf
end

function print_array(a::Vector{T}) where T<:Real
    # print 1d array
    println("Printing array with 1 dimensions:")
    for j = 1:size(a,1)
        println(a[j])
    end
    println()
end

function print_array(a::Matrix{T}) where T<:Real
    # print 2d array
    println("Printing array with 2 dimensions:")
    for j = 1:size(a,1)
        println(a[j,:])
    end
    println()
end

function quality_check_vf(mat)
	# check and print quality of given vf matrix
	println("Quality check of viewfactors:")
	controlvf = sum(mat,dims = 2)
	vmin = minimum(controlvf)
	vmax = maximum(controlvf)
	vsum = sum(controlvf)
	n = size(mat,1)
	vmean = vsum / n
	println("    Maximum: ", vmax, " / 1.00")
	println("    Minimum: ", vmin, " / 1.00")
	println("    Mean: ", vmean, " / 1.00")
end

function vf_elements_to_parts(mym::Mesh3D, vfmat)
	# Elements to Part viewfactors
	n_parts = size(mym.elements2parts,1)

	partarea = zeros(n_parts)
	partvf = zeros(n_parts,n_parts)

	for i = 1:n_parts
		e11 = mym.elements2parts[i,3]
		e12 = mym.elements2parts[i,4]
		
		partarea[i] = sum(mym.area[e11:e12,1])
		
		for j = 1:n_parts
			e21 = mym.elements2parts[j,3]
			e22 = mym.elements2parts[j,4]
			
			vfsplit = vfmat[e11:e12,e21:e22]
			areasplit = zeros(mym.elements2parts[i,2],mym.elements2parts[j,2])
			areasplit[:,1] = mym.area[e11:e12,1]
			for k = 2:mym.elements2parts[j,2]
				areasplit[:,k] = areasplit[:,1] # repeat probieren
			end
			vfsplit .*= areasplit
			partvf[i,j] = sum(vfsplit) ./ partarea[i]

		end
	end

	controlvf = sum(partvf,dims = 2)
	print_array(controlvf)

	return partvf
end