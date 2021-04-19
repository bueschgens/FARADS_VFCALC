function existing_vf!(vfmat, mym::VecMesh3D)
	# check for existing viewfactors
	n_elements = mym.nelements
    println("starting check for existing vf with ",Threads.nthreads()," threads")
    progress = Progress(n_elements, dt=1, barglyphs=BarGlyphs("[|| ]"), barlen=50)
    for i1 = 1:1:n_elements
        next!(progress)
        #for i2 = (i1+1):1:n_elements
        Threads.@threads for i2 = (i1+1):1:n_elements
        #Threads.@threads for i2 = shuffle!(collect((i1+1):n_elements))
            isexisting = checkElementsExistingVF(i1, i2, mym)
            #println(i1, " --> ", i2, " existing vf --> ", isexisting)
            if isexisting
                # vf existing
                vfmat[i1, i2] = 1
                vfmat[i2, i1] = 1
            else
                # no existing vf
                vfmat[i1, i2] = 0
                vfmat[i2, i1] = 0
            end
        end
    end 
    println("check for existing vf done")
end

function existing_vf_2elem(mym::VecMesh3D, i1, i2)
	# run the viewfactor calculation
    isexisting = checkElementsExistingVF(i1, i2, mym)
    println(i1, " --> ", i2, ": existing vf --> ", isexisting)
    return isexisting
end

@inline function checkElementsExistingVF(i1::T2, i2::T2, mym::VecMesh3D)::Bool where T2<:Integer
    # check if viewfactor is existing between two elements
    nvec1 = @SVector [mym.nvec[i1], mym.nvec[i1+1*mym.nelements], mym.nvec[i1+2*mym.nelements]]
    nvec2 = @SVector [mym.nvec[i2], mym.nvec[i2+1*mym.nelements], mym.nvec[i2+2*mym.nelements]]
    com1 = @SVector [mym.com[i1], mym.com[i1+1*mym.nelements], mym.com[i1+2*mym.nelements]]
    com2 = @SVector [mym.com[i2], mym.com[i2+1*mym.nelements], mym.com[i2+2*mym.nelements]]
    phi1, phi2, length_s = calcVFmath(nvec1, nvec2, com1, com2)
    EPSILON = 1E-12
    exist = true
    if phi1 <= EPSILON || phi2 <= EPSILON
        # no viewfactor existing
        exist = false
    end
    return exist
end

@inline function calcVFmath(nvec1, nvec2, com1, com2)
    # do math for viewfactor calculation
    comvec = @SVector [com2[1] - com1[1], com2[2] - com1[2], com2[3] - com1[3]]
    s = sqrt(comvec[1] * comvec[1] + comvec[2] * comvec[2] + comvec[3] * comvec[3])
    skalarProd1 = comvec[1] * nvec1[1] + comvec[2] * nvec1[2] + comvec[3] * nvec1[3]
    skalarProd2 = (-1)*comvec[1] * nvec2[1] + (-1)*comvec[2] * nvec2[2] + (-1)*comvec[3] * nvec2[3]
    n1 = sqrt(nvec1[1] * nvec1[1] + nvec1[2] * nvec1[2] + nvec1[3] * nvec1[3])
    n2 = sqrt(nvec2[1] * nvec2[1] + nvec2[2] * nvec2[2] + nvec2[3] * nvec2[3])
    cos1 = skalarProd1 / (s * n1)
    cos2 = skalarProd2 / (s * n2)
    # catch error: cos only between -1 and 1 possible
    EPSILON2 = 1E-8
    if cos1 < -1 || cos1 > 1
        if cos1 > (-1 - EPSILON2)
            cos1 = -1.0
        elseif cos1 < (1 + EPSILON2)
            cos1 = 1.0
        else 
            error("vf calc cos1")
        end
    end
    if cos2 < -1 || cos2 > 1
        if cos2 > (-1 - EPSILON2)
            cos2 = -1.0
        elseif cos2 < (1 + EPSILON2)
            cos2 = 1.0
        else 
            error("vf calc cos2")
        end
    end
    return cos1, cos2, s
end


