module FARADS_VFCALC

    using StaticArrays
    using ProgressMeter
    using Random
    using LinearAlgebra

    using FARADS_GEOM
    using FARADS_MESHING
    
    include("./existing_fct.jl")
    include("./bucket_types.jl")
    include("./circle_fct.jl")
    include("./bucket_fct.jl")
    include("./blocking_fct.jl")
    include("./viewfactor_fct.jl")
    include("./vfpost_fct.jl")

    include("./blocking_fct_new.jl")

    include("./ray_tri_intersect_fct.jl")

    export EmptyBuckets, Circles, OccBuckets

    export existing_vf!

    export buckets_main, info_bucket_occupation, make_buckets_vector
    # export create_buckets, get_circlepoints
    # export check_bucket_and_parts, check_bucket_and_elements
    export get_elements_inside_bucket

    export blocking_vf!

    export calculate_vf!

    export vf_elements_to_faces
    export quality_check_vf
    export print_array

    # for testing
    export existing_vf_2elem, blocking_vf_2elem

    # new blocking
    export blocking_vf_2elem_new
    export blocking_vf_2elem_bucket
    export blocking_vf_2elem_elem

    export MoellerTrumbore


end