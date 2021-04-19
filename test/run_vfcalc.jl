
using FARADS_GEOM
using FARADS_MESHING
using FARADS_PLOT
using FARADS_VFCALC

c1 = Cylinder(1.0, 2.0, [0.0,0.0,0.0])
c2 = Cylinder(0.5, 2.0, [0.0,0.0,0.0])

p1 = discretisation(c1, [5,20,30])
reverse_nvec_of_faces!(p1)
delete_face_of_part!(p1, 3)
delete_face_of_part!(p1, 1)

p2 = discretisation(c2, [5,15,30])
delete_face_of_part!(p2, 3)
delete_face_of_part!(p2, 1)

m = compose_mesh([p1, p2])
information_mesh(m)

# plot_mesh_parts(m, shownvec = true)
# plot_mesh_faces(m, shownvec = false)

mvec = make_mesh_vector(m)
vfmat = Array{Float64,2}(undef,mvec.nelements,mvec.nelements)

existing_vf!(vfmat, mvec)
# plot_existing_vf(m, vfmat, 1)

b = buckets_main(m, [0.1,0.1,0.1])
info_bucket_occupation(b)
bvec = make_buckets_vector(b)

blocking_vf!(vfmat, bvec, mvec)
# plot_existing_vf(m, vfmat, 1)

calculate_vf!(vfmat, 0, mvec)

vffaces = vf_elements_to_faces(m, vfmat)
