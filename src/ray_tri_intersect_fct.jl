
function MoellerTrumbore(v0, v1, v2, orig, dir)
	hit = 0
	EPSILON = 1E-10
	edge1 = @SVector [v1[1]-v0[1], v1[2]-v0[2], v1[3]-v0[3]]
	edge2 = @SVector [v2[1]-v0[1], v2[2]-v0[2], v2[3]-v0[3]]
	pvec = @SVector [(dir[2]*edge2[3]) - (dir[3]*edge2[2]),
						(dir[3]*edge2[1]) - (dir[1]*edge2[3]),
						(dir[1]*edge2[2]) - (dir[2]*edge2[1])]

	det = edge1[1]*pvec[1] + edge1[2]*pvec[2] + edge1[3]*pvec[3]
	# println("    starting MoellerTrumbore algorithm")
	if det < EPSILON
		# println("        no intersection: det < EPSILON")
	else 
		tvec = @SVector [orig[1]-v0[1], orig[2]-v0[2], orig[3]-v0[3]]
		u = tvec[1]*pvec[1] + tvec[2]*pvec[2] + tvec[3]*pvec[3]
		if u < 0 || u > det
			# println("        no intersection: u wrong")
		else
			qvec = @SVector [(tvec[2]*edge1[3]) - (tvec[3]*edge1[2]),
						(tvec[3]*edge1[1]) - (tvec[1]*edge1[3]),
						(tvec[1]*edge1[2]) - (tvec[2]*edge1[1])]

			v = dir[1]*qvec[1] + dir[2]*qvec[2] + dir[3]*qvec[3]
			if v < 0 || (u + v) > det
				# println("        no intersection: v wrong")
			else
				t = edge2[1]*qvec[1] + edge2[2]*qvec[2] + edge2[3]*qvec[3]
				inv_det = 1 / det
				t = t * inv_det
				u = u * inv_det
				v = v * inv_det
				# println("        last else check with t = ", t)
				if (t >= 0 && t < 1) 
					is = @SVector [orig[1] + t * dir[1], 
									orig[2] + t * dir[2], 
									orig[3] + t * dir[3]]
					println("        RT-Intersection at: ", is)
					hit = 1
				end
			end
		end
	end
	return hit
end


