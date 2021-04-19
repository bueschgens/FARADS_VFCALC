abstract type AbstractBuckets end

struct EmptyBuckets{T1<:AbstractFloat, T2<:Integer} <: AbstractBuckets
    nodes::Array{T1,2}
    volumes::Array{T2,2}
    delta::Vector{T1}
    n_buckets_dir::Vector{T2}
end

struct Circles{T1<:AbstractFloat, T2<:Integer} <: AbstractBuckets
    nodes::Array{T1,2}
    elements::Array{T2,2}
end

struct OccBuckets{T1<:AbstractFloat, T2<:Integer} <: AbstractBuckets
    nodes::Array{T1,2}
    volumes::Array{T2,2}
    occupied::Vector{T2}
    buckets2parts::Vector{Vector{T2}}
    buckets2elements::Vector{Vector{T2}}
    delta::Vector{T1}
    n_XYZ::Vector{T2}
end

struct VecOccBuckets{T1<:AbstractFloat, T2<:Integer} <: AbstractBuckets
    nodes::Vector{T1}
    volumes::Vector{T2}
    occupied::Vector{T2}
    buckets2elements_cellvec::Vector{T2}
    buckets2elements_assign::Vector{T2}
    delta::Vector{T1}
    n_XYZ::Vector{T2}
    nbuckets::T2
end