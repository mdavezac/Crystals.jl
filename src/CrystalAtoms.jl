module CrystalAtoms
using Crystals.Structures: Crystal, position_for_crystal
export eachatom
import DataFrames

""" Wrapper around a row/atom in a Crystal for iteration """
immutable CrystalAtom{PARENT <: Crystal}
    """ Parent crystal """
    parent::PARENT
    """ Index in parent """
    index::Int64
end

Base.getindex(a::CrystalAtom, c::Any) = getindex(a.parent, a.index, c)
Base.getindex(a::CrystalAtom, c::Symbol, d::Integer) = getindex(a.parent, a.index, c, d)
Base.setindex!(a::CrystalAtom, v::Any, c::Any) = setindex!(a.parent, v, a.index, c)
function Base.setindex!(atom::CrystalAtom, v::Any, c::Symbol, d::Integer)
    setindex!(atom.parent, v, atom.index, c, d)
end
Base.names(a::CrystalAtom) = names(a.parent)

""" Wraps crystal for iteration purposes """
immutable AtomicIterator{PARENT <: Crystal}
    """ Parent crystal structure """
    parent::PARENT
end

""" Iterator over each atom in the crystal """
eachatom(crystal::Crystal) = AtomicIterator(crystal)
""" Iterator over each atom in the crystal """
DataFrames.eachrow(crystal::Crystal) = eachatom(crystal)

Base.start(itr::AtomicIterator) = 1
Base.done(itr::AtomicIterator, i::Integer) = i > size(itr.parent, 1)
Base.next(itr::AtomicIterator, i::Integer) = (CrystalAtom(itr.parent, i), i + 1)
Base.size(itr::AtomicIterator) = (length(itr.parent), )
Base.length(itr::AtomicIterator) = size(itr.parent, 1)
Base.getindex(itr::AtomicIterator, i::Any) = CrystalAtom(itr.parent, i)
Base.eltype{T <: Crystal}(::Type{AtomicIterator{T}}) = CrystalAtom{T}
end
