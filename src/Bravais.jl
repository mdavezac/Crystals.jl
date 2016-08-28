" Creates a BCC lattice with a single site "
bcc(T::Type=Float64; scale=1) = Crystal(
    T[-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5],
    scale,
    tposition=T[0 0 0;],
    specie=['A']
)

" Creates an FCC lattice with a single site "
fcc(T::Type=Float64; scale=1) = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0],
    scale,
    tposition=T[0 0 0],
    specie=['A']
)
