" Rock-salt lattice "
rock_salt(T::Type=Float64; scale=1) = Crystal(
    eye(T, 3),
    scale,
    tposition=T[0 0 0; 0.5 0.5 0.5],
    specie=['A', 'B']
)

" Zinc-blende lattice "
zinc_blende(T::Type=Float64; scale=1) = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0],
    scale,
    tposition=T[0 0 0; 0.25 0.25 0.25],
    specie=['A', 'B']
)

" Diamond lattice "
diamond(T::Type=Float64; scale=1) = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0],
    scale,
    tposition=T[0 0 0; 0.25 0.25 0.25],
    specie=['A', 'A']
)

" Wurtzite lattice "
wurtzite(T::Type=Float64; scale=1) = Crystal(
    T[0.5 0.5 0; -sqrt(3)/2 sqrt(3)/2 0; 0 0 1],
    scale,
    tposition=T[0.5  sqrt(3)/6 0;
                0.5 -sqrt(3)/6 0.5;
                0.5  sqrt(3)/6 0.25;
                0.5 -sqrt(3)/6 0.75],
    specie=['A', 'A', 'B', 'B']
)
