" Rock-salt lattice "
rock_salt(T::Type=Float64; unit=u"nm") = Crystal(
    eye(T, 3) * unit,
    position=T[0.5, 0.5, 0.5] * unit,
    position=T[0, 0, 0] * unit,
    species=['A', 'B']
)

" Zinc-blende lattice "
zinc_blende(T::Type=Float64; unit=u"nm") = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0] * unit,
    tposition=T[0 0 0; 0.25 0.25 0.25] * unit,
    species=['A', 'B']
)

" Diamond lattice "
diamond(T::Type=Float64; unit=u"nm") = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0] * unit,
    tposition=T[0 0 0; 0.25 0.25 0.25] * unit
)

" Wurtzite lattice "
wurtzite(T::Type=Float64; unit=u"nm") = Crystal(
    T[0.5 0.5 0; -sqrt(3)/2 sqrt(3)/2 0; 0 0 1] * unit,
    tposition=T[0.5  sqrt(3)/6 0;
                0.5 -sqrt(3)/6 0.5;
                0.5  sqrt(3)/6 0.25;
                0.5 -sqrt(3)/6 0.75] * unit,
    species=['A', 'A', 'B', 'B']
)
