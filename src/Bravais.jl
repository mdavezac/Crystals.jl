""" Creates a body-centered lattice with a single site """
bcc(T::Type=Float64; unit=u"nm") = Crystal(
    T[-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5] * unit,
    position=T[0, 0, 0] * unit
)

""" Creates an face centered lattice with a single site """
fcc(T::Type=Float64; unit=u"nm") = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0] * unit,
    position=T[0, 0, 0] * unit
)
