""" Creates a BCC lattice with a single site """
bcc(T::Type=Float64; unita=u"nm") = Crystal(
    T[-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5] * unit,
    position=T[0, 0, 0] * unit
)

""" Creates an FCC lattice with a single site """
fcc(T::Type=Float64; unit=u"nm") = Crystal(
    T[0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0] * unit,
    position=T[0, 0, 0] * unit
)
