facts("Periodic images") do
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]

  for i in 1:10
    position = cell * rand(3)
    image = position + cell * rand(Int8, (3, ))
    @fact is_periodic(position, image, cell) --> true
    @fact is_periodic(position, image + rand(3) * 1e-4, cell) --> false

    positions = cell * rand((3, 10))
    images = positions + cell * rand(Int8, (3, 10))
    positions[:, [1, 3, 4]] += rand((3, 3))
    @fact is_periodic(positions, images, cell)[2] --> true
    @fact is_periodic(positions, images, cell)[5:end] --> all
    @fact is_periodic(positions, images, cell)[[1, 3, 4]] --> not(any)
  end
end

facts("Fold back into cell") do
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]

  for i in 1:10
    position = cell * (rand(3) + rand(Int8, (3, )))
    folded = into_cell(position, cell)
    @fact is_periodic(position, folded, cell) --> true
    @fact 0 .≤ inv(cell) * folded .< 1 --> all

    position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
    folded = into_cell(position, cell)
    @fact is_periodic(position, folded, cell) --> all
    @fact 0 .≤ inv(cell) * folded .< 1 --> all
  end
end

facts("Fold back around origin") do
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]

  for i in 1:10
    position = cell * (rand(3) + rand(Int8, (3, )))
    folded = origin_centered(position, cell)
    @fact is_periodic(position, folded, cell) --> true
    @fact -0.5 .< inv(cell) * folded .≤ 0.5 --> all

    position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
    folded = origin_centered(position, cell)
    @fact is_periodic(position, folded, cell) --> all
    @fact -0.5 .< inv(cell) * folded .≤ 0.5 --> all
  end
end
