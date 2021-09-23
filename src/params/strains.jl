struct StrainData
  outside_multiplier::Float32
  household_multiplier::Float32
end

const StrainTable = SVector{3, StrainData}

function make_strains(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5)::StrainTable
  SA[
    StrainData(base_multiplier, base_multiplier),
    StrainData(base_multiplier * british_multiplier, base_multiplier * british_multiplier),
    StrainData(base_multiplier* delta_multiplier, base_multiplier * delta_multiplier)
  ]
end

getdata(table::StrainTable, strain::StrainKind) = table[UInt8(strain)]
