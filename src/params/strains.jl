struct StrainData
  constant_kernel_param::Float32
  household_kernel_param::Float32
end

const StrainTable = SVector{3, StrainData}

function make_strains(constant_kernel_param::Real, household_kernel_param; british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5)::StrainTable
  SA[
    StrainData(constant_kernel_param, household_kernel_param),
    StrainData(constant_kernel_param * british_multiplier, household_kernel_param * british_multiplier),
    StrainData(constant_kernel_param * delta_multiplier, household_kernel_param * delta_multiplier)
  ]
end

getdata(table::StrainTable, strain::StrainKind) = table[UInt8(strain)]