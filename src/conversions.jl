module Conversions

export bar_to_Pa, GPa_to_Pa, kbar_to_GPa, kbar_to_MPa, Bohr_to_A
export eV_to_J, Ry_to_eV, A_to_m, Gyr_to_sec

const bar_to_Pa = 1e5
const GPa_to_Pa = 1e9
const kbar_to_GPa = 0.1
const kbar_to_MPa = 100
const eV_to_J = 1.602176565e-19
const Ry_to_eV = 13.605662285137
const A_to_m = 1e-10
const Bohr_to_A = 0.52918
const Gyr_to_sec = 1e9 * 365 * 24 * 60 * 60

end # module
