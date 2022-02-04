"""
    Util

Author: David James, davidabraham@ucla.edu
Date: 20211225
Notes: extra functions that don't fit into other libraries

Contains:
- get_date
- logspace
"""
module Util

using Dates

export get_date, logspace

"""
    get_date()

Generates a string of the current date and timestime in the following
format: YYYYMMDD-HHMMSS.sss
"""
function get_date()

    d = Dates.format(now(), dateformat"Ymmdd-HHMMSS.sss")

    return d

end

"""
    logspace(a::Real, b::Real, n::Integer)

Generates a range of length `n` in logspace going from start `a` to
end `b`.
"""
function logspace(a::Real, b::Real, n::Integer)

    return 10 .^ range(log10(a), log10(b), length=n)

end

end # module
