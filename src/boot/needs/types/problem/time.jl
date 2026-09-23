abstract type AbstractTime end

export Time

struct Time{T1<:Integer,T2<:Real}
    t  ::Vector{T2}
    te ::T2
    tg ::T2
    tep::T2
end
@adapt_struct Time
