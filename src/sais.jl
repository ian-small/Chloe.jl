#=
 * sais
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 =#

struct Alignment <: AbstractVector{UInt8}
    seq1::LongDNA{2} # sequence * sequence[1:end -1]
    lenseq1::Int32 # length of original sequence, not length of seq1 
    seq2::LongDNA{2} # sequence * sequence[1:end -1]
    lenseq2::Int32 # length of original sequence, not length of seq2 
end

function Alignment(seq1::CircularSequence, seq2::CircularSequence)
    return Alignment(seq1.sequence, seq1.length, seq2.sequence, seq2.length)
end

Base.size(a::Alignment) = (length(a.seq1) + 1 + length(a.seq2), )

@inline fulcrum(a::Alignment) = length(a.seq1) + 1

@inline function Base.getindex(a::Alignment, i::Integer)
    f = fulcrum(a)
    if i < f; return reinterpret(UInt8,a.seq1[i]); end
    if i > f; return reinterpret(UInt8,a.seq2[i - f]); end
    return reinterpret(UInt8,DNA_Gap)
end

Base.setindex!(a::Alignment, v, i::Integer) = value #don't do anything


struct IntVector <: AbstractVector{Integer}
    vec::Array{Int,1}
    off::Int
end
Base.size(v::IntVector) = (length(v.vec)-v.off,)
Base.getindex(v::IntVector, key) = v.vec[v.off + key]
function Base.setindex!(v::IntVector, value, key)
    v.vec[v.off+key] = value
    return value
end

function getcounts(T::AbstractVector{<:Integer}, C::IntVector, n::Int, k::Int)
    for i = 1:k
        C[i] = 0
    end
    for i = 1:n
        C[T[i]+1] += 1
    end
end

function getbuckets(C::IntVector, B::IntVector, k::Int, isend::Bool)
    s = 0
    if isend != false
        for i = 1:k
            s += C[i]
            B[i] = s
        end
    else
        for i = 1:k
            s += C[i]
            B[i] = s - C[i]
        end
    end
end

function sais(T::AbstractVector{<:Integer}, SA::Vector{Int32}, fs::Int, n::Int, k::Int)::Vector{Int32}
    pidx = 0
    flags = 0
    C = IntVector(zeros(Int, k), 0)
    if k <= fs
        B = IntVector(SA, n + fs - k)
        flags = 1
    else
        B = IntVector(zeros(Int, k), 0)
        flags = 3
    end
    
    # stage 1
    getcounts(T, C, n, k)
    getbuckets(C, B, k, true)
    for i = 1:n
        SA[i] = 0
    end
    b = -1
    i = j = n
    m = 0
    c0 = c1 = T[n]
    i -= 1
    while 1 <= i && ((c0 = T[i]) >= c1)
        c1 = c0
        i -= 1
    end
    while 1 <= i
        c1 = c0
        i -= 1
        while 1 <= i && ((c0 = T[i]) <= c1)
            c1 = c0
            i -= 1
        end
        if 1 <= i
            0 <= b && (SA[b+1] = j)
            b = (B[c1+1] -= 1)
            j = i - 1
            m += 1
            c1 = c0
            i -= 1
            while 1 <= i && ((c0 = T[i]) >= c1)
                c1 = c0
                i -= 1
            end
        end
    end
    if 1 < m
        LMSsort(T, SA, C, B, n, k)
        name = LMSpostproc(T, SA, n, m)
    elseif m == 1
        SA[b+1] = j + 1
        name = 1
    else
        name = 0
    end
    # stage 2
    if name < m
        newfs = n + fs - 2m
        if flags & (1 | 4 | 8) == 0
            if (k + name) <= newfs
                newfs -= k
            else
                flags |= 8
            end
        end
        j = 2m + newfs
        for i = (m+(n>>1)):-1:(m+1)
            if SA[i] != 0
                SA[j] = SA[i] - 1
                j -= 1
            end
        end
        RA = IntVector(SA, m + newfs)
        sais(RA, SA, newfs, m, name)

        i = n
        j = 2m
        c0 = c1 = T[n]
        while 1 <= (i -= 1) && ((c0 = T[i]) >= c1)
            c1 = c0
        end
        while 1 <= i
            c1 = c0
            while 1 <= (i -= 1) && ((c0 = T[i]) <= c1)
                c1 = c0
            end
            if 1 <= i
                SA[j] = i
                j -= 1
                c1 = c0
                while 1 <= (i -= 1) && ((c0 = T[i]) >= c1)
                    c1 = c0
                end
            end
        end
        for i = 1:m
            SA[i] = SA[m + SA[i] + 1]
        end
        if flags & 4 != 0
            C = B = IntVector(zeros(Int, k), 0)
        end
        if flags & 2 != 0
            B = IntVector(zeros(Int, k), 0)
        end
    end
    # stage 3
    flags & 8 != 0 && getcounts(T, C, n, k)
    if 1 < m
        getbuckets(C, B, k, true)
        i = m - 1
        j = n
        p = SA[m]
        c1 = T[p+1]
        while true
            c0 = c1
            q = B[c0+1]
            while q < j
                j -= 1
                SA[j+1] = 0
            end
            while true
                j -= 1
                SA[j+1] = p
                i -= 1
                i < 0 && break
                p = SA[i+1]
                c1 = T[p+1]
                c1 != c0 && break
            end
            i < 0 && break
        end
        while 0 < j
            j -= 1
            SA[j+1] = 0
        end
    end
    induceSA(T, SA, C, B, n, k)
    return SA
end

function LMSsort(
    T::AbstractVector{<:Integer},
    SA::Vector{Int32},
    C::IntVector,
    B::IntVector,
    n::Int,
    k::Int,
)
    C == B && getcounts(T, C, n, k)
    getbuckets(C, B, k, false)
    j = n - 1
    c1 = T[j+1]
    b = B[c1+1]
    j -= 1
    SA[b+1] = T[j+1] < c1 ? ~j : j
    b += 1
    for i = 1:n
        if 0 < (j = SA[i])
            if (c0 = T[j+1]) != c1
                B[c1+1] = b
                c1 = c0
                b = B[c1+1]
            end
            j -= 1
            SA[b+1] = T[j+1] < c1 ? ~j : j
            b += 1
            SA[i] = 0
        elseif j < 0
            SA[i] = ~j
        end
    end
    C == B && getcounts(T, C, n, k)
    getbuckets(C, B, k, true)
    c1 = 0
    b = B[c1+1]
    for i = n:-1:1
        if 0 < (j = SA[i])
            c0 = T[j+1]
            if c0 != c1
                B[c1+1] = b
                c1 = c0
                b = B[c1+1]
            end
            j -= 1
            b -= 1
            SA[b+1] = T[j+1] > c1 ? ~(j + 1) : j
            SA[i] = 0
        end
    end
end

function LMSpostproc(T::AbstractVector{<:Integer}, SA::Vector{Int32}, n::Int, m::Int)
    i = 1
    while (p = SA[i]) < 0
        SA[i] = ~p
        i += 1
    end
    if i - 1 < m
        j = i
        i += 1
        while true
            if (p = SA[i]) < 0
                SA[j] = ~p
                j += 1
                SA[i] = 0
                j - 1 == m && break
            end
            i += 1
        end
    end

    i = j = n
    c0 = c1 = T[n]
    while 1 <= (i -= 1) && ((c0 = T[i]) >= c1)
        c1 = c0
    end
    while 1 <= i
        c1 = c0
        while 1 <= (i -= 1) && ((c0 = T[i]) <= c1)
            c1 = c0
        end
        if 1 <= i
            SA[m + (i >> 1) + 1] = j - i
            j = i + 1
            c1 = c0
            while 1 <= (i -= 1) && ((c0 = T[i]) >= c1)
                c1 = c0
            end
        end
    end
    name = 0
    q = n
    qlen = 0
    for i = 1:m
        p = SA[i]
        plen = SA[m + (p >> 1) + 1]
        diff = true
        if plen == qlen && (q + plen < n)
            j = 0
            while j < plen && T[p+j+1] == T[q+j+1]
                j += 1
            end
            j == plen && (diff = false)
        end
        if diff != false
            name += 1
            q = p
            qlen = plen
        end
        SA[m + (p >> 1) + 1] = name
    end
    return name
end

function induceSA(
    T::AbstractVector{<:Integer},
    SA::Vector{Int32},
    C::IntVector,
    B::IntVector,
    n::Int,
    k::Int,
)
    C == B && getcounts(T, C, n, k)
    getbuckets(C, B, k, false)
    j = n - 1
    c1 = T[j+1]
    b = B[c1+1]
    SA[b+1] = 0 < j && T[j] < c1 ? ~j : j
    b += 1
    for i = 1:n
        j = SA[i]
        SA[i] = ~j
        if 0 < j
            j -= 1
            if (c0 = T[j+1]) != c1
                B[c1+1] = b
                c1 = c0
                b = B[c1+1]
            end
            SA[b+1] = 0 < j && T[j] < c1 ? ~j : j
            b += 1
        end
    end
    C == B && getcounts(T, C, n, k)
    getbuckets(C, B, k, true)
    c1 = 0
    b = B[c1+1]
    for i = n:-1:1
        if 0 < (j = SA[i])
            j -= 1
            c0 = T[j+1]
            if c0 != c1
                B[c1+1] = b
                c1 = c0
                b = B[c1+1]
            end
            b -= 1
            SA[b+1] = j == 0 || T[j] > c1 ? ~j : j
        else
            SA[i] = ~j
        end
    end
end

function lcparray(a::Alignment, sa::Vector{Int32}, ra::Vector{Int32})
    n = length(ra)
    lcparr = similar(ra)
    h = 0
    for i in 1:n
        if ra[i] == 1
            continue
        end
        j = sa[ra[i]-1]
        maxh = n - max(i, j)
        while h <= maxh && a[i+h] == a[j+h]
            h += 1
        end
        lcparr[ra[i]] = h
        h = max(h-1, 0)
    end
    lcparr[1] = 0
    lcparr
end