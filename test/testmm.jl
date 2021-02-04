
using Test
import Test: @test

include("../src/UInt8Utf8.jl")
include("../src/MMappedString.jl")

import .MappedString: ASCII, MMappedString

function gen_ascii(n)
    rand(UInt8(0):UInt8(127), n)
end

function gen_latin1(n)
    rand(UInt8(128):UInt8(255), n)
end

function gen_unicode(n)
    v = rand(Char, n)
    String(v)
end


function rrange(m::AbstractString)
    rrange(length(m))
end
function rrange(m::Int)
    l, h = rand(1:m, 2)
    if l > h
        h, l = l, h
    end
    return l, h
end

function randindex(m::Int, n::Int)
    rand(1:m, n)
end

function run_ascii(n=2000)

    m = MMappedString{ASCII}("")
    @test length(m) == 0

    v = gen_ascii(n)
    m = MMappedString{ASCII}(v)

    s = String(m)

    @test length(m) == length(s)

    @test Base.isvalid(s)
    @test Base.isvalid(m)

    @test m == s
    @test s == m

    for (t1, t2) in zip(s, m)
        @test t1 === t2
    end
    @test length(m) == length(s)
    @test ncodeunits(m) == ncodeunits(s)

    l1 = ncodeunits(s) + 1
    @test thisind(m, 0) == 0 && thisind(s, 0) == 0
    @test thisind(m, l1) == l1
    @test thisind(s, l1) == l1
    @test firstindex(m) == firstindex(s)
    @test lastindex(m) == lastindex(s)

    for i in randindex(length(m) - 1, 100)
        # only for ascii
        @test thisind(m, i) == i
        @test getindex(m, i) == getindex(s, i)
        @test nextind(m, i) == nextind(s, i)
        @test prevind(m, i) == prevind(s, i)
        @test thisind(m, i) == thisind(s, i)
        @test firstindex(m) == firstindex(s)
        @test lastindex(m) == lastindex(s)
        @test isvalid(m, i)
    end
    # test ranges, SubString
    for p in 1:1000
        l, h = rrange(m)

        s1 = SubString(s, l, h)
        m1 = SubString(m, l, h)

        @test m[l:h] == m1
        @test s1 == m1
        @test s[l:h] == m[l:h]
        @test cmp(s1, m1) == 0

        for (t1, t2) in zip(s1, m1)
            @test t1 == t2
        end
        @test findfirst(m, m1) == findfirst(s, s1)
        @test findlast(m, m1) == findlast(s, s1)
        @test occursin(m, m1) === occursin(s, s1)
    end

    m2 = MMappedString{ASCII}(s)
    @test m2 == m && m2 == s

    for p in 1:1000
        l, h = rrange(100)
        v1, v2 = gen_ascii(l), gen_ascii(h)
        m1, m2 = MMappedString{ASCII}(v1), MMappedString{ASCII}(v2)
        s1, s2 = String(m1), String(m2)

        @test ( m1 < m2 && s1 < s2 ) || (m1 > m2 && s1 > s2) || (m1 == m2 && s1 == s2)
    end


end

function run_unicode(n=2000)

    m = MMappedString("")
    @test length(m) == 0

    s = gen_unicode(n)
    m = MMappedString(s)
    
    not_ascii = any(map(x -> x >= 0x80, m.ptr))

    if not_ascii
        @test !Base.isvalid(MMappedString{ASCII}(m.ptr))
    end 

    @test Base.isvalid(s) && Base.isvalid(m)
    @test m == s
    @test cmp(m, s) == 0
    @test length(m) == length(s)
    @test ncodeunits(m) == ncodeunits(s)
    l1 = ncodeunits(s) + 1
    @test thisind(m, 0) == 0 && thisind(s, 0) == 0
    @test thisind(m, l1) == l1
    @test thisind(s, l1) == l1
    @test firstindex(m) == firstindex(s)
    @test lastindex(m) == lastindex(s)

    for i in randindex(length(m), 100)
        idx = thisind(s, i)
        @test getindex(m, idx) == getindex(s, idx)
        @test nextind(m, i) == nextind(s, i)
        @test prevind(m, i) == prevind(s, i)
        @test thisind(m, i) == thisind(s, i)
        @test isvalid(m, idx)
        @test codeunit(m, i) == codeunit(s, i)

    end
    
    s2 = String(m)
    @test s2 == m && s2 == s
    
    for (t1, t2) in zip(s, m)
        @test t1 === t2
    end
    
    for p in 1:1000
        l, h = rrange(m)
        # get valid indexes
        l, h = thisind(s, l), thisind(s, h)
        @test length(m, l, h) == length(s, l, h)

        s1 = SubString(s, l, h)
        m1 = SubString(m, l, h)
        @test s1 == m1
        @test cmp(s1, m1) == 0

        @test s[l:h] == m[l:h]
        @test m[l:h] == m1
        for (t1, t2) in zip(s1, m1)
            @test t1 == t2
        end
        @test findfirst(m, m1) == findfirst(s, s1)
        @test findlast(m, m1) == findlast(s, s1)
        @test occursin(m, m1) === occursin(s, s1)
    end

    for p in 1:1000
        l, h = rrange(100)
        s1, s2 = gen_unicode(l), gen_unicode(h)
        m1, m2 = MMappedString(s1), MMappedString(s2)
        t1 = s1 < m1 
        c1 = cmp(s1, m1)
        @test t1 && c1 < 0 || !t1 && c1 >= 0


        @test cmp(s1, m1) == 0 && s1 == m1

        @test ( m1 < m2 && s1 < s2 ) || (m1 > m2 && s1 > s2) || (m1 == m2 && s1 == s2)
    end
    
end

function run_latin1(n=100)
    m = MMappedString{ASCII}(gen_latin1(n))
    s = m[1:end] # MMappedString[n:m] coughs up Characters
    @test isvalid(s)
    @test s == String(s)
    # MMappendString can hold latin1 bytes 
    # but they don't play well with Strings
    @test m != s # bitwise comparison latin1 != utf-8
    @test m == String(m) # String(m) is now bitwise the same but **invalid* as utf8
    @test !isvalid(String(m))
    @test m == SubString(m, 1, length(m)) # substring is fine
    @test s == SubString(m, 1, length(m)) # bitwise comparison?
    @test length(m) == length(s)
    # comparison character by charater is fine
    for (a, b) in zip(m, s)
        @test a == b
    end
end
function run_cmp(n=100)
    sl = [gen_unicode(rand(1:n)) for i in 1:20]
    ml = sl .|> MMappedString
    for i in 1:length(sl)
        for j in 1:length(sl)
            @test cmp(sl[i], sl[j]) == cmp(ml[i], ml[j])
            @test cmp(sl[i], ml[j]) == cmp(sl[i], sl[j])
        end
    end
end

function runall()
    for i in 1:50
        run_ascii()
        run_unicode()
        run_latin1()
        run_cmp()
    end
end

runall() 