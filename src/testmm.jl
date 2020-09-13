include("MMappedString.jl")

function gen_ascii(n)
    v = rand(0:127, n)
    MMappedString{ASCII}(v)
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
function run_ascii(n=2000)
    m = gen_ascii(n)

    s = String(m)

    @assert Base.isvalid(s)
    @assert Base.isvalid(m)

    @assert m == s
    @assert s == m

    for (t1, t2) in zip(s, m)
        @assert t1 === t2
    end
    for p in 1:1000
        l, h = rrange(m)

        s1 = SubString(s, l, h)
        m1 = SubString(m, l, h)
        @assert s1 == m1

        @assert s[l:h] == m[l:h]
        for (t1, t2) in zip(s1, m1)
            @assert t1 == t2
        end
    end

    m2 = MMappedString{ASCII}(s)
    @assert m2 == m && m2 == s

    for p in 1:1000
        l, h = rrange(1000)
        m1, m2 = gen_ascii(l), gen_ascii(h)
        s1, s2 = String(m1), String(m2)

        @assert ( m1 < m2 && s1 < s2 ) || (m1 > m2 && s1 > s2) || (m1 == m2 && s1 == s2)
    end
    println("OK $n $(length(m))")

end

function run_unicode(n=2000)
    s = gen_unicode(n)
    m = MMappedString(s)
    
    not_ascii = any(map(x -> x >= 0x80, m.ptr))

    if not_ascii
        @assert !Base.isvalid(MMappedString{ASCII}(m.ptr))
    end 

    @assert Base.isvalid(s) && Base.isvalid(m)
    @assert m == s
    @assert length(m) == length(s)
    
    s2 = String(m)
    @assert s2 == m && s2 == s
    
    for (t1, t2) in zip(s, m)
        @assert t1 === t2
    end
    
    for p in 1:1000
        l, h = rrange(m)
        l, h = thisind(s, l), thisind(s, h)

        s1 = SubString(s, l, h)
        m1 = SubString(m, l, h)
        @assert s1 == m1

        @assert s[l:h] == m[l:h]
        for (t1, t2) in zip(s1, m1)
            @assert t1 == t2
        end
    end

    for p in 1:1000
        l, h = rrange(1000)
        s1, s2 = gen_unicode(l), gen_unicode(h)
        m1, m2 = MMappedString(s1), MMappedString(s2)

        @assert ( m1 < m2 && s1 < s2 ) || (m1 > m2 && s1 > s2) || (m1 == m2 && s1 == s2)
    end

    println("OK $n $(length(m))")
    
end


for i in 1:50
    run_ascii()
    run_unicode()
end
