import Mmap
# test if String copies memory
open("test.txt", "w") do f
    write(f, "abcdefghijk\n")
end
# open readonly
f = open("test.txt")
size = filesize(f)
v = Mmap.mmap(f)
@assert length(v) === size
s = String(v)
# length of v is now zero
@assert length(v) === 0
# length is 0 but the data is still there!
for i in 1:size
    vs = unsafe_load(pointer(s), i)
    vp = unsafe_load(pointer(v), i)
    @assert vs == vp
end
println("string start...")
print(s)
# you can write to a string....
unsafe_store!(pointer(s), 66, 1)
vp = unsafe_load(pointer(s), 1)
@assert vp == 66
@assert s[1] == Char(66)
println("string now...")
print(s)
# can't write to the vector though: ReadOnlyMemoryError
# must be different address from String
unsafe_store!(pointer(v), 66, 1)
