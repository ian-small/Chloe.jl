

import IntervalTrees: IntervalTree, AbstractInterval, Interval
import Base

struct FeatureInterval <: AbstractInterval{Int32}
    feature::Feature
end



# Base.convert(::Type{FeatureInterval}, f::Feature) = FeatureInterval(f)

toIntervalTree(v::AbstractArray{Feature}) = IntervalTree{Int32,FeatureInterval}(sort(map(FeatureInterval, v)))

Base.first(f::FeatureInterval) = f.feature.start
Base.last(f::FeatureInterval) = f.feature.start + f.feature.length - one(Int32)
Base.intersect(i::IntervalTree{Int32,FeatureInterval}, lo::Integer, hi::Integer) = intersect(i, Interval{Int32}(lo, hi))

