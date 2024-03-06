
struct ChloeBoost
    coding_xgb_model::XGBoost.Booster
    noncoding_xgb_model::XGBoost.Booster
    xgb_coding_model::XGBoost.Booster
    ChloeBoost() = new(
        XGBoost.load(XGBoost.Booster, joinpath(@__DIR__, "coding_xgb.model")),
        XGBoost.load(XGBoost.Booster, joinpath(@__DIR__, "noncoding_xgb.model")),
        XGBoost.load(XGBoost.Booster, joinpath(@__DIR__, "xgb.coding.model")))
end

# function show_handles(msg)
#     boost = get_boost()
#     mm = "$(msg)[tid=$(Threads.threadid())]: $(boost.coding_xgb_model.handle) and $(boost.noncoding_xgb_model.handle) and $(boost.xgb_coding_model.handle)"
#     @info mm
#     # println(mm)
# end

function get_boost()::ChloeBoost
    # if length(boost) == 0
    #     push!(boost, ChloeBoost())
    # end
    __boost[1]
end

const __boost::Vector{ChloeBoost} = []

# see Julia module initialisation
# https://docs.julialang.org/en/v1/manual/modules/#Module-initialization-and-precompilation
# for when __init__() is called.

# originally we had "global" variables like:

# const xgb_coding_model = XGBoost.Booster(XGBoost.DMatrix[], model_file=joinpath(@__DIR__, "xgb.coding.model"))

# but were getting errors in Chloe when using from a module viz:

#     nested task error: XGBoostError: (caller: XGBoosterPredictFromDMatrix)
#     [14:45:22] /workspace/srcdir/xgboost/src/c_api/c_api.cc:1059: Booster has not been initialized or has already been disposed.

# probably from precompilation not retaining the xgboost C-code handles....
function __init__()
    # @info "__init__ called"
    push!(__boost, ChloeBoost())
end
