function check_orientation(annotations::Vector{SFF_Model}, refstrands::Dict{String,Char})::Bool
    strandmatches = 0
    for genemodel in annotations
        if haskey(refstrands, genemodel.gene) && refstrands[genemodel.gene] == '+'
            strandmatches += 1
        end
    end
    strandmatches > length(annotations) / 2 ? true : false
end

function rotate!(
    target::FwdRev{CircularSequence},
    result::ChloeAnnotation,
    fstart::Integer
)::Tuple{FwdRev{CircularSequence},ChloeAnnotation}
    gl = result.target_length
    rstart = gl - fstart + 2
    new_target = FwdRev(rotate(target.forward, fstart), rotate(target.reverse, rstart))
    for gene in result.annotation.forward, f in gene.features
        f.feature.start = mod1(f.feature.start - fstart + 1, gl)
    end
    for gene in result.annotation.reverse, f in gene.features
        f.feature.start = mod1(f.feature.start - rstart + 1, gl)
    end
    new_target, result
end

function flip!(
    target::FwdRev{CircularSequence},
    result::ChloeAnnotation,
    fliprange::UnitRange{<:Integer}
)::Tuple{FwdRev{CircularSequence},ChloeAnnotation}
    gl = result.target_length
    if length(fliprange) == gl
        # flip whole sequence and its annotations
        new_annotations = FwdRev(result.annotation.reverse, result.annotation.forward)
        for a in new_annotations.forward
            a.strand = '+'
        end
        for a in new_annotations.reverse
            a.strand = '-'
        end
        target = FwdRev(target.reverse, target.forward)
        result = ChloeAnnotation(result.target_id, result.target_length, result.coverages, new_annotations)
    else
        # handle f and r seqs independently in case they differ because of editing
        new_fseq =
            LongDNA{2}() *
            target.forward[1:fliprange.start-1] *
            BioSequences.reverse_complement(target.forward[fliprange]) *
            target.forward[fliprange.stop+1:gl]
        new_fseq_mask = CircularMask(
            vcat(
                target.forward.mask[1:fliprange.start-1],
                reverse(target.forward.mask[fliprange]),
                target.forward.mask[fliprange.stop+1:gl]
            )
        )
        new_rseq =
            LongDNA{2}() *
            target.reverse[1:rc(fliprange.stop + 1, gl)] *
            BioSequences.reverse_complement(target.reverse[rc(fliprange, gl)]) *
            target.reverse[rc(fliprange.start - 1, gl):gl]
        new_rseq_mask = CircularMask(
            vcat(
                target.reverse.mask[1:rc(fliprange.stop + 1, gl)],
                reverse(target.reverse.mask[rc(fliprange, gl)]),
                target.reverse.mask[rc(fliprange.start - 1, gl):gl]
            )
        )
        target = FwdRev(CircularSequence(new_fseq, new_fseq_mask), CircularSequence(new_rseq, new_rseq_mask))
        flipped_forward_annotations =
            filter(x -> length(circularintersect(gene_span(x), fliprange, gl)) > 0, result.annotation.forward) #should correctly handle genes that start or finish in flipped region
        for a in flipped_forward_annotations
            a.strand = '-'
            for f in a.features
                f.feature.start = mod1(rc(fliprange.stop, gl) + f.feature.start - fliprange.start, gl)
            end
        end
        flipped_reverse_annotations =
            filter(x -> length(circularintersect(gene_span(x), rc(fliprange, gl), gl)) > 0, result.annotation.reverse)
        for a in flipped_reverse_annotations
            a.strand = '+'
            for f in a.features
                f.feature.start = mod1(fliprange.start + f.feature.start - rc(fliprange.stop, gl), gl)
            end
        end
        new_fwd_annotations = append!(
            filter(x -> x ∉ flipped_forward_annotations, result.annotation.forward),
            flipped_reverse_annotations
        )
        new_rev_annotations = append!(
            filter(x -> x ∉ flipped_reverse_annotations, result.annotation.reverse),
            flipped_forward_annotations
        )
        result = ChloeAnnotation(
            result.target_id,
            result.target_length,
            result.coverages,
            FwdRev(new_fwd_annotations, new_rev_annotations)
        )
    end
    target, result
end

function IR_range(annotation::Vector{SFF_Model})::Union{Nothing,UnitRange{Int}}
    IR_idx = findfirst(x -> startswith(x.gene, "IR"), annotation)
    if isnothing(IR_idx)
        return nothing
    end
    IR = first(annotation[IR_idx].features).feature
    IR.start:IR.start+IR.length-1
end

function transform!(
    target::FwdRev{CircularSequence},
    result::ChloeAnnotation,
    templates::Dict{String,FeatureTemplate}
)::Tuple{FwdRev{CircularSequence},ChloeAnnotation}
    gl = result.target_length
    refstrands = Dict{String,Char}()
    for t in values(templates)
        refstrands[first(split(t.path, "/"))] = t.reference_strand
    end
    IRa_range = IR_range(result.annotation.forward)
    IRb_range = IR_range(result.annotation.reverse)
    if isnothing(IRa_range)
        # no IRs so check orientation of whole genome
        if ~check_orientation(result.annotation.forward, refstrands)
            # flip whole sequence and its annotations
            new_annotations = FwdRev(result.annotation.reverse, result.annotation.forward)
            for a in new_annotations.forward
                a.strand = '+'
            end
            for a in new_annotations.reverse
                a.strand = '-'
            end
            target = FwdRev(target.reverse, target.forward)
            result = ChloeAnnotation(result.target_id, result.target_length, result.coverages, new_annotations)
        end
        ## rotate IR-less cp genomes to end with ndhF
        ndhFidx = findfirst(x -> x.gene == "ndhF", result.annotation.forward)
        if ~isnothing(ndhFidx)
            ndhF = result.annotation.forward[ndhFidx].features[1].feature
            target, result = rotate!(target, result, ndhF.start + ndhF.length)
        end
        ## rotate nuclear rDNA repeats to start with 18S rRNA gene
        rrn18idx = findfirst(x -> x.gene == "18SRRNA", result.annotation.forward)
        if ~isnothing(rrn18idx)
            fstart = result.annotation.forward[rrn18idx].features[1].feature.start
            target, result = rotate!(target, result, fstart)
        end
    else
        # In the target genome, pick the IR that has the most genes on the same strand as in the templates and make this IR1
        IRa_annotations =
            filter(x -> length(circularintersect(gene_span(x), IRa_range, gl)) > 0, result.annotation.forward)
        if check_orientation(IRa_annotations, refstrands)
            # IRa is in correct orientation, so rotate to start from the nucleotide following IRb; no flip required
            target, result = rotate!(target, result, gl - IRb_range.start + 2)
        else
            # IRb is in correct orientation, so rotate to end of IRa, then flip the IR annotations
            target, result = rotate!(target, result, IRa_range.stop + 1)
            #flip IR annotations
            IRa_feature =
                first(
                    result.annotation.forward[findfirst(
                        x -> startswith(x.gene, "IR"),
                        result.annotation.forward
                    )].features
                ).feature
            IRb_feature =
                first(
                    result.annotation.reverse[findfirst(
                        x -> startswith(x.gene, "IR"),
                        result.annotation.reverse
                    )].features
                ).feature
            tmp = IRa_feature.start
            IRa_feature.start = rc(IRb_feature.start + IRb_feature.length - 1, gl)
            IRb_feature.start = rc(tmp + IRa_feature.length - 1, gl)
        end
        # reload ranges as they may have altered
        IRa_range = IR_range(result.annotation.forward)
        IRb_range = IR_range(result.annotation.reverse)
        # orient LSC to maximise strand agreement with the templates
        LSC_range = 1:IRa_range.start-1
        LSC_annotations =
            filter(x -> length(circularintersect(gene_span(x), LSC_range, gl)) > 0, result.annotation.forward)
        if ~check_orientation(LSC_annotations, refstrands)
            target, result = flip!(target, result, LSC_range)
        end
        # orient SSC to maximise strand agreement with the templates
        SSC_range = IRa_range.stop+1:rc(IRb_range.stop + 1, gl)
        SSC_annotations =
            filter(x -> length(circularintersect(gene_span(x), SSC_range, gl)) > 0, result.annotation.forward)
        if ~check_orientation(SSC_annotations, refstrands)
            target, result = flip!(target, result, SSC_range)
        end
        #rename IRs?
    end
    target, result
end
