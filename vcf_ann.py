def decode(ann_string):
"""
Decodes VCF ANN field standard 1.0
http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf

input: annotation string (either with or without ANN= prefix)
output: dictionary representing annotations
    absent elements in the ANN field will not have
    a corresponding key in the dict
"""

    if not ann_string: return dict()

    # if prefixed with ANN=, trim this
    if ann_string[0:4] == "ANN=": ann_string = ann_string[4:]

    # specs page 2:
    # Multiple “effects / consequences” are separated by comma.
    ann_list = ann_string.split(',')

    # TODO: set up loop. for now examine first only
    eff_string = ann_list[0]

    # specs page 1:
    # Data fields are encoded separated by pipe sign "|";
    # the order of fields is written in the VCF header.
    eff = eff_string.split('|')

    result = dict()

    # TODO: read field order from header. FOr now, define statically
    field_list = [  "allele",
                    "annotation",
                    "putative_impact",
                    "gene_name",
                    "gene_id",
                    "feature_type",
                    "feature_id",
                    "transcript_biotype",
                    "rank_total",
                    "hgvs.c",
                    "hgvs.p",
                    "cdna_position",
                    "cds_position",
                    "protein_position",
                    "distance_to_feature",
                    "errors_warnings_info" ]

    fields = dict()
    for i,j in enumerate(field_list):
        fields[j] = i


    for i,field in enumerate(eff):
	result[ field ] = eff[i]

    # Done, return object
    print result
    return result
