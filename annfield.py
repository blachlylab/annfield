import unittest

class TestAnnDecode(unittest.TestCase):
    # Test annotations
    ANN1 = "C|stop_gained|HIGH|BTK|ENSG...|transcript|featid|Coding|7|c.123T>C|p.L456R|123|234|456|-100|W1"
    ANN2 = "T|histone_binding_site|LOW|BTK|ENSG...|transcript|featid|Coding|7/10|c.123T>C|p.L456R|123/1000|234/900|456/500|-100|WARNING_REF_DOES_NOT_MATCH_GENOME"
    ANN3 = "A|intron_variant&nc_transcript_variant|MEDIUM|||transcript|featid|Coding|7/10|c.123T>C|p.L456R|123|234|456|-100|W1"
    ANN_multiple = ANN1 + ", " + ANN2
    ANN_weq = "ANN=" + ANN1

    # Results
    RES1 = {'distance_to_feature': '-100', 'gene_name': 'BTK', 'feature_type': 'transcript', 'protein_position': '456', 'errors_warnings_info': 'W1', 'feature_id': 'featid', 'rank_total': '7', 'gene_id': 'ENSG...', 'hgvs_p': 'p.L456R', 'transcript_biotype': 'Coding', 'cds_position': '234', 'cdna_position': '123', 'effect': 'stop_gained', 'impact': 'HIGH', 'allele': 'C', 'hgvs_c': 'c.123T>C'}
    RES2 = {'impact': 'LOW', 'errors_warnings_info': 'WARNING_REF_DOES_NOT_MATCH_GENOME', 'cdna_position': '123/1000', 'hgvs_p': 'p.L456R', 'distance_to_feature': '-100', 'rank_total': '7/10', 'protein_position': '456/500', 'feature_id': 'featid', 'cds_position': '234/900', 'effect': 'histone_binding_site', 'allele': 'T', 'transcript_biotype': 'Coding', 'gene_id': 'ENSG...', 'feature_type': 'transcript', 'hgvs_c': 'c.123T>C', 'gene_name': 'BTK'}
    RES3 = {'feature_id': 'featid', 'protein_position': '456', 'gene_id': '', 'rank_total': '7/10', 'transcript_biotype': 'Coding', 'hgvs_c': 'c.123T>C', 'distance_to_feature': '-100', 'feature_type': 'transcript', 'hgvs_p': 'p.L456R', 'errors_warnings_info': 'W1', 'cds_position': '234', 'effect': 'intron_variant&nc_transcript_variant', 'impact': 'MEDIUM', 'cdna_position': '123', 'gene_name': '', 'allele': 'A'}
    RES3A = {'feature_id': 'featid', 'protein_position': '456', 'gene_id': '', 'rank_total': '7/10', 'transcript_biotype': 'Coding', 'hgvs_c': 'c.123T>C', 'distance_to_feature': '-100', 'feature_type': 'transcript', 'hgvs_p': 'p.L456R', 'errors_warnings_info': 'W1', 'cds_position': '234', 'effect': 'intron_variant', 'impact': 'MEDIUM', 'cdna_position': '123', 'gene_name': '', 'allele': 'A'}
    RES3B = {'feature_id': 'featid', 'protein_position': '456', 'gene_id': '', 'rank_total': '7/10', 'transcript_biotype': 'Coding', 'hgvs_c': 'c.123T>C', 'distance_to_feature': '-100', 'feature_type': 'transcript', 'hgvs_p': 'p.L456R', 'errors_warnings_info': 'W1', 'cds_position': '234', 'effect': 'nc_transcript_variant', 'impact': 'MEDIUM', 'cdna_position': '123', 'gene_name': '', 'allele': 'A'}

    def test_get_ann_empty(self):
        "Test decoding of empty string"
        res = [x for x in decode("")]
        self.assertEqual(res, [])

    def test_decode_basic(self):
        "Test basic decoding on single effect, single annotation"
        res = [x for x in decode(self.ANN1)]
        self.assertEqual(res, [self.RES1])
    
    def test_decode_range(self):
        "Test pos/max in rank_total and xxx_position fields"
        res = [x for x in decode(self.ANN2)]
        self.assertEqual(res, [self.RES2])

    def test_decode_multiple_conseq(self):
        "Test splitting on '&' in the consequence field"
        res = [x for x in decode(self.ANN3)]
        self.assertEqual(res, [self.RES3A, self.RES3B])

    def test_decode_multiple_effects(self):
        "Test splitting on ','"
        res = [x for x in decode(self.ANN_multiple)]
        self.assertEqual(res, [self.RES1, self.RES2] )

    def test_decode_weq(self):
        "Test decoding when prefixed with ANN="
        res = [x for x in decode(self.ANN_weq)]
        self.assertEqual(res, [self.RES1])

        
# TODO: read field order from header. For now, define statically
# snake_case is 20% easier to read than camelCase
# http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&tp=&arnumber=5521745
field_list = [  "allele",
                "effect",           # aka annotation
                "impact",           # aka putative_impact
                "gene_name",
                "gene_id",
                "feature_type",
                "feature_id",
                "transcript_biotype",
                "rank_total",
                "hgvs_c",           # previously hgvs.c
                "hgvs_p",           # previously hgvs.p
                "cdna_position",    # should cdna_len be a sep field?
                "cds_position",     # likewise
                "protein_position", # likewise
                "distance_to_feature",
                "errors_warnings_info" ]

def decode(ann_string):
    """
    Decodes one or more annotations from VCF ANN field standard 1.0
    http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf

    input: annotation string (either with or without ANN= prefix)
    output: (Generator)
        one or more dictionaries representing annotations
        absent elements in the ANN field will not have
        a corresponding key in the dict
    """

    # if empty string
    if not ann_string: return

    # if prefixed with ANN=, trim this
    if ann_string[0:4] == "ANN=": ann_string = ann_string[4:]

    # specs page 2:
    # Multiple effects / consequences are separated by comma.
    ann_list = [x.strip() for x in ann_string.split(',')] # remove whitespace around comma
    for eff_string in ann_list:

        # specs page 1:
        # Data fields are encoded separated by pipe sign "|";
        # the order of fields is written in the VCF header.
        eff = eff_string.split('|')

        result = dict()

        # TODO: read field order from header. For now, define statically
        fields = dict()
        for i,j in enumerate(field_list):
            fields[j] = i


        for i,value in enumerate(eff):
            result[ field_list[i] ] = value

        # Look for '&' in the annotation field
        # if present, split into two separate objects
        if '&' in result['effect']:
            # Count number of '&' conjunctions
            annotations = result['effect'].split('&')
            for i in range(len(annotations)):
                copy = dict(result)
                copy['effect'] = annotations[i]
                yield copy
        # No compound effect annotation ('&')
        # So just yield result
        else:
            yield result

    

