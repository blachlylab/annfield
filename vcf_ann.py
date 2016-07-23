import unittest

class TestAnnDecode(unittest.TestCase):
    # Test annotations
    ANN1 = "C|stop_gained|HIGH|BTK|ENSG...|transcript|featid|Coding|7|c.123T>C|p.L456R|123|234|456|-100|W1"
    ANN2 = "T|histone_binding_site|LOW|BTK|ENSG...|transcript|featid|Coding|7/10|c.123T>C|p.L456R|123/1000|234/900|456/500|-100|WARNING_REF_DOES_NOT_MATCH_GENOME"
    ANN3 = "A|intron_variant&nc_transcript_variant|MEDIUM|||transcript|featid|Coding|7/10|c.123T>C|p.L456R|123|234|456|-100|W1"
    ANN_multiple = ANN1 + ", " + ANN3
    ANN_weq = "ANN=" + ANN1

    # Results
    RES1 = {'distance_to_feature': '-100', 'gene_name': 'BTK', 'feature_type': 'transcript', 'protein_position': '456', 'errors_warnings_info': 'W1', 'feature_id': 'featid', 'rank_total': '7', 'gene_id': 'ENSG...', 'hgvs.p': 'p.L456R', 'transcript_biotype': 'Coding', 'cds_position': '234', 'cdna_position': '123', 'annotation': 'stop_gained', 'putative_impact': 'HIGH', 'allele': 'C', 'hgvs.c': 'c.123T>C'}
    RES2 = {'putative_impact': 'LOW', 'errors_warnings_info': 'WARNING_REF_DOES_NOT_MATCH_GENOME', 'cdna_position': '123/1000', 'hgvs.p': 'p.L456R', 'distance_to_feature': '-100', 'rank_total': '7/10', 'protein_position': '456/500', 'feature_id': 'featid', 'cds_position': '234/900', 'annotation': 'histone_binding_site', 'allele': 'T', 'transcript_biotype': 'Coding', 'gene_id': 'ENSG...', 'feature_type': 'transcript', 'hgvs.c': 'c.123T>C', 'gene_name': 'BTK'}
    RES3 = {'feature_id': 'featid', 'protein_position': '456', 'gene_id': '', 'rank_total': '7/10', 'transcript_biotype': 'Coding', 'hgvs.c': 'c.123T>C', 'distance_to_feature': '-100', 'feature_type': 'transcript', 'hgvs.p': 'p.L456R', 'errors_warnings_info': 'W1', 'cds_position': '234', 'annotation': 'intron_variant&nc_transcript_variant', 'putative_impact': 'MEDIUM', 'cdna_position': '123', 'gene_name': '', 'allele': 'A'}

    def test_get_ann_empty(self):
        self.assertEqual(get_annotations(""), [])

    def test_decode_basic(self):
        self.assertEqual(decode(self.ANN1), self.RES1)
    
    def test_decode_range(self):
        self.assertEqual(decode(self.ANN2), self.RES2)

    def test_decode_multiple_conseq(self):
        self.assertEqual(decode(self.ANN3), self.RES3)

    def test_decode_multiple_effects(self):
        self.assertEqual(get_annotations(self.ANN_multiple), [self.RES1, self.RES3] )

    def test_decode_weq(self):
        self.assertEqual(get_annotations(self.ANN_weq), [self.RES1])

        

def get_annotations(ann_string):
    """
    Retreives one or more annotations from VCF ANN field standard 1.0
    http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf

    input: annotation string (either with or without ANN= prefix)
    output: LIST of dictionaries representing annotations
        absent elements in the ANN field will not have
        a corresponding key in the dict
    """

    results = list()
    if not ann_string: return results

    # if prefixed with ANN=, trim this
    if ann_string[0:4] == "ANN=": ann_string = ann_string[4:]

    # specs page 2:
    # Multiple effects / consequences are separated by comma.
    ann_list = [x.strip() for x in ann_string.split(',')] # remove whitespace around comma
    for eff_string in ann_list:
        results.append( decode(eff_string) )

    return results


def decode(eff_string):
    """
    Decodes a single VCF ANN annotation fieldset

    called by: get_annotations
    input: string (previously split on ',')
    output: dictionary representing annotations
        absent elements in the ANN field will not have
        a corresponding key in the dict
    """


    # specs page 1:
    # Data fields are encoded separated by pipe sign "|";
    # the order of fields is written in the VCF header.
    eff = eff_string.split('|')

    result = dict()

    # TODO: read field order from header. For now, define statically
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


    for i,value in enumerate(eff):
        result[ field_list[i] ] = value

    # Done, return object
    #print(result)
    return result
