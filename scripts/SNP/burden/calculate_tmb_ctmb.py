#!/mnt/bin/TMB/2.0.0/python3
"""Calculate tissue TMB only (frameshift burden included), using VEP VCF or hotfilter TSV."""

import argparse
import gzip
import logging
import os
import pandas as pd

all_consequences = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
    "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant","splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant","splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
    "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant","upstream_gene_variant", "downstream_gene_variant",
    "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_variant", "regulatory_region_ablation", "regulatory_region_amplification",
    "intergenic_variant", "sequence_variant"
}

# Inclusion list for TMB: protein-impacting or core splice consequences
include_consequences = {
    'transcript_ablation',
    'splice_acceptor_variant', 'splice_donor_variant',
    'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
    "feature_elongation", "feature_truncation",
    'inframe_insertion', 'inframe_deletion',
    'missense_variant', 'protein_altering_variant'
}

def calculate_vep(infile:str, tumor_sample:str, vaf:float=0.05):
    """Parse VEP-annotated VCF and count TMB/indel using VAF and consequence filters.

    Consequence policy (per Ensembl/VEP predicted consequences: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences):
        * Keep vep High and moderate consequence (protein-impacting and core splice acceptor/donor consequences).
        * Exclude non-coding / UTR / synonymous / distal splice-region / regulatory / upstream/downstream / TFBS / promoter, etc., as defined in `excluded_consequences`.
        * Frameshift variants additionally contribute to indel burden.

    Logic:
        * Read VCF (gz or plain), skipping meta header lines.
        * For each ALT allele, collect CSQ entries matching that allele (VEP v115 format).
        * Apply consequence filters above; fetch AF from tumor sample column (allele-specific AF in FORMAT).
        * Count as TMB if AF >= vaf and passes filters; frameshift adds to indel burden.

    Args:
        infile (str): Path to VEP-annotated VCF file (supports .vcf or .vcf.gz).
        vaf (float): Minimum variant allele frequency threshold.
        tumor_sample (str): Tumor sample column name in the VCF (used for AF extraction).

    Returns:
        tuple[int, int]: (tmb_count, indel_count)
    """
    # calculate the skip rows for VCF header; support gzipped input
    skiprows = 0
    opener = gzip.open if infile.endswith('.gz') else open
    with opener(infile, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                skiprows += 1
            else:
                break
    df = pd.read_csv(infile, skiprows=skiprows, sep='\t', keep_default_na=False, dtype=str, compression='infer')

    if tumor_sample not in df.columns:
        raise ValueError(f"Sample column '{tumor_sample}' not found in VCF; available columns: {df.columns.tolist()}")

    def _get_sample_af(format_field, sample_field, alt_idx):
        fmt_keys = format_field.split(':')
        if 'AF' not in fmt_keys:
            return None
        af_idx = fmt_keys.index('AF')
        sample_parts = sample_field.split(':')
        if af_idx >= len(sample_parts):
            return None
        af_values = sample_parts[af_idx].split(',')
        if alt_idx >= len(af_values):
            return None
        try:
            return float(af_values[alt_idx])
        except ValueError:
            return None

    def _parse_csq(info_field):
        for item in info_field.split(';'):
            if item.startswith('CSQ='):
                return item.split('=', 1)[1].split(',')
        return []

    def _allele_match(csq_allele, alt, ref):
        # VEP may encode deletions as '-' in CSQ
        if csq_allele == alt:
            return True
        if csq_allele == '-' and len(alt) < len(ref):
            return True
        return False
    # Only keep autosomes; include human (1-22/chr1-22) and mouse (1-19/chr1-19)
    autosomes = {str(i) for i in range(1, 23)} | {f'chr{i}' for i in range(1, 23)}
    autosomes |= {str(i) for i in range(1, 20)} | {f'chr{i}' for i in range(1, 20)}

    tmb, indel = 0, 0
    for _, r in df.iterrows():
        chrom = str(r['#CHROM'])
        if chrom not in autosomes:
            continue
        csq_entries = _parse_csq(r['INFO'])
        if not csq_entries:
            logging.warning('Pos({}:{}) missing CSQ annotations'.format(r['#CHROM'], r['POS']))
            continue

        alts = str(r['ALT']).split(',')
        for alt_idx, alt in enumerate(alts):
            alt_csq = [c for c in csq_entries if _allele_match(c.split('|', 1)[0], alt, str(r['REF']))]
            if not alt_csq:
                continue

            consequences = []
            frameshift_flag = False
            for entry in alt_csq:
                fields = entry.split('|')
                if len(fields) < 2:
                    continue
                cons = fields[1]
                tokens = cons.split('&') if cons else []
                consequences.extend(tokens)
                if 'frameshift_variant' in tokens:
                    frameshift_flag = True

            if not consequences:
                continue


            # Inclusion policy: require at least one token in include_consequences
            if not any(token in include_consequences for token in consequences):
                continue

            af = _get_sample_af(r['FORMAT'], r[tumor_sample], alt_idx)
            if af is None or af < vaf:
                continue

            tmb += 1
            if frameshift_flag:
                indel += 1
    return tmb, indel


def calculate_filter(infile, vaf):
    """Parse hotfilter TSV and count TMB/indel using VAF and functional filters.

    Logic:
        * Read hotfilter TSV.
        * Keep variants with AutoInterpStatus in Y/H/HN, Primary == Y, caseAF >= vaf.
        * Exclude splice/intron/promoter consequences from TMB; frameshift contributes to indel burden.

    Args:
        infile (str): Path to hotfilter TSV file.
        vaf (float): Minimum variant allele frequency threshold.

    Returns:
        tuple[int, int]: (tmb_count, indel_count)
    """
    df = pd.read_csv(infile, sep='\t', keep_default_na=False, dtype=str)
    tmb, indel = 0, 0
    for _, r in df.iterrows():
        auto, caseaf, func, prim = r['AutoInterpStatus'], float(r['caseAF']), r['Function'], r['Primary']
        if auto in ['Y', 'H', 'HN'] and caseaf >= vaf and not (
                func == 'splice-3' or func == 'splice-5' or func == 'intron' or func == 'promoter') and prim == 'Y':
            tmb += 1
        if auto in ['Y', 'H', 'HN'] and caseaf >= vaf and func == 'frameshift' and prim == 'Y':
            indel += 1
    return tmb, indel


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate tissue TMB (no cTMB).')
    parser.add_argument('--vep', help='input VEP annotated VCF')
    parser.add_argument('--filter', help='input hotfilter TSV')
    parser.add_argument('--sample', required=True, help='sample name')
    parser.add_argument('-s', '--size', type=float, default=32, help='cds size(MB), default 32MB')
    parser.add_argument('-v', '--vaf', help='min VAF filter', default=0.05, type=float)
    parser.add_argument('-o', '--oud', required=True, help='output directory')
    args = parser.parse_args()

    if not args.vep and not args.filter:
        raise ValueError('Either --vep or --filter must be provided to calculate TMB')

    logging.basicConfig(level=logging.NOTSET)
    vaf = args.vaf

    if args.vep:
        logging.info('use VEP result to calculate TMB')
        tmb, indel = calculate_vep(args.vep, vaf, args.sample)
    else:
        logging.info('use hotfilter result to calculate TMB')
        tmb, indel = calculate_filter(args.filter, vaf)

    tmbf = os.path.join(args.oud, args.sample + '_tmb.txt')
    headers = ['SampleID', 'TMB', 'TMB(Mutations/Mb)', 'Indel_Burden', 'Indel_Burden(Mutations/Mb)']
    with open(tmbf, 'w') as outfile:
        outfile.write('\t'.join(headers) + '\n')
        outfile.write('\t'.join(map(str, [args.sample, tmb, '{:.2f}'.format(tmb * 1. / args.size), indel,
                                          '{:.2f}'.format(indel * 1. / args.size)])) + '\n')
