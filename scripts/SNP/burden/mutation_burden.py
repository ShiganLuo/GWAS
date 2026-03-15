import os
import pysam
import pandas as pd
import bisect


def compute_nonoverlap_exon_length(gtf_file: str, main_chroms: list, exclude_sex: bool = True) -> int:
    """Return merged exon length (bp) from GTF; can drop sex chromosomes if requested."""
    chrom_set = set(main_chroms)
    if exclude_sex:
        sex_chroms = {'X', 'Y', 'chrX', 'chrY'}
        chrom_set = {c for c in chrom_set if c not in sex_chroms}

    exon_dict = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5 or fields[2] != 'exon':
                continue
            chrom = fields[0]
            if chrom not in chrom_set:
                continue
            start = int(fields[3])
            end = int(fields[4])
            exon_dict.setdefault(chrom, []).append((start, end))

    def merge(intervals):
        if not intervals:
            return []
        intervals.sort(key=lambda x: x[0])
        merged = [intervals[0]]
        for st, ed in intervals[1:]:
            pst, ped = merged[-1]
            if st <= ped:
                merged[-1] = (pst, max(ped, ed))
            else:
                merged.append((st, ed))
        return merged

    total_len = 0
    for chrom in exon_dict:
        merged = merge(exon_dict[chrom])
        total_len += sum(ed - st + 1 for st, ed in merged)
    return total_len


def compute_nonoverlap_cds_length(gtf_file: str, main_chroms: list, exclude_sex: bool = True) -> int:
    """Return merged CDS length (bp) from GTF; can drop sex chromosomes if requested."""
    chrom_set = set(main_chroms)
    if exclude_sex:
        sex_chroms = {'X', 'Y', 'chrX', 'chrY'}
        chrom_set = {c for c in chrom_set if c not in sex_chroms}

    cds_dict = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5 or fields[2] != 'CDS':
                continue
            chrom = fields[0]
            if chrom not in chrom_set:
                continue
            start = int(fields[3])
            end = int(fields[4])
            cds_dict.setdefault(chrom, []).append((start, end))

    def merge(intervals):
        if not intervals:
            return []
        intervals.sort(key=lambda x: x[0])
        merged = [intervals[0]]
        for st, ed in intervals[1:]:
            pst, ped = merged[-1]
            if st <= ped:
                merged[-1] = (pst, max(ped, ed))
            else:
                merged.append((st, ed))
        return merged

    total_len = 0
    for chrom in cds_dict:
        merged = merge(cds_dict[chrom])
        total_len += sum(ed - st + 1 for st, ed in merged)
    return total_len

def compute_exonic_mutation_burden_from_vcf_dir(
        vcf_dir: str,
        gtf_file: str,
        qual_filter: float = 30,
        main_chroms: list = ["chr" + str(i) for i in range(1,23)] + ["chrX","chrY","chrM"],
        chunksize: int = 100000
    ) -> pd.DataFrame:
    """
    Functions: compute exonic mutation burden from a directory of VCF files, using GTF annotation to define exons.
    Parameters:
        - vcf_dir: vcf files directory, can contain .vcf, .vcf.gz or _SNV.txt.gz files
        - gtf_file: GTF annotation file for exon definitions
        - qual_filter: minimum QUAL score to consider a mutation (default 30)
        - main_chroms: list of main chromosomes to consider (default human chr1-22,X,Y,M)
        - chunksize: number of lines to read at a time for non-standard VCF text files (default 100000)
    Returns:
        - DataFrame with columns: sample, mutation_count, TMB
    """

    # ---------------------------------------------------------------
    # 1) 构建外显子区间，并预处理二分查找所需结构
    # ---------------------------------------------------------------
    def build_exon_intervals(gtf_file, main_chroms):
        exon_dict = {}
        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "exon":
                    continue
                chrom = fields[0]
                if chrom not in main_chroms:
                    continue
                start = int(fields[3])
                end = int(fields[4])
                exon_dict.setdefault(chrom, []).append((start, end))

        # 合并区间
        def merge(intervals):
            if not intervals:
                return []
            intervals.sort(key=lambda x: x[0])
            merged = [intervals[0]]
            for st, ed in intervals[1:]:
                pst, ped = merged[-1]
                if st <= ped:
                    merged[-1] = (pst, max(ped, ed))
                else:
                    merged.append((st, ed))
            return merged

        for c in exon_dict:
            merged = merge(exon_dict[c])
            # 生成便于二分查找的起点和终点列表
            starts = [s for s, e in merged]
            ends = [e for s, e in merged]
            exon_dict[c] = {"intervals": merged, "starts": starts, "ends": ends}

        return exon_dict

    exon_dict = build_exon_intervals(gtf_file, main_chroms)
    total_exon_length = sum(e - s + 1 for chrom in exon_dict for s, e in exon_dict[chrom]["intervals"])

    # ---------------------------------------------------------------
    # 2) 判断 POS 是否落在外显子区间（用二分法）
    # ---------------------------------------------------------------
    def is_pos_in_exon(chrom, pos):
        if chrom not in exon_dict:
            return False
        starts = exon_dict[chrom]["starts"]
        ends = exon_dict[chrom]["ends"]
        idx = bisect.bisect_right(starts, pos) - 1
        if idx >= 0 and pos <= ends[idx]:
            return True
        return False

    # ---------------------------------------------------------------
    # 3) 计算外显子突变数
    # ---------------------------------------------------------------
    def count_exonic_mutations(vcf_file, qual_filter, main_chroms, chunksize=100000):
        count = 0
        try:
            # 标准VCF
            vcf = pysam.VariantFile(vcf_file)
            vcf_chroms = set(vcf.header.contigs)
            valid_chroms = set(main_chroms) & vcf_chroms
            for chrom in valid_chroms:
                for rec in vcf.fetch(chrom):
                    if rec.qual is not None and rec.qual < qual_filter:
                        continue
                    if is_pos_in_exon(chrom, rec.pos):
                        count += 1
        except:
            # 非标准VCF文本，分块读取
            if vcf_file.endswith(".gz"):
                file_like = vcf_file
                compression = "gzip"
            else:
                file_like = vcf_file
                compression = None

            reader = pd.read_csv(
                file_like,
                sep=None,
                engine="python",
                compression=compression,
                chunksize=chunksize
            )

            required_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

            for chunk in reader:
                if not all(c in chunk.columns for c in required_cols):
                    raise ValueError(f"{vcf_file} 缺少必要列: {required_cols}")
                chunk = chunk[chunk["CHROM"].isin(main_chroms)]
                # chunk = chunk[chunk["QUAL"].notna() & (chunk["QUAL"] >= qual_filter)]

                for chrom, group in chunk.groupby("CHROM"):
                    for pos in group["POS"]:
                        if is_pos_in_exon(chrom, pos):
                            count += 1

        return count

    # ---------------------------------------------------------------
    # 4) 批量处理目录
    # ---------------------------------------------------------------
    results = []
    vcf_files = [f for f in os.listdir(vcf_dir)
                 if f.endswith(".vcf") or f.endswith(".vcf.gz") or f.endswith("_SNV.txt.gz")]

    for vcf_name in vcf_files:
        vcf_path = os.path.join(vcf_dir, vcf_name)
        mut_count = count_exonic_mutations(
            vcf_path, qual_filter, main_chroms, chunksize
        )
        mut_burden = mut_count / (total_exon_length / 1e6)
        results.append({
            "sample": vcf_name.replace(".vcf.gz","").replace(".vcf","").replace("_SNV.txt.gz",""),
            "mutation_count": mut_count,
            "TMB": mut_burden
        })

    df = pd.DataFrame(results)
    return df



if __name__ == "__main__":
    main_chroms_mouse = ["chr" + str(i) for i in range(1,20)] + ["chrX","chrY","chrM"]
    gtf = "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf"
    exon_length = compute_nonoverlap_exon_length(gtf, main_chroms_mouse, exclude_sex=True)
    print(f"Total non-overlapping exon length (excluding sex chromosomes): {exon_length} bp")
    cds_length = compute_nonoverlap_cds_length(gtf, main_chroms_mouse, exclude_sex=True)
    print(f"Total non-overlapping CDS length (excluding sex chromosomes): {cds_length} bp")