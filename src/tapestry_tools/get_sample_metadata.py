import polars as pl 

def get_sampleid_generation_gender_age_methylation(): 
    return (
        pl
        # https://docs.google.com/spreadsheets/d/1JrfEwe9fdD2HdCKd7IN4aTGkVrXJ8rBdF4gN6lLwbPs/edit 
        .read_csv('/scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/CEPH K1463 - All samples - DNA ages - deidentified.csv')
        .rename({
            '': 'id',
            'Generation?': 'generation',
            'Gender:': 'gender',
            'Age at blood draw': 'age_at_blood_draw'
        })
        .select([
            'id',
            'generation',
            'gender',
            'age_at_blood_draw'
        ])
        .with_columns(
            pl
            .format(
                "/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs/{}.dna-methylation.founder-phased.all_cpgs.sorted.bed.gz", 
                pl.col("id")
            ).alias('founder_phased_methylation')
        )
    )
