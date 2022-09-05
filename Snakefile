import pandas as pd
import sourmash
import os

# create wildcard variables for workflow
metadata = pd.read_csv("metadata.txt", sep = "\t")            # read in metadata as a pandas dataframe
MTX = metadata['mtx_run_accession'].unique().tolist()         # make run accession in mtx col into list
MGX = metadata['mgx_run_accession'].unique().tolist()         # make run accession in mgx col into list
RUN_ACCESSIONS = MGX + MTX                                    # combine mtx and mgx into one list with all run accessions
SAMPLES = metadata['sample_name'].unique().tolist()           # make a list of sample names
KSIZES = [31]                                                 # create a list of k-mer sizes for the workflow

MTX_MINUS_MGX = [x + '-minus-' + y for x, y in zip(MTX, MGX)]
print(MTX_MINUS_MGX)

rule all:
    input:
        expand("outputs/sourmash_sketch_subtract/{sample}_k{ksize}.sig", sample = SAMPLES, ksize = KSIZES),
        expand("outputs/sourmash_sketch_subtract/{mtx_minus_mgx}-k{ksize}.sig", mtx_minus_mgx = MTX_MINUS_MGX, ksize = KSIZES)

rule sourmash_sketch:
    output: "outputs/sourmash_sketch/{run_accession}.sig"
    shell:'''
    fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {wildcards.run_accession} |
        sourmash sketch dna -p k=21,k=31,k=51,scaled=200,abund --name {wildcards.run_accession} -o {output} -
    '''

# This is super slow. It basically re-implements sourmash sig subtract but removing some loops i didn't need
rule sig_subtract_in_python_api:
    input:
        metadata = 'metadata.txt',
        sigs = expand('outputs/sourmash_sketch/{run_accession}.sig', run_accession = RUN_ACCESSIONS)
    output:
        sigs = expand("outputs/sourmash_sketch_subtract/{sample}_k{{ksize}}.sig", sample = SAMPLES)
    benchmark: "benchmark_api_full.{ksize}.txt"
    run:
        ksize = wildcards.ksize
        # read in metadata dataframe to derive sample pairs
        metadata2 = pd.read_csv("metadata.txt", sep = "\t") # read in metadata as a pandas dataframe
        metadata2 = metadata2.reset_index()  # make sure indexes pair with number of rows
        for index, row in metadata2.iterrows():
            # grab the run accessions for a given paired metagenome and metatranscriptome
            mtx_run_accession = row['mtx_run_accession']
            mgx_run_accession = row['mgx_run_accession']
            # using the run assession, create paths of the signatures
            mtx_sigfile = os.path.join('outputs/sourmash_sketch', mtx_run_accession + '.sig')
            mgx_sigfile = os.path.join('outputs/sourmash_sketch', mgx_run_accession + '.sig')
                
            # read in mtx signature and grab both the minhash object and the hashes in that object
            mtx_sigobj = sourmash.load_one_signature(mtx_sigfile, ksize=ksize, select_moltype='DNA')
            mtx_mh = mtx_sigobj.minhash
    
            # turn the mtx hashes into a set
            mtx_subtract_mins = set(mtx_mh.hashes)
            
            # read in mgx signature
            mgx_sigobj = sourmash.load_one_signature(mgx_sigfile, ksize=ksize, select_moltype='DNA')
            
            # do the subtraction
            mtx_subtract_mins -= set(mgx_sigobj.minhash.hashes)
            
            # make a new signature that only contains the mtx hashes that were not in the mgx
            mtx_subtract_mh = mtx_sigobj.minhash.copy_and_clear().flatten()
            mtx_subtract_mh.add_many(mtx_subtract_mins)
    
            # re-establish hash abundances from mtx sig
            # re-read in mtx sig just in case abundances were removed in place
            abund_sig = sourmash.load_one_signature(mtx_sigfile, ksize=ksize, select_moltype='DNA')
            mtx_subtract_mh = mtx_subtract_mh.inflate(abund_sig.minhash)
            
            # create a new sig object that can be written to file
            mtx_subtract_sigobj = sourmash.SourmashSignature(mtx_subtract_mh)
    
            # create a name that reflects the signature contents
            name = row['sample_name']
            mtx_subtract_sigobj.name = name
    
            # create output sig file name
            sig_filename = os.path.join('outputs/sourmash_sketch_subtract/' + name + "_k" + ksize + ".sig")
            # write sig to file
            with sourmash.sourmash_args.FileOutput(sig_filename, 'wt') as fp:
                sourmash.save_signatures([mtx_subtract_sigobj], fp=fp)


# This is also slow -- commented out because it produces the same output as the above rule
#rule sig_subtract_in_python_by_shell:
#    input:
#        metadata = 'inputs/metadata-paired-mgx-mtx.tsv',
#        sigs = expand('outputs/sourmash_sketch/{run_accession}.sig', run_accession = RUN_ACCESSIONS)
#    output:
#        sigs = expand("outputs/sourmash_sketch_subtract/{sample}_k{{ksize}}.sig", sample = SAMPLES)
#    benchmark: "benchmark_api_shell.txt"
#    run:
#        ksize = wildcards.ksize
#        # read in metadata dataframe to derive sample pairs
#        metadata2 = pd.read_csv("metadata.tsv", sep = "\t") # read in metadata as a pandas dataframe
#        metadata2 = metadata2.reset_index()  # make sure indexes pair with number of rows
#        for index, row in metadata2.iterrows():
#            # grab the run accessions for a given paired metagenome and metatranscriptome
#            mtx_run_accession = row['mtx_run_accession']
#            mgx_run_accession = row['mgx_run_accession']
#            # using the run assession, create paths of the signatures
#            mtx_sigfile = os.path.join('outputs/sourmash_sketch', mtx_run_accession + '.sig')
#            mgx_sigfile = os.path.join('outputs/sourmash_sketch', mgx_run_accession + '.sig')
#            # create a name that reflects the signature contents
#            name = row['sample_name']
#            out_filename = os.path.join('outputs/sourmash_sketch_subtract/' + name + "_k" + ksize + ".sig")
#            print("subtracting!")
#            shell("sourmash sig subtract -k {ksize} -o {out_filename} -A {mtx_sigfile} {mtx_sigfile} {mgx_sigfile}")


# This is not slow
rule sig_subtract_by_cli:
    input:
        mtx_sig = 'outputs/sourmash_sketch/{mtx_run_accession}.sig',
        mgx_sig = 'outputs/sourmash_sketch/{mgx_run_accession}.sig',
    output: "outputs/sourmash_sketch_subtract/{mtx_run_accession}-minus-{mgx_run_accession}-k{ksize}.sig"
    benchmark: "benchmark_cli_full.{mtx_run_accession}.{mgx_run_accession}.{ksize}.txt"
    shell:'''
    sourmash sig subtract -k {wildcards.ksize} -o {output} -A {input.mtx_sig} {input.mtx_sig} {input.mgx_sig}
    '''
