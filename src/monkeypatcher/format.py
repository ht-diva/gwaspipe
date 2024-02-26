import pandas as pd
import os
import gzip

from gwaslab.io_to_formats import _check_indel, _adjust_position, _output_bed_like, _process_vcf_header, _bgzip_tabix_md5sum, _configure_output_cols_and_args, _configure_ssf_meta, md5sum_file, calculate_md5sum_file
from gwaslab.bd_common_data import get_formats_list
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_format_dict
from gwaslab.g_Log import Log
from gwaslab.io_preformat_input import print_format_info


def tofmt(sumstats,
          meta,
          path=None,
          suffix=None,
          fmt=None,
          cols=[],
          xymt_number=False,
          xymt=["X","Y","MT"],
          build="19",
          chr_prefix="",
          ssfmeta=False,
          md5sum=False,
          bgzip=False,
          tabix=False,
          tabix_indexargs={},
          verbose=True,
          no_status=False,
          log=Log(),
          to_csvargs=None):
    
    print("+++++++++++++CUSTOM TOFMT+++++++++++++++++")
    
    if to_csvargs is None:
        to_csvargs=dict()
    
    if fmt in ["ssf"]: 
        xymt_number=True
        if "SNPID" in sumstats.columns:
            log.write(' -Replacing SNPID separator from ":" to "_"...')
            sumstats["SNPID"] = sumstats["SNPID"].str.replace(":","_")
    log.write(" -Start outputting sumstats in "+fmt+" format...")
    
    if "CHR" in sumstats.columns:
        if xymt_number is False and pd.api.types.is_integer_dtype(sumstats["CHR"]):
            sumstats["CHR"]= sumstats["CHR"].map(get_number_to_chr(xymt=xymt,prefix=chr_prefix))
        elif chr_prefix is not None:
            sumstats["CHR"]= chr_prefix + sumstats["CHR"].astype("string")

    ####################################################################################################################
    if fmt=="bed":
        # bed-like format, 0-based, 
        # first 3 columns : chromosome, start, end
        # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
        is_snp,is_insert,is_delete = _check_indel(sumstats,log,verbose)
        log.write(" -formatting to 0-based bed-like file...")
        log.write(" -format description: {}".format("https://genome.ucsc.edu/FAQ/FAQformat.html#format1"))
        
        sumstats = _adjust_position(sumstats, fmt, is_snp, is_insert, is_delete, log, verbose )

        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"] + cols

        _output_bed_like(sumstats,  path, "bed", suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose)
    ####################################################################################################################   
    elif fmt=="vep":
        # bed-like format, 1-based
        # first 6 columns : chromosome, start, end, allele, strand, identifier
        # https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html
        
        is_snp,is_insert,is_delete = _check_indel(sumstats,log,verbose)
            
        log.write(" -formatting to 1-based bed-like file (for vep)...")
        log.write(" -format description: {}".format("http://asia.ensembl.org/info/docs/tools/vep/vep_formats.html"))
        sumstats = _adjust_position(sumstats, fmt, is_snp, is_insert, is_delete , log, verbose)
        
        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"]+ cols

        _output_bed_like(sumstats,  path,"vep", suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose)

    ####################################################################################################################
    elif fmt=="annovar":
        # bed-like format, 1-based, 
        # first 3 columns : Chromosome ("chr" prefix is optional), Start, End, Reference Allelel, Alternative Allele
        # https://annovar.openbioinformatics.org/en/latest/user-guide/input/
        is_snp,is_insert,is_delete = _check_indel(sumstats,log,verbose)

        log.write(" -formatting to 1-based bed-like file...")
        log.write(" -format description: {}".format("https://annovar.openbioinformatics.org/en/latest/user-guide/input/"))
        
        sumstats = _adjust_position(sumstats, fmt, is_snp, is_insert, is_delete, log, verbose )

        ouput_cols=["CHR","START","END","NEA_out","EA_out","SNPID"]+ cols
        
        _output_bed_like(sumstats, path, fmt, suffix, ouput_cols,to_csvargs,bgzip, tabix, tabix_indexargs, md5sum, log, verbose)
    
    ####################################################################################################################       
    elif fmt=="vcf":
        # GWAS-VCF
        log.write(" -"+fmt+" format will be loaded...",verbose=verbose)
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True, skip_meta_records=["format_fixed_header","format_contig_19","format_contig_38"])
        
        # determine which ID to use
        if "rsID" in sumstats.columns:
            rename_dictionary["rsID"]="ID"
        else:
            rename_dictionary["SNPID"]="ID" 
        
        # get the columns to output
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        ouput_cols = ouput_cols +["STATUS"]+ cols
        sumstats = sumstats[ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        
        # replace : with _
        sumstats["ID"] = sumstats["ID"].str.replace(":","_")
        
        # process Allele frequency data
        if "AF" in sumstats.columns:
            sumstats["INFO"] = "AF="+sumstats["AF"].astype("string")
        else:
            sumstats["INFO"] = "."

        # sumstats columns placed in vcf-FORMAT column     
        output_format=[]
        for i in sumstats.columns:
            if i in meta_data["format_format"]:
                output_format.append(i)  

        # determine path
        path = path + "."+suffix
        
        
        vcf_header =  _process_vcf_header(sumstats, meta, meta_data, build, log, verbose)

        log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)
        # output header
        with open(path,"w") as file:
            file.write(vcf_header)
        
        with open(path,"a") as file:
            log.write(" -Output columns:"," ".join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]]))
            file.write("\t".join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]])+"\n")
            log.write(" -Outputing data...")
            QUAL="."
            FILTER="PASS"
            for index,row in sumstats.iterrows():
                CHROM=str(row["#CHROM"])
                POS=str(row["POS"])
                ID=str(row["ID"])
                REF=str(row["REF"])
                ALT=str(row["ALT"])
                INFO=str(row["INFO"])
                FORMAT=":".join(output_format)
                DATA=":".join(row[output_format].astype("string"))
                file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, DATA))
        _bgzip_tabix_md5sum(path, fmt, bgzip, md5sum, tabix, tabix_indexargs, log, verbose)
    
    ####################################################################################################################
    elif fmt in get_formats_list():      
        # tabular 
        log.write(" -"+fmt+" format will be loaded...",verbose=verbose)
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True)
        
        ymal_path = path + "."+suffix+".tsv-meta.ymal"
        path = path + "."+suffix+".tsv.gz"
        log.write(" -Output path:",path, verbose=verbose) 
        
        sumstats,to_csvargs = _configure_output_cols_and_args(sumstats, rename_dictionary, cols, no_status, path, meta_data, to_csvargs, log, verbose)
        
        log.write(" -Writing sumstats to: {}...".format(path),verbose=verbose)

        try:
            fast_to_csv(sumstats, path, compress=True)
        except:
            sumstats.to_csv(path, index=None,**to_csvargs)

        if md5sum == True: 
            md5_value = md5sum_file(path,log,verbose)
        else:
            md5_value = calculate_md5sum_file(path)
        
        ## update ssf-style meta data and export to yaml file
        _configure_ssf_meta(sumstats, fmt, ssfmeta, meta, meta_data, path, md5_value, ymal_path, log, verbose)
        
        return sumstats  
    

def fast_to_csv(dataframe, path, compress=True):
        if path.endswith(".gz"):
            tsv_path = path[:-3]
        else:
            tsv_path = path

        # np.savetext() is faster than df.to_csv, however it loops through the rows of X and formats each row individually:
        # https://github.com/numpy/numpy/blob/d35cd07ea997f033b2d89d349734c61f5de54b0d/numpy/lib/npyio.py#L1613
        # We can speed up the process building the whole format string and then appling the formatting in one single call
        with open(tsv_path, 'w') as f:
            f.write(' '.join(dataframe.columns) + '\n') # header
            fmt = ' '.join(['%s']*dataframe.shape[1]) # build formatting for one single row
            fmt = '\n'.join([fmt]*dataframe.shape[0]) # add newline and replicate the formatting for all rows
            data = fmt % tuple(dataframe.to_numpy().ravel()) # flatten the array and then apply formatting
            f.write(data + "\n")

        # gzip tsv
        if compress:
            with open(tsv_path, 'rb') as file:
                bindata = bytearray(file.read())
                with gzip.open(tsv_path+".gz", "wb", compresslevel=1) as f:
                    f.write(bindata)

            # delete tsv
            os.remove(tsv_path)