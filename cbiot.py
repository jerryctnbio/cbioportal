"""This little script extracts all the defined variations in cbioportal.

It tallies all the unique number of samples for each variation.
Each variation shall have defined 'chr', 'pos_start' and 'pos_end'.

The output shall have 'Gene', 'Chr', 'Start', 'End', 'Ref', 'Var'
, the unique number of samples

The output is in a table format, i.e., in '.xlsx' format (Excel).
"""

__version='R1.0.0.0'


import os, sys, time
import argparse
from pathlib import Path

from bravado.client import SwaggerClient
import bravado.response
import bravado.http_future

import pandas as pd


def main():
    
    title="This is a script to extract variations from cbioportal"
    usage=title+"\n\n"
    usage+="The input is a gene name or a list of gene names by -g\n"
    usage+="e.g. -g TP53 EGFR or -g TP53 -g EGFR.\n"
    usage+="or a file with a column of gene names.\n\n"
    usage+="If any gene name is given, the input file, if given, is"
    usage+=" ignored.\n\n"
    usage+="The output defaults to Excel file with all variation details,\n"
    usage+="such as gene, chr, pos_start, pos_end, reference_allele\n"
    usage+=", variation allele and unique number of samples.\n"
    usage+="But it can be in csv format by -csv\n\n"
    usage+="If the output file name is not provided, it is given as below:\n"
    usage+="If the input file is given, it will mirror that with '.out'\n"
    usage+="added; otherwise it is named as 'your.cbioportal.out'.\n\n"
    usage+="The output defaults to Excel (.xslx) format. You need to install\n"
    usage+="the following modules for it to work:\n\n"
    usage+="            xlwt openpyxl xlsxwriter xlrd\n\n"
    
    argParser=argparse.ArgumentParser(description=title)
    
    argParser.add_argument("-g", help="input gene name[s]"
                               , action="append", nargs='+')
    
    msg="input file with a column of gene names"
    argParser.add_argument("-f", help=msg)
    argParser.add_argument("-o", help="output file name")

    msg="set output format to '.csv', otherwise '.xlsx'."
    argParser.add_argument("-csv", help=msg
                                 , action='store_true')
    argParser.add_argument("-t", '--timeout', default=300.0, type=float
                    , help="default timeout in seconds [300.0]")
    
    args=argParser.parse_args()
    
    if not args.g and not args.f:
        
        print("*****************************************")
        print("Please enter at least one gene name, see usage below:\n\n")
        print(f"{usage}\n")
        argParser.print_help()
        
        sys.exit(1)
    
    # get the gene names in a list
    gene=getGene(args)
    
    print("The input gene names are:")
    print(gene)

    timeout=args.timeout

    # get the variations
    variation=getVariation(gene, timeout)

    # get the output file name
    outFile=getOutFileName(args, args.csv)

    # save the output
    if args.csv:
        # save to .csv file
        variation.to_csv(outFile)
    else:
        # save to Excel file (.xlsx)
        variation.to_excel(outFile)


def getVariation(gene:list, timeout=300.0)->list:
    
    # connect to cbioportal
    url="https://www.cbioportal.org/api/api-docs"
    cbioportal=SwaggerClient.from_url(url, config={"validate_requests":False
                                                ,"validate_responses":False})    

    # get the gene entrez IDs.
    geneID, geneID_display=getGeneInfo(cbioportal, gene)
    
    print(f"The gene entrez IDs are:")
    print(geneID_display)

    # how many studies
    operation=cbioportal.Studies.getAllStudiesUsingGET
    args={}
    study=cbioportal_exe(operation, args, timeout)

    print(f"There are {len(study)} studies.")
    
    # how many cancer types
    operation=cbioportal.Cancer_Types.getAllCancerTypesUsingGET
    args={}
    cancer_type=cbioportal_exe(operation, args, timeout)

    print(f"There are {len(cancer_type)} cancer types.")
    
    # how many molecular profiles
    operation=cbioportal.Molecular_Profiles.getAllMolecularProfilesUsingGET
    args={}
    molecular_profile=cbioportal_exe(operation, args, timeout)

    print(f"There are {len(molecular_profile)} molecular profiles.")
    
    # how many sampleLists
    operation=cbioportal.Sample_Lists.getAllSampleListsUsingGET
    args={}
    sample_list=cbioportal_exe(operation, args, timeout)
    
    print(f"There are {len(sample_list)} sample lists.")

    # get all the samples
    operation=cbioportal.Samples.getSamplesByKeywordUsingGET
    args={}
    sample=cbioportal_exe(operation, args, timeout)
    
    print(f"There are {len(sample)} samples.")

    variationH=getVariationHardWay(cbioportal, molecular_profile
                                   , sample_list, geneID, timeout)

    variationE=getVariationEasyWay(cbioportal, molecular_profile
                                             , geneID, timeout)

    # merge the two DataFrames
    on_t=['Gene', 'Chr', 'Start', 'End', 'Ref', 'Var']
    
    variation_combined=pd.merge(variationH, variationE, on=on_t, how='outer')

    print(variation_combined)

    variation_combined['Diff']=variation_combined.UniqueSample_x \
                               -variation_combined.UniqueSample_y

    print(variation_combined)
    
    by_t=['Gene', 'Diff', 'Chr', 'Start', 'End', 'Ref', 'Var', 'UniqueSample_x']
    asc_t=[True, False, True, True, True, True, True, True]

    variation_combined.sort_values(by=by_t, inplace=True, ascending=asc_t)

    return variation_combined


def getVariationEasyWay(cbioportal, molecular_profile:list
                                             , geneID:list, timeout=300.0):
    """Get all the variations in one scoop.
    
    This method is easy, runs a little faster but may choke the system
    if there are many hits
    
    Parameters:
    cbioportal        : object --- handler to cbioportal server
    molecular_profile : list   --- list of objects for molecular profiles 
    geneID            : list   --- Entrez gene IDs
    timeout           : float  --- maximum time in seconds to retrieve data
                                   default:300
                                   
    Return:
    A pandas DataFrame. The columns are the necessary info for defining a
                        variant ['Chr', 'Start', 'End', 'Ref' and 'Var']
                        , plus ['Gene', 'UniqueSample'], and the number of
                        unique samples.
    """

    # get the molecular profile IDs    
    mpfID=[m.molecularProfileId for m in molecular_profile]
        
    filter_t={"entrezGeneIds":geneID, "molecularProfileIds":mpfID}
    
    operation=cbioportal.Mutations \
                       .fetchMutationsInMultipleMolecularProfilesUsingPOST
    args=dict(mutationMultipleStudyFilter=filter_t, projection='DETAILED')
    mutation_all_list=cbioportal_exe(operation, args, timeout)
        
    print(f"There are {len(mutation_all_list)} mutations returned.")
                    
    variation=getUniqueVariation(mutation_all_list)
            
    return variation  

    
def getVariationHardWay(cbioportal, molecular_profile:list
                         , sample_list:list, geneID:list, timeout=300.0):
    """Extract the variations in a batch mode.
    
    First generate all the combinations of molecular profile and sample
    list IDs based on study ID. Then loop through all the combinations    

    Parameters:
    cbioportal        : object --- handler to cbioportal server
    molecular_profile : list   --- list of objects for molecular profiles 
    sample_list       : list   --- list of objects for sample lists
    geneID            : list   --- Entrez gene IDs
    timeout           : float  --- maximum time in seconds to retrieve data
                                   default:300
                                   
    Return:
    A pandas DataFrame. The columns are the necessary info for defining a
                        variant ['Chr', 'Start', 'End', 'Ref' and 'Var']
                        , plus ['Gene', 'UniqueSample'], and the number of
                        unique samples.

    """
        
    # only use molecular profiles ending with '_mutations'
    # as others may crash. Reasons not known yet as of June 30, 2020.  
    mpf=pd.DataFrame([{'StudyID':m.studyId, 'Mol_ProID':m.molecularProfileId}
                                    for m in molecular_profile
                            if m.molecularProfileId.endswith("_mutations")])

    # only use sample list ending with '_all'. This should be enough
    sl=pd.DataFrame([{'StudyID':s.studyId, 'Sample_ListID':s.sampleListId}
                            for s in sample_list
                            if s.sampleListId.endswith("_all")])

    # get all the combinations of molecular profiles and sample lists
    mpf_sl=pd.merge(mpf, sl)
    
    msg_t=" combinations of molecular profile and sample list."
    print(f"There are {len(mpf_sl)}{msg_t}")  
    
    # now loop through all the combinations
    mutation_all_list=[]
    for i in range(len(mpf_sl)):

        print(f"doing {i+1} now.")
        
        mpfID=mpf_sl.iloc[i].Mol_ProID
        sampleListID=mpf_sl.iloc[i].Sample_ListID

        print(f"Molecular profile ID is {mpfID}.")
        print(f"Sample list ID is {sampleListID}.")
        
        operation=cbioportal.Mutations \
                            .fetchMutationsInMolecularProfileUsingPOST
        filter_t={"entrezGeneIds":geneID, "sampleListId":sampleListID}

        args=dict(molecularProfileId=mpfID, mutationFilter=filter_t
                                          , projection='DETAILED')
        mutation=cbioportal_exe(operation, args, timeout)

        print(f"There are {len(mutation)} mutations returned.")
        
        mutation_all_list+=mutation
        
        # a short delay so not to choke the system
        time.sleep(0.01)

    variation=getUniqueVariation(mutation_all_list)
            
    return variation  
        

def getUniqueVariation(mutation_all_before:list)->pd.DataFrame:
    """Get the unique variations.
    
    The key here is to use the attribute 'uniqueSampleKey' of the mutation
    object.
    
    Parameters:
    mutation_all_before   :  list --- a list of mutations  
    
    Return:
    A pandas DataFrame  The columns are the necessary info for defining a
                        variant ['Chr', 'Start', 'End', 'Ref' and 'Var']
                        , plus ['Gene', 'UniqueSample'], and the number of unique
                        samples.

    """

    # keep only the minimum info to simplify the mutations,
    # and filter out the ones without defined coordinates
    mutation_all_list=[{'Chr':m.chr
                      , 'Start':m.startPosition
                      , 'End':m.endPosition
                      , 'Ref':m.referenceAllele
                      , 'Var':m.variantAllele
#                      , 'StudyID':m.studyId
#                      , 'MolProfID':m.molecularProfileId
#                      , 'SampleID':m.sampleId
                      , 'UniqueSample':m.uniqueSampleKey
                      , 'Gene':m.gene.hugoGeneSymbol}
                for m in mutation_all_before 
                  if (m.startPosition !=-1 and m.endPosition !=-1)]
    
    # turn the list into a pandas DataFrame        
    mutation_all_pd=pd.DataFrame(mutation_all_list)

    print(mutation_all_pd)
    
    # group by gene, mutation identity
    groupByItem=['Gene', 'Chr', 'Start', 'End', 'Ref', 'Var']
        
    mutation_by_gene_mutationID=mutation_all_pd \
                               .groupby(groupByItem)['UniqueSample'] \
                               .nunique().to_frame()
    
    print(mutation_by_gene_mutationID)
        
    nTotal=mutation_by_gene_mutationID.UniqueSample.sum()
    print(f"There are {nTotal} unique samples.")
    
    return mutation_by_gene_mutationID    
    

def getGeneInfo(cbioportal, gene:list, timeout=300.0):
    """Get the Entrez gene ID for all the genes."""
    
    operation=cbioportal.Genes.fetchGenesUsingPOST
    args_t={'geneIdType':'HUGO_GENE_SYMBOL', 'geneIds':gene}
    geneInfo=cbioportal_exe(operation, args_t, timeout)

    geneInfo_dict={}
    for g in geneInfo:        
        geneInfo_dict[g.hugoGeneSymbol]=g.entrezGeneId
        
    geneID=[]
    for g in gene:
        
        try:
            geneID.append(geneInfo_dict[g.upper()])

        except KeyError:
            
            print(f"*** Gene '{g}' does not have a valid Entrez gene ID. ***")
            
            sys.exit(1)

    # The info displayed to the terminal, matching the input order
    geneID_display=["{} ==> {}"
                     .format(g, geneInfo_dict[g.upper()]) for g in gene]
                     
    return (geneID, geneID_display)

    
def cbioportal_exe(operation, args, timeout=300.0):
    """wrapper function to all cbioportal API request.
    
    Parameters:
    operation  :  object      --- callable function to access cbioportal
    args       :  dictionary  --- arguments to the operation 
    timeout    :  float       --- maximum time in seconds to retrieve data
                                  default:300

    Note:
    If it takes longer than the 'timeout' time, the call will fail. A
    message will say so. You may want to increase the time.   
    
    The response is returned only if the status code is '200'.
    
    Return:  a list of objects
    """
    
    try:
        cbio_response=operation(**args).response(timeout=timeout)
    
    except bravado.http_future.BravadoTimeoutError:
        print(f"Request timed out!")
        print(f"Please increase the timeout by '-t' or '--timeout'.")

        sys.exit(1)

    except Exception as e:
        
        msg="Error occurred extracting data from cbioportal!\n"
        print(f"{msg}")

        sys.exit(1)

    # check if the request worked    
    if isinstance(cbio_response, bravado.response.BravadoResponse) and \
                   (cbio_response.incoming_response.status_code==200):
        cbio_result=cbio_response.result
    else:
        cbio_result=''

    return cbio_result
    
            
def getOutFileName(args, csv:bool):
    """Get the output file name."""

    ext='.csv' if csv else '.xlsx'

    if args.o and len(args.o)>0:
        
        parent=Path(args.o).parent
        stem=Path(args.o).stem
    
        return Path(parent/(stem+ext))

    # if not given, generate as follows:
    if args.f and len(args.f)>0:

        parent=Path(args.f).parent
        stem=Path(args.f).stem
        
        stem=stem+".out"
        outFile=Path(parent/(stem+ext))
        
    else:
        outFile=Path("your.cbioportal.out"+ext)

    return outFile


def getGene(args)->list:
    """Get all the gene names in a list."""
    
    gene_t=[]
    if args.g:
        
        # input file is ignored
        args.f=""
        
        # combine the lists
        for g in args.g:
            gene_t+=g
        
    if args.f and len(args.f)>0:
 
        try:
            with open(args.f, 'r') as fh:
                f_st=fh.read()
        except Exception as e:
            
            print(f"*** Error reading input file '{args.f}'! ***")

            sys.exit(1)
            
        gene_t=f_st.rstrip("\n").split("\n")

    return gene_t


if __name__ == '__main__':
    """ script execution point"""
    
    main()
    
