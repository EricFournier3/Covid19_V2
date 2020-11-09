#!/usr/local/bin/python
# vim:tabstop=4 shiftwidth=4 softtabstop=4 expandtab

'''
Created on 11 mars 2013

modified by Eric Fournier 2020-10-02

@author: sam
@author: Eric Fournier
'''

import inspect
import os , sys , argparse, re
import glob
import json
from subprocess import Popen, PIPE
from time import time
from datetime import datetime

ARGS            = {}       # Command line arguments dictionary
PROCESSING_DIR  = ""
TECHNO_LIST     = []

LRUNHEADER      = []
LPLATEHEADER    = []
# Status 
    # RAW   = no processing (or no analysis)
    # PROC  = Processed but not referenced 
    # REF   = Processed and Referenced



def parseoptions( ):

        parser = argparse.ArgumentParser( description=" " )
        parser.add_argument( '-d',  '--dir',  help="Path to Beluga COVID_full_processing directory",   required=True )
        parser.add_argument( '-p',  '--fplate', help="File with list of plates and samples" )
        parser.add_argument( '-o',  '--repo', help="Path to Repository output", required=True  )

        #Eric Fournier 2020-10-02 add here
        parser.add_argument('-v','--verbose',help="Print output message in logfile",action='store_true')

        global ARGS
        ARGS = parser.parse_args()

def main() :
    t1    = time()
    print ("\nBegin @ " + str( datetime.now() ))

    parseoptions( )
   
    global logfile
    global runMessageFile
    global missingConsensusLog
    global consensusListFile
    global missingConsPercNLog
    global hiddenSampleFile
    global LRUNHEADER
    global ILLUMINA_OLD_RUNS
    global LPLATEHEADER
    global PROCESSING_DIR   # full_processing dir
    global TECHNO_LIST
    #Eric Fournier 2020-10-02 add here
    global TRACE_DIR
    global verbose
    verbose = ARGS.verbose

    #Eric Fournier 2020-10-02 add here
    if not os.path.isdir(ARGS.repo):
        os.makedirs(ARGS.repo)
 
    #Eric Fournier 2020-10-02 change here
    fjson   = open("runParam_v2.json")
    pjson   = json.load( fjson )
    LRUNHEADER      = pjson[ 'LRUNHEADER' ] 
    LPLATEHEADER    = pjson[ 'LPLATEHEADER']
    TECHNO_LIST     = pjson[ 'TECHNO_LIST' ]
    PROCESSING_DIR  = ARGS.dir
    TRACE_DIR = pjson['TRACE_DIR']
    ILLUMINA_OLD_RUNS = pjson['ILLUMINA_OLD_RUNS']
    fjson.close()

    cwd = os.getcwd()
    
    hiddenSampleFile_path = os.path.join(TRACE_DIR,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_hiddenSamples.list")
    hiddenSampleFile = open(hiddenSampleFile_path,'w')

    logfile_path = os.path.join(TRACE_DIR,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') +".log")
    logfile = open(logfile_path,'w')
    
    runMessageFile_path = os.path.join(TRACE_DIR,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_runMessage" +".txt") 
    runMessageFile = open(runMessageFile_path,'w')

    missingConsensusLog_path = os.path.join(TRACE_DIR,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_MissingConsensus" +".txt")
    missingConsensusLog = open(missingConsensusLog_path,'w')

    consensusListFile_path = os.path.join(TRACE_DIR,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_consensusList" +".list")
    consensusListFile = open(consensusListFile_path,'w')
    consensusListFile.write("SAMPLE\tSTATUS\tPATH\tTECHNO\tPERC_N\tRUN_NAME\n")

    missingConsPercNLog_path = os.path.join(TRACE_DIR,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_MissingConsPercN" +".txt")
    missingConsPercNLog = open(missingConsPercNLog_path,'w')


    if verbose:
        logfile.write("***************** in updateRunFile ********************\n")
        
    drun    = updateRunFile( )
    # Read list of plates ; This list has to be provided by LSPQ
    dplate  = {}    # key : plate ; value : list of samples
    if ARGS.fplate :
        if verbose:
            logfile.write("\n***************** in readPlateFile ********************\n")

        dplate = readPlateFile( ARGS.fplate )
    if verbose:
        logfile.write("\n***************** in updateRepo ********************\n")

    updateRepo( drun, dplate,  ARGS.repo )
    os.chdir( cwd )

    if verbose:
        logfile.write("\n***************** in writeRunFile *******************\n")

    writeRunFile( drun )

    dt = time()-t1
    print (" ---- time(s) = " + str(dt) + "\n")

def updateRepo( drun, dplate, repodir ) :
    cwd = os.getcwd()   
    # key : samplekey ; value : list of col
    dsample = {}

    # For each run
    for path,lrun in drun.items() :
        os.chdir( cwd )
            
        rundir  = os.path.join( ARGS.dir , path )
        techno  = lrun[1]
        date    = lrun[0]
        run     = lrun[2]
        status  = lrun[4]
        if verbose:
            logfile.write("\n     Update for:\n\t\t" + rundir + "\n\t\t" + date + "\n\t\t" + run + "\n\t\t" + run + "\n\t\t" + status + "\n")

        lsampletot  = lrun[6]
        lrmess      = []        # List of run messages

        # if run is  already processed (not only raw data)
        if not status == "PROC" :
            #print("\tStatus = " + status + ", run not processed" )
            pass

        # Run processed
        else :
            
            # Processing directory
            procdir = ""
            # nanopore path
            if techno == "nanopore" :
                dtmp    = ""
                procdirtmp  = os.path.join( rundir , "analysis" )
                try :
                    ld          = os.listdir( procdirtmp ) # Sometimes, more than one analysis dir
                except : 
                    #print ("\tERR: permission denied for directory = " + procdirtmp )
                    mess = "\tERR: permission denied for directory = " + procdirtmp
                    #Eric Fournier 2020-10-15 desactive cat line b existe pas
                    #line[ 4 ] = "ERR"
                    lrun[5] += ";" + mess
                    
                else :
                    dtmp        = ld[0]     # By default, choose first one
                    if len(ld) > 1 :    # if more than one ...
                        for d in ld :
                            if 'nanopolish_800x' in d :     # ... take the one named nanopolish...
                                dtmp = d
                    procdir = os.path.join( procdirtmp , dtmp )
        
            # Illumina and mgi path
            else :
                if techno == "illumina":
                    if run  in ILLUMINA_OLD_RUNS:
                        procdir = os.path.join( rundir , "alignment")
                    else:
                        procdir = os.path.join( rundir , "consensus")
                else:
                    procdir = os.path.join( rundir , "alignment")

            # Processing directory
            if len( procdir ) == 0 or not  os.path.isdir( procdir ) :
                mess = "ERR: Processing directory not found"
                #print("\t" + mess)
                lrun[5] += ";" + mess

            else :
                
                # Check if we have permission on processing directory
                try :
                    os.chdir( procdir )
                    # filter only directories
                    lfproc      = os.listdir( "." )
                
                #except PermissionError:
                except :
                    mess = "ERR: permission denied for directory = " + procdir 
                    #print("\t" + mess)
                    lrun[5] += ";" + mess

                #Permission to read ok
                else :
                    # Find list of samples processed in run
                    lsample    = []
                    for f in lfproc :
                        if os.path.isdir( f ) and not "multiqc" in f :
                            lsample.append( f )
                    lrun[7] = lsample
                    # Nb of samples in run
                    #print( "\t%i samples total, %i samples to process, %s" %(len(lsampletot), len(lsample), procdir))
                    s = 0
                    
                    # For all the samples in processing directory
                    for sname in lsample :
                        s += 1
                        
                        # Relative path to processing dir => NO : abs path
                        # relativepath    = os.path.join( "../../", procdir)     # Because we add the "plate" level
                        sampledirsrc    = os.path.join( procdir , sname )
                        
                        # Sample info
                        #!!!!! SNAME : nom de l'echantillon utilise par McGill ; 
                        #!!!!! SAMPLE : echantillon selon nommage LSPQ
                        plate,sample       = plateForSample( sname, dplate )
                        #print "PLATE " + plate + " sample " + sample
                        #print("%i dir=%s ; sample=%s ; "%(s, sample,sname) )

                        # On suppose que samplekey est unique
                        # car le meme echantillon n'est jamais sequence 2x avec la meme techno le meme jour
                        samplekey   = getSampleKey( plate , sample , techno , date )
                        # Name of plate used by McGill
                        plateMcGill = getMcGillPlate( run , techno )

                        #If dir for plate does not exist
                        platedirdest    = os.path.join( repodir , plate )
                        #Eric Fournier 2020-11-09 on desactive les 3 lignes suivantes
                        #ltmp            = glob.glob( platedirdest + "*" ) # Check if exists
                        #if len( ltmp ) == 0 :
                        #    os.mkdir( platedirdest )
                        #Eric Fournier 2020-11-09 on ajoute
                        if not os.path.isdir(platedirdest):
                            os.mkdir(platedirdest)
                        

                        # if repo for sample does not exist
                        sampledirdest   = os.path.join( platedirdest , samplekey ) 
                        # check if sampledir already exists
                        ltmp        = glob.glob( sampledirdest + "*" )

                        # No sample dir
                        if len( ltmp ) == 0 : 
                            # Create repository
                            os.mkdir( sampledirdest )
                           
                            consensus_access = ""
                            save_consensus_path = False
                            # Link to consensus seq
                            if techno == "nanopore" :
                                lconsensus  = glob.glob( os.path.join( sampledirsrc ,  "*consensus.nanopore.*.fasta" ) )
                            elif techno == "mgi":
                                lconsensus  = glob.glob( os.path.join( sampledirsrc ,  "*consensus.gisaid_renamed.fa" ) )
                            else :
                                if run  in ILLUMINA_OLD_RUNS:
                                    lconsensus  = glob.glob( os.path.join( sampledirsrc ,  "*consensus.gisaid_renamed.fa" ) )
                                else:
                                    lconsensus  = glob.glob( os.path.join( sampledirsrc ,  "*.consensus.illumina.pass.fasta" ) )
                                    if len(lconsensus) == 0:
                                        lconsensus  = glob.glob( os.path.join( sampledirsrc ,  "*.consensus.illumina.flag.fasta" ) )
                                        if len(lconsensus) == 0:
                                            lconsensus  = glob.glob( os.path.join( sampledirsrc ,  "*.consensus.illumina.rej.fasta" ) )
                            if len( lconsensus ) > 0 :
                                save_consensus_path = True
                                os.symlink( lconsensus[0] , os.path.join( sampledirdest , sample + ".consensus.fasta" ) )
                            else:
                                if(os.access(sampledirsrc,os.R_OK)):
                                    consensus_access = "ACCESS_TRUE"
                                else:
                                    consensus_access = "ACCESS_FALSE"
                                missingConsensusLog.write("Missing consensus for " + sample + " in " + sampledirdest + " -> " + sampledirsrc + "  " + consensus_access +  "\n\n")

                            # Link to bam file (primer trimmed)
                            lbam        = glob.glob( os.path.join( sampledirsrc ,  "*primer*.bam" ) )
                            if len( lbam ) > 0 :
                                os.symlink( lbam[0] , os.path.join( sampledirdest , sample + ".primertrim.bam" ) )
                            lbami        = glob.glob( os.path.join( sampledirsrc ,  "*primer*.bam.bai" ) )
                            if len( lbami ) > 0 :
                                os.symlink( lbami[0] , os.path.join( sampledirdest , sample + ".primertrim.bam.bai" ) )
                          
                            # Link to vcf major variants
                            if techno == "nanopore" :
                                lvcfmaj     = glob.glob( os.path.join( sampledirsrc ,  "*.pass.vcf" ) )
                            else :
                                lvcfmaj     = glob.glob( os.path.join( sampledirsrc ,  "*.primerTrim.vcf" ) )
                            if len(lvcfmaj) > 0 :
                                os.symlink( lvcfmaj[0] , os.path.join( sampledirdest , sample + ".major.vcf" ) )
                            
                            # Link to vcf minor variants
                            if techno == "nanopore" :
                                lvcfmin     = glob.glob( os.path.join( sampledirsrc ,  "*.merged.vcf" ) )
                            else :
                                # MISSING for illumina and mgi
                                lvcfmin     = []
                            if len(lvcfmin) > 0 :
                                os.symlink( lvcfmin[0] , os.path.join( sampledirdest , sample + ".minor.vcf" ) )

                            # Link to raw fastq file
                            if techno == "nanopore" :
                                lfastqraw   =  glob.glob( os.path.join( sampledirsrc ,  "*barcode*.fastq" ) )
                                if len( lfastqraw ) > 0 :
                                    os.symlink( lfastqraw[0] , os.path.join( sampledirdest , sample + ".raw.fastq" ) )
                            else :
                                # Eric Fournier 2020-09-21 change here
                                if techno == "illumina":
                                    lfastqraw   =  glob.glob( os.path.join( rundir ,  "data/" + sname + "*.fastq.gz" ) ) 
                                elif techno == "mgi":
                                    lfastqraw   =  glob.glob( os.path.join( rundir ,  "data/L0*/" + sname + "*.fastq.gz" ) ) 
                                for raw in lfastqraw :
                                    #Eric Fournier 2020-09-21 change here
                                    if "_R1_" in raw  or "R1." in raw:
                                        os.symlink( raw , os.path.join( sampledirdest , sample + ".raw_R1.fastq.gz" ) )
                                    elif "_R2_" in raw or "R2." in raw: 
                                        os.symlink( raw , os.path.join( sampledirdest , sample + ".raw_R2.fastq.gz" ) )
                            
                            # Link to clean fastq file
                            if techno == "nanopore" :
                                lfastqclean     =  glob.glob( os.path.join( sampledirsrc ,  "*clean.fastq" ) )
                                if len( lfastqclean ) > 0 :
                                    os.symlink( lfastqclean[0] , os.path.join( sampledirdest , sample + ".clean.fastq" ) )
                            else :
                                lfastqclean     = glob.glob( os.path.join( rundir ,  "trim/" + sname + "/*.fastq.gz" ) )
                                for clean in lfastqclean :
                                    #Eric Fournier 2020-09-21 change here
                                    if ".pair1." in clean:
                                        os.symlink( clean , os.path.join( sampledirdest , sample + ".clean_R1.fastq.gz" ) )
                                    #Eric Fournier 2020-09-21 change here
                                    elif ".pair2." in clean:
                                        os.symlink( clean , os.path.join( sampledirdest , sample + ".clean_R2.fastq.gz" ) )

                            # Link to metrics
                            dmetrics    = {}
                            if techno == "nanopore" :
                                lmetrics    = glob.glob( os.path.join( sampledirsrc , "*metrics.csv" ) )
                            else :
                                lmetrics    = glob.glob( os.path.join( rundir , "metrics/metrics.csv" ) )
                            if len( lmetrics ) > 0 :
                                os.symlink( lmetrics[0] , os.path.join( sampledirdest , sample + ".metrics.csv" ) )
                                dmetrics    = parseMetrics( lmetrics[0], sname ) 
                            
                            hidden_sample = False
                            if str(sample).startswith('.'):
                                #print("Hidden sample " + sample)
                                hidden_sample = True
                                hiddenSampleFile.write(sample + "\t" + sampledirsrc + "\t" + techno + "\n")
                            else:
                                # Write JSON file
                                dsamplejson = {}
                                dsamplejson[ "key" ]    = samplekey 
                                dsamplejson[ "sample"]  = sample
                                dsamplejson[ "samplemcgill"]  = sname
                                dsamplejson[ "plate"]   = plate
                                dsamplejson[ "platemcgill"] = plateMcGill
                                dsamplejson[ "techno"]  = techno 
                                dsamplejson[ "rundate"] = date
                                dsamplejson[ "runpath"] = rundir
                                dsamplejson[ "sampledir"] = sampledirdest
                                for key, val in dmetrics.iteritems() :
                                    if key == "cons.per.N":
                                        #Eric Fournier 2020-10-12 pour nanopore => cons.perc.N et pour illumina et MGI => cons.per.N
                                        key = "cons.perc.N"
                                    dsamplejson[ key ] = val
                                setQCStatus( dsamplejson )
                                dsamplejson[ "qccuration"] = dsamplejson["qcstatus"]
                                dsamplejson[ "qccurationmess"] = ""

                                with open( os.path.join( sampledirdest , sample + ".json" ) ,  'w') as outfile:
                                    json.dump( dsamplejson, outfile )

                                if dsamplejson["qcstatus"] in ["MISSING_METRICS_HEADER","MISSING_CONS_PERC_N"] and save_consensus_path:
                                    missingConsPercNLog.write(sample + "\t" + sampledirsrc + "\t" + techno + "\t" + dsamplejson["qcstatus"] + "\n")
                                    consensusListFile.write(str(sample).split('_')[0] + "\t" + dsamplejson["qcstatus"] + "\t" + "NA"  + "\t" + techno + "\t" + "NA" + "\t"+ run + "\n")

                                if dsamplejson["qcstatus"] in ["PASS","FLAG","REJ"] and save_consensus_path:
                                    consensusListFile.write(str(sample).split('_')[0] + "\t" + dsamplejson["qcstatus"] + "\t" + lconsensus[0]  + "\t" + techno + "\t" + dsamplejson["cons.perc.N"] + "\t"+ run + "\n")

                                if not save_consensus_path: 
                                    consensusListFile.write(str(sample).split('_')[0] + "\t" + "NA" + "\t" + "NA"  + "\t" + techno + "\t" + "NA" + "\t"+ run + "\n")

def setQCStatus( dsamplejson ) :
    qc      = "PASS"
    lmess   = []
   
    if not "cons.perc.N" in dsamplejson and not "bam.perc.50x" in dsamplejson :
        qc  = "MISSING_METRICS_HEADER"
        lmess.append( "No metrics to do QC" )

    #TODO bug a corriger : car si on entre dans le except le qc demeure a PASS Eric Fournier 2020-10-12
    if "cons.perc.N" in dsamplejson :
        try :
            consperN = float( dsamplejson["cons.perc.N"] )
        except :
            qc = "MISSING_CONS_PERC_N"
            lmess.append( "No cons.per.N to do QC" )
        else :
            if consperN > 5.0 :
                qc      = "REJ"
                lmess.append( "Percentage of N > 5%" )
            elif consperN > 1.0 :
                qc = "FLAG"
                lmess.append( "Percentage of N > 1% and < 5%" )

    if qc == "PASS" and "bam.perc.50x" in dsamplejson :
        try :
            bamperc50x = float( dsamplejson["bam.perc.50x"] )
        except : 
            lmess.append( "No bam.perc.50x to do QC" )
        else:
            if bamperc50x < 90.0 :
                qc = "FLAG"
                lmess.append( "Regions covered at 50X or more < 90%" )

    dsamplejson["qcstatus"] = qc  
    dsamplejson["qcmess"]   = ""
    if len( lmess ) > 0 :
        dsamplejson["qcmess"] = ",".join( lmess )

    return qc        


def getSampleKey( plate , sample , techno , date ) :
    return plate + "." + sample + "." + techno + "." + date 

def parseSampleKey( samplekey ) :
    return samplekey.split( '.' )

def getMcGillPlate( run , techno ) :
    ltmp    = run.split( '_' )
    if techno == "nanopore" and len( ltmp ) == 3 :
        return ltmp[ 0 ]
    if (techno == "illumina" or techno == "mgi" ) and len( ltmp ) == 4 :
        return ltmp[ 2 ]
    return ""
        
def plateForSample( sname, dplate ) :
    lplateforsample = []
    lplate  = sorted( dplate.keys() )
    sample  = sname     # SNAME: McGill name ; SAMPLE : LSPQ name
        
    for plate in lplate :
        if sample in dplate[ plate ] :
            lplateforsample.append(plate)

    # One plate found
    if len( lplateforsample ) ==  1 :
        return lplateforsample[0], sample
    
    # More than one plate
    elif len( lplateforsample ) > 1 :
        print( "**: More than one plate for sample = " + sample )
        return ','.join( lplateforsample ),sample

    # No plate found
    elif len( lplateforsample ) == 0 : 
        lsname  = sname.split('_')
        if len( lsname ) > 1 :
            snametmp        = '_'.join(lsname[:-1])
            plate,sample    = plateForSample( snametmp, dplate )
            if not plate  == "NOTFOUND" :
                #print( "Trimmed sample name " + sname + " => " + sample +" found in plate " + plate)
                return plate,sample
            else :
                pass
                #print( "ERROR: no plate found for sample = " + sname )
                return "NOTFOUND",sname
        else :
            #print( "ERROR: no plate found for sample = " + sname )
            return "NOTFOUND",sname 

def updateRunFile(  ) :
    drun        = {}

    # nanopore
    nanoprocdir = "analysis"
    # Illumina
    illuprocdir = "alignment"
    # mgi
    mgiprocdir  = "alignment"

    # For each technology (higher level folder)
    for tech in TECHNO_LIST :
        if verbose:
            logfile.write("\n    Parsing " + tech  + "\n")
        
        basedir = tech
        techdir = os.path.join( ARGS.dir , basedir )
        # Run name position
        i       = 0
        if tech == "nanopore" :
            i       = -1

        lrun    = os.listdir( techdir )
        # For each run of the  technology
        for run in lrun :
            if verbose:
                logfile.write("\n         run " + str(run) + "\n")

            date    = run.split('_')[i] 
            path    = os.path.join( basedir , run )   # Starting from COVID_full_processing/
            rundir  = os.path.join( techdir , run )   # System path
            lmess   = []

            # Not able to parse run name
            if not date :
                line = ["" , tech , run , path , "ERR", "",  [] , []]
                lmess.append( "ERR: Run name parsing fail" )

            # Run name parsable
            else :

                line    = [date , tech , run,  path, "", "", [] , [] ]
                drun[ path ] = line 
                
                # Find list of samples
                lsamples = [ ]
                if tech == "illumina" or tech == "mgi" :
                    fsamplepath = os.path.join( rundir, run + ".csv" )
                    try :
                        f = open(fsamplepath, "r")
                    except :

                        lmess.append( "ERR: csv file " + run + ".csv" +" with list of samples not found" )

                    else:
                        lines = f.readlines()
                        f.close()
                        for line1 in lines :
                            lcol    = line1.split(',') 
                            if len(lcol) > 0 :
                                lsamples.append(lcol[0])
                # !!!!!! NANOPORE MISSING A WAY OF HAVING TOTAL LIST OF SAMPLES
                elif tech == "nanopore" :
                    lmess.append("ERR: Not possible to get list of samples for nanopore" )

                # Name of the processing folder depends on technology
                if tech == "illumina" :
                    p    = illuprocdir
                elif tech == "nanopore" :
                    p    = nanoprocdir 
                else :
                    p    = mgiprocdir

                procdir = os.path.join( rundir , p )
                line[ 4 ] = "RAW"
                # If we don't have a processing folder => only raw data
                if os.path.isdir( procdir ) and  os.access( procdir , os.X_OK) :

                    try:
                        lproc   = os.listdir ( procdir )
                    #except PermissionError:
                    except :
                        lmess.append("ERROR: permission denied for directory = " + procdir )
                        line[ 4 ] =  "ERR"
                    else :
                        # If processing folder not empty, => data analysed
                        if len( lproc ) > 0 :
                            line[ 4 ] = "PROC"  
            line[ 5 ] = ";".join( lmess )
            # Report messages
            runMessageFile.write("        ###################### Message from " + inspect.stack()[0][3] + " ###################### \n")
            if len( lmess ) > 0 :
                runMessageFile.write(path + " : \t\n")
                for m in lmess :
                    runMessageFile.write("\t" + m + "\n")
            else:
                pass
                runMessageFile.write(path + ": OK\n")
            runMessageFile.write("\n\n")
    return drun

def readPlateFile( fplate ) :
    global LPLATEHEADER
    fplatestr   = open( fplate ) 
    dplate      = {}        # key : plate ; value : list of samples
    lplate      = []

    for line in fplatestr :
        lcol = line.split("\t")

        if len(lcol) < 3 :
            continue

        #Eric Fournier 2020-10-02 modif here, convert to str to avoid warning message
        if str(lcol[0].strip()) == str(LPLATEHEADER[0]) :
            continue

        plate = lcol[2].strip()
        
        if len(plate) > 0 :
            if not plate in lplate :
                lplate.append( plate )
                dplate[ plate ] = []
            dplate[ plate ].append( lcol[0].strip() )

    for plate in lplate :
        if verbose:
            logfile.write("\n         " + plate + " : %i samples"%(len( dplate[ plate ] )) + "\n")
            
    fplatestr.close()
    return dplate

def parseMetrics( fmetrics, samplename ) :
    dmet    = {}
    try :
        fmetstr = open( fmetrics )
    except Exception as e:
        print ("ERROR with " + fmetrics ,  e )
    else :
        lheader = []

        for line in fmetstr :
            lcol = line.split( "," )

            if len(lcol) < 10 :
                continue
            
            if lcol[0] == "sample":
                lheader = lcol        

            if samplename in lcol[0] :
                for i in range( 1, len(lcol) ) :
                    try :
                        val = "%.2f" %float( lcol[i] )
                    except ValueError :
                        val = lcol[i]
                    dmet[ lheader[i] ] = val
                break
    return dmet


def writeRunFile( drun ) :
    fn  = "listRuns_" + datetime.today().strftime('%Y%m%d') + ".tsv"

    #Eric Fournier 2020-10-02 change here
    f   = open( os.path.join(TRACE_DIR,fn) , "w" )

    global LRUNHEADER
    f.write( "\t".join( LRUNHEADER ) + "\n" )

    lpath = sorted( drun.keys() ) 
    for p in lpath :
        nbsampletot     = len( drun[ p ][6] )
        drun[ p ][6]    = "%i"%(nbsampletot)
        nbsamplerep     = len( drun[ p ][7] )
        drun[ p ][7]    = "%i"%(nbsamplerep)
        f.write( "\t".join( drun[ p ] ) + "\n" )
    
    print "End : >>>>>>>   New list of run file in " + os.path.join(TRACE_DIR,fn) + "\n"
    
    if verbose:
        logfile.write("\n   Save runs list to " + os.path.join(TRACE_DIR,fn) + "\n\n********************** End ********************\n")

    f.close()

if __name__ == "__main__":
    try:
        main()
        logfile.close()
        runMessageFile.close()
        missingConsensusLog.close()
        consensusListFile.close()
        missingConsPercNLog.close()
        hiddenSampleFile.close()
    except KeyboardInterrupt:
        print  "Program canceled by user..."
