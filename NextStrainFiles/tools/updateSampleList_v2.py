#!/usr/local/bin/python
# vim:tabstop=4 shiftwidth=4 softtabstop=4 expandtab

'''
Created on 11 mars 2013

Modified by Eric Fournier 2020-10-02

@author: sam
@author: Eric Fournier
'''

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
        """ Docstring 
            .... """

        parser = argparse.ArgumentParser( description=" " )
        parser.add_argument( '-r',  '--repodir',  help="Path to Repository directory",   required=True )
        parser.add_argument( '-p',  '--fplate', help="File with list of plates and samples" )
        parser.add_argument( '-t',  '--tracedir', help="Path to Trace output", required=True  )

        #Eric Fournier 2020-10-02 add here
        parser.add_argument( '-v','--verbose',help="Print output message in logfile",action='store_true')

        global ARGS 
        ARGS = parser.parse_args()

def main() :
    t1    = time()
    print ("\nBegin @ " + str( datetime.now() ))
                                                                
    parseoptions( )
    
    global LPLATEHEADER
    global LSAMPLEHEADER
   
    #Eric Fournier 2020-10-02 add here
    global verbose 
    global logfile
    global missing_json_listfile

    verbose = ARGS.verbose

    #Eric Fournier 2020-10-02 add here
    if not os.path.isdir(ARGS.tracedir):
        os.makedirs(ARGS.tracedir)

    # Load json param
    fjson   = open("runParam_v2.json")
    pjson   = json.load( fjson )
    LPLATEHEADER     = pjson[ 'LPLATEHEADER']
    LSAMPLEHEADER    = pjson[ 'LSAMPLEHEADER' ]
    
    # Current working directory
    cwd = os.getcwd()

    #Eric Fournier 2020-10-02 add here
    logfile_path = os.path.join(ARGS.tracedir,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') +".log")
    logfile = open(logfile_path,'w')
    missing_json_listfile_path =  os.path.join(ARGS.tracedir,os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') +"_missing_json.txt")
    missing_json_listfile = open(missing_json_listfile_path,'w')

    # Read list of plates ; This list has to be provided by LSPQ
    dplate  = {}    # key : plate ; value : list of samples
    if ARGS.fplate :
        if verbose:
            logfile.write("******************** in readPlateFile ************************\n")
        dplate = readPlateFile( ARGS.fplate )

    if verbose:
        logfile.write("\n******************** in updateTrace ************************\n")
    dsampleALL = updateTrace( dplate, ARGS.repodir, ARGS.tracedir )
    
    os.chdir( cwd )

    fjson.close()
    dt = time()-t1
    print (" ---- time(s) = " + str(dt) + " // @SaM")


def updateTrace( dplate, repodir, tracedir ) :
    cwd         = os.getcwd()   
    dsampleALL  = {}

    bckupdir    = os.path.join( tracedir , "BCKUP" )
    ltmp        = glob.glob( bckupdir )

    if len( ltmp ) == 0 :
         os.mkdir( bckupdir )

    lplateLSPQ = sorted( dplate.keys() )

    for plate in lplateLSPQ :
        if verbose:
            logfile.write("\n       in plate : " + plate + "\n")

        # Create plate list of sample file if does not exist
        dsample     = {}
        fsample     = plate + "." + datetime.today().strftime('%Y%m%d') + ".list"
        samplepath  = os.path.join( tracedir , fsample )

        # if exist, read and create bckup file
        if len( glob.glob( samplepath ) ) :
            dsample         = readSampleFile( samplepath )
            os.rename( samplepath , os.path.join( bckupdir , fsample ) )

        # For all samples in LSPQ list of samples, search if exist in repo
        lsamplerepovisit    = []    # List of repo path visited
        platerepodir        = os.path.join( repodir , plate ) 
        for sample in dplate[ plate ] :
            if verbose:
                logfile.write("\n              sample is " + sample + "\n")
            lsamplerepo =  glob.glob( os.path.join( platerepodir , "*" + sample + "*"  ) )
            
            for samplerepo in lsamplerepo :
                lsamplerepovisit.append( samplerepo )
                samplekey = os.path.basename( samplerepo )
                if samplekey in dsample :
                    pass
                    if verbose:
                        logfile.write("\n              !!!!  samplekey " + samplekey + " already in trace file\n")
                
                # New sample not in trace file
                else :
                    addSampleToTrace( plate, samplekey , platerepodir , dsample, "NEW" )
        
        # For all sample in repo (and not LSPQ), add info in trace
        for samplerepo in glob.glob( os.path.join( platerepodir , "*"  ) ) :
            if not samplerepo in lsamplerepovisit :
                samplekey = os.path.basename( samplerepo )
                addSampleToTrace( plate, samplekey , platerepodir , dsample, "Sample not in LSPQ plate" )
        
        # Write sample trace file
        writeSampleFile( samplepath, dsample, False )
        dsampleALL.update( dsample )

    # for plates not in LSPQ file
    lplaterepo  = glob.glob( os.path.join( repodir , "*" ) )
    dsample     = {}
    fsample     = "NOTFOUND" + "." + datetime.today().strftime('%Y%m%d') + ".list"
    samplepath  = os.path.join( tracedir , fsample )
    # if exist, read and create bckup file
    if len( glob.glob( samplepath ) ) :
        os.rename( samplepath , os.path.join( bckupdir , fsample ) )

    for platerepodir in lplaterepo :
        platerepo   = os.path.basename( platerepodir )
        # For plates not in LSPQ list of plates
        if not platerepo in lplateLSPQ :
            # for all samples
            lsamplerepo =  glob.glob( os.path.join( platerepodir , "*"  ) )
            for samplerepo in lsamplerepo :
                samplekey = os.path.basename( samplerepo )
                if not samplekey in dsampleALL :
                    # add line in trace file
                    addSampleToTrace( "NOTFOUND", samplekey , platerepodir , dsample, "Plate not in LSPQ" )
    # Sample list with ALL
    dsampleALL.update( dsample )
    fALL    = "AllSample" +  "." + datetime.today().strftime('%Y%m%d') + ".list"
    writeSampleFile( fALL , dsampleALL, True )

    fstat   = "AllSample" +  "." + datetime.today().strftime('%Y%m%d') + ".stat"
    writeSampleStat( fstat , dsampleALL )

    return dsampleALL

def addSampleToTrace( plate, samplekey , platerepodir , dsample, mess ):

    lsample = [""] * 10

    lsample[0]   = plate
    lsample[1]   = samplekey
    samplerepodir   = os.path.join( platerepodir , samplekey  )
    #Eric Fournier 2020-10-08 add try except
    try:
        fjson           = glob.glob( os.path.join( samplerepodir , "*.json"   ))[0]
        fjsonstr        = open( fjson )
        djson           = json.load( fjsonstr )
    
        lsample[2]   = djson['sample']
        lsample[3]   = djson['techno']
        lsample[4]   = djson['rundate']
        lsample[5]   = djson["qcstatus"]
        if djson["qcstatus"] == "PASS" :
            lsample[6]  = "PASS"
        else :
            lsample[6]  = "UNK"

        lsample[7]  = ""
        lsample[8]  = djson['sampledir']
        lsample[9]  = mess

        dsample[ samplekey ] = lsample
    except:
        missing_json_listfile.write(samplekey + "\n")


def readPlateFile( fplate ) :
    global LPLATEHEADER
    fplatestr   = open( fplate ) 
    dplate      = {}        # key : plate ; value : list of samples

    for line in fplatestr :
        lcol = line.split("\t")

        if len(lcol) < 3 :
            continue

        #Eric Fournier 2020-10-02 convert to str
        if str(lcol[0].strip()) == str(LPLATEHEADER[0]) :
            continue

        plate = lcol[2].strip()
        if len(plate) > 0 :
            if not plate in dplate :
                dplate[ plate ] = []
            dplate[ plate ].append( lcol[0].strip() )

    for plate,lsample in dplate.items() :
        if verbose:
            logfile.write("\n" + plate + " : %i samples"%(len(lsample)))
    fplatestr.close()
    return dplate

# "LSAMPLEHEADER":["Plate[0]", "SampleKEY[1]","Sample[2]" , "Techno[3]" , "RunDate[4]" , "QCStatus{PASS,FLAG,REJ}[5]" , "CurrationStatus{UNK,PASS,REJ}[6]", "AnalysisStatus{SKIP,KEEP}[7]", "REPOPATH[8]" ],

def readSampleFile( fsample ) :
    global LSAMPLEHEADER
    fsamplestr  = open( fsample )
    dsample = {} # key : sample key ; value : list of columns

    #print ("\nParse list of samples")
    for line in fsamplestr : 
        lcol = line.split("\t")

        if len(lcol) < 9 :
            continue

        #Eric Fournier 2020-10-02 convert to str
        if str(lcol[0].strip()) == str(LSAMPLEHEADER[0]) :
            continue

        samplekey = lcol[1].strip()
        if len(samplekey) > 0 :
            if samplekey in dsample :
                print ("ERROR, sample key duplicated : " + samplekey)
            dsample[ samplekey ] = lcol
    
    fsamplestr.close()
    return dsample

def writeSampleFile( fsample, dsample, ismess ) :
    #Eric Fournier 2020-10-03 add here
    out_sample_file =  os.path.join(ARGS.tracedir,fsample)

    #Eric Fournier 2020-10-03 change here
    f   = open( out_sample_file , "w" ) 

    if verbose and ismess:
        logfile.write("\n******************** in writeSampleFile **********************\n\n")
        logfile.write("     write samples to " + out_sample_file + "\n\n")

    global LSAMPLEHEADER
    f.write( "\t".join( LSAMPLEHEADER ) + "\n" )

    lsample = sorted( dsample.keys() )
    for sample in lsample :
        if ismess :
            lsampleval = dsample[ sample ]
        else :
            lsampleval = dsample[ sample ][0:8]
        f.write( "\t".join( lsampleval ) + "\n" )
    f.close()


def writeSampleStat( fstat, dsampleALL ) :
    lsamplekeys = sorted( dsampleALL.keys() )
    dplate      = {}

    #Eric Fournier 2020-10-02 add here
    out_stat_file =  os.path.join(ARGS.tracedir,fstat)

    if verbose:
        logfile.write("\n******************** in writeSampleStat **********************\n\n")
        logfile.write("     write stat to " + out_stat_file)

    # List of plates
    for samplekey in lsamplekeys :
        plate   = dsampleALL[ samplekey ][0]
        if not plate in dplate :
            dplate[ plate ] = []
        dplate[ plate ].append( samplekey )

    #
    lplate = sorted( dplate.keys() )
    lqc         = ["PASS", "FLAG", "REJ", "UNK"]
    dtotsample  = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }
    duniqsample = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }
    dtotdate    = {}
    dtotplate   = {}
    dtottechno  = {}


    for plate in lplate :

        # On ne tiend pas compte des samples dans NOTFOUND
        if plate == "NOTFOUND" :
            continue

        dtotplate[ plate ]  = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] } 

        # Samplelist
        for samplekey in dplate[ plate ] :
            sample  = dsampleALL[ samplekey ][2]
            date    = dsampleALL[ samplekey ][4][:6]
            techno  = dsampleALL[ samplekey ][3]
            qc      = dsampleALL[ samplekey ][5]

            if not qc in lqc :
                qc = "UNK"
            
            dtotsample[qc].append( sample )
            if not date in dtotdate :
                dtotdate[ date ] = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }
            dtotdate[ date ][qc].append(sample)
            if not plate in dtotplate :
                dtotplate[ plate ] = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }
            dtotplate[ plate ][qc].append(sample)
            if not techno in dtottechno :
                dtottechno[ techno ] = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }
            dtottechno[ techno ][qc].append(sample) 

    duniqsample = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }
    duniqdate   = {}
    duniqplate  = {}
    duniqtechno = {}
    
    # All samples

    # PASS uniq
    duniqsample["PASS"] = list(dict.fromkeys(dtotsample["PASS"]))
    
    # FLAG uniq
    lflag   = list(dict.fromkeys(dtotsample["FLAG"]))
    for sample in lflag :
        if not sample in duniqsample["PASS"] :
            duniqsample["FLAG"].append( sample )
           
    # REJ uniq
    lrej    = list(dict.fromkeys(dtotsample["REJ"]))
    for sample in lrej :
        if not sample in duniqsample["FLAG"] and not sample in duniqsample["PASS"] :
            duniqsample["REJ"].append( sample )
           
    # UNK uniq
    lunk    = list(dict.fromkeys(dtotsample["UNK"]))
    for sample in lunk :
        if not sample in duniqsample["REJ"] and not sample in duniqsample["FLAG"] and not sample in duniqsample["PASS"] :
            duniqsample["UNK"].append( sample )
        
    dtotsample["TOTAL"] = len(dtotsample["PASS"]) + len(dtotsample["FLAG"]) + len(dtotsample["REJ"]) + len(dtotsample["UNK"])
    duniqsample["TOTAL"] = len(duniqsample["PASS"]) + len(duniqsample["FLAG"]) + len(duniqsample["REJ"]) + len(duniqsample["UNK"]) 

    # Samples group by techno / date / plate
    dtot    = [ dtottechno , dtotdate, dtotplate ] 
    duniq   = [ duniqtechno, duniqdate, duniqplate ]
    for i in range(0,3) :
        lkeys = sorted(dtot[i].keys()) # key = techno / date / plate
        for key in lkeys :
            duniq[i][key]   = {"PASS":[], "FLAG":[], "REJ":[], "UNK": [] }

            # PASS
            duniq[i][key]["PASS"] = list(dict.fromkeys(dtot[i][key]["PASS"]))
            
            # FLAG
            lflag   = list(dict.fromkeys(dtot[i][key]["FLAG"]))
            for sample in lflag :
                if not sample in duniq[i][key]["PASS"] :
                    duniq[i][key]["FLAG"].append( sample )
            
            # REJ
            lrej    = list(dict.fromkeys(dtot[i][key]["REJ"]))
            for sample in lrej :
                if not sample in duniq[i][key]["FLAG"] and not sample in duniq[i][key]["PASS"] :
                    duniq[i][key]["REJ"].append( sample )
            
            # UNK
            lunk    = list(dict.fromkeys(dtot[i][key]["UNK"]))
            for sample in lunk :
                if not sample in duniq[i][key]["REJ"] and not sample in duniq[i][key]["FLAG"] and not sample in duniq[i][key]["PASS"] :
                    duniq[i][key]["UNK"].append( sample )
       
            dtot[i][key]["TOTAL"] = len(dtot[i][key]["PASS"]) + len(dtot[i][key]["FLAG"]) + len(dtot[i][key]["REJ"]) + len(dtot[i][key]["UNK"])
            duniq[i][key]["TOTAL"] = len(duniq[i][key]["PASS"]) + len(duniq[i][key]["FLAG"]) + len(duniq[i][key]["REJ"]) + len(duniq[i][key]["UNK"]) 
    #Eric Fournier 2020-10-02 change here
    f   = open( out_stat_file , "w" )
    
    # Header
    f.write( "Sample\t" + "\t".join( lqc ) + "\tTotal\n" )

    # All Sample
    linetot     = "Total Sample" 
    lineuniq    = "Uniq Sample"
    for qc in lqc :
        linetot     += "\t%i" % len(dtotsample[qc])
        lineuniq    += "\t%i" % len(duniqsample[qc])
    f.write( linetot    + "\t%i\n" % dtotsample["TOTAL"] )
    f.write( lineuniq   + "\t%i\n" % duniqsample["TOTAL"] )

    # By techno / date / plate
    for i in range(0,3)  :
        lkeys = sorted(dtot[i].keys())  # key = techno / date / plate
        for key in lkeys :
            linetot     = "Total " + key
            lineuniq    = "Uniq " + key
            for qc in lqc :
                linetot     += "\t%i" % len(dtot[i][key][qc])
                lineuniq    += "\t%i" % len(duniq[i][key][qc])
            f.write( linetot    + "\t%i\n" % dtot[i][key]["TOTAL"] )
            f.write( lineuniq   + "\t%i\n" % duniq[i][key]["TOTAL"] )
    f.close()


if __name__ == "__main__":
    try:
        main()
        logfile.close()
    except KeyboardInterrupt:
        print "Program canceled by user..."
