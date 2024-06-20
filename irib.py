from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
import os
import sys
import subprocess
import time

import numpy as np

# This pipeline does: 1: FITLD, INDXR, 2: listr, prtan, plots CL1, makes a VPLOT, 3: POSSMs all srcs, Gnuplots all POSSMs, finds the peak chan for each source
# Run it using:
#    ParselTongue Irib_pipe.py x
# where x is the tmask to start from (defaults to 1 if not sepcified)



############ FITS file

# Get the fits file from the current directory
#fits = subprocess.check_output("ls *IDI* | cut -d'.' -f1 | uniq", shell=True).strip()
#fits = fits.decode()
#print("Fits =", fits)



# Or set the fits file explicitly
#fits = 'irib5.IDI'

#INNAME=''; INCLASS=''; INSEQ=0; INDISK=0; INTY=''; ALLD


def get_tsys(antab_file):
    tsys = dict()
    with open(antab_file, "r") as antab:
         antab_data = antab.readlines()
         sources = []
         for line in antab_data:
             if "source=" in line:
                 src = line.split("source=")[1].replace("\n", "").upper()
                 if src not in sources:
                     sources.append(src)
                     
                 if src not in tsys:
                     
                     tsys[src] = {str(ch):[] for ch in range(0, 16)}
                     
             elif line.startswith("!") or "INDEX=" in line or "TSYS IR FT" in line or "/" in line:
                continue
                
         for line in antab_data:
             if "source=" in line:
                 pass
             elif line.startswith("!") or "INDEX=" in line or "TSYS IR FT" in line or "/" in line:
                continue           
             else:
                 t_sys = line.split(" ")
                 for ch in range(2, 18):
                    for src in sources:
                        if "TSYS IB FT" not in line and "using" not in line:
                            tsys[src][str(ch-2)].append(float(t_sys[ch].replace("\n", "")))
                
    return tsys



exper_name = sys.argv[1]
directory = "/mnt/VLBI/correlations/sfxc/" + exper_name

############# Antab file
DATAC = os.environ['DATAC']

ib_antab_file = DATAC + exper_name.lower() + "ib.antabfs"
ir_antab_file = DATAC + exper_name.lower() + "ir.antabfs"
tsys_ib = get_tsys(ib_antab_file)
tsys_ir = get_tsys(ir_antab_file)

sources_ir = list(tsys_ir.keys())
sources_ib = list(tsys_ib.keys())
channels = list(tsys_ib[sources_ib[0]].keys())

mean_tsys_ib = 0
count_tsys_ib = 0
mean_tsys_ir = 0
count_tsys_ir = 0

for source in sources_ib:
    for ch in channels:
        count_tsys_ib += len(tsys_ib[source][ch])
        mean_tsys_ib += np.sum(tsys_ib[source][ch])

for source in sources_ir:
    for ch in channels:        
        count_tsys_ir += len(tsys_ir[source][ch])
        mean_tsys_ir += np.sum(tsys_ir[source][ch])

mean_tsys_ib /= count_tsys_ib
mean_tsys_ir /= count_tsys_ir


change_ib = False
if 20 > mean_tsys_ib < 100: 
    switchpol_ib = "cd " + directory +  "  &&  sed -e 's/L[1-9]/X/g' -e 's/R\\([1-9]\\)/L\\1|R\\1/g' " + ib_antab_file + " > ib.editedantab"
    os.system(switchpol_ib)
    change_ib = True

change_ir = False       
if 20 > mean_tsys_ir < 100: 
    switchpol_ir = "cd " + directory +  " &&  sed -e 's/L[1-9]/X/g' -e 's/R\\([1-9]\\)/L\\1|R\\1/g' " + ir_antab_file + " > ir.editedantab"
    os.system(switchpol_ir)
    change_ir = True

if change_ib and change_ir:
    cat = "cd " + directory + " &&  cat *editedantab > ./antab"
    os.system(cat)
    
elif change_ib and not change_ir:
    cat = "cd " + directory + " &&  cat  ib.editedantab " + ir_antab_file + " > ./antab"
    os.system(cat)

elif not change_ib and change_ir:
    cat = "cd " + directory + " &&  cat  ir.editedantab " + ib_antab_file + " > ./antab"
    os.system(cat)
 
else:
    cat = "cd " + directory + " &&  cat " + ir_antab_file + " " + ib_antab_file + " > ./antab"
    os.system(cat)   

############# Flag file
cat = "cd " + directory + " &&  cat *uvflgfs > ./uvflg"
os.system(cat)

refant = 1

#################### This is where you set the restfrequency of the maser
tfreq = 6668519200
freq_h = float(int(tfreq*1e-6))*1e6
freq_l = tfreq - freq_h

# Get the Tmask either from command or from this script
if len(sys.argv) > 2:
    tmask =int(sys.argv[2])
else:
    tmask = 1

# The AIPS ID
AIPS.userno = int(time.time()/100000)
print("  Doing from Tmask:", tmask)
print("  AIPS ID:", AIPS.userno)
print("  Refant:", refant)
print("  FreqRef:", tfreq, freq_h, freq_l)

# Open a file called "essentials" where some useful info gets printed
essentials = open(directory + "/" + "essentials.txt", "w")

##### Print some info at the start
print('DATAA:' + exper_name + '_contiuum.IDI')
print('DATAB:' + exper_name + '_line.IDI')

##### Define the different data types we will use
contdata = AIPSUVData('CONT' + exper_name, 'UVDATA',1,1)
linedata = AIPSUVData('LINE' + exper_name, 'UVDATA',1,1)


if tmask <= 1:
    if contdata.exists():
        print("Zapping old data")
        contdata.zap()
    if linedata.exists():
        print("Zapping old data")
        linedata.zap()
    fitld = AIPSTask('FITLD')
    fitld.clint = 0.166
    fitld.doconcat = 1
    fitld.ncount = 2
    datain = 'DATAA:' + exper_name + '_contiuum.IDI'
    fitld.datain = datain
    fitld.outdata = contdata
    print("Fits is:", datain)
    fitld.go()
    ##
    datain = 'DATAB:' + exper_name + '_line.IDI'
    fitld.datain = datain
    fitld.outdata = linedata
    print("Fits is:", datain)
    fitld.go()
    #########Index data
    indxr = AIPSTask('indxr')
    indxr.cparm[3] = 0.166
    indxr.indata = contdata
    indxr.go()
    ########
    indxr.indata = linedata
    indxr.go()
    ########

##### Define list of targets in the full project source list (Currently just three but there will be 30 or so when the source list is finalised by Artis. Call this 'all_masers'
##### Then define a list called 'masers' which are the sources that are in both the 'all_masers' list and also in the uvdata in this particlar fits file.
##### Then define a list of all the sources that are in the uvdata but are not masers. These will be the calibrators, call this 'conts'.
all_masers = ['W3OH',
              'G111',
              'G196.454-1.6',
              'G188.946+0.8',
              'S255',
              'G009.621+0.1',
              'M_17',
              'M_17_SWex',
              'G022.356+0.0',
              'G023.706-0.1',
              'G24.33+0.14',
              'G028.304-0.3',
              'G029.955-0.0',
              'G030.224-0.1',
              'G030.378-0.2',
              'G030.817-0.0',
              'W43',
              'G033.641-0.2',
              'G034.244+0.1',
              'G32.744-0.07',
              'G043.795-0.1',
              'G045.472+0.1',
              'G049.599-0.2',
              'G075.782+0.3',
              'G069.539-0.9',
              'G78.122',
              'DR21(OH)N',
              'W75N',
              'G085.410+0.0',
              'CEPA',
              'L1206',
              'G107.298+5.6']
masers = set(linedata.sources) & set(all_masers)
conts = list(set(linedata.sources).difference(masers))

# Or if you want to explicitly set the maser and cont list:
# masers = ['G85.411', 'CEPA']
# conts = ['3C345', '3C454.3']

essentials.write("Maser sources are " + str(masers) + "\n")
essentials.write("Calibrator sources are " + str(conts) + "\n")

############################## Listr and Prtan
if tmask <= 2:
    essentials.write("Reached tmask 2 \n")
    if os.path.exists('cont_listr.txt'):
        print("  Removing old listr.txt")
        os.remove('cont_listr.txt')
    if os.path.exists('line_listr.txt'):
        print("  Removing old listr.txt")
        os.remove('line_listr.txt')

    listr = AIPSTask('LISTR')
    listr.optype = 'SCAN'
    listr.indata = contdata
    listr.outprint = 'DATAC:cont_listr.txt'
    listr.go()
    listr.indata = linedata
    listr.outprint = 'DATAC:line_listr.txt'
    listr.go()
    #######
    if os.path.exists('prtan.txt'):
        print("  Removing old prtan.txt")
        os.remove('prtan.txt')

    prtan = AIPSTask('PRTAN')
    prtan.indata = contdata
    prtan.outprint = 'DATAC:prtan.txt'
    prtan.go()
    #######

    while contdata.table_highver('AIPS PL') > 0:
        print("  Removing old PL tables")
        contdata.zap_table('AIPS PL', -1)

    if os.path.exists('VPLOT.eps'):
        print("  Removing old VPLOT.eps")
        os.remove('VPLOT.eps')

    #######
    snplt = AIPSTask('SNPLT')
    snplt.indata = contdata
    snplt.inext = 'CL'
    snplt.inver = 1
    snplt.dotv = -1
    snplt.nplots = 8
    snplt.optype = 'AMP'
    snplt.go()
    #######
    if os.path.exists('CL1.eps'):
        print("removing old CL1.eps")
        os.remove('CL1.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:CL1.eps'
    lwpla.go()
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

    ###### Plotting CL1 VPLOT
    vplot = AIPSTask('vplot')
    vplot.indata = contdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    # vplot.avgchan = 1
    # vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    # vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (contdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######

    if os.path.exists('VPLOT_CL1.eps'):
        print("removing old VPLOT_CL1.eps")
        os.remove('VPLOT_CL1.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:VPLOT_CL1.eps'
    lwpla.go()
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

################ Do ACCOR, ANTAB, APCAL
if tmask <= 3:
    essentials.write("Reached tmask 3 \n")
    # Start by exdesting any old SN tables
    while contdata.table_highver('AIPS SN') > 0:
        contdata.zap_table('AIPS SN', 0)
    #######
    accor = AIPSTask('accor')
    accor.indata = contdata
    accor.solint = 1
    accor.go()
    #######
    clcal = AIPSTask('CLCAL')
    clcal.indata = contdata
    clcal.snver = 1
    clcal.gainver = 1
    clcal.gainuse = 2
    clcal.interpol = 'self'
    clcal.go()
    ####### ANTAB
    # Start by exdesting any old GC tables
    while contdata.table_highver('AIPS TY') > 0:
        contdata.zap_table('AIPS TY', 0)
    while contdata.table_highver('AIPS GC') > 0:
        contdata.zap_table('AIPS GC', 0)
    antab = AIPSTask('ANTAB')
    antab.indata = contdata
    antab.calin = 'DATAC:antab'
    antab.go()
    #######
    apcal = AIPSTask('APCAL')
    apcal.indata = contdata
    apcal.dofit[1] = -1
    print(apcal.dofit)
    apcal.go()
    #######
    clcal = AIPSTask('CLCAL')
    clcal.indata = contdata
    clcal.snver = 2
    clcal.gainver = 2
    clcal.gainuse = 3
    clcal.interpol = 'self'
    clcal.go()
    #######
    #######  Plot the CL table
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    snplt = AIPSTask('SNPLT')
    snplt.indata = contdata
    snplt.inext = 'CL'
    snplt.inver = 3
    snplt.dotv = -1
    snplt.nplots = 8
    snplt.optype = 'AMP'
    snplt.go()
    #######
    if os.path.exists('CL3.eps'):
        print("removing old CL3.eps")
        os.remove('CL3.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:CL3.eps'
    lwpla.go()

    #######
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    vplot = AIPSTask('vplot')
    vplot.indata = contdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    # vplot.avgchan = 1
    # vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    # vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (contdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.docal = 1
    vplot.gainuse = 3
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######
    if os.path.exists('VPLOT_CL3.eps'):
        print("removing old VPLOT_CL3.eps")
        os.remove('VPLOT_CL3.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = 2
    lwpla.outfile = 'DATAD:VPLOT_CL3.eps'
    lwpla.go()

################ Do SETJY BPASS and POSSSM
if tmask <= 4:
    essentials.write("Reached tmask 4 \n")
    ############ Flag using the flag file provided
    uvflg = AIPSTask('uvflg')
    uvflg.indata = contdata
    uvflg.intext = 'DATAC:uvflg'
    uvflg.intext = 'DATAC:uvflg'
    uvflg.opcode = 'flag'
    uvflg.go()
    ######### Flag the band edges (10%)
    ########## Have the Pipeline determine the number of channels, and calculate the 10% band edge
    CONT_edgeCh = int((contdata.header['naxis'][2]) / 10)
    print (' Number of chans, and the 10% of it', (contdata.header['naxis'][2]), CONT_edgeCh)
    ###
    LINE_edgeCh = int((linedata.header['naxis'][2]) / 10)
    print (' Number of chans, and the 10% of it', (linedata.header['naxis'][2]), LINE_edgeCh)
    ###
    uvflg = AIPSTask('uvflg')
    uvflg.opcode = 'flag'
    #
    uvflg.indata = contdata
    uvflg.bchan = 1
    uvflg.echan = CONT_edgeCh
    uvflg.go()
    uvflg.bchan = (contdata.header['naxis'][2]) - CONT_edgeCh
    uvflg.echan = (contdata.header['naxis'][2])
    uvflg.go()
    #
    uvflg.indata = linedata
    uvflg.bchan = 1
    uvflg.echan = LINE_edgeCh
    uvflg.go()
    uvflg.bchan = (linedata.header['naxis'][2]) - LINE_edgeCh
    uvflg.echan = (linedata.header['naxis'][2])
    uvflg.go()
    #########
    setjy = AIPSTask('setjy')
    setjy.indata = linedata
    # setjy.sources[1:] = ''
    setjy.restfreq[1:] = [freq_h, freq_l]
    setjy.veltyp = 'LSR'
    setjy.optype = 'VCAL'
    setjy.veldef = 'RADIO'
    setjy.go()
    setjy.indata = contdata
    setjy.go()
    ########## POSSM all IFs and Stokes for all srcs, no cal, no flag, no BP
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.flagver = -1
        possm.nplots = 2
        possm.dotv = -1
        # possm.doband = 1
        possm.docal = -1
        possm.aparm[7] = 1
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_CROSS_CL1.eps'):
            print("removing old POSSM_CROSS_CL1.eps")
            os.remove('POSSM_CROSS_CL1.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'DATAD:POSSM_CROSS_CL1.eps'
        lwpla.go()
    ########## POSSM all IFs and Stokes for all srcs, no cal, no flag, with BP
    while contdata.table_highver('AIPS BP') > 0:
        contdata.zap_table('AIPS BP', 0)
    bpass = AIPSTask('BPASS')
    bpass.indata = contdata
    bpass.bpassprm[1] = 1
    bpass.docal = 1
    bpass.gainuse = 3
    bpass.bpassprm[1] = 1
    # bpass.bpassprm[10] = 1
    # bpass.bpassprm[5] = 1
    bpass.calsour[1:] = conts
    bpass.go()
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.flagver = -1
        possm.nplots = 4
        possm.dotv = -1
        possm.doband = -1
        possm.docal = -1
        possm.aparm[7] = 1
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        possm.aparm[8] = 1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_AUTO_BP.eps'):
            print("removing old POSSM_AUTO_BP.eps")
            os.remove('POSSM_AUTO_BP.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'DATAD:POSSM_AUTO_BP.eps'
        lwpla.go()
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.flagver = 1
        possm.nplots = 4
        possm.dotv = -1
        possm.doband = -1
        possm.docal = 1
        possm.gainuse = 3
        possm.aparm[7] = 0
        possm.aparm[9] = 1
        possm.aparm[8] = 1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_AUTO_CL3.eps'):
            print("removing old POSSM_AUTO_CL3.eps")
            os.remove('POSSM_AUTO_CL3.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'DATAD:POSSM_AUTO_CL3.eps'
        lwpla.go()
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources[1:] = [source]
        possm.stokes = 'HALF'
        possm.aparm[1] = -1
        possm.flagver = 1
        possm.nplots = 2
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 3
        possm.aparm[7] = 1
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_CROSS_CL3_BP.eps'):
            print("removing old POSSM_CROSS_CL3_BP.eps")
            os.remove('POSSM_CROSS_CL3_BP.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = 12
        lwpla.outfile = 'DATAD:POSSM_CROSS_CL3_BP.eps'
        lwpla.go()

################ Do FRING (ALIGN IFs SO THEY CAN BE AVGD) zero Rates
if tmask <= 5:
    essentials.write("Reached tmask 5 \n")
    print("Tmaks 5")
    while contdata.table_highver('AIPS SN') > 2:
        contdata.zap_table('AIPS SN', 0)
    fring = AIPSTask('FRING')
    fring.indata = contdata
    fring.calsour[1:] = conts
    fring.aparm = AIPSList([2, -1, 0, 0, 0, 0, 3, 0, 0])
    fring.dparm = AIPSList([1, 200, 300, 1, 0, 0, 0, 1])
    fring.doband = -1
    fring.refant = 1
    fring.docal = 1
    fring.gainuse = 3
    fring.solint = 10
    # Makes SN3
    fring.go()
    #######
    while contdata.table_highver('AIPS CL') > 3:
        contdata.zap_table('AIPS CL', 0)
    clcal = AIPSTask('CLCAL')
    clcal.indata = contdata
    clcal.snver = 3
    clcal.gainver = 3
    clcal.gainuse = 4
    # clcal.sources = AIPSList(['3C345'])
    clcal.interpol = 'SIMP'
    clcal.bparm = AIPSList([0, 0, 0, 1])
    clcal.dobtween = 1
    clcal.go()
    #######  Plot the CL table
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    snplt = AIPSTask('SNPLT')
    snplt.indata = contdata
    snplt.inext = 'CL'
    snplt.inver = 4
    snplt.dotv = -1
    snplt.nplots = 8
    snplt.optype = 'DELA'
    snplt.go()
    #######
    if os.path.exists('CL4.eps'):
        print("removing old CL4.eps")
        os.remove('CL4.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:CL4.eps'
    lwpla.go()

    ###### Plotting CL4 VPLOT
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

    vplot = AIPSTask('vplot')
    vplot.indata = contdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    # vplot.avgchan = 1
    # vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    # vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (contdata.header['naxis'][2])
    vplot.docal = 1
    vplot.gainuse = 4
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######

    if os.path.exists('VPLOT_CL4.eps'):
        print("removing old VPLOT_CL4.eps")
        os.remove('VPLOT_CL4.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:VPLOT_CL4.eps'
    lwpla.go()
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

################ Do FRING (GET SINGLE SOL FOR ALL IFS)
if tmask <= 6:
    essentials.write("Reached tmask 6 \n")
    print("Tmaks 6")
    while contdata.table_highver('AIPS SN') > 3:
        contdata.zap_table('AIPS SN', 0)
    fring = AIPSTask('FRING')
    fring.indata = contdata
    fring.calsour[1:] = conts
    fring.aparm = AIPSList([2, -1, 0, 0, 1, 0, 3, 0, 0])
    fring.dparm = AIPSList([1, 200, 300, 1, 0, 0, 0, 0])
    fring.doband = -1
    fring.refant = 1
    fring.docal = 1
    fring.gainuse = 4
    fring.solint = 1
    # fring.solsub = 4
    fring.go()
    #######
    while contdata.table_highver('AIPS CL') > 4:
        contdata.zap_table('AIPS CL', 0)
    clcal = AIPSTask('CLCAL')
    clcal.indata = contdata
    clcal.snver = 4
    clcal.gainver = 4
    clcal.gainuse = 5
    # clcal.sources = AIPSList(['3C345'])
    clcal.interpol = 'AMBG'
    clcal.bparm = AIPSList([0, 0, 0, 1])
    clcal.dobtween = 1
    clcal.go()

    ### Plot CL5, kill PL files first
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

    if os.path.exists('CL5.eps'):
        print("removing old CL5.eps")
        os.remove('CL5.eps')

    #######
    snplt = AIPSTask('SNPLT')
    snplt.indata = contdata
    snplt.inext = 'CL'
    snplt.inver = 5
    snplt.do3col = 1
    snplt.opcode = 'ALIF'
    snplt.dotv = -1
    snplt.nplots = 4
    snplt.optype = 'PHAS'
    snplt.pixrange = AIPSList([-180, 180])
    snplt.go()
    #######

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:CL5.eps'
    lwpla.go()
    #### Finish plotting CL7

    ###### Plotting CL5 VPLOT
    vplot = AIPSTask('vplot')
    vplot.indata = contdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    # vplot.avgchan = 1
    # vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    # vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (contdata.header['naxis'][2])
    vplot.docal = 1
    vplot.gainuse = 5
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######

    if os.path.exists('VPLOT_CL5.eps'):
        print("removing old VPLOT_CL5.eps")
        os.remove('VPLOT_CL5.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:VPLOT_CL5.eps'
    lwpla.go()

    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    #### Plot CL5 POSSM
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        print('doing source' + source)
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.stokes = 'HALF'
        possm.aparm[1] = -1
        possm.flagver = 1
        possm.nplots = 2
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 5
        possm.aparm[7] = 1
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.factor = 0.1
        possm.go()

    if os.path.exists('POSSM_CROSS_CL5.eps'):
        print("removing old POSSM_CROSS_CL5.eps")
        os.remove('POSSM_CROSS_CL5.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:POSSM_CROSS_CL5.eps'
    lwpla.go()

###################  FRING the maser for rate sols
# First we need to run POSSM on the maser in the calibrated CONT data set to get the masIF
# Then copy calibration to the LINE data set and run FRING there, and POSSM again for nice zoomed spectra

if tmask <= 7:
    essentials.write("Reached tmask 7 \n")
    print("Tmask 7")
    while contdata.table_highver('AIPS CL') > 5:
        contdata.zap_table('AIPS CL', 0)

    ########## Superarray
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in masers:
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.nplots = 0
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 5
        possm.aparm[7] = 1
        possm.aparm[1] = -1
        possm.outtext = './' + source + '.TXT'
        sa_spec = './' + source + '.TXT'
        if os.path.exists(sa_spec):
            print("removing old plot", sa_spec)
            os.remove(sa_spec)
        possm.go()
        gnuplot = "gnuplot -e \"set term postscript eps; set out '" + source + "_SA.eps'; plot '" + sa_spec + "' u 5:6 w lines\""
        os.system(gnuplot)
        if os.path.exists('Spectra'):
            mv = "mv " + sa_spec + " Spectra"
            os.system(mv)
            mv = "mv " + source + "_SA.eps Spectra"
            os.system(mv)
        else:
            os.system('mkdir Spectra')
            mv = "mv " + sa_spec + " Spectra"
            os.system(mv)
            mv = "mv " + source + "_SA.eps Spectra"
            os.system(mv)

        #  for source in ['IRAS18056']:
        print("   we look for maschan")
        cmd = "sed -E -n '/LL|RR/p' Spectra/" + source + ".TXT | sort -k6 -n | tail -n 1"
        maspeak = subprocess.check_output(cmd, shell=True)
        print(cmd)
        print("Maspeak =")
        print(maspeak)
        mas = maspeak.split()
        maschan = int(mas[0])
        print("   Source and maser ch:", source, " ", maschan)
        # Write in the maser peak channels into Essentials.txt
        essentials.write("(CONT) Maser peak channels are \n")
        essentials.write(source + " " + str(maschan))
        essentials.write('\n')

        ####### Zoomed possm in the CONT data set
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.nplots = 0
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 5
        possm.aparm[7] = 1
        possm.aparm[1] = -1
        possm.stokes = 'I'
        if maschan > 50:
            possm.bchan = maschan - 50
        else:
            possm.bchan = 1
        if maschan + 50 > contdata.header['naxis'][2]:
            possm.echan = contdata.header['naxis'][2]
        else:
            possm.echan = maschan + 50
        possm.bif = int(mas[1])
        possm.eif = int(mas[1])
        possm.outtext = './' + source + '.zoom.TXT'
        sa_spec = './' + source + '.zoom.TXT'
        if os.path.exists(sa_spec):
            print("removing old plot", sa_spec)
            os.remove(sa_spec)
        possm.go()

        # Gnuplot the TXT spectra
        gnuplot = "gnuplot -e \"set term postscript eps; set out '" + source + "_SA_zoom.eps'; set xlabel 'Velocity [km/s]'; set ylabel 'Flux density [Jy]'; set yrange [:]; plot '" + sa_spec + "' u 5:6 w lines title 'IF " + str(
            possm.bif) + "' \""
        os.system(gnuplot)
        if os.path.exists('Spectra'):
            mv = "mv " + source + "_SA_zoom.eps Spectra"
            os.system(mv)
            mv = "mv " + source + ".zoom.TXT Spectra"
            os.system(mv)
        else:
            os.system('mkdir Spectra')
            mv = "mv " + source + "_SA_zoom.eps Spectra"
            os.system(mv)
            mv = "mv " + source + ".zoom.TXT Spectra"
            os.system(mv)

        ### Copy CONT calibration in masIF to allIFs then copy to linedata
        clcop = AIPSTask('CLCOP')
        clcop.indata = contdata
        clcop.inext = 'CL'
        clcop.invers = 5
        clcop.optype = 'AVER'
        clcop.fparm = AIPSList([int(mas[1])])
        clcop.vparm = AIPSList([1, 2, 4, 5, 6, 7, 8])
        clcop.go()
        ### Makes CL6
        while linedata.table_highver('AIPS CL') > 1:
            linedata.zap_table('AIPS CL', 0)
        tacop = AIPSTask('TACOP')
        tacop.indata = contdata
        tacop.outdata = linedata
        tacop.inext = 'CL'
        tacop.ncount = 1
        tacop.invers = 6
        tacop.go()
        ### Makes SN6

        ########## Superarray
        while linedata.table_highver('AIPS PL') > 0:
            linedata.zap_table('AIPS PL', 0)
        possm = AIPSTask('possm')
        possm.indata = linedata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.nplots = 0
        possm.dotv = -1
        possm.doband = -1
        possm.docal = 1
        possm.gainuse = 2
        possm.aparm[7] = 2
        possm.aparm[1] = -1
        possm.outtext = './' + source + '_LINE.TXT'
        sa_spec = './' + source + '_LINE.TXT'
        if os.path.exists(sa_spec):
            print("removing old plot", sa_spec)
            os.remove(sa_spec)
        possm.go()
        gnuplot = "gnuplot -e \"set term postscript eps; set out '" + source + "_SA_LINE.eps'; plot '" + sa_spec + "' u 5:6 w lines\""
        os.system(gnuplot)
        if os.path.exists('Spectra'):
            mv = "mv " + sa_spec + " Spectra"
            os.system(mv)
            mv = "mv " + source + "_SA_LINE.eps Spectra"
            os.system(mv)
        else:
            os.system('mkdir Spectra')
            mv = "mv " + sa_spec + " Spectra"
            os.system(mv)
            mv = "mv " + source + "_SA_LINE.eps Spectra"
            os.system(mv)

        #  for source in ['IRAS18056']:
        print("   we look for maschan")
        cmd = "sed -E -n '/LL|RR/p' Spectra/" + source + "_LINE.TXT | sort -k6 -n | tail -n 1"
        maspeak = subprocess.check_output(cmd, shell=True)
        print(cmd)
        print("Maspeak =")
        print(maspeak)
        mas = maspeak.split()
        maschan = int(mas[0])
        masflux = mas[5]
        masflux = masflux.decode()
        print("   Source and maser ch:", source, " ", maschan)
        print("   Max Flux:", source, " ", masflux)
        # Write in the maser peak channels into Essentials.txt
        essentials.write("(LINE) Maser peak channels are \n")
        essentials.write(source + " " + str(maschan))
        essentials.write('\n')
        essentials.write("Maser peak flux is \n")
        essentials.write(source + " " + masflux)
        essentials.write('\n')

        ####### Zoomed possm in the LINE data set
        possm = AIPSTask('possm')
        possm.indata = linedata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.nplots = 0
        possm.dotv = -1
        possm.doband = -1
        possm.docal = 1
        possm.gainuse = 2
        possm.aparm[7] = 2
        possm.aparm[1] = -1
        possm.stokes = 'I'
        if maschan > 100:
            possm.bchan = maschan - 100
        else:
            possm.bchan = 1
        if maschan + 100 > linedata.header['naxis'][2]:
            possm.echan = linedata.header['naxis'][2]
        else:
            possm.echan = maschan + 100
        possm.outtext = './' + source + '.LINE_zoom.TXT'
        sa_spec = './' + source + '.LINE_zoom.TXT'
        if os.path.exists(sa_spec):
            print("removing old plot", sa_spec)
            os.remove(sa_spec)
        possm.go()

        # Gnuplot the TXT spectra
        sed = "sed -n '13,33p' < " + source + ".LINE_zoom.TXT  | awk '{cnt+=$6; num+=1} END{print cnt/num}'"
        first10 = subprocess.check_output(sed, shell=True)
        first10 = first10.decode()
        print("First ten", first10)
        essentials.write("Maser first 10 channel flux is \n")
        essentials.write(source + " " + str(first10))

        essentials.write("Maser flux - baseline is \n")
        realflux = float(masflux) - float(first10)
        essentials.write(source + " " + str(realflux))
        essentials.write('\n')

        essentials.write("Maser SNR is \n")
        snr = float(masflux) / float(first10)
        essentials.write(source + " " + str(snr))
        essentials.write('\n')
        essentials.write('\n')

        gnuplot = "gnuplot -e \"set term postscript eps; set out '" + source + "_LINE_zoom.eps'; set xlabel 'Velocity [km/s]'; set ylabel 'Flux density [Jy]'; set yrange [-1:]; plot '" + sa_spec + "' u 5:(\$6-" + first10 + ") every ::14 w lines title 'Source:" + source + "' \""
        os.system(gnuplot)

        if os.path.exists('Spectra'):
            mv = "mv " + source + "_LINE_zoom.eps Spectra"
            os.system(mv)
            mv = "mv " + source + ".LINE_zoom.TXT Spectra"
            os.system(mv)
        else:
            os.system('mkdir Spectra')
            mv = "mv " + source + "_LINE_zoom.eps Spectra"
            os.system(mv)
            mv = "mv " + source + ".LINE_zoom.TXT Spectra"
            os.system(mv)

        ##### Dirty eorkaround whic will hopefully be elegantified later in PT
        splatdata = AIPSUVData(source, 'SPLAT', 1, 1)
        if splatdata.exists():
            print("   Zapping old SPLAT data ")
            splatdata.zap()
        splat = AIPSTask('SPLAT')
        splat.indata = linedata
        splat.outdata = splatdata
        splat.source[0:] = AIPSList([source])
        splat.aparm[1] = -1
        splat.doband = -1
        splat.docalib = 1
        splat.gainuse = 2
        splat.outdisk = 1
        splat.go()

        ############# Maser FRING
        # Need to figure out how to deal with SN numbers here. Maybe make, SNCOP, use, delete
        while splatdata.table_highver('AIPS SN') > 0:
            splatdata.zap_table('AIPS SN', 0)
        fring = AIPSTask('FRING')
        fring.indata = splatdata
        fring.calsour = AIPSList([source])
        fring.aparm = AIPSList([2, 0, 3, 0, 0, 0, 1, 0, 0])
        fring.dparm = AIPSList([0, -1, 200, 0, 0, 0, 0, 0])
        fring.doband = -1
        fring.refant = 1
        # fring.docal = 1
        # fring.gainuse = 6
        fring.docal = -1
        # fring.gainuse = 6
        # fring.bif = 3
        # fring.eif = 3
        fring.bchan = maschan
        fring.echan = maschan
        fring.solint = 0.333
        fring.go()
        ####### makes SN1

        ### Copy SN5 masIF sols to all IFs
        while linedata.table_highver('AIPS SN') > 0:
            linedata.zap_table('AIPS SN', 0)
        while contdata.table_highver('AIPS SN') > 4:
            contdata.zap_table('AIPS SN', 0)
        # tacop = AIPSTask('TACOP')
        # tacop.indata = splatdata
        # tacop.inext = 'SN'
        # tacop.invers = 1
        # tacop.ncount = 1
        # tacop.outdata = contdata
        # tacop.go()
        # tacop.outdata = linedata
        # tacop.go()

        ## NO! need to copy CONT SN5 single IF to all 8 IFs
        sncop = AIPSTask('SNCOP')
        sncop.indata = splatdata
        sncop.outdata = contdata
        sncop.invers = 1
        sncop.go()

        # Done with SPLAT data now so zap it
        if splatdata.exists():
            print("   Zapping old SPLAT data ")
        #   splatdata.zap()
        ### Need to copy LCP to RCP here (temporarily)
        # clcop = AIPSTask('CLCOP')
        # clcop.indata = contdata
        # clcop.inext = 'SN'
        # clcop.invers = 5
        # clcop.optype = 'R2L'
        # clcop.fparm = AIPSList([ 7])
        # clcop.vparm = AIPSList([ 1,2,4,5,6,7,8])
        # clcop.go()

        ### Need to copy IF not 1 or 3 phases into the others because 1 and 3 are bugs
        # clcop = AIPSTask('CLCOP')
        # clcop.indata = contdata
        # clcop.inext = 'SN'
        # clcop.invers = 6
        # clcop.optype = 'AVER'
        # clcop.fparm = AIPSList([ 7])
        # clcop.vparm = AIPSList([ 1,2,3,4,5,6,7,8])
        # clcop.go()

        # Makes SN5 in the CONT data, apply it since this is being done source by source
        clcal = AIPSTask('CLCAL')
        clcal.indata = contdata
        clcal.snver = 5
        clcal.gainver = 5
        clcal.gainuse = 7
        clcal.sources = AIPSList([source])
        clcal.calsour = AIPSList([source])
        clcal.interpol = 'AMBG'
        clcal.opcode = 'CALP'
        clcal.go()

        ### Now we're done with the SPLIT data just kill it
        # if splatdata.exists():
        #   print("   Zapping old SPLAT data ")
        #   splatdata.zap()

        ### Plot SN5, kill PL files first
        while contdata.table_highver('AIPS PL') > 0:
            contdata.zap_table('AIPS PL', 0)

        if os.path.exists('SN5_' + source + '.eps'):
            print("removing old SN5.eps")
            os.remove('SN5_' + source + '.eps')

        #######
        snplt = AIPSTask('SNPLT')
        snplt.indata = contdata
        snplt.inext = 'SN'
        snplt.inver = 5
        snplt.dotv = -1
        snplt.nplots = 4
        snplt.do3col = 1
        snplt.opcode = 'ALIF'
        snplt.optype = 'PHAS'
        snplt.pixrange = AIPSList([-180, 180])
        snplt.go()
        #######

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = contdata.table_highver('AIPS PL')
        lwpla.outfile = 'DATAD:SN5_' + source + '.eps'
        lwpla.go()
        #### Finish plotting SN5 for each source

    ### Plot CL7, kill PL files first
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

    if os.path.exists('CL7.eps'):
        print("removing old CL7.eps")
        os.remove('CL7.eps')

    #######
    snplt = AIPSTask('SNPLT')
    snplt.indata = contdata
    snplt.inext = 'CL'
    snplt.inver = 7
    snplt.do3col = 1
    snplt.opcode = 'ALIF'
    snplt.dotv = -1
    snplt.nplots = 4
    snplt.optype = 'PHAS'
    snplt.pixrange = AIPSList([-180, 180])
    snplt.go()
    #######

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:CL7.eps'
    lwpla.go()
    #### Finish plotting CL7

    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    vplot = AIPSTask('vplot')
    vplot.indata = contdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    # vplot.avgchan = 1
    # vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    # vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (contdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.docal = 1
    vplot.gainuse = 7
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######
    if os.path.exists('VPLOT_CL7.eps'):
        print("removing old VPLOT_CL7.eps")
        os.remove('VPLOT_CL7.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:VPLOT_CL7.eps'
    lwpla.go()

    #### Plot CL6 POSSM
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        print('doing source ' + source)
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.stokes = 'HALF'
        possm.aparm[1] = -1
        possm.flagver = 1
        possm.nplots = 2
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 7
        possm.aparm[7] = 2
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.factor = 0.1
        possm.go()

        if os.path.exists('POSSM_CROSS_CL7.eps'):
            print("removing old POSSM_CROSS_CL7.eps")
            os.remove('POSSM_CROSS_CL7.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = contdata.table_highver('AIPS PL')
        lwpla.outfile = 'DATAD:POSSM_CROSS_CL7.eps'
        lwpla.go()

if tmask <= 8:
    ########## FRING all sources: long solint, IFs indivi
    essentials.write("Reached tmask 8 \n")
    while contdata.table_highver('AIPS SN') > 6:
        contdata.zap_table('AIPS SN', 0)
    fring = AIPSTask('FRING')
    fring.indata = contdata
    # fring.calsour = AIPSList([ contdata.sources])
    fring.aparm = AIPSList([2, -1, 0, 0, 0, 0, 3, 0, 0])
    fring.dparm = AIPSList([1, 200, 300, 1, 0, 0, 0, 1])
    fring.doband = -1
    fring.refant = 1
    fring.docal = 1
    fring.gainuse = 7
    # fring.bif = int(mas[1])
    # fring.eif = int(mas[1])
    # fring.bchan = maschan
    # fring.echan = maschan
    fring.solint = 10
    fring.go()
    ####### makes SN6

    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)

    #######
    snplt = AIPSTask('SNPLT')
    snplt.indata = contdata
    snplt.inext = 'SN'
    snplt.inver = 6
    snplt.dotv = -1
    snplt.nplots = 8
    snplt.optype = 'DELA'
    snplt.go()
    #######

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:SN6.eps'
    lwpla.go()

if tmask <= 9:
    while contdata.table_highver('AIPS CL') > 7:
        contdata.zap_table('AIPS CL', 0)

    essentials.write("Reached tmask 9 \n")
    clcal = AIPSTask('CLCAL')
    clcal.indata = contdata
    clcal.snver = 6
    clcal.gainver = 7
    clcal.gainuse = 8
    # clcal.sources =
    # clcal.calsour = AIPSList([source])
    clcal.interpol = 'SELF'
    clcal.go()

    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    vplot = AIPSTask('vplot')
    vplot.indata = contdata
    vplot.solint = 0.166
    vplot.bparm[2] = -1
    # vplot.avgchan = 1
    # vplot.avgif = 1
    vplot.dotv = -1
    vplot.refant = refant
    # vplot.nplots = 12
    # Also want to make bchan echan more automated
    vplot.bchan = 1
    vplot.echan = (contdata.header['naxis'][2])
    vplot.crowded = 3
    vplot.do3col = 1
    vplot.docal = 1
    vplot.gainuse = 8
    vplot.stokes = 'RR'
    vplot.go()
    vplot.stokes = 'LL'
    vplot.go()
    print("bchan echan is", vplot.bchan, vplot.echan)
    #######
    if os.path.exists('VPLOT_CL7.eps'):
        print("removing old VPLOT_CL7.eps")
        os.remove('VPLOT_CL7.eps')

    lwpla = AIPSTask('lwpla')
    lwpla.indata = contdata
    lwpla.plver = 1
    lwpla.inver = contdata.table_highver('AIPS PL')
    lwpla.outfile = 'DATAD:VPLOT_CL8.eps'
    lwpla.go()

    #### Plot CL8 POSSM
    while contdata.table_highver('AIPS PL') > 0:
        contdata.zap_table('AIPS PL', 0)
    for source in contdata.sources:
        possm = AIPSTask('possm')
        possm.indata = contdata
        possm.solint = 666
        possm.sources = AIPSList([source])
        possm.stokes = 'HALF'
        possm.aparm[1] = -1
        possm.flagver = 1
        possm.nplots = 2
        possm.dotv = -1
        possm.doband = 1
        possm.docal = 1
        possm.gainuse = 8
        possm.aparm[7] = 2
        possm.aparm[9] = 1
        possm.aparm[1] = -1
        # Fixed scale for phase
        possm.aparm[2] = 1
        possm.aparm[5] = -180
        possm.aparm[6] = 180
        possm.go()

        if os.path.exists('POSSM_CROSS_CL8.eps'):
            print("removing old POSSM_CROSS_CL8.eps")
            os.remove('POSSM_CROSS_CL8.eps')

        lwpla = AIPSTask('lwpla')
        lwpla.indata = contdata
        lwpla.plver = 1
        lwpla.inver = contdata.table_highver('AIPS PL')
        lwpla.outfile = 'DATAD:POSSM_CROSS_CL8.eps'
        lwpla.go()

if tmask <= 10:
    essentials.write("Reached tmask 10 \n")
    for source in contdata.sources:
        try:
            print('\n')
            print(source)
            #'''
            
            if os.path.exists('UVFIT_' + source + '.txt'):
                print("  Removing old UVFIT results text file")
                os.remove('UVFIT_' + source + '.txt')
            uvfit = AIPSTask('uvfit')
            uvfit.indata = contdata
            uvfit.srcname = AIPSList([source]).pop(1)
            uvfit.docal = 1
            uvfit.gainuse = 8
            uvfit.stokes = 'I'
            uvfit.prtlev = 1
            # uvfit.dopos = -1
            uvfit.solmode = ''
            uvfit.fitout = 'DATAD:UVFIT_' + source + '.txt'
            uvfit.go()
            #'''
        except:
            print("error")      

'''
#### Concatenate all the CL VPLOTS into one
if os.path.exists('VPLOTS_ALL.pdf'):
    print("removing old VPLOTS_ALL.pdf")
    os.remove('VPLOTS_ALL.pdf')

gs = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=VPLOTS_ALL.pdf VPLOT*.eps"

os.system(gs)
###########################################

#### Concatenate all the CL POSSMS into one
if os.path.exists('POSSM_CROSS_ALL.pdf'):
    print("removing old POSSM_CROSS_ALL.pdf")
    os.remove('POSSM_CROSS_ALL.pdf')

gs = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=POSSM_CROSS_ALL.pdf POSSM_CROSS*.eps"

os.system(gs)
###########################################


############ extract just the flux values from the results for essentials
essentials.write("Final Flux and error values are \n")
essentials.close()

awk = "awk '$3 ~ /\*/ {print FILENAME, $5,$6} FNR==2 && $3 !~ /\*/ {print FILENAME, $7,$8}' ./UVFIT_* cat >> essentials.txt"
os.system(awk)
'''
