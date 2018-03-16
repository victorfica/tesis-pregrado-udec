# -*- coding: utf-8 -*-
"""Functions to parse DSSP output. For DSSP, see U{http://swift.cmbi.ru.nl/gv/dssp/}.

The DSSP codes for secondary structure used here are:

    - H        Alpha helix (4-12)
    - B        Isolated beta-bridge residue
    - E        Strand
    - G        3-10 helix
    - I        pi helix
    - T        Turn
    - S        Bend
    - -        None
"""

import re

class DSSPData:
  def __init__(self):
      self.num      = []
      self.resnum   = []
      self.inscode  = []
      self.chain    = []
      self.aa       = []
      self.struct   = []
      self.bp1      = []
      self.bp2      = []
      self.pacc     = []
      self.acc      = []
      self.h_oinj   = None
      self.h_pb     = None
      self.h_apb    = None
      self.h_oinin5 = None
      self.h_oinin4 = None
      self.h_oinin3 = None
      self.h_oinin2 = None
      self.h_oinin1 = None
      self.h_oinip0 = None
      self.h_oinip2 = None
      self.h_oinip3 = None
      self.h_oinip4 = None
      self.h_oinip5 = None
      self.h_nho1   = []
      self.h_ohn1   = []
      self.h_nho2   = []
      self.h_ohn2   = []
      self.tco      = []
      self.kappa    = []
      self.alpha    = []
      self.phi      = []
      self.psi      = []
      self.xca      = []
      self.yca      = []
      self.zca      = []

  def parseDSSP(self, file):
    with open(file,"r") as in_handle1:

      lines = in_handle1.readlines()

      self.pacc.append(   float(lines[7].split()[0])  )
    # percentage/nmbr Hbonds types, hbonds in parallel/antiparallel bridges in DSSP output 
      self.h_oinj   =   lines[8].split()[0:2]         
      self.h_pb     =   lines[9].split()[0:2]
      self.h_apb    =   lines[10].split()[0:2]
      self.h_oinin5 =   lines[11].split()[0:2]
      self.h_oinin4 =   lines[12].split()[0:2]
      self.h_oinin3 =   lines[13].split()[0:2]
      self.h_oinin2 =   lines[14].split()[0:2]
      self.h_oinin1 =   lines[15].split()[0:2]
      self.h_oinip0 =   lines[16].split()[0:2]
      self.h_oinip1 =   lines[17].split()[0:2]
      self.h_oinip2 =   lines[18].split()[0:2]
      self.h_oinip3 =   lines[19].split()[0:2]
      self.h_oinip4 =   lines[20].split()[0:2]
      self.h_oinip5 =   lines[21].split()[0:2]



    with open(file, 'r') as input_handle:

    
      line_num = 0
      start=False
      for line in input_handle:
    
        if( re.search('#', line) ):
          start=True
          continue

        if( start ):
          self.num.append(     int(line[0:5].strip()) )
          self.resnum.append(  line[5:10].strip() )
          self.inscode.append( line[10:11].strip() )
          self.chain.append(   line[11:12].strip() )
          self.aa.append(      line[12:14].strip() )
          self.struct.append(  line[16] )
          self.bp1.append(     line[25:29].strip() )
          self.bp2.append(     line[29:34].strip() )
          self.acc.append(     int(line[34:38].strip()) )
          self.h_nho1.append(  line[38:50].strip() )
          self.h_ohn1.append(  line[50:61].strip() )
          self.h_nho2.append(  line[61:72].strip() )
          self.h_ohn2.append(  line[72:83].strip() )
          self.tco.append(     float(line[83:91].strip()) )
          self.kappa.append(   float(line[91:97].strip()) )
          self.alpha.append(   float(line[97:103].strip()) )
          self.phi.append(     float(line[103:109].strip()) )
          self.psi.append(     float(line[109:115].strip()) )
          self.xca.append(     float(line[115:122].strip()) )
          self.yca.append(     float(line[122:129].strip()) )
          self.zca.append(     float(line[129:136].strip()) )

    

  def getResnums(self):
    return self.resnum
  def getInsCode(self):
    return self.inscode
  def getChain(self):
    return self.chain
  def getAAs(self):
    return self.aa
  def getSecStruc(self):
    ss = ["-" if x==" " else x for x in self.struct]
    return ss 
  def getBP1(self):
    return self.bp1
  def getBP2(self):
    return self.bp2
  def getACC(self):
    return self.acc
  def getPACC(self):
    return self.pacc
  def getH_OINJ(self):
    h_oinj = [float(i) for i in self.h_oinj]
    return h_oinj
  def getH_Pbridges(self):
    h_pb = [float(i) for i in self.h_pb]
    return h_pb
  def getH_APbridges(self):
    h_apb = [float(i) for i in self.h_apb]
    return h_apb
  def getH_OININ5(self):
    h_oinin5 = [float(i) for i in self.h_oinin5]
    return h_oinin5
  def getH_OININ4(self):
    h_oinin4 = [float(i) for i in self.h_oinin4]
    return h_oinin4
  def getH_OININ3(self):
    h_oinin3 = [float(i) for i in self.h_oinin3]
    return h_oinin3
  def getH_OININ2(self):
    h_oinin2 = [float(i) for i in self.h_oinin2]
    return h_oinin2
  def getH_OININ1(self):
    h_oinin1 = [float(i) for i in self.h_oinin1]
    return h_oinin1
  def getH_OINIP0(self):
    h_oinip0 = [float(i) for i in self.h_oinip0]
    return h_oinip0
  def getH_OINIP1(self):
    h_oinip1 = [float(i) for i in self.h_oinip1]
    return h_oinip1
  def getH_OINIP2(self):
    h_oinip2 = [float(i) for i in self.h_oinip2]
    return h_oinip2
  def getH_OINIP3(self):
    h_oinip3 = [float(i) for i in self.h_oinip3]
    return h_oinip3
  def getH_OINIP4(self):
    h_oinip4 = [float(i) for i in self.h_oinip4]
    return h_oinip4 
  def getH_OINIP5(self):
    h_oinip5 = [float(i) for i in self.h_oinip5]
    return h_oinip5   
  def getH_NHO1(self):
    return self.h_nho1
  def getH_NHO2(self):
    return self.h_nho2
  def getH_OHN1(self):
    return self.h_ohn1
  def getH_OHN2(self):
    return self.h_ohn2
  def getTCO(self):
    return self.tco
  def getKAPPA(self):
    return self.kappa
  def getALPHA(self):
    return self.alpha
  def getPHI(self):
    return self.phi
  def getPSI(self):
    return self.psi
  def getX(self):
    return self.xca
  def getY(self):
    return self.yca
  def getZ(self):
    return self.zca



